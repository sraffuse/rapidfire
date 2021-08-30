
# A script to test processing

library(rapidfire)
library(ggplot2)

# Can we successfully process two full months?
dt1 <- as.Date("2020-09-01")
dt2 <- as.Date("2020-10-31")

# Get and prep AirNow data
an_ws <- get_airnow_daterange(dt1, dt2, "CA")
an <- recast_monitors(an_ws)
an_vg <- create_airnow_variograms(an)
an_data <- split_airnow_data(an)
an_ok <- krige_airnow(an_data, an_vg)

# Get and prep PurpleAir data
pa <- get_purpleair_daterange(dt1, dt2, "CA")
pas <- purpleair_spatial(pa)
pa_clean <- purpleair_clean_spatial_outliers(pas)
pa_vg <- create_purpleair_variograms(pa_clean, cutoff = 100000)

# Krige purple air data for all locations in the an data set
pa_ok <- krige_purpleair(pa_clean, an, pa_vg)

compare <- pa_ok@data %>%
  mutate(PM25_purp = exp(purp.pred))

ggplot(compare, aes(x = PM25_log, y = purp.pred)) + geom_point() +
  geom_abline(slope = 1)
ggplot(compare, aes(x = PM25, y = PM25_purp)) + geom_point()

# NARR Data
# only need to run this one once
# get_narr(2020, out_path = "./data/NARR")
narr <- narr_at_airnow(an)


# Bluesky Data
# only need to run this once
# grab_bluesky(dt1, dt2, output_path = "./data/bluesky/NAM/", model = "NAM-3km")
bluesky <- bluesky_at_airnow(an)

ggplot(bluesky, aes(x = PM25_log, y = bluesky_log, color = Day)) + geom_point()
ggplot(bluesky, aes(x = PM25, y = PM25_bluesky, color = Day)) + geom_point() +
  scale_y_continuous(limits = c(0, 200), oob = scales::squish) +
  scale_x_continuous(limits = c(0, 200), oob = scales::squish)


# Now need to put altogether in prep for RF modeling
train <- an_data$train@data %>%
  mutate(Split = "Train")
test <- an_data$test@data %>%
  mutate(Split = "Test")

# These have the measured PM25 - the target variable
model_in <- bind_rows(train, test)

# All the inputs
# AirNow Kriged
ank_in <- an_ok@data %>%
  select(monitorID, Day, PM25_log_ANK = var1.pred)
# PurpleAir Kriged
pak_in <- pa_ok@data %>%
  select(monitorID, Day, PM25_log_PAK = purp.pred)
# NARR
narr_in <- narr %>%
  select(-PM25, -Hours, -PM25_log)
# Bluesky
bluesky_in <- bluesky %>%
  select(monitorID, Day, PM25_bluesky)

model_in <- model_in %>%
  left_join(ank_in, by = c("monitorID", "Day")) %>%
  left_join(pak_in, by = c("monitorID", "Day")) %>%
  left_join(narr_in, by = c("monitorID", "Day")) %>%
  left_join(bluesky_in, by = c("monitorID", "Day"))

# Because of two layers of modeling, I have test and training flipped here - now
# use the test data to train the RF model

library(randomForest)
library(caret)

set.seed(1977)
train_control <- trainControl(method = "cv", number = 10)
tune_grid <- data.frame(mtry = c(2,3,4,5))

model_in_test <- filter(model_in, Split == "Test")

mod1.cv <- train(PM25_log ~ PM25_log_ANK + PM25_log_PAK + PM25_bluesky +
                   air.2m + uwnd.10m + vwnd.10m + rhum.2m + apcp + hpbl,
                 data = model_in_test, tuneGrid = tune_grid, do.trace = 100,
                 ntree = 500, trControl = train_control, method = "rf",
                 importance = TRUE)
print(mod1.cv)
# save these model results
saveRDS(mod1.cv, "./data/models/mod1cv.RDS")

model_in_test$PM25_log_RF <- predict(mod1.cv$finalModel, model_in_test)

ggplot(model_in_test, aes(x = PM25_log, y = PM25_log_RF)) + geom_point()

model_in_test$PM25_RF <- exp(model_in_test$PM25_log_RF)
ggplot(model_in_test, aes(x = PM25, y = PM25_RF, color = Day)) + geom_point()

varImpPlot(mod1.cv$finalModel)
# PurpleAir adds a lot, the met data just a little, the model almost none


