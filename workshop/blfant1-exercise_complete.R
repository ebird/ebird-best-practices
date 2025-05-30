library(auk)
library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(mccf1)
library(ranger)
library(readr)
library(scam)
library(sf)
library(terra)
library(tidyr)

# set random number seed for reproducibility
set.seed(1)

# EXERCISE: Use the skills you learned in this workshop to train a relative
# abundance hurdle model for Black-faced Antthrush in Panama using eBird data
# from May-July of 2014-2023.


# Import eBird data ----

# checklist data
f_sed <- "data-raw/ebd_PA_blfant1_smp_relOct-2023_sampling.txt"
checklists <- read_sampling(f_sed)
# observation data
f_ebd <- "data-raw/ebd_PA_blfant1_smp_relOct-2023.txt"
observations <- read_ebd(f_ebd)


# Filter to region and season ----

# complete checklists from may-july of the last 10 years
# filter the checklist data
checklists <- checklists |>
  filter(all_species_reported,
         between(year(observation_date), 2014, 2023),
         between(month(observation_date), 5, 7))
# filter the observation data
observations <- observations |>
  filter(all_species_reported,
         between(year(observation_date), 2014, 2023),
         between(month(observation_date), 5, 7))

# subset to region boundary polygon to remove offshore checklists
# convert checklist locations to points geometries
checklists_sf <- checklists |>
  select(checklist_id, latitude, longitude) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
# boundary of study region, buffered by 1 km
study_region_buffered <- read_sf("data/gis-data.gpkg", layer = "ne_countries") |>
  filter(country_code == "PA") |>
  st_transform(crs = st_crs(checklists_sf)) |>
  st_buffer(dist = 1000)
# spatially subset the checklists to those in the study region
in_region <- checklists_sf[study_region_buffered, ]
# join to checklists and observations to remove checklists outside region
checklists <- semi_join(checklists, in_region, by = "checklist_id")
observations <- semi_join(observations, in_region, by = "checklist_id")

# remove observations without matching checklists
observations <- semi_join(observations, checklists, by = "checklist_id")


# Zero-fill eBird data ----

# combine checklist and observation data to produce detection/non-detection data
zf <- auk_zerofill(observations, checklists, collapse = TRUE)
# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}
# transform effort variables
# 1. converet counts to integer and "X" to NA
# 2. set distance to 0 for stationary checklists
# 3. convert duration to hours
# 4. create speed variable
# 5. convert time to hours since midnight
# 6. split date into year and day of year
zf <- zf |>
  mutate(
    # convert count to integer and X to NA
    # ignore the warning "NAs introduced by coercion"
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for stationary counts
    effort_distance_km = if_else(protocol_type == "Stationary",
                                 0, effort_distance_km),
    # convert duration to hours
    effort_hours = duration_minutes / 60,
    # speed km/h
    effort_speed_kmph = effort_distance_km / effort_hours,
    # convert time to decimal hours since midnight
    hours_of_day = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )


# Apply effort filters ----

# traveling or stationary counts with fewer than 10 observers
# duration <= 6 h, length <= 10 km, speed <= 100km/h
zf_filtered <- zf |>
  filter(protocol_type %in% c("Stationary", "Traveling"),
         effort_hours <= 6,
         effort_distance_km <= 10,
         effort_speed_kmph <= 100,
         number_observers <= 10)


# Test-train split ----

# split checklists into 20/80 test/train
zf_filtered$type <- if_else(runif(nrow(zf_filtered)) <= 0.8, "train", "test")

# subset to only those columns we need
checklists <- zf_filtered |>
  select(checklist_id, observer_id, type,
         observation_count, species_observed,
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         hours_of_day,
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers)
# save processed ebird data
write_csv(checklists, "data/checklists-zf_blfant1_may-jul_pa.csv", na = "")


# Observation map ----

# load gis data
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") |>
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") |>
  st_geometry()
study_region <- read_sf("data/gis-data.gpkg", "ne_countries") |>
  filter(country_code == "PA") |>
  st_geometry()

# prepare ebird data for mapping
checklists_sf <- checklists |>
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  select(species_observed)

# map
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
plot(st_geometry(checklists_sf), col = NA, border = NA)
# contextual gis data
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)
plot(study_region, col = "#e6e6e6", border = NA, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
# ebird observations
# not observed
plot(filter(checklists_sf, !species_observed),
     pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
     add = TRUE)
# observed
plot(filter(checklists_sf, species_observed),
     pch = 19, cex = 0.3, col = alpha("#4daf4a", 1),
     add = TRUE)
# legend
legend("bottomright", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("eBird checklist", "Detection"),
       pch = 19)
box()


# Environmental variables ----

# environmental variables: landcover and elevation
env_vars <- read_csv("data/environmental-variables_checklists_jan-dec_pa.csv")
# combine ebird and environmental data
checklists_env <- inner_join(checklists, env_vars, by = "checklist_id")


# Spatiotemporal subsampling ----

# sample one checklist per 3km x 3km x 1 week grid for each year
# sample detection/non-detection independently
checklists_ss <- grid_sample_stratified(checklists_env,
                                        obs_column = "species_observed",
                                        sample_by = "type")


# Encounter rate ----

# filter to training data, select only the columns to be used in the model
checklists_train <- checklists_ss |>
  filter(type == "train") |>
  select(species_observed, observation_count,
         year, day_of_year, hours_of_day,
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers,
         starts_with("pland_"),
         starts_with("ed_"),
         starts_with("elevation_"))

# calculate detection frequency
detection_freq <- mean(checklists_train$species_observed)

# train a balanced random forest to classify detection/non-detection
# remove observation_count prior to training model
train_er <- select(checklists_train, -observation_count)
er_model <- ranger(formula =  as.factor(species_observed) ~ .,
                   data = train_er,
                   importance = "impurity",
                   probability = TRUE,
                   replace = TRUE,
                   sample.fraction = c(detection_freq, detection_freq))


# ├ Calibration ----

# predicted encounter rate based on out of bag samples
er_pred <- er_model$predictions[, 2]
# observed detection, converted back from factor
det_obs <- as.integer(checklists_train$species_observed)
# construct a data frame to train the scam model
obs_pred <- data.frame(obs = det_obs, pred = er_pred)

# train calibration model
calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"),
                          gamma = 2,
                          data = obs_pred)

# calibration plot
# group the predicted encounter rate into bins of width 0.02
# then calculate the mean observed encounter rates in each bin
er_breaks <- seq(0, 1, by = 0.02)
mean_er <- obs_pred |>
  mutate(er_bin = cut(pred, breaks = er_breaks, include.lowest = TRUE)) |>
  group_by(er_bin) |>
  summarise(n_checklists = n(),
            pred = mean(pred),
            obs = mean(obs),
            .groups = "drop")
# make predictions from the calibration model
calibration_curve <- data.frame(pred = er_breaks)
cal_pred <- predict(calibration_model, calibration_curve, type = "response")
calibration_curve$calibrated <- cal_pred
# compared binned mean encounter rates to calibration model
ggplot(calibration_curve) +
  aes(x = pred, y = calibrated) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(color = "blue") +
  geom_point(data = mean_er,
             aes(x = pred, y = obs),
             size = 2, alpha = 0.6,
             show.legend = FALSE) +
  labs(x = "Estimated encounter rate",
       y = "Observed encounter rate",
       title = "Calibration model") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))


# ├ Thresholding ----

# mcc and fscore calculation for various thresholds
mcc_f1 <- mccf1(
  # observed detection/non-detection
  response = obs_pred$obs,
  # predicted encounter rate from random forest
  predictor = obs_pred$pred)
# identify best threshold
mcc_f1_summary <- summary(mcc_f1)
threshold <- mcc_f1_summary$best_threshold[1]


# ├ Prediction ----

# load prediction grid
pred_grid <- read_csv("data/environmental-variables_prediction-grid_pa.csv")
# raster template for the grid
r <- rast("data/prediction-grid_pa.tif")
# get the coordinate reference system of the prediction grid
crs <- st_crs(r)

# add standardized effort covariates to prediction grid
# 6am on 2023-06-15
# 2 km, 1 hour traveling checklist with 1 observer
pred_grid_eff <- pred_grid |>
  mutate(observation_date = ymd("2023-06-15"),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         hours_of_day = 6,
         effort_distance_km = 2,
         effort_hours = 1,
         effort_speed_kmph = 2,
         number_observers = 1)

# estimate encounter rate
pred_er <- predict(er_model, data = pred_grid_eff, type = "response")
pred_er <- pred_er$predictions[, 2]
# define range-boundary
pred_binary <- as.integer(pred_er > threshold)
# apply calibration
pred_er_cal <- predict(calibration_model,
                       data.frame(pred = pred_er),
                       type = "response") |>
  as.numeric()
# constrain to 0-1
pred_er_cal[pred_er_cal < 0] <- 0
pred_er_cal[pred_er_cal > 1] <- 1
# combine predictions with coordinates from prediction grid
predictions <- data.frame(cell_id = pred_grid_eff$cell_id,
                          x = pred_grid_eff$x,
                          y = pred_grid_eff$y,
                          in_range = pred_binary,
                          encounter_rate = pred_er_cal)

# rasterize predictions
r_pred <- predictions |>
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs) |>
  select(in_range, encounter_rate) |>
  # rasterize
  rasterize(r, field = c("in_range", "encounter_rate"))


# ├ Mapping ----

# project gis data to coordinate reference system of the prediction grid
study_region_proj<- st_transform(study_region, crs = crs)
ne_land_proj <- st_transform(ne_land, crs = crs)
ne_country_lines_proj <- st_transform(ne_country_lines, crs = crs)

# range boundary
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region_proj, col = NA, border = NA)
plot(ne_land_proj, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# convert binary prediction to categorical
r_range <- as.factor(r_pred[["in_range"]])
plot(r_range, col = c("#e6e6e6", "forestgreen"),
     maxpixels = ncell(r_range),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_country_lines_proj, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region_proj, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()


# encounter rate
par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region_proj, col = NA, border = NA)
plot(ne_land_proj, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks
brks <- global(r_pred[["encounter_rate"]], fun = quantile,
               probs = seq(0, 1, 0.1), na.rm = TRUE) |>
  as.numeric() |>
  unique()
# label the bottom, middle, and top value
lbls <- round(c(0, median(brks), max(brks)), 2)
# ebird status and trends color palette
pal <- ebirdst_palettes(length(brks) - 1)
plot(r_pred[["encounter_rate"]],
     col = pal, breaks = brks,
     maxpixels = ncell(r_pred),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_country_lines_proj, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region_proj, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Encounter Rate"
image.plot(zlim = c(0, 1), legend.only = TRUE,
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.04, 0.07),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))


# Count model ----

# subset to only observed or predicted detections
train_count <- checklists_train
train_count$pred_er <- er_model$predictions[, 2]
train_count <- train_count |>
  filter(!is.na(observation_count),
         observation_count > 0 | pred_er > threshold) |>
  select(-species_observed, -pred_er)

# add predicted encounter rate as an additional covariate
predicted_er <- predict(er_model, data = train_count, type = "response")
predicted_er <- predicted_er$predictions[, 2]
train_count$predicted_er <- predicted_er

# train a random forest to estimate expected count
count_model <- ranger(formula = observation_count ~ .,
                      data = train_count,
                      importance = "impurity",
                      replace = TRUE)


# ├ Prediction ----

# add predicted encounter rate required for count estimates
pred_grid_eff$predicted_er <- pred_er
# estimate count
pred_count <- predict(count_model, data = pred_grid_eff, type = "response")
pred_count <- pred_count$predictions
# combine with all other predictions
predictions$count <- pred_count
# relative abundance = encounter_rate * count
predictions$abundance <- predictions$encounter_rate * predictions$count

# rasterize
layers <- c("in_range", "encounter_rate", "count", "abundance")
r_pred <- predictions |>
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs) |>
  select(all_of(layers)) |>
  # rasterize
  rasterize(r, field = layers)


# ├ Mapping ----

# map of relative abundance
par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region_proj, col = NA, border = NA)
plot(ne_land_proj, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks, excluding zeros
brks <- ifel(r_pred[["abundance"]] > 0, r_pred[["abundance"]], NA) |>
  global(fun = quantile,
         probs = seq(0, 1, 0.1), na.rm = TRUE) |>
  as.numeric() |>
  unique()
# label the bottom, middle, and top value
lbls <- round(c(min(brks), median(brks), max(brks)), 2)
# ebird status and trends color palette
pal <- ebirdst_palettes(length(brks) - 1)
plot(r_pred[["abundance"]],
     col = c("#e6e6e6", pal), breaks = c(0, brks),
     maxpixels = ncell(r_pred),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_country_lines_proj, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region_proj, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Relative Abundance"
image.plot(zlim = c(0, 1), legend.only = TRUE,
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.04, 0.07),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))


# Assessment ----

# ├ Predict to test data ----

# get the test set held out from training
# only consider checklists with counts
checklists_test <- checklists_ss |>
  filter(type == "test", !is.na(observation_count)) |>
  mutate(species_observed = as.integer(species_observed))
# estimate encounter rate
pred_er <- predict(er_model, data = checklists_test, type = "response")
pred_er <- pred_er$predictions[, 2]
# convert predictions to binary (presence/absence) using the threshold
pred_binary <- as.integer(pred_er > threshold)
# calibrate
pred_calibrated <- predict(calibration_model,
                           newdata = data.frame(pred = pred_er),
                           type = "response") |>
  as.numeric()
# constrain probabilities to 0-1
pred_calibrated[pred_calibrated < 0] <- 0
pred_calibrated[pred_calibrated > 1] <- 1

# add predicted encounter rate required for count estimates
checklists_test$predicted_er <- pred_er
# estimate count
pred_count <- predict(count_model, data = checklists_test, type = "response")
pred_count <- pred_count$predictions
# relative abundance is the product of encounter rate and count
pred_abundance <- pred_calibrated * pred_count

# combine observations and estimates
obs_pred_test <- data.frame(
  id = seq_along(pred_abundance),
  # actual detection/non-detection
  obs_detected = as.integer(checklists_test$species_observed),
  obs_count = checklists_test$observation_count,
  # model estimates
  pred_binary = pred_binary,
  pred_er = pred_calibrated,
  pred_count = pred_count,
  pred_abundance = pred_abundance
)


# ├ Encounter rate PPMs ----

# mean squared error (mse)
mse <- mean((obs_pred_test$obs_detected - obs_pred_test$pred_er)^2)

# precision-recall auc
em <- precrec::evalmod(scores = obs_pred_test$pred_binary,
                       labels = obs_pred_test$obs_detected)
pr_auc <- precrec::auc(em) |>
  filter(curvetypes == "PRC") |>
  pull(aucs)

# calculate metrics for binary prediction: sensitivity, specificity
pa_metrics <- obs_pred_test |>
  select(id, obs_detected, pred_binary) |>
  PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)

# combine ppms together
er_ppms <- data.frame(
  mse = mse,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity,
  pr_auc = pr_auc
)


# ├ Count PPMs ----

# subset to only those checklists where detect occurred
detections_test <- filter(obs_pred_test, obs_detected > 0)

# count metrics, based only on checklists where detect occurred
count_spearman <- cor(detections_test$pred_count,
                      detections_test$obs_count,
                      method = "spearman")
log_count_pearson <- cor(log(detections_test$pred_count + 1),
                         log(detections_test$obs_count + 1),
                         method = "pearson")

# abundance metrics, based on all checklists
abundance_spearman <- cor(obs_pred_test$pred_abundance,
                          obs_pred_test$obs_count,
                          method = "spearman")
log_abundance_pearson <- cor(log(obs_pred_test$pred_abundance + 1),
                             log(obs_pred_test$obs_count + 1),
                             method = "pearson")

# combine ppms together
count_abd_ppms <- data.frame(
  count_spearman = count_spearman,
  log_count_pearson = log_count_pearson,
  abundance_spearman = abundance_spearman,
  log_abundance_pearson = log_abundance_pearson
)


# Habitat associations ----

# ├ Predictor importance ----

# extract partial dependence from the random forest model objects
# encounter rate
pi_er <- er_model$variable.importance
pi_er <- data.frame(predictor = names(pi_er), importance = pi_er) |>
  arrange(desc(importance))
# count
pi_count <- count_model$variable.importance
pi_count <- data.frame(predictor = names(pi_count), importance = pi_count) |>
  arrange(desc(importance))
# plot predictor importance for top 20 encounter rate predictors
gg_er <- ggplot(head(pi_er, 20)) +
  aes(x = reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL,
       y = "Predictor Importance (Gini Index)",
       title = "Predictor importance for encounter rate model") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", linewidth = 0.5))
# plot predictor importance for top 20 count predictors
gg_count <- ggplot(head(pi_count, 20)) +
  aes(x = reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL,
       y = "Predictor Importance (Gini Index)",
       title = "Predictor importance for count model") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", linewidth = 0.5))
grid.arrange(gg_er, gg_count, nrow = 1)


# ├ Partial dependence ----

# function to calculate partial dependence for a given predictor
calculate_pd <- function(predictor,
                         er_model, calibration_model, count_model,
                         data, x_res = 25, n = 1000) {
  # create prediction grid using quantiles
  x_grid <- quantile(data[[predictor]],
                     probs = seq(from = 0, to = 1, length = x_res),
                     na.rm = TRUE)
  # remove duplicates
  x_grid <- x_grid[!duplicated(signif(x_grid, 8))]
  x_grid <- unname(unique(x_grid))
  grid <- data.frame(predictor = predictor, x = x_grid)
  names(grid) <- c("predictor", predictor)

  # subsample training data
  n <- min(n, nrow(data))
  data <- data[sample(seq.int(nrow(data)), size = n, replace = FALSE), ]

  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)

  # estimate encounter rate
  pred_er <- predict(er_model, data = grid)
  # estimate count
  grid$predicted_er <- pred_er$predictions[, 2]
  pred_count <- predict(count_model, data = grid)

  # summarize
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "x")
  pd$encounter_rate <- pred_er$predictions[, 2]
  pd$count <- pred_count$predictions
  pd <- dplyr::group_by(pd, predictor, x)
  pd <- dplyr::summarise(pd,
                         encounter_rate = mean(encounter_rate, na.rm = TRUE),
                         count = mean(count, na.rm = TRUE),
                         .groups = "drop")

  # calibrate
  pd$encounter_rate <- predict(calibration_model,
                               newdata = data.frame(pred = pd$encounter_rate),
                               type = "response")
  pd$encounter_rate <- as.numeric(pd$encounter_rate)
  # constrain to 0-1
  pd$encounter_rate[pd$encounter_rate < 0] <- 0
  pd$encounter_rate[pd$encounter_rate > 1] <- 1

  # calculate relative abundance
  pd$abundance <- pd$encounter_rate * pd$count

  return(pd)
}

# make a partial dependence plot of hours_of_day
# assess whether we choose the optimal time for the standard checklist
pd_duration <- calculate_pd("hours_of_day",
                            er_model = er_model,
                            calibration_model = calibration_model,
                            count_model = count_model,
                            data = checklists_train)
pd_duration <- pivot_longer(pd_duration,
                            cols = c(encounter_rate, count, abundance))
ggplot(pd_duration) +
  aes(x = x, y = value) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 24, by = 1)) +
  facet_wrap(~ name, ncol = 1, scales = "free") +
  labs(x = "Checklist start time [hours since midnight]",
       y = "Encounter Rate")

# calculate abundance partial dependence for each of the top 6 predictors
pd <- NULL
for (predictor in head(pi_er$predictor)) {
  pd <- calculate_pd(predictor,
                     er_model = er_model,
                     calibration_model = calibration_model,
                     count_model = count_model,
                     data = checklists_train) |>
    bind_rows(pd)
}
# plot partial dependence
ggplot(pd) +
  aes(x = x, y = abundance) +
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor, levels = rev(unique(predictor))),
             ncol = 2, scales = "free") +
  labs(x = NULL, y = "Relative Abundance") +
  theme_minimal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))
