library(auk)
library(arrow)
library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(glue)
library(gridExtra)
library(lubridate)
library(ranger)
library(rnaturalearth)
library(scam)
library(sf)
library(terra)
library(tidyr)

# set random number seed for reproducibility
set.seed(1)

# define species, region, and season
species <- "taibar1"
species_name <- ebird_species(species, type = "common")


# eBird data ----

# checklist data
checklists_all <- read_parquet("data/erd_checklists_1km_tw.paquet")
glimpse(checklists_all)

# observation data
observations_all <- read_parquet("data/erd_observations_tw.parquet")
glimpse(observations_all)
# notice that some observations have counts of NA, these indicate the species
# was observed but no count was reported, i.e. an "X" on the checklist
filter(observations_all, is.na(observation_count))

# join checklists and observations and "zero-fill"
observations <- observations_all |>
  filter(species_code == species) |>
  # add a binary variable for detection
  mutate(species_observed = TRUE)
zf <- left_join(checklists_all, observations, by = "checklist_id") |>
  mutate(species_observed = coalesce(species_observed, FALSE),
         observation_count = ifelse(species_observed, observation_count, 0))
# number of detections
sum(zf$species_observed)


# ├ Apply effort filters ----

# this dataset has been pre-filtered for 3 km status and trends modeling:
# - traveling or stationary checklists on land
# - fewer than 50 observers
# - duration between 2 minutes and 8 hours
# - distance traveled <= 10 km and speed <= 100km/h

# histogram of remaining variation in distance traveled
ggplot(zf) +
  aes(x = effort_distance_km) +
  geom_histogram(binwidth = 0.5,
                 aes(y = after_stat(count / sum(count)))) +
  scale_y_continuous(limits = c(0, NA), labels = scales::label_percent()) +
  labs(x = "Distance traveled [km]",
       y = "% of eBird checklists",
       title = "Distribution of distance traveled")

# apply effort filters suitable for 1 km modeling
# also filter to observations from the main island from the last 10 years
region <- c(lat_min = 21.7, lat_max = 25.4,
            lon_min = 119.1, lon_max = 122.4)
zf_filtered <- zf |>
  filter(between(latitude, region["lat_min"], region["lat_max"]),
         between(longitude, region["lon_min"], region["lon_max"]),
         between(year, 2015, 2024),
         effort_hours <= 8,
         effort_distance_km <= 1.5,
         effort_speed_kmph <= 100,
         number_observers <= 10)

# number of checklists before and after filtering
nrow(zf_filtered) / nrow(zf)


# ├ Test-train split ----

# split checklists into 20/80 test/train
zf_filtered$type <- if_else(runif(nrow(zf_filtered)) <= 0.8, "train", "test")
table(zf_filtered$type) / nrow(zf_filtered)


# ├ Mapping ----

# country boundaries for context
countries <- ne_countries(scale = 10, continent = "Asia") |>
  select(name)

# prepare ebird data for mapping
checklists_sf <- zf_filtered |>
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  select(species_observed)

# map
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
plot(st_geometry(checklists_sf),
     #main = glue("{species_name} eBird observations\n 2015-2024"),
     col = NA, border = NA)
# country boundaries
plot(countries, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# ebird observations
# not observed
plot(filter(checklists_sf, !species_observed),
     pch = 19, cex = 0.05, col = alpha("#555555", 0.1),
     add = TRUE)
# observed
plot(filter(checklists_sf, species_observed),
     pch = 19, cex = 0.15, col = alpha("#4daf4a", 0.5),
     add = TRUE)
# legend
legend("topleft", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("eBird checklist", "Detections"),
       pch = 19)
box()


# ├ Spatiotemporal subsampling ----

# sample one checklist per 1 km x 1 km x 1 week grid for each year
# sample detection/non-detection independently
# stratify by type (test/train)
checklists_ss <- grid_sample_stratified(zf_filtered,
                                        res = c(1000, 1000, 7),
                                        obs_column = "species_observed",
                                        by_year = TRUE,
                                        case_control = TRUE,
                                        sample_by = "type")

# explore impact of sampling
# original data
nrow(zf_filtered)
count(zf_filtered, species_observed) |>
  mutate(percent = n / sum(n))
# after sampling
nrow(checklists_ss)
count(checklists_ss, species_observed) |>
  mutate(percent = n / sum(n))


# Prediction surface ----

# load the prediction surface environmental variables
pred_grid <- read_parquet("data/prediction-grid_1km_tw.paquet") |>
  filter(between(latitude, region["lat_min"], region["lat_max"]),
         between(longitude, region["lon_min"], region["lon_max"]))

# load the raster template for the grid
r <- rast("data/prediction-grid_1km_tw.tif") |>
  rast()

# EXERCISE: explore the environmental variables. For details visit
# https://docs.google.com/spreadsheets/d/1-_fmNGraFkkiLdycrQFYOz_j7OPBhG0hv4xKA_HPtaE

# insert elevation values into the raster
elevation <- pred_grid |>
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs(r)) |>
  # rasterize points
  rasterize(r, field = "elevation_30m_median")

# make a map of elevation
plot(elevation,
     axes = FALSE, box = FALSE, col = viridis(10),
     main = "Elevation [m]")


# Relative abundance model ----

# filter to training data, select only the columns to be used in the model
checklists_train <- checklists_ss |>
  filter(type == "train") |>
  select(species_observed, observation_count,
         year, day_of_year, hours_of_day,
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers,
         eastness_90m_median, eastness_90m_sd,
         northness_90m_median, northness_90m_sd,
         elevation_30m_median, elevation_30m_sd,
         starts_with("mcd"),
         astwbd_c1_ed, astwbd_c1_pland,
         astwbd_c2_ed, astwbd_c2_pland,
         gsw_c2_pland, gsw_c2_ed,
         starts_with("lake"),
         starts_with("river"),
         starts_with("road"),
         ntl_mean, ntl_sd)

# calculate detection frequency
detection_freq <- mean(checklists_train$species_observed)
n_detections <- sum(checklists_train$species_observed)


# ├ Range ----

# train a balanced binary classification forest for range boundary
# remove observation_count prior to training model
train_er <- select(checklists_train, -observation_count)
train_er$species_observed <- factor(train_er$species_observed)
rf_range <- ranger(
  formula = species_observed ~ .,
  data = train_er,
  importance = "impurity",
  sample.fraction = c(detection_freq, detection_freq)
)


# ├ Encounter rate ----

# train a balanced probability forest for encounter rate
rf_er <- ranger(
  formula = species_observed ~ .,
  data = train_er,
  importance = "impurity",
  probability = TRUE,
  sample.fraction = c(detection_freq, detection_freq),
  min.node.size = ceiling(n_detections * 0.02)
)


# ├ Calibration ----

# predicted encounter rate based on out of bag samples
er_pred <- rf_er$predictions[, "TRUE"]
# observed detection, converted back from factor
det_obs <- as.integer(checklists_train$species_observed == "TRUE")
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


# ├ Count ----

# hurdle: subset to detections only for the count model
train_count <- checklists_train |>
  filter(!is.na(observation_count), species_observed) |>
  select(-species_observed)

# plugin: add predicted encounter rate as an additional covariate
predicted_er <- predict(rf_er, data = train_count, type = "response")
predicted_er <- predicted_er$predictions[, "TRUE"]
train_count$predicted_er <- predicted_er

# train a regression forest for count
rf_count <- ranger(
  formula = observation_count ~ .,
  data = train_count,
  importance = "impurity"
)


# Prediction ----

# add standardized effort covariates to prediction grid
# 6am on 2024-05-15
# 2 km, 1 hour traveling checklist with 1 observer
pred_grid_eff <- pred_grid |>
  mutate(observation_date = ymd("2024-05-15"),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         hours_of_day = 6,
         effort_distance_km = 2,
         effort_hours = 1,
         effort_speed_kmph = 2,
         number_observers = 1)

# estimate range
pred_range <- predict(rf_range, data = pred_grid_eff, type = "response")
pred_range <- as.integer(pred_range$predictions == "TRUE")
# estimate encounter rate
pred_er <- predict(rf_er, data = pred_grid_eff, type = "response")
pred_er <- pred_er$predictions[, "TRUE"]
# apply calibration
pred_er_cal <- predict(calibration_model,
                       data.frame(pred = pred_er),
                       type = "response") |>
  as.numeric()
# constrain to 0-1
pred_er_cal[pred_er_cal < 0] <- 0
pred_er_cal[pred_er_cal > 1] <- 1
# estimate count
pred_grid_eff$predicted_er <- pred_er
pred_count <- predict(rf_count, data = pred_grid_eff, type = "response")
pred_count <- pred_count$predictions

# combine predictions with cell id and coordinates from prediction grid
predictions <- data.frame(cell_id = pred_grid_eff$cell_id,
                          x = pred_grid_eff$x,
                          y = pred_grid_eff$y,
                          range = pred_range,
                          encounter_rate = pred_er_cal,
                          count = pred_count,
                          abundance = pred_range * pred_er_cal * pred_count)


# ├ Rasterize predictions ----

# insert predictions into the raster template
r_pred <- predictions |>
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = st_crs(r)) |>
  # rasterize
  rasterize(r, field = c("range", "encounter_rate", "count", "abundance"))


# ├ Mapping ----

countries_proj <- st_transform(countries, st_crs(r))
tw_boundary <- filter(countries_proj, name == "Taiwan") |>
  st_geometry()

# map encounter rate
par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(tw_boundary, col = NA, border = NA)
plot(countries_proj, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks, excluding zeros
brks <- ifel(r_pred[["encounter_rate"]] > 0, r_pred[["encounter_rate"]], NA) |>
  global(fun = quantile,
         probs = seq(0, 1, 0.1), na.rm = TRUE) |>
  as.numeric() |>
  unique()
# label the bottom, middle, and top value
lbls <- round(c(min(brks), median(brks), max(brks)), 3)
# ebird status and trends color palette
pal <- ebirdst_palettes(length(brks) - 1)
plot(r_pred[["encounter_rate"]],
     col = c("#e6e6e6", pal), breaks = c(0, brks),
     maxpixels = ncell(r_pred),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- glue("{species_name} encounter rate (Oct 2024)")
image.plot(zlim = c(0, 1), legend.only = TRUE,
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.15, 0.85, 0.05, 0.08),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0.1))

# map relative abundance in range
par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(tw_boundary, col = NA, border = NA)
plot(countries_proj, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

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
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- glue("{species_name} relative abundance (Oct 2024)")
image.plot(zlim = c(0, 1), legend.only = TRUE,
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.15, 0.85, 0.05, 0.08),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0.1))


# Assessment ----

# ├ Predict to test data ----

# get the test set held out from training
# only consider checklists with counts
checklists_test <- checklists_ss |>
  filter(type == "test", !is.na(observation_count)) |>
  mutate(species_observed = as.integer(species_observed))

# estimate range
pred_range <- predict(rf_range, data = checklists_test, type = "response")
pred_range <- as.integer(pred_range$predictions == "TRUE")
# estimate encounter rate
pred_er <- predict(rf_er, data = checklists_test, type = "response")
pred_er <- pred_er$predictions[, "TRUE"]
# apply calibration
pred_er_cal <- predict(calibration_model,
                       data.frame(pred = pred_er),
                       type = "response") |>
  as.numeric()
# constrain to 0-1
pred_er_cal[pred_er_cal < 0] <- 0
pred_er_cal[pred_er_cal > 1] <- 1
# add predicted encounter rate required for count estimates
checklists_test$predicted_er <- pred_er
# estimate count
pred_count <- predict(rf_count, data = checklists_test, type = "response")
pred_count <- pred_count$predictions

# combine observations and estimates
obs_pred_test <- data.frame(
  checklist_id = checklists_test$checklist_id,
  # actual detection/non-detection
  obs_detected = checklists_test$species_observed,
  obs_count = checklists_test$observation_count,
  # model estimates
  pred_binary = pred_range,
  pred_er = pred_er_cal,
  pred_count = pred_count,
  pred_abundance = pred_er_cal * pred_count
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
  select(checklist_id, obs_detected, pred_binary) |>
  PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)

# combine ppms together
er_ppms <- data.frame(
  mse = mse,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity,
  pr_auc = pr_auc
)
print(er_ppms)


# ├ Count PPMs ----

# subset to only those checklists where detect occurred
detections_test <- filter(obs_pred_test, obs_detected > 0)

# count metrics, based only on checklists where detect occurred
count_spearman <- cor(detections_test$pred_count,
                      detections_test$obs_count,
                      method = "spearman")
log_count_pearson <- cor(log(detections_test$pred_count),
                         log(detections_test$obs_count),
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
print(count_abd_ppms)


# Habitat associations ----

# ├ Predictor importance ----

# extract partial dependence from the encounter rate random forest
pi_er <- rf_er$variable.importance
pi_er <- data.frame(predictor = names(pi_er), importance = pi_er) |>
  # scale so importances sum to 1
  mutate(importance = importance / sum(importance)) |>
  arrange(desc(importance))
# plot predictor importance for top 15 encounter rate predictors
ggplot(head(pi_er, 15)) +
  aes(x = reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL,
       y = "Predictor Importance",
       title = "Predictor importance for encounter rate model") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", linewidth = 0.5))


# ├ Partial dependence ----

# function to calculate partial dependence for a given predictor
calculate_pd <- function(predictor,
                         er_model, calibration_model,
                         data, x_res = 25, n = 1000) {
  # create prediction grid using quantiles based on detections
  detections <- data[data$species_observed, ]
  x_grid <- quantile(detections[[predictor]],
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
  pred_er <- predict(er_model, data = grid, type = "response")
  pred_er <- pred_er$predictions[, "TRUE"]
  # apply calibration
  pred_er_cal <- predict(calibration_model,
                         data.frame(pred = pred_er),
                         type = "response") |>
    as.numeric()
  # constrain to 0-1
  pred_er_cal[pred_er_cal < 0] <- 0
  pred_er_cal[pred_er_cal > 1] <- 1

  # summarize
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "x")
  pd$encounter_rate <- pred_er_cal
  pd <- dplyr::group_by(pd, predictor, x)
  pd <- dplyr::summarise(pd,
                         encounter_rate = mean(encounter_rate, na.rm = TRUE),
                         .groups = "drop")

  return(pd)
}

# calculate partial dependence for elevation
pd_elevation <- calculate_pd("elevation_30m_median",
                             er_model = rf_er,
                             calibration_model = calibration_model,
                             data = checklists_train)
# plot partial dependence
ggplot(pd_elevation) +
  aes(x = x, y = encounter_rate) +
  geom_line() +
  geom_point() +
  labs(x = "Elevation [m]", y = "Encounter Rate")

# calculate abundance partial dependence for each of the top 6 predictors
pd <- NULL
for (predictor in head(pi_er$predictor)) {
  pd <- calculate_pd(predictor,
                     er_model = rf_er,
                     calibration_model = calibration_model,
                     data = checklists_train) |>
    bind_rows(pd)
}
# plot partial dependence
ggplot(pd) +
  aes(x = x, y = encounter_rate) +
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor, levels = rev(unique(predictor))),
             ncol = 2, scales = "free") +
  labs(x = NULL, y = "Encouter Rate") +
  theme_minimal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

# EXERCISE: Produce the partial dependence plot for the checklist start time.
# How can this be used to help choose an optimal time of day for the
# standardized checklist that we predicted to prior to mapping?
pd_time <- calculate_pd("hours_of_day",
                        er_model = rf_er,
                        calibration_model = calibration_model,
                        data = checklists_train)
ggplot(pd_time) +
  aes(x = x, y = encounter_rate) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 24, by = 1)) +
  labs(x = "Checklist start time [hours since midnight]",
       y = "Encounter Rate")
