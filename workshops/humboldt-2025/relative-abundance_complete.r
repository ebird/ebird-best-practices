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
library(readr)
library(scam)
library(sf)
library(terra)
library(tidyr)

# set random number seed for reproducibility
set.seed(1)


# Load data ----

species <- "pnbfin1"
species_name <- ebird_species(species, type = "scientific")

# environmental variables
env_vars <- read_parquet("data/environmental-variables_checklists_co.parquet")
# zero-filled ebird data
checklists <- read_csv(glue("data/checklists-zf_{species}_co.csv"))
# combine ebird data and environmental variables
checklists_env <- inner_join(checklists, env_vars,
                             by = c("locality_id", "year"))

# prediction grid
pred_grid <- read_parquet("data/prediction-grid_co.parquet")

# convert mountain variable to factor
checklists_env$mountain <- factor(checklists_env$mountain)
# use the same levels in the prediction grid
pred_grid$mountain <- factor(pred_grid$mountain,
                             levels = levels(checklists_env$mountain))
pred_grid <- filter(pred_grid, !is.na(mountain))

# raster template for the grid
r <- rast("data/prediction-grid_co.tif") |>
  rast()
# get the coordinate reference system of the prediction grid
crs <- st_crs(r)

# load gis data for making maps
land <- read_sf("data/gis-data.gpkg", "land") |>
  st_transform(crs = crs) |>
  st_geometry()
country_lines <- read_sf("data/gis-data.gpkg", "country_lines") |>
  st_transform(crs = crs) |>
  st_geometry()
region_boundary <- read_sf("data/gis-data.gpkg", "region") |>
  st_transform(crs = crs) |>
  st_geometry()

# subset prediction grid to region
pred_grid_sf <- st_as_sf(pred_grid,
                         coords = c("x", "y"),
                         crs = crs,
                         remove = FALSE)
pred_grid <- pred_grid_sf[region_boundary, ] |>
  st_drop_geometry()


# Spatiotemporal subsampling ----

# sample one checklist per 3km x 3km x 1 week grid for each year
# sample detection/non-detection independently
# stratify by factor variables: mountain and type
checklists_ss <- grid_sample_stratified(checklists_env,
                                        obs_column = "species_observed",
                                        by_year = TRUE,
                                        case_control = TRUE,
                                        sample_by = c("type", "mountain"))

# original data
nrow(checklists_env)
count(checklists_env, species_observed) |>
  mutate(percent = n / sum(n))
# after sampling
nrow(checklists_ss)
count(checklists_ss, species_observed) |>
  mutate(percent = n / sum(n))


# Relative abundance model ----

# filter to training data, select only the columns to be used in the model
checklists_train <- checklists_ss |>
  filter(type == "train") |>
  select(species_observed, observation_count,
         year, day_of_year, hours_of_day,
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers,
         elevation_30m_median, elevation_30m_sd,
         northness_90m_median, northness_90m_sd,
         eastness_90m_median, eastness_90m_sd,
         mountain,
         starts_with("landcover"),
         starts_with("road"))
glimpse(checklists_train)

# calculate detection frequency and number of detections
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
  respect.unordered.factors = TRUE,
  replace = TRUE,
  sample.fraction = c(detection_freq, detection_freq),
  min.node.size = 1
)


# ├ Encounter rate ----

# train a balanced probability forest for encounter rate
rf_er <- ranger(
  formula = species_observed ~ .,
  data = train_er,
  importance = "impurity",
  respect.unordered.factors = TRUE,
  probability = TRUE,
  replace = TRUE,
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
  importance = "impurity",
  respect.unordered.factors = TRUE
)


# Prediction ----

# add standardized effort covariates to prediction grid
# 6:30am on 2024-10-15
# 2 km, 1 hour traveling checklist with 1 observer
pred_grid_eff <- pred_grid |>
  mutate(observation_date = ymd("2024-10-15"),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         hours_of_day = 6.5,
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
                          abundance = pred_er_cal * pred_count)
# abundance within range
predictions$abd_range <- predictions$abundance * predictions$range


# ├ Rasterize predictions ----

# insert predictions into the raster template
r_pred <- predictions |>
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs) |>
  # rasterize
  rasterize(r, field = c("range", "encounter_rate", "count",
                         "abundance", "abd_range"))


# ├ Mapping ----

# map encounter rate
par(mar = c(3, 0.25, 0.25, 0.25))
# set up plot area
plot(region_boundary, col = NA, border = NA)
plot(land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks
brks <- global(r_pred[["encounter_rate"]], fun = quantile,
               probs = seq(0, 1, 0.1), na.rm = TRUE) |>
  as.numeric() |>
  unique()
# label the bottom, middle, and top value
lbls <- round(c(min(brks), median(brks), max(brks)), 2)
# ebird status and trends color palette
pal <- ebirdst_palettes(length(brks) - 1)
plot(r_pred[["encounter_rate"]],
     col = pal, breaks = brks,
     maxpixels = ncell(r_pred),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(region_boundary, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- glue("{species_name} encounter rate (Oct 2024)")
image.plot(zlim = c(0, 1), legend.only = TRUE,
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.09, 0.12),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0.1))


# map relative abundance in range
par(mar = c(3, 0.25, 0.25, 0.25))
# set up plot area
plot(region_boundary, col = NA, border = NA)
plot(land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks, excluding zeros
brks <- ifel(r_pred[["abd_range"]] > 0, r_pred[["abd_range"]], NA) |>
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
plot(country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(region_boundary, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- glue("{species_name} relative abundance (Oct 2024)")
image.plot(zlim = c(0, 1), legend.only = TRUE,
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.09, 0.12),
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

# extract predictor importance from the random forest model objects
# encounter rate
pi_er <- rf_er$variable.importance
pi_er <- data.frame(predictor = names(pi_er), importance = pi_er) |>
  # scale so importances sum to 1
  mutate(importance = importance / sum(importance)) |>
  arrange(desc(importance))
# count
pi_count <- rf_count$variable.importance
pi_count <- data.frame(predictor = names(pi_count), importance = pi_count) |>
  # scale so importances sum to 1
  mutate(importance = importance / sum(importance)) |>
  arrange(desc(importance))
# plot predictor importance for top 20 encounter rate predictors
gg_er <- ggplot(head(pi_er, 20)) +
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
# plot predictor importance for top 20 count predictors
gg_count <- ggplot(head(pi_count, 20)) +
  aes(x = reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL,
       y = "Predictor Importance",
       title = "Predictor importance for count model") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", linewidth = 0.5))
grid.arrange(gg_er, gg_count, nrow = 1)


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
# exclude mountain since it is categorical
top6 <- head(setdiff(pi_er$predictor, "mountain"))
pd <- NULL
for (predictor in top6) {
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
