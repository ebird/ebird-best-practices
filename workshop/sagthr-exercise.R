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
# abundance hurdle model for Sage Thrasher in Wyoming in June-July using eBird
# data from 2014-2023.


# Import eBird data ----

# download checklist and observation data via the custom download form
# move both text files into the data-raw/ subdirectory of your rstudio project
# load both datasets



# Filter to region and season ----

# complete checklists in june-july from the last 10 years



# Zero-fill eBird data ----

# combine checklist and observation data to produce detection/non-detection data

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


# Apply effort filters ----

# traveling or stationary counts with fewer than 10 observers
# duration <= 6 h, length <= 10 km, speed <= 100km/h



# Test-train split ----

# split checklists into 20/80 test/train


# subset to only those columns we need



# Observation map ----

# load gis data
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") |>
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") |>
  st_geometry()
ne_state_lines <- read_sf("data/gis-data.gpkg", "ne_state_lines") |>
  st_geometry()
study_region <- read_sf("data/gis-data.gpkg", "ne_states") |>
  filter(state_code == "US-WY") |>
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
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
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
env_vars <- read_csv("data/environmental-variables_checklists_jun-jul_us-wy.csv")
# combine ebird and environmental data



# Spatiotemporal subsampling ----

# sample one checklist per 3km x 3km x 1 week grid for each year
# sample detection/non-detection independently



# Encounter rate ----

# train a balanced random forest to classify detection/non-detection
# remove observation_count prior to training model



# ├ Calibration ----

# train calibration model


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

# calculate threshold using mcc-f1


# ├ Prediction ----

# load prediction grid
pred_grid <- read_csv("data/environmental-variables_prediction-grid_us-wy.csv")
# raster template for the grid
r <- rast("data/prediction-grid_us-wy.tif")
# get the coordinate reference system of the prediction grid
crs <- st_crs(r)

# add standardized effort covariates to prediction grid
# 6am on 2023-07-01
# 2 km, 1 hour traveling checklist with 1 observer
pred_grid_eff <- pred_grid |>
  mutate(observation_date = ymd("2023-07-01"),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         hours_of_day = 6,
         effort_distance_km = 2,
         effort_hours = 1,
         effort_speed_kmph = 2,
         number_observers = 1)

# estimate calibrated encounter rate and range

# combine predictions with coordinates from prediction grid
predictions <- data.frame(cell_id = pred_grid_eff$cell_id,
                          x = pred_grid_eff$x,
                          y = pred_grid_eff$y,
                          in_range = pred_binary,
                          encounter_rate = pred_er_cal)

# rasterize predictions



# ├ Mapping ----

# project gis data to coordinate reference system of the prediction grid
study_region_proj<- st_transform(study_region, crs = crs)
ne_land_proj <- st_transform(ne_land, crs = crs)
ne_country_lines_proj <- st_transform(ne_country_lines, crs = crs)
ne_state_lines_proj <- st_transform(ne_state_lines, crs = crs)

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
plot(ne_state_lines_proj, col = "#ffffff", lwd = 0.75, add = TRUE)
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
plot(ne_state_lines_proj, col = "#ffffff", lwd = 0.75, add = TRUE)
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

# train a random forest to estimate expected count



# ├ Prediction ----

# add estiamted count and relative abundance to predictions data frame


# rasterize



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
plot(ne_state_lines_proj, col = "#ffffff", lwd = 0.75, add = TRUE)
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

# estimate encounter rate, range boundary, count, and relative abundance
# for test dataset

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




# ├ Count PPMs ----




# Habitat associations ----

# ├ Predictor importance ----

# extract partial dependence from the random forest model objects
# encounter rate

# count

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


# calculate abundance partial dependence for each of the top 6 predictors

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
