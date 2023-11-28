library(ebirdst)
library(fields)
library(gridExtra)
library(mccf1)
library(ranger)
library(scam)
library(sf)
library(terra)
library(verification)
library(tidyverse)

# set random number seed to insure fully repeatable results
set.seed(1)


# Load data ----

# environmental variables: landcover and elevation

# zero-filled ebird data

# combine ebird and environmental data

# prediction grid

# raster template for the grid

# get the coordinate reference system of the prediction grid


# load gis data for making maps
study_region <- read_sf("data/gis-data.gpkg", "ne_states") %>%
  filter(state_code == "US-GA") %>%
  st_transform(crs = crs) %>%
  st_geometry()
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>%
  st_transform(crs = crs) %>%
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") %>%
  st_transform(crs = crs) %>%
  st_geometry()
ne_state_lines <- read_sf("data/gis-data.gpkg", "ne_state_lines") %>%
  st_transform(crs = crs) %>%
  st_geometry()


# Spatiotemporal subsampling ----

# sample one checklist per 3km x 3km x 1 week grid for each year
# sample detection/non-detection independently


# EXERCISE: Compare the full set of eBird checklists to the set of checklists
# remaining after subsampling. What was the change in sampled size and how did
# the subsampling impact the prevalence of detections compared to
# non-detections?


# Random forests ----

# filter to training data, select only the columns to be used in the model

# calculate detection frequency

# train a balanced random forest to classify detection/non-detection


# ├ Calibration ----

# predicted encounter rate based on out of bag samples

# observed detection, converted back from factor

# construct a data frame to train the scam model

# fit calibration model

# prepare to make a calibration plot
# group the predicted encounter rate into bins of width 0.02
# then calculate the mean observed encounter rates in each bin
er_breaks <- seq(0, 1, by = 0.02)

# make predictions from the calibration model

# compared binned mean encounter rates to calibration model
ggplot(cal_pred) +
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

# identify best threshold


# ├ Assessment ----

# get the test set held out from training

# predict on test data using random forest model

# convert predictions to binary (presence/absence) using the threshold

# calibrate

# combine observations and estimates

# mean squared error (mse)
mse <- mean((obs_pred_test$obs - obs_pred_test$pred_calibrated)^2, na.rm = TRUE)

# spearman correlation, based on in range observations only
spearman <- cor(obs_pred_test$pred_calibrated[obs_pred_test$pred_binary > 0],
                obs_pred_test$obs[obs_pred_test$pred_binary > 0],
                method = "spearman")

# precision-recall auc
em <- precrec::evalmod(scores = obs_pred_test$pred_binary,
                       labels = obs_pred_test$obs)
pr_auc <- precrec::auc(em) %>%
  filter(curvetypes == "PRC") %>%
  pull(aucs)

# calculate metrics for binary prediction: kappa, sensitivity, specificity
pa_metrics <- obs_pred_test %>%
  select(id, obs, pred_binary) %>%
  PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)

# mcc and f1
mcc_f1 <- calculate_mcc_f1(obs_pred_test$obs, obs_pred_test$pred_binary)

# combine metrics together
data.frame(
  mse = mse,
  spearman = spearman,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity,
  kappa = pa_metrics$Kappa,
  pr_auc = pr_auc,
  mcc = mcc_f1$mcc,
  f1 = mcc_f1$f1
)


# Habitat associations ----

# ├ Predictor importance ----

# extract partial dependence from the random forest model object

# plot predictor importance for top 20 predictors
ggplot(head(pi, 20)) +
  aes(x = fct_reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL,
       y = "Predictor Importance (Gini Index)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", linewidth = 0.5))


# ├ Partial dependence ----

# function to calculate partial dependence for a given predictor
calculate_pd <- function(predictor, er_model, calibration_model,
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

  # predict encounter rate
  p <- predict(er_model, data = grid)

  # summarize
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "x")
  pd$encounter_rate <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, predictor, x)
  pd <- dplyr::summarise(pd,
                         encounter_rate = mean(encounter_rate, na.rm = TRUE),
                         .groups = "drop")

  # calibrate
  pd$encounter_rate <- predict(calibration_model,
                               newdata = data.frame(pred = pd$encounter_rate),
                               type = "response")
  pd$encounter_rate <- as.numeric(pd$encounter_rate)
  # constrain to 0-1
  pd$encounter_rate[pd$encounter_rate < 0] <- 0
  pd$encounter_rate[pd$encounter_rate > 1] <- 1

  return(pd)
}

# calculate partial dependence for each of the top 6 predictors

# plot partial dependence
ggplot(pd) +
  aes(x = x, y = encounter_rate) +
  geom_line() +
  geom_point() +
  facet_wrap(~ as_factor(predictor), ncol = 2, scales = "free") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

# EXERCISE: Examine the predictor importance data to identify an the next most
# important percent landcover variables after pland_c04_deciduous_broadleaf.
# Make a partial dependence plot for this variable.


# Prediction ----

# ├ Standardized effort variables ----

# identify peak time of day for observing this species
# partial dependence for hours_of_day

# partial dependence plot
ggplot(pd_time) +
  aes(x = x, y = encounter_rate) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(x = "Hours since midnight",
       y = "Encounter rate",
       title = "Observation start time partial dependence")
# identify time maximizing encounter rate


# add standardized effort variables to prediction grid


# ├ Model estimates ----

# predict encounter rate

# define range-boundary

# apply calibration

# constrain to 0-1

# combine predictions with coordinates from prediction grid

# rasterize predictions


# Mapping ----

# ├ Range ----

par(mar = c(0.25, 0.25, 1.25, 0.25))
# set up plot area
plot(study_region,
     main = "Wood Thrush Range (June 2023)",
     col = NA, border = NA)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# convert binary prediction to categorical
r_range <- as.factor(r_pred[["in_range"]])
plot(r_range, col = c("#e6e6e6", "forestgreen"),
     maxpixels = ncell(r_range),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()


# ├ Encounter rate ----

par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region, col = NA, border = NA)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks
brks <- global(r_pred[["encounter_rate"]], fun = quantile,
               probs = seq(0, 1, 0.1), na.rm = TRUE) %>%
  as.numeric() %>%
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
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Wood Thrush Encounter Rate (June 2023)"
image.plot(zlim = c(0, 1), legend.only = TRUE,
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.03, 0.06),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))


# ├ Encounter rate within range ----

# within range encounter rate


# map
par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region, col = NA, border = NA)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks, excluding zeros
brks <- ifel(er_range > 0, er_range, NA) %>%
  global(fun = quantile,
         probs = seq(0, 1, 0.1), na.rm = TRUE) %>%
  as.numeric() %>%
  unique()
# label the bottom, middle, and top value
lbls <- round(c(min(brks), median(brks), max(brks)), 2)
# ebird status and trends color palette
pal <- ebirdst_palettes(length(brks) - 1)
plot(er_range,
     col = c("#e6e6e6", pal), breaks = c(0, brks),
     maxpixels = ncell(r_pred),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Wood Thrush Encounter Rate within Range (June 2023)"
image.plot(zlim = c(0, 1), legend.only = TRUE,
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.03, 0.06),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))
