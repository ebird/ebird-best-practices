library(ebirdst)
library(fields)
library(gridExtra)
library(mccf1)
library(ranger)
library(scam)
library(sf)
library(terra)
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


# filter to training data, select only the columns to be used in the model
checklists_train <- checklists_ss %>%
  filter(type == "train") %>%
  # select only the columns to be used in the model
  select(species_observed, observation_count,
         year, day_of_year, hours_of_day,
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers,
         starts_with("pland_"),
         starts_with("ed_"),
         starts_with("elevation_"))


# Hurdle model ----

# ├ Encounter rate ----

# calculate detection frequency for the balance random forest

# train a random forest model for encounter rate
# remove observation_count prior to training model


# select the mcc-f1 optimizing occurrence threshold


# train calibration model


# ├ Count ----

# subset to only observed or predicted detections


# add predicted encounter rate as an additional covariate


# train a random forest to estimate expected count


# ├ Assessment ----

# get the test set held out from training

# estimate encounter rate

# convert predictions to binary (presence/absence) using the threshold

# calibrate

# add predicted encounter rate required for count estimates

# combine predictions

# estimate count

# combine observations and estimates


# get the test set held out from training
# subset to only those checklists where detection is predicted
detections_test <- filter(obs_pred_test, pred_binary > 0)

# count metrics
count_spearman <- cor(detections_test$pred_count,
                      detections_test$obs_count,
                      method = "spearman")
log_count_pearson <- cor(log(detections_test$pred_count + 1),
                         log(detections_test$obs_count + 1),
                         method = "pearson")

# abundance metrics
abundance_spearman <- cor(detections_test$pred_abundance,
                          detections_test$obs_count,
                          method = "spearman")
log_abundance_pearson <- cor(log(detections_test$pred_abundance + 1),
                             log(detections_test$obs_count + 1),
                             method = "pearson")

# combine metrics together
data.frame(
  count_spearman = count_spearman,
  log_count_pearson = log_count_pearson,
  abundance_spearman = abundance_spearman,
  log_abundance_pearson = log_abundance_pearson
)


# Prediction ----

# add standardized effort variables to prediction grid

# encounter rate estimate

# estimate encounter rate

# define range-boundary

# apply calibration

# constrain to 0-1

# add predicted encounter rate required for count estimates

# estimate count


# combine predictions with coordinates from prediction grid
predictions <- data.frame(cell_id = pred_grid_eff$cell_id,
                          x = pred_grid_eff$x,
                          y = pred_grid_eff$y,
                          in_range = pred_binary,
                          encounter_rate = pred_er_cal,
                          count = pred_count)
# add relative abundance estimate

# rasterize
layers <- c("in_range", "encounter_rate", "count", "abundance")


# in range abundance


par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region, col = NA, border = NA)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks, excluding zeros
brks <- ifel(r_plot > 0, r_plot, NA) %>%
  global(fun = quantile,
         probs = seq(0, 1, 0.1), na.rm = TRUE) %>%
  as.numeric() %>%
  unique()
# label the bottom, middle, and top value
lbls <- round(c(min(brks), median(brks), max(brks)), 2)
# ebird status and trends color palette
pal <- ebirdst_palettes(length(brks) - 1)
plot(r_plot,
     col = c("#e6e6e6", pal), breaks = c(0, brks),
     maxpixels = ncell(r_plot),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Wood Thrush Relative Abundance (June 2023)"
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
