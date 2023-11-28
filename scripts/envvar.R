library(dplyr)
library(readr)
library(sf)
library(terra)
library(viridis)


# environmental variables ----

# load and explore the environmenal variable dataset
f_habitat <- "data/environmental-variables_checklists_jun_us-ga.csv"



# prediction grid ----

# load the prediction grid environmental variables
pred_grid <- read_csv("data/environmental-variables_prediction-grid_us-ga.csv")
r <- rast("data/prediction-grid_us-ga.tif") %>%
  # this second rast() call removes all the values from the raster template
  rast()

# insert deciduous broadleaf forest % landcover into the raster
forest_cover <- env_variables_ps %>%
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = laea_crs) %>%
  # rasterize points
  rasterize(r, field = "pland_c04_deciduous_broadleaf")

