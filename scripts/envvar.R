library(dplyr)
library(readr)
library(sf)
library(terra)
library(viridis)


# environmental variables ----

# load and explore the environmental variable dataset
f_habitat <- "data/environmental-variables_checklists_jun_us-ga.csv"
habitat <- read_csv(f_habitat)
glimpse(habitat)


# prediction surface ----

# load the prediction surface environmental variables
pred_grid <- read_csv("data/environmental-variables_prediction-grid_us-ga.csv")
# load the raster template for the grid
r <- rast("data/prediction-grid_us-ga.tif") %>%
  # this second rast() call removes all the values from the raster template
  rast()

# insert deciduous broadleaf forest % landcover into the raster
forest_cover <- pred_grid %>%
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs(r)) %>%
  # rasterize points
  rasterize(r, field = "pland_c04_deciduous_broadleaf")

# make a map of deciduous broadleaf forest cover
plot(forest_cover,
     axes = FALSE, box = FALSE, col = viridis(10),
     main = "Deciduous Broadleaf Forest (% cover)")

# EXERCISE: Choose one of the other environmental variables and make a map. Does
# the spatial pattern match what you know about the region?

