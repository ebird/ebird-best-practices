library(arrow)
library(auk)
library(dplyr)
library(ggplot2)
library(glue)
library(lubridate)
library(readr)
library(sf)
library(terra)
library(viridis)


# Import eBird data ----

species <- "pnbfin1"
species_name <- ebird_species(species, type = "scientific")

# checklist data
f_sed <- glue("data-raw/ebd_CO_{species}_smp_relAug-2025_sampling.txt")


# observation data
f_ebd <- glue("data-raw/ebd_CO_{species}_smp_relAug-2025.txt")


# EXERCISE: Take some time to explore the variables in these datasets. If
# you're unsure about any of the variables, consult the metadata document that
# came with the data download ("eBird_Basic_Dataset_Metadata_v1.15.pdf").


# ├ Shared checklists ----

# import checklist data without collapsing shared checklists
checklists_shared <- read_sampling(f_sed, unique = FALSE)
# identify shared checklists

# collapse shared checklists



# ├ Taxonomic rollup ----

# import one of the auk example datasets without rolling up taxonomy
obs_ex <- system.file("extdata/ebd-rollup-ex.txt", package = "auk") |>
  read_ebd(rollup = FALSE)
# rollup taxonomy

# identify the taxonomic categories present in each dataset

# yellow-rumped warbler observations prior to rollup

# yellow-rumped warbler observations after rollup



# Filter to region and season ----

# filter to complete checklists from the last 10 years


# subset to andinas region boundary polygon
# convert checklist locations to points geometries

# boundary of andinas region
region_boundary <- read_sf("data/gis-data.gpkg", layer = "region") |>
  st_transform(crs = st_crs(checklists_sf))
# spatially subset the checklists to those in the study region


# remove observations without matching checklists
# this applies the same filters to observations that were applied to checklists



# Zero-fill eBird data ----

# combine checklist and observation data to produce detection/non-detection data

# function to convert observation time to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}
# transform effort variables
# 1. convert counts to integer and "X" to NA
# 2. set distance to 0 for stationary checklists
# 3. convert duration to hours
# 4. create speed variable
# 5. convert time to hours since midnight
# 6. split date into year and day of year


# Apply effort filters ----

# traveling or stationary counts with fewer than 10 observers
# duration <= 8 h and >= 2 minutes, length <= 10 km, speed <= 100km/h


# EXERCISE: Pick one of the four effort variables we filtered on above and
# explore how much variation remains.



# Test-train split ----

# split checklists into 20/80 test/train


# subset to only those columns we need
checklists <- zf_filtered |>
  select(checklist_id, observer_id, type,
         observation_count, species_observed,
         state_code, locality_id, latitude, longitude,
         observation_type,
         observation_date, year, day_of_year, hours_of_day,
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers)
# save dataset for use in next lesson
write_csv(checklists, glue("data/checklists-zf_{species}_co.csv"), na = "")


# Mapping ----

# load gis data
land <- read_sf("data/gis-data.gpkg", "land") |>
  st_geometry()
country_lines <- read_sf("data/gis-data.gpkg", "country_lines") |>
  st_geometry()
region_boundary <- read_sf("data/gis-data.gpkg", "region") |>
  st_geometry()

# prepare ebird data for mapping
checklists_sf <- checklists |>
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  select(species_observed)

# map
par(mar = c(0.25, 0.25, 4, 0.25))
# set up plot area
plot(st_geometry(checklists_sf),
     main = glue("{species_name} eBird observations\n 2015-2024"),
     col = NA, border = NA)
# contextual gis data
plot(land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)
plot(region_boundary, col = "#e6e6e6", border = NA, add = TRUE)
plot(country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
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


# Environmental variables ----

# load and explore the environmental variable dataset
f_habitat <- "data/environmental-variables_checklists_co.parquet"



# Prediction surface ----

# load the prediction surface environmental variables

# load the raster template for the grid


# insert evergreen broadleaf forest % landcover into the raster


# make a map of deciduous broadleaf forest cover
plot(forest_cover,
     axes = FALSE, box = FALSE, col = viridis(10),
     main = "Evergreen Broadleaf Forest (% cover)")

# EXERCISE: Choose one of the other environmental variables and make a map. Does
# the spatial pattern match what you know about the region?
