library(dplyr)
library(exactextractr)
library(fs)
library(ggplot2)
library(lubridate)
library(scico)
library(sf)
library(terra)
library(ebirdst)
extract <- terra::extract


# Data access ----

# get an ebirdst access key https://ebird.org/st/request
# store it using set_ebirdst_access_key("XXXXXXXXX")
# run the following to test the access key
ebirdst_download_status("taibar1",
                        pattern = "abundance_median_27km_2023",
                        force = TRUE)

# copy pre-downloaded data to ebirdst data download location
# ONLY TO BE USED FOR WORKSHOP, DELETE THIS AFTER YOU RUN IT!
ebirdst_dir <- path(ebirdst_data_dir(), "2023")
dir_create(ebirdst_dir)
dirs <- dir_ls("ebirdst-data/2023/")
dir_copy(dirs, path(ebirdst_dir, basename(dirs)), overwrite = TRUE)


# Status results ----

# ebirdst_runs contains all species for which status data products exist


# EXERCISE: Look up one migratory species and one resident species of interest
# to you. Identify the seasonal dates and the review quality ratings.



# ├ Downloading status data ----

# download data for Brown-headed Thrush



# ├ Loading raster estimates ----

# load weekly relative abundance

# identify weeks corresponding to each layer


# load confidence intervals



# EXERCISE: Try loading the weekly proportion of population cube at 27 km
# resolution.


# load seasonal relative abundance

# identify seasons corresponding to each layer


# load full year maximum relative abundance



# ├ working with raster data ----

# load weekly and seasonal relative abundance at 9km resolution


# subset to a single week or season


# subset to all weeks in december and average


# make a map


# crop to the island of Taiwan
bounding_box <- read_sf("data/gis-data.gpkg", layer = "bounding_box") |>
  st_transform(crs = st_crs(abd_dec))

# map



# Application: proportions of population within regions ----
# GOAL: estimate the proportion of the population of Taiwan Barbet within each
# of the national parks of Taiwan.

# download just the 3km seasonal proportion of population


# seasonal proportion of population


# national park boundaries from world database of protected areas
np_boundaries <- read_sf("data/gis-data.gpkg", layer = "national_parks") |>
  st_transform(st_crs(prop_pop_seasonal))

# sum proportion of population within national parks



# Application: migration chronology ----
# GOAL: generate a migration chronology with uncertainty showing the weekly
# proportion of population in Taiwan for Brown-headed Thrush.

# load the median weekly relative abundance and lower/upper confidence limits


# boundary of Taiwan
taiwan_boundary <- read_sf("data/gis-data.gpkg", layer = "countries") |>
  filter(name == "Taiwan") |>
  st_transform(st_crs(abd_median))

# extract values within region and calculate total abundance


# transform to data frame format with rows corresponding to weeks
chronology <- data.frame(week = as.Date(names(abd_median)),
                         median = as.numeric(abd_median_region),
                         lower = as.numeric(abd_lower_region),
                         upper = as.numeric(abd_upper_region))

# convert to percent of population by dividing by total abundance globally


# plot chronology
ggplot(chronology) +
  aes(x = week, y = median) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(x = "Week",
       y = "Proportion of Population in Taiwan",
       title = "Migration chronology for Brown-headed Thrush") +
  coord_cartesian(ylim = c(0, 1))

# EXERCISE: notice that there is a portion of the migration chronology missing
# in September-October. What do you think is happening here? Look at the weekly
# relative abundance maps for this period online for a hint.
# https://science.ebird.org/en/status-and-trends/species/brhthr1/abundance-map-weekly?week=38


# Application: areas of importance ----
# GOAL: identify areas highest importance for a set of 6 migratory shorebirds
# within Taiwan during the non-breeding season.

# six shorebirds
# Long-toed Stint (lotsti)
# Red Knot (redkno)
# Red-necked Stint (rensti)
# Little Ringed Plover (lirplo)
# Ruddy Turnstone (rudtur)
# Black-bellied Plover (bkbplo)
species_list <- c("lotsti", "redkno", "rensti", "lirplo", "rudtur", "bkbplo")

# boundary of Taiwan
taiwan_boundary <- read_sf("data/gis-data.gpkg", layer = "countries") |>
  filter(name == "Taiwan") |>
  st_transform(st_crs(abd_median))
# buffer by a single grid cell to ensure we capture estimates just offshore
taiwan_boundary <- st_buffer(taiwan_boundary, dist = 3000)


# ├ Richness ----

# calculate richness at 3km resolution for this group

# sum across species to calculate richness


# map richness



# ├ Importance ----

# calculate mean proportion of population at 3km resolution for this group
prop_pop <- list()
for (species in species_list) {
  # download seasonal proportion of population at 3km
  ebirdst_download_status(species,
                          pattern = "proportion-population_seasonal_mean_3km")

  # load non-breeding season proportion of population
  pp <- load_raster(species,
                    product = "proportion-population",
                    period = "seasonal") |>
    subset("nonbreeding")
  # crop and mask to region
  prop_pop[[species]] <- mask(crop(pp, taiwan_boundary), taiwan_boundary)
}
# take mean across species


# map areas of importance


# drop zeros

# convert to percentile ranks
convert_to_pct_rank <- function(x) {
  rank(x, na.last = "keep", ties.method = "max") / length(na.omit(x))
}


# make a simple map
plot(importance, main = "Important areas for migratory shorebirds", axes = FALSE)
plot(taiwan_boundary, col = "grey", axes = FALSE, add = TRUE)
plot(importance, axes = FALSE, legend = FALSE, add = TRUE)


# Application: assessing model performance ----
# GOAL: visualize the spatial pattern of predictive performance for Brown-headed
# Thrush during the non-breeding season.

# download the spatial predictive performance metrics (ppms)


# load the raster for proportion of bernoulli deviance explained


# seasonal dates

# subset to non-breeding weeks and average


# mask to within range

# convert to range

# mask

# trim to remove areas with no data


# define a diverging color ramp
ppm_cols <- rev(scico(20, palette = "vik"))
# get the maximum value


# map
land <- read_sf("data/gis-data.gpkg", layer = "land") |>
  st_transform(crs = st_crs(bern_dev_masked))
plot(bern_dev_masked,
     range = c(-max_val, max_val),
     plg = list(at = seq(-max_val, max_val, length.out = 5)),
     col = ppm_cols,
     axes = FALSE, box = TRUE)
plot(land, col = "grey", border = FALSE, axes = FALSE, add = TRUE)
plot(bern_dev_masked,
     range = c(-max_val, max_val),
     plg = list(at = seq(-max_val, max_val, length.out = 5)),
     col = ppm_cols,
     axes = FALSE, box = TRUE, legend = FALSE, add = TRUE)
