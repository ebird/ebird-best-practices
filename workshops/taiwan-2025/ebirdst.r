library(dplyr)
library(fs)
library(ggplot2)
library(lubridate)
library(scico)
library(sf)
library(terra)
library(ebirdst)
extract <- terra::extract


# Status results ----

# ebirdst_runs contains all species for which status data products exist
View(ebirdst_runs)

# EXERCISE: Look up one migratory species and one resident species of interest
# to you. Identify the seasonal dates and the review quality ratings.
ebirdst_runs |>
  filter(common_name %in% c("Taiwan Barbet", "Brown-headed Thrush")) |>
  View()


# ├ Downloading status data ----

# download data for Brown-headed Thrush
ebirdst_download_status("Brown-headed Thrush")


# ├ Loading raster estimates ----

# load weekly relative abundance
abd_median <- load_raster("brhthr1", product = "abundance")
# identify weeks corresponding to each layer
names(abd_median)

# load confidence intervals
abd_lower <- load_raster("brhthr1", product = "abundance", metric = "lower")
abd_upper <- load_raster("brhthr1", product = "abundance", metric = "upper")


# EXERCISE: Try loading the weekly proportion of population cube at 27 km
# resolution.
prop_pop <- load_raster("brhthr1", product = "proportion-population",
                        resolution = "27km")

# load seasonal relative abundance
abd_seasonal_mean <- load_raster("brhthr1",
                                 product = "abundance",
                                 period = "seasonal",
                                 metric = "mean")
# identify seasons corresponding to each layer
names(abd_seasonal_mean)

# load full year maximum relative abundance
abd_yr_max <- load_raster("brhthr1",
                          product = "abundance",
                          period = "full-year",
                          metric = "max")


# ├ working with raster data ----

# load weekly and seasonal relative abundance at 9km resolution
abd_weekly <- load_raster("brhthr1", product = "abundance", resolution = "9km")
abd_seasonal <- load_raster("brhthr1", product = "abundance",
                            period = "seasonal", resolution = "9km")

# subset to a single week or season
abd_median[["2023-12-27"]]
abd_seasonal[["nonbreeding"]]

# subset to all weeks in december and average
in_dec <- month(names(abd_median)) == 12
abd_dec <- mean(abd_weekly[[in_dec]], na.rm = TRUE)

# make a map
plot(abd_dec)

# crop to the island of Taiwan
bounding_box <- read_sf("data/gis-data.gpkg", layer = "bounding_box") |>
  st_transform(crs = st_crs(abd_dec))
abd_dec_tw <- crop(abd_dec, bounding_box)
# map
plot(abd_dec_tw, axes = FALSE)


# Application: proportions of population within regions ----
# GOAL: estimate the proportion of the population of Taiwan Barbet within each
# of the national parks of Taiwan.

# download just the 3km seasonal proportion of population
ebirdst_download_status("taibar1",
                        pattern = "proportion-population_seasonal_mean_3km")

# seasonal proportion of population
prop_pop_seasonal <- load_raster("taibar1",
                                 product = "proportion-population",
                                 period = "seasonal")

# national park boundaries from world database of protected areas
np_boundaries <- read_sf("data/gis-data.gpkg", layer = "national_parks") |>
  st_transform(st_crs(prop_pop_seasonal))

# sum proportion of population within national parks
prop_pop_np <- extract(prop_pop_seasonal,
                       np_boundaries,
                       fun = "sum", na.rm = TRUE, weights = TRUE,
                       bind = TRUE) |>
  as.data.frame() |>
  arrange(desc(resident))
print(prop_pop_np)


# Application: migration chronology ----
# GOAL: generate a migration chronology with uncertainty showing the weekly
# proportion of population in Taiwan for Brown-headed Thrush.

# load the median weekly relative abundance and lower/upper confidence limits
abd_median <- load_raster("brhthr1", product = "abundance", metric = "median")
abd_lower <- load_raster("brhthr1", product = "abundance", metric = "lower")
abd_upper <- load_raster("brhthr1", product = "abundance", metric = "upper")

# boundary of Taiwan
taiwan_boundary <- read_sf("data/gis-data.gpkg", layer = "countries") |>
  filter(name == "Taiwan") |>
  st_transform(st_crs(abd_median))

# extract values within region and calculate total abundance
abd_median_region <- extract(abd_median, taiwan_boundary,
                             fun = "sum", na.rm = TRUE, ID = FALSE)
abd_lower_region <- extract(abd_lower, taiwan_boundary,
                            fun = "sum", na.rm = TRUE, ID = FALSE)
abd_upper_region <- extract(abd_upper, taiwan_boundary,
                            fun = "sum", na.rm = TRUE, ID = FALSE)

# transform to data frame format with rows corresponding to weeks
chronology <- data.frame(week = as.Date(names(abd_median)),
                         median = as.numeric(abd_median_region),
                         lower = as.numeric(abd_lower_region),
                         upper = as.numeric(abd_upper_region))

# convert to percent of population by dividing by total abundance globally
total_abd <- global(abd_median, sum, na.rm = TRUE)
total_abd <- total_abd$sum
chronology$median <- chronology$median / total_abd
chronology$lower <- chronology$lower / total_abd
chronology$upper <- chronology$upper / total_abd

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
range_rasters <- list()
for (species in species_list) {
  # download seasonal abundance at 3km
  ebirdst_download_status(species, pattern = "abundance_seasonal_mean_3km")

  # load non-breeding season relative abundance
  abd <- load_raster(species, period = "seasonal") |>
    subset("nonbreeding")
  # crop and mask to region
  abd_masked <- mask(crop(abd, taiwan_boundary), taiwan_boundary)
  # convert to binary, presence-absence
  range_rasters[[species]] <- abd_masked > 0
}
# sum across species to calculate richness
richness <- sum(rast(range_rasters), na.rm = TRUE)

# map richness
plot(richness, axes = FALSE)


# ├ Importance ----

# calculate mean proportion of population at 3km resolution for this group
prop_pop <- list()
for (species in species_list) {
  # download seasonal proportion of population at 3km
  # ebirdst_download_status(species,
  #                         pattern = "proportion-population_seasonal_mean_3km")

  # load non-breeding season proportion of population
  pp <- load_raster(species,
                    product = "proportion-population",
                    period = "seasonal") |>
    subset("nonbreeding")
  # crop and mask to region
  prop_pop[[species]] <- mask(crop(pp, taiwan_boundary), taiwan_boundary)
}
# take mean across species
importance <- mean(rast(prop_pop), na.rm = TRUE)

# map areas of importance
plot(importance, axes = FALSE)

# drop zeros
importance <- ifel(importance == 0, NA, importance)
# convert to percentile ranks
convert_to_pct_rank <- function(x) {
  rank(x, na.last = "keep", ties.method = "max") / length(na.omit(x))
}
importance <- app(importance, convert_to_pct_rank)

# make a simple map
plot(importance, main = "Important areas for migratory shorebirds", axes = FALSE)
plot(taiwan_boundary, col = "grey", axes = FALSE, add = TRUE)
plot(importance, axes = FALSE, legend = FALSE, add = TRUE)


# Application: assessing model performance ----
# GOAL: visualize the spatial pattern of predictive performance for Brown-headed
# Thrush during the non-breeding season.

# download the spatial predictive performance metrics (ppms)
ebirdst_download_status("brhthr1", download_ppms = TRUE)

# load the raster for proportion of bernoulli deviance explained
bern_dev <- load_ppm("brhthr1", ppm = "occ_bernoulli_dev")

# seasonal dates
seasonal_dates <- ebirdst_runs |>
  filter(species_code == "brhthr1") |>
  select(start = nonbreeding_start, end = nonbreeding_end)
weeks <- paste0("2023-", names(bern_dev))

# subset to non-breeding weeks and average
in_nonbreeding <- weeks >= seasonal_dates$start | weeks <= seasonal_dates$end
bern_dev_nonbreeding <- bern_dev[[in_nonbreeding]] |>
  mean(na.rm = TRUE)

# mask to within range
abd_nonbreeding <- load_raster("brhthr1",
                               product = "abundance",
                               period = "seasonal",
                               resolution = "27km") |>
  subset("nonbreeding")
# convert relative abundance to range
range_nonbreeding <- ifel(abd_nonbreeding > 0, 1, NA)
# mask
bern_dev_masked <- mask(bern_dev_nonbreeding, range_nonbreeding)
# trim to remove areas with no data
bern_dev_masked <- trim(bern_dev_masked)

# define a diverging color ramp
ppm_cols <- rev(scico(20, palette = "vik"))
# get the maximum value
max_val <- global(abs(bern_dev_masked), fun = max, na.rm = TRUE) |>
  as.numeric()

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
