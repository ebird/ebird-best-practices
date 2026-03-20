library(dplyr)
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
ebirdst_download_status("babwar",
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
View(ebirdst_runs)

# EXERCISE: Look up one migratory species and one resident species of interest
# to you. Identify the seasonal dates and the review quality ratings.
ebirdst_runs |>
  filter(common_name %in% c("Bay-breasted Warbler", "Bay-headed Tanager")) |>
  View()


# ├ Downloading status data ----

# download data for Bay-breasted Warbler (Setophaga castanea)
ebirdst_download_status("Bay-breasted Warbler")


# ├ Loading raster estimates ----

# load weekly relative abundance
abd_median <- load_raster("babwar", product = "abundance")
# identify weeks corresponding to each layer
names(abd_median)

# load confidence intervals
abd_lower <- load_raster("babwar", product = "abundance", metric = "lower")
abd_upper <- load_raster("babwar", product = "abundance", metric = "upper")


# EXERCISE: Try loading the weekly proportion of population cube at 27 km
# resolution.
prop_pop <- load_raster("babwar", product = "proportion-population",
                        resolution = "27km")

# load seasonal relative abundance
abd_seasonal_mean <- load_raster("babwar",
                                 product = "abundance",
                                 period = "seasonal",
                                 metric = "mean")
# identify seasons corresponding to each layer
names(abd_seasonal_mean)

# load full year maximum relative abundance
abd_yr_max <- load_raster("babwar",
                          product = "abundance",
                          period = "full-year",
                          metric = "max")


# ├ working with raster data ----

# load weekly and seasonal relative abundance at 9km resolution
abd_weekly <- load_raster("babwar", product = "abundance", resolution = "9km")
abd_seasonal <- load_raster("babwar", product = "abundance",
                            period = "seasonal", resolution = "9km")


# subset to a single week or season
abd_median[["2023-12-27"]]
abd_seasonal[["nonbreeding"]]

# subset to all weeks in december and average
in_dec <- month(names(abd_median)) == 12
abd_dec <- mean(abd_weekly[[in_dec]], na.rm = TRUE)

# make a map
plot(abd_dec)

# crop and mask to Colombia
boundary <- read_sf("data/gis-data.gpkg", layer = "countries") |>
  filter(name == "Colombia") |>
  st_transform(crs = st_crs(abd_dec))
abd_dec_co <- abd_dec |>
  crop(boundary) |>
  mask(boundary)
# map
plot(abd_dec_co, axes = FALSE)
plot(st_geometry(boundary), add = TRUE)


# Application: proportions of population within regions ----
# GOAL: estimate the proportion of the non-breeding population of Bay-breasted
# Warbler (Setophaga castanea) within the departments of Colombia.

# seasonal proportion of population
prop_pop_seasonal <- load_raster("babwar",
                                 product = "proportion-population",
                                 period = "seasonal")

# department boundaries
departments <- read_sf("data/gis-data.gpkg", layer = "departments") |>
  st_transform(st_crs(prop_pop_seasonal))

# sum proportion of population within departments
department_prop_pop <- extract(prop_pop_seasonal[["nonbreeding"]],
                               departments,
                               fun = "sum", na.rm = TRUE, weights = TRUE,
                               bind = TRUE) |>
  as.data.frame() |>
  # sort in descending order or non-breeding proportion of population
  arrange(desc(nonbreeding))
head(department_prop_pop)


# Application: migration chronology ----
# GOAL: generate a migration chronology with uncertainty showing the weekly
# proportion of population in Colombia for Bay-breasted Warbler (Setophaga
# castanea).

# load the median weekly relative abundance and lower/upper confidence limits
abd_median <- load_raster("babwar", product = "abundance", metric = "median")
abd_lower <- load_raster("babwar", product = "abundance", metric = "lower")
abd_upper <- load_raster("babwar", product = "abundance", metric = "upper")

# Colombia boundary
region_boundary <- read_sf("data/gis-data.gpkg", layer = "countries") |>
  filter(name == "Colombia") |>
  st_transform(crs = st_crs(abd_median))

# extract values within region and calculate total abundance
abd_median_region <- extract(abd_median, region_boundary,
                             fun = "sum", na.rm = TRUE, ID = FALSE)
abd_lower_region <- extract(abd_lower, region_boundary,
                            fun = "sum", na.rm = TRUE, ID = FALSE)
abd_upper_region <- extract(abd_upper, region_boundary,
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
       y = "Proportion of Population in Colombia",
       title = "Migration chronology for Setophaga castanea") +
  coord_cartesian(ylim = c(0, 1))


# Application: areas of importance ----
# GOAL: identify areas highest importance for a set of 6 migratory warblers
# within the Andinas region of Colombia during the non-breeding season.

# six migratory warblers
# Blackburnian Warbler (Setophaga fusca)
# Cerulean Warbler (Setophaga cerulea)
# Mourning Warbler (Geothlypis philadelphia)
# Bay-breasted Warbler (Setophaga castanea)
# Canada Warbler (Cardellina canadensis)
# Blackpoll Warbler (Setophaga striata)
species_list <- c("bkbwar", "cerwar", "mouwar", "babwar", "canwar", "bkpwar")

# Andinas region boundary
region_boundary <- read_sf("data/gis-data.gpkg", layer = "region") |>
  st_transform(st_crs(abd_median))


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
  abd_masked <- mask(crop(abd, region_boundary), region_boundary)
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
  ebirdst_download_status(species,
                          pattern = "proportion-population_seasonal_mean_3km")

  # load non-breeding season proportion of population
  pp <- load_raster(species,
                    product = "proportion-population",
                    period = "seasonal") |>
    subset("nonbreeding")
  # crop and mask to region
  prop_pop[[species]] <- mask(crop(pp, region_boundary), region_boundary)
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
plot(importance, main = "Important areas for migratory warblers", axes = FALSE)
plot(region_boundary, col = "grey", axes = FALSE, add = TRUE)
plot(importance, axes = FALSE, legend = FALSE, add = TRUE)


# Application: assessing model performance ----
# GOAL: visualize the spatial pattern of predictive performance for Cerulean
# Warbler (Setophaga cerulea) during migration.

# download the spatial predictive performance metrics (ppms)
ebirdst_download_status("cerwar",
                        download_ppms = TRUE,
                        download_abundance = FALSE)

# load the raster for proportion of bernoulli deviance explained
bern_dev <- load_ppm("cerwar", ppm = "occ_bernoulli_dev")

# map the week of april 5
bern_dev_apr <- trim(bern_dev[["04-05"]])
# define a diverging color ramp
ppm_cols <- rev(scico(20, palette = "vik"))
# get the maximum value
max_val <- global(abs(bern_dev_apr), fun = max, na.rm = TRUE) |>
  as.numeric()
# map
plot(bern_dev_apr,
     range = c(-max_val, max_val),
     plg = list(at = seq(-max_val, max_val, length.out = 5)),
     col = ppm_cols,
     axes = FALSE, box = TRUE)
