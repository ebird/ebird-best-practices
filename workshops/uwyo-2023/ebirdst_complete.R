library(dplyr)
library(exactextractr)
library(fs)
library(ggplot2)
library(lubridate)
library(sf)
library(terra)
library(ebirdst)


# copy pre-downloaded data to ebirdst data download location
# ONLY TO BE USED FOR WORKSHOP, DELETE THIS AFTER YOU RUN IT!
ebirdst_dir <- path(ebirdst_data_dir(), "2022")
dir_create(ebirdst_dir)
dirs <- dir_ls("ebirdst-data/")
dir_copy(dirs, path(ebirdst_dir, basename(dirs)))


# data access ----

# get an ebirdst access key https://ebird.org/st/request
# store it using set_ebirdst_access_key("XXXXXXXXX")
# run the following to test the access key
ebirdst_download_status("sagthr",
                        pattern = "abundance_median_27km_2022",
                        force = TRUE)


# status species ----

# ebirdst_runs contains all species for which status data products exist
View(ebirdst_runs)

# EXERCISE: Look up a species of interest to you. Identify the seasonal dates
# and the review quality ratings.
ebirdst_runs |>
  filter(common_name == "Gray-crowned Rosy-Finch") |>
  glimpse()


# downloading status data ----

# download sage thrasher data
ebirdst_download_status("Sage Thrasher")


# loading data ----

# ├ weekly raster estimates ----

# load sage thrasher relative abundance
abd_median <- load_raster("sagthr", product = "abundance")
# identify weeks corresponding to each layer
names(abd_median)

# load confidence intervals
abd_lower <- load_raster("sagthr", product = "abundance", metric = "lower")
abd_upper <- load_raster("sagthr", product = "abundance", metric = "upper")


# EXERCISE: Try loading the weekly proportion of population cube at 27 km
# resolution.
prop_pop <- load_raster("sagthr", product = "proportion-population",
                        resolution = "27km")


# ├ seasonal raster estimates ----

# load seasonal relative abundance
abd_seasonal_mean <- load_raster("sagthr",
                                 product = "abundance",
                                 period = "seasonal",
                                 metric = "mean")
# identify seasons corresponding to each layer
names(abd_seasonal_mean)

# load full year maximum relative abundance
abd_yr_max <- load_raster("sagthr",
                          product = "abundance",
                          period = "full-year",
                          metric = "max")


# working with raster data ----

# load weekly and seasonal relative abundance at 9km resolution
abd_weekly <- load_raster("sagthr", product = "abundance", resolution = "9km")
abd_seasonal <- load_raster("sagthr", product = "abundance",
                            period = "seasonal", resolution = "9km")


# subset to a single week or season
abd_median[["2022-06-07"]]
abd_seasonal[["breeding"]]

# subset to all weeks in june and average
in_june <- month(names(abd_median)) == 6
abd_june <- mean(abd_weekly[[in_june]], na.rm = TRUE)

# make a map
plot(abd_june)

# crop and mask to wyoming
boundary <- read_sf("data/gis-data.gpkg", layer = "ne_states") |>
  filter(state_code == "US-WY") |>
  st_transform(crs = st_crs(abd_june))
abd_june_wy <- abd_june |>
  crop(boundary) |>
  mask(boundary)
# map
plot(abd_june_wy, axes = FALSE)
plot(st_geometry(boundary), add = TRUE)


# application: trajectories ----

# load proportion of population cubes for sage thrasher at 27km
prop_pop_sagthr <- load_raster("sagthr", product = "proportion-population",
                               period = "weekly",
                               resolution = "27km")

# calculate weekly proportion in wyoming
prop_pop_sagthr_wy <- exact_extract(prop_pop_sagthr, boundary, fun = "sum")
prop_pop_sagthr_wy <- data.frame(week = ymd(names(prop_pop_sagthr)),
                                 prop_pop = as.numeric(prop_pop_sagthr_wy))

# plot the trajectories
ggplot(prop_pop_sagthr_wy, aes(x = week, y = prop_pop)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Week",
       y = "% of population",
       title = "Weekly % of population trajectory in Wyoming") +
  theme(legend.position = "bottom")

# compare the trajectories of the following species
species_list <- c("Sage Thrasher", "Sagebrush Sparrow", "Brewer's Sparrow")
trajectories <- NULL
for (species in species_list) {
  ebirdst_download_status(species)
  prop_pop <- load_raster(species, product = "proportion-population",
                          period = "weekly", resolution = "27km")
  traj <- exact_extract(prop_pop, boundary, fun = "sum")
  traj <- data.frame(common_name = species,
                     week = ymd(names(traj)),
                     prop_pop = as.numeric(traj))
  trajectories <- bind_rows(trajectories, traj)
}

# plot the trajectories
ggplot(trajectories, aes(x = week, y = prop_pop, color = common_name)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Week",
       y = "% of population",
       title = "Weekly % of population trajectory in Wyoming",
       color = NULL) +
  theme(legend.position = "bottom")

# EXERCISE: Pick a species and region of interest to you and generate a
# trajectory for either abundance or proportion of population.


# trends species ----

# filter ebirdst_runs to only those species with trends
trends_runs <- ebirdst_runs %>%
  filter(has_trends) %>%
  select(species_code, common_name,
         trends_season, trends_region,
         trends_start_year, trends_end_year,
         trends_start_date, trends_end_date,
         rsquared, beta0)
glimpse(trends_runs)


# download trends data ----

# trends data for sage thrasher
trends_sagthr <- load_trends("Sage Thrasher")


# applications ----

# ├ regional trend ----

# convert trends estimates to spatial format
trends_sagthr_sf <- trends_sagthr |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# wyoming boundary
boundary <- read_sf("data/gis-data.gpkg", layer = "ne_states") |>
  filter(state_code == "US-WY") |>
  st_transform(crs = st_crs(trends_sagthr_sf))

# subset cells to only those falling within wyoming
trends_sagthr_wy <- trends_sagthr_sf[boundary, ]
trends_sagthr_wy <- st_drop_geometry(trends_sagthr_wy)

# abundance weighted trend in wyoming
sum(trends_sagthr_wy$abd * trends_sagthr_wy$abd_ppy) / sum(trends_sagthr_wy$abd)


# ├ regional trend with uncertainty ----

# load fold-level trend data
trends_sagthr_folds <- load_trends("sagthr", fold_estimates = TRUE)

# subset to wyoming
trends_sagthr_folds_sf <- trends_sagthr_folds |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
trends_sagthr_folds_wy <- trends_sagthr_folds_sf[boundary, ] |>
  st_drop_geometry()

# calculate regional trend for each fold
region_trends_folds <- trends_sagthr_folds_wy |>
  group_by(fold) |>
  summarize(abd_ppy = sum(abd * abd_ppy) / sum(abd))

# calculate median and 80% confidence intervals
quantile(region_trends_folds$abd_ppy, probs = c(0.1, 0.5, 0.9))
