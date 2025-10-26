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


# EXERCISE: Look up a species of interest to you. Identify the seasonal dates
# and the review quality ratings.



# downloading status data ----

# download sage thrasher data



# loading data ----

# ├ weekly raster estimates ----

# load sage thrasher relative abundance

# identify weeks corresponding to each layer


# load confidence intervals



# EXERCISE: Try loading the weekly proportion of population cube at 27 km
# resolution.



# ├ seasonal raster estimates ----

# load seasonal relative abundance

# identify seasons corresponding to each layer


# load full year maximum relative abundance



# working with raster data ----

# load weekly and seasonal relative abundance at 9km resolution



# subset to a single week or season


# subset to all weeks in june and average


# make a map


# crop and mask to wyoming
boundary <- read_sf("data/gis-data.gpkg", layer = "ne_states") |>
  filter(state_code == "US-WY") |>
  st_transform(crs = st_crs(abd_june))

# map



# application: trajectories ----

# load proportion of population cubes for sage thrasher at 27km


# calculate weekly proportion in wyoming


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



# download trends data ----

# trends data for sage thrasher



# applications ----

# ├ regional trend ----

# convert trends estimates to spatial format


# wyoming boundary
boundary <- read_sf("data/gis-data.gpkg", layer = "ne_states") |>
  filter(state_code == "US-WY") |>
  st_transform(crs = st_crs(trends_sagthr_sf))

# subset cells to only those falling within wyoming


# abundance weighted trend in wyoming



# ├ regional trend with uncertainty ----

# load fold-level trend data


# subset to wyoming


# calculate regional trend for each fold


# calculate median and 80% confidence intervals

