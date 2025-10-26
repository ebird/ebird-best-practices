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
checklists_all <- read_sampling(f_sed)
glimpse(checklists_all)

# observation data
f_ebd <- glue("data-raw/ebd_CO_{species}_smp_relAug-2025.txt")
observations_all <- read_ebd(f_ebd)
glimpse(observations_all)

# EXERCISE: Take some time to explore the variables in these datasets. If
# you're unsure about any of the variables, consult the metadata document that
# came with the data download ("eBird_Basic_Dataset_Metadata_v1.15.pdf").


# ├ Shared checklists ----

# import checklist data without collapsing shared checklists
checklists_shared <- read_sampling(f_sed, unique = FALSE)
# identify shared checklists
checklists_shared |>
  filter(!is.na(group_identifier)) |>
  arrange(group_identifier) |>
  select(sampling_event_identifier, group_identifier)
# collapse shared checklists
checklists_unique <- auk_unique(checklists_shared, checklists_only = TRUE)
nrow(checklists_shared)
nrow(checklists_unique)


# ├ Taxonomic rollup ----

# import one of the auk example datasets without rolling up taxonomy
obs_ex <- system.file("extdata/ebd-rollup-ex.txt", package = "auk") |>
  read_ebd(rollup = FALSE)
# rollup taxonomy
obs_ex_rollup <- auk_rollup(obs_ex)
# identify the taxonomic categories present in each dataset
unique(obs_ex$category)
unique(obs_ex_rollup$category)
# yellow-rumped warbler observations prior to rollup
obs_ex |>
  filter(common_name == "Yellow-rumped Warbler") |>
  select(checklist_id, category, common_name, subspecies_common_name,
         observation_count)
# yellow-rumped warbler observations after rollup
obs_ex_rollup |>
  filter(common_name == "Yellow-rumped Warbler") |>
  select(checklist_id, category, common_name, observation_count)


# Filter to region and season ----

# filter to complete checklists from the last 10 years
checklists <- checklists_all |>
  filter(all_species_reported,
         between(year(observation_date), 2015, 2024))

# subset to andinas region boundary polygon
# convert checklist locations to points geometries
checklists_sf <- st_as_sf(checklists,
                          coords = c("longitude", "latitude"),
                          crs = 4326,
                          remove = FALSE)
# boundary of andinas region
region_boundary <- read_sf("data/gis-data.gpkg", layer = "region") |>
  st_transform(crs = st_crs(checklists_sf))
# spatially subset the checklists to those in the study region
checklists <- checklists_sf[region_boundary, ] |>
  st_drop_geometry()

# remove observations without matching checklists
# this applies the same filters to observations that were applied to checklists
observations <- semi_join(observations_all, checklists, by = "checklist_id")


# Zero-fill eBird data ----

# combine checklist and observation data to produce detection/non-detection data
zf <- auk_zerofill(observations, checklists, collapse = TRUE)
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
zf <- zf |>
  mutate(
    # convert count to integer and X to NA
    # ignore the warning "NAs introduced by coercion"
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for stationary counts
    effort_distance_km = if_else(observation_type == "Stationary",
                                 0, effort_distance_km),
    # convert duration to hours
    effort_hours = duration_minutes / 60,
    # speed km/h
    effort_speed_kmph = effort_distance_km / effort_hours,
    # convert time to decimal hours since midnight
    hours_of_day = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )


# Apply effort filters ----

# traveling or stationary counts with fewer than 10 observers
# duration <= 8 h and >= 2 minutes, length <= 10 km, speed <= 100km/h
zf_filtered <- zf |>
  filter(observation_type %in% c("Stationary", "Traveling"),
         !is.na(effort_hours), effort_hours >= 0.032, effort_hours <= 8,
         !is.na(effort_distance_km), effort_distance_km <= 10,
         effort_speed_kmph <= 100,
         number_observers <= 10)

# EXERCISE: Pick one of the four effort variables we filtered on above and
# explore how much variation remains.
ggplot(zf) +
  aes(x = effort_hours) +
  geom_histogram(binwidth = 0.5,
                 aes(y = after_stat(count / sum(count)))) +
  scale_y_continuous(limits = c(0, NA), labels = scales::label_percent()) +
  labs(x = "Duration [hours]",
       y = "% of eBird checklists",
       title = "Distribution of eBird checklist duration",
       subtitle = "Before effort filtering")
ggplot(zf_filtered) +
  aes(x = effort_hours) +
  geom_histogram(binwidth = 0.5,
                 aes(y = after_stat(count / sum(count)))) +
  scale_y_continuous(limits = c(0, NA), labels = scales::label_percent()) +
  labs(x = "Duration [hours]",
       y = "% of eBird checklists",
       title = "Distribution of eBird checklist duration",
       subtitle = "After effort filtering")


# Test-train split ----

# split checklists into 20/80 test/train
zf_filtered$type <- if_else(runif(nrow(zf_filtered)) <= 0.8, "train", "test")
table(zf_filtered$type) / nrow(zf_filtered)

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
habitat <- read_parquet(f_habitat)
glimpse(habitat)


# Prediction surface ----

# load the prediction surface environmental variables
pred_grid <- read_parquet("data/prediction-grid_co.parquet")
# load the raster template for the grid
r <- rast("data/prediction-grid_co.tif")

# insert evergreen broadleaf forest % landcover into the raster
forest_cover <- pred_grid |>
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs(r)) |>
  # rasterize points
  rasterize(r, field = "landcover_evergreen_broadleaf_pland")

# make a map of deciduous broadleaf forest cover
plot(forest_cover,
     axes = FALSE, box = FALSE, col = viridis(10),
     main = "Evergreen Broadleaf Forest (% cover)")

# EXERCISE: Choose one of the other environmental variables and make a map. Does
# the spatial pattern match what you know about the region?
