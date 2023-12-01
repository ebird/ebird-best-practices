library(auk)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(readr)
library(sf)


# Import eBird data ----

# checklist data
f_sed <- "data-raw/ebd_US-GA_woothr_smp_relOct-2023_sampling.txt"
checklists <- read_sampling(f_sed)
glimpse(checklists)

# observation data
f_ebd <- "data-raw/ebd_US-GA_woothr_smp_relOct-2023.txt"
observations <- read_ebd(f_ebd)
glimpse(observations)

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

# complete checklists in june from the last 10 years
# filter the checklist data
checklists <- checklists |>
  filter(all_species_reported,
         between(year(observation_date), 2014, 2023),
         month(observation_date) == 6)
# filter the observation data
observations <- observations |>
  filter(all_species_reported,
         between(year(observation_date), 2014, 2023),
         month(observation_date) == 6)


# Zero-fill eBird data ----

# combine checklist and observation data to produce detection/non-detection data
zf <- auk_zerofill(observations, checklists, collapse = TRUE)
# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}
# transform effort variables
# 1. converet counts to integer and "X" to NA
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
    effort_distance_km = if_else(protocol_type == "Stationary",
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
# duration <= 6 h, length <= 10 km, speed <= 100km/h
zf_filtered <- zf |>
  filter(protocol_type %in% c("Stationary", "Traveling"),
         effort_hours <= 6,
         effort_distance_km <= 10,
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
       title = "Distribution of eBird checklist duration")
ggplot(zf_filtered) +
  aes(x = effort_hours) +
  geom_histogram(binwidth = 0.5,
                 aes(y = after_stat(count / sum(count)))) +
  scale_y_continuous(limits = c(0, NA), labels = scales::label_percent()) +
  labs(x = "Duration [hours]",
       y = "% of eBird checklists",
       title = "Distribution of eBird checklist duration")


# Test-train split ----

# split checklists into 20/80 test/train
zf_filtered$type <- if_else(runif(nrow(zf_filtered)) <= 0.8, "train", "test")
table(zf_filtered$type) / nrow(zf_filtered)

# subset to only those columns we need
checklists <- zf_filtered |>
  select(checklist_id, observer_id, type,
         observation_count, species_observed,
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         hours_of_day,
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers)
# save dataset for use in next lesson
write_csv(checklists, "data/checklists-zf_woothr_jun_us-ga.csv", na = "")


# Mapping ----

# load gis data
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") |>
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") |>
  st_geometry()
ne_state_lines <- read_sf("data/gis-data.gpkg", "ne_state_lines") |>
  st_geometry()
study_region <- read_sf("data/gis-data.gpkg", "ne_states") |>
  filter(state_code == "US-GA") |>
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
     main = "Wood Thrush eBird Observations\nJune 2014-2023",
     col = NA, border = NA)
# contextual gis data
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)
plot(study_region, col = "#e6e6e6", border = NA, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
# ebird observations
# not observed
plot(filter(checklists_sf, !species_observed),
     pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
     add = TRUE)
# observed
plot(filter(checklists_sf, species_observed),
     pch = 19, cex = 0.3, col = alpha("#4daf4a", 1),
     add = TRUE)
# legend
legend("bottomright", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("eBird checklist", "Wood Thrush sighting"),
       pch = 19)
box()
