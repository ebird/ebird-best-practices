library(doParallel)
library(dplyr)
library(exactextractr)
library(foreach)
library(glue)
library(landscapemetrics)
library(readr)
library(sf)
library(stringr)
library(terra)
library(tidyr)
library(units)
library(viridis)
registerDoParallel(14)

f_checklists <- "data/checklists-zf_larbun_jun-jul_us-wy.csv"
time_period <- "jun-jul"
region_code <- "us-wy"

# load and inspect the landcover data
lc_classes <- read_csv("data-raw/mcd12q1_umd_classes.csv")
landcover <- rast(glue("data-raw/landcover_mcd12q1_umd_{region_code}_2014-2022.tif"))

# ebird checklist locations
checklists <- read_csv(f_checklists) %>%
  # landcover data not availble for the full period, so we use the closest year
  mutate(year_lc = as.character(pmin(year, 2022)))

# generate circular neighborhoods for all checkists
pts <- checklists %>%
  # identify unique location/year combinations
  distinct(locality_id, year_lc, latitude, longitude) %>%
  # generate a 3 km neighborhoods
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
centroid <- pts %>%
  st_combine() %>%
  st_centroid() %>%
  st_coordinates() %>%
  round(1) %>%
  data.frame()
buffers <- st_buffer(pts, dist = set_units(1.5, "km"))

# calculate landscape metrics
lsm <- foreach (i = seq_len(nrow(buffers)), .combine = bind_rows) %dopar% {
  buffer_i <- st_transform(buffers[i, ], crs = crs(landcover))
  year <- as.character(buffer_i$year_lc)

  # crop and mask landcover raster so all values outside buffer are missing
  crop(landcover[[year]], buffer_i) %>%
    mask(buffer_i) %>%
    # calcualte landscape metrics
    calculate_lsm(level = "class", metric = c("pland", "ed")) %>%
    # add variables to uniquely identify each point
    mutate(locality_id = buffer_i$locality_id,
           year_lc = buffer_i$year_lc) %>%
    select(locality_id, year_lc, class, metric, value)
}
# transform to wide format
lsm_wide <- lsm %>%
  # fill missing classes with zeros
  complete(nesting(locality_id, year_lc),
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) %>%
  # bring in more descriptive names
  inner_join(select(lc_classes, class, label), by = "class") %>%
  # transform from long to wide format
  pivot_wider(values_from = value,
              names_from = c(class, label, metric),
              names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{label}",
              names_sort = TRUE) %>%
  arrange(locality_id, year_lc)

# elevation raster
elevation <- rast(glue("data-raw/elevation_gmted_1km_{region_code}.tif"))

# mean and standard deviation within each circular neighborhood
elev_buffer <- exact_extract(elevation, buffers, fun = c("mean", "stdev"),
                             progress = FALSE) %>%
  # add variables to uniquely identify each point
  mutate(locality_id = buffers$locality_id, year_lc = buffers$year_lc) %>%
  select(locality_id, year_lc,
         elevation_mean = mean,
         elevation_sd = stdev)

# combine elevation and landcover
env_variables <- inner_join(elev_buffer, lsm_wide,
                            by = c("locality_id", "year_lc"))

# attach and expand to checklists
env_variables <- checklists %>%
  select(checklist_id, locality_id, year_lc) %>%
  inner_join(env_variables, by = c("locality_id", "year_lc")) %>%
  select(-locality_id, -year_lc)

# save to csv, dropping any rows with missing variables
write_csv(drop_na(env_variables),
          glue("data/environmental-variables_checklists_{time_period}_{region_code}.csv"),
          na = "")


# prediction grid ----

# lambert's azimuthal equal area projection
laea_crs <- glue("+proj=laea +lat_0={centroid$Y} +lon_0={centroid$X}") %>%
  st_crs()

# study region
study_region <- read_sf("data/gis-data.gpkg", layer = "ne_states") %>%
  filter(state_code == str_to_upper(region_code)) %>%
  st_transform(crs = laea_crs)

# create a raster template covering the region with 3 km resolution
r <- rast(study_region, res = c(3000, 3000))

# fill the raster with 1s inside the study region
r <- rasterize(study_region, r, values = 1) %>%
  setNames("study_region")

# save for later use
r <- writeRaster(r, glue("data/prediction-grid_{region_code}.tif"),
                 overwrite = TRUE,
                 gdal = "COMPRESS=DEFLATE")

# generate neighborhoods for the prediction grid cell centers
buffers_pg <- as.data.frame(r, cells = TRUE, xy = TRUE) %>%
  select(cell_id = cell, x, y) %>%
  st_as_sf(coords = c("x", "y"), crs = laea_crs, remove = FALSE) %>%
  st_transform(crs = 4326) %>%
  st_buffer(set_units(3, "km"))

# estimate landscape metrics for each cell in the prediction grid
lsm_pg <- foreach (i = seq_len(nrow(buffers_pg)), .combine = bind_rows) %dopar% {
  buffer_i <- st_transform(buffers_pg[i, ], crs = crs(landcover))

  # crop and mask landcover raster so all values outside buffer are missing
  crop(landcover[["2022"]], buffer_i) %>%
    mask(buffer_i) %>%
    # calcualte landscape metrics
    calculate_lsm(level = "class", metric = c("pland", "ed")) %>%
    # add variable to uniquely identify each point
    mutate(cell_id = buffer_i$cell_id) %>%
    select(cell_id, class, metric, value)
}
# transform to wide format
lsm_wide_pg <- lsm_pg %>%
  # fill missing classes with zeros
  complete(cell_id,
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) %>%
  # bring in more descriptive names
  inner_join(select(lc_classes, class, label), by = "class") %>%
  # transform from long to wide format
  pivot_wider(values_from = value,
              names_from = c(class, label, metric),
              names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{label}",
              names_sort = TRUE,
              values_fill = 0) %>%
  arrange(cell_id)

elev_buffer_pg <- exact_extract(elevation, buffers_pg,
                                fun = c("mean", "stdev"),
                                progress = FALSE) %>%
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_pg$cell_id) %>%
  select(cell_id, elevation_mean = mean, elevation_sd = stdev)

# combine landcover and elevation
env_variables_pg <- inner_join(elev_buffer_pg, lsm_wide_pg, by = "cell_id")

# attach the xy coordinates of the cell centers
env_variables_pg <- buffers_pg %>%
  st_drop_geometry() %>%
  select(cell_id, x, y) %>%
  inner_join(env_variables_pg, by = "cell_id")

# save as csv, dropping any rows with missing variables
write_csv(drop_na(env_variables_pg),
          glue("data/environmental-variables_prediction-grid_{region_code}.csv"),
          na = "")
