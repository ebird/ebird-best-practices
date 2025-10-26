library(arrow)
library(dplyr)
library(fs)
library(jsonlite)
library(rnaturalearth)
library(sf)
library(terra)
library(tidyr)

base_dir <- "/ocean/projects/deb200005p/sligocki/"
erd_dir <- path(base_dir, "erd_prep", "erd2024")
inputs_dir <- path(erd_dir, "covar_assignment", "inputs")
outputs_dir <- "/ocean/projects/deb200005p/mstrimas/tw-erd-srd-1km/"
dir_create(outputs_dir)

cov_def <- read_json("~/stem_hwf/erd/assign_covars/covariate-definitions.json",
                     simplifyVector = TRUE)

crs <- crs("epsg:8857")
start_year <- 2005
end_year <- 2024

# region boundary
region <- ne_countries(scale = 10, country = "Taiwan") |>
  st_set_precision(1e6) |>
  st_make_valid() |>
  st_union() |>
  st_buffer(10000, max_cells = 1e5)
region_proj <- st_transform(region, crs) |> vect()
bb <- st_bbox(region) |>
  as.list()


# 1 km raster template ----

f_tif <- path(outputs_dir, "prediction-grid_1km_tw.tif")
srd3km <- path(base_dir, "stem_hwf", "ERD2024_inputs", "srd_templates",
              "srd_3km_mask_land.tif") |>
  rast()
srd1km <- srd3km |>
  disagg(fact = 3, method = "near") |>
  crop(region_proj) |>
  mask(region_proj,
       filename = f_tif, overwrite = TRUE,
       datatype = "INT1U", NAflag = 0L,
       gdal = c("COMPRESS=DEFLATE",
                "TILED=YES",
                "COPY_SRC_OVERVIEWS=YES"))


# srd cell coordinates ----

f_pq <- path(outputs_dir, "srd_1km_tw_cell-coordinates.parquet")
srd_pts <- as.data.frame(srd1km,
                         cells = TRUE, xy = TRUE,
                         na.rm = TRUE) |>
  select(srd_id = cell, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = crs, remove = FALSE) |>
  st_transform(crs = 4326)
coords <- st_coordinates(srd_pts) |>
  as.data.frame() |>
  transmute(latitude = Y, longitude = X)
coords$srd_id_3km <- srd_pts[, c("x", "y")] |>
  st_drop_geometry() |>
  cellFromXY(srd3km, xy = _)
srd_pts <- bind_cols(srd_pts, coords) |>
  st_drop_geometry() |>
  mutate(year = end_year) |>
  select(srd_id, srd_id_3km, year, x, y, latitude, longitude) |>
  arrange(srd_id)
write_parquet(srd_pts, f_pq)

# inputs for feature assignment
srd_static_pq <- path(inputs_dir, "srd_3km_tw_static.parquet")
srd_pts |>
  mutate(row_id = srd_id, latitude, longitude, .keep = "none") |>
  write_parquet(srd_static_pq)
srd_year_pq <- path(inputs_dir, "srd_3km_tw_year.parquet")
srd_pts |>
  mutate(row_id = srd_id, year, latitude, longitude, .keep = "none") |>
  write_parquet(srd_year_pq)


# erd checklists ----

checklists <- path(erd_dir, "erd_checklists.parquet") |>
  open_dataset() |>
  filter(between(latitude, bb$ymin, bb$ymax),
         between(longitude, bb$xmin, bb$xmax),
         between(year, start_year, end_year),
         country_code == "TW",
         cci_es >= 0.9, effort_distance_km < 10,
         is_primary_observer, all_obs_reported) |>
  select(checklist_id, observer_id,
         loc_id, latitude, longitude,
         obs_dt, year, day_of_year, hours_of_day, solar_noon_diff,
         effort_hours = effort_hrs,
         effort_distance_km, effort_speed_kmph,
         number_observers = num_observers) |>
  collect()
f_pq <- path(outputs_dir, "erd_checklists_tw.parquet")
write_parquet(checklists, f_pq)

# inputs for feature assignment
# static
erd_static_pq <- path(inputs_dir, "erd_tw_static.parquet")
erd_static <- checklists |>
  group_by(row_id = loc_id) |>
  summarize(latitude_sd = sd(latitude), longitude_sd = sd(longitude),
            latitude = mean(latitude), longitude = mean(longitude),
            .groups = "drop")
problems <- erd_static |>
  filter(latitude_sd > 0.001 | longitude_sd > 0.001)
stopifnot(nrow(problems) == 0)
erd_static |>
  select(row_id, latitude, longitude) |>
  write_parquet(erd_static_pq)
# yearly
erd_year_pq <- path(inputs_dir, "erd_tw_year.parquet")
erd_year <- checklists |>
  group_by(row_id = loc_id, year) |>
  summarize(latitude_sd = sd(latitude), longitude_sd = sd(longitude),
            latitude = mean(latitude), longitude = mean(longitude),
            .groups = "drop")
problems <- erd_year |>
  filter(latitude_sd > 0.001 | longitude_sd > 0.001)
stopifnot(nrow(problems) == 0)
erd_year |>
  select(row_id, year, latitude, longitude) |>
  write_parquet(erd_year_pq)


# erd observations ----

obs_pq <- path(erd_dir, "erd_obs.parquet") |>
  open_dataset()
obs <- NULL
id <- min(checklists$checklist_id) |> as.numeric()
batch_size <- 5e7
while (id < max(checklists$checklist_id)) {
  message(id)
  obs <- obs_pq |>
    filter(between(checklist_id, id, id + batch_size)) |>
    collect() |>
    semi_join(checklists, by = "checklist_id") |>
    filter(only_presence_reported | obs_count > 0,
           valid, !only_slash_reported) |>
    mutate(observation_count = ifelse(only_presence_reported, NA, obs_count)) |>
    select(checklist_id, species_code, observation_count) |>
    bind_rows(obs)
  id <- id + batch_size + 1
}
f_pq <- path(outputs_dir, "erd_observations_tw.parquet")
write_parquet(obs, f_pq)

# transform to wide format
obs_wide <- obs |>
  mutate(obs_count = case_when(
    only_slash_reported == 1 | valid == 0 ~ -1L,
    only_presence_reported == 1 ~ -2L,
    .default = as.integer(obs_count)
  )) |>
  select(checklist_id, species_code, obs_count) |>
  group_by(species_code) |>
  filter(n() > 100) |>
  ungroup() |>
  arrange(species_code, checklist_id) |>
  pivot_wider(names_from = "species_code",
              values_from = "obs_count",
              values_fill = 0) |>
  mutate(across(-checklist_id, ~ ifelse(.x == -2, NA_integer_, .x)))
f_pq <- path(outputs_dir, "erd_tw_observations_wide.parquet")
write_parquet(obs_wide, f_pq)


# run assignments on bridges2 ----

# crop elevation, run locally
# sinu_crs <- crs("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")
# e <- project(r, sinu_crs) |> ext() |> as.list()
# elev_dir <- "/ocean/projects/deb200005p/sligocki/erd_prep/covariates/elevation"
# glue("gdal_translate -projwin {e$xmin} {e$ymax} {e$xmax} {e$ymin} ",
#      "{elev_dir}/astgtm_sinu_full.tif {elev_dir}/astgtm_sinu.tif")

# generate covariate list, run locally
# jsonlite::read_json("erd/assign_covars/covariate-definitions.json") |>
#   names() |>
#   purrr::discard(stringr::str_detect, pattern = "evi|ocean|weather") |>
#   sort() |>
#   writeLines("covariates.txt")

# interact -n 64
# module load anaconda3
# conda activate conda-env
# for x in $(cat covariates.txt); do
#   python3 schedule_assign_covars.py 2024 erd $x 4 --label=erd_tw --num-procs=8 --no-canary
# done
# sq | grep covar | awk '{print $1}' | sort | uniq | join_str ,
# for x in $(cat covariates.txt); do
#   python3 schedule_assign_covars.py 2024 srd_3km $x 1 --label=srd_3km_tw --num-procs=8 --no-canary
# done
# sq | grep covar | awk '{print $1}' | sort | uniq | join_str ,

# combine shards
# interact -n 64
# module load anaconda3
# conda activate conda-env
# for x in $(cat covariates.txt); do
#   python3 finish_covar.py 2024 erd_tw $x 4
# done

# for x in $(cat covariates.txt); do
#   python3 finish_covar.py 2024 srd_3km_tw $x 1
# done


# merge features into srd ----

srd <- path(outputs_dir, "srd_1km_tw_cell-coordinates.parquet") |>
  read_parquet()
srd_ass_dir <- path(erd_dir, "covar_assignment", "srd_3km_tw")
covs <- dir_ls(srd_ass_dir) |> basename()
for (cov in covs) {
  message(cov)
  assignment <- path(srd_ass_dir, cov, "assignment.parquet") |>
    read_parquet() |>
    rename(srd_id = row_id)
  if (cov_def[[cov]][["time_res"]] == "static") {
    srd <- left_join(srd, assignment, by = "srd_id")
  } else {
    srd <- left_join(srd, assignment, by = c("srd_id", "year"))
  }
}
stopifnot(complete.cases(srd))
srd <- srd |>
  filter(is_on_land) |>
  select(cell_id = srd_id, x, y, latitude, longitude,
         eastness_90m_median, eastness_90m_sd,
         northness_90m_median, northness_90m_sd,
         elevation_30m_median, elevation_30m_sd,
         astwbd_c1_ed, astwbd_c1_pland,
         astwbd_c2_ed, astwbd_c2_pland,
         gsw_c2_pland, gsw_c2_ed,
         ntl_mean, ntl_sd,
         road_density_motorway, road_density_primary,
         road_density_secondary, road_density_tertiary,
         road_density_local, road_density_nonmotorized,
         mcd12q1_lccs1_c1_ed, mcd12q1_lccs1_c1_pland,
         mcd12q1_lccs1_c2_ed, mcd12q1_lccs1_c2_pland,
         mcd12q1_lccs1_c11_ed, mcd12q1_lccs1_c11_pland,
         mcd12q1_lccs1_c12_ed, mcd12q1_lccs1_c12_pland,
         mcd12q1_lccs1_c13_ed, mcd12q1_lccs1_c13_pland,
         mcd12q1_lccs1_c14_ed, mcd12q1_lccs1_c14_pland,
         mcd12q1_lccs1_c15_ed, mcd12q1_lccs1_c15_pland,
         mcd12q1_lccs1_c16_ed, mcd12q1_lccs1_c16_pland,
         mcd12q1_lccs1_c21_ed, mcd12q1_lccs1_c21_pland,
         mcd12q1_lccs1_c22_ed, mcd12q1_lccs1_c22_pland,
         mcd12q1_lccs1_c31_ed, mcd12q1_lccs1_c31_pland,
         mcd12q1_lccs1_c32_ed, mcd12q1_lccs1_c32_pland,
         mcd12q1_lccs1_c41_ed, mcd12q1_lccs1_c41_pland,
         mcd12q1_lccs1_c42_ed, mcd12q1_lccs1_c42_pland,
         mcd12q1_lccs1_c43_ed, mcd12q1_lccs1_c43_pland,
         mcd12q1_lccs1_c255_ed, mcd12q1_lccs1_c255_pland,
         mcd12q1_lccs2_c25_ed, mcd12q1_lccs2_c25_pland,
         mcd12q1_lccs2_c35_ed, mcd12q1_lccs2_c35_pland,
         mcd12q1_lccs2_c36_ed, mcd12q1_lccs2_c36_pland,
         mcd12q1_lccs3_c27_ed, mcd12q1_lccs3_c27_pland,
         mcd12q1_lccs3_c50_ed, mcd12q1_lccs3_c50_pland,
         mcd12q1_lccs3_c51_ed, mcd12q1_lccs3_c51_pland,
         lake_pland,
         lake_dor_mean,  lake_dor_sd,
         lake_volume_mean, lake_volume_sd,
         lake_depth_mean, lake_depth_sd,
         river_dor_mean, river_dor_sd,
         river_area_mean, river_area_sd,
         river_volume_mean, river_volume_sd,
         river_logqavg_mean, river_logqavg_sd,
         river_logqvar_mean, river_logqvar_sd,
         river_power_mean, river_power_sd, river_density)
write_parquet(srd, path(outputs_dir, "prediction-grid_1km_tw.paquet"))


# merge features into erd ----

erd <- path(outputs_dir, "erd_checklists_tw.parquet") |>
  read_parquet()
erd_ass_dir <- path(erd_dir, "covar_assignment", "erd_tw")
covs <- dir_ls(erd_ass_dir) |> basename()
for (cov in covs) {
  message(cov)
  assignment <- path(erd_ass_dir, cov, "assignment.parquet") |>
    read_parquet() |>
    rename(loc_id = row_id)
  if (cov_def[[cov]][["time_res"]] == "static") {
    erd <- left_join(erd, assignment, by = "loc_id")
  } else {
    erd <- left_join(erd, assignment, by = c("loc_id", "year"))
  }
}
stopifnot(complete.cases(erd))
erd <- erd |>
  filter(is_on_land) |>
  select(checklist_id, observer_id,
         loc_id, latitude, longitude,
         obs_dt, year, day_of_year, hours_of_day, solar_noon_diff,
         effort_hours, effort_distance_km, effort_speed_kmph, number_observers,
         eastness_90m_median, eastness_90m_sd,
         northness_90m_median, northness_90m_sd,
         elevation_30m_median, elevation_30m_sd,
         astwbd_c1_ed, astwbd_c1_pland,
         astwbd_c2_ed, astwbd_c2_pland,
         gsw_c2_pland, gsw_c2_ed,
         ntl_mean, ntl_sd,
         road_density_motorway, road_density_primary,
         road_density_secondary, road_density_tertiary,
         road_density_local, road_density_nonmotorized,
         mcd12q1_lccs1_c1_ed, mcd12q1_lccs1_c1_pland,
         mcd12q1_lccs1_c2_ed, mcd12q1_lccs1_c2_pland,
         mcd12q1_lccs1_c11_ed, mcd12q1_lccs1_c11_pland,
         mcd12q1_lccs1_c12_ed, mcd12q1_lccs1_c12_pland,
         mcd12q1_lccs1_c13_ed, mcd12q1_lccs1_c13_pland,
         mcd12q1_lccs1_c14_ed, mcd12q1_lccs1_c14_pland,
         mcd12q1_lccs1_c15_ed, mcd12q1_lccs1_c15_pland,
         mcd12q1_lccs1_c16_ed, mcd12q1_lccs1_c16_pland,
         mcd12q1_lccs1_c21_ed, mcd12q1_lccs1_c21_pland,
         mcd12q1_lccs1_c22_ed, mcd12q1_lccs1_c22_pland,
         mcd12q1_lccs1_c31_ed, mcd12q1_lccs1_c31_pland,
         mcd12q1_lccs1_c32_ed, mcd12q1_lccs1_c32_pland,
         mcd12q1_lccs1_c41_ed, mcd12q1_lccs1_c41_pland,
         mcd12q1_lccs1_c42_ed, mcd12q1_lccs1_c42_pland,
         mcd12q1_lccs1_c43_ed, mcd12q1_lccs1_c43_pland,
         mcd12q1_lccs1_c255_ed, mcd12q1_lccs1_c255_pland,
         mcd12q1_lccs2_c25_ed, mcd12q1_lccs2_c25_pland,
         mcd12q1_lccs2_c35_ed, mcd12q1_lccs2_c35_pland,
         mcd12q1_lccs2_c36_ed, mcd12q1_lccs2_c36_pland,
         mcd12q1_lccs3_c27_ed, mcd12q1_lccs3_c27_pland,
         mcd12q1_lccs3_c50_ed, mcd12q1_lccs3_c50_pland,
         mcd12q1_lccs3_c51_ed, mcd12q1_lccs3_c51_pland,
         lake_pland,
         lake_dor_mean,  lake_dor_sd,
         lake_volume_mean, lake_volume_sd,
         lake_depth_mean, lake_depth_sd,
         river_dor_mean, river_dor_sd,
         river_area_mean, river_area_sd,
         river_volume_mean, river_volume_sd,
         river_logqavg_mean, river_logqavg_sd,
         river_logqvar_mean, river_logqvar_sd,
         river_power_mean, river_power_sd, river_density)
write_parquet(erd, path(outputs_dir, "erd_checklists_1km_tw.paquet"))

# filter obs
f_obs <- path(outputs_dir, "erd_observations_tw.parquet")
read_parquet(f_obs) |>
  semi_join(erd, by = "checklist_id") |>
  write_parquet(f_obs)


# gis data ----

library(dplyr)
library(fs)
library(rnaturalearth)
library(sf)
library(smoothr)
library(wdpar)

f_gpkg <- "data/gis-data.gpkg"
if (file_exists(f_gpkg)) file_delete(f_gpkg)

# bounding box
bb <- st_bbox(c(xmin = 119.1, xmax = 122.4, ymin = 21.7, ymax = 25.4), crs = 4326) |>
  st_as_sfc() |>
  densify(n = 1000) |>
  write_sf(dsn = f_gpkg, layer = "bounding_box")

clip_region <- st_bbox(c(xmin = 30, xmax = 179.9, ymin = -35, ymax = 70),
                       crs = 4326) |>
  st_as_sfc() |>
  densify(max_distance = 10000)

# land boundary
ne_land <- ne_download(scale = 10, category = "physical",
                       type = "land",
                       returnclass = "sf") |>
  st_combine() |>
  st_make_valid() |>
  st_intersection(clip_region) |>
  st_as_sf() |>
  mutate(id = "land")
write_sf(ne_land, dsn = f_gpkg, layer = "land")

# country lines
ne_country_lines <- ne_download(scale = 10, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") |>
  transmute(id = row_number()) |>
  st_combine() |>
  st_make_valid() |>
  st_intersection(clip_region) |>
  st_as_sf() |>
  mutate(id = "country_lines")
write_sf(ne_country_lines, dsn = f_gpkg, layer = "country_lines")

# country
ne_country <- ne_download(scale = 10, category = "cultural",
                          type = "admin_0_countries",
                          returnclass = "sf") |>
  filter(NAME == "Taiwan") |>
  select(name = NAME)
write_sf(ne_country, dsn = f_gpkg, layer = "countries")

# protected areas for Taiwan
tw_pa <- wdpa_fetch("Taiwan", wait = TRUE,
                    download_dir = rappdirs::user_data_dir("wdpar"))
tw_pa_clean <- wdpa_clean(tw_pa) |>
  filter(MARINE == "terrestrial", DESIG_ENG == "National Park") |>
  select(name = NAME)
write_sf(tw_pa_clean, dsn = f_gpkg, layer = "national_parks")
