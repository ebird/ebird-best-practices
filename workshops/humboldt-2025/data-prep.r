library(arrow)
library(auk)
library(dplyr)
library(fs)
library(glue)
library(lubridate)
library(readr)
library(rnaturalearth)
library(sf)
library(stringr)
library(terra)
setwd(here::here("workshop/humboldt-2025/"))

# prepare data for covariate assignment ----

# ebird
checklists <- read_sampling("data-raw/ebd_CO_smp_relAug-2025_sampling.txt") |>
  filter(all_species_reported,
         observation_type %in% c("Traveling", "Stationary")) |>
  group_by(row_id = locality_id, year = year(observation_date)) |>
  summarize(latitude = first(latitude), longitude = first(longitude),
            .groups = "drop")

# prediction grid
# colombia centered albers equal area conic
aeac_crs <- st_crs("+proj=aea +lat_1=-4 +lat_2=8 +lat_0=2 +lon_0=-73")

# study region
study_region <- ne_countries(scale = 10, country = "Colombia")
# clip out islands
clip_bb <- st_bbox(c(xmin = -79.4, ymin = -4.2, xmax = -66.9, ymax = 12.6),
                   crs = 4326) |>
  st_as_sfc()
study_region <- st_intersection(study_region, clip_bb) |>
  st_transform(crs = aeac_crs) |>
  st_geometry() |>
  vect()

# create a raster template covering the region with 3 km resolution
r <- study_region |>
  buffer(3000) |>
  rast(res = c(3000, 3000))

# fill the raster with 1s inside the study region
r <- rasterize(study_region, r, touches = TRUE) |>
  setNames("study_region")

# save for later use
r <- writeRaster(r, filename = "data/prediction-grid_co.tif",
                 overwrite = TRUE,
                 gdal = "COMPRESS=DEFLATE")

# cell coordinates
cell_coordinates <- as.data.frame(r, cells = TRUE, xy = TRUE) |>
  select(cell_id = cell, x, y) |>
  st_as_sf(coords = c("x", "y"), crs = st_crs(r), remove = FALSE) |>
  st_transform(crs = 4326)
cell_coordinates <- st_coordinates(cell_coordinates) |>
  as.data.frame() |>
  select(latitude = Y, longitude = X) |>
  bind_cols(cell_coordinates) |>
  mutate(row_id = paste0("PG", cell_id), year = 2024) |>
  select(cell_id, row_id, year, x, y, latitude, longitude)
write_csv(cell_coordinates, "data-raw/prediction-grid_cell-coordinates_co.csv")

# combine and save for covariate assignment
combined <- bind_rows(
  checklists,
  select(cell_coordinates, row_id, year, latitude, longitude)
)
write_parquet(combined, "data-raw/erd_co_year.parquet")
combined <- checklists |>
  group_by(row_id) |>
  summarize(latitude = first(latitude), longitude = first(longitude)) |>
  bind_rows(select(cell_coordinates, row_id, latitude, longitude))
write_parquet(combined, "data-raw/erd_co_static.parquet")


# assign covariates ----

# crop elevation
# sinu_crs <- crs("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")
# e <- project(r, sinu_crs) |> ext() |> as.list()
# elev_dir <- "/ocean/projects/deb200005p/sligocki/erd_prep/covariates/elevation"
# glue("gdal_translate -projwin {e$xmin} {e$ymax} {e$xmax} {e$ymin} ",
#      "{elev_dir}/astgtm_sinu_full.tif {elev_dir}/astgtm_sinu.tif")

# assign
# interact -n 64
# module load anaconda3
# conda activate conda-env
# for x in elevation slope_1km mountains roads ntl mcd12q1_lccs1; do
#   python3 schedule_assign_covars.py 2024 erd $x 4 --label=erd_co --num-procs=8 --no-canary
# done
# sq | grep covar | awk '{p rint $1}' | sort | uniq | join_str ,

# combine shards
# interact -n 64
# module load anaconda3
# conda activate conda-env
# for x in elevation slope_1km mountains roads ntl mcd12q1_lccs1; do
#   python3 finish_covar.py 2024 erd_co $x 4
# done

# combine
base_dir <- "/ocean/projects/deb200005p/sligocki/erd_prep/erd2024/covar_assignment"
inputs_dir <- path(base_dir, "inputs")
assign_dir <- path(base_dir, "erd_co")
outputs_dir <- "/ocean/projects/deb200005p/mstrimas"

erd_srd <- path(inputs_dir, "erd_co_year.parquet") |>
  read_parquet() |>
  select(row_id, year)

covs <- dir_ls(assign_dir) |> basename()
covs <- c("elevation", "mcd12q1_lccs1", "mountains",
          "ntl", "roads", "slope_1km", "slope_90m")
for (cov in covs) {
  message(cov)
  assignment <- path(assign_dir, cov, "assignment.parquet") |>
    read_parquet()
  if ("year" %in% names(assignment)) {
    erd_srd <- left_join(erd_srd, assignment, by = c("row_id", "year"))
  } else {
    erd_srd <- left_join(erd_srd, assignment, by = "row_id")
  }
}
stopifnot(complete.cases(erd_srd))
write_parquet(erd_srd, path(outputs_dir, "erd-srd_3km_co.paquet"))


# merge features into erd/srd ----

assignments <- read_parquet("data-raw/erd-srd_3km_co.paquet") |>
  select(row_id, year,
         starts_with("elevation"),
         starts_with("northness_90m"),
         starts_with("eastness_90m"),
         mountain,
         starts_with("mcd12q1"),
         starts_with("road")) |>
  select(-mcd12q1_lccs1_diversity)

# swap in descriptive names
lc <- read_csv("data-raw/mcd12q1_lccs1_classes.csv")
old_names <- c(glue("mcd12q1_lccs1_c{lc$class}_pland"),
               glue("mcd12q1_lccs1_c{lc$class}_ed"))
new_names <- c(glue("landcover_{lc$label}_pland"),
               glue("landcover_{lc$label}_ed"))
lookup <- match(names(assignments), old_names)
index <- which(!is.na(lookup))
lookup <- lookup[!is.na(lookup)]
names(assignments)[index] <- new_names[lookup]

# erd
assignments |>
  filter(str_starts(row_id, "L")) |>
  rename(locality_id = row_id) |>
  write_parquet("data/environmental-variables_checklists_co.parquet")

# srd
srd <- read_csv("data-raw/prediction-grid_cell-coordinates_co.csv")
assignments |>
  filter(str_starts(row_id, "PG")) |>
  inner_join(srd, y = _, by = c("row_id", "year")) |>
  select(-row_id, -year) |>
  write_parquet("data/prediction-grid_co.parquet")


# gis data ----

f_gpkg <- "data/gis-data.gpkg"
if (file_exists(f_gpkg)) file_delete(f_gpkg)

# focal region
region <- read_sf("data-raw/RegionFis_Col_WGS84_Andinas.gpkg") |>
  st_transform(crs = 4326) |>
  mutate(region = "Andinas") |>
  select(region) |>
  write_sf(dsn = f_gpkg, layer = "region")

clip_region <- st_geometry(st_buffer(region, 2e6))

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
  filter(ISO_A2 == "CO") |>
  select(name = NAME)
write_sf(ne_country, dsn = f_gpkg, layer = "countries")

# states
ne_state <- ne_download(scale = 10, category = "cultural",
                        type = "admin_1_states_provinces",
                        returnclass = "sf") |>
  filter(iso_a2 == "CO") |>
  select(code = iso_3166_2, name = name)
write_sf(ne_state, dsn = f_gpkg, layer = "departments")
