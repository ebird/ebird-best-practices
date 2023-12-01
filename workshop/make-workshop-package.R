library(fs)
library(glue)
library(usethis)

name <- "ebird-workshop_uwyo-2023"

# create rstudio project
proj_dir <- path(tempdir(), name)
dir_create(proj_dir)
create_project(proj_dir, rstudio = TRUE, open = FALSE)
dir_delete(path(proj_dir, "R"))

# copy scripts to the project
scripts <- c("ebird.R", "envvar.R", "relative-abundance.R",
             "sagthr-exercise.R")
file_copy(path("workshop", scripts), proj_dir)

# copy data
dir_create(path(proj_dir, "data-raw"))
dir_create(path(proj_dir, "data"))
dir_create(path(proj_dir, "ebird-downloads"))
# main content raw data
files <- c("elevation_gmted_1km_us-ga.tif",
           "landcover_mcd12q1_umd_us-ga_2014-2022.tif",
           "mcd12q1_umd_classes.csv")
file_copy(path("data-raw", files), path(proj_dir, "data-raw"))
# main content data
file_copy(dir_ls("data"), path(proj_dir, "data"))
# ebird downloads
file_copy("data-figures-prep/ebd_US-GA_woothr_smp_relOct-2023.zip",
          path(proj_dir, "ebird-downloads"))

# exercise data
files <- c("environmental-variables_checklists_jun-jul_us-wy.csv",
           "environmental-variables_prediction-grid_us-wy.csv",
           "prediction-grid_us-wy.tif",
           "environmental-variables_checklists_jan-dec_pa.csv",
           "environmental-variables_prediction-grid_pa.csv",
           "prediction-grid_pa.tif")
file_copy(path("workshop", "data", files),
          path(proj_dir, "data"))
# ebird downloads
files <- c("ebd_US-WY_larbun_smp_relOct-2023.zip",
           "ebd_US-WY_sagthr_smp_relOct-2023.zip",
           "ebd_PA_blfant1_smp_relOct-2023.zip")
file_copy(path("workshop", "ebird-downloads", files),
          path(proj_dir, "ebird-downloads"))

# zip files
oldwd <- setwd(dirname(proj_dir))
system(glue("zip -r {path(oldwd, 'workshop', paste0(name, '.zip'))} ",
            "{basename(proj_dir)}"))
setwd(oldwd)
dir_delete(proj_dir)
