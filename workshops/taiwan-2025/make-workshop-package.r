library(fs)
library(glue)
library(usethis)

label <- "taiwan-2025"
name <- glue("ebird-workshop_{label}")

# create rstudio project
proj_dir <- path(tempdir(), name)
dir_create(proj_dir)
create_project(proj_dir, rstudio = TRUE, open = FALSE)
dir_delete(path(proj_dir, "R"))

# copy scripts to the project
scripts <- c("ebirdst.r", "ebirdst_complete.r",
             "relative-abundance.r",  "relative-abundance_complete.r")
file_copy(path("workshops", label, scripts), proj_dir)

# copy data
for (dir in c("data", "ebirdst-data")) {
  dir_copy(path("workshops", label, dir), proj_dir)
}

# zip files
oldwd <- setwd(dirname(proj_dir))
system(glue("zip -r {path(oldwd, 'workshops', paste0(name, '.zip'))} ",
            "{basename(proj_dir)}"))
setwd(oldwd)
dir_delete(proj_dir)
