library(tidyverse)
library(minioclient)

# one-time setup (safe to re-run - config persists in minioclient's config dir)
install_mc()
mc_alias_set(
  "aw",
  endpoint = "s3-us-west-2.amazonaws.com",
  access_key = "",
  secret_key = ""
)

aw_url <- "https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6v73/normals/Normal_1991_2020_bioclim.zip"
aw_mc_path <- file.path("aw", sub("^https?://[^/]+/", "", aw_url))

aw_zip_path <- here::here("data", "climate", basename(aw_url))
fs::dir_create(dirname(aw_zip_path))

if (!file.exists(aw_zip_path)) {
  message(glue::glue("Syncing {basename(aw_url)} via mc..."))
  mc_cp(aw_mc_path, aw_zip_path)
} else {
  message(glue::glue("{basename(aw_url)} already downloaded. Skipping."))
}


aw_uz_fold <- tools::file_path_sans_ext(aw_zip_path)
fs::dir_create(aw_uz_fold)

all_names <- utils::unzip(aw_zip_path, list = TRUE)$Name
tif_files <- grep("\\.tif$", all_names, value = TRUE)

existing <- fs::dir_ls(aw_uz_fold, type = "file")
if (length(existing) >= length(tif_files) && length(tif_files) > 0) {
  message(glue::glue("{basename(aw_uz_fold)} already extracted. Skipping."))
} else {
  message(glue::glue(
    "Extracting {length(tif_files)} raster(s) to {aw_uz_fold}..."
  ))
  utils::unzip(
    aw_zip_path,
    files = tif_files,
    exdir = aw_uz_fold,
    junkpaths = TRUE
  )
}

message("Done: current climate normals downloaded and extracted.")
