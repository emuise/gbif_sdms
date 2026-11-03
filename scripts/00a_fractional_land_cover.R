library(tidyverse)
library(terra)
library(curl)

covariate_path <- here::here("data", "covariates")
fs::dir_create(dirname(covariate_path))
scratch <- here::here("scratch")
fs::dir_create(scratch)


# this function downloads a zipped file to scratch, and then unzips it into the covariates folder
# it only extracts rasters (tif files)
download_unzip <- function(url) {
  zip_path <- here::here(scratch, basename(url))
  fs::dir_create(dirname(zip_path))

  if (!file.exists(zip_path)) {
    message(glue::glue("Downloading {basename(url)}..."))
    part_path <- paste0(zip_path, ".part")
    res <- curl::multi_download(
      url,
      destfile = part_path,
      resume = TRUE,
      progress = TRUE
    )
    stopifnot(res$success, res$status_code %in% c(200, 206))
    file.rename(part_path, zip_path)
  } else {
    message(glue::glue("{basename(url)} already downloaded. Skipping."))
  }

  uz_fold <- here::here(
    covariate_path,
    url %>% basename() %>% tools::file_path_sans_ext()
  )
  fs::dir_create(uz_fold)

  all_names <- utils::unzip(zip_path, list = TRUE)$Name
  tif_files <- grep("\\.tif$", all_names, value = TRUE)

  existing <- fs::dir_ls(uz_fold, type = "file")
  if (length(existing) >= length(tif_files) && length(tif_files) > 0) {
    message(glue::glue("{basename(uz_fold)} already extracted. Skipping."))
  } else {
    message(glue::glue(
      "Extracting {length(tif_files)} raster(s) to {uz_fold}..."
    ))
    utils::unzip(
      zip_path,
      files = tif_files,
      exdir = uz_fold,
      junkpaths = TRUE
    )
  }

  message(glue::glue("Done: {basename(url)}"))

  fs::dir_ls(here::here(uz_fold))
}

# snap raster to project all other rasters to
snap <- here::here(covariate_path, "Normal_1991_2020_bioclim") %>%
  fs::dir_ls() %>%
  head(1) %>%
  rast()

# land cover
# from here https://www.cec.org/north-american-environmental-atlas/land-cover-30m-2020/

cec_classes <- c(
  "conifers",
  "taiga",
  "evergreen",
  "tropical",
  "deciduous",
  "mixed",
  "tropical_shrubland",
  "temperate_shrubland",
  "tropical_grassland",
  "temperate_grassland",
  "polar_shrubland",
  "polar_grassland",
  "polar_barren",
  "wetland",
  "cropland",
  "barren",
  "urban",
  "water",
  "snow"
)

url <- "https://www.cec.org/files/atlas_layers/1_terrestrial_ecosystems/1_01_0_land_cover_2020_30m/land_cover_2020v2_30m_tif.zip"

lc_30 <- download_unzip(url) %>%
  rast()

tmp <- Sys.getenv("SLURM_TMPDIR")
if (nzchar(tmp)) {
  terraOptions(tempdir = tmp)
}
terraOptions(threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

out_dir <- here::here(covariate_path, "landcover")
fs::dir_create(out_dir)
# ok this makes a fractional cover metric for all of the 19 classes
k <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
message(cec_classes[k])

savepath <- here::here(out_dir, glue::glue("lc_frac_{cec_classes[k]}.tif"))
partpath <- paste0(tools::file_path_sans_ext(savepath), ".part.tif")

if(file.exists(partpath)) {
  fs::file_move(partpath, savepath)
}

if (file.exists(savepath)) {
  message("already done, skipping")
} else {
  rcl <- cbind(0:19, as.numeric(0:19 == k))
  rcl[1, 2] <- NA

  message(Sys.time(), " classify start")
  r <- classify(lc_30, rcl)

  message(Sys.time(), " project start")
  project(
    r, snap,
    method = "average", threads = TRUE,
    filename = partpath,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
  )
  message(Sys.time(), " project done")

  fs::file_move(partpath, savepath)
}
message("done")

fs::dir_ls(here::here("data", "covariates", "landcover")) %>%
  map(\(x) {
    fs::file_move(x, str_remove(x, ".part"))
  })
