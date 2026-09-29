library(tidyverse)
library(minioclient)
library(terra)
library(curl)
library(whitebox)
whitebox::install_whitebox()

# one-time setup (safe to re-run - config persists in minioclient's config dir)
install_mc()
mc_alias_set(
  "aw",
  endpoint = "s3-us-west-2.amazonaws.com",
  access_key = "",
  secret_key = ""
)

covariate_path <- here::here("data", "covariates")
dir.create(dirname(covariate_path))

# climate

download_unzip_bucket <- function(aw_url) {
  aw_mc_path <- file.path("aw", sub("^https?://[^/]+/", "", aw_url))

  aw_zip_path <- here::here(covariate_path, basename(aw_url))
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

  message(glue::glue("Done: {basename(aw_url)}"))
}
# these links are from here: https://adaptwest.databasin.org/pages/adaptwest-climatena/
download_unzip_bucket(
  "https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6v73/normals/Normal_1991_2020_bioclim.zip"
)
download_unzip_bucket(
  "https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6v73/ensembles/ensemble_13GCMs_ssp585_2081_2100_bioclim.zip"
)


# snap raster to project all other rasters to
snap <- here::here(covariate_path, "Normal_1991_2020_bioclim") %>%
  fs::dir_ls() %>%
  head(1) %>%
  rast()

# make it a snap raster of 0s where values exist
snap <- classify(snap, cbind(-Inf, Inf, 0))

snap_bound <- ext(snap) %>%
  vect(crs = snap)

snap_poly <- as.polygons(snap, dissolve = T)

# human footprint index
hfi_path <- fs::dir_create(covariate_path, "hfi")

hfi_download <- function(url) {
  dest <- file.path(hfi_path, basename(url))
  if (file.exists(dest)) {
    message(glue::glue("Already done: {basename(url)}"))
    return(dest)
  }
  scratchfile <- here::here("scratch", basename(url))
  if (!file.exists(scratchfile)) {
    res <- multi_download(
      url,
      destfile = scratchfile,
      resume = TRUE,
      progress = TRUE
    )
    stopifnot(res$success, res$status_code %in% c(200, 206))
  }
  hfi <- rast(scratchfile)
  message("Downloaded")

  snap_proj_file <- here::here("scratch", "snap_hfi.tif")

  if (!file.exists(snap_proj_file)) {
    ply <- snap_poly %>%
      project(hfi)

    r_mask <- rasterize(ply, hfi) %>%
      trim() %>%
      writeRaster(snap_proj_file, overwrite = T)
  }

  snap_proj <- rast(snap_proj_file)

  message("Aligning to climate rasters: crop mask cover project")
  hfi <- hfi %>%
    crop(snap_proj, mask = T) %>%
    # trim() %>%
    cover(snap_proj) %>%
    project(snap, method = "mean", threads = T)

  writeRaster(hfi, dest)
  message(glue::glue("Done: {basename(url)}"))
  return(dest)
}
# url from here https://wcshumanfootprint.org/data-access
# single years
hfi_download(
  "https://storage.googleapis.com/hii-export/2020-01-01/hii_2020-01-01.tif"
)

# all years
glue::glue(
  "https://storage.googleapis.com/hii-export/{c(2015:2020)}-01-01/hii_{c(2015:2020)}-01-01.tif"
) %>%
  map(hfi_download)

# IUCN range maps
# manually download from https://www.iucnredlist.org/resources/spatial-data-download
# unzip into here::here("data", "covariates", "iucn")
# then run this section

iucn_loc <- here::here("data", "covariates", "iucn")

iucn_folds <- fs::dir_ls(here::here(iucn_loc, "raw"))

map(iucn_folds, \(fold) {
  savename <- glue::glue(
    here::here(iucn_loc, fold %>% basename() %>% str_to_lower()),
    ".gpkg"
  )

  if (file.exists(savename)) {
    return(savename)
  }

  print(basename(fold))

  rangemaps <- fs::dir_ls(fold, regexp = ".shp$") %>%
    vect()

  rm_crop <- rangemaps %>%
    crop(
      snap_poly %>%
        project(rangemaps)
    ) %>%
    project(snap_poly)

  writeVector(rm_crop, savename)
})

# DEM and other derived data
download_elev <- function() {
  elev_dir <- here::here(covariate_path, "elev")

  files <- c("dem", "twi", "tri")

  if(all(file.exists(here::here(elev_dir, glue::glue("{files}.tif"))))) {
    message("All DEM based files already downloaded and aligned")
    return(fs::dir_ls(elev_dir))
  }
  fs::dir_create(elev_dir)
  elev <- geodata::elevation_global(0.5, "scatch")

  elev_snap_loc <- here::here("scratch", "snap_elev.tif")

  if (!file.exists(elev_snap_loc)) {
    snap_elev <- snap_poly %>%
      project(elev) %>%
      rasterize(elev) %>%
      trim() %>%
      writeRaster(elev_snap_loc)
  }

  snap_elev <- rast(elev_snap_loc)

  elev_a_loc <- here::here(elev_dir, "dem.tif")

  if (!file.exists(elev_a_loc)) {
    elev_aligned <- elev %>%
      crop(snap_elev, mask = T) %>%
      project(snap, threads = T, method = "average") %>%
      writeRaster(elev_a_loc)
  }

  elev_aligned <- rast(elev_a_loc)

  tri_loc <- here::here(elev_dir, "tri.tif")

  if (!file.exists(tri_loc)) {
    tri <- terrain(elev_aligned, v = "TRI") %>%
      writeRaster(tri_loc)
  }

  twi_loc <- here::here(elev_dir, "twi.tif")

  if (!file.exists(twi_loc)) {
    # breach depressions (carves through barriers rather than filling/flattening)
    breached_file <- here::here("scratch", "elev_breached.tif")
    wbt_breach_depressions_least_cost(
      dem = elev_a_loc,
      output = breached_file,
      dist = 5, # search radius (cells) for the least-cost breach path; increase if pits persist
      fill = TRUE # fill any depressions that couldn't be breached within `dist`
    )

    slope_file <- here::here("scratch", "slope.tif")
    wbt_slope(breached_file, slope_file)

    sca_file <- here::here("scratch", "sca.tif")
    wbt_d_inf_flow_accumulation(
      breached_file,
      sca_file,
      out_type = "Specific Contributing Area"
    )

    wbt_wetness_index(sca_file, slope_file, twi_loc)
  }
  twi <- rast(twi_loc)

  message("Downloaded and aligned elev (climateNA), calculated TRI (terra), and TWI (whiteboxtools)")
  return(fs::dir_ls(elev_dir))
}

download_elev()
