# spec_tbl is a tibble with columns kingdom, species (one row per taxon)

calc_group_intensity <- function(spec_tbl, counts = cc) {
  message("calculating group intensity")
  counts %>%
    semi_join(spec_tbl, by = c("kingdom", "species")) %>%
    group_by(cell) %>%
    summarize(n = sum(n)) %>%
    collect() %>%
    mutate(n = as.numeric(n))
}

intensity_to_raster <- function(intensity_df, template = snap, template_vals = snap_vals) {
  message("rasterizing")
  r <- template
  values(r) <- NA
  r[!is.nan(template_vals)] <- 0
  r[intensity_df$cell] <- intensity_df$n
  r
}

smooth_intensity <- function(intensity_r, w = 31) {
  message("smoothing")
  focal(intensity_r, w = w, fun = "mean", na.rm = TRUE, na.policy = "omit")
}

pipe_save_raster <- function(r, path) {
  names(r) <- path %>% basename() %>% tools::file_path_sans_ext()
  writeRaster(r, path, overwrite = TRUE)
  r
}

calc_smooth_intensity <- function(
  spec_tbl,
  name,
  w = 31,
  dir = here::here(covariate_path, "sample_intensity")
) {
  message("--- ", name, " ---")
  fs::dir_create(dir)

  raw_path        <- here::here(dir, glue::glue("{name}-raw.tif"))
  smooth_path     <- here::here(dir, glue::glue("{name}-smooth_{w}.tif"))
  smooth_l1p_path <- here::here(dir, glue::glue("{name}-smooth_{w}-l1p.tif"))

  if (all(file.exists(raw_path, smooth_path, smooth_l1p_path))) {
    message(glue::glue("{name}_w-{w} already processed"))
    return(rast(smooth_l1p_path))
  }

  calc_group_intensity(spec_tbl) %>%
    intensity_to_raster() %>%
    pipe_save_raster(raw_path) %>%
    smooth_intensity(w = w) %>%
    pipe_save_raster(smooth_path) %>%
    log1p() %>%
    pipe_save_raster(smooth_l1p_path)
}