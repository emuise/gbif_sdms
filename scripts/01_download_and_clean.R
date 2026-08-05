library(gbifdb)
library(arrow)
library(tidyverse)
library(countrycode)
library(geodata)
library(terra)
library(tidyterra)
library(sf) # fast intersect

# set paths and make folder structure
raw_gbif_loc <- here::here("data", "gbif_global")

scratch_dir <- here::here("scratch", "noram_process")
clean_dir <- here::here("data", "gbif_noram")

fs::dir_create(scratch_dir)
fs::dir_create(clean_dir)

run_download <- function(max_retries = Inf) {
  download_flag <- here::here("flags", "download_done.txt")

  if (fs::file_exists(download_flag)) {
    message("Download already marked as complete. Skipping.")
    return(invisible(NULL))
  }

  # minioclient::install_mc()
  attempt <- 1
  while (!fs::file_exists(download_flag) && attempt <= max_retries) {
    message(sprintf("Starting GBIF download (Attempt %d)...", attempt))

    tryCatch(
      {
        gbifdb::gbif_download(dir = raw_gbif_loc)

        fs::dir_create(dirname(download_flag))
        fs::file_create(download_flag)
        message("Download complete.")
      },
      error = function(e) {
        message(sprintf(
          "Download failed on attempt %d: %s",
          attempt,
          e$message
        ))
        message("Retrying in 5 seconds...")
        Sys.sleep(5)
      }
    )

    attempt <- attempt + 1
  }
}

run_download(max_retries = 30)

if (!file.exists(here::here("flags", "download_done.txt"))) {
  gbif_files <- fs::dir_ls(raw_gbif_loc, recurse = TRUE, type = "file")
  fs::file_delete(gbif_files[is.na(as.numeric(basename(gbif_files)))])
}

noram <- vect(
  "https://github.com/gbif/continents/raw/master/continent_cookie_cutter.gpkg"
) %>%
  filter(continent_part == "north_america_east")

world <- geodata::world(path = here::here("data", "shps")) %>%
  project(noram)

na_countries <- intersect(world, noram)
na_bbox <- ext(na_countries)

na_cc <- na_countries %>% pull(GID_0)
na_cc_iso2c <- countrycode(na_cc, origin = "iso3c", destination = "iso2c")

ymin <- na_bbox$ymin
ymax <- na_bbox$ymax
xmin <- na_bbox$xmin
xmax <- na_bbox$xmax

# save these for downstream scripts / provenance
saveRDS(
  list(
    na_cc_iso2c = na_cc_iso2c,
    ymin = ymin,
    ymax = ymax,
    xmin = xmin,
    xmax = xmax
  ),
  here::here("data", "noram_extent.rds")
)

# standard GBIF filters for geospatial issues
geo_issues <- c(
  "COORDINATE_REPROJECTION_FAILED",
  "COORDINATE_REPROJECTION_SUSPICIOUS",
  "COORDINATE_UNCERTAINTY_METERS_INVALID",
  "PRESUMED_NEGATED_LATITUDE",
  "FOOTPRINT_WKT_MISMATCH",
  "FOOTPRINT_WKT_INVALID",
  "COUNTRY_COORDINATE_MISMATCH",
  "COORDINATE_PRECISION_INVALID",
  "PRESUMED_NEGATED_LONGITUDE",
  "CONTINENT_COUNTRY_MISMATCH",
  "CONTINENT_COORDINATE_MISMATCH",
  "PRESUMED_SWAPPED_COORDINATE"
)

run_split_grouped <- function() {
  split_flag <- here::here("flags", "split_done.txt")

  if (fs::file_exists(split_flag)) {
    message("Split already marked as complete. Skipping.")
    return(invisible(NULL))
  }

  dones <- fs::dir_ls(scratch_dir)
  files <- fs::dir_ls(raw_gbif_loc, recurse = TRUE, type = "file")
  files_left <- files[!(basename(files) %in% basename(dones))]

  glob_max <- files_left %>% basename() %>% as.numeric() %>% max()

  groups <- split(files_left, ceiling(seq_along(files_left) / 10))

  walk(groups, \(group) {
    tempfiles <- here::here("flags", "noram_process", basename(group))

    remaining <- group[!fs::file_exists(tempfiles)]

    if (length(remaining) == 0) {
      return()
    }

    nums <- remaining %>%
      basename() %>%
      as.numeric()

    template <- glue::glue("{min(nums)}-{max(nums)}")

    message(glue::glue(
      "\n=== Processing Files [{template}] (Max: {glob_max}) ==="
    ))

    df <- remaining %>%
      open_dataset()

    message(glue::glue(
      "  │ Raw records loaded : {format(nrow(df), big.mark = ',')}"
    ))

    filtered_df <- df %>%
      filter(
        occurrencestatus == "PRESENT",
        !basisofrecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"),
        !is.na(species),
        !is.na(decimallatitude),
        !is.na(decimallongitude),
        countrycode %in% na_cc_iso2c | is.na(countrycode),
        decimallatitude < ymax,
        decimallatitude > ymin,
        decimallongitude < xmax,
        decimallongitude > xmin
      ) %>%
      collect()

    message(glue::glue(
      "  │ Post spatial filter: {format(nrow(filtered_df), big.mark = ',')}"
    ))

    filtered_df <- filtered_df %>%
      filter(!purrr::map_lgl(issue, ~ any(.x$array_element %in% geo_issues)))

    message(glue::glue(
      "  │ Post issue filter  : {format(nrow(filtered_df), big.mark = ',')}"
    ))

    if (nrow(filtered_df) > 0) {
      message("  └─ Saving Parquet... ")
      write_dataset(
        filtered_df,
        path = clean_dir,
        format = "parquet",
        existing_data_behavior = "overwrite",
        max_partitions = 20000,
        max_rows_per_file = 1000000,
        basename_template = paste0("part_", template, "_{i}.parquet")
      )
    } else {
      message("  └─ No rows remaining. Skipping write.")
    }

    fs::dir_create(dirname(tempfiles)[[1]])
    fs::file_create(tempfiles)

    rm(filtered_df)
    gc(full = TRUE)
    gc(full = TRUE)
  })

  fs::dir_create(dirname(split_flag))
  fs::file_create(split_flag)
}

run_split_grouped()

if (fs::dir_exists(raw_gbif_loc)) {
  fs::dir_delete(raw_gbif_loc)
}

message("Done: cleaned, filtered parquet written to ", clean_dir)
