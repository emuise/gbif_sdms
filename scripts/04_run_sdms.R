library(arrow)
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)



species_dir <- here::here("data", "bc_species")
species_list <- readRDS(here::here(species_dir, "species_list.rds"))

# ---------------------------------------------------------------------------
# Full run: loop over all BC species and save predictions
# ---------------------------------------------------------------------------
# Uncomment to run for every species in species_list. Wrapped in tryCatch
# since some species will fail (too few points, non-convergence, etc).

# sdm_out_dir <- here::here("data", "sdm_predictions")
# fs::dir_create(sdm_out_dir)
#
# walk(species_list, \(spec) {
#   outfile <- here::here(sdm_out_dir, glue::glue("{make.names(spec)}.tif"))
#   if (file.exists(outfile)) return()
#
#   res <- tryCatch(
#     run_sdm(spec, cc, covariates, psabs_tib),
#     error = function(e) {
#       message(glue::glue("Failed {spec}: {e$message}"))
#       NULL
#     }
#   )
#
#   if (!is.null(res)) {
#     writeRaster(res$prediction, outfile, overwrite = TRUE)
#   }
# })
