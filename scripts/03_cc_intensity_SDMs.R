library(arrow)
library(tidyverse)
library(terra)
library(sf)
library(bcmaps)
library(gbm3)
library(maxnet)

# set paths and mkdir

clean_dir <- here::here("data", "gbif_noram")
cellcount_dir <- here::here("data", "cellcounts")
fs::dir_create(cellcount_dir)

species_list <- here::here("data", "bc_species", "species_list.rds") %>%
  read_rds()

c_files <- fs::dir_ls(clean_dir, recurse = TRUE, type = "file")

bcb <- bcmaps::bc_bound_hres() %>%
  vect()

# this is downloaded in 00
# i have double checked that all climate layers use the same mask
bioclim_fold <- here::here("data", "climate", "Normal_1991_2020_bioclim")

snap <- bioclim_fold %>%
  fs::dir_ls(type = "file") %>%
  head(1) %>%
  rast()

snap_vals <- terra::values(snap)[, 1]

# keep birds since 1990, keep all other observations since 1970

walk(seq_along(c_files), \(n) {
  x <- c_files[[n]]
  savename <- here::here(cellcount_dir, basename(x))
  if (file.exists(savename)) {
    return()
  }

  message(glue::glue("[{n}/{length(c_files)}]: {basename(x)}"))

  df <- open_dataset(x) %>%
    filter(
      (class == "Aves" & year >= 1990) |
        ((is.na(class) | class != "Aves") & year >= 1970)
    ) %>%
    select(species, decimallatitude, decimallongitude) %>%
    collect()

  if (nrow(df) == 0) {
    # probably should save an empty file
    return()
  }

  pts <- terra::vect(
    df,
    geom = c("decimallongitude", "decimallatitude"),
    crs = "epsg:4326"
  ) %>%
    project(snap)

  df$cell <- terra::cellFromXY(snap, terra::crds(pts))

  cc <- df %>%
    filter(!is.na(cell)) %>% # this removes anything outside the raster bounds
    mutate(snap_val = snap_vals[cell]) %>%
    filter(!is.na(snap_val)) %>% # this only keeps values that have matching climate values
    select(-snap_val) %>%
    count(species, cell, name = "n") %>%
    arrow::write_parquet(savename)
})
# 35 and 1494 evidently have nothing in them
# 35 is all birds prior to 1982
# 1494 is some plants prior to 1940

cc <- arrow::open_dataset(cellcount_dir)

cellcounts <- cc %>%
  group_by(species, cell) %>%
  summarise(n = sum(n))

sample_intensity <- cc %>%
  group_by(cell) %>%
  summarize(n = sum(n)) %>%
  collect()

r <- snap
values(r) <- NA

r[!is.nan(snap_vals)] <- 0
r[sample_intensity$cell] <- sample_intensity$n

smooth_r <- focal(r, w = 15, fun = "mean", na.rm = TRUE, na.policy = "omit")
names(smooth_r) <- "sample_intensity"

l1p_smooth_r <- log1p(smooth_r)

plot(
  l1p_smooth_r,
  main = "log(mean sampling intensity) (30 km focal window)"
)

bc_smooth <- smooth_r %>%
  crop(
    bcb %>%
      project(smooth_r),
    mask = TRUE
  ) %>%
  project(bcb)

plot(log(bc_smooth + 1))

vars <- c(
  "MAP", # mean annual precipitation (mm)
  "DD_0", # chilling degree days (Degree days below 0 °C)
  "PAS", # precipitation as snow (mm)
  "CMD", # Hargreave's climatic moisture index
  "DD18" # warming degree days above 18 °C
)

bioclim <- here::here(
  bioclim_fold,
  glue::glue("Normal_1991_2020_{vars}.tif")
) %>%
  rast()

names(bioclim) <- vars

covariates <- c(bioclim, l1p_smooth_r)

set.seed(42)

psabs <- spatSample(
  covariates,
  10000,
  values = FALSE,
  na.rm = TRUE,
  cells = TRUE
) %>%
  as.numeric()

psabs_tib <- tibble(
  species = "Pseudo absence",
  cell = psabs,
  n = 0
)


run_sdm <- function(
  spec,
  cellcounts,
  covariates,
  psabs_tib,
  n_presence = 5000
) {
  spec_df <- cellcounts %>%
    filter(species == spec) %>%
    collect()

  if (nrow(spec_df) == 0) {
    message(glue::glue("Skipping {spec}: no cell records."))
    return(NULL)
  }

  spec_psab <- spec_df %>%
    slice_sample(n = min(n_presence, nrow(spec_df))) %>%
    bind_rows(psabs_tib) %>%
    ungroup()

  coords <- xyFromCell(covariates, spec_psab$cell)

  pa <- spec_psab %>%
    mutate(
      x = coords[, 1],
      y = coords[, 2],
      covariates[cell]
    ) %>%
    arrange(species)

  for_model <- pa %>%
    mutate(pa = as.numeric(n > 0)) %>%
    select(names(covariates), pa) %>%
    drop_na()

  if (n_distinct(for_model$pa) < 2) {
    message(glue::glue(
      "Skipping {spec}: presence/absence not both present after NA drop."
    ))
    return(NULL)
  }
  # # testing glm
  # model <- gbm::gbm(pa ~ ., data = for_model, family = binomial())

  # # maxent

  model <- maxnet(
    p = for_model$pa,
    data = for_model %>% select(-pa)
  )

  # # brt
  train_params <- training_params(
    num_trees = 2000,
    shrinkage = 0.01,
    interaction_depth = 3,
    bag_fraction = 0.5,
    num_train = round(0.8 * nrow(for_model)), # or nrow() for no held-out split
    min_num_obs_in_node = 10
  )

  model <- gbmt(
    pa ~ .,
    data = for_model,
    distribution = gbm_dist("Bernoulli"),
    train_params = train_params,
    cv_folds = 5,
    is_verbose = FALSE
  )

  best_iter <- gbmt_performance(model, method = "cv")

  pred_covariates <- covariates
  pred_covariates[["sample_intensity"]] <- max(
    pa$sample_intensity,
    na.rm = TRUE
  )

  prediction <- predict(
    pred_covariates,
    model,
    n.trees = best_iter,
    type = "response"
  )

  # for maxent
  prediction <- predict(pred_covariates, model, type = "cloglog", na.rm = T)
  # for uhhhhh glm
  prediction <- predict(pred_covariates, model, type = "response")

  list(species = spec, model = model, prediction = prediction, pa = pa)
}

# ---------------------------------------------------------------------------
# Example: single species (matches original script's demo)
# ---------------------------------------------------------------------------

spec <- "Gulo gulo"
result <- run_sdm(spec, cc, covariates, psabs_tib)

if (!is.null(result)) {
  plot(result$prediction, main = spec)
  plot(
    result$pa %>%
      filter(n > 0) %>%
      vect(geom = c("x", "y")),
    add = TRUE
  )
}
