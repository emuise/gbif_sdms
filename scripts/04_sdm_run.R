library(arrow)
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(bcmaps)
library(gbm3)
library(maxnet)
library(viridis)
library(taxize)
muiselib2::muiselib_setup()

# set paths and mkdir

covariate_path <- here::here("data", "covariates")
clean_dir <- here::here("data", "gbif_noram")
cellcount_dir <- here::here("data", "cellcounts")

species_list <- here::here("data", "bc_species", "species_list.rds") %>%
  read_rds()

taxa_join_loc <- here::here("data", "taxa_join.parquet")

if (!file.exists(taxa_join_loc)) {
  taxa_join <- here::here("data", "spatial_counts") %>%
    arrow::open_dataset() %>%
    distinct(kingdom, phylum, class, order, family, genus, species) %>%
    mutate(
      si_group = case_when(
        family %in%
          c(
            "Papilionidae",
            "Hesperiidae",
            "Pieridae",
            "Lycaenidae",
            "Riodinidae",
            "Nymphalidae",
            "Hedylidae"
          ) ~ "butterflies",
        class == "Aves" ~ "birds",
        class == "Mammalia" ~ "mammals",
        class %in%
          c(
            "Reptilia",
            "Testudines",
            "Squamata",
            "Crocodylia",
            "Amphibia"
          ) ~ "reptiles_amphibians",
        kingdom == "Plantae" ~ "plants",
        kingdom == "Fungi" ~ "fungi",
        (class == "Insecta" &
          !order %in%
            c("Ephemeroptera", "Plecoptera", "Trichoptera", "Megaloptera")) |
          class %in% c("Arachnida", "Diplopoda", "Chilopoda", "Collembola") |
          order %in% c("Stylommatophora", "Crassiclitellata") ~ "invertebrates",
        TRUE ~ "all_species"
      )
    )
  
  write_parquet(taxa_join, taxa_join_loc)
}
taxa_join <- read_parquet(taxa_join_loc)

# this is the cellcounts index
# within the sdm function, we identify all files with each
cc_index <- here::here("data", "cellcounts_species_index.parquet") %>%
  read_parquet()

# covariate stac
bc_vars <- c(
  "MAP", # mean annual precipitation (mm)
  "DD_0", # chilling degree days (Degree days below 0 °C)
  "PAS", # precipitation as snow (mm)
  "CMD", # Hargreave's climatic moisture index
  "DD18" # warming degree days above 18 °C
)

bioclim <- here::here(
  covariate_path,
  "Normal_1991_2020_bioclim",
  glue::glue("Normal_1991_2020_{bc_vars}.tif")
) %>%
  rast()

names(bioclim) <- bc_vars

lc <- here::here(
  covariate_path,
  "landcover_input"
) %>%
  fs::dir_ls() %>%
  rast()

topo <- here::here(
  covariate_path,
  "elev"
) %>%
  fs::dir_ls() %>%
  rast()

names(topo) <- c("elev", "tri", "twi")

soil <- here::here(covariate_path, "soil") %>%
  fs::dir_ls(regexp = ".tif$") %>%
  rast()

hfi <- here::here(covariate_path, "hfi", "hii_2020-01-01.tif") %>%
  rast()

set.seed(42)

sample_intensity <- fs::dir_ls(
  here::here(covariate_path, "sample_intensity"),
  regexp = "-l1p.tif$"
) %>%
  rast()
# these are log 1 plus sample intensities
names(sample_intensity) <- str_remove(names(sample_intensity), "-smooth_.*")

covariates <- c(bioclim, lc, topo, soil, hfi, sample_intensity$all_species)

base_covariate_names <- names(covariates)[!(names(covariates) %in% names(sample_intensity))]

# generate and save pseudoabscences

psabs_fold <- here::here("data", "psuedoabsences")
fs::dir_create(psabs_fold)

psabs_tib_loc <- here::here(psabs_fold, "pseudoabsences.csv")

if (!file.exists(psabs_tib_loc)) {
  psabs <- spatSample(
    covariates,
    20000,
    values = FALSE,
    na.rm = TRUE,
    cells = TRUE
  ) %>%
    as.numeric()

  psabs_tib <- tibble(
    species = "Pseudoabsence",
    cell = psabs,
    n = 0
  )

  coords <- xyFromCell(covariates, psabs_tib$cell)

  psabs_tib <- psabs_tib %>%
    mutate(
      x = coords[, 1],
      y = coords[, 2],
      covariates[cell]
    ) %>%
    drop_na() %>%
    slice_sample(n = 10000)

  write_csv(psabs_tib, psabs_tib_loc)
}

psabs_tib <- read_csv(psabs_tib_loc)

# one weighted target-group background per si group, cached
make_tgb <- function(group, n = 10000) {
  
  savename <- here::here(psabs_fold, glue::glue("tgb_{group}.csv"))
  if (file.exists(savename)) {
    return(read_csv(savename, show_col_types = FALSE))
  }
  message(glue::glue("making target group background for {group}"))

  tg_raw <- here::here(covariate_path, "sample_intensity", glue::glue("{group}-raw.tif")) %>%
    rast()

  cells <- spatSample(tg_raw, n * 2, method = "weights", cells = TRUE, na.rm = TRUE)[, "cell"] %>%
    unique()

  coords <- xyFromCell(covariates, cells)

  tgb <- tibble(species = "Pseudoabsence", cell = cells, n = 0) %>%
    mutate(x = coords[, 1], y = coords[, 2], covariates[cell]) %>%
    drop_na() %>%
    slice_sample(n = n)

  write_csv(tgb, savename)
  tgb
}

groups <- fs::dir_ls(here::here(covariate_path, "sample_intensity"), regexp = "-raw.tif$") %>%
  basename() %>%
  str_remove("-raw.tif")

walk(groups, make_tgb)

spec <- "Turdus migratorius"

run_sdm <- function(
  spec,
  covariates,
  psabs_tib,
  n_presence = 5000
) {
  taxa_info <- taxa_join %>%
    filter(species == spec)

  si_group <- taxa_info %>%
    pull(si_group)

  if (!(si_group %in% names(sample_intensity))) {
    message(
      "no sample intensity group generated, can't run some models, skipping"
    )
    return()
  }

  n_filter <- 1
  if (si_group == "birds") {
    n_filter <- 2
  }

  spec_df <- cc_index %>%
    filter(species == spec) %>%
    pull(file) %>%
    open_dataset() %>%
    filter(species == spec) %>%
    filter(n >= n_filter) %>%
    select(species, cell, n) %>%
    collect() %>%
    # double max presence to account for nas
    slice_sample(n = min(n_presence * 2, nrow(.)))


  if (nrow(spec_df) == 0) {
    message(glue::glue("Skipping {spec}: no cell records."))
    return(NULL)
  }

  coords <- xyFromCell(covariates, spec_df$cell)

  spec_psab <- spec_df %>%
    mutate(
      x = coords[, 1],
      y = coords[, 2],
      covariates[cell]
    ) %>%
    # drop nas, resample
    drop_na() %>% 
    slice_sample(n = min(n_presence, nrow(spec_df))) %>%
    bind_rows(psabs_tib) %>%
    ungroup()

  select_vars <- c("pa", base_covariate_names, si_group, "all_species")

  for_model <- spec_psab %>%
    mutate(pa = as.numeric(n > 0)) %>%
    select(names(covariates), pa) %>%
    drop_na() %>%
    select(all_of(select_vars))

  if (n_distinct(for_model$pa) < 2) {
    message(glue::glue(
      "Skipping {spec}: presence/absence not both present after NA drop."
    ))
    return(NULL)
  }
  # # testing glm
  model <- glm(pa ~ ., data = for_model, family = binomial())

  # # maxent

  # model <- maxnet(
  #   p = for_model$pa,
  #   data = for_model %>% select(-pa)
  # )

  # # brt
  # train_params <- training_params(
  #   num_trees = 2000,
  #   shrinkage = 0.01,
  #   interaction_depth = 3,
  #   bag_fraction = 0.5,
  #   num_train = round(0.8 * nrow(for_model)), # or nrow() for no held-out split
  #   min_num_obs_in_node = 10
  # )

  # model <- gbmt(
  #   pa ~ .,
  #   data = for_model,
  #   distribution = gbm_dist("Bernoulli"),
  #   train_params = train_params,
  #   cv_folds = 5,
  #   is_verbose = FALSE
  # )

  # best_iter <- gbmt_performance(model, method = "cv")

  pred_covariates <- covariates
  pred_covariates[["sample_intensity"]] <- max(
    pa$sample_intensity,
    na.rm = TRUE
  )

  # prediction <- predict(
  #   pred_covariates,
  #   model,
  #   n.trees = best_iter,
  #   type = "response"
  # )

  # for maxent
  # prediction <- predict(pred_covariates, model, type = "cloglog", na.rm = T)
  # for uhhhhh glm
  prediction <- predict(pred_covariates, model, type = "response")

  list(species = spec, model = model, prediction = prediction, pa = pa)
}


# ---------------------------------------------------------------------------
# Example: single species (matches original script's demo)
# ---------------------------------------------------------------------------

spec <- "Phrynosoma douglasii"
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

suit_vect <- result$pa %>%
  mutate(result$prediction[cell]) %>%
  vect(geom = c("x", "y"), crs = result$prediction)


ggplot() +
  geom_spatraster(data = result$prediction) +
  geom_spatvector(
    data = suit_vect %>%
      filter(n > 0),
    aes(fill = lyr1),
    shape = 21
  ) +
  scale_fill_viridis_c(
    na.value = NA,
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("Low", "High")
  ) +

  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.2),
    panel.grid = element_line(colour = "grey70"),
    panel.background = element_rect(fill = "grey90"),
    legend.background = element_rect(fill = "grey90", colour = "grey70")
  ) +
  labs(title = result$species, fill = "Suitability")

ggsave("wolv_suit.png", height = 5, width = 5, dpi = 600)

ambi <- 2 * abs(result$prediction - 0.5)

ambi_vect <- result$pa %>%
  mutate(ambi[cell]) %>%
  vect(geom = c("x", "y"), crs = result$prediction)

ggplot() +
  geom_spatraster(data = ambi) +
  geom_spatvector(
    data = ambi_vect %>%
      filter(n > 0),
    aes(fill = lyr1),
    shape = 21
  ) +
  scale_fill_viridis_c(
    na.value = NA,
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("High", "Low")
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.2),
    panel.grid = element_line(colour = "grey70"),
    panel.background = element_rect(fill = "grey90"),
    legend.background = element_rect(fill = "grey90", colour = "grey70")
  ) +
  labs(title = "Gulo gulo", fill = "Ambiguity")

ggsave("wolv_amb.png", height = 5, width = 5, dpi = 600)


bind_rows(
  ambi_vect %>% as_tibble() %>% mutate(name = "ambiguity"),
  suit_vect %>% as_tibble() %>% mutate(name = "suitability")
) %>%
  ggplot() +
  geom_density(aes(x = lyr1, fill = species), alpha = 0.5) +
  facet_wrap(~ fct_rev(name), ncol = 1) +
  theme(legend.position = "bottom") +
  labs(x = NULL, fill = NULL)

ggsave("suit-ambi.png", height = 6, width = 4, dpi = 600)
