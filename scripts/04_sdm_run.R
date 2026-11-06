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

taxa_join <- here::here("data", "spatial_counts") %>%
  arrow::open_dataset() %>%
  distinct(kingdom, phylum, class, order, family, genus, species) %>%
  collect()

cc <- arrow::open_dataset(cellcount_dir)

cellcounts <- cc %>%
  group_by(species, cell) %>%
  summarise(n = sum(n))

# covariate stac
bc_vars <- c(
  "MAP", # mean annual precipitation (mm)
  "DD_0", # chilling degree days (Degree days below 0 °C)
  "PAS", # precipitation as snow (mm)
  "CMD", # Hargreave's climatic moisture index
  "DD18" # warming degree days above 18 °C
)

bioclim <- here::here(
  bioclim_fold,
  glue::glue("Normal_1991_2020_{bc_vars}.tif")
) %>%
  rast()

names(bioclim) <- bc_vars

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
