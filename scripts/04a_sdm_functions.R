library(tidymodels)
library(tidysdm) # maxent() parsnip model (maxnet engine), tss_max()
library(ranger)
library(xgboost)
library(maxnet)

# sdm fitting functions, sourced in 04
#
# methods follow:
# - eckert et al. 2023 (nat comms, 30x30, pollock lab): brt + maxent,
#   5000 presences / 10000 background, sample effort covariate set to max at
#   prediction, maxent hinge-only, brt tc = 4 / bag = 0.6
# - wightman et al. 2025 (biorxiv, pollock lab): bias covariate correction
#   (30 km gbif density, fixed at max for prediction), brt lr = 0.08,
#   70/30 train/test, auroc + auprc + tss
# - valavi et al. 2022 (ecol monogr), cited by both: down-sampled random forest
#   (equal presence/background per tree, replace = TRUE, 1000 trees)
#
# for_model is the tibble built in run_sdm(): pa (0/1) + covariate columns,
# including the bias covariates (si_group column and all_species)

# ---------------------------------------------------------------------------
# data prep
# ---------------------------------------------------------------------------

# tidysdm/yardstick convention: factor outcome, presence as the first level
prep_for_model <- function(for_model) {
  for_model %>%
    mutate(
      pa = if_else(pa == 1, "presence", "pseudoabs") %>%
        factor(levels = c("presence", "pseudoabs"))
    )
}

sdm_recipe <- function(for_model) {
  recipe(pa ~ ., data = for_model) %>%
    step_zv(all_predictors()) %>%
    step_dummy(all_nominal_predictors()) # xgboost needs numeric predictors
}

# ---------------------------------------------------------------------------
# model specs
# ---------------------------------------------------------------------------

# down-sampled rf (valavi et al. 2022): each tree sees n_pres presences and
# n_pres background points, sampled with replacement
# ranger class-wise sample.fraction is ordered by factor level (presence first)
rf_spec <- function(for_model, trees = 1000, n_threads = 1) {
  n_pres <- sum(for_model$pa == "presence")
  n_tot <- nrow(for_model)

  rand_forest(trees = trees) %>%
    set_engine(
      "ranger",
      sample.fraction = c(n_pres, n_pres) / n_tot,
      replace = TRUE,
      num.threads = n_threads
    ) %>%
    set_mode("classification")
}

# brt: tree complexity 4, bag fraction 0.6 (eckert, wightman)
# lr 0.08 from wightman; eckert used 0.001 (needs far more trees)
# gbm.step picks n trees by cv; here xgboost early-stops on a 20% holdout
brt_spec <- function(
  learn_rate = 0.08,
  tree_depth = 4,
  bag_fraction = 0.6,
  max_trees = 5000,
  stop_iter = 50,
  n_threads = 1
) {
  boost_tree(
    trees = max_trees,
    tree_depth = tree_depth,
    learn_rate = learn_rate,
    sample_size = bag_fraction,
    stop_iter = stop_iter
  ) %>%
    set_engine("xgboost", validation = 0.2, nthread = n_threads) %>%
    set_mode("classification")
}

# maxent: hinge features only (eckert, wightman), default regularization
maxent_spec <- function(feature_classes = "h", regularization_multiplier = 1) {
  maxent(
    feature_classes = feature_classes,
    regularization_multiplier = regularization_multiplier
  ) %>%
    set_engine("maxnet")
}

sdm_specs <- function(for_model, n_threads = 1) {
  list(
    rf = rf_spec(for_model, n_threads = n_threads),
    brt = brt_spec(n_threads = n_threads),
    maxent = maxent_spec()
  )
}

# ---------------------------------------------------------------------------
# fit + evaluate
# ---------------------------------------------------------------------------

sdm_metrics <- metric_set(roc_auc, pr_auc, tss_max)

fit_one <- function(spec, for_model) {
  workflow() %>%
    add_recipe(sdm_recipe(for_model)) %>%
    add_model(spec) %>%
    fit(data = for_model)
}

# 70/30 split (wightman) stratified on pa, evaluate on the test set,
# then refit on all data for prediction
# returns a tibble: one row per model with metrics and the full-data fit
fit_sdms <- function(
  for_model,
  prop = 0.7,
  refit_full = TRUE,
  n_threads = 1,
  seed = 42
) {
  set.seed(seed)

  for_model <- prep_for_model(for_model)

  split <- initial_split(for_model, prop = prop, strata = pa)
  train <- training(split)
  test <- testing(split)

  # rf sample.fraction depends on the data it is fit to, so build specs per set
  train_specs <- sdm_specs(train, n_threads)

  imap(train_specs, \(spec, name) {
    message(glue::glue("fitting {name}"))

    train_fit <- fit_one(spec, train)

    metrics <- train_fit %>%
      predict(test, type = "prob") %>%
      bind_cols(test %>% select(pa)) %>%
      sdm_metrics(truth = pa, .pred_presence)

    full_fit <- train_fit
    if (refit_full) {
      full_fit <- fit_one(sdm_specs(for_model, n_threads)[[name]], for_model)
    }

    tibble(
      model = name,
      metrics = list(metrics),
      fit = list(full_fit)
    )
  }) %>%
    bind_rows()
}

# ---------------------------------------------------------------------------
# prediction
# ---------------------------------------------------------------------------

# bias covariates are fixed at their max in the training data everywhere
# (eckert, wightman), so predictions reflect "fully sampled" conditions
bias_max_values <- function(for_model, bias_vars) {
  for_model %>%
    summarise(across(all_of(bias_vars), \(x) max(x, na.rm = TRUE))) %>%
    as.list()
}

predict_sdm <- function(fit, covariates, bias_values) {
  pred_vars <- fit %>%
    extract_preprocessor() %>%
    pluck("var_info") %>%
    filter(role == "predictor") %>%
    pull(variable)

  terra::predict(
    covariates[[pred_vars]],
    fit,
    fun = \(model, data, ...) {
      for (v in names(bias_values)) {
        data[[v]] <- bias_values[[v]]
      }
      predict(model, data, type = "prob")$.pred_presence
    },
    na.rm = TRUE
  )
}

# predicts every model in a fit_sdms() result, returns a named multi-layer raster
predict_sdms <- function(sdm_fits, covariates, for_model, bias_vars) {
  bias_values <- bias_max_values(for_model, bias_vars)

  preds <- map(sdm_fits$fit, \(f) predict_sdm(f, covariates, bias_values)) %>%
    terra::rast()

  names(preds) <- sdm_fits$model
  preds
}

# weighted mean across models
# eckert et al. used brt 0.7 / maxent 0.3; default here is equal weights
ensemble_sdm <- function(preds, weights = c(rf = 1, brt = 1, maxent = 1)) {
  weights <- weights[names(preds)]
  weights <- weights / sum(weights)

  terra::weighted.mean(preds, w = weights)
}
