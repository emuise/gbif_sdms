# DEFUNCT SAMPLING INTESNITY MAP CODE


ggplot() +
  geom_spatraster(data = l1p_smooth_r) +
  scale_fill_viridis_c(
    na.value = NA,
    limits = c(0, minmax(l1p_smooth_r)[2]),
    breaks = c(0, minmax(l1p_smooth_r)[2]),
    labels = c("Low", "High")
  ) +
  labs(
    title = "log(mean sampling intensity + 1) (30 km focal window)",
    fill = NULL
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.2),
    panel.grid = element_line(colour = "grey70"),
    panel.background = element_rect(fill = "grey90"),
    legend.background = element_rect(fill = "grey90", colour = "grey70")
  )

ggsave("sampling_intensity.png")

bc_smooth <- smooth_r %>%
  crop(
    bcb %>%
      project(smooth_r),
    mask = TRUE
  ) %>%
  project(bcb)

plot(log(bc_smooth + 1))

# species accumulation curves
bc_cells <- rasterize(bcb %>% project(smooth_r), smooth_r) %>%
  cells()

species_accumulation_curve <- function(spec_list) {
  subset_counts <- cellcounts %>%
    filter(species %in% spec_list) %>%
    filter(cell %in% bc_cells) %>%
    collect() %>%
    ungroup()

  # --- Chao2 ---
  incidence <- subset_counts %>%
    distinct(species, cell) %>%
    count(species, name = "n_i")

  Q1 <- sum(incidence$n_i == 1)
  Q2 <- sum(incidence$n_i == 2)
  S_obs <- nrow(incidence)

  chao2 <- S_obs + (Q1 * (Q1 - 1)) / (2 * (Q2 + 1))

  # --- accumulation curve (unchanged) ---
  iters <- map_dfr(
    1:100,
    \(x) {
      slice_sample(subset_counts, prop = 1) %>%
        mutate(
          is_new = !duplicated(species),
          cum_species = cumsum(is_new),
          cum_records = row_number(),
          cum_n = cumsum(n)
        ) %>%
        select(starts_with("cum")) %>%
        mutate(x)
    },
    .progress = T
  )

  smry_iters <- iters %>%
    group_by(cum_records) %>%
    summarize(
      mean_richness = mean(cum_species),
      sd_richness = sd(cum_species),
      lower_ci = mean_richness - 1.96 * sd_richness,
      upper_ci = mean_richness + 1.96 * sd_richness
    )

  p <- ggplot(smry_iters, aes(x = cum_records, y = mean_richness)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.5) +
    geom_line() +
    geom_hline(yintercept = chao2, linetype = "dashed", color = "firebrick") +
    geom_hline(yintercept = S_obs, linetype = "dotted", color = "grey40") +
    annotate(
      "text",
      x = max(smry_iters$cum_records),
      y = chao2,
      label = glue::glue("Chao2 \u2248 {round(chao2, 1)}"),
      hjust = 1,
      vjust = -0.5,
      color = "firebrick",
      size = 3
    ) +
    labs(
      x = "Cumulative records",
      y = "Species richness",
      caption = glue::glue(
        "S_obs = {S_obs}  |  Chao2 = {round(chao2, 1)}  |  ",
        "Est. undetected = {round(chao2 - S_obs, 1)}  |  ",
        "Q1 = {Q1}, Q2 = {Q2}"
      )
    )

  p

  list(
    plot = p,
    summary = smry_iters,
    chao2 = chao2,
    S_obs = S_obs,
    Q1 = Q1,
    Q2 = Q2,
    incidence = incidence
  )
}

egg <- species_accumulation_curve(
  taxa_join %>% filter(class == "Aves") %>% pull(species)
)
