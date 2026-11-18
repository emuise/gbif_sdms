library(arrow)
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(bcmaps)
muiselib2::muiselib_setup()

# functions to get grouped intensities (uses cc, snap, snap_vals at call time)
source(here::here("scripts", "03a_calc_group_intensity_functions.R"))

# set paths and mkdir

covariate_path <- here::here("data", "covariates")
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
bioclim_fold <- here::here(covariate_path, "Normal_1991_2020_bioclim")

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
      !is.na(species),
      (class == "Aves" & year >= 1990) |
        ((is.na(class) | class != "Aves") & year >= 1970)
    ) %>%
    select(
      kingdom, phylum, order, class, genus, species,
      decimallatitude, decimallongitude
    ) %>%
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

  df %>%
    filter(!is.na(cell)) %>% # removes anything outside the raster bounds
    mutate(snap_val = snap_vals[cell]) %>%
    filter(!is.na(snap_val)) %>% # only keeps cells with climate values
    select(-snap_val) %>%
    count(kingdom, phylum, order, class, genus, species, cell, name = "n") %>%
    arrow::write_parquet(savename)
})
# 35 and 1494 evidently have nothing in them
# 35 is all birds prior to 1982
# 1494 is some plants prior to 1940


# per-phylum sorted files for fast SDM lookups

out_dir <- here::here("data", "cellcounts_by_phylum")
fs::dir_create(out_dir)

cc <- open_dataset(cellcount_dir)

groups <- cc %>%
  distinct(kingdom, phylum) %>%
  collect()

safe <- \(x) {
  if (is.na(x)) "unknown" else str_replace_all(x, "[^A-Za-z0-9]", "_")
}

filter_rank <- function(q, col, val) {
  if (is.na(val)) {
    q %>% filter(is.na({{ col }}))
  } else {
    q %>% filter({{ col }} == val)
  }
}

walk(seq_len(nrow(groups)), \(i) {
  k <- groups$kingdom[i]
  p <- groups$phylum[i]

  cc %>%
    filter_rank(kingdom, k) %>%
    filter_rank(phylum, p) %>%
    group_by(kingdom, phylum, species, cell) %>% # drops class/order/genus
    summarise(n = sum(n)) %>%
    collect() %>%
    arrange(species, cell) %>%
    write_parquet(
      here::here(out_dir, glue::glue("{safe(k)}__{safe(p)}.parquet")),
      chunk_size = 50000
    )
})

species_index <- cc %>%
  distinct(kingdom, phylum, species) %>%
  collect() %>%
  mutate(
    file = here::here(
      out_dir,
      glue::glue("{map_chr(kingdom, safe)}__{map_chr(phylum, safe)}.parquet")
    )
  )

write_parquet(
  species_index,
  here::here("data", "cellcounts_species_index.parquet")
)


# sampling intensity surfaces for taxonomic groups


taxa <- cc %>%
  distinct(kingdom, class, order, species) %>%
  collect() %>%
  filter(!is.na(kingdom))

iter_groups <- list(
  all_species         = taxa,
  mammals             = filter(taxa, class == "Mammalia"),
  birds               = filter(taxa, class == "Aves"),
  plants              = filter(taxa, kingdom == "Plantae"),
  reptiles_amphibians = filter(taxa, class %in% c("Reptilia", "Amphibia")),
  butterflies         = filter(taxa, order == "Lepidoptera")
) %>%
  map(\(x) distinct(x, kingdom, species))

imap(iter_groups, \(spp, name) {
  calc_smooth_intensity(spp, name)
})