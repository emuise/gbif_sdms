library(arrow)
library(tidyverse)
library(terra)
library(sf)
library(bcmaps)

# set paths and make them

clean_dir <- here::here("data", "gbif_noram")
spatcount_dir <- here::here("data", "spatial_counts")
fs::dir_create(spatcount_dir)

# operate on condensed files (fast to list since already batched in script 1)
c_files <- fs::dir_ls(clean_dir, recurse = TRUE, type = "file")


bcb <- bcmaps::bc_bound_hres() %>%
  vect()

bcb_sf <- bcb %>%
  st_as_sf()

bcb_bbox <- bcb %>%
  project("epsg:4326") %>%
  ext()

# ---------------------------------------------------------------------------
# Per-file in-BC / not-in-BC species counts
# ---------------------------------------------------------------------------
# birds are going to be unable to process all at once
# i.e. geese have > 100 million observations
# every other family is /reasonably large/ and can be done at once
# so this is done on a per-file basis, condensed at the end

walk(
  seq_along(c_files),
  \(n) {
    x <- c_files[n]
    savename <- here::here(
      spatcount_dir,
      glue::glue("{basename(x)}")
    )
    if (file.exists(savename)) {
      return()
    }

    message(glue::glue("[{n}/{length(c_files)}]: {basename(x)}"))

    df <- arrow::open_dataset(x) %>%
      select(
        kingdom,
        phylum,
        class,
        order,
        family,
        genus,
        species,
        decimallatitude,
        decimallongitude
      ) %>%
      collect()

    # index of if points are in bounding box
    in_bbox <- df$decimallatitude <= bcb_bbox$ymax &
      df$decimallatitude >= bcb_bbox$ymin &
      df$decimallongitude <= bcb_bbox$xmax &
      df$decimallongitude >= bcb_bbox$xmin

    # static integer column
    df$in_bc <- 0L

    # if any are in the bounding box, do a spatial intersect to see if they are in the BC polygon
    if (any(in_bbox)) {
      # sf is fastest for these intersections, keeping the polygon in its normal crs is also fastest
      # so we project the points from wgs84 into bc albers for speed
      # we don't need to reproject out or anything
      spat <- vect(
        df[in_bbox, ],
        geom = c("decimallongitude", "decimallatitude"),
        crs = "epsg:4326"
      ) %>%
        project(bcb) %>%
        st_as_sf()

      inters <- st_intersects(spat, bcb_sf, sparse = TRUE)

      # if number of intersects greater than 0, it is within the bounding box of bc
      df$in_bc[in_bbox] <- as.integer(lengths(inters) > 0)
    }

    # group by summarize if it is inside or outside bc
    inout <- df %>%
      group_by(kingdom, phylum, class, order, family, genus, species) %>%
      summarize(
        n_bc = sum(in_bc == 1L),
        n_notbc = sum(in_bc == 0L),
        n_total = n(),
        .groups = "drop"
      ) %>%
      mutate(file = x)

    arrow::write_parquet(inout, savename)
  }
)

# ---------------------------------------------------------------------------
# Aggregate across files, keep species with at least one BC record
# ---------------------------------------------------------------------------

counts_wfiles <- arrow::open_dataset(spatcount_dir)

all_counts <- counts_wfiles %>%
  group_by(kingdom, phylum, class, order, family, genus, species) %>%
  summarize(
    n_bc = sum(n_bc),
    n_notbc = sum(n_notbc),
    n_total = sum(n_total)
  )

bc_counts <- all_counts %>%
  filter(n_bc > 0) %>%
  collect() %>%
  arrange(kingdom, phylum, class, order, family, genus, species)

species_list <- bc_counts %>%
  pull(species)

out_dir <- here::here("data", "bc_species")
fs::dir_create(out_dir)

arrow::write_parquet(bc_counts, here::here(out_dir, "bc_counts.parquet"))
saveRDS(species_list, here::here(out_dir, "species_list.rds"))

message(
  glue::glue(
    "Done: {length(species_list)} species with at least one BC record."
  )
)
