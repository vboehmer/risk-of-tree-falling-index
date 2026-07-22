library(sf)
library(dplyr)
library(purrr)
library(future)

#Functions ---------------------------------------------------------
# Process R34_count, R34_USA_SSHS_sum, R64_count, R64_USA_SSHS_sum in chunks
process_chunks <- function(baseline_sf, stormtracks_sf, chunk_size = 10, checkpoint_file = "checkpoint_34_64.rds") {
  
  # Determine starting chunk
  if(file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    start_chunk <- checkpoint$last_chunk + 1
    results <- checkpoint$results
    message("Resuming from chunk ", start_chunk)
  } else {
    start_chunk <- 1
    results <- list()
  }
  
  n_chunks <- ceiling(nrow(baseline_sf) / chunk_size)
  
  for(i in start_chunk:n_chunks) {
    message("Processing chunk ", i, " of ", n_chunks)
    
    idx <- ((i-1)*chunk_size + 1) : min(i*chunk_size, nrow(baseline_sf))
    chunk <- baseline_sf[idx, ]
    
    #### R34 ####
    # First, union R34 polygons by SID so each storm has one merged buffer
    storm_r34 <- stormtracks_sf |>
      group_by(SID) |>
      summarise(R34_geom = st_union(R34_geom)) |>
      st_as_sf()
    
    # Compute intersections between baseline segments and R34 storm polygons
    # Each row in the result corresponds to one road–storm intersection
    intersections_r34 <- st_intersects(chunk, storm_r34)
    
    # Count how many unique storm SIDs intersect each road segment
    chunk$R34_count <- lengths(intersections_r34)
    
    # We'll only keep relevant columns from stormtracks_swio
    stormtracks_r34 <- stormtracks_sf |>
      select(SID, USA_SSHS, R34_geom) |>
      st_set_geometry("R34_geom")
    
    
    # For R34_USA_SSHS_sum
    inter_r34_sf <- st_intersection(chunk %>% select(id), stormtracks_r34) |> select(-geometry)
    
    # Summarise maximum USA_SSHS per (id, SID) and then sum per id
    tmp_r34 <- inter_r34_sf %>%
      st_drop_geometry() %>%
      select(id, SID, USA_SSHS) %>%
      group_by(id, SID) %>%
      slice_max(order_by = USA_SSHS, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      group_by(id) %>%
      summarise(R34_USA_SSHS_sum = sum(USA_SSHS, na.rm = TRUE)) %>%
      ungroup()
    
    chunk <- chunk %>% left_join(tmp_r34, by = "id")
    
    #### R64 ####
    # First, union R64 polygons by SID so each storm has one merged buffer
    storm_r64 <- stormtracks_sf |>
      group_by(SID) |>
      summarise(R64_geom = st_union(R64_geom)) |>
      st_as_sf()
    
    # Compute intersections between baseline segments and R64 storm polygons
    # Each row in the result corresponds to one road–storm intersection
    intersections_r64 <- st_intersects(chunk, storm_r64)
    
    # Count how many unique storm SIDs intersect each road segment
    chunk$R64_count <- lengths(intersections_r64)
    
    # We'll only keep relevant columns from stormtracks_swio
    stormtracks_r64 <- stormtracks_sf |>
      select(SID, USA_SSHS, R64_geom) |>
      st_set_geometry("R64_geom")
    
    
    # For R64_USA_SSHS_sum
    inter_r64_sf <- st_intersection(chunk %>% select(id), stormtracks_r64) |> select(-geometry)
    
    # Summarise maximum USA_SSHS per (id, SID) and then sum per id
    tmp_r64 <- inter_r64_sf %>%
      st_drop_geometry() %>%
      select(id, SID, USA_SSHS) %>%
      group_by(id, SID) %>%
      slice_max(order_by = USA_SSHS, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      group_by(id) %>%
      summarise(R64_USA_SSHS_sum = sum(USA_SSHS, na.rm = TRUE)) %>%
      ungroup()
    
    chunk <- chunk %>% left_join(tmp_r64, by = "id")
    
    # Store chunk and save checkpoint
    results[[i]] <- chunk
    saveRDS(list(last_chunk = i, results = results), checkpoint_file)
  }
  
  # Combine all chunks
  baseline_sf_processed <- bind_rows(results)
  
  return(baseline_sf_processed)
}
# Compute mean distance to storm tracks in chunks
compute_mean_distance_chunks <- function(baseline_sf, stormtracks_sf, chunk_size = 10, checkpoint_file = "distance_checkpoint.rds") {
  
  # Determine starting chunk
  if(file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    start_chunk <- checkpoint$last_chunk + 1
    results <- checkpoint$results
    message("Resuming from chunk ", start_chunk)
  } else {
    start_chunk <- 1
    results <- list()
  }
  
  n_chunks <- ceiling(nrow(baseline_sf) / chunk_size)
  
  for(i in start_chunk:n_chunks) {
    message("Processing chunk ", i, " of ", n_chunks)
    
    idx <- ((i-1)*chunk_size + 1) : min(i*chunk_size, nrow(baseline_sf))
    chunk <- baseline_sf[idx, ]
    
    #### R34 ####
    # First, union R34 polygons by SID so each storm has one merged buffer
    storm_r34 <- stormtracks_sf |>
      group_by(SID) |>
      summarise(R34_geom = st_union(R34_geom),
                geometry = st_union(geometry)) |>
      st_as_sf()
    
    # Compute intersections between baseline segments and R34 storm polygons
    # Each row in the result corresponds to one road–storm intersection
    intersections_r34 <- st_intersects(chunk, storm_r34$R34_geom)
    
    chunk$R34_mean_distance_m <- sapply(seq_len(nrow(chunk)), function(i) {
      # storm indices intersecting this road segment
      storm_ids <- intersections_r34[[i]]
      
      if (length(storm_ids) == 0) {
        # no intersecting storms
        return(NA_real_)
      }
      
      # compute mean distance between this segment and all intersecting storms
      distances <- st_distance(chunk[i, ], storm_r34$geometry[storm_ids, ])
      mean(as.numeric(distances))
    })
    
    # Store chunk and save checkpoint
    results[[i]] <- chunk
    saveRDS(list(last_chunk = i, results = results), checkpoint_file)
  }
  
  # Combine all chunks
  baseline_sf_processed <- bind_rows(results)
  
  return(baseline_sf_processed)
}

# Min-max normalization helper
min_max_norm <- function(x) {
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    return(rep(0, length(x))) # avoid divide-by-zero
  }
  (x - rng[1]) / (rng[2] - rng[1])
}
# Load data --------------------------------------------------------------
stormtracks <- st_read("data/IBTrACS.since1980.list.v04r01.lines.shp")
baseline_sf <- st_read("data/mdg_baseline_roads.gpkg")

# Data Preparation ---------------------------------------------------
stormtracks_utm <- st_transform(stormtracks, st_crs(baseline_sf))

# Filter storms with wind >= 1
stormtracks_subset <- stormtracks_utm %>%
  filter(USA_SSHS >= 1) %>%
  rowwise() %>%
  mutate(
    R34_buffer = mean(c(USA_R34_NE, USA_R34_NW, USA_R34_SE, USA_R34_SW), na.rm = TRUE),
    R64_buffer = mean(c(USA_R64_NE, USA_R64_NW, USA_R64_SE, USA_R64_SW), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(USA_SSHS) %>%
  mutate(
    R34_buffer = ifelse(is.na(R34_buffer), mean(R34_buffer, na.rm = TRUE), R34_buffer),
    R64_buffer = ifelse(is.na(R64_buffer), mean(R64_buffer, na.rm = TRUE), R64_buffer)
  ) %>%
  ungroup() %>%
  mutate(
    R34_buffer_km = R34_buffer * 1.852,
    R64_buffer_km = R64_buffer * 1.852
  )

# Keep only valid geometries
geom_points <- sapply(st_geometry(stormtracks_subset), function(g) nrow(st_coordinates(g)))
stormtracks_clean <- stormtracks_subset[geom_points > 1, ] %>%
  filter(st_geometry_type(geometry) %in% c("LINESTRING", "MULTILINESTRING"))

# SWIO bounding box
swio_bbox <- st_bbox(c(xmin = 30, ymin = -40, xmax = 90, ymax = 10), crs = st_crs(4326))
swio_poly_utm <- st_transform(st_as_sfc(swio_bbox), st_crs(stormtracks_utm))

stormtracks_swio <- stormtracks_clean %>%
  st_filter(swio_poly_utm) %>%
  select(
    SID, SEASON, USA_SSHS,
    USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW,
    USA_R64_NE, USA_R64_SE, USA_R64_SW, USA_R64_NW,
    R34_buffer, R34_buffer_km, R64_buffer, R64_buffer_km,
    geometry
  ) %>%
  rowwise() %>%
  mutate(
    R34_geom = st_buffer(geometry, R34_buffer_km * 1000),
    R64_geom = st_buffer(geometry, R64_buffer_km * 1000)
  ) %>%
  ungroup()

# Processing -------------------------------------------------------
# Run chunk processing
baseline_sf <- process_chunks(baseline_sf, stormtracks_swio, chunk_size = 1000)
# Compute mean distances
baseline_sf <- compute_mean_distance_chunks(baseline_sf, stormtracks_swio, chunk_size = 1000)

# Follow-up Work --------------------------------------------------------
# Replace NAs with 0 
baseline_sf[] <- lapply(baseline_sf, function(x) {
  if (is.numeric(x)) x[is.na(x)] <- 0
  x
})

# Columns to normalize
cols_to_normalize <- c("R34_count", "R34_USA_SSHS_sum", 
                       "R64_count", "R64_USA_SSHS_sum", 
                       "R34_mean_distance_m")

# Apply function
baseline_sf <- baseline_sf %>%
  mutate(across(all_of(cols_to_normalize), min_max_norm, .names = "{.col}_norm"))

# Final calculation of risk --------------------------------------
# invert distance to get high values -> more risk
baseline_sf$distance_risk <- 1 - baseline_sf$R34_mean_distance_m_norm

# Compute cyclone_risk_score
baseline_sf <- baseline_sf %>%
  mutate(cyclone_risk_score = (
    0.2*R34_count_norm +
      0.2*R34_USA_SSHS_sum_norm +
      0.25*R64_count_norm +
      0.25*R64_USA_SSHS_sum_norm +
      0.1*(1 - R34_mean_distance_m_norm)
  ))

# Save final output -------------------------------------------------------
st_write(baseline_sf, "data/mdg_roads_cyclone_risk.gpkg", delete_dsn = TRUE)
