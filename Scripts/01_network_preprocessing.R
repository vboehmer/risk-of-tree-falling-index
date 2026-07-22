# ===========================================================================
# 01_network_preprocessing.R
# Downloads OSM road data for Madagascar via ohsome API, builds a
# cleaned sfnetwork, and exports baseline road segments (100 m each).
#
# Input:  data/20240108_MAD_CRM_Warehouses_updated.gpkg
#         data/mdg_admin2.shp
# Output: data/cleaned_mdg_roads.gpkg
#         data/mdg_baseline_roads.gpkg
# ===========================================================================

library(ohsome)
library(dplyr)
library(sf)
library(tidyverse)
library(sfnetworks)
library(tidygraph)
library(igraph)

# FUNCTIONS -----------------------------------------------------------------

snap_line_point <- function(sf_lines, moving_id, target_id,
                            moving_end = "end", target_end = "end") {
  line_move <- sf_lines |> filter(osmId == moving_id)
  line_target <- sf_lines |> filter(osmId == target_id)

  coords_move <- st_coordinates(line_move)[,1:2]
  coords_target <- st_coordinates(line_target)[,1:2]

  move_idx <- if (moving_end == "end") nrow(coords_move) else 1
  target_idx <- if (target_end == "end") nrow(coords_target) else 1

  coords_move[move_idx, ] <- coords_target[target_idx, ]

  new_line <- st_linestring(as.matrix(coords_move))

  sf_lines <- sf_lines |>
    mutate(
      geometry = case_when(
        osmId == moving_id ~ st_sfc(new_line, crs = st_crs(sf_lines)),
        TRUE ~ geometry
      )
    ) |>
    st_as_sf()

  return(sf_lines)
}

split_line <- function(line, segment_length = 50) {
  total_length <- st_length(line)
  distances <- seq(0, as.numeric(total_length), by = segment_length)

  if (tail(distances, 1) < as.numeric(total_length)) {
    distances <- c(distances, as.numeric(total_length))
  }

  pieces <- st_line_sample(line, sample = distances / as.numeric(total_length))
  coords <- st_coordinates(pieces)

  segments <- list()
  for (i in seq_len(nrow(coords) - 1)) {
    segments[[i]] <- st_linestring(rbind(coords[i, c("X", "Y")], coords[i+1, c("X", "Y")]))
  }

  sfc <- st_sfc(segments, crs = st_crs(line))
  lengths <- sapply(sfc, st_length)

  st_sf(
    id = seq_along(segments),
    length_m = as.numeric(lengths),
    geometry = sfc
  )
}

# GET OSM DATA ---------------------------------------------------------------

warehouses <- st_read("data/20240108_MAD_CRM_Warehouses_updated.gpkg",
                       layer = "20240108_MAD_CRM_Warehouses_updated") |>
  st_transform(29738) |>
  st_make_valid()

mdg_boundary <- st_read("data/mdg_admin2.shp") |> st_union() |> st_make_valid()
bb_mdg <- st_bbox(mdg_boundary)

mdg_roads_raw <- ohsome_elements_geometry(
  boundary = bb_mdg,
  filter = "highway=primary or highway=secondary or highway=tertiary or highway=tertiary_link or route=ferry or bridge=yes or ford=yes",
  time = "2025-08-01",
  properties = "tags",
  clipGeometry = FALSE
) |>
  ohsome_post()

mdg_roads <- mdg_roads_raw |>
  mutate(
    type = case_when(
      highway == "primary"       ~ "primary",
      highway == "secondary"     ~ "secondary",
      highway == "tertiary"      ~ "tertiary",
      highway == "tertiary_link" ~ "tertiary",
      route   == "ferry"         ~ "ferry",
      bridge  == "yes"           ~ "bridge",
      ford    == "yes"           ~ "ford",
      TRUE ~ "other")) |>
  rename(osmId = '@osmId', snapshotTimestamp = '@snapshotTimestamp') |>
  select(osmId, type, snapshotTimestamp)

# NETWORK PREPROCESSING ------------------------------------------------------

mdg_roads <- snap_line_point(mdg_roads,
                             moving_id = "way/765455429",
                             target_id = "way/1395745576",
                             moving_end = "end", target_end = "start")

mdg_roads <- snap_line_point(mdg_roads,
                             moving_id = "way/199885144",
                             target_id = "way/199885154",
                             moving_end = "start", target_end = "end")

mdg_roads <- snap_line_point(mdg_roads,
                             moving_id = "way/685977021",
                             target_id = "way/685977020",
                             moving_end = "end", target_end = "end")

mdg_roads <- snap_line_point(mdg_roads,
                             moving_id = "way/44889718",
                             target_id = "way/53763258",
                             moving_end = "end", target_end = "start")

mdg_roads <- snap_line_point(mdg_roads, "way/593009607", "way/507252620")

mdg_roads <- st_as_sf(mdg_roads)

## Convert to sfnetwork ------------------------------------------------------

weighting_profile <- c(
  primary = 65, secondary = 60, tertiary = 50,
  ferry = 20, bridge = 30, ford = 10)

mdg_roads <- mdg_roads |>
  filter(st_geometry_type(geometry) %in% c("LINESTRING", "MULTILINESTRING")) |>
  st_cast("LINESTRING") |>
  st_make_valid()

mdg_net <- as_sfnetwork(mdg_roads, directed = FALSE) |>
  st_transform(29738) |>
  activate("edges") |>
  mutate(speed_km_h = weighting_profile[type]) |>
  mutate(edge_length = edge_length()) |>
  mutate(speed_m_s = speed_km_h * 1000 / 3600) |>
  mutate(travel_time_s = as.double((edge_length / speed_m_s)))

## Simplify network ----------------------------------------------------------

simple_mdg_net <- mdg_net |>
  activate("edges") |>
  arrange(edge_length()) |>
  filter(!edge_is_multiple())

## Subdivide edges -----------------------------------------------------------

subdivision_mdg_net <- convert(simple_mdg_net, to_spatial_subdivision)

subdivision_mdg_net <- subdivision_mdg_net |>
  activate("edges") |>
  filter(!is.na(travel_time_s))

## Largest connected component -----------------------------------------------

comp <- igraph::components(as_tbl_graph(subdivision_mdg_net))
largest_comp <- which.max(comp$csize)

subdivision_mdg_net <- subdivision_mdg_net |>
  activate("nodes") |>
  mutate(component = comp$membership)

filtered_subdivision_mdg_net <- subdivision_mdg_net |>
  filter(component == largest_comp)

## Export cleaned roads -------------------------------------------------------

cleaned_mdg_roads <- filtered_subdivision_mdg_net |>
  activate("edges") |>
  as_tibble() |>
  st_as_sf() |>
  select(osmId, type, speed_km_h, edge_length, speed_m_s, travel_time_s)

st_write(cleaned_mdg_roads, "data/cleaned_mdg_roads.gpkg",
         layer = "cleaned_mdg_roads", delete_dsn = TRUE)

# BASELINE ROAD DATASET (100 m segments) -------------------------------------

edges <- st_as_sf(cleaned_mdg_roads)
merged <- st_line_merge(st_union(edges$geom))
lines <- st_cast(merged, "LINESTRING")

segments_list <- lapply(seq_along(lines), function(i) {
  segs <- split_line(lines[i], segment_length = 100)
  segs |> mutate(id = paste0(i, "_", id))
})

baseline_sf <- do.call(rbind, segments_list)

st_write(baseline_sf, "data/mdg_baseline_roads.gpkg",
         layer = "mdg_baseline_roads", delete_dsn = TRUE)
