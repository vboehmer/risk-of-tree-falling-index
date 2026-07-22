# ===========================================================================
# 02_centrality_analysis.R
# Computes population-weighted betweenness centrality and travel-time
# impact (Delta-T) for Madagascar's road network relative to CRM warehouses.
#
# Requires: Scripts/01_network_preprocessing.R (run first)
# Input:    data/cleaned_mdg_roads.gpkg
#           data/mdg_baseline_roads.gpkg
#           data/20240108_MAD_CRM_Warehouses_updated.gpkg
#           data/mdg_pop_2025_CN_100m_R2025A_v1.tif
# Output:   data/mdg_roads_nearest.gpkg
#           data/mdg_road_centrality.gpkg
# ===========================================================================

library(dplyr)
library(sf)
library(tidyverse)
library(sfnetworks)
library(tidygraph)
library(igraph)
library(terra)
library(purrr)
library(parallel)

# FUNCTIONS -----------------------------------------------------------------

min_max_norm <- function(x) {
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

compute_edge_delta <- function(eid) {
  g_mod <- delete_edges(g, eid)
  new_dist <- distances(g_mod, v = warehouse_nodes, to = demand_nodes_unique,
                        weights = E(g_mod)$weight)

  new_mean <- sapply(seq_along(demand_nodes_unique), function(i) {
    mean(new_dist[nearest_3[, i], i], na.rm = TRUE)
  })

  new_mean[!is.finite(new_mean)] <- finite_max

  mean((new_mean[demand_to_unique_idx] - base_min_all) * pop_weights,
       na.rm = TRUE) / mean(pop_weights)
}

# LOAD PROCESSED DATA --------------------------------------------------------

cleaned_mdg_roads <- st_read("data/cleaned_mdg_roads.gpkg",
                             layer = "cleaned_mdg_roads")

warehouses <- st_read("data/20240108_MAD_CRM_Warehouses_updated.gpkg",
                       layer = "20240108_MAD_CRM_Warehouses_updated") |>
  st_transform(29738) |>
  st_make_valid()

baseline_sf <- st_read("data/mdg_baseline_roads.gpkg",
                        layer = "mdg_baseline_roads")

filtered_subdivision_mdg_net <- as_sfnetwork(cleaned_mdg_roads, directed = FALSE)

# CLOSEST WAREHOUSE ANALYSIS ------------------------------------------------

warehouses <- st_transform(warehouses, st_crs(baseline_sf))

nodes_sf <- filtered_subdivision_mdg_net |> activate("nodes") |> st_as_sf()
warehouses_nodes <- st_nearest_feature(warehouses, nodes_sf)

net_igraph <- as.igraph(filtered_subdivision_mdg_net, directed = FALSE)

dist_matrix <- igraph::distances(
  net_igraph,
  v = warehouses_nodes,
  to = V(net_igraph),
  weights = E(net_igraph)$travel_time_s
)

edges_sf <- filtered_subdivision_mdg_net |> activate("edges") |> st_as_sf()
from_nodes <- edges_sf$from

nearest_warehouse <- apply(dist_matrix[, from_nodes, drop = FALSE], 2, which.min)
edges_sf <- edges_sf |> mutate(nearest_warehouse_index = nearest_warehouse)

st_write(edges_sf, "data/mdg_roads_nearest.gpkg", delete_dsn = TRUE)

# BETWEENNESS CENTRALITY -----------------------------------------------------

net <- filtered_subdivision_mdg_net

ghspop <- rast('data/mdg_pop_2025_CN_100m_R2025A_v1.tif')
ghspop_proj <- project(ghspop, st_crs(net)$wkt)

hull <- net |>
  activate("nodes") |>
  st_geometry() |>
  st_combine() |>
  st_convex_hull()

ghspop_crop <- crop(ghspop_proj, hull, mask = TRUE)

set.seed(128)
pop_points <- as.points(ghspop_crop, values = TRUE)
n_samples <- 1000

cell_ids <- sample(seq_len(nrow(pop_points)), size = n_samples,
                   prob = pop_points$mdg_pop_2025_CN_100m_R2025A_v1,
                   replace = TRUE)

demand_sites <- pop_points[cell_ids, ] |> st_as_sf()
demand_sites <- demand_sites |> rename(pop_value = mdg_pop_2025_CN_100m_R2025A_v1)
demand_nodes <- st_nearest_feature(demand_sites, st_as_sf(net, "nodes"))

g <- as.igraph(net)
warehouse_nodes <- st_nearest_feature(warehouses, st_as_sf(net, "nodes"))

## 1 nearest warehouse -------------------------------------------------------

edge_btwn_scores <- numeric(ecount(g))

for (dest_id in demand_nodes) {
  dists <- distances(g, v = dest_id, to = warehouse_nodes, weights = E(g)$travel_time_s)
  nearest_wh <- warehouse_nodes[which.min(dists)]

  sp <- shortest_paths(g, from = nearest_wh, to = dest_id,
                       weights = E(g)$travel_time_s, output = "epath")$epath[[1]]

  if (length(sp) > 0) {
    edge_btwn_scores[sp] <- edge_btwn_scores[sp] + 1
  }
}

net <- net |> activate("edges") |> mutate(betweenness_wh = edge_btwn_scores)

## 3 nearest warehouses ------------------------------------------------------

edge_btwn_scores <- numeric(ecount(g))

for (dest_id in demand_nodes) {
  dists <- distances(g, v = dest_id, to = warehouse_nodes, weights = E(g)$travel_time_s)
  nearest_whs <- warehouse_nodes[order(dists)[1:3]]

  for (wh in nearest_whs) {
    sp <- shortest_paths(g, from = wh, to = dest_id,
                         weights = E(g)$travel_time_s, output = "epath")$epath[[1]]
    if (length(sp) > 0) {
      edge_btwn_scores[sp] <- edge_btwn_scores[sp] + 1
    }
  }
}

net <- net |> activate("edges") |> mutate(betweenness_wh3 = edge_btwn_scores)

## Population-weighted (3 nearest warehouses) ---------------------------------

edge_btwn_scores <- numeric(ecount(g))

for (i in seq_along(demand_nodes)) {
  dest_id <- demand_nodes[i]
  demand_pop <- demand_sites$pop_value[i]
  if (is.na(demand_pop) || demand_pop <= 0) next

  dists <- distances(g, v = dest_id, to = warehouse_nodes, weights = E(g)$travel_time_s)
  nearest_whs <- warehouse_nodes[order(dists)[1:3]]

  for (wh in nearest_whs) {
    sp <- shortest_paths(g, from = wh, to = dest_id,
                         weights = E(g)$travel_time_s, output = "epath")$epath[[1]]
    if (length(sp) > 0) {
      edge_btwn_scores[sp] <- edge_btwn_scores[sp] + demand_pop
    }
  }
}

net <- net |> activate("edges") |> mutate(betweenness_popweighted = edge_btwn_scores)

## Population-weighted (all warehouses) --------------------------------------

edge_btwn_scores <- numeric(ecount(g))

for (i in seq_along(demand_nodes)) {
  dest_id <- demand_nodes[i]
  demand_pop <- demand_sites$pop_value[i]
  if (is.na(demand_pop) || demand_pop <= 0) next

  for (wh in warehouse_nodes) {
    sp <- shortest_paths(g, from = wh, to = dest_id,
                         weights = E(g)$travel_time_s, output = "epath")$epath[[1]]
    if (length(sp) > 0) {
      edge_btwn_scores[sp] <- edge_btwn_scores[sp] + demand_pop
    }
  }
}

net <- net |> activate("edges") |> mutate(betweenness_popweighted_allwh = edge_btwn_scores)

# TRAVEL-TIME IMPACT (Delta-T) -----------------------------------------------

g <- as.igraph(net)

edges_sf <- filtered_subdivision_mdg_net |> activate("edges") |> st_as_sf()

E(g)$weight <- edges_sf$travel_time_s
E(g)$edge_id <- seq_len(ecount(g))

warehouse_nodes <- st_nearest_feature(warehouses, st_as_sf(net, "nodes"))
demand_nodes_unique <- unique(demand_nodes)
demand_to_unique_idx <- match(demand_nodes, demand_nodes_unique)
pop_weights <- replace(demand_sites$pop_value, is.na(demand_sites$pop_value), 0)

base_dist <- distances(g, v = warehouse_nodes, to = demand_nodes_unique,
                       weights = E(g)$weight)
finite_max <- max(base_dist[is.finite(base_dist)], na.rm = TRUE)

nearest_3 <- apply(base_dist, 2, function(x) order(x)[1:3])

base_min_unique <- apply(base_dist, 2, function(x) mean(sort(x)[1:3], na.rm = TRUE))
base_min_unique[!is.finite(base_min_unique)] <- finite_max
base_min_all <- base_min_unique[demand_to_unique_idx]

ncores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(ncores)
clusterExport(cl, c("g", "warehouse_nodes", "demand_nodes_unique",
                     "demand_to_unique_idx", "base_min_all", "pop_weights",
                     "finite_max", "nearest_3", "compute_edge_delta"),
              envir = environment())
clusterEvalQ(cl, library(igraph))

delta_times <- parSapply(cl, seq_len(ecount(g)), compute_edge_delta)
stopCluster(cl)

net <- net |> activate("edges") |> mutate(delta_time_s = delta_times)

# JOIN TO BASELINE & EXPORT --------------------------------------------------

centrality_sf <- net |>
  activate("edges") |>
  as_tibble() |>
  st_as_sf() |>
  select(-c(from, to))

centrality_sf <- centrality_sf |>
  mutate(
    betweenness_wh_norm = min_max_norm(betweenness_wh),
    betweenness_wh3_norm = min_max_norm(betweenness_wh3),
    betweenness_popweighted_norm = min_max_norm(betweenness_popweighted),
    betweenness_popweighted_allwh_norm = min_max_norm(betweenness_popweighted_allwh),
    centrality_score = (
      betweenness_wh_norm *
        betweenness_wh3_norm *
        betweenness_popweighted_norm *
        betweenness_popweighted_allwh_norm
    )^(1/4)
  )

baseline_joined <- baseline_sf |>
  st_join(
    centrality_sf |>
      select(betweenness_wh_norm, betweenness_wh3_norm,
             betweenness_popweighted_norm, betweenness_popweighted_allwh_norm,
             centrality_score),
    join = st_nearest_feature
  )

st_write(baseline_joined, "data/mdg_road_centrality.gpkg",
         layer = "mdg_road_centrality", delete_dsn = TRUE)
