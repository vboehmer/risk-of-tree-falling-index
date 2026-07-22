library(sf)
library(dplyr)
library(classInt)
library(scales)

# load data ---------------------------------------------
baseline_sf <- st_read("data/mdg_baseline_roads.gpkg")
tree_height <- st_read("data/Output/mdg_baseline_canopy_height.gpkg")
cyclones <- st_read("data/mdg_roads_cyclone_risk.gpkg")
centrality <- st_read("data/mdg_road_centrality.gpkg")


# Preparation ---------------------------------------------
tree_height <- tree_height |> select(id, THS_risk) |> st_drop_geometry()
cyclones <- cyclones |> filter(!is.na(id)) |> distinct(id, .keep_all = TRUE) |> select(id, cyclone_risk_score) |> st_drop_geometry()
centrality <- centrality |> select(id, centrality_score) |> st_drop_geometry()

#join all indicators
rtfi_sf <- baseline_sf |> 
  left_join(tree_height, by="id") |> 
  left_join(cyclones, by="id") |> 
  left_join(centrality, by="id")

# Normalization -----------------------------------------
## centrality score normalization  ---------------------

# test different methods
df <- rtfi_sf |> st_drop_geometry() # for convenience

# Summary function returning 6 values
summary_methods <- function(x) {
  c(
    min = min(x, na.rm = TRUE),
    Q1 = quantile(x, .25, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    Q3 = quantile(x, .75, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

# 1. Original min-max normalization
df$norm_original <- rescale(df$centrality_score)

# 2. Log transform + min-max normalization
df$log1 <- log1p(df$centrality_score)
df$norm_log1 <- rescale(df$log1)

# 3. Double-log transform
df$log2 <- log1p(log1p(df$centrality_score + 1e-12))
df$norm_log2 <- rescale(df$log2)

# 4. Winsorization @ 99th percentile
p99 <- quantile(df$centrality_score, 0.99, na.rm = TRUE)
df$wins <- pmin(df$centrality_score, p99)
df$norm_wins <- rescale(df$wins)

# 5. Quantile normalization (5–95)
q05 <- quantile(df$centrality_score, 0.05, na.rm = TRUE)
q95 <- quantile(df$centrality_score, 0.95, na.rm = TRUE)

norm_q <- (df$centrality_score - q05) / (q95 - q05)
df$norm_quantile <- pmin(pmax(norm_q, 0), 1)

# 6. Rank percentile normalization
df$norm_rank <- rank(df$centrality_score, ties.method = "average") / nrow(df)

# Build comparison table
methods <- rbind(
  summary_methods(df$norm_original),
  summary_methods(df$norm_log1),
  summary_methods(df$norm_log2),
  summary_methods(df$norm_wins),
  summary_methods(df$norm_quantile),
  summary_methods(df$norm_rank)
)

result <- data.frame(
  method = c("original", "log1", "log2", "winsorized_99", "quantile_5_95", "rank"),
  round(methods, 5)
)



# used Quantile normalization (5–95)
q05 <- quantile(rtfi_sf$centrality_score, 0.05, na.rm = TRUE)
q95 <- quantile(rtfi_sf$centrality_score, 0.95, na.rm = TRUE)

norm_q <- (rtfi_sf$centrality_score - q05) / (q95 - q05)
rtfi_sf$centrality_score_norm <- pmin(pmax(norm_q, 0), 1)

# THS and SHS normalization (min-max) ------------------------------
rtfi_sf <- rtfi_sf |> 
  mutate(
    THS_risk_norm = (THS_risk - min(THS_risk, na.rm = TRUE)) /
      (max(THS_risk, na.rm = TRUE) - min(THS_risk, na.rm = TRUE)),
    
    cyclone_risk_score_norm = (cyclone_risk_score - min(cyclone_risk_score, na.rm = TRUE)) /
      (max(cyclone_risk_score, na.rm = TRUE) - min(cyclone_risk_score, na.rm = TRUE))
  )

# Risk of Tree Falling Index (RTFI) ------------------------------

# calculate rtfi (geometric mean; RTFI = 0 when no canopy detected)
rtfi_sf <- rtfi_sf |>
  mutate(
    THS_risk_norm_f = ifelse(THS_risk_norm < 1e-12, 0, THS_risk_norm),
    centrality_score_norm_f = ifelse(centrality_score_norm < 1e-12, 0, centrality_score_norm),
    cyclone_risk_score_norm_f = ifelse(cyclone_risk_score_norm < 1e-12, 0, cyclone_risk_score_norm),

    rtfi_index =
      if_else(
        THS_risk_norm_f == 0,
        0,
        (THS_risk_norm_f * centrality_score_norm_f * cyclone_risk_score_norm_f)^(1/3)
      )
  )

st_write(rtfi_sf, "data/Output/mdg_baseline_rtfi.gpkg", delete_dsn = TRUE)



# THS quality assessment ------------------------------
library(sf)
library(igraph)

set.seed(Sys.time())

# Ensure valid geometries
rtfi_sf <- st_make_valid(rtfi_sf)

touch_list <- st_touches(rtfi_sf)

# Build edge list
edges <- do.call(
  rbind,
  lapply(seq_along(touch_list), function(i) {
    if (length(touch_list[[i]]) > 0) {
      cbind(i, touch_list[[i]])
    }
  })
)

# Create graph
g <- graph_from_edgelist(edges, directed = FALSE)

# ---- Function to sample 20 connected segments ----
sample_connected_20 <- function() {
  repeat {
    # Random starting node
    start <- sample(V(g), 1)
    
    bfs_nodes <- bfs(
      g,
      root = start,
      neimode = "all"
    )$order
    
    bfs_nodes <- bfs_nodes[!is.na(bfs_nodes)]
    
    # Only accept if we can reach at least 20 connected segments
    if (length(bfs_nodes) >= 20) {
      sel <- bfs_nodes[1:20]
      return(rtfi_sf[as.integer(sel), ])
    }
  }
}

# ---- Draw 6 samples ----
samples <- lapply(1:6, function(i) {
  s <- sample_connected_20()
  s$sample_id <- i
  s
})

# Combine
rtfi_samples_sf <- do.call(rbind, samples)

rtfi_samples_sf

st_write(rtfi_samples_sf, "data/Output/mdg_baseline_canopy_height_qualitycheck.gpkg", delete_dsn = TRUE)

