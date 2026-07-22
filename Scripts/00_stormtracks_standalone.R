
# Packages ------------------------------------------------
library(sf)
library(dplyr)
library(purrr)

# Preparation ------------------------------------------------------

stormtracks <- st_read("data/IBTrACS.since1980.list.v04r01.lines.shp")
baseline_sf <- st_read("data/mdg_baseline_roads.gpkg")

#baseline_sf <- baseline_sf |> head(100)

stormtracks_utm <- st_transform(stormtracks, st_crs(baseline_sf))

 
stormtracks_subset <- stormtracks_utm |> 
  filter(USA_SSHS >= 1)

# Keep geometry for buffering
stormtracks_subset <- stormtracks_utm |>
  filter(USA_SSHS >= 1) |>
  rowwise() |>
  # Compute mean R34
  mutate(R34_buffer = mean(c(USA_R34_NE, USA_R34_NW, USA_R34_SE, USA_R34_SW), na.rm = TRUE)) |>
  ungroup() |>
  group_by(USA_SSHS) |>
  mutate(R34_buffer = ifelse(is.na(R34_buffer), mean(R34_buffer, na.rm = TRUE), R34_buffer)) |>
  ungroup() |>
  mutate(R34_buffer_km = R34_buffer * 1.852) |>
  rowwise() |>
  # Compute mean R64
  mutate(R64_buffer = mean(c(USA_R64_NE, USA_R64_NW, USA_R64_SE, USA_R64_SW), na.rm = TRUE)) |>
  ungroup() |>
  group_by(USA_SSHS) |>
  mutate(R64_buffer = ifelse(is.na(R64_buffer), mean(R64_buffer, na.rm = TRUE), R64_buffer)) |>
  ungroup() |>
  mutate(R64_buffer_km = R64_buffer * 1.852)


# Count number of points per geometry
geom_points <- sapply(st_geometry(stormtracks_subset), function(g) {
  coords <- st_coordinates(g)
  nrow(coords)
})

# Keep only geometries with >1 point
stormtracks_clean <- stormtracks_subset[geom_points > 1, ]

# Now filter by LINESTRING/MULTILINESTRING types
stormtracks_clean <- stormtracks_clean |>
  filter(st_geometry_type(geometry) %in% c("LINESTRING", "MULTILINESTRING"))


# Define SWIO bounding box in WGS84
swio_bbox <- st_bbox(c(xmin = 30, ymin = -40, xmax = 90, ymax = 10), crs = st_crs(4326))

# Convert bbox to polygon
swio_poly <- st_as_sfc(swio_bbox)

# Transform to your stormtracks CRS (UTM zone 38S)
swio_poly_utm <- st_transform(swio_poly, st_crs(stormtracks_utm))

stormtracks_swio <- stormtracks_clean |>
  st_filter(swio_poly_utm) |>
  select(
    SID, SEASON, USA_SSHS, WMO_WIND,
    USA_R34_NE, USA_R34_SE, USA_R34_SW, USA_R34_NW,
    USA_R64_NE, USA_R64_SE, USA_R64_SW, USA_R64_NW,
    R34_buffer, R34_buffer_km, R64_buffer, R64_buffer_km,
    geometry
  ) 


# Create buffers row-wise using map2 to avoid errors
stormtracks_swio <- stormtracks_swio |>
  rowwise() |>
  mutate(
    R34_geom = st_buffer(geometry, R34_buffer_km * 1000),
    R64_geom = st_buffer(geometry, R64_buffer_km * 1000)
  ) |>
  ungroup()


# R34 ---------------------------------------------------------------------
# First, union R34 polygons by SID so each storm has one merged buffer
stormtracks_r34_union <- stormtracks_swio |>
  group_by(SID) |>
  summarise(R34_geom = st_union(R34_geom)) |>
  st_as_sf()

# Compute intersections between baseline segments and R34 storm polygons
# Each row in the result corresponds to one road–storm intersection
intersections <- st_intersects(baseline_sf, stormtracks_r34_union)

# Count how many unique storm SIDs intersect each road segment
baseline_sf$R34_count <- lengths(intersections)

# Make sure R34_geom exists in the sf object
stormtracks_swio <- stormtracks_swio |>
  st_set_geometry("R34_geom")

# We'll only keep relevant columns from stormtracks_swio
stormtracks_r34 <- stormtracks_swio |>
  select(SID, USA_SSHS, R34_geom) |>
  st_as_sf(sf_column_name = "R34_geom")

# Spatial intersection (returns only overlapping features)
intersections <- st_intersection(
  baseline_sf |> select(id),
  stormtracks_r34
)

# Clean up and keep relevant fields
result_df <- intersections |>
  st_drop_geometry() |>
  select(id, SID, USA_SSHS)


# If multiple (id, SID) combinations exist, keep the one with the highest USA_SSHS
result_df <- result_df |>
  group_by(id, SID) |>
  slice_max(order_by = USA_SSHS, n = 1, with_ties = FALSE) |>
  ungroup()

# Now sum USA_SSHS per road segment (id)
df_sum <- result_df |>
  group_by(id) |>
  summarise(R34_USA_SSHS_sum = sum(USA_SSHS, na.rm = TRUE)) |>
  ungroup()

# join back to baseline_sf
baseline_sf <- baseline_sf |>
  left_join(df_sum, by = "id") 


## Compute mean distance between baseline segments and storm track -------------
# We'll use `result_df` which has `id` (road) and `SID` (storm)
# First, join storm geometries to result_df
result_with_geom <- result_df  |>
  left_join(stormtracks_swio |> select(SID, geometry), by = "SID") 

# only keep necessary columns
storms_lines <- result_with_geom |> select(id, SID, geometry)

# Compute the minimum distance between each storm line and baseline with the same id
storms_with_distance <- storms_lines |>
  rowwise() |>  # operate row by row
  mutate(
    road_distance = st_distance(
      geometry, 
      baseline_sf$geom[baseline_sf$id == id], 
      by_element = TRUE
    )
  ) |>
  ungroup()

# If multiple entries per (id, SID), keep the one with the minimum distance
storms_unique <- storms_with_distance |>
  group_by(id, SID) |>
  slice_min(road_distance, n = 1, with_ties = FALSE) |>
  ungroup()

# Now compute mean distance per road segment (id)
storms_mean_distance <- storms_unique |>
  st_drop_geometry() |>
  group_by(id) |>
  summarise(mean_distance = mean(road_distance, na.rm = TRUE)) |>
  ungroup()

# Join mean distance back to baseline_sf
baseline_sf <- baseline_sf |>
  left_join(storms_mean_distance, by = "id")


# R64 ---------------------------------------------------------------------
# First, union R64 polygons by SID so each storm has one merged buffer
stormtracks_r64_union <- stormtracks_swio |>
  group_by(SID) |>
  summarise(R64_geom = st_union(R64_geom)) |>
  st_as_sf()

# Compute intersections between baseline segments and R64 storm polygons
# Each row in the result corresponds to one road–storm intersection
intersections <- st_intersects(baseline_sf, stormtracks_r64_union)

# Count how many unique storm SIDs intersect each road segment
baseline_sf$R64_count <- lengths(intersections)

# Make sure R64_geom exists in the sf object
stormtracks_swio <- stormtracks_swio |>
  st_set_geometry("R64_geom")

# We'll only keep relevant columns from stormtracks_swio
stormtracks_r64 <- stormtracks_swio |>
  select(SID, USA_SSHS, R64_geom) |>
  st_as_sf(sf_column_name = "R64_geom")

# Spatial intersection (returns only overlapping features)
intersections <- st_intersection(
  baseline_sf |> select(id),
  stormtracks_r64
)

# Clean up and keep relevant fields
result_df <- intersections |>
  st_drop_geometry() |>
  select(id, SID, USA_SSHS)


# If multiple (id, SID) combinations exist, keep the one with the highest USA_SSHS
result_df <- result_df |>
  group_by(id, SID) |>
  slice_max(order_by = USA_SSHS, n = 1, with_ties = FALSE) |>
  ungroup()

# Now sum USA_SSHS per road segment (id)
df_sum <- result_df |>
  group_by(id) |>
  summarise(R64_USA_SSHS_sum = sum(USA_SSHS, na.rm = TRUE)) |>
  ungroup()

# join back to baseline_sf
baseline_sf <- baseline_sf |>
  left_join(df_sum, by = "id") 

#replace NAs with 0
baseline_sf$R64_USA_SSHS_sum[is.na(baseline_sf$R64_USA_SSHS_sum)] <- 0


# Ensure your geometry column is correctly set
st_geometry(stormtracks_union_R34) <- "R34_geom"
stormtracks_union_R34
# Plot with mapview
mapview::mapview(stormtracks_swio, zcol = "USA_SSHS", layer.name = "R34 Buffers", alpha.regions = 0.3)
