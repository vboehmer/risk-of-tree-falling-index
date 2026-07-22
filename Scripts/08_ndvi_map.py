# ===========================================================================
# ndvi_map.py
# Downloads MODIS NDVI spatial standard-deviation tiles for Madagascar
# from Google Earth Engine and merges them into a single GeoTIFF.
#
# Output: data/ndvi_madagascar_2023_merged.tif
# ===========================================================================

import ee
import geemap
import rasterio
from rasterio.merge import merge
import glob
import os

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------
try:
    ee.Initialize(project="ee-aa-automatization")
except Exception:
    ee.Authenticate()
    ee.Initialize(project="ee-aa-automatization")

# ------------------------------------------------------------------
# 2. Compute NDVI std-dev
# ------------------------------------------------------------------
madagascar = ee.Geometry.Rectangle([43.2, -26.0, 52, -11])

ndvi = (
    ee.ImageCollection("MODIS/061/MOD13Q1")
    .select("NDVI")
    .filterBounds(madagascar)
    .filterDate("2023-11-01", "2024-03-31")
    .mean()
    .multiply(0.0001)
)

kernel = ee.Kernel.circle(radius=1000, units="meters")
ndvi_std = ndvi.reduceNeighborhood(reducer=ee.Reducer.stdDev(), kernel=kernel).clip(
    madagascar
)

print("NDVI spatial std-dev image ready.")

# ------------------------------------------------------------------
# 3. Export tiles
# ------------------------------------------------------------------
tiles = geemap.fishnet(madagascar, h_spacing=0.5, v_spacing=0.5, delta=0)
tile_count = tiles.size().getInfo()
tiles_list = tiles.toList(tile_count)

print("Number of tiles:", tile_count)

export_dir = "data/tmp_ndvi_tiles"
os.makedirs(export_dir, exist_ok=True)

for i in range(tile_count):
    tile = ee.Feature(tiles_list.get(i))
    region = tile.geometry()
    filename = f"{export_dir}/ndvi_std_tile_{i}.tif"
    print("Exporting:", filename)
    try:
        geemap.ee_export_image(
            ndvi_std, filename=filename, region=region, scale=250, crs="EPSG:4326"
        )
    except Exception as e:
        print("Tile", i, "failed:", e)

print("Tile export done.")

# ------------------------------------------------------------------
# 4. Merge tiles
# ------------------------------------------------------------------
tile_files = sorted(glob.glob(f"{export_dir}/*.tif"))
src_files = [rasterio.open(t) for t in tile_files]

mosaic, out_transform = merge(src_files)

out_meta = src_files[0].meta.copy()
out_meta.update(
    {
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_transform,
    }
)

merged_path = "data/ndvi_madagascar_2023_merged.tif"

with rasterio.open(merged_path, "w", **out_meta) as dest:
    dest.write(mosaic)

print(f"Merged file saved: {merged_path}")

# Clean up temporary tiles
for f in tile_files:
    os.remove(f)
os.rmdir(export_dir)
