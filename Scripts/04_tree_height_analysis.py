import ee
import geemap
import geopandas as gpd
import pandas as pd
import numpy as np
import os
import pickle
import math

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------
# Authenticate with Earth Engine (only needed once per machine).
# If running locally, call:  earthengine authenticate
try:
    ee.Initialize(project="ee-aa-automatization")
except Exception:
    ee.Authenticate()
    ee.Initialize(project="ee-aa-automatization")
results_file = "processed_chunks.pkl"
chunk_size = 50
buffer_m = 10  # max buffer distance (m)
distance_zones = [2, 5, 10]  # sub-buffers for 1m canopy weighting
distance_weights = {"2m": 1.0, "5m": 0.6, "10m": 0.3}

# ------------------------------------------------------------------
# 2. Input data
# ------------------------------------------------------------------
gdf = gpd.read_file("data/mdg_baseline_roads.gpkg").to_crs(29738)
n = len(gdf)
total_chunks = math.ceil(n / chunk_size)
max_chunks = None

# ------------------------------------------------------------------
# 3. Load datasets
# ------------------------------------------------------------------
canopy_10m = ee.Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1").rename("chm")
canopy_1m = ee.ImageCollection(
    "projects/sat-io/open-datasets/facebook/meta-canopy-height"
).median()


def make_s2_composite(
    start_date, end_date, region, cloud_prob_thresh=40, percentile=75
):
    """Create Sentinel-2 NDVI composite (pre-cyclone season)."""
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
        .filterDate(start_date, end_date)
        .filterBounds(region)
    )
    clouds = (
        ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
        .filterDate(start_date, end_date)
        .filterBounds(region)
    )
    joined = ee.Join.saveFirst("cloud_mask").apply(
        primary=s2,
        secondary=clouds,
        condition=ee.Filter.equals(leftField="system:index", rightField="system:index"),
    )

    def apply_cloud_mask(img):
        cloud_mask = ee.Image(img.get("cloud_mask")).select("probability")
        mask = cloud_mask.lt(cloud_prob_thresh)
        return ee.Image(img).updateMask(mask)

    s2_masked = ee.ImageCollection(joined).map(apply_cloud_mask).select(["B4", "B8"])
    comp = (
        s2_masked.median()
        if percentile == 50
        else s2_masked.reduce(ee.Reducer.percentile([percentile]))
        .select([f"B4_p{percentile}", f"B8_p{percentile}"])
        .rename(["B4", "B8"])
    )
    return comp


# ------------------------------------------------------------------
# 4. Resume progress if checkpoint exists
# ------------------------------------------------------------------
results, start_idx, chunk_count = [], 0, 0
if os.path.exists(results_file):
    with open(results_file, "rb") as f:
        saved = pickle.load(f)
    results = saved["results"]
    start_idx = saved["last_idx"]
    chunk_count = saved["chunk_count"]
    print(
        f"Resuming from chunk starting at index {start_idx} ({chunk_count}/{total_chunks})"
    )

# ------------------------------------------------------------------
# 5. Process chunks
# ------------------------------------------------------------------
while start_idx < n and (max_chunks is None or chunk_count < max_chunks):
    end_idx = min(start_idx + chunk_size, n)
    gdf_chunk = gdf.iloc[start_idx:end_idx]
    print(
        f"Processing chunk {start_idx}-{end_idx} (# {chunk_count + 1}/{total_chunks}) ..."
    )

    try:
        fc = geemap.geopandas_to_ee(gdf_chunk)
        fc_buffered = fc.map(lambda f: f.buffer(buffer_m))

        # ------------------------- 5.1 Canopy 10m -------------------------
        reducer = (
            ee.Reducer.mean()
            .combine(ee.Reducer.max(), sharedInputs=True)
            .combine(ee.Reducer.stdDev(), sharedInputs=True)
        )

        stats_fc_10m = canopy_10m.reduceRegions(
            collection=fc_buffered, reducer=reducer, scale=10
        )
        tmp_csv_10m = f"temp_canopy10m_{start_idx}.csv"
        geemap.ee_to_csv(stats_fc_10m, filename=tmp_csv_10m)
        df_chunk_10m = pd.read_csv(tmp_csv_10m).rename(
            columns={
                "mean": "canopy10m_mean",
                "max": "canopy10m_max",
                "stdDev": "canopy10m_std",
            }
        )
        os.remove(tmp_csv_10m)

        # ------------------------- 5.2 Canopy 1m multi-distance -------------------------
        results_by_zone = []
        for dist in distance_zones:
            fc_zone = fc.map(lambda f: f.buffer(dist))
            stats_fc = canopy_1m.reduceRegions(
                collection=fc_zone, reducer=reducer, scale=1
            )

            tmp_csv = f"temp_canopy1m_{dist}m_{start_idx}.csv"
            geemap.ee_to_csv(stats_fc, filename=tmp_csv)
            df_zone = pd.read_csv(tmp_csv)
            os.remove(tmp_csv)

            df_zone = df_zone.drop(
                columns=[
                    c
                    for c in ["system:index", "length_m", "geometry"]
                    if c in df_zone.columns
                ]
            )

            df_zone = df_zone.rename(
                columns={
                    "mean": f"canopy_mean_{dist}m",
                    "max": f"canopy_max_{dist}m",
                    "stdDev": f"canopy_std_{dist}m",
                }
            )
            results_by_zone.append(df_zone)

        # Merge distance zones
        df_chunk_1m = results_by_zone[0]
        for df_z in results_by_zone[1:]:
            df_chunk_1m = (
                df_chunk_1m.set_index("id")
                .join(df_z.set_index("id"), how="outer")
                .reset_index()
            )

        # ------------------------- 5.3 Apply distance weighting -------------------------
        total_w = sum(distance_weights.values())

        def weighted_sum(row, metric):
            return (
                sum(
                    [
                        row.get(f"canopy_{metric}_{d}", 0) * w
                        for d, w in distance_weights.items()
                    ]
                )
                / total_w
            )

        for metric in ["mean", "max", "std"]:
            df_chunk_1m[f"canopy1m_{metric}_w"] = df_chunk_1m.apply(
                lambda r: weighted_sum(r, metric), axis=1
            )

        # ------------------------- 5.4 NDVI composite -------------------------
        start_date, end_date = "2024-11-01", "2025-04-01"
        region = fc_buffered.geometry().bounds()
        s2_comp = make_s2_composite(start_date, end_date, region)
        ndvi = s2_comp.expression(
            "(NIR - RED) / (NIR + RED)",
            {"NIR": s2_comp.select("B8"), "RED": s2_comp.select("B4")},
        ).rename("ndvi")

        ndvi_stats = ndvi.reduceRegions(
            collection=fc_buffered,
            reducer=ee.Reducer.mean().combine(ee.Reducer.stdDev(), sharedInputs=True),
            scale=10,
        )
        tmp_csv_ndvi = f"temp_ndvi_{start_idx}.csv"
        geemap.ee_to_csv(ndvi_stats, filename=tmp_csv_ndvi)
        df_ndvi = pd.read_csv(tmp_csv_ndvi).rename(
            columns={"mean": "ndvi_mean", "stdDev": "ndvi_std"}
        )
        os.remove(tmp_csv_ndvi)

        # ------------------------- 5.5 Merge all results -------------------------
        for df in [df_chunk_1m, df_ndvi]:
            for col in ["system:index", "length_m", "geometry"]:
                df.drop(columns=[col], inplace=True, errors="ignore")

        df_chunk = (
            df_chunk_10m.set_index("id")
            .join(df_chunk_1m.set_index("id"), how="left")
            .join(df_ndvi.set_index("id"), how="left")
        )

        results.append(df_chunk)

        # ------------------------- 5.6 Save checkpoint -------------------------
        with open(results_file, "wb") as f:
            pickle.dump(
                {
                    "results": results,
                    "last_idx": end_idx,
                    "chunk_count": chunk_count + 1,
                },
                f,
            )

    except Exception as e:
        print(f"Error processing chunk {start_idx}-{end_idx}: {e}")

    start_idx = end_idx
    chunk_count += 1

print(f"Processed {chunk_count}/{total_chunks} chunks.")
print(f"Intermediate results saved to {results_file}")

# ------------------------------------------------------------------
# 6. Combine all chunks
# ------------------------------------------------------------------
if not results:
    raise RuntimeError("No canopy or NDVI stats computed.")
df_all = pd.concat(results)
df_all.drop(columns=["length_m", "geometry"], inplace=True, errors="ignore")
gdf = gdf.set_index("id").join(df_all, how="left", rsuffix="_new").fillna(0)


# ------------------------------------------------------------------
# 7. Compute Tree Height Risk (THS)
# ------------------------------------------------------------------
def norm_series(s):
    return (
        (s - s.min()) / (s.max() - s.min())
        if s.max() > s.min()
        else pd.Series(0.0, index=s.index)
    )


mask_tree = gdf["canopy10m_max"] >= 5.0
gdf["tree_flag"] = mask_tree
gdf_tree = gdf[gdf["tree_flag"]].copy()

if not gdf_tree.empty:
    # Normalize
    gdf_tree["canopy10m_mean_norm"] = norm_series(gdf_tree["canopy10m_mean"])
    gdf_tree["canopy1m_mean_norm"] = norm_series(gdf_tree["canopy1m_mean_w"])
    gdf_tree["canopy1m_max_norm"] = norm_series(gdf_tree["canopy1m_max_w"])
    gdf_tree["canopy1m_sd_norm"] = norm_series(gdf_tree["canopy1m_std_w"])
    gdf_tree["ndvi_sd_norm"] = norm_series(gdf_tree["ndvi_std"])

    w = {
        "canopy10m_mean_norm": 0.15,
        "canopy1m_mean_norm": 0.2,
        "canopy1m_max_norm": 0.25,
        "canopy1m_sd_norm": 0.25,
        "ndvi_sd_norm": 0.15,
    }
    tot = sum(w.values())

    gdf_tree["THS_risk"] = (
        gdf_tree["canopy10m_mean_norm"] * w["canopy10m_mean_norm"]
        + gdf_tree["canopy1m_mean_norm"] * w["canopy1m_mean_norm"]
        + gdf_tree["canopy1m_max_norm"] * w["canopy1m_max_norm"]
        + gdf_tree["canopy1m_sd_norm"] * w["canopy1m_sd_norm"]
        + gdf_tree["ndvi_sd_norm"] * w["ndvi_sd_norm"]
    ) / tot

else:
    gdf_tree["THS_risk"] = pd.Series(dtype=float)

# Merge back
gdf["THS_risk"] = 0.0
gdf.loc[gdf_tree.index, "THS_risk"] = gdf_tree["THS_risk"]

# ------------------------------------------------------------------
# 8. Export
# ------------------------------------------------------------------
output_file = "data/Output/mdg_baseline_canopy_height.gpkg"
os.makedirs(os.path.dirname(output_file), exist_ok=True)
gdf.reset_index().to_file(output_file, driver="GPKG", overwrite=True)

print(f"Tree Height Risk Score (THS_risk) saved to {output_file}")
