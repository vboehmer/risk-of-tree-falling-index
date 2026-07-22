import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import geopandas as gpd

gdf = gpd.read_file("data/Output/mdg_baseline_rtfi.gpkg").to_crs(29738)


x = gdf.centrality_score_norm
y = gdf.cyclone_risk_score_norm
z = gdf.THS_risk_norm

# Select top 100 by rtfi_index
top1000 = gdf.nlargest(1000, "rtfi_index")

fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(projection="3d")

# Base scatter (all points)
ax.scatter(x, y, z, linewidths=0.2, alpha=0.7, edgecolor="k", s=20, c=z)

# Overlay top 100 points in red
ax.scatter(
    top1000.centrality_score_norm,
    top1000.cyclone_risk_score_norm,
    top1000.THS_risk_norm,
    linewidths=0.8,
    alpha=0.7,
    edgecolor="red",
    s=20,
    c=top1000.THS_risk_norm,
    label="Top 1000 RTFI values",
)

ax.set_xlabel("CS", fontweight="bold")
ax.set_ylabel("STS", fontweight="bold")
ax.set_zlabel("THS", rotation=90, fontweight="bold")

ax.invert_yaxis()

ax.set_title("3D Distribution of RTFI Indicators", fontsize=14, fontweight="bold")

ax.title.set_y(0.85)

ax.legend(loc="center right", bbox_to_anchor=(0.98, 0.9), frameon=True)

# Export the figure
plt.savefig(
    "Figures/3d_scatter.png",
    dpi=300,
)

plt.show()
