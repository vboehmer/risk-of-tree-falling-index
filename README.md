# Risk of Tree Falling Index (RTFI) for Madagascar

A composite spatial index that estimates the risk of trees falling onto Madagascar's primary, secondary, and tertiary road network during tropical cyclones. The index combines three indicators: canopy characteristics near roads, cyclone exposure, and road network centrality relative to health commodity warehouses.

This repository contains all processing scripts, output figures, and the thesis document.

## Study Area

Madagascar. The analysis covers all classified roads (primary, secondary, tertiary, plus ferries, bridges, and fords) obtained from OpenStreetMap via the ohsome API.

## Methodology

The RTFI is computed as the geometric mean of three normalised scores assigned to each 100-m road segment:

| Indicator | Abbreviation | What it captures |
|-----------|-------------|------------------|
| Tree Height Score | THS | Canopy structure within 10 m of the road (GEDI 10 m, Meta 1 m canopy height, Sentinel-2 NDVI) |
| Storm Tracks Score | STS | Cyclone exposure (track count, cumulative Saffir-Simpson intensity, distance to storm centre) |
| Centrality Score | CS | Road importance for warehouse-to-population access (population-weighted betweenness, travel-time impact) |

```
RTFI = (THS × STS × CS)^(1/3)
```

All input indicators are normalised to [0, 1] before combination. A road segment receives THS = 0 (and therefore RTFI = 0) when no canopy is detected within its 10-m buffer.

## Repository Structure

```
├── Scripts/
│   ├── 01_network_preprocessing.R      # OSM download, network cleaning, baseline segments
│   ├── 02_centrality_analysis.R        # Betweenness centrality, travel-time impact (ΔT)
│   ├── 03_cyclone_risk.R               # Cyclone track buffering, storm score
│   ├── 04_tree_height_analysis.py      # GEE canopy height + NDVI extraction, THS
│   ├── 05_rtfi_calculation.R           # Combine indicators into RTFI
│   ├── 06_visualization.R              # Histogram plots for all indicators
│   ├── 07_3d_plot.py                   # 3D scatter of RTFI indicators
│   └── 08_ndvi_map.py                  # NDVI spatial std-dev tile export (GEE)
├── Figures/                             # All output maps, histograms, workflow diagrams
├── data/
│   ├── README.md                        # Data sources, download links, file descriptions
│   └── Output/                          # Intermediate geopackages (generated, git-ignored)
├── Literature/                          # Reference papers (git-ignored, see data/README.md)
├── Masterarbeit.Rproj                   # RStudio project file
├── master_thesis.bib                    # Bibliography
├── requirements.R                       # R package installer
├── requirements.txt                     # Python package installer
└── README.md
```

## Prerequisites

- **R** (>= 4.2) with packages listed in `requirements.R`
- **Python** (>= 3.9) with packages listed in `requirements.txt`
- A **Google Earth Engine** account (for scripts 04 and 08)
- An **ohsome API** key or public access (script 01 uses the ohsome R package)

## Replication

1. Clone the repository and set the working directory to the repo root.

2. Install R and Python dependencies:
   ```r
   source("requirements.R")
   ```
   ```bash
   pip install -r requirements.txt
   ```

3. Download input data as described in `data/README.md`.

4. Run scripts in order:
   ```
   01_network_preprocessing.R
   02_centrality_analysis.R
   03_cyclone_risk.R
   04_tree_height_analysis.py    # requires GEE authentication
   05_rtfi_calculation.R
   06_visualization.R            # optional
   07_3d_plot.py                 # optional
   ```

   Scripts 06 and 07 produce figures only. All other scripts write geopackages to `data/` or `data/Output/`.

## Data Sources

| Dataset | Source |
|---------|--------|
| Road network | OpenStreetMap via [ohsome API](https://docs.ohsome.org/ohsome-api/) |
| Warehouse locations | CRM Madagascar |
| Cyclone tracks | [IBTrACS](https://www.ncei.noaa.gov/products/international-best-track-archive) v04r01 |
| Canopy height (10 m) | [ETH Global Canopy Height 2020](https://doi.org/10.1038/s41586-023-06787-x) |
| Canopy height (1 m) | [Meta / Open Buildings](https://dataFORDEM.com/dataset/1431e721-1e0f-46e0-a751-836d7c46b945) |
| NDVI | MODIS MOD13Q1 (250 m) and Sentinel-2 (10 m) |
| Population | [WorldPop](https://www.worldpop.org/) constrained 100 m |
| Admin boundaries | [DIVA-GIS](https://diva-gis.org/gdata) |

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use this work, please cite the thesis:

> Boehmer, A. (2025). *Risk of Tree Falling Index for Madagascar's road network during tropical cyclones.* Master's thesis, University of Göttingen.

The full thesis is available for download: [thesis PDF](https://heibox.uni-heidelberg.de/f/56ee9f61270b45bbbe92/?dl=1).

A `CITATION.cff` file is included for machine-readable citation metadata.
