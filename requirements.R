# ===========================================================================
# requirements.R
# Install all R packages needed for the RTFI analysis.
# Run once: source("requirements.R")
# ===========================================================================

required_packages <- c(
  # Spatial data
  "sf",            # Simple features — reading/writing/vector ops
  "terra",         # Raster data
  "sfnetworks",    # Tidy spatial networks
  "tidygraph",     # Tidy interface to igraph

  # Data manipulation
  "dplyr",
  "tidyverse",
  "purrr",
  "scales",

  # Network analysis
  "igraph",        # Graph algorithms

  # OSM data
  "ohsome",        # ohsome API client

  # Parallel computing
  "parallel",
  "future",

  # Visualisation
  "patchwork",     # Combining ggplot2 plots
  "classInt",      # Classification intervals (used in RTFI)
  "ggspatial"      # Spatial annotations on ggplot (if used)
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    message("Already installed: ", pkg)
  }
}
