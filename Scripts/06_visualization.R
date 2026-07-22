library(tidyverse)
library(sf)
library(patchwork)


# Load dataset
rtfi_sf <- st_read("data/Output/mdg_baseline_rtfi.gpkg")

# Centrality Score -----------------------------------------------

# Calculate mean and median for both variables
mean_cs <- mean(rtfi_sf$centrality_score, na.rm = TRUE)
median_cs <- median(rtfi_sf$centrality_score, na.rm = TRUE)

mean_cs_norm <- mean(rtfi_sf$centrality_score_norm, na.rm = TRUE)
median_cs_norm <- median(rtfi_sf$centrality_score_norm, na.rm = TRUE)

# Plot 1: original centrality score
p1 <- ggplot(rtfi_sf, aes(x = centrality_score)) +
  geom_histogram(bins = 40) +
  geom_vline(aes(xintercept = mean_cs), color = "darkgrey", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median_cs), color = "darkgrey", linetype = "dotted", size = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  labs(title = "Centrality Score", x = "CS", y = "Frequency")

# Plot 2: quantile normalized
p2 <- ggplot(rtfi_sf, aes(x = centrality_score_norm)) +
  geom_histogram(bins = 40) +
  geom_vline(aes(xintercept = mean_cs_norm), color = "darkgrey", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median_cs_norm), color = "darkgrey", linetype = "dotted", size = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  labs(title = "Centrality Score (Quantile Normalized)", x = "CS", y = "Frequency")

# Combine plots side by side
cs_hist <- p1 + p2 

# Save as PNG
ggsave(
  filename = "Figures/centrality_histograms.png",
  plot = cs_hist,
  width = 12,    # adjust size as needed
  height = 6,
  dpi = 600
)


# Tree Height Score -----------------------------------------------

# Calculate mean and median for both variables
mean_ths <- mean(rtfi_sf$THS_risk, na.rm = TRUE)
median_ths <- median(rtfi_sf$THS_risk, na.rm = TRUE)

mean_ths_norm <- mean(rtfi_sf$THS_risk_norm, na.rm = TRUE)
median_ths_norm <- median(rtfi_sf$THS_risk_norm, na.rm = TRUE)

# Plot 1: original ths score
p1 <- ggplot(rtfi_sf, aes(x = THS_risk)) +
  geom_histogram(bins = 40) +
  geom_vline(aes(xintercept = mean_ths), color = "darkgrey", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median_ths), color = "darkgrey", linetype = "dotted", size = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  labs(title = "Tree Height Score", x = "THS", y = "Frequency")

# Plot 2: normalized
p2 <- ggplot(rtfi_sf, aes(x = THS_risk_norm)) +
  geom_histogram(bins = 40) +
  geom_vline(aes(xintercept = mean_ths_norm), color = "darkgrey", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median_ths_norm), color = "darkgrey", linetype = "dotted", size = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  labs(title = "Tree Height Score (Min-Max Normalized)", x = "THS", y = "Frequency")
# Combine plots side by side
ths_hist <- p1 + p2 

# Save as PNG
ggsave(
  filename = "Figures/treeheight_histograms.png",
  plot = ths_hist,
  width = 12,    # adjust size as needed
  height = 6,
  dpi = 600
)



# Storm Tracks Score -----------------------------------------------

# Calculate mean and median for both variables
mean_sts <- mean(rtfi_sf$cyclone_risk_score, na.rm = TRUE)
median_sts <- median(rtfi_sf$cyclone_risk_score, na.rm = TRUE)

mean_sts_norm <- mean(rtfi_sf$cyclone_risk_score_norm, na.rm = TRUE)
median_sts_norm <- median(rtfi_sf$cyclone_risk_score_norm, na.rm = TRUE)

# Plot 1: original sts score
p1 <- ggplot(rtfi_sf, aes(x = cyclone_risk_score)) +
  geom_histogram(bins = 40) +
  geom_vline(aes(xintercept = mean_sts), color = "darkgrey", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median_sts), color = "darkgrey", linetype = "dotted", size = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  labs(title = "Storm Tracks Score", x = "STS", y = "Frequency")

# Plot 2: normalized
p2 <- ggplot(rtfi_sf, aes(x = cyclone_risk_score_norm)) +
  geom_histogram(bins = 40) +
  geom_vline(aes(xintercept = mean_sts_norm), color = "darkgrey", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median_sts_norm), color = "darkgrey", linetype = "dotted", size = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  labs(title = "Storm Tracks Score (Min-Max Normalized)", x = "STS", y = "Frequency")
# Combine plots side by side
sts_hist <- p1 + p2 

# Save as PNG
ggsave(
  filename = "Figures/stormtracks_histograms.png",
  plot = sts_hist,
  width = 12,    # adjust size as needed
  height = 6,
  dpi = 600
)



#  RTFI -----------------------------------------------

# Calculate mean and median for both variables
mean_rtfi <- mean(rtfi_sf$rtfi_index, na.rm = TRUE)
median_rtfi <- median(rtfi_sf$rtfi_index, na.rm = TRUE)

mean_sts_norm <- mean(rtfi_sf$cyclone_risk_score_norm, na.rm = TRUE)
median_sts_norm <- median(rtfi_sf$cyclone_risk_score_norm, na.rm = TRUE)

# Plot 1: original sts score
p1 <- ggplot(rtfi_sf, aes(x = rtfi_index)) +
  geom_histogram(bins = 40) +
  geom_vline(aes(xintercept = mean_rtfi), color = "darkgrey", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median_rtfi), color = "darkgrey", linetype = "dotted", size = 1) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  labs(title = "Risk of Tree Falling Index", x = "RTFI", y = "Frequency")


# Save as PNG
ggsave(
  filename = "Figures/rtfi_histograms.png",
  plot = p1,
  width = 8,    # adjust size as needed
  height = 6,
  dpi = 600
)









