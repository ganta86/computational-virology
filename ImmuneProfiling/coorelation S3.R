# Load necessary libraries
library(patchwork)
library(ggplot2)
library(pheatmap)
library(readxl)
library(reshape2)
library(ggpubr)
library(matrixStats)
library(GeneNet)
library(RColorBrewer)
library(caret)
library(dplyr)
library(ggraph)
library(igraph)
library(corrr)
library(tidyr)
library(stringr)
library(tidyverse)
library(scales)
library(rstatix)
library(writexl)
library(viridis)
library(vip)
library(Hmisc)
library(extrafont)
library(ggthemes)
library(ggbeeswarm)
library(ggsignif)
library(cowplot)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)

# Load the dataset
COMPILED_DATA_EBNA1b <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
                                   sheet = "Data_2024")
COMPILED_DATA_EBNA1b <- as.data.frame(COMPILED_DATA_EBNA1b)

# Define the list of DRB1_1501_Positive participants
DRB1_1501_Positive_participants <- c("1580", "1586", "1611", "1620", "1652", "1654", "1665", "1555", "1559", "1602", "1617", "1624", "1645", "1653", "1664", "1666", "1679", "1681", "1684", "ES346", "ES357", "ES389", "ES406", "ES411", "ES418", "ES459", "ES372", "ES396", "ES452", "ES488")

# Categorize DRB status and correct Sero_Pos group
COMPILED_DATA_EBNA1b <- COMPILED_DATA_EBNA1b %>%
  mutate(DRB = ifelse(participants %in% DRB1_1501_Positive_participants, "DRB1_1501_Positive", "DRB1_1501_Negative"))

# Filter for combined Sero_Pos category
COMPILED_DATA_sero_pos <- filter(COMPILED_DATA_EBNA1b, EBV_infection == "Sero_Pos (n = 50)")

# Scatter plot for Sero_Pos (Combined Category)
A <- ggplot(COMPILED_DATA_sero_pos, aes(x = IgG1_EBNA1_460_641, y = IgG1_EBNA1_377_459)) +
  geom_point(aes(color = DRB)) +
  scale_color_manual(values = c("DRB1_1501_Positive" = "dodgerblue3", "DRB1_1501_Negative" = "red")) +
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10') +
  coord_cartesian(xlim = c(100, 1e7), ylim = c(100, 1e7)) +
  labs(x = "IgG1_EBNA1_460_641", y = "IgG1_EBNA1_377_459", title = "Sero_Pos (n = 50)") +
  geom_smooth(method = "lm", col = "darkgrey", fill = "grey") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 15, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.y = element_blank(),   # Remove y-axis labels
    title = element_text(size = 12),
    legend.position = "right"
  ) +
  stat_cor(method = "spearman", label.x = 3, label.y = 6, size = 3)

# Filter for 6 months category
COMPILED_DATA_6M <- filter(COMPILED_DATA_EBNA1b, EBV_infection == "6 months (n = 30)")

# Scatter plot for 6 months with labeled y-axis
B <- ggplot(COMPILED_DATA_6M, aes(x = IgG1_EBNA1_460_641, y = IgG1_EBNA1_377_459)) +
  geom_point(aes(color = DRB)) +
  scale_color_manual(values = c("DRB1_1501_Positive" = "dodgerblue3", "DRB1_1501_Negative" = "red")) +
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10') +
  coord_cartesian(xlim = c(100, 1e7), ylim = c(100, 1e7)) +
  labs(x = "IgG1_EBNA1_460_641", y = "IgG1_EBNA1_377_459", title = "6 months (n = 30)") +
  geom_smooth(method = "lm", col = "darkgrey", fill = "grey") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 15, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(size = 12),
    legend.position = "none"  # Remove legend to avoid redundancy
  ) +
  stat_cor(method = "spearman", label.x = 3, label.y = 6, size = 3)

# Filter for 1 year category
COMPILED_DATA_1Y <- filter(COMPILED_DATA_EBNA1b, EBV_infection == "1 year (n = 67)")

# Scatter plot for 1 year
C <- ggplot(COMPILED_DATA_1Y, aes(x = IgG1_EBNA1_460_641, y = IgG1_EBNA1_377_459)) +
  geom_point(aes(color = DRB)) +
  scale_color_manual(values = c("DRB1_1501_Positive" = "dodgerblue3", "DRB1_1501_Negative" = "red", "NEGATIVE" = "black")) +
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10') +
  coord_cartesian(xlim = c(100, 1e7), ylim = c(100, 1e7)) +
  labs(x = "IgG1_EBNA1_460_641", y = "IgG1_EBNA1_377_459", title = "1 year (n = 67)") +
  geom_smooth(method = "lm", col = "darkgrey", fill = "grey") + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 15, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12, hjust = 0.5),
    legend.position = "none"
  ) +
  stat_cor(method = "spearman", label.x = 3, label.y = 6, size = 3)

# Combine the plots
combined_plot <- B | C | A

# Display combined plot
combined_plot

# Define save path (ensure directory exists)
#save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/plots/"
#ggsave(paste0(save_path, "377_vs_393.tiff"), 
      # plot = combined_plot, 
      # device = "tiff", 
      # dpi = 800, 
       #width = 10, 
       #height = 4, 
      # units = "in")
