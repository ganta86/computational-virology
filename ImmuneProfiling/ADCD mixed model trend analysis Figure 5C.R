# Load necessary libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)  # Added for read_excel function

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
library(readxl)
library(stringr)
library(tidyverse)
library("ggplot2")
library("ggpubr")
library("ggsci")
library(patchwork)
library(scales)
library(rstatix)
library(stringr)
library(writexl)
library(pheatmap)
library("viridis")  
library(vip)
library(Hmisc)
library(extrafont)
library(ggthemes)
library(RColorBrewer)
library(readxl)
library(ggbeeswarm)
# Required libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(ez)  # For repeated measures ANOVA
library(lme4)  # For mixed-effects model
library(lmerTest) 
library(emmeans)
# Load the data
COMPILED_DATA_EBNA1b <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
                                   sheet = "Data_2024")

# List of DRB1_1501_Positive participants
DRB1_1501_Positive_participants <- c("1580", "1586", "1611", "1620", "1652", "1654", "1665", "1555", "1559", "1602", "1617", "1624", "1645", "1653", "1664", "1666", "1679", "1681", "1684", "ES346", "ES357", "ES389", "ES406", "ES411", "ES418", "ES459", "ES372", "ES396", "ES452", "ES488")

# Add DRB1_1501_Positive_status to the data
COMPILED_DATA_EBNA1b <- COMPILED_DATA_EBNA1b %>%
  mutate(DRB1_1501_Positive_status = case_when(
    participants %in% DRB1_1501_Positive_participants ~ "DRB1_1501_Positive",
    EBV_infection == "Sero_Negative" ~ "Sero_Negative",
    EBV_infection == "PBS" ~ "PBS",
    TRUE ~ "DRB1_1501_Negative"
  ))

# Filter the data to remove unwanted groups
COMPILED_DATA_EBNA1b <- COMPILED_DATA_EBNA1b %>%
  filter(!`SAMPLE ID` %in% c("LLOD", "Average Negative", "STD dev","NA")) %>%
  filter(!EBV_infection %in% c("NA", "SC (n = 21)"))

# Define dodge width for consistent plots
dodge_width <- 1

# Define the order and colors for the plots
order <- c("Acute (n = 97)", "6 weeks (n = 67)", "6 months (n = 30)", "1 year (n = 67)", "Sero_Pos (n = 30)","Sero_Pos (n = 20)_HX", "Sero_Negative")
colorlist <- c("DRB1_1501_Positive" = "dodgerblue3", "DRB1_1501_Negative" = "red", "Sero_Negative" = "black")

# Set the correct order of EBV_infection factor levels
COMPILED_DATA_EBNA1b$EBV_infection <- factor(COMPILED_DATA_EBNA1b$EBV_infection, levels = order)

# List of variables to plot
variables <- c( "C3_EBNA1_377_459", 
                "C3_EBNA1_365_420_B958", 
                "C3_EBNA1_393_448_B958", "C3_EBNA1_393_448_AG876")

COMPILED_DATA_EBNA1b_filtered <- COMPILED_DATA_EBNA1b %>%
  filter(!EBV_infection %in% c("Sero_Negative", "Sero_Pos (n = 30)", "Sero_Pos (n = 20)_HX"))

# Load necessary libraries
library(ggplot2)
library(patchwork)  # For combining plots in a row

# Define the four variables to plot
variables <- c("C3_EBNA1_377_459", 
               "C3_EBNA1_365_420_B958", 
               "C3_EBNA1_393_448_B958", 
               "C3_EBNA1_393_448_AG876")

# Function to generate trend plots with straight lines (no smoothing)
# Function to generate trend plots with straight lines (no smoothing) and uniform y-axis
# Function to generate trend plots with straight lines (no smoothing) and remove NA values
generate_trend_plot <- function(variable_name, show_legend = FALSE) {
  plot <- COMPILED_DATA_EBNA1b_filtered %>%
    filter(!is.na(EBV_infection), !is.na(DRB1_1501_Positive_status)) %>%  # Remove rows with NA in key variables
    group_by(EBV_infection, DRB1_1501_Positive_status) %>%
    summarise(mean_value = mean(get(variable_name), na.rm = TRUE),
              se_value = sd(get(variable_name), na.rm = TRUE) / sqrt(n())) %>%
    ggplot(aes(x = EBV_infection, y = mean_value, group = DRB1_1501_Positive_status, 
               color = DRB1_1501_Positive_status)) +
    
    # Add straight lines to connect points
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2) +
    
    # Set uniform y-axis limits across all plots
    scale_y_continuous(trans = "log10", limits = c(5e3, 5e5)) +
    scale_color_manual(values = colorlist) +
    theme_classic() +
    labs(title = variable_name, x = "EBV Infection Stage", y = 'Mean Fluorescence Units (log10)') +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      plot.title = element_text(size = 9),  # Adjust title size
      plot.margin = margin(20, 20, 20, 20)  # Increase margins to avoid clipping
    )
  
  # Remove legend if not required
  if (!show_legend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}

# Generate the three plots
plot1 <- generate_trend_plot("C3_EBNA1_377_459", show_legend = FALSE)
plot2 <- generate_trend_plot("C3_EBNA1_365_420_B958", show_legend = FALSE)
plot3 <- generate_trend_plot("C3_EBNA1_393_448_B958", show_legend = TRUE)

# Combine the three plots in a single row using patchwork
combined_plot <- plot1 | plot2 | plot3

# Adjust the overall output size
combined_plot <- combined_plot + plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Display the combined plot
print(combined_plot)

save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/manuscript/Manuscript figures/final manuscript/PDF figures/"

ggsave(paste0(save_path, "C3trend_plots.pdf"), 
       plot = combined_plot, 
       device = "pdf", 
       dpi = 1000, 
       width = 12, 
       height = 4, 
       units = "in")
      dev.off()  # Close the device
# ==============================================================
# Mixed Model Test Results (Only Mixed-Effects Model)
# ==============================================================

# Ensure the required libraries are loaded
# Ensure the required libraries are loaded
library(lme4)
library(lmerTest)
library(emmeans)  # For post-estimation contrasts

# Loop through each variable for mixed-effects model
for (variable in variables) {
  
  # Create a formula dynamically for the mixed model
  formula_mixed <- as.formula(paste(variable, "~ EBV_infection * DRB1_1501_Positive_status + (1 | participants)", sep=""))
  
  # Run the mixed model
  mixed_model <- lmer(formula_mixed, data = COMPILED_DATA_EBNA1b_filtered)
  
  # Print the summary of the mixed model for each variable
  print(paste("Mixed Model Results for:", variable))
  print(summary(mixed_model))
  
  # Perform post-estimation contrasts using emmeans
  emmeans_results <- emmeans(mixed_model, ~ EBV_infection | DRB1_1501_Positive_status)
  
  # Post-hoc comparisons between EBV_infection stages within each DRB1_1501_Positive_status group
  contrasts <- contrast(emmeans_results, method = "pairwise", adjust = "tukey")
  
  # Print the contrast results for each variable
  print(paste("Post-estimation contrasts for:", variable))
  print(contrasts)
}
