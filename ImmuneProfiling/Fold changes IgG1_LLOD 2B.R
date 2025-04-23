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

# Load the data
COMPILED_DATA_EBNA1b <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
                                   sheet = "Data_2024")

# List of DRB1_1501_Positive participants
DRB1_1501_Positive_participants <- c("1580", "1586", "1611", "1620", "1652", "1654", "1665", "1555", 
                                     "1559", "1602", "1617", "1624", "1645", "1653", "1664", "1666", 
                                     "1679", "1681", "1684", "ES346", "ES357", "ES389", "ES406", 
                                     "ES411", "ES418", "ES459", "ES372", "ES396", "ES452", "ES488")

# Filter the data for Acute and 1-year time points
acute_data <- COMPILED_DATA_EBNA1b %>% filter(EBV_infection == "Acute (n = 97)")
year1_data <- COMPILED_DATA_EBNA1b %>% filter(EBV_infection == "1 year (n = 67)")

# Merge Acute and 1-year data based on participant IDs
merged_data <- merge(acute_data, year1_data, by = "participants", suffixes = c("_acute", "_year1"))

# Add DRB1_1501_Positive_status column based on participant ID
merged_data <- merged_data %>%
  mutate(DRB1_1501_Positive_status_acute = ifelse(participants %in% DRB1_1501_Positive_participants, 
                                                  "DRB1_1501_Positive", "DRB1_1501_Negative"))

# List of variables to analyze
variables <- c("IgG1_EBNA1_377_459", "IgG1_EBNA1_365_420_B958", "IgG1_EBNA1_393_448_B958")

# Extract LLOD values from the 387th row
LLOD_values <- COMPILED_DATA_EBNA1b[387, variables]

# Function to calculate individual fold changes for each participant with LLOD applied
calculate_individual_fold_change <- function(row, variable, LLOD) {
  value_acute <- row[[paste0(variable, "_acute")]]
  value_year1 <- row[[paste0(variable, "_year1")]]
  
  # Apply the LLOD threshold for acute values, set value to LLOD if below the threshold
  if (!is.na(LLOD) && value_acute < LLOD) {
    value_acute <- LLOD  # Set to the threshold (LLOD) if below the threshold
  }
  
  # Avoid division by zero and missing values
  if (is.na(value_acute) | value_acute == 0) {
    return(NA)
  } else {
    fold_change <- value_year1 / value_acute
    return(fold_change)
  }
}

# Calculate fold changes for each participant for each variable using the suffixed columns
fold_changes <- merged_data %>%
  rowwise() %>%
  mutate(across(all_of(paste0(variables, "_acute")), 
                ~ calculate_individual_fold_change(cur_data(), str_remove(cur_column(), "_acute"), 
                                                   LLOD_values[[str_remove(cur_column(), "_acute")]]), 
                .names = "FC_{col}"))

# Reshape the fold change data for statistical testing and table output
fold_changes_long <- fold_changes %>%
  pivot_longer(cols = starts_with("FC_"), names_to = "variable", values_to = "fold_change") %>%
  mutate(variable = str_replace(variable, "FC_", ""))

# Function to perform Mann-Whitney U test
perform_mann_whitney_fold_change <- function(data, var) {
  subset_data <- data %>% filter(variable == var)
  
  # Ensure there are at least two groups to compare
  if (length(unique(subset_data$DRB1_1501_Positive_status_acute)) < 2) {
    return(NA)
  }
  
  mw_test <- wilcox.test(fold_change ~ DRB1_1501_Positive_status_acute, data = subset_data)
  return(mw_test$p.value)
}

# Calculate p-values for each variable individually
fold_changes_long <- fold_changes_long %>%
  group_by(variable) %>%
  mutate(p_value = perform_mann_whitney_fold_change(., variable))

# Apply FDR correction on p-values
fold_changes_long <- fold_changes_long %>%
  ungroup() %>%
  mutate(p_value_fdr = p.adjust(p_value, method = "BH"),
         Significance = ifelse(p_value_fdr < 0.001, "***",
                               ifelse(p_value_fdr < 0.01, "**",
                                      ifelse(p_value_fdr < 0.05, "*", "ns"))))

# Output the results as a table of fold changes and p-values before and after FDR correction
stats_with_fold_changes <- fold_changes_long %>%
  select(participants, variable, fold_change, p_value, p_value_fdr, Significance, DRB1_1501_Positive_status_acute)

# Display the p-values before and after FDR correction
print(stats_with_fold_changes)

# Save the results (including fold changes and p-values) to a CSV
write.csv(stats_with_fold_changes, "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/fold_changes_and_pvalues.csv")

# Plot using the pre-calculated p-values

# Your custom labels for facet titles
custom_labels <- c(
  "IgG1_EBNA1_377_459_acute" = "aa 377-459",
  "IgG1_EBNA1_365_420_B958_acute" = "aa 365-420",
  "IgG1_EBNA1_393_448_B958_acute" = "aa 393-448"
)

# Summarize the data to get one entry per variable for p-values and significance
fold_changes_summary <- fold_changes_long %>%
  group_by(variable) %>%
  dplyr::summarise(
    max_fold_change = max(fold_change, na.rm = TRUE),
    p_value_fdr = first(p_value_fdr),
    Significance = first(Significance)
  )

# Plot the fold changes
p <- ggplot(fold_changes_long, aes(x = DRB1_1501_Positive_status_acute, y = fold_change, fill = DRB1_1501_Positive_status_acute)) +
  geom_boxplot(outlier.shape = NA) +    # Boxplot without outliers
  geom_jitter(width = 0.2, size = 1) +  # Add jittered points to avoid overplotting
  
  # Use facet_wrap with labeller to apply custom labels
  facet_wrap(~ variable, scales = "free", nrow = 1, 
             labeller = labeller(variable = custom_labels)) +
  
  # Add the p-values and significance stars using the summarized data
  geom_text(data = fold_changes_summary, 
            aes(x = 1.5, y = max_fold_change * 1.1, label = Significance), 
            size = 5, color = "black", fontface = "bold", inherit.aes = FALSE) +  # Turn off inheriting aesthetics from ggplot
  geom_text(data = fold_changes_summary, 
            aes(x = 1.5, y = max_fold_change * 1.05, label = paste0("P = ", round(p_value_fdr, 3))), 
            size = 3, color = "black", fontface = "bold", inherit.aes = FALSE) +  # Turn off inheriting aesthetics from ggplot
  
  # Apply log scale to y-axis and customize limits
  scale_y_log10(limits = c(1, max(fold_changes_long$fold_change, na.rm = TRUE) * 1.2)) +
  
  # Customize x-axis labels and fill colors
  scale_x_discrete(labels = c("Negative", "Positive")) +
  scale_fill_manual(values = c("DRB1_1501_Positive" = "dodgerblue3", "DRB1_1501_Negative" = "red")) +
  
  # Customize the overall plot appearance
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "plain"),
        legend.position = "none")  # This removes the legend

# Print the plot
print(p)
# Define the path where you want to save the plot
save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/manuscript/Manuscript figures/final manuscript/PDF figures"

# Use ggsave to save your plot as PDF with high resolution
ggsave(paste0(save_path, "2b modified.pdf"), 
       plot = p, 
       device = "pdf", 
       dpi = 1000,   # DPI doesn't have a significant effect for vector formats like PDF, but it's good to include
       width = 8, 
       height = 4, 
       units = "in")
dev.off()  # Close the device