# Load necessary libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(reshape2)
library(patchwork)
library(rstatix)
library(knitr)

# Load your data

all_blocking_hi_cryab <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/cross reactive antigens luminex/pure cross antigens/blocking/blocking by beads hi cryab/all blocking hi cryab.xlsx", 
                                    sheet = "DATA")
all_blocking_hi_cryab <- as.data.frame(all_blocking_hi_cryab)

# Define variable pairs for block treated and untreated
variable_pairs <- list(
  list(block_vars = "CRYAB_EBNAspiked", untreated_vars = "CRYAB_FLU_E6_spiked", label = "Anti-CRYAB IgG reactivity"),
  list(block_vars = "FLU_EBNAspiked", untreated_vars = "FLU_FLU_E6_spiked", label = "Anti-HA flu IgG reactivity"),
  list(block_vars = "VCA_EBNAspiked", untreated_vars = "VCA_FLU_E6_spiked", label = "Anti-VCA IgG reactivity")
)

# Filter the dataset to include all groups (Convalescence, Sero_Pos, Negative, PBS if required)
all_blocking_hi_cryab <- all_blocking_hi_cryab %>%
  filter(EBV_infection %in% c("Convalescence (n=55)", "Sero_Pos (n = 25)", "Negative", "PBS"))

# Ensure the correct order of EBV_infection factor levels
all_blocking_hi_cryab$EBV_infection <- factor(
  all_blocking_hi_cryab$EBV_infection, 
  levels = c("Convalescence (n=55)", "Sero_Pos (n = 25)", "Negative", "PBS")
)

# Function to calculate percentage reduction with proper handling for smaller untreated values
calculate_percentage_reduction <- function(blocked_vals, untreated_vals) {
  mean_blocked <- mean(blocked_vals, na.rm = TRUE)
  mean_untreated <- mean(untreated_vals, na.rm = TRUE)
  
  if (mean_untreated == 0) {
    return("Undefined")  # Prevent division by zero
  } else if (mean_blocked >= mean_untreated) {
    reduction <- ((mean_blocked - mean_untreated) / mean_blocked) * 100
  } else {
    reduction <- ((mean_untreated - mean_blocked) / mean_untreated) * 100
  }
  
  return(sprintf("%.1f%%", reduction))
}

# Initialize an empty data frame to store statistical results
stat_results <- data.frame(
  Variable = character(),
  Group = character(),
  P_value = numeric(),
  Significance = character(),
  Reduction = character(),
  stringsAsFactors = FALSE
)

# Function to perform paired tests within each group and create violin plots
create_group_comparison_plot <- function(data, block_var, untreated_var, title) {
  
  # Convert data to long format
  data_long <- melt(data, measure.vars = c(block_var, untreated_var), variable.name = "variable", value.name = "value")
  
  # Add treatment information and convert the values to numeric
  data_long <- data_long %>%
    mutate(treatment = ifelse(variable == block_var, "EBNA1 peptides", "HA-flu peptides"),
           value = as.numeric(value)) %>%
    filter(!is.na(value))  # Remove NA values
  
  # Group the data by infection status (Convalescence, Sero_Pos, Negative, PBS)
  for (group in levels(data_long$EBV_infection)) {
    
    group_data <- data_long %>% filter(EBV_infection == group)
    
    # Perform paired Wilcoxon test, use exact = FALSE to handle zeros
    blocked_vals <- group_data$value[group_data$treatment == "EBNA1 peptides"]
    untreated_vals <- group_data$value[group_data$treatment == "HA-flu peptides"]
    
    wilcox_test_result <- wilcox.test(blocked_vals, untreated_vals, paired = TRUE, exact = FALSE)
    p_value <- wilcox_test_result$p.value
    significance <- case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
    
    # Calculate percentage reduction
    reduction <- calculate_percentage_reduction(blocked_vals, untreated_vals)
    
    # Store the results in the stat_results data frame
    stat_results <<- rbind(stat_results, data.frame(
      Variable = title,
      Group = group,
      P_value = p_value,
      Significance = significance,
      Reduction = reduction
    ))
  }
  
  # Create the violin plot with color fill based on treatment (EBNA1 vs HA-flu)
  plot <- ggplot(data_long, aes(x = EBV_infection, y = value, fill = treatment)) +
    geom_violin(size = 0.5, alpha = 0.5, position = position_dodge(width = 1), 
                width = 1, scale = "width", adjust = 1) +
    geom_boxplot(aes(color = treatment), width = 0.2, outlier.shape = NA, fill = "white", 
                 position = position_dodge(width = 1), show.legend = FALSE) +
    scale_y_continuous(trans = "log10", limits = c(10, 100000)) +
    theme_classic() +
    labs(title = title, x = 'EBV infection status', y = 'Mean Fluorescence Units', fill = 'Treatment') +  
    theme(aspect.ratio = 0.75, axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
    scale_fill_manual(values = c("EBNA1 peptides" = "violetred4", "HA-flu peptides" = "grey30")) +
    # Add significance stars for Convalescence, Sero_Pos, Negative, and PBS
    geom_text(aes(x = "Convalescence (n=55)", y = 1e5, 
                  label = stat_results$Significance[stat_results$Group == "Convalescence (n=55)" & stat_results$Variable == title]), 
              color = "black", size = 4) +
    geom_text(aes(x = "Sero_Pos (n = 25)", y = 1e5, 
                  label = stat_results$Significance[stat_results$Group == "Sero_Pos (n = 25)" & stat_results$Variable == title]), 
              color = "black", size = 4) +
    geom_text(aes(x = "Negative", y = 1e5, 
                  label = stat_results$Significance[stat_results$Group == "Negative" & stat_results$Variable == title]), 
              color = "black", size = 4) +
    geom_text(aes(x = "PBS", y = 1e5, 
                  label = stat_results$Significance[stat_results$Group == "PBS" & stat_results$Variable == title]), 
              color = "black", size = 4)
  
  return(plot)
}

# Loop through the variable pairs and generate the plots with group comparisons
plots <- lapply(variable_pairs, function(vars) {
  create_group_comparison_plot(all_blocking_hi_cryab, vars$block_vars, vars$untreated_vars, vars$label)
})

# Combine all the plots into a single layout (3 columns)
combined_plot <- wrap_plots(plots, ncol = 3)

# Print the combined violin plot
print(combined_plot)

# Display the statistical results as a table
kable(stat_results, caption = "Statistical Comparisons and Percentage Reductions by Group")
