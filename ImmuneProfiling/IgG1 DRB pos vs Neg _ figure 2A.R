# Required libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)

# Load the data
COMPILED_DATA_EBNA1b <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
                                   sheet = "Data_2024")

# List of DRB1_1501_Positive participants
DRB1_1501_Positive_participants <- c("1580", "1586", "1611", "1620", "1652", "1654", "1665", "1555", "1559", "1602", "1617", "1624", "1645", "1653", "1664", "1666", "1679", "1681", "1684", "ES346", "ES357", "ES389", "ES406", "ES411", "ES418", "ES459", "ES372", "ES396", "ES452", "ES488")

# Define the list of variables to analyze
variables <- c("IgG1_EBNA1_421_476_B958", "IgG1_EBNA1_589_641_B958", "IgG1_EBNA1_589_641_AG876", 
               "IgG1_EBNA1_1_56_GD1", "IgG1_EBNA1_1_56_B958", "IgG1_EBNA1_421_476_AG876", 
               "IgG1_EBNA1_449_504_B958", "IgG1_EBNA1_1_90", "IgG1_EBNA1_377_459", 
               "IgG1_EBNA1_460_641", "IgG1_EBNA1_29_84_GD1", "IgG1_EBNA1_365_420_B958", 
               "IgG1_EBNA1_393_448_B958", "IgG1_EBNA1_393_448_AG876", "IgG1_GP350_R", 
               "IgG1_VCA_p23", "IgG1_Early_Ag", "IgG1_GH", "IgG1_GP350_V")

custom_labels <- c(
  "IgG1_EBNA1_421_476_B958" = "aa 421-476_B958",
  "IgG1_EBNA1_589_641_B958" = "aa 589-641_B958",
  "IgG1_EBNA1_589_641_AG876" = "aa 589-641_AG876",
  "IgG1_EBNA1_1_56_GD1" = "aa 1-56_GD1",
  "IgG1_EBNA1_1_56_B958" = "aa 1-56_B958",
  "IgG1_EBNA1_421_476_AG876" = "aa 421-476_AG876",
  "IgG1_EBNA1_449_504_B958" = "aa 449-504_B958",
  "IgG1_EBNA1_1_90" = "aa 1-90",
  "IgG1_EBNA1_377_459" = "aa 377-459",
  "IgG1_EBNA1_460_641" = "aa 460-641",
  "IgG1_EBNA1_29_84_GD1" = "aa 29-84_GD1",
  "IgG1_EBNA1_365_420_B958" = "aa 365-420",
  "IgG1_EBNA1_393_448_B958" = "aa 393-448",
  "IgG1_EBNA1_393_448_AG876" = "aa 393-448_AG876",
  "IgG1_GP350_R" = "GP350_R",
  "IgG1_VCA_p23" = "VCA_p23",
  "IgG1_Early_Ag" = "Early_Ag",
  "IgG1_GH" = "GH",
  "IgG1_GP350_V" = "GP350_V"
)

# Add DRB1_1501_Positive_status to the data
COMPILED_DATA_EBNA1b <- COMPILED_DATA_EBNA1b %>%
  mutate(DRB1_1501_Positive_status = case_when(
    participants %in% DRB1_1501_Positive_participants ~ "DRB1_1501_Positive",
    EBV_infection == "Sero_Negative" ~ "Sero_Negative",
    EBV_infection == "PBS" ~ "PBS",
    TRUE ~ "DRB1_1501_Negative"
  ))

# Combine and update the 'EBV_infection' column for the new groups
COMPILED_DATA_EBNA1b <- COMPILED_DATA_EBNA1b %>%
  mutate(EBV_infection = case_when(
    EBV_infection == "6 months (n = 30)" ~ "6 months (n = 30)",
    EBV_infection == "1 year (n = 67)" ~ "1 year (n = 67)",
    EBV_infection == "Sero_Pos (n = 30)" ~ "Sero_Pos (n = 50)",
    EBV_infection == "Sero_Pos (n = 20)_HX" ~ "Sero_Pos (n = 50)",
    TRUE ~ EBV_infection  # Leave the rest unchanged
  )) %>%
  # Remove unwanted groups in a single step
  filter(!EBV_infection %in% c("6 weeks (n = 67)", "Sero_Negative", "SC (n = 21)", "NA")) %>%
  filter(!`SAMPLE ID` %in% c("LLOD", "Average Sero_Negative", "STD dev")) 

# Ensure the correct order of EBV_infection factor levels
COMPILED_DATA_EBNA1b$EBV_infection <- factor(COMPILED_DATA_EBNA1b$EBV_infection, 
                                             levels = c("Acute (n = 97)", "6 months (n = 30)", "1 year (n = 67)", "Sero_Pos (n = 50)"))

# Define the order and colors for the plots
colorlist <- c("DRB1_1501_Positive" = "dodgerblue3", 
               "DRB1_1501_Negative" = "red", 
               "Sero_Negative" = "grey")
border_colorlist <- c("DRB1_1501_Positive" = "dodgerblue4", 
                      "DRB1_1501_Negative" = "darkred", 
                      "Sero_Negative" = "darkgrey")

# Function to generate violin plots with significance stars using precomputed p-values
generate_plot_with_significance <- function(variable_name) {
  # Extract p-values for the variable
  p_values <- mann_whitney_results_df %>% filter(Variable == variable_name)
  
  # Map p-values to significance stars
  acute_signif <- p_value_to_significance(p_values$Acute_pvalue)
  months6_signif <- p_value_to_significance(p_values$Months6_pvalue)
  year1_signif <- p_value_to_significance(p_values$Year1_pvalue)
  
  # Generate the plot
  plot <- ggplot(COMPILED_DATA_EBNA1b, aes(x = EBV_infection, y = !!sym(variable_name), fill = DRB1_1501_Positive_status)) +
    geom_violin(scale = "width", adjust = 1.5, color = "black", width = 0.9, position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(color = DRB1_1501_Positive_status), 
                 width = 0.2, show.legend = FALSE, outlier.shape = NA, 
                 fill = "white", position = position_dodge(width = 0.8)) +
    scale_y_continuous(trans = "log10", limits = c(1e3, 1e7)) +
    scale_fill_manual(values = colorlist) + 
    scale_color_manual(values = border_colorlist) + 
    theme_classic() +
    labs(
      title = custom_labels[[variable_name]], 
      x = 'EBV infection status', 
      y = 'Mean Fluorescence Units'
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_text(margin = margin(t = 15))) +
    annotate("text", x = 1, y = 1e6, label = acute_signif, size = 5, color = "black", fontface = "bold") +
    annotate("text", x = 2, y = 1e6, label = months6_signif, size = 5, color = "black", fontface = "bold") +
    annotate ("text", x = 3, y = 1e6, label = year1_signif, size = 5, color = "black", fontface = "bold")
  
  return(plot)
}

# Function to perform Mann-Whitney U test (Wilcoxon rank-sum test)
perform_mann_whitney <- function(data, time_point, variable) {
  # Filter data for the given time point
  subset_data <- data %>%
    filter(EBV_infection == time_point) %>%
    select(DRB1_1501_Positive_status, !!sym(variable)) %>%
    filter(!is.na(!!sym(variable)))  # Remove NA values
  
  # Perform Mann-Whitney U test
  mw_test_result <- wilcox.test(
    formula = as.formula(paste0(variable, " ~ DRB1_1501_Positive_status")),
    data = subset_data
  )
  return(mw_test_result$p.value)
}

# Helper function to convert p-values to significance stars
p_value_to_significance <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Apply the Mann-Whitney U test for each variable for all three time points
mann_whitney_results <- lapply(variables, function(var) {
  acute_pvalue <- perform_mann_whitney(COMPILED_DATA_EBNA1b, "Acute (n = 97)", var)
  months6_pvalue <- perform_mann_whitney(COMPILED_DATA_EBNA1b, "6 months (n = 30)", var)
  year1_pvalue <- perform_mann_whitney(COMPILED_DATA_EBNA1b, "1 year (n = 67)", var)
  data.frame(Variable = var, Acute_pvalue = acute_pvalue, Months6_pvalue = months6_pvalue, Year1_pvalue = year1_pvalue)
})

# Combine the results into a data frame
mann_whitney_results_df <- do.call(rbind, mann_whitney_results)

# Generate individual plots for all variables; legend only for the last plot
plot_list <- lapply(seq_along(variables), function(i) {
  generate_plot_with_significance(variables[i]) + 
    theme(legend.position = ifelse(i == length(variables), "right", "none"))
})

# Specify the indices of the plots you want to combine
plot9 <- plot_list[[12]]
plot12 <- plot_list[[9]]
plot13 <- plot_list[[13]]  # This plot will have the legend

# Combine the plots 9, 12, and 13 into one plot with only the legend in plot13
combined_plot <- (plot12 | plot9 | plot13)

# Print the combined plot
print(combined_plot)
# Combine the results into a data frame
mann_whitney_results_df <- do.call(rbind, mann_whitney_results)

# Print the results
print(mann_whitney_results_df)


# Define the path to save the plot
save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/manuscript/Manuscript figures/"

# Save the combined plot
ggsave(paste0(save_path, "IgG1_DRB_Ebna1_no_threshold.tiff"), 
       plot = combined_plot, 
       device = "tiff", 
       dpi = 450, 
       width = 8, 
       height = 4, 
       units = "in")
dev.off()  # Close the device

  