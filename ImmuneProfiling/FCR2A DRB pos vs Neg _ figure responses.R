# Required libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)

# Define the p_value_to_significance function
p_value_to_significance <- function(p_value) {
  if (is.na(p_value)) {
    return("")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Load the data
COMPILED_DATA_EBNA1b <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
                                   sheet = "Data_2024")

# List of DRB1_1501_Positive participants
DRB1_1501_Positive_participants <- c("1580", "1586", "1611", "1620", "1652", "1654", "1665", "1555", "1559", "1602", "1617", "1624", "1645", "1653", "1664", "1666", "1679", "1681", "1684", "ES346", "ES357", "ES389", "ES406", "ES411", "ES418", "ES459", "ES372", "ES396", "ES452", "ES488")

# Define the list of variables to analyze
variables <- c("FcγRIIa_EBNA1_421_476_B958", "FcγRIIa_EBNA1_589_641_B958", "FcγRIIa_EBNA1_589_641_AG876", 
               "FcγRIIa_EBNA1_1_56_GD1", "FcγRIIa_EBNA1_1_56_B958", "FcγRIIa_EBNA1_421_476_AG876", 
               "FcγRIIa_EBNA1_449_504_B958", "FcγRIIa_EBNA1_1_90", "FcγRIIa_EBNA1_377_459", 
               "FcγRIIa_EBNA1_460_641", "FcγRIIa_EBNA1_29_84_GD1", "FcγRIIa_EBNA1_365_420_B958", 
               "FcγRIIa_EBNA1_393_448_B958", "FcγRIIa_EBNA1_393_448_AG876", "FcγRIIa_GP350_R", 
               "FcγRIIa_VCA_p23", "FcγRIIa_Early_Ag", "FcγRIIa_GH", "FcγRIIa_GP350_V")

custom_labels <- c(
  "FcγRIIa_EBNA1_421_476_B958" = "aa 421-476_B958",
  "FcγRIIa_EBNA1_589_641_B958" = "aa 589-641_B958",
  "FcγRIIa_EBNA1_589_641_AG876" = "aa 589-641_AG876",
  "FcγRIIa_EBNA1_1_56_GD1" = "aa 1-56_GD1",
  "FcγRIIa_EBNA1_1_56_B958" = "aa 1-56_B958",
  "FcγRIIa_EBNA1_421_476_AG876" = "aa 421-476_AG876",
  "FcγRIIa_EBNA1_449_504_B958" = "aa 449-504_B958",
  "FcγRIIa_EBNA1_1_90" = "aa 1-90",
  "FcγRIIa_EBNA1_377_459" = "aa 377-459",
  "FcγRIIa_EBNA1_460_641" = "aa 460-641",
  "FcγRIIa_EBNA1_29_84_GD1" = "aa 29-84_GD1",
  "FcγRIIa_EBNA1_365_420_B958" = "aa 365-420",
  "FcγRIIa_EBNA1_393_448_B958" = "aa 393-448",
  "FcγRIIa_EBNA1_393_448_AG876" = "aa 393-448_AG876",
  "FcγRIIa_GP350_R" = "GP350_R",
  "FcγRIIa_VCA_p23" = "VCA_p23",
  "FcγRIIa_Early_Ag" = "Early_Ag",
  "FcγRIIa_GH" = "GH",
  "FcγRIIa_GP350_V" = "GP350_V"
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
# Function to apply FDR correction
apply_fdr_correction <- function(p_values) {
  p.adjust(p_values, method = "BH")  # Benjamini-Hochberg correction
}

# Function to perform Mann-Whitney U test and return FDR-corrected p-values
perform_mann_whitney_with_fdr <- function(data, variables, time_points) {
  results <- list()
  
  for (time_point in time_points) {
    p_values <- sapply(variables, function(var) {
      subset_data <- data %>%
        filter(EBV_infection == time_point) %>%
        select(DRB1_1501_Positive_status, !!sym(var)) %>%
        filter(!is.na(!!sym(var)))  # Remove NA values
      
      if (nrow(subset_data) > 1) {
        test_result <- wilcox.test(
          formula = as.formula(paste0(var, " ~ DRB1_1501_Positive_status")),
          data = subset_data
        )
        return(test_result$p.value)
      } else {
        return(NA)  # If insufficient data, return NA
      }
    })
    
    # Apply FDR correction to p-values
    fdr_corrected <- apply_fdr_correction(p_values)
    results[[time_point]] <- data.frame(
      Variable = variables,
      Raw_pvalue = p_values,
      FDR_pvalue = fdr_corrected
    )
  }
  
  return(results)
}

# Define the time points to test
time_points <- c("Acute (n = 97)", "6 months (n = 30)", "1 year (n = 67)")

# Perform Mann-Whitney U test with FDR correction for all variables and time points
mann_whitney_fdr_results <- perform_mann_whitney_with_fdr(COMPILED_DATA_EBNA1b, variables, time_points)

# Combine results into a single data frame for easier access
mann_whitney_fdr_combined <- bind_rows(lapply(names(mann_whitney_fdr_results), function(tp) {
  data <- mann_whitney_fdr_results[[tp]]
  data$Time_point <- tp
  return(data)
}))

# Print the combined results
print(mann_whitney_fdr_combined)

# Generate violin plots with FDR-corrected significance stars
generate_plot_with_fdr <- function(variable_name) {
  # Extract FDR-corrected p-values for the variable
  fdr_data <- mann_whitney_fdr_combined %>%
    filter(Variable == variable_name)
  
  # Prepare significance stars based on FDR p-values
  acute_signif <- p_value_to_significance(fdr_data %>% filter(Time_point == "Acute (n = 97)") %>% pull(FDR_pvalue))
  months6_signif <- p_value_to_significance(fdr_data %>% filter(Time_point == "6 months (n = 30)") %>% pull(FDR_pvalue))
  year1_signif <- p_value_to_significance(fdr_data %>% filter(Time_point == "1 year (n = 67)") %>% pull(FDR_pvalue))
  
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
          axis.title.x = element_text(margin = margin(t = 15)))
  
  # Add FDR-corrected significance stars
  plot <- plot +
    annotate("text", x = 1, y = 1e7, label = acute_signif, size = 5, color = "black", fontface = "bold") +
    annotate("text", x = 2, y = 1e7, label = months6_signif, size = 5, color = "black", fontface = "bold") +
    annotate("text", x = 3, y = 1e7, label = year1_signif, size = 5, color = "black", fontface = "bold")+
    annotate("text", x = 4, y = 1e7, label = year1_signif, size = 5, color = "black", fontface = "bold")
  return(plot)
}

# Generate plots for all variables using FDR correction
plot_list <- lapply(seq_along(variables), function(i) {
  generate_plot_with_fdr(variables[i]) + 
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
# Open the combined results in the RStudio viewer
View(mann_whitney_fdr_combined)

# Save the plot
save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/manuscript/Manuscript figures/"
ggsave(paste0(save_path, "FcγRIIa_DRB_Ebna1_with_FDR.tiff"), 
       plot = combined_plot, 
       device = "tiff", 
       dpi = 450, 
       width = 8, 
       height = 4, 
       units = "in")
dev.off()  # Close the device