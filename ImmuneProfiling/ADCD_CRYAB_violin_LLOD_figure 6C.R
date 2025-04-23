# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# Load the dataset
COMPILED_DATA_EBNA1b <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
                                   sheet = "Data_2024")

# Ensure the data is a data frame
COMPILED_DATA_EBNA1b <- as.data.frame(COMPILED_DATA_EBNA1b)

# Define common theme for the plots
labssizes <- theme(axis.title = element_text(size = 14, face = 'bold', color = 'black'),
                   axis.text = element_text(size = 12, color = 'black'),
                   axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
                   plot.title = element_text(size = 14, face = 'bold', hjust = 0.5, vjust = 0.5, color = 'black'))

# Define groups and colors
groups <- scale_x_discrete(limits = c("Acute (n = 97)", "6 weeks (n = 67)", "6 months (n = 30)", 
                                      "1 year (n = 67)", "Sero_Pos (n = 50)", "Sero_Negative"))
mycolor <- scale_color_manual(breaks = c("Acute (n = 97)", "6 weeks (n = 67)", "6 months (n = 30)", 
                                         "1 year (n = 67)", "Sero_Pos (n = 50)", "Sero_Negative"), 
                              values = c("magenta4", "dodgerblue3", "chocolate3", "darkgreen", "blue", "black"))
myfill <- scale_fill_manual(breaks = c("Acute (n = 97)", "6 weeks (n = 67)", "6 months (n = 30)", 
                                       "1 year (n = 67)", "Sero_Pos (n = 50)", "Sero_Negative"), 
                            values = c("magenta4", "dodgerblue3", "chocolate3", "darkgreen", "blue", "black"))

# Create a named vector for custom labels (updated with the new variables)
# Create a named vector for custom labels
custom_labels <- c(
  "Total_IgG_ANO2" = "ANO2",
  "Total_IgG_GlialCAM" = "GlialCAM",
  "Total_IgG_CRYAB" = "CRYAB",
  "Total_IgG_ANO2_urea" = "ANO2",
  "Total_IgG_GlialCAM_urea" = "GlialCAM",
  "Total_IgG_CRYAB_urea" = "CRYAB",
  "ADCD_ANO2" = "ANO2",
  "ADCD_GlialCAM" = "GlialCAM",
  "ADCD_CRYAB" = "CRYAB"
)

# List of variables to test
variables_to_test <- c("Total_IgG_ANO2", "Total_IgG_GlialCAM", "Total_IgG_CRYAB", 
                       "Total_IgG_ANO2_urea", "Total_IgG_GlialCAM_urea", "Total_IgG_CRYAB_urea", 
                       "ADCD_ANO2", "ADCD_GlialCAM", "ADCD_CRYAB")

# Extract LLOD values (assuming they're stored in the 387th row)
LLOD_values <- COMPILED_DATA_EBNA1b[387, variables_to_test]

# Function to apply the LLOD threshold
apply_llod_threshold <- function(row, variable, LLOD) {
  value <- row[[variable]]  # Access the value in the specific column
  # Ensure value and LLOD are not NA before comparing
  if (!is.na(value) && !is.na(LLOD) && value < LLOD) {
    value <- LLOD  # Replace with LLOD value if the value is below the threshold
  }
  return(value)
}

# Function to apply LLOD threshold to all acute samples
apply_llod_to_acute <- function(data, variables, LLOD_values) {
  for (variable in variables) {
    LLOD_value <- as.numeric(LLOD_values[[variable]])  # Get the LLOD for the variable
    # Apply threshold only if LLOD is not NA
    if (!is.na(LLOD_value)) {
      data <- data %>%
        rowwise() %>%
        mutate(!!variable := apply_llod_threshold(cur_data(), variable, LLOD_value)) %>%
        ungroup()
    } else {
      message(paste("Skipping variable", variable, "due to missing LLOD value"))
    }
  }
  return(data)
}

# Now apply the LLOD threshold to a copy of the dataset
COMPILED_DATA_EBNA1b_with_LLOD <- apply_llod_to_acute(COMPILED_DATA_EBNA1b, variables_to_test, LLOD_values)


# Initialize a data frame to store all p-values (paired and independent)
all_p_values <- data.frame()

# Define your independent comparisons here as a list of named pairs
independent_comparisons <- list(
  c("1 year (n = 67)", "Sero_Negative"),
  c("Sero_Pos (n = 50)", "Sero_Negative"),
  c("Sero_Pos (n = 50)", "1 year (n = 67)"),
  c("Acute (n = 97)", "Sero_Pos (n = 50)"),
  c("Acute (n = 97)", "Sero_Negative")
)

# Function to perform paired tests, including Acute vs 6 months
perform_paired_tests_with_LLOD <- function(variables_to_test) {
  # Filter relevant data for Acute, 1 year, and 6 months (visit 7)
  paired_data <- COMPILED_DATA_EBNA1b_with_LLOD %>%
    filter((EBV_infection == "Acute (n = 97)" & visit %in% c("v1", "v2")) | 
             (EBV_infection == "1 year (n = 67)" & visit == "v8") |
             (EBV_infection == "6 months (n = 30)" & visit == "v7")) %>%
    pivot_wider(
      id_cols = participants,  # Use participants to identify pairs
      names_from = EBV_infection,
      values_from = variables_to_test
    )
  
  # Iterate over each variable and perform paired Wilcoxon test
  for (variable in variables_to_test) {
    # Define column names for Acute, 1 year, and 6 months
    acute_col <- paste0(variable, "_Acute (n = 97)")
    one_year_col <- paste0(variable, "_1 year (n = 67)")
    six_months_col <- paste0(variable, "_6 months (n = 30)")
    
    # Perform paired test for Acute vs 1 year
    complete_pairs_1year <- paired_data %>%
      filter(!is.na(!!sym(acute_col)) & !is.na(!!sym(one_year_col)))
    
    if (nrow(complete_pairs_1year) == 67) {
      paired_test <- wilcox.test(
        complete_pairs_1year[[acute_col]], 
        complete_pairs_1year[[one_year_col]], 
        paired = TRUE
      )
      all_p_values <<- rbind(all_p_values, data.frame(
        response_var = variable,
        comparison = "Acute_vs_1Year",
        p_value = paired_test$p.value,
        group1 = "Acute (n = 97)",
        group2 = "1 year (n = 67)",
        test_type = "Paired"
      ))
    }
    
    # Perform paired test for Acute vs 6 months (visit 7)
    complete_pairs_6months <- paired_data %>%
      filter(!is.na(!!sym(acute_col)) & !is.na(!!sym(six_months_col)))
    
    if (nrow(complete_pairs_6months) == 30) {
      paired_test_6months <- wilcox.test(
        complete_pairs_6months[[acute_col]], 
        complete_pairs_6months[[six_months_col]], 
        paired = TRUE
      )
      all_p_values <<- rbind(all_p_values, data.frame(
        response_var = variable,
        comparison = "Acute_vs_6Months",
        p_value = paired_test_6months$p.value,
        group1 = "Acute (n = 97)",
        group2 = "6 months (n = 30)",
        test_type = "Paired"
      ))
    }
  }
}


# Function to perform independent tests
perform_independent_tests_with_LLOD <- function(variables_to_test) {
  for (variable in variables_to_test) {
    for (comparison in independent_comparisons) {
      subset_data <- COMPILED_DATA_EBNA1b_with_LLOD %>% filter(EBV_infection %in% comparison)
      independent_test <- wilcox.test(get(variable) ~ EBV_infection, data = subset_data)
      all_p_values <<- rbind(all_p_values, data.frame(
        response_var = variable,
        comparison = paste(comparison, collapse = "_vs_"),
        p_value = independent_test$p.value,
        group1 = comparison[1],
        group2 = comparison[2],
        test_type = "Independent"
      ))
    }
  }
}

# Perform paired and independent tests
perform_paired_tests_with_LLOD(variables_to_test)
perform_independent_tests_with_LLOD(variables_to_test)

# Apply FDR correction for all p-values
all_p_values$fdr_corrected <- p.adjust(all_p_values$p_value, method = "fdr")
all_p_values$significance <- ifelse(all_p_values$fdr_corrected < 0.0001, "****",
                                    ifelse(all_p_values$fdr_corrected < 0.001, "***",
                                           ifelse(all_p_values$fdr_corrected < 0.01, "**",
                                                  ifelse(all_p_values$fdr_corrected < 0.05, "*", "ns"))))

# Filter out non-significant results and those with a single star ("*")
significant_p_values <- all_p_values %>%
  filter(significance != "ns" & significance != "*")

# Function to create plots
create_plot <- function(df, response_var, group_var, stats_results) {
  ggplot(df, aes_string(x = group_var, y = response_var)) +
    geom_violin(aes(fill = !!sym(group_var)), color = NA, alpha = 0.6, scale = "width", adjust = 1.5, show.legend = FALSE) +
    geom_boxplot(aes(color = !!sym(group_var)), width = 0.2, outlier.shape = NA, show.legend = FALSE) +
    scale_y_continuous(trans = "log10") +
    coord_cartesian(ylim = c(1e3, 1e8)) +
    theme_classic() +
    labs(
      title = custom_labels[[response_var]],
      x = "EBV infection status",
      y = 'C3 Mean Fluorescence Units'
    ) +
    labssizes + groups + mycolor + myfill +
    # Only add p-value annotations if there are significant results
    if (nrow(stats_results) > 0) stat_pvalue_manual(data = stats_results, label = "significance", y.position = 6.5, step.increase = 0.05, size = 6) else NULL +
    theme(legend.position = "none")
}

# Define common parameters for plots
group_var <- "EBV_infection"
response_vars <- variables_to_test

# Loop through response variables to create plots
plots <- lapply(response_vars, function(response_var) {
  stats_results <- significant_p_values %>% filter(response_var == !!response_var)
  create_plot(COMPILED_DATA_EBNA1b, response_var, group_var, stats_results)
})

# Print original p-values and FDR-corrected p-values
print(all_p_values[, c("response_var", "comparison", "p_value", "fdr_corrected", "significance")])

# Combine specific plots using patchwork
# Combine specific plots using patchwork
# Combine specific plots using patchwork
final_plot1 <- plots[[1]] | plots[[2]] | plots[[3]] 
final_plot2 <- plots[[4]] | plots[[5]] | plots[[6]]
final_plot3 <- plots[[7]] | plots[[8]] | plots[[9]] 

# Display plots
print(final_plot1)
print(final_plot2)
print(final_plot3)
# Print original p-values and FDR-corrected p-values
print(all_p_values[, c("response_var", "comparison", "p_value", "fdr_corrected", "significance")])


# Define the path where you want to save the plot
save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/manuscript/Manuscript figures/final manuscript/PDF figures"

# Use ggsave to save your plot as PDF with high resolution
ggsave(paste0(save_path, "Cross_adcd.pdf"), 
       plot = final_plot3, 
       device = "pdf", 
       dpi = 1000,   # DPI doesn't have a significant effect for vector formats like PDF, but it's good to include
       width = 8, 
       height = 4, 
       units = "in")
dev.off()  # Close the device


##
# If you want to display it in a nicer table format (requires knitr package)
library(knitr)

# Display the results as a kable table for better visualization
kable(all_p_values, caption = "Statistical Comparisons and p-values")


# Function to calculate percentage of positive samples above LLOD for each group and variable
calculate_percentage_positive <- function(data, variables, LLOD_values, group_var = "EBV_infection") {
  results <- data.frame()
  
  for (variable in variables) {
    LLOD_value <- as.numeric(LLOD_values[[variable]])  # Get the LLOD for the variable
    
    # Calculate for each group in the data
    group_stats <- data %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        total_samples = n(),
        positive_count = sum(!!sym(variable) > LLOD_value, na.rm = TRUE),
        percentage_positive = (positive_count / total_samples) * 100
      ) %>%
      mutate(variable = variable)  # Add the variable name to the results
    
    results <- rbind(results, group_stats)
  }
  
  return(results)
}

# Calculate the percentage of positivity for each group and variable
percentage_positive_results <- calculate_percentage_positive(COMPILED_DATA_EBNA1b, variables_to_test, LLOD_values)

# Display the results
print(percentage_positive_results)

library(knitr)

# Display the results in a kable table for better visualization
kable(percentage_positive_results, caption = "Percentage of Samples Above LLOD by Group and Variable")

