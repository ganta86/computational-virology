# Load necessary libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(ggpubr)
library(tidyr)
library(patchwork)
library(rstatix)
library(reshape2)

# Load data
COMPILED_DATA_EBNA1b <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
                                   sheet = "Data_2024")
COMPILED_DATA_EBNA1b <- as.data.frame(COMPILED_DATA_EBNA1b)

# Define variables for urea treated and untreated
urea_vars <- c("Total_IgG_ANO2_urea", "Total_IgG_GlialCAM_urea", "Total_IgG_CRYAB_urea")
untreated_vars <- c("Total_IgG_ANO2", "Total_IgG_GlialCAM","Total_IgG_CRYAB")

# Convert to long format using melt
data_long <- melt(COMPILED_DATA_EBNA1b, measure.vars = c(urea_vars, untreated_vars), variable.name = "variable", value.name = "value")

# Add treatment information
data_long <- data_long %>%
  mutate(treatment = ifelse(variable %in% urea_vars, "urea", "no_urea"),
         variable = str_replace(variable, "_urea", ""),
         value = as.numeric(value))

# Remove non-numeric values
data_long <- data_long %>% filter(!is.na(value))

# Define common themes and limits
labssizes <- theme(axis.title = element_text(size = 14, face = 'bold'),
                   axis.text = element_text(size = 10),
                   axis.text.x = element_text(angle = 45, hjust = 1),
                   plot.title = element_text(size = 14, face = 'bold', hjust = 0.5, vjust = 0.5))

groups <- scale_x_discrete(limits = c("Acute (n = 97)", "6 weeks (n = 67)", "6 months (n = 30)", "1 year (n = 67)", "Sero_Pos (n = 50)", "Sero_Negative"))

# Extended color palette to cover all categories
mycolor <- scale_color_manual(breaks = c("urea", "no_urea"), values = c("seagreen", "burlywood"))

myfill <- scale_fill_manual(breaks = c("urea", "no_urea"), values = c("seagreen", "burlywood"))

# Function to create violin plot
create_violin_plot <- function(data, variable, title) {
  filtered_data <- data %>% filter(variable == !!variable)
  ggplot(filtered_data, aes(x = EBV_infection, y = value, fill = treatment)) +
    geom_violin(size = 0.5, alpha = 0.5, position = position_dodge(width = 0.9), width = 0.8, scale = "width", adjust = 1.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(trans = "log10", limits = c(1e1, 1e7)) +
    theme_classic() +
    theme(aspect.ratio = 1.2) + # Adjust aspect ratio
    labs(title = title, x = 'EBV_infection', y = 'Mean Fluorescence Units', fill = 'Treatment') +
    labssizes + groups + mycolor + myfill
}

# Create violin plots for each variable
plots <- list()
for (var in unique(data_long$variable)) {
  title <- paste(var)
  plots[[var]] <- create_violin_plot(data_long, var, title)
}

# Split plots into 4 groups for 4 individual plots, each with 3 columns
plot_groups <- split(plots, ceiling(seq_along(plots)/3))

# Combine plots into 4 individual plots, each with 3 columns
combined_plots <- lapply(plot_groups, function(group) {
  wrap_plots(group, ncol = 3)
})

# Display the combined plots
for (plot in combined_plots) {
  print(plot)
}
save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/manuscript/Manuscript figures/final manuscript/PDF figures/"

ggsave(paste0(save_path, "IgG_cross_urea.pdf"), 
       plot = plot, 
       device = "pdf", 
       dpi = 1000, 
       width = 12, 
       height = 4, 
       units = "in")
dev.off()  # Close the device
