# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)  
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
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(dplyr)
library(readxl)
library(patchwork)

# Load the dataset
COMPILED_DATA_EBNA1b <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
                                   sheet = "Data_2024")

COMPILED_DATA_EBNA1b <- as.data.frame(COMPILED_DATA_EBNA1b)

# Define common theme for the plots
labssizes <- theme(axis.title = element_text(size = 12, face = 'bold'),
                   axis.text = element_text(size = 10),
                   axis.text.x = element_text(angle = 45, hjust = 1),
                   plot.title = element_text(size = 12, face = 'bold', hjust = 0.5, vjust = 0.5))

# Define groups and colors
groups <- scale_x_discrete(limits = c("Acute (n = 97)", "6 weeks (n = 67)", "6 months (n = 30)", "1 year (n = 67)", "Sero_Pos (n = 50)", "NEGATIVE"))
mycolor <- scale_color_manual(breaks = c("Acute (n = 97)", "6 weeks (n = 67)", "6 months (n = 30)", "1 year (n = 67)", "Sero_Pos (n = 50)", "NEGATIVE"), 
                              values = c("magenta4", "dodgerblue3", "chocolate3", "darkgreen", "blue", "yellow4", "black"))
myfill <- scale_fill_manual(breaks = c("Acute (n = 97)", "6 weeks (n = 67)", "6 months (n = 30)", "1 year (n = 67)", "Sero_Pos (n = 50)", "NEGATIVE"), 
                            values = c("magenta4", "dodgerblue3", "chocolate3", "darkgreen", "blue", "yellow4", "black"))

# Function to perform ANOVA, adjust p-values, and filter significant results
perform_anova <- function(df, response_var, group_var) {
  tukey_result <- aov(as.formula(paste(response_var, "~", group_var)), data = df) %>% tukey_hsd()
  tukey_result$p.adj <- p.adjust(tukey_result$p.adj, method = "fdr")
  
  tukey_result$p.adj.signif <- ifelse(tukey_result$p.adj < 0.0001, "****",
                                      ifelse(tukey_result$p.adj < 0.001, "***", ""))
  
  tukey_result <- tukey_result %>% 
    filter(p.adj < 0.05 & p.adj.signif != "") %>%
    filter(
      (group1 == "Acute (n = 97)" & group2 %in% c("1 year (n = 67)", "Sero_Pos (n = 50)", "NEGATIVE")) |
        (group2 == "Acute (n = 97)" & group1 %in% c("1 year (n = 67)", "Sero_Pos (n = 50)", "NEGATIVE")) |
        (group1 == "Sero_Pos (n = 50)" & group2 == "NEGATIVE") |
        (group2 == "Sero_Pos (n = 50)" & group1 == "NEGATIVE") |
        (group1 == "1 year (n = 67)" & group2 == "NEGATIVE") |
        (group2 == "1 year (n = 67)" & group1 == "NEGATIVE")
    )
  
  return(tukey_result)
}
custom_labels <- c(
  "ADCD_ANO2" = "ANO2",
  "ADCD_GlialCAM" = "GlialCAM",
  "ADCD_CRYAB" = "CRYAB"
)

# Function to create ggplot for each response variable
create_plot <- function(df, response_var, group_var, tukey_result) {
  ggplot(df, aes_string(x = group_var, y = response_var)) +
    geom_violin(aes(fill = !!sym(group_var)), color = NA, alpha = 0.6, scale = "width", adjust = 1.5) +
    geom_boxplot(aes(color = !!sym(group_var)), width = 0.2, outlier.shape = NA, show.legend = FALSE) +
    scale_y_continuous(trans = "log10") +
    coord_cartesian(ylim = c(3000, 1e8)) +
    theme_classic() +
    # Use custom labels for the plot title
    labs(title = custom_labels[[response_var]], x = 'EBV infection status', y = 'C3 Mean Fluorescence Units', fill = 'Treatment') +  # Change made here
    labssizes + groups + mycolor + myfill +
    stat_pvalue_manual(data = tukey_result, label = "p.adj.signif", y.position = 6, step.increase = 0.07, size = 6) +
    theme(legend.position = "none")
}


# Define common parameters
group_var <- "EBV_infection"
response_vars <- c("ADCD_ANO2", "ADCD_GlialCAM", "ADCD_CRYAB")

# Loop through each response variable, perform analysis, and create plots
plots <- lapply(response_vars, function(response_var) {
  tukey_result <- perform_anova(COMPILED_DATA_EBNA1b, response_var, group_var)
  create_plot(COMPILED_DATA_EBNA1b, response_var, group_var, tukey_result)
})

# Combine specific plots using patchwork
final_plot3 <- plots[[1]] | plots[[2]] | plots[[3]] 

# Display plots (adjust according to your needs)
print(final_plot3)
#save
save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/manuscript/Manuscript figures/"

ggsave(paste0(save_path, "ADCD_cross.tiff"), 
       plot = final_plot3,  # Use `p` since this is the plot object
       device = "tiff", 
       dpi = 800, 
       width = 10, 
       height = 4, 
       units = "in")
