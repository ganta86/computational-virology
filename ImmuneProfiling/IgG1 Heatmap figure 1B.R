library(ropls)
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

#Data_mean <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/viroscan/lUMNEX/Luminex 2024/COMPILED_DATA_EBNA1b.xlsx", 
#sheet = "Data_2024")

Data_mean <- read_excel("C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/EBV_DATA.xlsx", 
                        sheet = "IgG1")
Data_mean <- Data_mean %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)),
         across(where(is.character), ~ replace_na(., "")))

# Alternatively, if you want to replace NAs with 0 in numeric columns and "NA" in character columns
Data_mean <- Data_mean %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)),
         across(where(is.character), ~ replace_na(., "NA")))
Data_mean[ , 7:ncol(Data_mean)] <- lapply(Data_mean[ , 7:ncol(Data_mean)], as.numeric)

# This will create a vector of row names by pasting components needed
row_names <- paste (Data_mean$participants, 
                    Data_mean$`SAMPLE ID`,
                    Data_mean$EBV_infection,
                    sep = "_")



Data_mean <- Data_mean[, -c(1:6)]

desired_col_order <- c("1537_1537 v1_Acute (n = 97)",
                       "1531_1531 v1_Acute (n = 97)",
                       "1538_1538 v1_Acute (n = 97)",
                       "1536_1536 v1_Acute (n = 97)",
                       "1539_1539 v1_Acute (n = 97)",
                       "1527_1527 v1_Acute (n = 97)",
                       "1568_1568 v1_Acute (n = 97)",
                       "1582_1582 v1_Acute (n = 97)",
                       "1541_1541 v1_Acute (n = 97)",
                       "1545_1545 v1_Acute (n = 97)",
                       "1558_1558 v1_Acute (n = 97)",
                       "1572_1572 v1_Acute (n = 97)",
                       "1575_1575 v2_Acute (n = 97)",
                       "1585_1585 v1_Acute (n = 97)",
                       "1550_1550 v1_Acute (n = 97)",
                       "1564_1564 v1_Acute (n = 97)",
                       "1578_1578 v1_Acute (n = 97)",
                       "1542_1542 v1_Acute (n = 97)",
                       "1569_1569 v1_Acute (n = 97)",
                       "1583_1583 v1_Acute (n = 97)",
                       "1547_1547 v1_Acute (n = 97)",
                       "1559_1559 v1_Acute (n = 97)",
                       "1586_1586 v1_Acute (n = 97)",
                       "1576_1576 v1_Acute (n = 97)",
                       "1587_1587 v1_Acute (n = 97)",
                       "1555_1555 v1_Acute (n = 97)",
                       "1567_1567 v1_Acute (n = 97)",
                       "1580_1580 v1_Acute (n = 97)",
                       "1543_1543 v1_Acute (n = 97)",
                       "1570_1570 v1_Acute (n = 97)",
                       "1573_1573 v1_Acute (n = 97)",
                       "1584_1584 v1_Acute (n = 97)",
                       "1548_1548 v1_Acute (n = 97)",
                       "1556_1556 v1_Acute (n = 97)",
                       "1563_1563 v1_Acute (n = 97)",
                       "1577_1577 v1_Acute (n = 97)",
                       "1588_1588 v1_Acute (n = 97)",
                       "1591_1591 v1_Acute (n = 97)",
                       "1598_1598 v1_Acute (n = 97)",
                       "1612_1612 v1_Acute (n = 97)",
                       "1624_1624 v2_Acute (n = 97)",
                       "1603_1603 v2_Acute (n = 97)",
                       "1618_1618 v1_Acute (n = 97)",
                       "1627_1627 v1_Acute (n = 97)",
                       "1589_1589 v1_Acute (n = 97)",
                       "1593_1593 v2_Acute (n = 97)",
                       "1607_1607 v1_Acute (n = 97)",
                       "1622_1622 v1_Acute (n = 97)",
                       "1633_1633 v1_Acute (n = 97)",
                       "1600_1600 v1_Acute (n = 97)",
                       "1614_1614 v1_Acute (n = 97)",
                       "1625_1625 v1_Acute (n = 97)",
                       "1604_1604 v1_Acute (n = 97)",
                       "1620_1620 v1_Acute (n = 97)",
                       "1630_1630 v1_Acute (n = 97)",
                       "1590_1590 v1_Acute (n = 97)",
                       "1594_1594 v1_Acute (n = 97)",
                       "1611_1611 v1_Acute (n = 97)",
                       "1623_1623 v1_Acute (n = 97)",
                       "1634_1634 v1_Acute (n = 97)",
                       "1602_1602 v1_Acute (n = 97)",
                       "1617_1617 v1_Acute (n = 97)",
                       "1626_1626 v1_Acute (n = 97)",
                       "1605_1605 v1_Acute (n = 97)",
                       "1621_1621 v1_Acute (n = 97)",
                       "1632_1632 v2_Acute (n = 97)",
                       "1635_1635 v1_Acute (n = 97)",
                       "1649_1649 v1_Acute (n = 97)",
                       "1659_1659 v1_Acute (n = 97)",
                       "1665_1665 v1_Acute (n = 97)",
                       "1679_1679 v1_Acute (n = 97)",
                       "1639_1639 v1_Acute (n = 97)",
                       "1653_1653 v1_Acute (n = 97)",
                       "1670_1670 v1_Acute (n = 97)",
                       "1685_1685 v1_Acute (n = 97)",
                       "1646_1646 v1_Acute (n = 97)",
                       "1657_1657 v1_Acute (n = 97)",
                       "1675_1675 v1_Acute (n = 97)",
                       "1636_1636 v1_Acute (n = 97)",
                       "1651_1651 v2_Acute (n = 97)",
                       "1663_1663 v1_Acute (n = 97)",
                       "1666_1666 v1_Acute (n = 97)",
                       "1681_1681 v1_Acute (n = 97)",
                       "1641_1641 v1_Acute (n = 97)",
                       "1654_1654 v1_Acute (n = 97)",
                       "1671_1671 v1_Acute (n = 97)",
                       "1648_1648 v1_Acute (n = 97)",
                       "1658_1658 v1_Acute (n = 97)",
                       "1664_1664 v1_Acute (n = 97)",
                       "1676_1676 v1_Acute (n = 97)",
                       "1637_1637 v1_Acute (n = 97)",
                       "1652_1652 v1_Acute (n = 97)",
                       "1667_1667 v1_Acute (n = 97)",
                       "1684_1684 v1_Acute (n = 97)",
                       "1645_1645 v1_Acute (n = 97)",
                       "1655_1655 v2_Acute (n = 97)",
                       "1674_1674 v1_Acute (n = 97)",
                       "1527_1527 v6_6 weeks (n = 67)",
                       "1537_1537 v5_6 weeks (n = 67)",
                       "1531_1531 v6_6 weeks (n = 67)",
                       "1548_1548 v6_6 weeks (n = 67)",
                       "1563_1563 v6_6 weeks (n = 67)",
                       "1577_1577 v6_6 weeks (n = 67)",
                       "1568_1568 v6_6 weeks (n = 67)",
                       "1582_1582 v6_6 weeks (n = 67)",
                       "1545_1545 v6_6 weeks (n = 67)",
                       "1558_1558 v6_6 weeks (n = 67)",
                       "1575_1575 v6_6 weeks (n = 67)",
                       "1585_1585 v6_6 weeks (n = 67)",
                       "1550_1550 v6_6 weeks (n = 67)",
                       "1564_1564 v5_6 weeks (n = 67)",
                       "1542_1542 v6_6 weeks (n = 67)",
                       "1569_1569 v5_6 weeks (n = 67)",
                       "1559_1559 v6_6 weeks (n = 67)",
                       "1576_1576 v5_6 weeks (n = 67)",
                       "1587_1587 v5_6 weeks (n = 67)",
                       "1567_1567 v6_6 weeks (n = 67)",
                       "1543_1543 v6_6 weeks (n = 67)",
                       "1570_1570 v6_6 weeks (n = 67)",
                       "1573_1573 v5_6 weeks (n = 67)",
                       "1588_1588 v6_6 weeks (n = 67)",
                       "1632_1632 v6_6 weeks (n = 67)",
                       "1624_1624 v6_6 weeks (n = 67)",
                       "1618_1618 v5_6 weeks (n = 67)",
                       "1607_1607 v6_6 weeks (n = 67)",
                       "1622_1622 v6_6 weeks (n = 67)",
                       "1633_1633 v6_6 weeks (n = 67)",
                       "1600_1600 v6_6 weeks (n = 67)",
                       "1614_1614 v6_6 weeks (n = 67)",
                       "1625_1625 v5_6 weeks (n = 67)",
                       "1620_1620 v6_6 weeks (n = 67)",
                       "1630_1630 v6_6 weeks (n = 67)",
                       "1590_1590 v6_6 weeks (n = 67)",
                       "1594_1594 v5_6 weeks (n = 67)",
                       "1611_1611 v5_6 weeks (n = 67)",
                       "1623_1623 v5_6 weeks (n = 67)",
                       "1634_1634 v6_6 weeks (n = 67)",
                       "1617_1617 v6_6 weeks (n = 67)",
                       "1626_1626 v6_6 weeks (n = 67)",
                       "1645_1645 v6_6 weeks (n = 67)",
                       "1655_1655 v6_6 weeks (n = 67)",
                       "1674_1674 v6_6 weeks (n = 67)",
                       "1635_1635 v6_6 weeks (n = 67)",
                       "1649_1649 v6_6 weeks (n = 67)",
                       "1665_1665 v6_6 weeks (n = 67)",
                       "1679_1679 v6_6 weeks (n = 67)",
                       "1639_1639 v6_6 weeks (n = 67)",
                       "1653_1653 v5_6 weeks (n = 67)",
                       "1670_1670 v6_6 weeks (n = 67)",
                       "1646_1646 v6_6 weeks (n = 67)",
                       "1657_1657 v6_6 weeks (n = 67)",
                       "1675_1675 v6_6 weeks (n = 67)",
                       "1636_1636 v6_6 weeks (n = 67)",
                       "1651_1651 v6_6 weeks (n = 67)",
                       "1666_1666 v5_6 weeks (n = 67)",
                       "1681_1681 v6_6 weeks (n = 67)",
                       "1641_1641 v6_6 weeks (n = 67)",
                       "1648_1648 v6_6 weeks (n = 67)",
                       "1658_1658 v6_6 weeks (n = 67)",
                       "1664_1664 v6_6 weeks (n = 67)",
                       "1676_1676 v5_6 weeks (n = 67)",
                       "1637_1637 v6_6 weeks (n = 67)",
                       "1652_1652 v5_6 weeks (n = 67)",
                       "1667_1667 v6_6 weeks (n = 67)",
                       "1531_1531 v7_6 months (n = 30)",
                       "1538_1538 v7_6 months (n = 30)",
                       "1563_1563 v7_6 months (n = 30)",
                       "1577_1577 v7_6 months (n = 30)",
                       "1568_1568 v7_6 months (n = 30)",
                       "1582_1582 v7_6 months (n = 30)",
                       "1564_1564 v7_6 months (n = 30)",
                       "1578_1578 v7_6 months (n = 30)",
                       "1602_1602 v7_6 months (n = 30)",
                       "1621_1621 v7_6 months (n = 30)",
                       "1612_1612 v7_6 months (n = 30)",
                       "1624_1624 v7_6 months (n = 30)",
                       "1627_1627 v7_6 months (n = 30)",
                       "1589_1589 v7_6 months (n = 30)",
                       "1667_1667 v7_6 months (n = 30)",
                       "1684_1684 v7_6 months (n = 30)",
                       "1655_1655 v7_6 months (n = 30)",
                       "1674_1674 v7_6 months (n = 30)",
                       "1649_1649 v7_6 months (n = 30)",
                       "1665_1665 v7_6 months (n = 30)",
                       "1679_1679 v7_6 months (n = 30)",
                       "1670_1670 v7_6 months (n = 30)",
                       "1685_1685 v7_6 months (n = 30)",
                       "1675_1675 v7_6 months (n = 30)",
                       "1666_1666 v7_6 months (n = 30)",
                       "1681_1681 v7_6 months (n = 30)",
                       "1671_1671 v7_6 months (n = 30)",
                       "1648_1648 v7_6 months (n = 30)",
                       "1664_1664 v7_6 months (n = 30)",
                       "1676_1676 v7_6 months (n = 30)",
                       "1527_1527 v8_1 year (n = 67)",
                       "1537_1537 v8_1 year (n = 67)",
                       "1536_1536 v8_1 year (n = 67)",
                       "1539_1539 v8_1 year (n = 67)",
                       "1543_1543 v8_1 year (n = 67)",
                       "1556_1556 v8_1 year (n = 67)",
                       "1570_1570 v8_1 year (n = 67)",
                       "1573_1573 v8_1 year (n = 67)",
                       "1584_1584 v8_1 year (n = 67)",
                       "1548_1548 v8_1 year (n = 67)",
                       "1541_1541 v8_1 year (n = 67)",
                       "1545_1545 v8_1 year (n = 67)",
                       "1558_1558 v8_1 year (n = 67)",
                       "1572_1572 v8_1 year (n = 67)",
                       "1575_1575 v8_1 year (n = 67)",
                       "1585_1585 v8_1 year (n = 67)",
                       "1550_1550 v8_1 year (n = 67)",
                       "1542_1542 v8_1 year (n = 67)",
                       "1569_1569 v8_1 year (n = 67)",
                       "1586_1586 v8_1 year (n = 67)",
                       "1583_1583 v8_1 year (n = 67)",
                       "1547_1547 v8_1 year (n = 67)",
                       "1555_1555 v8_1 year (n = 67)",
                       "1559_1559 v8_1 year (n = 67)",
                       "1576_1576 v8_1 year (n = 67)",
                       "1587_1587 v8_1 year (n = 67)",
                       "1567_1567 v8_1 year (n = 67)",
                       "1580_1580 v8_1 year (n = 67)",
                       "1617_1617 v8_1 year (n = 67)",
                       "1626_1626 v8_1 year (n = 67)",
                       "1588_1588 v8_1 year (n = 67)",
                       "1591_1591 v8_1 year (n = 67)",
                       "1605_1605 v8_1 year (n = 67)",
                       "1632_1632 v8_1 year (n = 67)",
                       "1598_1598 v8_1 year (n = 67)",
                       "1603_1603 v8_1 year (n = 67)",
                       "1618_1618 v8_1 year (n = 67)",
                       "1593_1593 v8_1 year (n = 67)",
                       "1607_1607 v8_1 year (n = 67)",
                       "1622_1622 v8_1 year (n = 67)",
                       "1633_1633 v8_1 year (n = 67)",
                       "1600_1600 v8_1 year (n = 67)",
                       "1614_1614 v8_1 year (n = 67)",
                       "1625_1625 v8_1 year (n = 67)",
                       "1604_1604 v8_1 year (n = 67)",
                       "1620_1620 v8_1 year (n = 67)",
                       "1630_1630 v8_1 year (n = 67)",
                       "1590_1590 v8_1 year (n = 67)",
                       "1594_1594 v8_1 year (n = 67)",
                       "1611_1611 v8_1 year (n = 67)",
                       "1623_1623 v8_1 year (n = 67)",
                       "1634_1634 v8_1 year (n = 67)",
                       "1637_1637 v8_1 year (n = 67)",
                       "1652_1652 v8_1 year (n = 67)",
                       "1645_1645 v8_1 year (n = 67)",
                       "1635_1635 v8_1 year (n = 67)",
                       "1659_1659 v8_1 year (n = 67)",
                       "1639_1639 v8_1 year (n = 67)",
                       "1653_1653 v8_1 year (n = 67)",
                       "1646_1646 v8_1 year (n = 67)",
                       "1657_1657 v8_1 year (n = 67)",
                       "1663_1663 v8_1 year (n = 67)",
                       "1636_1636 v8_1 year (n = 67)",
                       "1651_1651 v8_1 year (n = 67)",
                       "1641_1641 v8_1 year (n = 67)",
                       "1654_1654 v8_1 year (n = 67)",
                       "1658_1658 v8_1 year (n = 67)",
                       "ES375_ES375_Sero_Pos (n = 50)",
                       "ES450_ES450_Sero_Pos (n = 50)",
                       "ES587_ES587_Sero_Pos (n = 50)",
                       "ES363_ES363_Sero_Pos (n = 50)",
                       "ES411_ES411_Sero_Pos (n = 50)",
                       "ES432_ES432_Sero_Pos (n = 50)",
                       "ES379_ES379_Sero_Pos (n = 50)",
                       "ES451_ES451_Sero_Pos (n = 50)",
                       "ES344_ES344_Sero_Pos (n = 50)",
                       "ES365_ES365_Sero_Pos (n = 50)",
                       "ES414_ES414_Sero_Pos (n = 50)",
                       "ES437_ES437_Sero_Pos (n = 50)",
                       "ES382_ES382_Sero_Pos (n = 50)",
                       "ES452_ES452_Sero_Pos (n = 50)",
                       "ES346_ES346_Sero_Pos (n = 50)",
                       "ES366_ES366_Sero_Pos (n = 50)",
                       "ES415_ES415_Sero_Pos (n = 50)",
                       "ES457_ES457_Sero_Pos (n = 50)",
                       "ES393_ES393_Sero_Pos (n = 50)",
                       "ES455_ES455_Sero_Pos (n = 50)",
                       "ES350_ES350_Sero_Pos (n = 50)",
                       "ES373_ES373_Sero_Pos (n = 50)",
                       "ES418_ES418_Sero_Pos (n = 50)",
                       "ES459_ES459_Sero_Pos (n = 50)",
                       "ES396_ES396_Sero_Pos (n = 50)",
                       "ES487_ES487_Sero_Pos (n = 50)",
                       "ES352_ES352_Sero_Pos (n = 50)",
                       "ES374_ES374_Sero_Pos (n = 50)",
                       "ES420_ES420_Sero_Pos (n = 50)",
                       "ES494_ES494_Sero_Pos (n = 50)",
                       "ES356_ES356_Sero_Pos (n = 50)",
                       "ES433_ES433_Sero_Pos (n = 50)",
                       "ES488_ES488_Sero_Pos (n = 50)",
                       "ES353_ES353_Sero_Pos (n = 50)",
                       "ES389_ES389_Sero_Pos (n = 50)",
                       "ES424_ES424_Sero_Pos (n = 50)",
                       "ES499_ES499_Sero_Pos (n = 50)",
                       "ES371_ES371_Sero_Pos (n = 50)",
                       "ES438_ES438_Sero_Pos (n = 50)",
                       "ES489_ES489_Sero_Pos (n = 50)",
                       "ES357_ES357_Sero_Pos (n = 50)",
                       "ES394_ES394_Sero_Pos (n = 50)",
                       "ES425_ES425_Sero_Pos (n = 50)",
                       "ES501_ES501_Sero_Pos (n = 50)",
                       "ES372_ES372_Sero_Pos (n = 50)",
                       "ES440_ES440_Sero_Pos (n = 50)",
                       "ES556_ES556_Sero_Pos (n = 50)",
                       "ES360_ES360_Sero_Pos (n = 50)",
                       "ES406_ES406_Sero_Pos (n = 50)",
                       "ES426_ES426_Sero_Pos (n = 50)",
                       "N_ES644_NEGATIVE",
                       "N_ES640_NEGATIVE",
                       "N_ES628_NEGATIVE",
                       "N_ES624_NEGATIVE",
                       "N_ES622_NEGATIVE",
                       "N_ES617_NEGATIVE",
                       "N_ES612_NEGATIVE",
                       "N_ES607_NEGATIVE",
                       "N_ES603_NEGATIVE",
                       "N_ES593_NEGATIVE",
                       "N_ES644_NEGATIVE",
                       "N_ES640_NEGATIVE",
                       "N_ES628_NEGATIVE",
                       "N_ES624_NEGATIVE",
                       "N_ES622_NEGATIVE"
)


#Center and scale the numeric columns
Data_mean_scale <- scale(Data_mean, center = TRUE, scale = TRUE)

#choose this line for plot heatmap for raw MFI
#Data_mean_scale = Data_mean

Data_mean_scale = as.matrix (Data_mean_scale)

rownames(Data_mean_scale) = row_names

Data_mean_scale = t(Data_mean_scale)

Data_mean_scale = Data_mean_scale[, desired_col_order]

colnames(Data_mean_scale) <- make.unique(colnames(Data_mean_scale))
EBVstatus <- gsub(".*\\d+_(.*)", "\\1", desired_col_order)

# Create the dataframe
df <- data.frame(
  EBVstatus = EBVstatus
)

# Print the dataframe
print(df)
row.names(df) <- colnames(Data_mean_scale)
# Generate colors for two lines of annotations

group_colors <- list(
  EBVstatus = c( "Acute (n = 97)" = "magenta4",
                 "6 weeks (n = 67)" = "dodgerblue3", 
                 "6 months (n = 30)" = "chocolate3", 
                 "1 year (n = 67)" = "darkgreen", 
                 "Sero_Pos (n = 50)" = "blue",
                 "NEGATIVE" = "black")
)

# Function to create heatmap for a specific range of rows
create_heatmap <- function(data, row_indices, title) {
  pheatmap(data[row_indices, ],
           border_color = "black",    
           cluster_rows = FALSE,  
           cluster_cols = FALSE,  
           color = cividis(100),  # Use cividis color palette,  
           annotation_col = df,  
           annotation_colors = group_colors,  
           annotation_legend = TRUE,  
           breaks = seq(0, 1, length.out = 100),  
           show_colnames = FALSE,  
           fontsize_row = 12,  
           fontsize = 12,  
           cellheight = 12,  
           main = title
  )
}
# Flatten the data to calculate the global mean and standard deviation
global_mean <- mean(as.matrix(Data_mean), na.rm = TRUE)
global_sd <- sd(as.matrix(Data_mean), na.rm = TRUE)

# Apply Z-score normalization across the entire dataset
Data_mean_global_zscore <- as.data.frame(lapply(Data_mean, function(x) (x - global_mean) / global_sd))

# Convert to matrix for heatmap generation
Data_mean_global_zscore <- as.matrix(Data_mean_global_zscore)
rownames(Data_mean_global_zscore) <- row_names
Data_mean_global_zscore <- t(Data_mean_global_zscore)
Data_mean_global_zscore <- Data_mean_global_zscore[, desired_col_order]

# Generate the heatmap for global Z-score normalization
# Step 1: Display the heatmap in the R console or RStudio viewer
create_heatmap(Data_mean_global_zscore, 1:11, "IgG1")

# Generate the heatmap as a ggplot object
heatmap_plot <- create_heatmap(Data_mean_global_zscore, 1:11, "IgG1")

# Define the output path with a filename and .pdf extension
save_path <- "C:/Users/gantak/OneDrive - UMass Chan Medical School/Desktop/EBV manuscript/data/manuscript/Manuscript figures/final manuscript/PDF figures/"

# Use ggsave to save the plot as PDF
ggsave(filename = paste0(save_path, "heatmap_EBNA_IgG1.pdf"), 
       plot = heatmap_plot, 
       device = "pdf", 
       dpi = 1000,  # DPI isn't critical for PDFs but can help with raster conversions
       width = 8, 
       height = 4, 
       units = "in")
dev.off()  # Close the device