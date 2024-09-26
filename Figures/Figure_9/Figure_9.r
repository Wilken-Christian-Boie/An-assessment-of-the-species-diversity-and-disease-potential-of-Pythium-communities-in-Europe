# loading packages
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)
library(patchwork)

set.seed(42)

color_codes <- c("NC" = "#FFFF99",
                 "G. attrantheridium" = "#1F78B4",
                 "G. ultimum var. ultimum" = "#FDBF6F")

# import data
Pathogenicity_test <- read.table("../Data/Data_Pathogenicity.txt", header = TRUE, sep = "\t", na.strings = c("", "NA"))

# germination
data_path_germination <- Pathogenicity_test
data <- ddply(data_path_germination,
              c("Cultivar", "Variety", "Pythium_species"),
              summarize,
              Germination = mean(Germination_rate * 100),
              SD = sd(Germination_rate * 100))

data$Pythium_species <- factor(data$Pythium_species, levels = c("NC", "G. attrantheridium", "G. ultimum var. ultimum"))
barplot_germination_rate <- ggplot(data, aes(x = Pythium_species, y = Germination, fill = Pythium_species)) +
                                geom_bar(stat = "summary", fun = "mean", color = "black", position = "dodge") +
								geom_point(data = data_path_germination, aes(x = Pythium_species, y = Germination_rate * 100, fill = Pythium_species)) +
                                scale_fill_manual(values = color_codes) +
                                geom_errorbar(aes(ymin = Germination - SD, ymax = Germination), position = position_dodge(width = 0.9), color = "black", linewidth = 0.5, width = 0.3) +
                                annotate("text", x = 2, y = data$Germination[2] + 5, size = 5, label = "", fontface = "bold", color = "black") +
                                annotate("text", x = 3, y = data$Germination[3] + 5, size = 5, label = "*", fontface = "bold", color = "black") +
                                scale_y_continuous(name = "Germination rate (%)", limits = c(0, 100), expand = c(0, 0)) + 
                                scale_x_discrete(name = "") +
                                labs(fill = "Species") +
                                theme_bw() +
                                theme(text = element_text(size = 16, color = "black"),
                                      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 12, face = "italic"),
                                      axis.text.y = element_text(size = 16, color = "black"),
                                      legend.text = element_text(size = 16, face = "italic"),
                                      axis.title.y = element_text(size = 16, color = "black"),
                                      legend.position = "none")

# shoot fresh-weight
data_path_shoot <- Pathogenicity_test
data <- ddply(data_path_shoot,
              c("Cultivar", "Variety", "Pythium_species", "Species"),
              summarize,
              mean_Shoot_weight_fresh_per_plant = mean(Shoot_weight_fresh_per_plant),
              SD = sd(Shoot_weight_fresh_per_plant))

data$Pythium_species <- factor(data$Pythium_species, levels = c("NC", "G. attrantheridium", "G. ultimum var. ultimum"))
barplot_shoot_freshweight <- ggplot(data, aes(x = Pythium_species, y = mean_Shoot_weight_fresh_per_plant, fill = Pythium_species)) +
                                geom_bar(stat = "summary", fun = "mean", color = "black", position = "dodge") +
								geom_point(data = data_path_shoot, aes(x = Pythium_species, y = Shoot_weight_fresh_per_plant, fill = Pythium_species)) +
                                scale_fill_manual(values = color_codes) +
                                geom_errorbar(aes(ymin = mean_Shoot_weight_fresh_per_plant - SD, ymax = mean_Shoot_weight_fresh_per_plant + SD),
                                              position = position_dodge(width = 0.9),
                                              color = "black",
                                              linewidth = 0.5,
                                              width = 0.3) +
                                annotate("text", x = 2, y = data$mean_Shoot_weight_fresh_per_plant[2] + data$SD[2] + 0.1, size = 4, label = "", fontface = "bold", color = "black") +
                                annotate("text", x = 3, y = 3.5, size = 5, label = "*", fontface = "bold", color = "black") +
                                scale_y_continuous(name = "Shoot fresh weight/plant (g)", limits = c(0, 6), breaks = seq(0, 6, 1), expand = c(0, 0)) +
                                scale_x_discrete(name = "") +
                                labs(fill = "Species") +
                                theme_bw() +
                                theme(text = element_text(size = 16, color = "black"),
                                      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 12, face = "italic"),
                                      axis.text.y = element_text(size = 16, color = "black"),
                                      legend.text = element_text(size = 16, face = "italic"),
                                      legend.position = "none",
                                      axis.title.y = element_text(size = 16, color = "black"))

# root dry-weight
data_path_root <- Pathogenicity_test
data <- ddply(data_path_root,
              c("Cultivar", "Variety", "Pythium_species", "Species"),
              summarize,
              mean_Root_weight_dry_per_plant = mean(Root_weight_dry_per_plant),
              SD = sd(Root_weight_dry_per_plant))

data$Pythium_species <- factor(data$Pythium_species, levels = c("NC", "G. attrantheridium", "G. ultimum var. ultimum"))
barplot_root_dryweight <- ggplot(data, aes(x = Pythium_species, y = mean_Root_weight_dry_per_plant, fill = Pythium_species)) +
                            geom_bar(stat = "summary", fun = "mean", color = "black", position = "dodge") +
							geom_point(data = data_path_root, aes(x = Pythium_species, y = Root_weight_dry_per_plant, fill = Pythium_species)) +
                            scale_fill_manual(values = color_codes) +
                            geom_errorbar(aes(ymin = mean_Root_weight_dry_per_plant - SD, ymax = mean_Root_weight_dry_per_plant + SD),
                                          position = position_dodge(width = 0.9),
                                          color = "black",
                                          linewidth = 0.5,
                                          width = 0.3) +
                            annotate("text", x = 2, y = data$mean_Root_weight_dry_per_plant[2] + data$SD[2] + 0.005, size = 4, label = "", fontface = "bold", color = "black") +
                            annotate("text", x = 3, y = 0.2, size = 5, label = "*", fontface = "bold", color = "black") +
                            scale_y_continuous(name = "Root dry weight/plant (g)", limits = c(0, 0.3), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), expand = c(0, 0)) +
                            scale_x_discrete(name = "") +
                            theme(legend.position = "none") +
                            labs(fill = "Species") +
                            theme_bw() +
                            theme(text = element_text(size = 16, color = "black"),
                                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 12, face = "italic"),
                                  axis.text.y = element_text(size = 16, color = "black"),
                                  legend.text = element_text(size = 16, face = "italic"),
                                  legend.position = "none",
                                  axis.title.y = element_text(size = 16, color = "black"))

pyth_pathogenicity <- ggarrange(barplot_germination_rate, barplot_shoot_freshweight, barplot_root_dryweight,
                                nrow = 1,
                                ncol = 3)

jpeg(file = "Figure_10.jpg", width = 4500, height = 4500, res = 600)
pyth_pathogenicity
dev.off()

pdf(file = "Figure_10.pdf", width = 9, height = 9)
pyth_pathogenicity
dev.off()
