#load_packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(ggpubr)
library(jpeg)
library(png)
library(grid)
library(cowplot)
library(patchwork)
library(figpatch)

######################
# Shoot-fresh-weight #
######################
Pathogenitätstest <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/Pathogenitätstest/Pathogenitätstest_Gewächshaus_Mais/Abbildungen/Mario/Pathogenitätstest_Gewächshaus_Dezember_2021.txt", header = TRUE, sep = "\t", na.strings = c("","NA"))


data_path <- Pathogenitätstest
data_path <- data_path %>% filter(Sorte == "Benedictio")
data_path <- data_path %>% filter(Pythiumstamm %in% c("NC", "P. attrantheridium", "P. ultimum"))
data <- ddply(data_path, c("Kultur", "Sorte", "Pythiumstamm", "Species"), summarize, mean_Sprossgewicht_pro_Pflanze_frisch = mean(Sprossgewicht_pro_Pflanze_frisch),SD = sd(Sprossgewicht_pro_Pflanze_frisch))


barplot_spross_frisch_pythiumstamm <- ggplot2::ggplot(data, aes(x = Pythiumstamm, y = mean_Sprossgewicht_pro_Pflanze_frisch, fill = Pythiumstamm)) +
                    ggplot2::geom_bar(stat = "summary", color ="black", position = "dodge") +
                    ggplot2::scale_fill_manual(values = c("NC" = "#ffff99", "P. attrantheridium" = "#1f78b4", "P. ultimum" = "#6a3d9a")) +
                    ggplot2::geom_errorbar(aes(ymin = mean_Sprossgewicht_pro_Pflanze_frisch-SD, ymax = mean_Sprossgewicht_pro_Pflanze_frisch+SD), position = position_dodge(width = 0.9), color = "black", size = 0.5, width = 0.3) +
                    ggplot2::annotate("text", x = 2, y = data$mean_Sprossgewicht_pro_Pflanze_frisch[2] + data$SD[2] + 0.1, size = 4, label = "", fontface = "bold", color = "black") +
                    ggplot2::annotate("text", x = 3, y = data$mean_Sprossgewicht_pro_Pflanze_frisch[3] + data$SD[3] + 0.1, size = 4, label = "**", fontface = "bold", color = "black") +   
                            ggplot2::scale_y_continuous(name = "shoot fresh weight/plant (g)", limits = c(0,6), breaks =c(0, 1, 2, 3, 4, 5, 6)) +
                            ggplot2::scale_x_discrete(name = "") +
                            ggplot2::theme(legend.position = "none") +
                            ggplot2::labs(fill = "Species") +
                            ggplot2::theme_bw() +
                            ggplot2::theme(text = element_text(size = 12, color = "black")) +
                            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")) +
                                        theme(axis.text.y = element_text(size = 12, color = "black"), 
                                        legend.text = element_text(size = 12, face = "italic"),
                                        axis.text.x = element_text(size = 12, face = "italic", color = "black"),
                                        axis.title.y = element_text(size = 12, color = "black"))


barplot_spross_frisch_pythiumstamm


###################
# Root-dry-weight #
###################
Pathogenitätstest <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/Pathogenitätstest/Pathogenitätstest_Gewächshaus_Mais/Abbildungen/Mario/Pathogenitätstest_Gewächshaus_Dezember_2021_Abbildungen_Wurzel.txt", header = TRUE, sep = "\t", na.strings = c("","NA"))


data_path_beiz <- Pathogenitätstest
data_path_beiz <- data_path_beiz %>% filter(Sorte == "Benedictio")
data_path_beiz <- data_path_beiz %>% filter(Pythiumstamm %in% c("NC", "P. attrantheridium", "P. ultimum"))
data <- ddply(data_path_beiz, c("Kultur", "Sorte", "Pythiumstamm", "Species"), summarize, Wurzelgewicht_pro_Pflanze = mean(Wurzelgewicht_pro_Pflanze_trocken),SD = sd(Wurzelgewicht_pro_Pflanze_trocken))


barplot_wurzel_trocken_sorte <- ggplot2::ggplot(data, aes(x = Pythiumstamm, y = Wurzelgewicht_pro_Pflanze, fill = Pythiumstamm)) +
                    ggplot2::geom_bar(stat = "summary", color ="black", position = "dodge") +
                    ggplot2::scale_fill_manual(values = c("NC" = "#ffff99", "P. attrantheridium" = "#1f78b4", "P. ultimum" = "#6a3d9a")) +
                    ggplot2::geom_errorbar(aes(ymin = Wurzelgewicht_pro_Pflanze-SD, ymax = Wurzelgewicht_pro_Pflanze+SD), position = position_dodge(width = 0.9), color = "black", size = 0.5, width = 0.3) +
                    ggplot2::annotate("text", x = 2, y = data$Wurzelgewicht_pro_Pflanze[2] + data$SD[2] + 0.005, size = 4, label = "", fontface = "bold", color = "black") +
                    ggplot2::annotate("text", x = 3, y = data$Wurzelgewicht_pro_Pflanze[3] + data$SD[3] + 0.005, size = 4, label = "**", fontface = "bold", color = "black") +   
                            ggplot2::scale_y_continuous(name = "root dry weight/plant (g)", limits = c(0,0.3), breaks =c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)) +
                            ggplot2::scale_x_discrete(name = "") +
                            ggplot2::theme(legend.position = "none") +
                            ggplot2::labs(fill = "Species") +
                            ggplot2::theme_bw() +
                            ggplot2::theme(text = element_text(size = 12, color = "black")) +
                            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")) +
                                        theme(axis.text.y = element_text(size = 12, color = "black"), 
                                        legend.text = element_text(size = 12, face = "italic"),
                                        axis.text.x = element_text(size = 12, face = "italic", color = "black"),
                                        axis.title.y = element_text(size = 12, color = "black"))


barplot_wurzel_trocken_sorte


###############
# Germination #
###############
Pathogenitätstest <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/Pathogenitätstest/Pathogenitätstest_Gewächshaus_Mais/Abbildungen/Mario/Pathogenitätstest_Gewächshaus_Dezember_2021.txt", header = TRUE, sep = "\t", na.strings = c("","NA"))
data_path <- Pathogenitätstest
data_path <- data_path %>% filter(Sorte == "Benedictio")
data_path <- data_path %>% filter(Pythiumstamm %in% c("NC", "P. attrantheridium", "P. ultimum"))


data <- ddply(data_path, c("Kultur", "Sorte", "Pythiumstamm"), summarize, Germination = mean(Keimung*100),SD = sd(Keimung*100))

barplot_keimung_sorte <- ggplot2::ggplot(data, aes(x = Pythiumstamm, y = Germination, fill = Pythiumstamm)) +
                                    ggplot2::geom_bar(stat = "summary", color ="black", position = "dodge") +
                                    ggplot2::scale_fill_manual(values = c("NC" = "#ffff99", "P. attrantheridium" = "#1f78b4", "P. ultimum" = "#6a3d9a")) +
                                    ggplot2::geom_errorbar(aes(ymin = Germination-SD, ymax = Germination), position = position_dodge(width = 0.9), color = "black", size = 0.5, width = 0.3) +
                                    ggplot2::annotate("text", x = 2, y = data$Germination[2] + 5, size = 4, label = "", fontface = "bold", color = "black") +
                                    ggplot2::annotate("text", x = 3, y = data$Germination[3] + 5, size = 4, label = "**", fontface = "bold", color = "black") +   
                                    ggplot2::scale_y_continuous(name = "germination capacity (%)", limits = c(0,100)) +
                                    ggplot2::scale_x_discrete(name = "") +
                                    ggplot2::theme(legend.position = "none") +
                                    ggplot2::labs(fill = "Species") +
                                    ggplot2::theme_bw() +
                                    ggplot2::theme(text = element_text(size = 12, color = "black")) +
                                    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")) +
                                                theme(axis.text.y = element_text(size = 12, color = "black"), 
                                                legend.text = element_text(size = 12, face = "italic"),
                                                axis.text.x = element_text(size = 12, face = "italic", color = "black"),
                                                axis.title.y = element_text(size = 12, color = "black"))

barplot_keimung_sorte



#####################
# Plot_all_together #
#####################
pyth_pathogenität <- ggarrange(barplot_keimung_sorte, barplot_spross_frisch_pythiumstamm, barplot_wurzel_trocken_sorte,
                        common.legend = TRUE,
                        legend = "bottom",
                      nrow = 1,
                      ncol = 3)
pyth_pathogenität



png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Pathogenitätstest/pyth_pathogenität.png", width = 933, height = 850, res = 175)
pyth_pathogenität
dev.off()

