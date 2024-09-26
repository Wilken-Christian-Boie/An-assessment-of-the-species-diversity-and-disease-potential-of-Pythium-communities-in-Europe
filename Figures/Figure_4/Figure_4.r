library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)

set.seed(42)

# import data
pythium_species <- read.table("../Data/Relative_abundance_top_10_Pythium_species.txt", header = TRUE, sep = "\t")
pyth_metadata <- read.table("../../../output/Metadata_pythium.txt", header = TRUE, sep = "\t")

# prepare data
pythium_species <- merge(pyth_metadata, pythium_species, by = "Id")
pythium_species <- tidyr::gather(data = pythium_species,
                                 key = Taxa,
                                 value = Abundance,
                                 Globisporangium.sylvaticum:Globisporangium.intermedium) %>%
                   na.omit()

pythium_species$Abundance <- as.numeric(pythium_species$Abundance) * 100
pythium_species$Year <- as.factor(pythium_species$Year)
pythium_species$Taxa <- gsub("Pythium.", "P. ", pythium_species$Taxa)
pythium_species$Taxa <- gsub("Globisporangium.", "G. ", pythium_species$Taxa)
pythium_species$Taxa <- gsub("_", " ", pythium_species$Taxa)

# plot stacked data per species
pythium_relative_year_species <- ggplot(pythium_species, aes(x = Year, y = Abundance, fill = Taxa)) +
                                    geom_bar(stat = "identity", position = "fill") +
                                    labs(y = "Relative abundance (%)", x = "") +
                                    scale_fill_brewer(palette = "Paired", name = "Species") +
                                    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
                                    theme_bw() +
                                    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
                                        axis.text.y = element_text(size = 12, color = "black"),
                                        axis.title.y = element_text(size = 12, color = "black"),
                                        legend.position = "right",
                                        legend.text = element_text(size = 12, face = "italic", color = "black"),
                                        legend.title = element_text(size = 12, color = "black"),
                                        legend.key.size = unit(0.5, "lines"))
pythium_relative_year_species

# helper function to get number of observations per group
n_fun <- function(x) return(data.frame(y = median(x), label = paste0("n = ",length(x))))

# plot relative abundance per species
pythium_boxplot_species_year <- ggplot(pythium_species, aes(x = Year, y = Abundance, fill = Taxa)) +
                                    geom_boxplot() +
									geom_jitter(width = 0.2, color = "gray30") +
                                    stat_summary(fun = "mean", shape = 13, size = 0.4, color = "grey20") +
                                    facet_wrap(. ~ Taxa, ncol = 5) +
                                    scale_fill_brewer(palette = "Paired", name = "Species") +
                                    labs(y = "Relative abundance (%)", x = "") +
                                    theme_bw() +
                                    theme(strip.background = element_rect(fill = "white"),
                                          strip.text.x = element_text(size = 9, color = "black", face = "bold.italic"),
                                          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
                                          axis.text.y = element_text(size = 12, color = "black"),
                                          axis.title.y = element_text(size = 12, color = "black"),
                                          legend.position = "right",
                                          legend.text = element_text(size = 12, face = "italic", color = "black"),
                                          legend.title = element_text(size = 12, color = "black"),
                                          legend.key.size = unit(0.5, "lines"))

pythium_relative_plot <- ggarrange(pythium_relative_year_species, pythium_boxplot_species_year,
                                   labels = c("A", "B"),
                                   common.legend = TRUE,
                                   widths = c(25, 75),
                                   legend = "bottom",
                                   ncol = 2)

jpeg(file = "Figure 3.jpg", width = 7000, height = 4000, res = 600)
pythium_relative_plot
dev.off()

pdf(file = "Figure 3.pdf", width = 14, height = 8)
pythium_relative_plot
dev.off()

