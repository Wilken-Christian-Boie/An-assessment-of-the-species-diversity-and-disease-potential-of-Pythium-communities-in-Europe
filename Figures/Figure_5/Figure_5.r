library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)

set.seed(42)

# countries comparison
# import data
pythium_species <- read.table("../Data/Relative_abundance_top_10_Pythium_species.txt", header = TRUE, sep = "\t")
pyth_metadata <- read.table("../../../output/Metadata_pythium.txt", header = TRUE, sep = "\t")

# prepare data
pythium_species <- merge(pyth_metadata, pythium_species, by = "Id")
pythium_species <- tidyr::gather(data = pythium_species,
                                 key = Taxa,
                                 value = Abundance,
                                 Globisporangium.sylvaticum:Globisporangium.intermedium)

pythium_species$Abundance <- as.numeric(pythium_species$Abundance) * 100
pythium_species$Taxa <- gsub("Pythium.", "P. ", pythium_species$Taxa)
pythium_species$Taxa <- gsub("Globisporangium.", "G. ", pythium_species$Taxa)
pythium_species$Taxa <- gsub("_", " ", pythium_species$Taxa)
df <- pythium_species %>% filter(Country != "Czech_republic")
countries <- c("Netherland", "Belgium", "Hungary", "Germany", "Spain", "Romania", "France", "Austria", "Switzerland", "Italy")

pythium_relative_country <- ggplot(df, aes(x = Country, y = Abundance, fill = Taxa)) +
                                geom_bar(stat = "identity", position = "fill") +
                                labs(y = "Relative abundance (%)", x = "") +
                                scale_fill_brewer(palette = "Paired", name = "Species") +
                                theme_bw() +
                                scale_y_continuous(expand = c(0, 0), labels = scales::percent_format()) +
                                scale_x_discrete(limits = countries) +
                                guides(fill = guide_legend(nrow = 3)) +
                                theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
                                      axis.text.y = element_text(size = 12, color = "black"),
                                      axis.title.y = element_text(size = 12, color = "black"),
                                      legend.position = "bottom",
                                      legend.text = element_text(size = 12, face = "italic", color = "black"),
                                      legend.title = element_text(size = 12, color = "black"),
                                      legend.key.size = unit(0.5, "lines"))

# number of species per soilweight
# import data
pyth_species_asvs <- read.table("../Data/ASV_table_pythium_species.txt", header = TRUE, sep = "\t")
pyth_metadata <- read.table("../../../output/Metadata_pythium.txt", header = TRUE, sep = "\t")

# prepare data
pyth_species_asvs <- pyth_species_asvs %>% mutate_if(is.numeric, ~1 * (. != 0))
pythium_matrix <- merge(pyth_metadata, pyth_species_asvs, by.x = "Id", by.y = "Id")
pythium_matrix <- pythium_matrix %>% filter(Country != "Czech_republic")

df <- tidyr::gather(data = pythium_matrix,
                    key = Taxa,
                    value = Abundance,
                    Globisporangium.sylvaticum:Globisporangium.terrestris)

df <- df %>%
        dplyr::filter(Abundance != "0") %>%
        dplyr::group_by(Year, Location, Country, Soilweight) %>%
        dplyr::summarize(n_taxa = n_distinct(Taxa),
                         mean_taxa_count = mean(n_taxa), .groups = "drop")

pythium_relative_soilweight <- ggplot(data = df, aes(x = Soilweight, y = mean_taxa_count, fill = Soilweight)) +
                                 geom_boxplot() +
						 	     geom_jitter(width = 0.2, color = "gray30") + 
                                 scale_fill_brewer(palette = "Paired", name = "weight class soil") +
                                 annotate("text", x = 2.01, y = 24, size = 3, label = "p = 0.089", color = "black") +
                                 labs(y = "Number of species", x = "") +
                                 scale_y_continuous(limits = c(0, 27)) +
                                 scale_x_discrete(limits = c("light", "medium", "heavy")) +
                                 theme_bw() +
                                 theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic"),
                                       axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
                                       axis.text.y = element_text(size = 12, color = "black"),
                                       axis.title.y = element_text(size = 12, color = "black"),
                                       legend.position = "none",
                                       legend.text = element_text(size = 12, color = "black"),
                                       legend.title = element_text(size = 12, color = "black"),
                                       legend.key.size = unit(0.5, "lines"))



pyth_species_asvs <- read.table("../Data/ASV_table_pythium_species.txt", header = TRUE, sep = "\t")
pyth_metadata <- read.table("../../../output/Metadata_pythium.txt", header = TRUE, sep = "\t")

# prepare data
pythium_matrix <- merge(pyth_metadata, pyth_species_asvs, by = "Id")

df <- tidyr::gather(data = pythium_matrix,
                    key = Taxa,
                    value = Abundance,
                    Globisporangium.sylvaticum:Globisporangium.terrestris) %>%
        dplyr::filter(Abundance != "0") %>%
        dplyr::group_by(Year, Location, Country, Soilweight) %>%
        dplyr::summarize(Abundance = sum(Abundance), .groups = "drop")

# plot data
pythium_relative_soiltype <- ggplot(data = df, aes(x = Soilweight, y = Abundance, fill = Soilweight)) +
                                 geom_boxplot() +
						 	     geom_jitter(width = 0.2, color = "gray30") + 
                                 scale_fill_brewer(palette = "Paired", name = "weight class soil") +
                                 labs(y = "Total reads", x = "") +
                                 scale_x_discrete(limits = c("light", "medium", "heavy")) +
								 scale_y_continuous(limits = c(0, 55000)) +
                                 theme_bw() +
                                 theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic"),
                                       axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, color = "black"),
                                       axis.text.y = element_text(size = 12, color = "black"),
                                       axis.title.y = element_text(size = 12, color = "black"),
                                       legend.position = "bottom",
                                       legend.text = element_text(size = 12, color = "black"),
                                       legend.title = element_text(size = 12, color = "black"),
                                       legend.key.size = unit(0.5, "lines"))

pythium_compare_soilweight <- ggarrange(pythium_relative_soiltype, pythium_relative_soilweight,
                                labels = c("B", "C"),
                                common.legend = TRUE,
                                align = "v",
                                legend = "bottom",
                                nrow = 2)

pythium_country_compare_soilweight <- ggarrange(pythium_relative_country, pythium_compare_soilweight,
                                        labels = c("A", "B"),
                                        legend = "bottom",
                                        widths = c(0.7, 0.3),
                                        ncol = 2)

jpeg(file = "Figure 5.jpg", width = 7000, height = 4000, res = 600)
pythium_country_compare_soilweight
dev.off()

pdf(file = "Figure 5.pdf", width = 14, height = 8)
pythium_country_compare_soilweight
dev.off()

