library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

set.seed(42)

# import data
tax_classes <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
pythium_count_matrix <- read.table("../../../output/ASVs_Pythium_counts.txt", header = TRUE, sep = "\t")
pyth_metadata <- read.table("../../../output/Metadata_pythium.txt", header = TRUE, sep = "\t")
pyth_metadata <- pyth_metadata[pyth_metadata$Compartment == "Soil", ]
top10_names <- read.table("../Data/top10_names_soil.txt", sep = "\t", head = FALSE, col.names = "Species")
new_taxonomy <- read.table("../Data/new_taxonomy_and_clades.txt", sep = "\t", header = TRUE) %>%
				dplyr::select(Genus_neu, Species)

pythium_taxa <- read.table("../../../output/ASVs_Pythium_taxonomy.txt", header = FALSE, sep = "\t", col.names = c("ASV", "taxonomy")) %>%
                  tidyr::separate(col = taxonomy, into = tax_classes, sep = " ") %>%
				  dplyr::left_join(new_taxonomy, by = "Species") %>%
                  tidyr::unite(col = "Organism", Genus_neu, Species, sep = " ") %>%
                  dplyr::select(ASV, Organism) %>%
				  dplyr::filter(Organism != "Pythium sp._nov")

pythium_taxa <- pythium_taxa %>% dplyr::filter(Organism %in% top10_names$Species)
pythium_taxa$Organism <- gsub("Pythium", "P. ", pythium_taxa$Organism)
pythium_taxa$Organism <- gsub("Globisporangium", "G. ", pythium_taxa$Organism)
pythium_taxa$Organism <- gsub("_", " ", pythium_taxa$Organism)

pythium_count_matrix <- tidyr::gather(data = pythium_count_matrix,
                                      key = Id,
                                      value = Count,
                                      year2019_bait_P1_adaptertrimmed_Pythium_R1.fastq:year2021_soil_P9_adaptertrimmed_Pythium_R1.fastq)

colnames(pythium_count_matrix)[1] <- "ASV"
pythium_count_taxa <- merge(pythium_count_matrix, pythium_taxa, by = "ASV")
pythium_species <- merge(pyth_metadata, pythium_count_taxa, by = "Id")
pythium_species$Year <- as.factor(pythium_species$Year)
pythium_species <- subset(pythium_species, Count >= 1)
pythium_species <- pythium_species %>% mutate_if(is.numeric, ~1 * (. > 0))

pythium_species_für_df1 <- pythium_species %>% group_by(ASV, Organism, Year) %>% dplyr::summarise(count = n(), .groups = "drop")
df1 <- pythium_species_für_df1 %>% group_by(Organism, Year) %>% dplyr::summarise(count = n(), .groups = "drop")
pythium_species_für_df2 <- pythium_species %>% group_by(ASV, Organism, Compartment) %>% dplyr::summarise(count = n(), .groups = "drop")
df2 <- pythium_species_für_df2 %>% group_by(Organism, Compartment) %>% dplyr::summarise(count = n(), .groups = "drop")

df2[df2 == "Soil"] <- "All years"
colnames(df2)[2] <- "Year"
df3 <- bind_rows(df1, df2)
df3$Organism <- gsub("Pythium", "P. ", df3$Organism)
df3$Organism <- gsub("Globisporangium", "G. ", df3$Organism)

# plot bubble plot
bubble_plot <- ggplot(df3, aes(x = Organism, y = Year, size = count, fill = Organism)) +
                geom_point(shape = 21) +
                scale_size(range = c(1, 30), name = "Number of ASVs", breaks = c(5, 10, 30, 50, 70)) +
                geom_hline(yintercept = 3.5, linewidth = 1) +
                scale_x_discrete(name = "") +
                geom_text(aes(label = ifelse(Year == "All years", count, ""), size = 4, vjust = 3.4), show.legend = FALSE) +
                theme_bw() +
                theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15, face = "italic", color = "black"),
                      axis.text.y = element_text(size = 15, color = "black"),
                      axis.title.y = element_text(size = 15, color = "black"),
                      legend.position = "right",
                      legend.text = element_text(face = "italic", size = 15, color = "black"),
                      legend.title = element_text(size = 15, color = "black"),
                      legend.key.size = unit(1, "lines"))

jpeg(file = "Figure 6.jpg", width = 7000, height = 5000, res = 600)
bubble_plot + scale_fill_brewer(palette = "Paired", name = "Species", guide = "none")
dev.off()

pdf(file = "Figure 6.pdf", width = 14, height = 10)
bubble_plot + scale_fill_brewer(palette = "Paired", name = "Species", guide = "none")
dev.off()

write.table(df3, file = "../Data/Number_ASV_per_year_and_species.txt", sep = "\t", row.names = FALSE, quote = FALSE)
