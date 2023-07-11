# Effect of Country and soilweight (soil-weight-class)
library(ggplot2)
library(xlsx)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)



########################
#       Plotten        #
########################
#####################
#Species_per_Country#
#####################
pythium_species <- read.xlsx("C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Relative_abundance/Relative_abundance_pythium_top_10.xlsx", sheetIndex=1)
pyth_metadata <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/output/metadata_pythium.txt", header = TRUE, sep = "\t")

pythium_species <- merge(pyth_metadata, pythium_species, by.x = "Id", by.y="Id")
pythium_species <- pythium_species %>% filter(Year !="2022")

pythium_species <- tidyr::gather(pythium_species, 
                                key = Taxa, 
                                value = Abundance,
                                s__sylvaticum:s__attrantheridium) 

pythium_species$Taxa <- str_replace(pythium_species$Taxa, "s__", "P. ")
pythium_species$Abundance <- as.numeric(pythium_species$Abundance)
pythium_species <- pythium_species %>% mutate(Abundance = Abundance*100)


pythium_species_soil <- pythium_species %>% filter(Compartment =="Soil")
df <- pythium_species_soil %>% filter(Country != "Czech_republic")
df <- df %>% mutate(Abundance = Abundance*100)
pyth_relativ_country_species_soil <- ggplot2::ggplot(df, aes(x=Country, y= Abundance, fill = Taxa)) +
                                        ggplot2::geom_bar(stat = "identity", position = "fill") +
                                        ggplot2::labs(y = "Relative abundance (%)", x = "") +
                                        ggplot2::scale_fill_brewer(palette = "Paired", name = "Species") +
                                        ggplot2::theme_bw() +
                                        ggplot2::scale_y_continuous(labels = scales::percent_format()) +
                                        ggplot2::scale_x_discrete(limits = c("Netherland", "Belgium", "Hungary", "Germany", "Spain",  "Romania", "France","Austria", "Switzerland", "Italy")) +
                                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
                                                    axis.text.y = element_text(size = 12, color = "black"),
                                                    axis.title.y = element_text(size = 12, color = "black"),
                                                    legend.position = "bottom",
                                                    legend.text = element_text(size = 12, face = "italic", color = "black"),
                                                    legend.title = element_text(size = 12, color = "black"),
                                                    legend.key.size = unit(0.5, "lines"))

pyth_relativ_country_species_soil
png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Relative_abundance/relative_abundance_pythium_soil_country.png", width = 5000, height = 2500, res = 500)
pyth_relativ_country_species_soil
dev.off()



###################################
#Species_per_sample_per_soilweight#
###################################
pyth_species_asvs <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/R/Phyloseq/Phyloseq_Objekte/Pythium/ASV_table_pythium_species.txt", header = TRUE, sep = "\t")
pyth_species_asvs <- pyth_species_asvs %>% mutate_if(is.numeric, ~1 * (. != 0))
pyth_metadata <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/output/metadata_pythium.txt", header = TRUE, sep = "\t")

pythium_matrix <- merge(pyth_metadata, pyth_species_asvs, by.x = "Id", by.y="Id")
pythium_matrix <- pythium_matrix %>% filter(Year != "2022")
pythium_matrix <- pythium_matrix %>% filter(Compartment == "Soil")
pythium_matrix <- pythium_matrix %>% filter(Country != "Czech_republic")

df <- tidyr::gather(pythium_matrix, 
                                key = Taxa, 
                                value = Abundance,
                                s__sylvaticum:s__paroecandrum)


df <- df %>% filter(Abundance != "0")
df <- df %>% group_by(Year, Location, Country, Soilweight) %>% summarise(n_taxa = n_distinct(Taxa))
df <- df %>% group_by(Location, Year, Country, Soilweight) %>% summarise(mean_taxa_count = mean(n_taxa))


pyth_relativ_Soilweight_species_soil <- ggplot2::ggplot(data = df, aes(x = Soilweight, y = mean_taxa_count, fill = Soilweight)) +
                                        ggplot2::geom_boxplot() +
                                        ggplot2::scale_fill_brewer(palette = "Paired", name = "weight class soil") +
                                        ggplot2::annotate("text", x = 2.01, y = 19.2, size = 7, label = "*", fontface = "bold", color = "black") +
                                        ggplot2::theme_bw() +
                                        ggplot2::labs(y = "Number of Species", x = "") +
                                        ggplot2::scale_x_discrete(limits = c("light", "medium", "heavy")) +
                                        ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
                                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
                                                        axis.text.y = element_text(size = 12, color = "black"),
                                                        axis.title.y = element_text(size = 12, color = "black"),
                                                        legend.position = "bottom",
                                                        legend.text = element_text(size = 12, color = "black"),
                                                        legend.title = element_text(size = 12, color = "black"),
                                                        legend.key.size = unit(0.5, "lines"))


pyth_relativ_Soilweight_species_soil




#######################################
#abundance_per_location_per_soilweight#
#######################################
pyth_species_asvs <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/R/Phyloseq/Phyloseq_Objekte/Pythium/ASV_table_pythium_species.txt", header = TRUE, sep = "\t")
pyth_metadata <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/output/metadata_pythium.txt", header = TRUE, sep = "\t")

pythium_matrix <- merge(pyth_metadata, pyth_species_asvs, by.x = "Id", by.y="Id")
pythium_matrix <- pythium_matrix %>% filter(Year != "2022")
pythium_matrix <- pythium_matrix %>% filter(Compartment == "Soil")
pythium_matrix <- pythium_matrix %>% filter(Country != "Czech_republic")

df <- tidyr::gather(pythium_matrix, 
                                key = Taxa, 
                                value = Abundance,
                                s__sylvaticum:s__paroecandrum)


df <- df %>% filter(Abundance != "0")
df <- df %>% group_by(Year, Location, Country, Soilweight) %>% summarise(Abundance = sum(Abundance))



pyth_relativ_soiltype_species_soil <- ggplot2::ggplot(data = df, aes(x = Soilweight, y = Abundance, fill = Soilweight)) +
                                        ggplot2::geom_boxplot() +
                                        ggplot2::scale_fill_brewer(palette = "Paired", name = "weight class soil") +
                                        ggplot2::theme_bw() +
                                        ggplot2::labs(y = "Abundance", x = "") +
                                        ggplot2::scale_x_discrete(limits = c("light", "medium", "heavy")) +
                                        ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
                                        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
                                                    axis.text.y = element_text(size = 12, color = "black"),
                                                    axis.title.y = element_text(size = 12, color = "black"),
                                                    legend.position = "bottom",
                                                    legend.text = element_text(size = 12, color = "black"),
                                                    legend.title = element_text(size = 12, color = "black"),
                                                    legend.key.size = unit(0.5, "lines"))


pyth_relativ_soiltype_species_soil





pyth_compar_soilweight <- ggarrange(pyth_relativ_soiltype_species_soil, pyth_relativ_Soilweight_species_soil, 
                                    labels =c("B", "C"),
                                    common.legend = TRUE,
                                    legend="bottom",
                                    nrow = 2)
pyth_compar_soilweight


pyth_country_compar_soilweight <- ggarrange(pyth_relativ_country_species_soil, pyth_compar_soilweight,
                                    labels =c("A", "B"),
                                    #common.legend = TRUE,
                                    legend="bottom",
                                    widths = c(0.6, 0.4), 
                                    ncol = 2)
pyth_country_compar_soilweight



png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Boxplot/pyth_relativ_soiltype.png", width = 9000, height = 4000, res = 600)
pyth_country_compar_soilweight
dev.off()




