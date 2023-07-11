# Jahres_Vergleich
library(ggplot2)
library(xlsx)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggstar)
library(stringr)


#########################
#import and prepare OTUs#
#########################
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
########################
#       Plotten        #
########################

#################
#Phylum_per_Year#
#################
pythium_species_soil <- pythium_species %>% filter(Compartment =="Soil")
pyth_relativ_year_species_soil <- ggplot2::ggplot(pythium_species_soil, aes(x=Year, y= Abundance, fill = Taxa)) +
                                    ggplot2::geom_bar(stat = "identity", position = "fill") +
                                    ggplot2::labs(y = "Relative abundance (%)", x = "") +
                                    ggplot2::scale_fill_brewer(palette = "Paired", name = "Species") +
                                    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
                                    ggplot2::theme_bw() +
                                    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
                                                axis.text.y = element_text(size = 12, color = "black"),
                                                    axis.title.y = element_text(size = 12, color = "black"),
                                                    legend.position = "right",
                                                    legend.text = element_text(size = 12, face = "italic", color = "black"),
                                                    legend.title = element_text(size = 12, color = "black"),
                                                    legend.key.size = unit(0.5, "lines"))

pyth_relativ_year_species_soil
png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Relative_abundance/pyth_relativ_year_species_soil..png", width = 3000, height = 4000, res = 600)
pyth_relativ_year_species_soil
dev.off()



#########################
#import and prepare OTUs#
#########################
pythium_species <- read.xlsx("C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Relative_abundance/Relative_abundance_pythium_top_10.xlsx", sheetIndex=1)
pyth_metadata <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/output/metadata_pythium.txt", header = TRUE, sep = "\t")


pythium_species <- merge(pyth_metadata, pythium_species, by.x = "Id", by.y="Id")
pythium_species <- pythium_species %>% filter(Year !="2022")
pythium_species <- pythium_species %>% filter(Compartment =="Soil")
pythium_species <- pythium_species %>% filter(Country !="Czech_republic")

pythium_species <- tidyr::gather(pythium_species, 
                                key = Taxa, 
                                value = Abundance,
                                s__sylvaticum:s__attrantheridium) 


pythium_species$Taxa <- str_replace(pythium_species$Taxa, "s__", "P. ")
pythium_species$Abundance <- as.numeric(pythium_species$Abundance)
pythium_species <- pythium_species %>% mutate(Abundance = Abundance*100)
pythium_species$Year <- as.factor(pythium_species$Year)


pythium_boxplot_species_year <- ggplot2::ggplot(pythium_species, aes(x = Year, y=Abundance, fill=Taxa))+
                                    ggplot2::geom_boxplot() +
                                    ggplot2::stat_summary(fun.y="mean", shape = 13, size = 0.4, color = "grey20") +
                                    ggplot2::facet_wrap(vars(Taxa), ncol = 5) +
                                    ggplot2::scale_fill_brewer(palette = "Paired", name = "Species") +
                                    ggplot2::labs(y = "Relative abundance (%)", x = "") +
                                    ggplot2::theme_bw() +
                                    ggplot2::theme(strip.background = element_rect(fill="white")) +
                                    ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
                                    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
                                                axis.text.y = element_text(size = 12, color = "black"),
                                                    axis.title.y = element_text(size = 12, color = "black"),
                                                    legend.position = "right",
                                                    legend.text = element_text(size = 12, face = "italic", color = "black"),
                                                    legend.title = element_text(size = 12, color = "black"),
                                                    legend.key.size = unit(0.5, "lines"))

pythium_boxplot_species_year






pyth_relative_abundance_barplot_boxplot <- ggarrange(pyth_relativ_year_species_soil, pythium_boxplot_species_year,
                                                labels =c("A", "B"),
                                                common.legend = TRUE,
                                                widths = c(25,75),
                                                legend="bottom", 
                                                ncol = 2)

pyth_relative_abundance_barplot_boxplot                                            


png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Relative_abundance/pyth_relativ_year_barplot_boxplot.png", width = 7000, height = 4000, res = 600)
pyth_relative_abundance_barplot_boxplot
dev.off()