# Relative_abundance_pythium
library(ggplot2)
library(xlsx)
library(dplyr)
library(tidyr)
library(ggpubr)


#########################
#import and prepare Data#
#########################
pythium_count_matrix <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/output/Pythium/ASVs_Pythium_counts.txt", , header = TRUE, sep = "\t")

pythium_count_matrix <- tidyr::gather(pythium_count_matrix, 
                                key = Id, 
                                value = Count,
                                pmOomyceten2_BD_S40_L001_adaptertrimmed_Pythium_2020_R1.fastq:year2021_soil_P9_adaptertrimmed_Pythium_2021_R1.fastq) 


colnames(pythium_count_matrix)[1] = "ASV"


pythium_taxa <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Bubble_plot/ASVs_taxonomy_blast.txt", , header = TRUE, sep = "\t")
pythium_taxa <- pythium_taxa %>% filter(Genus == "g__Pythium")
pythium_taxa <- pythium_taxa %>% filter(Species %in% c("s__intermedium", "s__aristosporum", "s__attrantheridium", "s__heterothallicum", "s__monospermum", "s__oligandrum", "s__sp", "s__sylvaticum", "s__ultimum", "s__rostratifingens"))
pythium_taxa$Species <- gsub("s__", "P. ", pythium_taxa$Species)
pythium_count_taxa <- merge(pythium_count_matrix, pythium_taxa, by.x = "ASV", by.y="ASV")



pyth_metadata <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/output/metadata_pythium.txt", header = TRUE, sep = "\t")
pyth_metadata <- pyth_metadata %>% filter(Year !="2022")

pythium_species <- merge(pyth_metadata, pythium_count_taxa, , by.x = "Id", by.y="Id")
pythium_species <- pythium_species %>% filter(Compartment == "Soil")
pythium_species$Year <- as.factor(pythium_species$Year)
pythium_species <- subset(pythium_species, Count >= 1)
pythium_species <- pythium_species %>% mutate_if(is.numeric, ~1 * (. > 0))


#########
#plotten#
#########

pythium_species_für_df1 <- pythium_species %>% group_by(ASV, Species, Year) %>% dplyr::summarise(count = n())
df1 <- pythium_species_für_df1 %>% group_by(Species, Year) %>% dplyr::summarise(count = n())


bubble_plot1 <- ggplot(df1, aes(x=Species, y=Year, size = count, fill = Species)) +
                geom_point(shape = 21) +
                scale_fill_brewer(palette = "Spectral", name = "Species") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0, size = 15, color = "black"),
                            axis.text.y = element_text(size = 12, color = "black"),
                                   axis.title.y = element_text(size = 12, color = "black"),
                                   legend.position = "right",
                                   legend.text = element_text(size = 12, color = "black"),
                                   legend.title = element_text(size = 12, color = "black"),
                                   legend.key.size = unit(0.5, "lines"))


bubble_plot1

pythium_species_für_df2 <- pythium_species %>% group_by(ASV, Species, Compartment) %>% dplyr::summarise(count = n())
df2 <- pythium_species_für_df2 %>% group_by(Species, Compartment) %>% dplyr::summarise(count = n())

bubble_plot2 <- ggplot(df2, aes(x=Species, y=Compartment, size = count, fill = Species)) +
                geom_point(shape = 21) +
                scale_fill_brewer(palette = "Spectral", name = "Species") +
                theme_bw() +
                theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.title.x = element_blank(),
                            axis.text.y = element_text(size = 12, color = "black"),
                                   axis.title.y = element_text(size = 12, color = "black"),
                                   legend.text = element_text(size = 12, color = "black"),
                                   legend.title = element_text(size = 12, color = "black"),
                                   legend.key.size = unit(0.5, "lines"))


bubble_plot2

pyth_bubble <- ggarrange(bubble_plot2, bubble_plot1,
                                            labels =c("A", "B"), 
                                            legend="bottom",
                                            heights =c(1,3), 
                                            ncol = 1)
pyth_bubble




#########################################################
#                    plot all in one                    #
#########################################################


df2[df2 == "Soil"] <- "All years"
colnames(df2)[2] = "Year"
df3 <- bind_rows(df1, df2)



bubble_plot3 <- ggplot(df3, aes(x=Species, y=Year, size = count, fill = Species)) +
                geom_point(shape = 21) +
                scale_size(range = c(1, 30), name = "Number of ASVs", breaks = c(5, 10, 30, 50, 70)) +
                scale_fill_brewer(palette = "Paired", name = "Species") +
                theme_bw() +
                theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                geom_hline(yintercept = 3.5, size = 1) +
                scale_x_discrete(name = "") +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25, size = 15, face = "italic", color = "black"),
                            axis.text.y = element_text(size = 15, color ="black"),
                                   axis.title.y = element_text(size = 15, color = "black"),
                                   legend.position = "right",
                                   legend.text = element_text(face = "italic", size = 15, color = "black"),
                                   legend.title = element_text(size = 15, color = "black"),
                                   legend.key.size = unit(1, "lines"))

bubble_plot3



png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Bubble_plot/Bubble_plot_soil.png", width = 5000, height = 3700, res = 400)
bubble_plot3
dev.off()


