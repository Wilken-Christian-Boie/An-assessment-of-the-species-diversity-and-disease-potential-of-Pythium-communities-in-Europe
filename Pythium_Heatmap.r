#rm(list=ls())
library(phyloseq)
library(gridGraphics)
library(dplyr)
library(ggplot2)
library(ape)
library(patchwork)
library(RColorBrewer)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(xlsx)
library(gridExtra)
library(grid)


set.seed(42)
setwd("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/scripts")

colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))


ImportOTUTable <- function(filepath){
    otu_import <- as.matrix(read.table(filepath, sep = "\t", header = TRUE))
    rownames(otu_import) <- paste0("ASV_",c(1:nrow(otu_import)))
    class(otu_import) <- "numeric"
    otu_import <- otu_table(otu_import, taxa_are_rows = TRUE)
    return(otu_import)
}
ImportMetadata <- function(filepath){
    sampledata <- read.table(filepath, sep = "\t", header = TRUE)
    colnames(sampledata) <- c("All", "Location", "Country", "Year","Town","Compartment","Soiltexture","Precrop","Tillage")
    rownames(sampledata) <- sampledata$All
    sampledata <- sampledata[grep("R1", rownames(sampledata)), ]
    return(sampledata)
}
ImportTaxonomy <- function(filepath){
    taxonomy <- read.table(filepath, header = FALSE)
    taxonomy <- taxonomy[, 2:ncol(taxonomy)]
    colnames(taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    rownames(taxonomy) <- paste0("ASV_",c(1:nrow(taxonomy)))
    taxonomy <- as.matrix(taxonomy[c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
    return(taxonomy)
}
ImportTree <- function(filepath){
    taxonomy_tree <- read.tree(filepath)
    return(taxonomy_tree)

}


otu_pythium <- ImportOTUTable("../output/Pythium/ASVs_Pythium_counts.txt")
metadata <- ImportMetadata("../output/metadata_pythium.txt")
tree_pythium <- ImportTree("../output/Pythium/ASVs_Pythium.nwk")
taxa_df_pythium <- ImportTaxonomy("../output/Pythium/ASVs_Pythium_taxonomy_blast.txt")
taxa_df_pythium <- gsub("s__", "P. ", taxa_df_pythium)

physeq_pythium = phyloseq(phyloseq::sample_data(metadata),
                          otu_pythium,
                          tree_pythium,
                          phyloseq::tax_table(taxa_df_pythium))
physeq_pythium


#############################################
#           Pythium Heatmap Soil            #
#############################################
    phylo_complex_heat <- function(tax_level, compartment, year){
        target_taxonomy <- physeq_pythium  %>%
                                phyloseq::tax_glom(taxrank = tax_level) %>%                                
                                phyloseq::subset_taxa(Genus %in% c("g__Pythium")) %>% 
                                phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%               
                                phyloseq::psmelt()                                                          
                            
        target_taxonomy <- target_taxonomy %>% filter(Compartment == compartment & Year == year)

        rel_abund <- target_taxonomy %>% mutate(organism = paste0(Species))  %>% select(Location, organism, Abundance)

        rel_abund_spread <- tidyr::spread(rel_abund, 
                                    key = organism, 
                                    value = Abundance) 
        rel_abund_spread[is.na(rel_abund_spread)] <- 0 
        cols_inorder <- unique(rel_abund_spread$Location[order(nchar(rel_abund_spread$Location), rel_abund_spread$Location)])



        
        rownames(rel_abund_spread) <- as.character(rel_abund_spread$Location)
        rel_abund_spread <- rel_abund_spread %>% dplyr::select(-c(Location)) 

        rel_abund_spread <- t(as.matrix(rel_abund_spread))

        class(rel_abund_spread) <- "numeric"
        
        col_fun = colorRamp2(c(100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,70,69,68,67,66,65,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,
                                49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0), 
                            c("#2F1E12","#301F12","#322013","#342113","#352214","#372314","#392415","#3A2516","#3C2616","#3E2717","#3F2817","#412918","#432A18","#442B19","#462C19",
                            "#482D1A","#4A2E1A","#4B2F1B","#4D301B","#4F311C","#50321C","#52331D","#54341D","#55351E","#57361E","#59371F","#5B381F","#5C3920","#5E3A20","#603B21","#623C21",
                            "#633D22","#653E22","#673F23","#684023","#6A4124","#6C4224","#6E4325","#6F4425","#714525","#734626","#754726","#774827","#784927","#7A4A28","#7C4B28","#7E4C28",
                            "#7F4D29","#814E29","#834F2A","#85502A","#87512B","#88522B","#8A532B","#8C542C","#8E552C","#90562D","#91572D","#93582D","#95592E","#975A2E","#995B2F","#9B5C2F",
                            "#9C5D2F","#9E5E30","#A05F30","#A26030","#A46131","#A66231","#A76332","#A96432","#AB6532","#AD6633","#AF6733","#B16833","#B36934","#B46A34","#B66B34","#B86C35",
                            "#BA6D35","#BE7039","#C3733D","#C77741","#CB7A45","#CF7D4A","#D2814E","#D68453","#D98758","#DD8B5D","#E08F62","#E39267","#E6966D","#E89A72","#EB9E78","#EDA27E",
                            "#EFA684","#F1AA8A","#F3AE91","#F5B397","#F7B79E","white"))

        col_fun2 = colorRamp2(c(0,20), c("white", "blue"))
        col_fun3 = colorRamp2(c(0,40), c("white", "darkgreen"))
        return(Heatmap(rel_abund_spread*100, 
                                column_title = paste0(compartment," ", year),
                                border = TRUE,
                                column_title_side = "top",
                                name = "Abundance (%)",
                                col = col_fun,
                                row_names_gp = gpar(fontface = "italic"),
                                row_order = order(as.numeric(gsub("P", "Pythium ", rownames(rel_abund_spread)))),
                                column_order = order(as.numeric(gsub("P", "", colnames(rel_abund_spread)))),
                                top_annotation = HeatmapAnnotation ("No. of Taxa" = rowSums(t(rel_abund_spread)!=0), col = list("No. of Taxa" = col_fun2), border = TRUE),
                                right_annotation = rowAnnotation("No. Loc." = colSums(t(rel_abund_spread)!=0), col = list("No. Loc." = col_fun3), border = TRUE)))
    }

png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Heatmap/Pythium/Abbildungen/2019/Pythium_Heatmap_soil_2019.png", width = 9000, height = 9000, res = 900)
phylo_complex_heat("Species", "Soil", "2019")
dev.off()

png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Heatmap/Pythium/Abbildungen/2020/Pythium_Heatmap_soil_2020.png", width = 9000, height = 9000, res = 900)
phylo_complex_heat("Species", "Soil", "2020")
dev.off()

png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Heatmap/Pythium/Abbildungen/2021/Pythium_Heatmap_soil_2021.png", width = 9000, height = 9000, res = 900)
phylo_complex_heat("Species", "Soil", "2021")
dev.off()


#################################################
#               Jahresvergleich                 #
#################################################
pythium_species <- read.xlsx("C:/Users/Wilken Boie/Desktop/Pythium/R/Heatmap/Pythium/Relative_abundanz_absolut.xlsx", sheetIndex=1)
pyth_metadata <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/output/metadata_pythium.txt", header = TRUE, sep = "\t")

pythium_species <- merge(pyth_metadata, pythium_species, by.x = "Id", by.y="Id")
pythium_species <- pythium_species %>% filter (Year != "2022")
pythium_species <- pythium_species %>% filter (Compartment == "Soil")


pythium_species <- pythium_species[,-c(1:3, 5:10)]

pythium_species <- tidyr::gather(pythium_species, 
                                key = Taxa, 
                                value = Abundance,
                                s__sylvaticum:s__paroecandrum) 


pythium_species_count <- pythium_species 
s__biforme <- c(2022, "s__biforme", 0.0000000000000000000001)
s__violae <- c(2022, "s__violae", 0.00000000000000000000001)
s__splendens <- c(2022, "s__splendens", 0.000000000000000001)
pythium_species_count <- rbind(pythium_species_count, s__biforme, s__violae, s__splendens)
pythium_species_count$Abundance <- as.numeric(pythium_species_count$Abundance)
pythium_species_count$Year <- as.integer(pythium_species_count$Year)

pythium_species_count <- pythium_species_count %>% filter (Abundance != "0")
pythium_species_count  <- pythium_species_count %>% group_by(Year, Taxa) %>% summarize(Abundance = n())
pythium_species_count$Abundance <- as.numeric(pythium_species_count$Abundance)
pythium_species_count <- tidyr::spread(pythium_species_count, 
                                Taxa, 
                                Abundance) 
pythium_species_count <- pythium_species_count %>% mutate_all(~replace_na(., 0))
pythium_species_count <- pythium_species_count [-4,]
rownames(pythium_species_count) <- as.character(pythium_species_count$Year)
pythium_species_count <- t(as.matrix(pythium_species_count))
pythium_species_count <- pythium_species_count[-1,]
species_totals <- rowSums(pythium_species_count)




pythium_species <- pythium_species %>% group_by(Year, Taxa) %>% summarise(Abundance=sum(Abundance)) %>% mutate(Abundance = Abundance / sum(Abundance))
cols_inorder <- unique(pythium_species$Taxa[order(nchar(pythium_species$Taxa), pythium_species$Taxa)])
pythium_species$Taxa <- str_replace(pythium_species$Taxa, "s__", "P. ")
pythium_species <- tidyr::spread(pythium_species, 
                                Taxa, 
                                Abundance) 
rownames(pythium_species) <- as.character(pythium_species$Year)
pythium_species <- t(as.matrix(pythium_species))
pythium_species <- pythium_species[-1,]
class(pythium_species) <- "numeric"


    col_fun = colorRamp2(c(100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,70,69,68,67,66,65,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,
                            49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0), 
                        c("#000608","#000A0D","#000E12","#001217","#00151C","#001921","#001C26","#001F2B","#002130","#002435","#00273A","#00293F","#002B44","#002D49","#002F4E",
                        "#003153","#003358","#00345D","#003663","#003768","#00386D","#003972","#003B77","#003C7C","#003C81","#003D86","#003E8B","#003F90","#003F95","#01409A","#02419F",
                        "#0341A4","#0441A9","#0642AE","#0742B3","#0943B8","#0B43BD","#0C43C1","#0E44C4","#1044C8","#1244CB","#1445CE","#1645D2","#1845D5","#1B46D8","#1D46DB","#1F47DE",
                        "#2247E1","#2447E3","#2748E6","#284AE8","#2A4CEB","#2B4FED","#2D51EF","#2E53F1","#3056F3","#3258F5","#345AF6","#355DF8","#375FFA","#3962FC","#3B64FD","#3D67FF",
                        "#3F6AFF","#416CFF","#446FFF","#4672FF","#4874FF","#4B77FF","#4D7AFF","#507DFF","#5480FF","#5882FF","#5C85FF","#6088FF","#648BFF","#678EFF","#6B91FF","#6F94FF",
                        "#7397FF","#779AFF","#7A9CFF","#7E9FFF","#82A2FF","#86A5FF","#8AA8FF","#8EABFF","#91AEFF","#95B1FF","#99B4FF","#9DB7FF","#A1BAFF","#A5BDFF","#A8C0FF","#ACC3FF",
                        "#B0C6FF","#B4C9FF","#B8CCFF","#BBCFFF","#BFD2FF","white"))



THRESH <- 42 #(P.intermedium) 
OFFSET <- 10
greater_top_ten_threshold <- which(unname(unlist(species_totals)) %in% c(THRESH:100))


png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Heatmap/Pythium/Abbildungen/Heatmap_compar_years.png", width = 4200, height = 9000, res = 900)
Heatmap(pythium_species*100, 
                            column_title = paste0("Year"),
                            border = TRUE,
                            column_title_side = "top",
                            name = "Abundance (%)",
                            col = col_fun,
                            row_names_gp = gpar(fontface = "italic"),
                            row_order = order(as.numeric(gsub("P", "Pythium", rownames(pythium_species)))),
                            column_order = order(as.numeric(gsub("Year", "", colnames(pythium_species)))),
                            right_annotation = rowAnnotation("No. Locations" = anno_barplot(species_totals,
                                                             border = TRUE,
                                                             ylim=c(0,100),
                                                             width = unit(2, "cm"),
                                                             gp = gpar(fill = ifelse(data.frame(species_totals)$species_totals >= THRESH, "darkslategray", "gray80")))))

                            decorate_annotation("No. Locations", grid.lines(unit(c(THRESH, THRESH), "native"), c(0, 1), gp = gpar(col = "red")))
                            

dev.off()

