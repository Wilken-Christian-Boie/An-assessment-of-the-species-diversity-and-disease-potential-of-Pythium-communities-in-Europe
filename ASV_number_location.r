#ASV_Location
#packages
library(phyloseq)
library(gridGraphics)
library(dplyr)
library(ggplot2)
library(ggtext)
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
library(cowplot)

#set seed for reproducible randomization
set.seed(42)
setwd("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/scripts")

colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))

#import functions for creating phyloseq object  
ImportOTUTable <- function(filepath){
    #delete starting "#OTU ID" from original OTU table (coming out of USEARCH) otherwise import will crush
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
    #import phylogenetic tree build with mafft and FastTree
    taxonomy_tree <- read.tree(filepath)
    return(taxonomy_tree)

}

# Finally import all necessary data for analysis
otu_pythium <- ImportOTUTable("../output/Pythium/ASVs_Pythium_counts.txt")
metadata <- ImportMetadata("../output/metadata_pythium.txt")
tree_pythium <- ImportTree("../output/Pythium/ASVs_Pythium.nwk")
taxa_df_pythium <- ImportTaxonomy("C:/Users/Wilken Boie/Desktop/Pythium/NGS/Workingstation/output/Pythium/ASVs_Pythium_taxonomy_blast.txt")
taxa_df_pythium <- gsub("s__", "P. ", taxa_df_pythium)
# TODO: edit sample names in "run_dada2.r" script -> no '.fastq' in name!
# TODO: remove first line of taxonomy file

# combine all imported data into phyloseq object
physeq_pythium = phyloseq(phyloseq::sample_data(metadata),
                          otu_pythium,
                          tree_pythium,
                          phyloseq::tax_table(taxa_df_pythium))
physeq_pythium



#############################################
#           Pythium prepare data            #
#############################################

relative_data <- function(tax_level){
    target_taxonomy <- physeq_pythium  %>%
                            phyloseq::subset_taxa(Genus %in% c("g__Pythium")) %>% # exlude positiv control -> Saprolegnia ferax
                            #phyloseq::tax_glom(taxrank = tax_level) %>%                                 # glomerate on specific level
                            phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%               # transform to relative counts
                            phyloseq::psmelt()                                                          # melt to regular data frame
                        
    rel_abund <- target_taxonomy %>% mutate(organism = paste0(Genus,"_",Species))  #%>% select(Location, organism, Abundance)

    return(rel_abund)

}
d <- relative_data() %>% filter(!is.na(Abundance))





##############################
#       P. aristosporum      #
##############################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. aristosporum")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_aristosporum <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")

#################################
#       P. attrantheridium      #
#################################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. attrantheridium")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_attrantheridium <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")

#################################
#       P. heterothallicum      #
#################################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. heterothallicum")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_heterothallicum <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")


#############################
#       P. intermedium      #
#############################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. intermedium")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_intermedium <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "10 most common ASVs") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")


#############################
#       P. monospermum      #
#############################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. monospermum")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_monospermum <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")

############################
#       P. oligandrum      #
############################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. oligandrum")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_oligandrum <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "10 most common ASVs") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")


#################################
#       P. rostratifingens      #
#################################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. rostratifingens")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_rostratifingens <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")

####################
#       P. sp      #
####################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. sp")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_sp <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")

############################
#       P. sylvaticum      #
############################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. sylvaticum")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)

P_sylvaticum <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")

#########################
#       P. ultimum      #
#########################
dat <- d %>% dplyr::filter(Compartment == "Soil", Species == "P. ultimum")
dat <- dat %>% filter (Abundance != "0")
dat <- dat [-c(3:20,22)]
dat <- dat %>% group_by(OTU, Species) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n =15, order_by =n)
dat <- dat %>% arrange(desc(n)) %>% head(10)




P_ultimum <- ggplot2::ggplot(dat, aes(x = reorder(OTU, -Sample), y=Sample, fill = Species)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_manual(values = c("P. aristosporum" = "#a6cee3", "P. attrantheridium" = "#1f78b4", "P. heterothallicum" = "#b2df8a", "P. intermedium" = "#33a02c", "P. monospermum" = "#fb9a99", 
                                                    "P. oligandrum" = "#e31a1c", "P. rostratifingens" = "#fdbf6f", "P. sp" = "#ff7f00", "P. sylvaticum" = "#cab2d6", "P. ultimum" = "#6a3d9a")) +
            ggplot2::coord_flip() +
            ggplot2::facet_wrap(~ Species) +
            ggplot2::labs(y = "Number of Locations", x = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = element_rect(fill="white")) +
            ggplot2::theme(strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black"),
                                                plot.title = element_textbox(hjust = 0.5, face = "italic"),
                                                axis.text.y = element_text(size = 8, color = "black"),
                                                axis.title.y = element_text(size = 8, color = "black"),
                                                axis.title.x = element_text(size = 8, color = "black"),
                                                        legend.position = "noen")


#########################
#   Plot_all_species    #
#########################
ASVs_of_top_10_abundance <- cowplot::plot_grid(P_aristosporum, P_attrantheridium, P_heterothallicum, P_intermedium, P_monospermum, 
                                                    P_oligandrum, P_rostratifingens, P_sp, P_sylvaticum, P_ultimum,
                                                    ncol = 5,
                                                    nrow = 2)

ASVs_of_top_10_abundance


png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/ASV Barplot/ASVs_of_top_10_abundance.png", width = 7000, height = 4000, res = 600)
ASVs_of_top_10_abundance
dev.off()