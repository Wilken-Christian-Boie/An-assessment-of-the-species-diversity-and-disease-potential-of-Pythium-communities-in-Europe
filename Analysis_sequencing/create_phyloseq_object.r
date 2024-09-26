#!/usr/bin/Rscript --vanilla
library(phyloseq)
library(ape)
library(dplyr)
library(tibble)

# load all data of previous pipeline
# 1. Count table
# 2. Metadata
# 3. Taxonomy
# 4. Phylogenetic tree

# helper functions to properly import all required parts
ImportOTUTable <- function(filepath) {
    otu_import <- as.matrix(read.table(filepath, sep = "\t", header = TRUE))
    rownames(otu_import) <- paste0("ASV_", seq_len(nrow(otu_import)))
    class(otu_import) <- "numeric"
    otu_import <- phyloseq::otu_table(otu_import, taxa_are_rows = TRUE)
    return(otu_import)
}

ImportMetadata <- function(filepath) {
    sampledata <- read.table(filepath, sep = "\t", header = TRUE)
	soil_data <- read.table("./figures/Data/nearest_locations.tsv", sep = "\t", head = TRUE)
    colnames(sampledata) <- c("All", "Location", "Country", "Year", "Compartment", "Soilweight")
    sampledata <- dplyr::left_join(sampledata, soil_data, by = c("Location" = "id"))
    rownames(sampledata) <- sampledata$All
    return(sampledata)
}

ImportTaxonomy <- function(filepath, new_tax) {
	# import annotated taxonomy
    taxonomy <- read.table(filepath, header = FALSE, col.names = c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

	# import actual taxonomy based of Nguyen et al. 2022
	# 'Whole genome sequencing and phylogenomic analysis show support for the splitting of genus Pythium', Mycologia, Volume 114, 2022, Issue 3
	new_taxonomy <- read.table(new_tax, header = TRUE, sep = "\t")
	taxonomy <- taxonomy %>%
					dplyr::left_join(new_taxonomy, by = "Species") %>%
					dplyr::select(-Genus, -Genus_alt) %>%
					dplyr::rename(Genus = Genus_neu) %>%
                    dplyr::mutate(Organism = paste(Genus, Species, sep = " ")) %>%
					tibble::column_to_rownames(var = "ASV")

    taxonomy <- as.matrix(taxonomy[c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Clade", "Organism")])
    return(taxonomy)
}

ImportTree <- function(filepath) {
    return(read.tree(filepath))
}

# read all parts
otu_pythium <- ImportOTUTable("../output/ASVs_Pythium_counts.txt")
metadata <- ImportMetadata("../output/Metadata_pythium.txt")
tree_pythium <- ImportTree("../output/ASVs_Pythium.nwk")
taxa_df_pythium <- ImportTaxonomy("../output/ASVs_Pythium_taxonomy.txt", "./figures/Data/new_taxonomy_and_clades.txt")

# build phyloseq object
physeq_pythium <- phyloseq(phyloseq::sample_data(metadata),
                           otu_pythium,
                           tree_pythium,
                           phyloseq::tax_table(taxa_df_pythium))

# focus on Pythium sensu latu
pythium_sensu_latu <- c("Globisporangium", "Elongisporangium", "Phytopythium", "Pilasporangium", "Pythium")
physeq_pythium <- phyloseq::subset_taxa(physeq_pythium, Genus %in% pythium_sensu_latu)

# save phyloseq object as reference to create plots
saveRDS(physeq_pythium, file = "./figures/Data/phyloseq_object.rds")
