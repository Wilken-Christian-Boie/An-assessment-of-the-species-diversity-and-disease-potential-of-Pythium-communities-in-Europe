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
library(gridExtra)
library(grid)

set.seed(42)

# load phyloseq object
physeq_pythium <- readRDS("../Data/phyloseq_object.rds")

# define clades and get row order of heatmap
clades <- read.table("../../figures/Data/new_taxonomy_and_clades.txt", head = TRUE, sep = "\t")
clades$tax <- paste(clades$Genus_neu, clades$Species, sep = " ")
clades$tax <- gsub("_", " ", clades$tax)
clades <- clades %>% dplyr::select("tax", "Clade") %>% dplyr::arrange(Clade) %>% column_to_rownames(var = "tax")

# Pythium heatmap soil for all locations
phylo_complex_heat <- function(compartment, year, output = "plot") {
    # prepare the data to get relative abundance values
    target_taxonomy <- physeq_pythium  %>%
                            phyloseq::tax_glom(taxrank = "Species") %>%
                            phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%
                            phyloseq::psmelt()

    # focus on specific year and compartment
    target_taxonomy <- target_taxonomy %>%
                        dplyr::filter(Compartment == compartment & Year == year) %>%
                   		dplyr::mutate(Organism = paste(Genus, Species, sep = " "))

    # select relevance data out of all metadata
    rel_abund <- target_taxonomy %>% dplyr::select(Location, Organism, Abundance)

    # change data format from long to wide format
    rel_abund_spread <- tidyr::spread(data = rel_abund,
                                      key = Organism,
                                      value = Abundance)
    rel_abund_spread[is.na(rel_abund_spread)] <- 0

    # define order of columns for heatmap plot
    cols_inorder <- unique(rel_abund_spread$Location[order(nchar(rel_abund_spread$Location), rel_abund_spread$Location)])

    # prepare data and adjust format for heatmap
    rownames(rel_abund_spread) <- as.character(rel_abund_spread$Location)
    rel_abund_spread <- rel_abund_spread %>% dplyr::select(-c(Location))
    rel_abund_spread <- t(as.matrix(rel_abund_spread))
	class(rel_abund_spread) <- "numeric"

    # define color ramps for heatmap
    col_fun <- colorRamp2(seq(100, 0, -1),
                          c("#2F1E12", "#301F12", "#322013", "#342113", "#352214", "#372314", "#392415", "#3A2516", "#3C2616", "#3E2717", "#3F2817", "#412918", "#432A18", "#442B19",
                            "#462C19", "#482D1A", "#4A2E1A", "#4B2F1B", "#4D301B", "#4F311C", "#50321C", "#52331D", "#54341D", "#55351E", "#57361E", "#59371F", "#5B381F", "#5C3920",
                            "#5E3A20", "#603B21", "#623C21", "#633D22", "#653E22", "#673F23", "#684023", "#6A4124", "#6C4224", "#6E4325", "#6F4425", "#714525", "#734626", "#754726",
                            "#774827", "#784927", "#7A4A28", "#7C4B28", "#7E4C28", "#7F4D29", "#814E29", "#834F2A", "#85502A", "#87512B", "#88522B", "#8A532B", "#8C542C", "#8E552C",
                            "#90562D", "#91572D", "#93582D", "#95592E", "#975A2E", "#995B2F", "#9B5C2F", "#9C5D2F", "#9E5E30", "#A05F30", "#A26030", "#A46131", "#A66231", "#A76332",
                            "#A96432", "#AB6532", "#AD6633", "#AF6733", "#B16833", "#B36934", "#B46A34", "#B66B34", "#B86C35", "#BA6D35", "#BE7039", "#C3733D", "#C77741", "#CB7A45",
                            "#CF7D4A", "#D2814E", "#D68453", "#D98758", "#DD8B5D", "#E08F62", "#E39267", "#E6966D", "#E89A72", "#EB9E78", "#EDA27E", "#EFA684", "#F1AA8A", "#F3AE91",
                            "#F5B397", "#F7B79E", "white"))

    col_fun2 <- colorRamp2(c(0, 20), c("white", "blue"))
    col_fun3 <- colorRamp2(c(0, 40), c("white", "darkgreen"))

    # check if heatmap or data is requested
    if (output == "data")
        return(data.frame(rel_abund_spread) %>%
                    tibble::rownames_to_column(var = "Species") %>%
                    dplyr::mutate(Year = as.factor(year)) %>%
                    tidyr::gather(key = "Location", value = "Abundance", -Year, -Species))

    if (output == "order") return(gsub("_", " ", rownames(rel_abund_spread)))
    if (output == "plot") {
		rownames(rel_abund_spread) <- gsub("_", " ", rownames(rel_abund_spread))
        return(Heatmap(rel_abund_spread * 100,
                        column_title = year,
                        border = TRUE,
                        column_title_side = "top",
                        name = "Abundance (%)",
                        col = col_fun,
                        row_names_side = "left",
                        row_names_gp = gpar(fontface = "italic", fontsize = 10),
                        column_title_gp = gpar(fontsize = 34),
                        #row_order = order(as.numeric(rownames(rel_abund_spread))),
                        row_order = rownames(clades),
                        column_order = order(as.numeric(colnames(rel_abund_spread))),
                        top_annotation = HeatmapAnnotation("No. of Taxa" = rowSums(t(rel_abund_spread) != 0),
                                                            col = list("No. of Taxa" = col_fun2),
                                                            border = TRUE),
                        right_annotation = rowAnnotation("No. Loc." = colSums(t(rel_abund_spread) != 0),
                                                        col = list("No. Loc." = col_fun3),
                                                        border = TRUE)))
        }
}

# build heatmaps of every year
heatmap_2019_plot <- phylo_complex_heat("Soil", "2019", output = "plot")
heatmap_2020_plot <- phylo_complex_heat("Soil", "2020", output = "plot")
heatmap_2021_plot <- phylo_complex_heat("Soil", "2021", output = "plot")

# extract data and build heatmaps summarizing years
year_2019_data_soil <- phylo_complex_heat("Soil", "2019", output = "data")
year_2020_data_soil <- phylo_complex_heat("Soil", "2020", output = "data")
year_2021_data_soil <- phylo_complex_heat("Soil", "2021", output = "data")

year_2019_data_bait <- phylo_complex_heat("Bait", "2019", output = "data")
year_2020_data_bait <- phylo_complex_heat("Bait", "2020", output = "data")
year_2021_data_bait <- phylo_complex_heat("Bait", "2021", output = "data")
heatmap_row_order <- phylo_complex_heat("Bait", "2021", output = "order")
heatmap_row_order2 <- phylo_complex_heat("Soil", "2021", output = "order")

###############################################
# Export tables as foundation of further plots#
###############################################

export_data <- function(df19, df20, df21, compartment) {
    # determine occurence per species over all years
    # in order to finally get the top 10 species
    species_occurrence_all_years <- rbind(df19, df20, df21) %>%
                                        dplyr::group_by(Species) %>%
                                        dplyr::summarize(Count = sum(Abundance != 0), .groups = "drop")

    species_occurrence_per_years <- rbind(df19, df20, df21) %>%
                                        dplyr::group_by(Year, Species) %>%
                                        dplyr::summarize(Count = sum(Abundance != 0), .groups = "drop")

    # get the top 10 species count and names
    top10_species_counts_all_year <- species_occurrence_all_years[order(species_occurrence_all_years$Count, decreasing = TRUE), ][1:11, ]
    top10_names <- top10_species_counts_all_year$Species
    write.table(top10_names, file = paste0("../Data/top10_names_", compartment, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(top10_species_counts_all_year, file = paste0("../Data/top10_species_counts_all_year_", compartment, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    # get the top 10 species count per year
    top10_species_counts_per_year <- species_occurrence_per_years %>% dplyr::group_by(Year) %>% dplyr::slice_max(order_by = Count, n = 11, with_ties = FALSE)
    write.table(top10_species_counts_per_year, file = paste0("../Data/top10_species_counts_per_year_", compartment, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    # get total number of species detected
    write.table(unique(species_occurrence_all_years$Species), file = paste0("../Data/total_species_names_", compartment, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    # get total number of species detected per year
    write.table(species_occurrence_per_years %>% dplyr::group_by(Year) %>% dplyr::summarize(Total_species = sum(Count != 0)),
                file = paste0("../Data/species_per_year_", compartment, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    # get species count per location and year
    write.table(rbind(df19, df20, df21) %>%
                    dplyr::group_by(Location, Year) %>%
                    dplyr::summarize(Species = sum(Abundance != 0)) %>%
                    dplyr::arrange(desc(Year), desc(Species)),
                file = paste0("../Data/species_count_per_location_and_year_", compartment, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

export_data(year_2019_data_soil, year_2020_data_soil, year_2021_data_soil, "soil")
export_data(year_2019_data_bait, year_2020_data_bait, year_2021_data_bait, "bait")


composite_heatmap <- function(df19, df20, df21, compartment = "soil", output = "plot", col_angle = 0) {
    top10_names <- read.table(paste0("../Data/top10_names_", compartment, ".txt"), sep = "\t", head = FALSE, col.names = "Species")$Species
	top10_names <- top10_names[!grepl("nov", top10_names)]

    # combine data of all years
    years_concatenated <- rbind(df19, df20, df21) %>%
                            dplyr::group_by(Species, Year) %>%
                            dplyr::summarize(Abundance = mean(replace_na(Abundance, 0)))

    species_occurrence_all_years <- rbind(df19, df20, df21) %>%
                                    dplyr::group_by(Species) %>%
                                    dplyr::summarize(Count = sum(Abundance != 0))
    # check if only data is requested
    if (output == "data") return(years_concatenated)

    # prepare data to plot heatmap
    years_concatenated_matrix <- tidyr::spread(data = years_concatenated,
                                               key = Species,
                                               value = Abundance) %>%
                                data.frame()

    # prepare data and adjust format for heatmap
    rownames(years_concatenated_matrix) <- as.character(years_concatenated_matrix$Year)
    colnames(years_concatenated_matrix) <- gsub("\\.", " ", colnames(years_concatenated_matrix))
    years_concatenated_matrix <- t(as.matrix(years_concatenated_matrix))[-1, ]
    class(years_concatenated_matrix) <- "numeric"

    col_fun <- colorRamp2(seq(100, 0, -1),
                          c("#000608", "#000A0D", "#000E12", "#001217", "#00151C", "#001921", "#001C26", "#001F2B", "#002130", "#002435", "#00273A", "#00293F", "#002B44", "#002D49", "#002F4E",
                            "#003153", "#003358", "#00345D", "#003663", "#003768", "#00386D", "#003972", "#003B77", "#003C7C", "#003C81", "#003D86", "#003E8B", "#003F90", "#003F95", "#01409A", "#02419F",
                            "#0341A4", "#0441A9", "#0642AE", "#0742B3", "#0943B8", "#0B43BD", "#0C43C1", "#0E44C4", "#1044C8", "#1244CB", "#1445CE", "#1645D2", "#1845D5", "#1B46D8", "#1D46DB", "#1F47DE",
                            "#2247E1", "#2447E3", "#2748E6", "#284AE8", "#2A4CEB", "#2B4FED", "#2D51EF", "#2E53F1", "#3056F3", "#3258F5", "#345AF6", "#355DF8", "#375FFA", "#3962FC", "#3B64FD", "#3D67FF",
                            "#3F6AFF", "#416CFF", "#446FFF", "#4672FF", "#4874FF", "#4B77FF", "#4D7AFF", "#507DFF", "#5480FF", "#5882FF", "#5C85FF", "#6088FF", "#648BFF", "#678EFF", "#6B91FF", "#6F94FF",
                            "#7397FF", "#779AFF", "#7A9CFF", "#7E9FFF", "#82A2FF", "#86A5FF", "#8AA8FF", "#8EABFF", "#91AEFF", "#95B1FF", "#99B4FF", "#9DB7FF", "#A1BAFF", "#A5BDFF", "#A8C0FF", "#ACC3FF",
                            "#B0C6FF", "#B4C9FF", "#B8CCFF", "#BBCFFF", "#BFD2FF", "white"))

	rownames(years_concatenated_matrix) <- gsub("_", " ", rownames(years_concatenated_matrix))
    heatmap_years_plot <- Heatmap(years_concatenated_matrix * 100,
                            column_title = paste0("Year"),
                            width = unit(1.5, "cm"),
                            border = TRUE,
                            column_title_side = "top",
                            name = "Abundance (%) ",
                            col = col_fun,
                            column_names_rot = col_angle,
                            row_names_side = "left",
                            row_names_gp = gpar(fontface = "italic"),
                            heatmap_legend_param = list(direction = "vertical"),
                            row_order = order(as.numeric(rownames(years_concatenated_matrix))),
                            column_order = order(as.numeric(colnames(years_concatenated_matrix))),
                            right_annotation = rowAnnotation("No. Loc" = anno_barplot(species_occurrence_all_years$Count,
                                                                border = TRUE,
                                                                width = unit(1.5, "cm"),
                                                                gp = gpar(fill = ifelse(species_occurrence_all_years$Species %in% top10_names,
                                                                                        "darkslategray",
                                                                                        "gray80")))))
    return(heatmap_years_plot)
}

# create composite heatmap (2019 + 2020 + 2021 + all years)
heat_soil_plot <- composite_heatmap(year_2019_data_soil, year_2020_data_soil, year_2021_data_soil, compartment = "soil")
heat_bait_plot <- composite_heatmap(year_2019_data_bait, year_2020_data_bait, year_2021_data_bait, compartment = "bait")

heat_soil_data <- composite_heatmap(year_2019_data_soil, year_2020_data_soil, year_2021_data_soil, compartment = "soil", output = "data")
heat_bait_data <- composite_heatmap(year_2019_data_bait, year_2020_data_bait, year_2021_data_bait, compartment = "bait", output = "data")


clades <- clades[match(heatmap_row_order2, rownames(clades)), ]
year_heatmap_soil <- heatmap_2019_plot + heatmap_2020_plot + heatmap_2021_plot +
					 rowAnnotation("Clade" = as.matrix(clades),
								   border = TRUE, 
								   col = list(Clade = c("A" = "#A6CEE3", "B" = "#1F78B4", "C" = "#B2DF8A", "D" = "#33A02C", "E" = "#FB9A99", "F" = "#E31A1C", "G" = "#FDBF6F", "I" = "#FF7F00", "J" = "#CAB2D6", "K" = "#6A3D9A", "NA" = "gray")))

jpeg(file = "Figure 2.jpg", width = 8000, height = 4000, res = 300)
draw(year_heatmap_soil, legend_grouping = "original", merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf(file = "Figure 2.pdf", width = 24, height = 13)
draw(year_heatmap_soil, legend_grouping = "original", merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# save composite heatmap of bait and soil samples
jpeg(file = "Figure 4.jpg", width = 8000, height = 4000, res = 300)
draw(heat_bait_plot + heat_soil_plot, legend_grouping = "original", merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

pdf(file = "Figure 4.pdf", width = 19, height = 14)
draw(heat_bait_plot + heat_soil_plot, legend_grouping = "original", merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# calculate abundance difference per species and all years considering bait and soil samples
difference_bait_soil_year <- rbind(heat_soil_data %>% dplyr::mutate(Compartment = "Soil"),
                                   heat_bait_data %>% dplyr::mutate(Compartment = "Bait")) %>%
                                dplyr::group_by(Species, Compartment) %>%
                                dplyr::summarize(mabund = sum(Abundance), .groups = "drop")

comparison_species_all_years_plot <- ggplot(difference_bait_soil_year, aes(x = factor(Species, levels = rev(heatmap_row_order)), y = mabund * 100, fill = Compartment)) +
                                        geom_bar(stat = "identity") +
                                        coord_flip() +
                                        scale_fill_manual(values = c("navajowhite", "#954535")) +
                                        labs(x = "Species", y = "Relative abundance all years (%)") +
                                        facet_wrap(. ~ Compartment) +
                                        theme_light(base_size = 12) +
                                        theme(legend.position = "none")

## prepare data to analyze top10 and all pythium species by relative and absolute abundance
top10_names <- read.table("../Data/top10_names_soil.txt", head = FALSE, sep = "\t", col.names = c("Species"))
metadata <- read.table("../../../output/Metadata_pythium.txt", head = TRUE, sep = "\t")

export_species_table <- function(relative = TRUE, top10 = TRUE, add_meta = FALSE) {
    species_abundance <- physeq_pythium %>%
                            phyloseq::subset_samples(Compartment == "Soil") %>%
                            phyloseq::tax_glom(taxrank = "Species")
    if (relative) {
        species_abundance <- species_abundance %>% phyloseq::transform_sample_counts(function(x) x / sum(x))
    }
    if (top10) {
		species_abundance <- species_abundance %>% phyloseq::subset_taxa(Species %in% sapply(strsplit(top10_names$Species, " "), `[`,2))
    }

    count_table <- data.frame(phyloseq::otu_table(species_abundance)) %>%
                        tibble::rownames_to_column(var = "ASV")

    asv_species_link <- data.frame(phyloseq::tax_table(species_abundance)) %>%
                            tibble::rownames_to_column(var = "ASV") %>%
                   			dplyr::mutate(Organism = paste(Genus, Species, sep = " ")) %>%
                            dplyr::select(ASV, Organism)

    count_table$ASV <- asv_species_link$Organism[match(count_table$ASV, asv_species_link$ASV)]
    colnames(count_table)[which(names(count_table) == "ASV")] <- "Id"
    count_table <- setNames(as.data.frame(t(count_table[-1])), count_table[[1]]) # maintain data structure during transpose
    count_table <- count_table %>% tibble::rownames_to_column(var = "Id")

    if (add_meta) {
        count_table <- count_table %>% dplyr::left_join(metadata, by = c("Id" = "Id"))
    }
    return(count_table)
}

write.table(export_species_table(relative = TRUE, top10 = TRUE), file = "../Data/Relative_abundance_top_10_Pythium_species.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(export_species_table(relative = FALSE, top10 = FALSE), file = "../Data/ASV_table_pythium_species.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(export_species_table(relative = TRUE, top10 = FALSE, add_meta = TRUE), file = "../Data/Pythium_matrix_relativ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
