library(phyloseq)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggtext)

set.seed(42)

# load phyloseq object
physeq_pythium <- readRDS("../Data/phyloseq_object.rds")

# prepare data
relative_data <- function(tax_level) {
    target_taxonomy <- physeq_pythium  %>%
                            phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%
                            phyloseq::psmelt()
    return(target_taxonomy)
}

data <- relative_data() %>% dplyr::filter(!is.na(Abundance))

color_codes <- c("G. attrantheridium" = "#1F78B4",
                 "G. heterothallicum" = "#B2DF8A",
                 "G. sylvaticum" = "#E31A1C",
                 "P. aff. hydnosporum" = "#FF7F00",
                 "P. monospermum" = "#6A3D9A",
                 "P. arrhenomanes" = "#CAB2D6",
                 "G. rostratifingens" = "#FB9A99",
                 "G. ultimum var. ultimum" = "#FDBF6F",
                 "G. intermedium" = "#33A02C",
                 "G. apiculatum" = "#A6CEE3")

plot_species <- function(species_name) {
    dat <- data %>% dplyr::filter(Compartment == "Soil", Organism == species_name)
    dat <- dat %>% dplyr::filter(Abundance != "0")
    dat <- dat[-c(4:15, 17)]
    dat <- dat %>% group_by(OTU, Organism) %>% count() %>% ungroup() %>% mutate(Sample = `n`) %>% top_n(15, n) %>% arrange(desc(n)) %>% slice_max(n = 15, order_by = n)
    dat <- dat %>% arrange(desc(n)) %>% head(10)
	
	dat$Organism <- gsub("Pythium", "P.", dat$Organism)
	dat$Organism <- gsub("Globisporangium", "G.", dat$Organism)
	dat$Organism <- gsub("_", " ", dat$Organism)
	write.table(dat, file = gsub(" ", "_", paste0("../Data/", species_name, "_ASV_distribution.txt", sep = "")), col.names = TRUE, row.names = FALSE, quote = FALSE)

    abund_plot <- ggplot(dat, aes(x = reorder(OTU, -Sample), y = Sample, fill = Organism)) +
                    geom_bar(stat = "identity") +
                    scale_fill_manual(values = color_codes) +
                    coord_flip() +
                    facet_wrap(~ Organism) +
                    labs(y = "Number of Locations", x = "10 most common ASVs") +
                    theme_bw() +
                    theme(strip.background = element_rect(fill = "white"),
                          strip.text.x = element_text(size = 8, color = "black", face = "bold.italic"),
                          plot.title = element_textbox(hjust = 0.5, face = "italic"),
                          axis.text.y = element_text(size = 8, color = "black"),
                          axis.title.y = element_text(size = 8, color = "black"),
                          axis.title.x = element_text(size = 8, color = "black"),
                          legend.position = "none")
    return(abund_plot)
}

jpeg(file = "Figure 7.jpg", width = 7000, height = 4000, res = 600)
(plot_species("Globisporangium apiculatum") |
 plot_species("Globisporangium attrantheridium") |
 plot_species("Globisporangium heterothallicum") |
 plot_species("Globisporangium intermedium") |
 plot_species("Globisporangium rostratifingens")) /
(plot_species("Globisporangium sylvaticum") |
 plot_species("Globisporangium ultimum_var._ultimum") |
 plot_species("Pythium aff._hydnosporum") |
 plot_species("Pythium arrhenomanes") |
 plot_species("Pythium monospermum"))
dev.off()

pdf(file = "Figure 7.pdf", width = 14, height = 8)
(plot_species("Globisporangium apiculatum") |
 plot_species("Globisporangium attrantheridium") |
 plot_species("Globisporangium heterothallicum") |
 plot_species("Globisporangium intermedium") |
 plot_species("Globisporangium rostratifingens")) /
(plot_species("Globisporangium sylvaticum") |
 plot_species("Globisporangium ultimum_var._ultimum") |
 plot_species("Pythium aff._hydnosporum") |
 plot_species("Pythium arrhenomanes") |
 plot_species("Pythium monospermum"))
dev.off()

