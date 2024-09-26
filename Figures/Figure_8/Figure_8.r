library(ggplot2)
library(dplyr)
library(patchwork)
library(phyloseq)
library(ape)
library(stringr)
library(tidyverse)

set.seed(42)

# load phyloseq object
physeq_pythium <- readRDS("../Data/phyloseq_object.rds")
physeq_pythium <- phyloseq::subset_samples(physeq_pythium, Compartment == "Soil")

# helper functions
stat_sign <- function(pval) {
    symbol <- ifelse(pval <= 0.001, "***",
                     ifelse(pval > 0.001 & pval <= 0.01, "**",
                     ifelse(pval > 0.01 & pval <= 0.05, "*",
                     ifelse(pval > 0.05 & pval <= 0.1, ".", ""))))
    return(symbol)
}

# calculate relative abundance of all species
target_taxonomy <- physeq_pythium  %>%
                    tax_glom(taxrank = "Species") %>%
                    transform_sample_counts(function(x) x / sum(x)) %>%
                    psmelt() %>%
					dplyr::mutate(Organism = paste(Genus, Species, sep = " "))

# calculate pearson correlation per species
# do not consider zero abundance of species
build_correlation_data <- function(which_species) {
    dat <- target_taxonomy %>% dplyr::filter(Abundance != 0, Organism == which_species)
    props <- c("mpH_CaCl2", "mOC")
    dat <- data.frame(do.call("rbind", (setNames(lapply(props, function(x) {
                stats <- tryCatch(cor.test(formula = as.formula(paste0("~ Abundance +", x)),
                                           data = dat,
                                           use = "complete.obs",
                                           method = "pearson",
                                           alternative = "two.sided"),
                                  error = function(e) NULL)
                return(data.frame(species = which_species,
                                  cor = ifelse(length(stats) != 0, round(as.numeric(stats$estimate), 3), NA),
                                  pval = ifelse(length(stats) != 0, round(as.numeric(stats$p.value), 3), NA)))
    }), props))))
    return(dat %>%
            rownames_to_column(var = "property") %>%
            dplyr::mutate(direction = ifelse(cor < 0, "NEG", "POS")))
}

# build correlation data frame
corr_data <- do.call("rbind", lapply(unique(target_taxonomy$Organism), function(x) build_correlation_data(x)))
corr_data$property <- as.factor(ifelse(corr_data$property == "mOC", "OC", "pH"))

# create plot of correlation data
corr_plot <- function(df) {
    correlation <- ggplot(df, aes(x = property, y = species, fill = cor)) +
                     geom_tile(color = "gray80", lwd = 0.25, linetype = 1) +
                     geom_text(aes(label = stat_sign(pval))) +
                     scale_fill_gradientn(colours = c("red", "white", "blue"), na.value = "grey98", limits = c(-1, 1)) +
                     scale_x_discrete(labels = c("mpH_CaCl2" = "pH", "mOC" = "OC")) +
                     theme_classic(base_size = 12) +
                     coord_flip() +
                     labs(x = "Soil property", y = "Species", fill = "Pearson correlation") +
                     theme(legend.position = "top",
                           legend.key.width = unit(1, "cm"),
                           axis.text.x = element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5))
    return(correlation)
}

# create plot of significant intercations of soil property and species
corr_plot_sign <- ggplot(corr_data[corr_data$pval <= 0.05, ] %>% na.omit(), aes(x = species, y = pval)) +
                    geom_point(aes(color = direction, size = abs(cor))) +
                    labs(y = "p-value", x = "Species", color = "Direction", size = "Correlation") +
                    coord_flip() +
                    scale_y_continuous(limits = c(0, 0.06)) +
                    scale_color_manual(values = c("firebrick", "cornflowerblue")) +
                    facet_grid(. ~ property) +
                    theme_light(base_size = 14) +
                    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1, vjust = 1),
                          axis.text.y = element_text(face = "italic", size = 12),
                          panel.grid.minor = element_blank(),
                          strip.text = element_text(colour = "black"))


ggsave(plot = corr_plot(corr_data), file = "correlation_plot_all.jpg", width = 10, height = 4.2, dpi = 600)
ggsave(plot = corr_plot_sign, file = "correlation_plot_sign.jpg", width = 6.5, height = 3.5, dpi = 600)

ggsave(plot = corr_plot(corr_data), file = "correlation_plot_all.pdf", width = 10, height = 4)
ggsave(plot = corr_plot_sign, file = "correlation_plot_sign.pdf", width = 6.5, height = 3.5)
