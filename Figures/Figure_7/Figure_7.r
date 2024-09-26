library(cluster)
library(reshape2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)

set.seed(42)

# import data
correlation <- read.table("../../../3_cooccurrence/correlation_fastspar.tsv", sep = "\t", head = FALSE)
pvalue <- read.table("../../../3_cooccurrence/pvalues_fastspar.tsv", sep = "\t", head = FALSE)
top10_names <- read.table("../Data/top10_names_soil.txt", sep = "\t", head = FALSE, col.names = "Species")
tax_names <- setNames(read.table("../../../3_cooccurrence/ASVs_Pythium_counts_summarized_name.tsv", head = FALSE, sep = "\t"), "taxonomy")

# prepare data
correlation <- as.matrix(correlation[, 2:ncol(correlation)])
rownames(correlation) <- tax_names$taxonomy
colnames(correlation) <- tax_names$taxonomy

pvalue <- as.matrix(pvalue[, 2:ncol(pvalue)])
class(correlation) <- "numeric"
class(pvalue) <- "numeric"

# select only significant correlation and reset self correlation (1)
correlation[pvalue < 0.05] <- 0
correlation[correlation == 1] <- 0

# helper functions
get_upper_tri <- function(cormat) {
        cormat[lower.tri(cormat)] <- NA
        return(cormat)
}

reorder_cormat <- function(cormat) {
        dd <- as.dist((1 - cormat) / 2)
        hc <- hclust(dd)
        cormat <- cormat[hc$order, hc$order]
        return(cormat)
}

top10 <- function(df) {
      return(df[df$Var1 %in% top10_names$Species & df$Var2 %in% top10_names$Species, ])
}

# rearrange correlation matrix as precondition for correct visualization
melted_cormat <- melt(get_upper_tri(reorder_cormat(correlation)), na.rm = TRUE)

# create correlation heatmap (upper triangle) of the correlation matrix
correlation_heatmap <- function(df, bsize = 6, text = FALSE, title_ = "", show_legend = TRUE, top10 = FALSE) {
      p <- ggplot(data = df, aes(Var2, Var1, fill = value)) +
            geom_tile(color = "gray80") +
            scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-0.16, 0.16), na.value = "white", space = "Lab", name = "Pearson\ncorrelation") +
            theme(axis.text.x = element_text(face = "italic", angle = 90, vjust = 0.5, size = bsize, hjust = 1),
                  axis.text.y = element_text(face = "italic", size = bsize)) +
            coord_fixed() +
            labs(title = title_) +
            scale_y_discrete(position = "right") +
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
                  axis.ticks = element_blank(),
                  legend.justification = c(1, 0),
                  legend.position = c(0.8, 0.15),
                  legend.direction = "horizontal") +
            guides(fill = guide_colorbar(barwidth = 15, barheight = 1.5,
                   title.position = "top", title.hjust = 0.5))
      if (text) p <- p + geom_text(aes(label = value))
      if (!show_legend) p <- p + theme(legend.position = "none")
	  if (top10) {
            p <- p + scale_y_discrete(position = "right", labels = function(x) gsub("Globisporangium", "G. ", gsub("Pythium", "P. ", x, fixed = TRUE), x, fixed = TRUE))
            p <- p + scale_x_discrete(labels = function(x) gsub("Globisporangium", "G. ", gsub("Pythium", "P. ", x, fixed = TRUE), x, fixed = TRUE))
	  }
      return(p)
}

# in order to summarize the data
melted_cormat$direction <- ifelse(melted_cormat$value > 0, "POS", "NEG")
cormat_summary <- melted_cormat %>%
                    dplyr::group_by(Var1, direction) %>%
                    dplyr::summarize(count = n(), .groups = "drop")

totals <- cormat_summary %>% dplyr::group_by(Var1) %>% dplyr::summarize(totals = sum(count))

# plot summary of all correlations of every species
summ <- ggplot(cormat_summary, aes(x = reorder(Var1, count), y = count, fill = direction)) +
            geom_bar(position = "stack", stat = "identity") +
            geom_text(aes(label = count), size = 2, position = "stack", hjust = 0.5, vjust = 1, color = "white") +
            geom_text(aes(label = ifelse(Var1 %in% top10_names$Species, "x", ""), y = 40), size = 4, hjust = 0.5, vjust = -2, color = "black") +
            scale_fill_manual(values = c("gold3", "darkslategray")) +
            scale_y_continuous(expand = c(0, 0)) +
            labs(x = "Species", y = "No. of significant correlations", fill = "") +
            theme_classic(base_size = 12) +
            theme(axis.text.y = element_text(size = 12, hjust = 0),
                  axis.text.x = element_text(face = "italic", size = 8, angle = 90, hjust = 1, vjust = 0.25),
                  legend.position = c(0.7, 0.9),
                  legend.direction = "horizontal",
                  plot.margin = unit(c(1, 1, 1, 1), "cm"))

# plot summarized data of correlation directions
cormat_summary2 <- cormat_summary %>% dplyr::group_by(direction) %>% summarize(count2 = sum(count))
summ2 <- ggplot(cormat_summary2, aes(x = direction, y = count2)) +
            geom_bar(stat = "identity", fill = c("gold3", "darkslategray")) +
            geom_text(aes(label = count2), size = 5, vjust = 0.5, hjust = 1, color = "white") +
            scale_fill_manual(values = c("gold3", "darkslategray")) +
            labs(y = "No. of significant correlations", x = "Direction") +
            coord_flip() +
            scale_y_continuous(expand = c(0, 0)) +
            theme_classic(base_size = 18)


ggsave(plot = correlation_heatmap(melted_cormat), file = "correlation_heatmap.jpg", width = 10, height = 10, dpi = 400)
ggsave(plot = correlation_heatmap(top10(melted_cormat), title = "Top 10", bsize = 28, text = FALSE, show_legend = FALSE, top10 = TRUE), file = "correlation_heatmap_top10.jpg", width = 10, height = 10, dpi = 400)
ggsave(plot = summ, file = "correlation_summary_direction.jpg", width = 10, height = 5.5, dpi = 400)
ggsave(plot = summ2, file = "correlation_summary_direction2.jpg", width = 6, height = 2, dpi = 400)

ggsave(plot = correlation_heatmap(melted_cormat), file = "correlation_heatmap.pdf", width = 10, height = 10)
ggsave(plot = correlation_heatmap(top10(melted_cormat), title = "Top 10", bsize = 28, text = FALSE, show_legend = FALSE, top10 = TRUE), file = "correlation_heatmap_top10.pdf", width = 10, height = 10)
ggsave(plot = summ, file = "correlation_summary_direction.pdf", width = 10, height = 5.5)
ggsave(plot = summ2, file = "correlation_summary_direction2.pdf", width = 6, height = 2)
