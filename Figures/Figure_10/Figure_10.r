library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(scales)
library(ggbreak)

set.seed(42)

# import and prepare data
qRT_PCR_data <- read.table("../Data/qRT_PCR_data.txt", header = TRUE, sep = "\t")
qRT_PCR_data$Target <- gsub("PDF1,2", "PDF1.2", qRT_PCR_data$Target)

dat <- plyr::ddply(qRT_PCR_data,
                   c("Target", "Species"),
                   summarize,
                   Value_mean = mean(Expression),
                   SD = sd(Expression))

# graphical properties of the plot
color_codes <- c("NC" = "#FFFF99",
                 "G. attrantheridium" = "#1F78B4",
                 "G. ultimum var. ultimum" = "#FDBF6F")

# create gene expression plot
dat$Species <- factor(dat$Species, levels = c("NC", "G. attrantheridium", "G. ultimum var. ultimum"))
qRT_PCR_data$Species <- factor(qRT_PCR_data$Species, levels = c("NC", "G. attrantheridium", "G. ultimum var. ultimum"))
gene_expression <- ggplot(data = dat, aes(x = Target, y = Value_mean, group = Species, fill = Species)) +
                    geom_bar(stat = "identity", color = "black", position = "dodge") +
					geom_point(data = qRT_PCR_data, aes(x = Target, y = Expression, group = Species, fill = Species), position = position_dodge(width = 0.9)) + 
                    geom_errorbar(aes(ymin = dat$Value_mean - dat$SD, ymax = dat$Value_mean + dat$SD), position = position_dodge(width = 0.9), color = "black", linewidth = 0.5, width = 0.3) +
                    scale_fill_manual(values = color_codes) +
                    scale_y_continuous(limits  = c(0, 120), expand = c(0, 0)) +
					scale_y_break(breaks = c(2.5, 20), scales = 0.4) + 
                    annotate("text", x = 3.3, y = 0.5, size = 4, label = "***", fontface = "bold", color = "black") +
                    annotate("text", x = 4.0, y = 25, size = 4, label = "***", fontface = "bold", color = "black") +
                    annotate("text", x = 4.3, y = 110, size = 4, label = "***", fontface = "bold", color = "black") +
                    labs(x = "", y = "Relative gene expression", fill = "Species") +
                    theme_bw() +
                    theme(legend.position = "bottom",
                          text = element_text(size = 12, color = "black"),
                          axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12, face = "italic", color = "black"),
                          axis.text.y = element_text(size = 12, color = "black"),
                          legend.text = element_text(size = 12, face = "italic"),
                          axis.title.y = element_text(size = 12, color = "black"))

jpeg(file = "Figure_11.jpg", width = 5000, height = 4000, res = 600)
gene_expression
dev.off()  

pdf(file = "Figure_11.pdf", width = 10, height = 8)
gene_expression
dev.off()
