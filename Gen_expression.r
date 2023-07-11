#############################
#           qRT-PCR         #
#############################

#load packages
library(ggplot2)
library(xlsx)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(plotly)
library(patchwork)
library(dplyr)
library(plyr)
library(scales)




##############################
#     With y-axis-break      #
##############################

metadata <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Gene_expression/qRT_PCR_data.txt", header = TRUE, sep = "\t")

dat <- ddply(metadata, c("Target", "Species"), summarize, Value_mean = mean(Value),SD = sd(Value))


#Function to transform data to y positions
trans <- function(x){pmin(x,4) + 0.03*pmax(x-3,0)}
yticks <- c(0, 1, 2, 3, 20, 40, 60, 80, 100)

#Transform the data onto the display scale
dat$mean_t <- trans(dat$Value_mean)
dat$sd_up_t <- trans(dat$Value_mean + dat$SD)
dat$sd_low_t <- pmax(trans(dat$Value_mean - dat$SD),0.1)


gen_expres <- ggplot2::ggplot(data = dat, aes(x = Target, y = mean_t, group = Species, fill = Species)) +
                ggplot2::geom_bar(stat = "identity", color = "black", position = "dodge") +
                ggplot2::geom_errorbar(aes(ymin = sd_low_t, ymax = sd_up_t), position = position_dodge(width = 0.9), color = "black", size = 0.5, width = 0.3) +
                ggplot2::scale_fill_manual(values = c("NC" = "#ffff99", "P. attrantheridium" = "#1f78b4", "P. ultimum" = "#6a3d9a")) +
                ggplot2::annotate("text", x = 3.3, y = dat$mean_t[9] + 0.2, size = 4, label = "***", fontface = "bold", color = "black") +
                ggplot2::annotate("text", x = 4, y = dat$mean_t[11] + 0.2, size = 4, label = "***", fontface = "bold", color = "black") +
                ggplot2::annotate("text", x = 4.3, y = dat$mean_t[12] + 0.3, size = 4, label = "***", fontface = "bold", color = "black") +
                ggplot2::theme_bw() +
                ggplot2::geom_rect(aes(xmin = 0, xmax = 5, ymin = 3.2, ymax = 4), fill = "white") +
                ggplot2::scale_y_continuous(limits  = c(0, NA), breaks = trans(yticks), labels = yticks) +
                ggplot2::labs(y = "relative gene expression") +
                ggplot2::labs(x = "") +
                ggplot2::theme(legend.position = "bottom") +
                ggplot2::labs(fill = "Species") +
                ggplot2::theme(text = element_text(size = 12, color = "black")) +
                ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")) +
                                  theme(axis.text.y = element_text(size = 12, color = "black"),
                                  legend.text = element_text(size = 12, face = "italic"),
                                  axis.text.x = element_text(size = 12, face = "italic", color = "black"),
                                  axis.title.y = element_text(size = 12, color = "black"))

gen_expres




png(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Pythium/Output/Gene_expression/gene_expression.png", width = 1000, height = 750, res = 150)
gen_expres
dev.off()