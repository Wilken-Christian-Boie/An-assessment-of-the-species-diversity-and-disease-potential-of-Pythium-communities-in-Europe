# get command line options
library(optparse)

option_list = list(
  make_option(c("-f", "--forward"), type = "character", default = NULL, help = "Path of forward read histogram"),
  make_option(c("-r", "--reverse"), type = "character", default = NULL, help = "Path of reverse read histogram"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Path to output folder"),
  make_option(c("--trimmf"), type = "integer", default = NULL, help = "Path to output folder"),
  make_option(c("--trimmr"), type = "integer", default = NULL, help = "Path to output folder")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$forward)){
  print_help(opt_parser)
  stop("At least one argument must be supplied\n", call.=FALSE)
}

library(ggplot2)
library(cowplot)

R1 <- read.table(opt$forward, sep = " ", head = FALSE) 
R2 <- read.table(opt$reverse, sep = " ", head = FALSE)

colnames(R1) <- c("bp", "count")
colnames(R2) <- c("bp", "count")

R1_plot <- ggplot(R1, aes(x = bp, y = log(count))) + 
            geom_bar(color="cornsilk3", fill = "cornsilk3", stat = 'identity') +
            geom_vline(xintercept = opt$trimmf, color = "firebrick") +
            labs(y = "log(count)", x = "Length (bp)") +
            scale_y_continuous(expand = c(0,0)) +
            scale_x_continuous(expand = c(0,0), limits = c(0,320)) +
            theme_classic(base_size = 16) +
            theme(legend.position = "none")

R2_plot <- ggplot(R2, aes(x = bp, y = log(count))) + 
            geom_bar(color="aquamarine4", fill = "aquamarine4", stat = 'identity') +
            geom_vline(xintercept = opt$trimmr, color = "firebrick") +
            labs(y = "log(count)", x = "Length (bp)") +
            scale_y_continuous(expand = c(0,0)) +
            scale_x_continuous(expand = c(0,0), limits = c(0,320)) +
            theme_classic(base_size = 16) +
            theme(legend.position = "none")

both_plots <- cowplot::plot_grid(R1_plot, R2_plot, ncol = 1, nrow = 2)
ggsave(opt$output, plot = both_plots, height = unit(4, "cm"), width = unit(5, "cm"))