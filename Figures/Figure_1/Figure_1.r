library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(sp)
library(ggplot2)
library(patchwork)

# import and prepare data
world <- ne_countries(scale = "medium", returnclass = "sf")

sites <- function(filepath) {
	return(read.table(filepath, header = TRUE, sep = '\t', colClasses = c("character", "character", "character", "numeric", "numeric"), na.strings = c("", "NA")))
}

# load geographic coordinates
sites_2019 <- sites("../Data/Coordinates_2019.txt")
sites_2020 <- sites("../Data/Coordinates_2020.txt")
sites_2021 <- sites("../Data/Coordinates_2021.txt")

map_border <- data.frame(x1 = -15, x2 = 32, y1 = 35, y2 = 60)

show_map <- function(sampling_sites, year) {
	ggplot(data = world) +
		geom_sf() +
		geom_point(data = sampling_sites, aes(y = longitude, x = latitude), color = "gray30", fill = "black", size = 3, shape = 21, alpha = 0.8) +                 
		xlab("Latitude") + ylab("Longitude") +
		labs(fill = "No. of areas", title = year) +
		annotation_scale(location = "bl", width_hint = 0.5) +
		annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering) +
		coord_sf(xlim = c(map_border$x1, map_border$x2), ylim = c(map_border$y1, map_border$y2)) + 
		scale_fill_grey(na.value = "gray80", position = "top", start = 0.8, end = 0.2) + 
		theme_classic(base_size = 14) +
		theme(panel.background = element_rect(fill = "lightsteelblue1"), axis.text = element_text(size = 14))
}

# plot and export figure
jpeg("../Figure_1/Figure_1.jpg", width = 6000, height = 3000, res = 600)
show_map(sites_2019, "2019") | show_map(sites_2020, "2020") | show_map(sites_2021, "2021")
dev.off()

pdf("../Figure_1/Figure_1.pdf", width = 12, height = 6)
show_map(sites_2019, "2019") | show_map(sites_2020, "2020") | show_map(sites_2021, "2021")
dev.off()
