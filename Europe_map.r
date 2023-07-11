#Sampling map
#load packages
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(sp)
library(viridisLite)
library(RColorBrewer)
library(viridis)
library(kableExtra)
library(ggplot2)
library(dplyr)
library(cowplot)





 #show areas
 world <- ne_countries(scale = "medium", returnclass = "sf")

 #2019
 places2019 <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/R/Karten/Koordinaten/Soil/places2019.txt", 
 header = TRUE, sep = '\t', colClasses = c("character", "character", "character", "numeric", "numeric"), na.strings = c("","NA"))
 metadata2019 <- places2019
 coordinates(places2019) <- ~longitude + latitude
 proj4string(places2019) <- CRS("+proj=longlat +datum=WGS84")

 #2020
places2020 <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/R/Karten/Koordinaten/Soil/places2020.txt",  
header = TRUE, sep = '\t', colClasses = c("character", "character", "character", "numeric", "numeric"), na.strings = c("","NA"))
metadata2020 <- places2020
coordinates(places2020) <- ~longitude + latitude
proj4string(places2020) <- CRS("+proj=longlat +datum=WGS84")

#2021
places2021 <- read.table("C:/Users/Wilken Boie/Desktop/Pythium/R/Karten/Koordinaten/Soil/places2021.txt", 
header = TRUE, sep = '\t', colClasses = c("character", "character", "character", "numeric", "numeric"), na.strings = c("","NA"))
metadata2021 <- places2021
coordinates(places2021) <- ~longitude + latitude
proj4string(places2021) <- CRS("+proj=longlat +datum=WGS84")

#create border data frame to select specific countries
#set border to exclusivly display region of interest
aoi <- function(area){
    showcountry <- data.frame(country = c("Europe","France", "Germany", "Hungary", "Spain", "Belgium", "Italy", "Romania", "Netherlands", "Austria", "Switzerland"),
    x1 = c(-15, -5.855218, 5.091968, 15.059340, -10.366753, 2.018984, 5.522602, 19.819668, 2.482055, 8.523945, 5.442530),
                              x2 = c(32, 9.789313, 15.883674, 23.630425, 4.710256, 6.742236, 19.188612, 30.462262, 7.692492, 17.411251, 10.773058),
                              y1 = c(35, 41.789743, 46.935093, 45.212388, 35.359605, 49.258207, 35.844183, 42.610391, 50.485291, 45.856503, 45.475358),
                              y2 = c(60, 51.986570, 55.643060, 49.042165, 44.357896, 51.689979, 47.484770, 48.636408, 53.881754, 49.295000, 47.920979))
                    
        return(showcountry %>% filter(country == area))
        }

area_per_country <- data.frame(table(metadata2019$country))
area_per_country$Freq <- as.factor(area_per_country$Freq)
with_data <- left_join(world, area_per_country, by = c("name_long" = "Var1"))

#################
#   Map_2019    #
#################
europe_map_2019 <- ggplot(data = world) +
                    geom_sf() +
                    geom_sf(data = with_data, fill = "gray90")+
                    geom_point(data = metadata2019, aes(y = longitude, x = latitude), color = "gray30", fill = "black", size = 3, shape = 21, alpha = 0.8) +                 
                    xlab("Latitude") + ylab("Longitude") +
                    labs(fill = "No. of areas", title = "2019") +
                    annotation_scale(location = "bl", width_hint = 0.5) +
                    annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering) +
                    coord_sf(xlim = c(aoi("Europe")$x1, aoi("Europe")$x2), ylim = c(aoi("Europe")$y1, aoi("Europe")$y2)) + 
                    scale_fill_grey(na.value = "gray80", position = "top", start = 0.8, end = 0.2) + 
                    theme_classic(base_size = 14) +
                    theme(panel.background = element_rect(fill = "lightsteelblue1"), axis.text = element_text(size = 14))

europe_map_2019


#################
#   Map_2020    #
#################
europe_map_2020 <- ggplot(data = world) +
                    geom_sf() +
                    geom_sf(data = with_data, fill = "gray90")+
                    geom_point(data = metadata2020, aes(y = longitude, x = latitude), color = "gray30", fill = "black", size = 3, shape = 21, alpha = 0.8) +             
                    xlab("Latitude") + ylab("Longitude") +
                    labs(fill = "No. of areas", title = "2020") +
                    annotation_scale(location = "bl", width_hint = 0.5) +
                    annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering) +
                    coord_sf(xlim = c(aoi("Europe")$x1, aoi("Europe")$x2), ylim = c(aoi("Europe")$y1, aoi("Europe")$y2)) + 
                    scale_fill_grey(na.value = "gray80", position = "top", start = 0.8, end = 0.2) + 
                    theme_classic(base_size = 14) +
                    theme(panel.background = element_rect(fill = "lightsteelblue1"), axis.text = element_text(size = 14))

europe_map_2020


#################
#   Map_2021    #
#################
europe_map_2021 <- ggplot(data = world) +
                    geom_sf() +
                    geom_sf(data = with_data, fill = "gray90")+
                    geom_point(data = metadata2021, aes(y = longitude, x = latitude), color = "gray30", fill = "black", size = 3, shape = 21, alpha = 0.8) +                
                    xlab("Latitude") + ylab("Longitude") +
                    labs(fill = "No. of areas", title = "2021") +
                    annotation_scale(location = "bl", width_hint = 0.5) +
                    annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering) +
                    coord_sf(xlim = c(aoi("Europe")$x1, aoi("Europe")$x2), ylim = c(aoi("Europe")$y1, aoi("Europe")$y2)) + 
                    scale_fill_grey(na.value = "gray80", position = "top", start = 0.8, end = 0.2) + 
                    theme_classic(base_size = 14) +
                    theme(panel.background = element_rect(fill = "lightsteelblue1"), axis.text = element_text(size = 14))

europe_map_2021


#######################
#   Plot_all_years    #
#######################
plot_grid(europe_map_2019, europe_map_2020, europe_map_2021, ncol = 3)
ggsave(file="C:/Users/Wilken Boie/Desktop/Pythium/R/Karten/Abbildungen/Soil/map_all_years.jpeg", width = 15, height = 7)
