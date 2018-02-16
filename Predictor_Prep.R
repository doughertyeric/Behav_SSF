library(raster)
library(lme4)
library(sf)
library(tidyverse)
library(adehabitatHR)
library(mapview)
library(fasterize)

#################################################

zebra09 <- read_csv("Zebra_Anthrax_2009_Cleaned.csv") %>% 
  dplyr::select(x,y,date,ID) %>% 
  dplyr::filter(!is.na(x)) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

zebra10 <- read_csv("Zebra_Anthrax_2010_Cleaned.csv") %>% 
  dplyr::select(x,y,date,ID) %>% 
  dplyr::filter(!is.na(x)) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

zebra.ext <- extent(c(extent(zebra10)@xmin - 5000,
                      extent(zebra10)@xmax + 5000,
                      extent(zebra09)@ymin - 5000,
                      extent(zebra10)@ymax + 5000))
  
#################################################

# Import vegetation polygons
veg <- st_read("ENP_Predictors/vegetation.shp", crs = "+init=epsg:4326")
veg <- st_transform(veg, "+init=epsg:32733")

for (i in 1:nrow(veg)) {
  if (as.character(veg$VEG_CLASS[i]) == "Steppe - L" || 
      as.character(veg$VEG_CLASS[i]) == "Steppe - M" || 
      as.character(veg$VEG_CLASS[i]) == "Steppe - VL") {
    veg$reclass[i] <- 'Steppe'
  } else if (as.character(veg$VEG_CLASS[i]) == "High tree savanna - H" || 
             as.character(veg$VEG_CLASS[i]) == "High tree savanna - M" || 
             as.character(veg$VEG_CLASS[i]) == "High tree savanna - L") {
    veg$reclass[i] <- 'High tree savanna'
  } else if (as.character(veg$VEG_CLASS[i]) == "Low tree savanna - H" || 
             as.character(veg$VEG_CLASS[i]) == "Low tree savanna - M" || 
             as.character(veg$VEG_CLASS[i]) == "Low tree savanna - L") {
    veg$reclass[i] <- 'Low tree savanna'
  } else if (as.character(veg$VEG_CLASS[i]) == "Shrub savanna - L" || 
             as.character(veg$VEG_CLASS[i]) == "Shrub savanna - M") {
    veg$reclass[i] <- 'Shrub savanna'
  } else if (as.character(veg$VEG_CLASS[i]) == "Grass savanna" || 
             as.character(veg$VEG_CLASS[i]) == "Grassland") {
    veg$reclass[i] <- 'Grassland'
  } else {
    veg$reclass[i] <- as.character(veg$VEG_CLASS[i])
  }
}

veg$reclass <- factor(veg$reclass, levels = c("Bare ground", "Steppe", "Grassland", 
                                              "Shrub savanna", "Low tree savanna", 'High tree savanna'))

r <- raster(extent(dist_roads))
res(r) <- c(30,30)
projection(r) <- "+init=epsg:32733"
veg_raster <- fasterize(veg, raster=r, field='reclass', fun='first')
writeRaster(veg_raster, "Vegetation_6Classes.tif", format='GTiff')

##################################################

file.list <- c("Etosha_East_20090322", "Etosha_East_20090423", "Etosha_East_20090509", "Etosha_East_20090525",
               "Etosha_East_20100205", "Etosha_East_20100410", "Etosha_East_20100512", "Etosha_East_20100528")

for (i in 1:length(file.list)) {
  RED <- raster(paste0("ENP_Predictors/", file.list[i], "/B4.TIF")) %>%
    projectRaster(veg_raster) %>% crop(zebra.ext) %>%
    resample(veg_raster, method='bilinear')
  NIR <- raster(paste0("ENP_Predictors/", file.list[i], "/B5.TIF")) %>%
    projectRaster(veg_raster) %>% crop(zebra.ext) %>%
    resample(veg_raster)
  NDVI <- (NIR - RED) / (NIR + RED)
  writeRaster(NDVI, paste0(file.list[i],"_NDVI"), format="GTiff", overwrite=TRUE)
  print(i)
}

Brightness.vals <- c(0.3037, 0.2793, 0.4343, 0.5585, 0.5082, 0.1863)
Greenness.vals <-	c(-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800)
Wetness.vals <-	c(0.1509, 0.1793, 0.3299, 0.3406, -0.7112, -0.4572)

for (i in 1:length(file.list)) {
  CH1 <- raster(paste0("ENP_Predictors/", file.list[i], "/B1.TIF")) %>%
    projectRaster(veg_raster) %>% crop(zebra.ext) %>%
    resample(veg_raster, method='bilinear')
  CH2 <- raster(paste0("ENP_Predictors/", file.list[i], "/B2.TIF")) %>%
    projectRaster(veg_raster) %>% crop(zebra.ext) %>%
    resample(veg_raster, method='bilinear')
  CH3 <- raster(paste0("ENP_Predictors/", file.list[i], "/B3.TIF")) %>%
    projectRaster(veg_raster) %>% crop(zebra.ext) %>%
    resample(veg_raster, method='bilinear')
  CH4 <- raster(paste0("ENP_Predictors/", file.list[i], "/B4.TIF")) %>%
    projectRaster(veg_raster) %>% crop(zebra.ext) %>%
    resample(veg_raster, method='bilinear')
  CH5 <- raster(paste0("ENP_Predictors/", file.list[i], "/B5.TIF")) %>%
    projectRaster(veg_raster) %>% crop(zebra.ext) %>%
    resample(veg_raster, method='bilinear')
  CH7 <- raster(paste0("ENP_Predictors/", file.list[i], "/B7.TIF")) %>%
    projectRaster(veg_raster) %>% crop(zebra.ext) %>%
    resample(veg_raster, method='bilinear')
  #Brightness <- (Brightness.vals[1] * CH1) + (Brightness.vals[2] * CH2) + (Brightness.vals[3] * CH3) + 
  #  (Brightness.vals[4] * CH4) + (Brightness.vals[5] * CH5) + (Brightness.vals[6] * CH7)
  #writeRaster(Brightness, paste0(file.list[i],"_Brightness"), format="GTiff", overwrite=TRUE)
  Greenness <- (Greenness.vals[1] * CH1) + (Greenness.vals[2] * CH2) + (Greenness.vals[3] * CH3) + 
    (Greenness.vals[4] * CH4) + (Greenness.vals[5] * CH5) + (Greenness.vals[6] * CH7)
  writeRaster(Greenness, paste0(file.list[i],"_Greenness"), format="GTiff", overwrite=TRUE)
  Wetness <- (Wetness.vals[1] * CH1) + (Wetness.vals[2] * CH2) + (Wetness.vals[3] * CH3) + 
    (Wetness.vals[4] * CH4) + (Wetness.vals[5] * CH5) + (Wetness.vals[6] * CH7)
  writeRaster(Wetness, paste0(file.list[i],"_Wetness"), format="GTiff", overwrite=TRUE)
  print(i)
}

################################################

# Import ENP fenceline to use as extent
ENP <- st_read("ENP_Predictors/enp fence poly.shp") %>% 
  st_union
ENP <- st_transform(ENP, "+init=epsg:32733")

ext.poly <- as(zebra.ext, 'SpatialPolygons')
#enp_grid <- st_make_grid(ext.poly, cellsize=c(30,30), what= "centers") #n = c(1266, 476)
#st_crs(enp_grid) = 32733
# convert to raster
#dist_raster <- raster(as(enp_grid, "Spatial"))

dist_raster <- raster(ext.poly, res=c(30,30))
projection(dist_raster) <- "+init=epsg:32733"
dist_raster[] <- 1
grid <- data.frame(coordinates(dist_raster)[!is.na(values(dist_raster)),]) %>%
  st_as_sf(., coords=1:2, crs="+init=epsg:32733")

#Import roads file, transform, and crop
roads <- st_read("ENP_Predictors/enp roads.shp", crs = "+init=epsg:4326") %>%
  filter(TYPE %in% c("Gravel", "Tar"))
roads <- st_transform(roads, "+init=epsg:32733")
#roads <- st_intersection(roads, st_set_crs(st_as_sf(ext.poly), st_crs(roads)))
#ggplot(roads) + geom_sf()

#Calculate distances from each raster cell center to nearest primary road
dist_road_grid <- st_distance(roads, grid)
dist_raster[] = apply(dist_road_grid, 2, min)
writeRaster(dist_raster, "Dist_PrimaryRoads.tif", format="GTiff", overwrite=TRUE)

#Import water file, transform, and crop
water <- st_read("ENP_Predictors/functional water.shp") %>% 
  st_transform("+init=epsg:32733")
#water <- st_intersection(water, st_set_crs(st_as_sf(ext.poly), st_crs(water)))
#ggplot(water) + geom_sf()

#Calculate distances from each raster cell center to nearest functional water source
dist_water_grid <- st_distance(water, grid)
dist_raster[] = apply(dist_water_grid, 2, min)
writeRaster(dist_raster, "Dist_Water.tif", format="GTiff")

# Import temporary files, flip, and save a new flipped version
dist_road <- raster("Dist_PrimaryRoads.tif") %>% flip(2)
writeRaster(dist_road, "Dist_PrimaryRoads_Flip.tif", format="GTiff", overwrite=TRUE)
dist_water <- raster("Dist_Water.tif") %>% flip(2)
writeRaster(dist_water, "Dist_Water_Flip.tif", format="GTiff", overwrite=TRUE)

###########################################

list.09 <- c("Etosha_East_20090322", "Etosha_East_20090423", "Etosha_East_20090509", "Etosha_East_20090525")
list.10 <- c("Etosha_East_20100205", "Etosha_East_20100410", "Etosha_East_20100512", "Etosha_East_20100528")

NDVI.rast1 <- raster(paste0('ENP_Predictors/', list.09[1], '_NDVI.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
NDVI.rast2 <- raster(paste0('ENP_Predictors/', list.09[2], '_NDVI.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
NDVI.rast3 <- raster(paste0('ENP_Predictors/', list.09[3], '_NDVI.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
NDVI.rast4 <- raster(paste0('ENP_Predictors/', list.09[4], '_NDVI.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
NDVI.09 <- mean(NDVI.rast1, NDVI.rast2, NDVI.rast3, NDVI.rast4)
NDVI.09.max <- max(NDVI.rast1, NDVI.rast2, NDVI.rast3, NDVI.rast4)
NDVI.09.min <- min(NDVI.rast1, NDVI.rast2, NDVI.rast3, NDVI.rast4)
NDVI.09.range <- range(NDVI.rast1, NDVI.rast2, NDVI.rast3, NDVI.rast4)
writeRaster(NDVI.09, "Mean_NDVI_2009.tif", format='GTiff')
writeRaster(NDVI.09.max, "Max_NDVI_2009.tif", format='GTiff')
writeRaster(NDVI.09.min, "Min_NDVI_2009.tif", format='GTiff')
writeRaster(NDVI.09.range, "Range_NDVI_2009.tif", format='GTiff')

NDVI.rast1 <- raster(paste0('ENP_Predictors/', list.10[1], '_NDVI.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
NDVI.rast2 <- raster(paste0('ENP_Predictors/', list.10[2], '_NDVI.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
NDVI.rast3 <- raster(paste0('ENP_Predictors/', list.10[3], '_NDVI.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
NDVI.rast4 <- raster(paste0('ENP_Predictors/', list.10[4], '_NDVI.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
NDVI.10 <- mean(NDVI.rast1, NDVI.rast2, NDVI.rast3, NDVI.rast4)
NDVI.10.max <- max(NDVI.rast1, NDVI.rast2, NDVI.rast3, NDVI.rast4)
NDVI.10.min <- min(NDVI.rast1, NDVI.rast2, NDVI.rast3, NDVI.rast4)
NDVI.10.range <- range(NDVI.rast1, NDVI.rast2, NDVI.rast3, NDVI.rast4)
writeRaster(NDVI.10, "Mean_NDVI_2010.tif", format='GTiff')
writeRaster(NDVI.10.max, "Max_NDVI_2010.tif", format='GTiff')
writeRaster(NDVI.10.min, "Min_NDVI_2010.tif", format='GTiff')
writeRaster(NDVI.10.range, "Range_NDVI_2010.tif", format='GTiff')

Greenness.rast1 <- raster(paste0('ENP_Predictors/', list.09[1], '_Greenness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Greenness.rast2 <- raster(paste0('ENP_Predictors/', list.09[2], '_Greenness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Greenness.rast3 <- raster(paste0('ENP_Predictors/', list.09[3], '_Greenness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Greenness.rast4 <- raster(paste0('ENP_Predictors/', list.09[4], '_Greenness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Greenness.09 <- mean(Greenness.rast1, Greenness.rast2, Greenness.rast3, Greenness.rast4)
writeRaster(Greenness.09, "Mean_Greenness_2009.tif", format='GTiff')

Greenness.rast1 <- raster(paste0('ENP_Predictors/', list.10[1], '_Greenness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Greenness.rast2 <- raster(paste0('ENP_Predictors/', list.10[2], '_Greenness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Greenness.rast3 <- raster(paste0('ENP_Predictors/', list.10[3], '_Greenness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Greenness.rast4 <- raster(paste0('ENP_Predictors/', list.10[4], '_Greenness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Greenness.10 <- mean(Greenness.rast1, Greenness.rast2, Greenness.rast3, Greenness.rast4)
writeRaster(Greenness.10, "Mean_Greenness_2010.tif", format='GTiff')

Wetness.rast1 <- raster(paste0('ENP_Predictors/', list.09[1], '_Wetness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Wetness.rast2 <- raster(paste0('ENP_Predictors/', list.09[2], '_Wetness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Wetness.rast3 <- raster(paste0('ENP_Predictors/', list.09[3], '_Wetness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Wetness.rast4 <- raster(paste0('ENP_Predictors/', list.09[4], '_Wetness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Wetness.09 <- mean(Wetness.rast1, Wetness.rast2, Wetness.rast3, Wetness.rast4)
writeRaster(Wetness.09, "Mean_Wetness_2009.tif", format='GTiff')

Wetness.rast1 <- raster(paste0('ENP_Predictors/', list.10[1], '_Wetness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Wetness.rast2 <- raster(paste0('ENP_Predictors/', list.10[2], '_Wetness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Wetness.rast3 <- raster(paste0('ENP_Predictors/', list.10[3], '_Wetness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Wetness.rast4 <- raster(paste0('ENP_Predictors/', list.10[4], '_Wetness.tif')) %>%
  projectRaster(veg_raster) %>% crop(zebra.ext)
Wetness.10 <- mean(Wetness.rast1, Wetness.rast2, Wetness.rast3, Wetness.rast4)
writeRaster(Wetness.10, "Mean_Wetness_2010.tif", format='GTiff')

################################################

dist_roads <- raster('ENP_Predictors/Dist_PrimaryRoads.tif')
dist_water <- raster('ENP_Predictors/Dist_Water.tif')
veg_raster <- raster('ENP_Predictors/Vegetation_Crop.tif')
NDVI_09 <- raster('ENP_Predictors/Mean_NDVI_2009.tif')
NDVI_10 <- raster('ENP_Predictors/Mean_NDVI_2010.tif')
Green_09 <- raster('ENP_Predictors/Mean_Greenness_2009.tif')
Green_10 <- raster('ENP_Predictors/Mean_Greenness_2010.tif')
Wet_09 <- raster('ENP_Predictors/Mean_Wetness_2009.tif')
Wet_10 <- raster('ENP_Predictors/Mean_Wetness_2010.tif')




