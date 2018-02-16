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

r <- raster(extent(fast_rast))
res(r) <- c(30,30)
projection(r) <- "+init=epsg:32733"
veg_raster <- fasterize(veg, raster=r, field='reclass', fun='first')
veg_df <- as.data.frame(veg_raster)

veg_df$prop_bare = 0
veg_df$prop_steppe = 0
veg_df$prop_grass = 0
veg_df$prop_shrub = 0
veg_df$prop_low = 0
veg_df$prop_high = 0

veg_df$prop_bare[veg_df$layer == 1] <- 1
veg_df$prop_steppe[veg_df$layer == 2] <- 1
veg_df$prop_grass[veg_df$layer == 3] <- 1
veg_df$prop_shrub[veg_df$layer == 4] <- 1
veg_df$prop_low[veg_df$layer == 5] <- 1
veg_df$prop_high[veg_df$layer == 6] <- 1

r <- raster(extent(fast_rast))
res(r) <- c(30,30)
projection(r) <- "+init=epsg:32733"
r[] <- 0
r <- setValues(r, veg_df$prop_bare)
writeRaster(r, "Prop_Bare.tif", format='GTiff')
r <- setValues(r, veg_df$prop_steppe)
writeRaster(r, "Prop_Steppe.tif", format='GTiff')
r <- setValues(r, veg_df$prop_grass)
writeRaster(r, "Prop_Grass.tif", format='GTiff')
r <- setValues(r, veg_df$prop_shrub)
writeRaster(r, "Prop_Shrub.tif", format='GTiff')
r <- setValues(r, veg_df$prop_low)
writeRaster(r, "Prop_Low.tif", format='GTiff')
r <- setValues(r, veg_df$prop_high)
writeRaster(r, "Prop_High.tif", format='GTiff')

