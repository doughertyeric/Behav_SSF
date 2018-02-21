anth09 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Anthrax_GARP/Rerun_1km/B.anthracis_2009_1km.tif")
anth10 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Anthrax_GARP/Rerun_1km/B.anthracis_2010_1km.tif")

#setwd("~/Box Sync/Dissertation/Behavioral_SSF/Figures")
#pdf("B.anthracis_2009_Prediction_1km.pdf", height=6, width=8)
#plot(anth09)
#dev.off()
#pdf("B.anthracis_2010_Prediction_1km.pdf", height=6, width=8)
#plot(anth10)
#dev.off()

# Create a binomial plot of the probability of anthrax presence given a relatively liberal threshold
anth09_df <- as.data.frame(anth09)
for (i in 1:nrow(anth09_df)) {
  if (!is.na(anth09_df[i,1])) {
    if (anth09_df[i,1] < 100) {
      anth09_df[i,1] <- 0 
    } else {
      anth09_df[i,1] <- 1
    }
  }
}
anth09_liberal <- setValues(anth09, anth09_df[,1])
plot(anth09_liberal)


# Create a similar binomial plot with a moderate estimate of the area with anthrax present
anth09_df <- as.data.frame(anth09)
for (i in 1:nrow(anth09_df)) {
  if (!is.na(anth09_df[i,1])) {
    if (anth09_df[i,1] < 500) {
      anth09_df[i,1] <- 0 
    } else {
      anth09_df[i,1] <- 1
    }
  }
}
anth09_moderate <- setValues(anth09, anth09_df[,1])
plot(anth09_moderate)

# Create a final binomial plot with the most conservative estimate
anth09_df <- as.data.frame(anth09)
for (i in 1:nrow(anth09_df)) {
  if (!is.na(anth09_df[i,1])) {
    if (anth09_df[i,1] < 750) {
      anth09_df[i,1] <- 0 
    } else {
      anth09_df[i,1] <- 1
    }
  }
}
anth09_conserve <- setValues(anth09, anth09_df[,1])
plot(anth09_conserve)

#######

# Create a binomial plot of the probability of anthrax presence given a relatively liberal threshold
anth10_df <- as.data.frame(anth10)
for (i in 1:nrow(anth10_df)) {
  if (!is.na(anth10_df[i,1])) {
    if (anth10_df[i,1] < 100) {
      anth10_df[i,1] <- 0 
    } else {
      anth10_df[i,1] <- 1
    }
  }
}
anth10_liberal <- setValues(anth10, anth10_df[,1])
plot(anth10_liberal)

# Create a similar binomial plot with a moderate estimate of the area with anthrax present
anth10_df <- as.data.frame(anth10)
for (i in 1:nrow(anth10_df)) {
  if (!is.na(anth10_df[i,1])) {
    if (anth10_df[i,1] < 500) {
      anth10_df[i,1] <- 0 
    } else {
      anth10_df[i,1] <- 1
    }
  }
}
anth10_moderate <- setValues(anth10, anth10_df[,1])
plot(anth10_moderate)

# Create a final binomial plot with the most conservative estimate
anth10_df <- as.data.frame(anth10)
for (i in 1:nrow(anth10_df)) {
  if (!is.na(anth10_df[i,1])) {
    if (anth10_df[i,1] < 750) {
      anth10_df[i,1] <- 0 
    } else {
      anth10_df[i,1] <- 1
    }
  }
}
anth10_conserve <- setValues(anth10, anth10_df[,1])
plot(anth10_conserve)

###################################################

all09 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Rerun_Results/Zebra_2009_All_Risk.tif")
v_all09 <- velox(all09)
forage09 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Rerun_Results/Zebra_2009_Forage_Risk.tif")
v_forage09 <- velox(forage09)
all10 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Rerun_Results/Zebra_2010_All_Risk.tif")
v_all10 <- velox(all10)
forage10 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Rerun_Results/Zebra_2010_Forage_Risk.tif")
v_forage10 <- velox(forage10)

anth09_lib_poly <- SpatialPolygons(list(rasterToPolygons(anth09_liberal, dissolve=TRUE)@polygons[[2]]))
anth09_lib_sf <- st_as_sf(anth09_lib_poly)
anth09_mod_poly <- SpatialPolygons(list(rasterToPolygons(anth09_moderate, dissolve=TRUE)@polygons[[2]]))
anth09_mod_sf <- st_as_sf(anth09_mod_poly)
anth09_cons_poly <- SpatialPolygons(list(rasterToPolygons(anth09_conserve, dissolve=TRUE)@polygons[[2]]))
anth09_cons_sf <- st_as_sf(anth09_cons_poly)

v_all09$extract(anth09_lib_sf, fun=sum) #689017.8
v_all09$extract(anth09_mod_sf, fun=sum) #248667.4
v_all09$extract(anth09_cons_sf, fun=sum) #154536.1

v_forage09$extract(anth09_lib_sf, fun=sum) #783595.3
v_forage09$extract(anth09_mod_sf, fun=sum) #264261.1
v_forage09$extract(anth09_cons_sf, fun=sum) #156968.1

anth10_lib_poly <- SpatialPolygons(list(rasterToPolygons(anth10_liberal, dissolve=TRUE)@polygons[[2]]))
anth10_lib_sf <- st_as_sf(anth10_lib_poly)
anth10_mod_poly <- SpatialPolygons(list(rasterToPolygons(anth10_moderate, dissolve=TRUE)@polygons[[2]]))
anth10_mod_sf <- st_as_sf(anth10_mod_poly)
anth10_cons_poly <- SpatialPolygons(list(rasterToPolygons(anth10_conserve, dissolve=TRUE)@polygons[[2]]))
anth10_cons_sf <- st_as_sf(anth10_cons_poly)

v_all10$extract(anth10_lib_sf, fun=sum) #591296.6
v_all10$extract(anth10_mod_sf, fun=sum) #245518.3
v_all10$extract(anth10_cons_sf, fun=sum) #153657.6

v_forage10$extract(anth10_lib_sf, fun=sum) #727306
v_forage10$extract(anth10_mod_sf, fun=sum) #262596.6
v_forage10$extract(anth10_cons_sf, fun=sum) #156903.1

#####

v_all09$extract(anth09_lib_sf, fun=mean) #0.357
v_all09$extract(anth09_mod_sf, fun=mean) #0.418
v_all09$extract(anth09_cons_sf, fun=mean) #0.465

v_forage09$extract(anth09_lib_sf, fun=mean) #0.406
v_forage09$extract(anth09_mod_sf, fun=mean) #0.445
v_forage09$extract(anth09_cons_sf, fun=mean) #0.472

v_all09$extract(anth09_lib_sf_out, fun=mean) #0.273
v_all09$extract(anth09_mod_sf_out, fun=mean) #0.283
v_all09$extract(anth09_cons_sf_out, fun=mean) #0.286

v_forage09$extract(anth09_lib_sf_out, fun=mean) #0.345
v_forage09$extract(anth09_mod_sf_out, fun=mean) #0.354
v_forage09$extract(anth09_cons_sf_out, fun=mean) #0.356


##############################################################

anth09_lib_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_liberal, dissolve=TRUE)@polygons[[1]]))
anth09_lib_sf_out <- st_as_sf(anth09_lib_poly_out)
anth09_mod_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_moderate, dissolve=TRUE)@polygons[[1]]))
anth09_mod_sf_out <- st_as_sf(anth09_mod_poly_out)
anth09_cons_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_conserve, dissolve=TRUE)@polygons[[1]]))
anth09_cons_sf_out <- st_as_sf(anth09_cons_poly_out)

v_all09$extract(anth09_lib_sf_out, fun=sum) #1632125
v_all09$extract(anth09_mod_sf_out, fun=sum) #2072475
v_all09$extract(anth09_cons_sf_out, fun=sum) #2166607

v_all09$extract(anth09_lib_sf, fun=sum)/v_all09$extract(anth09_lib_sf_out, fun=sum) #0.42216
v_all09$extract(anth09_mod_sf, fun=sum)/v_all09$extract(anth09_mod_sf_out, fun=sum) #0.1199857
v_all09$extract(anth09_cons_sf, fun=sum)/v_all09$extract(anth09_cons_sf_out, fun=sum) #0.07132631

v_forage09$extract(anth09_lib_sf_out, fun=sum) #2063128
v_forage09$extract(anth09_mod_sf_out, fun=sum) #2582462
v_forage09$extract(anth09_cons_sf_out, fun=sum) #2689755

v_forage09$extract(anth09_lib_sf, fun=sum)/v_forage09$extract(anth09_lib_sf_out, fun=sum) #0.3798094
v_forage09$extract(anth09_mod_sf, fun=sum)/v_forage09$extract(anth09_mod_sf_out, fun=sum) #0.1023291
v_forage09$extract(anth09_cons_sf, fun=sum)/v_forage09$extract(anth09_cons_sf_out, fun=sum) #0.05835779

anth10_lib_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_liberal, dissolve=TRUE)@polygons[[1]]))
anth10_lib_sf_out <- st_as_sf(anth10_lib_poly_out)
anth10_mod_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_moderate, dissolve=TRUE)@polygons[[1]]))
anth10_mod_sf_out <- st_as_sf(anth10_mod_poly_out)
anth10_cons_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_conserve, dissolve=TRUE)@polygons[[1]]))
anth10_cons_sf_out <- st_as_sf(anth10_cons_poly_out)

v_all10$extract(anth10_lib_sf_out, fun=sum) #604706.1
v_all10$extract(anth10_mod_sf_out, fun=sum) #950484.4
v_all10$extract(anth10_cons_sf_out, fun=sum) #1042345

v_all10$extract(anth10_lib_sf, fun=sum)/v_all10$extract(anth10_lib_sf_out, fun=sum) #0.9778246
v_all10$extract(anth10_mod_sf, fun=sum)/v_all10$extract(anth10_mod_sf_out, fun=sum) #0.2583086
v_all10$extract(anth10_cons_sf, fun=sum)/v_all10$extract(anth10_cons_sf_out, fun=sum) #0.1474153

v_forage10$extract(anth10_lib_sf_out, fun=sum) #1292533
v_forage10$extract(anth10_mod_sf_out, fun=sum) #1757242
v_forage10$extract(anth10_cons_sf_out, fun=sum) #1862936

v_forage10$extract(anth10_lib_sf, fun=sum)/v_forage10$extract(anth10_lib_sf_out, fun=sum) #0.5626983
v_forage10$extract(anth10_mod_sf, fun=sum)/v_forage10$extract(anth10_mod_sf_out, fun=sum) #0.1494368
v_forage10$extract(anth10_cons_sf, fun=sum)/v_forage10$extract(anth10_cons_sf_out, fun=sum) #0.08422357

#######

v_all10$extract(anth10_lib_sf, fun=mean) #0.283
v_all10$extract(anth10_mod_sf, fun=mean) #0.406
v_all10$extract(anth10_cons_sf, fun=mean) #0.456

v_forage10$extract(anth10_lib_sf, fun=mean) #0.348
v_forage10$extract(anth10_mod_sf, fun=mean) #0.434
v_forage10$extract(anth10_cons_sf, fun=mean) #0.466

v_all10$extract(anth10_lib_sf_out, fun=mean) #0.104
v_all10$extract(anth10_mod_sf_out, fun=mean) #0.130
v_all10$extract(anth10_cons_sf_out, fun=mean) #0.138

v_forage10$extract(anth10_lib_sf_out, fun=mean) #0.223
v_forage10$extract(anth10_mod_sf_out, fun=mean) #0.241
v_forage10$extract(anth10_cons_sf_out, fun=mean) #0.246

