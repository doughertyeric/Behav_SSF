anth09 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Anthrax_GARP/Rerun_1km_New/B.anthracis_2009_1km.tif")
anth10 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Anthrax_GARP/Rerun_1km_New/B.anthracis_2010_1km.tif")

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
#pdf('Anthrax_2009_Liberal.pdf', height=6, width=8)
#plot(anth09_liberal)
#dev.off()


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
#pdf('Anthrax_2009_Moderate.pdf', height=6, width=8)
#plot(anth09_moderate)
#dev.off()

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
#pdf('Anthrax_2009_Conservative.pdf', height=6, width=8)
#plot(anth09_conserve)
#dev.off()

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
#pdf('Anthrax_2010_Liberal.pdf', height=6, width=8)
#plot(anth10_liberal)
#dev.off()

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
#pdf('Anthrax_2010_Moderate.pdf', height=6, width=8)
#plot(anth10_moderate)
#dev.off()

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
#pdf('Anthrax_2010_Conservative.pdf', height=6, width=8)
#plot(anth10_conserve)
#dev.off()

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

v_all09$extract(anth09_lib_sf, fun=sum) #918501.3
v_all09$extract(anth09_mod_sf, fun=sum) #436484.1
v_all09$extract(anth09_cons_sf, fun=sum) #273227.8

v_forage09$extract(anth09_lib_sf, fun=sum) #1050767
v_forage09$extract(anth09_mod_sf, fun=sum) #458987.6
v_forage09$extract(anth09_cons_sf, fun=sum) #278090

anth10_lib_poly <- SpatialPolygons(list(rasterToPolygons(anth10_liberal, dissolve=TRUE)@polygons[[2]]))
anth10_lib_sf <- st_as_sf(anth10_lib_poly)
anth10_mod_poly <- SpatialPolygons(list(rasterToPolygons(anth10_moderate, dissolve=TRUE)@polygons[[2]]))
anth10_mod_sf <- st_as_sf(anth10_mod_poly)
anth10_cons_poly <- SpatialPolygons(list(rasterToPolygons(anth10_conserve, dissolve=TRUE)@polygons[[2]]))
anth10_cons_sf <- st_as_sf(anth10_cons_poly)

v_all10$extract(anth10_lib_sf, fun=sum) #642090.7
v_all10$extract(anth10_mod_sf, fun=sum) #286855.1
v_all10$extract(anth10_cons_sf, fun=sum) #202190.3

v_forage10$extract(anth10_lib_sf, fun=sum) #762175.8
v_forage10$extract(anth10_mod_sf, fun=sum) #289575.4
v_forage10$extract(anth10_cons_sf, fun=sum) #194998.5

##############################################################

anth09_lib_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_liberal, dissolve=TRUE)@polygons[[1]]))
anth09_lib_sf_out <- st_as_sf(anth09_lib_poly_out)
anth09_mod_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_moderate, dissolve=TRUE)@polygons[[1]]))
anth09_mod_sf_out <- st_as_sf(anth09_mod_poly_out)
anth09_cons_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_conserve, dissolve=TRUE)@polygons[[1]]))
anth09_cons_sf_out <- st_as_sf(anth09_cons_poly_out)

v_all09$extract(anth09_lib_sf_out, fun=sum) #1402641
v_all09$extract(anth09_mod_sf_out, fun=sum) #1884659
v_all09$extract(anth09_cons_sf_out, fun=sum) #2047915

v_all09$extract(anth09_lib_sf, fun=sum)/v_all09$extract(anth09_lib_sf_out, fun=sum) #0.6548368
v_all09$extract(anth09_mod_sf, fun=sum)/v_all09$extract(anth09_mod_sf_out, fun=sum) #0.2315985
v_all09$extract(anth09_cons_sf, fun=sum)/v_all09$extract(anth09_cons_sf_out, fun=sum) #0.1334175

v_forage09$extract(anth09_lib_sf_out, fun=sum) #1795956
v_forage09$extract(anth09_mod_sf_out, fun=sum) #2387735
v_forage09$extract(anth09_cons_sf_out, fun=sum) #2568633

v_forage09$extract(anth09_lib_sf, fun=sum)/v_forage09$extract(anth09_lib_sf_out, fun=sum) #0.585074
v_forage09$extract(anth09_mod_sf, fun=sum)/v_forage09$extract(anth09_mod_sf_out, fun=sum) #0.1922272
v_forage09$extract(anth09_cons_sf, fun=sum)/v_forage09$extract(anth09_cons_sf_out, fun=sum) #0.1082638

anth10_lib_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_liberal, dissolve=TRUE)@polygons[[1]]))
anth10_lib_sf_out <- st_as_sf(anth10_lib_poly_out)
anth10_mod_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_moderate, dissolve=TRUE)@polygons[[1]]))
anth10_mod_sf_out <- st_as_sf(anth10_mod_poly_out)
anth10_cons_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_conserve, dissolve=TRUE)@polygons[[1]]))
anth10_cons_sf_out <- st_as_sf(anth10_cons_poly_out)

v_all10$extract(anth10_lib_sf_out, fun=sum) #553912
v_all10$extract(anth10_mod_sf_out, fun=sum) #909147.6
v_all10$extract(anth10_cons_sf_out, fun=sum) #993812.4

v_all10$extract(anth10_lib_sf, fun=sum)/v_all10$extract(anth10_lib_sf_out, fun=sum) #1.159193
v_all10$extract(anth10_mod_sf, fun=sum)/v_all10$extract(anth10_mod_sf_out, fun=sum) #0.315521
v_all10$extract(anth10_cons_sf, fun=sum)/v_all10$extract(anth10_cons_sf_out, fun=sum) #0.2034491

v_forage10$extract(anth10_lib_sf_out, fun=sum) #1257663
v_forage10$extract(anth10_mod_sf_out, fun=sum) #1730263
v_forage10$extract(anth10_cons_sf_out, fun=sum) #1824840

v_forage10$extract(anth10_lib_sf, fun=sum)/v_forage10$extract(anth10_lib_sf_out, fun=sum) #0.6060254
v_forage10$extract(anth10_mod_sf, fun=sum)/v_forage10$extract(anth10_mod_sf_out, fun=sum) #0.1673592
v_forage10$extract(anth10_cons_sf, fun=sum)/v_forage10$extract(anth10_cons_sf_out, fun=sum) #0.1068579

#######################################################################
###################### Overlap Analysis (3 km) ########################
#######################################################################

anth09 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Anthrax_GARP/Rerun_3km/B.anthracis_2009_3km.tif")
anth10 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Anthrax_GARP/Rerun_3km/B.anthracis_2010_3km.tif")

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
# pdf('Anthrax_2009_Liberal_3k.pdf', height=6, width=8)
# plot(anth09_liberal)
# dev.off()


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
# pdf('Anthrax_2009_Moderate_3k.pdf', height=6, width=8)
# plot(anth09_moderate)
# dev.off()

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
# pdf('Anthrax_2009_Conservative_3k.pdf', height=6, width=8)
# plot(anth09_conserve)
# dev.off()

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
# pdf('Anthrax_2010_Liberal_3k.pdf', height=6, width=8)
# plot(anth10_liberal)
# dev.off()

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
# pdf('Anthrax_2010_Moderate_3k.pdf', height=6, width=8)
# plot(anth10_moderate)
# dev.off()

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
# pdf('Anthrax_2010_Conservative_3k.pdf', height=6, width=8)
# plot(anth10_conserve)
# dev.off()

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

v_all09$extract(anth09_lib_sf, fun=sum) #654251.5
v_all09$extract(anth09_mod_sf, fun=sum) #257675.7
v_all09$extract(anth09_cons_sf, fun=sum) #153164.8

v_forage09$extract(anth09_lib_sf, fun=sum) #742450.5
v_forage09$extract(anth09_mod_sf, fun=sum) #275962
v_forage09$extract(anth09_cons_sf, fun=sum) #155764.7

anth10_lib_poly <- SpatialPolygons(list(rasterToPolygons(anth10_liberal, dissolve=TRUE)@polygons[[2]]))
anth10_lib_sf <- st_as_sf(anth10_lib_poly)
anth10_mod_poly <- SpatialPolygons(list(rasterToPolygons(anth10_moderate, dissolve=TRUE)@polygons[[2]]))
anth10_mod_sf <- st_as_sf(anth10_mod_poly)
anth10_cons_poly <- SpatialPolygons(list(rasterToPolygons(anth10_conserve, dissolve=TRUE)@polygons[[2]]))
anth10_cons_sf <- st_as_sf(anth10_cons_poly)

v_all10$extract(anth10_lib_sf, fun=sum) #597350
v_all10$extract(anth10_mod_sf, fun=sum) #244124.2
v_all10$extract(anth10_cons_sf, fun=sum) #135443.6

v_forage10$extract(anth10_lib_sf, fun=sum) #734826.9
v_forage10$extract(anth10_mod_sf, fun=sum) #264268.7
v_forage10$extract(anth10_cons_sf, fun=sum) #143562.5


##############################################################

anth09_lib_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_liberal, dissolve=TRUE)@polygons[[1]]))
anth09_lib_sf_out <- st_as_sf(anth09_lib_poly_out)
anth09_mod_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_moderate, dissolve=TRUE)@polygons[[1]]))
anth09_mod_sf_out <- st_as_sf(anth09_mod_poly_out)
anth09_cons_poly_out <- SpatialPolygons(list(rasterToPolygons(anth09_conserve, dissolve=TRUE)@polygons[[1]]))
anth09_cons_sf_out <- st_as_sf(anth09_cons_poly_out)

v_all09$extract(anth09_lib_sf_out, fun=sum) #1676823
v_all09$extract(anth09_mod_sf_out, fun=sum) #2073399
v_all09$extract(anth09_cons_sf_out, fun=sum) #2177910

v_all09$extract(anth09_lib_sf, fun=sum)/v_all09$extract(anth09_lib_sf_out, fun=sum) #0.3901733
v_all09$extract(anth09_mod_sf, fun=sum)/v_all09$extract(anth09_mod_sf_out, fun=sum) #0.124277
v_all09$extract(anth09_cons_sf, fun=sum)/v_all09$extract(anth09_cons_sf_out, fun=sum) #0.07032651

v_forage09$extract(anth09_lib_sf_out, fun=sum) #2117500
v_forage09$extract(anth09_mod_sf_out, fun=sum) #2583988
v_forage09$extract(anth09_cons_sf_out, fun=sum) #2704186

v_forage09$extract(anth09_lib_sf, fun=sum)/v_forage09$extract(anth09_lib_sf_out, fun=sum) #0.350626
v_forage09$extract(anth09_mod_sf, fun=sum)/v_forage09$extract(anth09_mod_sf_out, fun=sum) #0.1067969
v_forage09$extract(anth09_cons_sf, fun=sum)/v_forage09$extract(anth09_cons_sf_out, fun=sum) #0.05760132

anth10_lib_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_liberal, dissolve=TRUE)@polygons[[1]]))
anth10_lib_sf_out <- st_as_sf(anth10_lib_poly_out)
anth10_mod_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_moderate, dissolve=TRUE)@polygons[[1]]))
anth10_mod_sf_out <- st_as_sf(anth10_mod_poly_out)
anth10_cons_poly_out <- SpatialPolygons(list(rasterToPolygons(anth10_conserve, dissolve=TRUE)@polygons[[1]]))
anth10_cons_sf_out <- st_as_sf(anth10_cons_poly_out)

v_all10$extract(anth10_lib_sf_out, fun=sum) #600670.9
v_all10$extract(anth10_mod_sf_out, fun=sum) #953896.8
v_all10$extract(anth10_cons_sf_out, fun=sum) #1062577

v_all10$extract(anth10_lib_sf, fun=sum)/v_all10$extract(anth10_lib_sf_out, fun=sum) #0.9944714
v_all10$extract(anth10_mod_sf, fun=sum)/v_all10$extract(anth10_mod_sf_out, fun=sum) #0.2559231
v_all10$extract(anth10_cons_sf, fun=sum)/v_all10$extract(anth10_cons_sf_out, fun=sum) #0.127467

v_forage10$extract(anth10_lib_sf_out, fun=sum) #1290121
v_forage10$extract(anth10_mod_sf_out, fun=sum) #1760679
v_forage10$extract(anth10_cons_sf_out, fun=sum) #1881385

v_forage10$extract(anth10_lib_sf, fun=sum)/v_forage10$extract(anth10_lib_sf_out, fun=sum) #0.5695798
v_forage10$extract(anth10_mod_sf, fun=sum)/v_forage10$extract(anth10_mod_sf_out, fun=sum) #0.1500948
v_forage10$extract(anth10_cons_sf, fun=sum)/v_forage10$extract(anth10_cons_sf_out, fun=sum) #0.07630683

