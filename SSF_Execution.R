library(raster)
library(tidyverse)
library(sf)
library(Epi)
library(survival)

setwd("~/Box Sync/Dissertation/Behavioral_SSF")

########################################################

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

zebra09 <- read_csv("Zebra_Anthrax_2009_Cleaned.csv") %>%
    dplyr::select(x,y,date,ID) %>%
    dplyr::filter(!is.na(x)) %>%
    mutate(row =  seq(1,nrow(.),1), dists = 0, weights = 1) %>%
    st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

zebra10 <- read_csv("Zebra_Anthrax_2010_Cleaned.csv") %>%
    dplyr::select(x,y,date,ID) %>%
    dplyr::filter(!is.na(x)) %>%
    mutate(row =  seq(1,nrow(.),1), dists = 0, weights = 1) %>%
    st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

all_zebra <- rbind(zebra09, zebra10)

zebra.ext <- extent(c(extent(zebra10)@xmin - 5000,
                      extent(zebra10)@xmax + 5000,
                      extent(zebra09)@ymin - 5000,
                      extent(zebra10)@ymax + 5000))

#########################################################

road_dens <- raster('ENP_Predictors/Road_Density.tif')
NDVI_09 <- raster('ENP_Predictors/Mean_NDVI_2009.tif') 
NDVI_10 <- raster('ENP_Predictors/Mean_NDVI_2010.tif') 
Green_09 <- raster('ENP_Predictors/Mean_Greenness_2009.tif')
Green_10 <- raster('ENP_Predictors/Mean_Greenness_2010.tif')
Wet_09 <- raster('ENP_Predictors/Mean_Wetness_2009.tif')
Wet_10 <- raster('ENP_Predictors/Mean_Wetness_2010.tif')
# veg <- raster('ENP_Predictors/Vegetation_Crop.tif')
# bare <- raster('ENP_Predictors/Prop_Bare.tif')
# steppe <- raster('ENP_Predictors/Prop_Steppe.tif')
# grass <- raster('ENP_Predictors/Prop_Grass.tif')
# shrub<- raster('ENP_Predictors/Prop_Shrub.tif')
# low <- raster('ENP_Predictors/Prop_Low.tif')
# high <- raster('ENP_Predictors/Prop_High.tif')

pred_stack_09 <- stack(NDVI_09, Green_09, Wet_09, road_dens)#, bare, steppe, grass, shrub, low, high)
pred_stack_10 <- stack(NDVI_10, Green_10, Wet_10, road_dens)#, bare, steppe, grass, shrub, low, high)

means09 <- cellStats(pred_stack_09, stat='mean', na.rm=TRUE)
sd09 <- cellStats(pred_stack_09, stat='sd', na.rm=TRUE)
means10 <- cellStats(pred_stack_10, stat='mean', na.rm=TRUE)
sd10 <- cellStats(pred_stack_10, stat='sd', na.rm=TRUE)

norm_stack_09 <- stack((Green_09 - means09[2]) / sd09[2],
                       (Wet_09 - means09[3]) / sd09[3],
                       (road_dens - means09[4])/ sd09[4])
#                       bare, steppe, grass, shrub, low, high)

norm_stack_10 <- stack((Green_10 - means10[2]) / sd10[2],
                       (Wet_10 - means10[3]) / sd10[3],
                       (road_dens - means10[4])/ sd10[4])
#                       bare, steppe, grass, shrub, low, high)


############################################################

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Original_Final")

cond.logit.all <- function(i) {
  if (i < 6) {
    all_final <- read_csv(paste0(name_list[i],"_All_Final.csv")) %>%
      dplyr::select("Green", "Wet", "Road_Dens",
                    "Green_Norm", "Wet_Norm", "Road_Dens_Norm", 
                    "ID", "binom", "fix")
  } else {
    all_final <- read_csv(paste0(name_list[i],"_All_Final.csv")) %>%
        dplyr::select("Green", "Wet", "Road_Dens",
                      "Green_Norm", "Wet_Norm", "Road_Dens_Norm", 
                      "ID", "binom", "fix")
  }
  assign(paste0(name_list[i],"_All"), all_final)
  assign(paste0(name_list[i],"_Corr"), cor(all_final[,1:4]))
  assign(paste0(name_list[i],"_Model"), clogistic(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm,
  #                                                  prop_bare + prop_steppe + prop_grass + 
  #                                                  prop_shrub + prop_low + prop_high,
                                                  data=all_final, strata=fix))
  assign(paste0(name_list[i],"_clogit"), clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm +
  #                                                 prop_bare + prop_steppe + prop_grass + 
  #                                                 prop_shrub + prop_low + prop_high + 
                                                    strata(fix),
                                                 data=all_final))
  assign(paste0(name_list[i],"_GLM"), glm(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm,
  #                                          prop_bare + prop_steppe + prop_grass + 
  #                                          prop_shrub + prop_low + prop_high,
                                          data=all_final, family=binomial(link='logit')))

  out_list <- list(get(paste0(name_list[i],"_All")),
                   get(paste0(name_list[i],"_Corr")),
                   get(paste0(name_list[i],"_Model")),
                   get(paste0(name_list[i],"_clogit")),
                   get(paste0(name_list[i],"_GLM")))
  return(out_list)
}

cond.logit.forage <- function(i) {
  if (i < 6) {
    forage_final <- read_csv(paste0(name_list[i],"_Foraging_Final.csv")) %>%
      dplyr::select("Green", "Wet", "Road_Dens",
                    "Green_Norm", "Wet_Norm", "Road_Dens_Norm", 
                    "ID", "binom", "fix")
  } else {
    forage_final <- read_csv(paste0(name_list[i],"_Foraging_Final.csv")) %>%
      dplyr::select("Green", "Wet", "Road_Dens",
                    "Green_Norm", "Wet_Norm", "Road_Dens_Norm", 
                    "ID", "binom", "fix")
  }
  all_pts <- rbind(used, avail)
  assign(paste0(name_list[i],"_All"), forage_final)
  assign(paste0(name_list[i],"_Corr"), cor(forage_final[,1:4]))
  assign(paste0(name_list[i],"_Model"), clogistic(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm,
  #                                                  prop_bare + prop_steppe + prop_grass + 
  #                                                  prop_shrub + prop_low + prop_high,
                                                  data=forage_final, strata=fix))
  assign(paste0(name_list[i],"_clogit"), clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm + 
  #                                                 prop_bare + prop_steppe + prop_grass + 
  #                                                 prop_shrub + prop_low + prop_high + 
                                                    strata(fix), 
                                                 data=forage_final))
  assign(paste0(name_list[i],"_GLM"), glm(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm,
  #                                          prop_bare + prop_steppe + prop_grass + 
  #                                          prop_shrub + prop_low + prop_high,
                                          data=forage_final, family=binomial(link='logit')))

  out_list <- list(get(paste0(name_list[i],"_All")),
                   get(paste0(name_list[i],"_Corr")),
                   get(paste0(name_list[i],"_Model")),
                   get(paste0(name_list[i],"_clogit")),
                   get(paste0(name_list[i],"_GLM")))
  return(out_list)
}

out1 <- cond.logit.all(1)
out2 <- cond.logit.all(2)
out3 <- cond.logit.all(3)
out4 <- cond.logit.all(4)
out5 <- cond.logit.all(5)
out6 <- cond.logit.all(6)
out7 <- cond.logit.all(7)
out8 <- cond.logit.all(8)
out9 <- cond.logit.all(9)
out10 <- cond.logit.all(10)
out11 <- cond.logit.all(11)

forage1 <- cond.logit.forage(1)
forage2 <- cond.logit.forage(2)
forage3 <- cond.logit.forage(3)
forage4 <- cond.logit.forage(4)
forage5 <- cond.logit.forage(5)
forage6 <- cond.logit.forage(6)
forage7 <- cond.logit.forage(7)
forage8 <- cond.logit.forage(8)
forage9 <- cond.logit.forage(9)
forage10 <- cond.logit.forage(10)
forage11 <- cond.logit.forage(11)

#################################################################
################# Epi clogistic Predictions #####################

pred_fun09 <- function(model, object) {
  coef <- model$coefficients
  out_vals <- ((coef[1]*object$Mean_Wetness_2009) + (coef[2]*object$Mean_Greenness_2009) + 
              (coef[3]*object$Road_Density))
  return(out_vals)
}

pred_fun10 <- function(model, object) {
  coef <- model$coefficients
  out_vals <- ((coef[1]*object$Mean_Wetness_2010) + (coef[2]*object$Mean_Greenness_2010) +
              (coef[3]*object$Road_Density))
  return(out_vals)
}

pred1 <- exp(raster::predict(object=norm_stack_09, model=out1[[3]], fun=pred_fun09))
pred2 <- exp(raster::predict(object=norm_stack_09, model=out2[[3]], fun=pred_fun09))
pred3 <- exp(raster::predict(object=norm_stack_09, model=out3[[3]], fun=pred_fun09))
pred4 <- exp(raster::predict(object=norm_stack_09, model=out4[[3]], fun=pred_fun09))
pred5 <- exp(raster::predict(object=norm_stack_09, model=out5[[3]], fun=pred_fun09))
pred6 <- exp(raster::predict(object=norm_stack_10, model=out6[[3]], fun=pred_fun10))
pred7 <- exp(raster::predict(object=norm_stack_10, model=out7[[3]], fun=pred_fun10))
pred8 <- exp(raster::predict(object=norm_stack_10, model=out8[[3]], fun=pred_fun10))
pred9 <- exp(raster::predict(object=norm_stack_10, model=out9[[3]], fun=pred_fun10))
pred10 <- exp(raster::predict(object=norm_stack_10, model=out10[[3]], fun=pred_fun10))
pred11 <- exp(raster::predict(object=norm_stack_10, model=out11[[3]], fun=pred_fun10))

pred1_forage <- exp(raster::predict(object=norm_stack_09, model=forage1[[3]], fun=pred_fun09))
pred2_forage <- exp(raster::predict(object=norm_stack_09, model=forage2[[3]], fun=pred_fun09))
pred3_forage <- exp(raster::predict(object=norm_stack_09, model=forage3[[3]], fun=pred_fun09))
pred4_forage <- exp(raster::predict(object=norm_stack_09, model=forage4[[3]], fun=pred_fun09))
pred5_forage <- exp(raster::predict(object=norm_stack_09, model=forage5[[3]], fun=pred_fun09))
pred6_forage <- exp(raster::predict(object=norm_stack_10, model=forage6[[3]], fun=pred_fun10))
pred7_forage <- exp(raster::predict(object=norm_stack_10, model=forage7[[3]], fun=pred_fun10))
pred8_forage <- exp(raster::predict(object=norm_stack_10, model=forage8[[3]], fun=pred_fun10))
pred9_forage <- exp(raster::predict(object=norm_stack_10, model=forage9[[3]], fun=pred_fun10))
pred10_forage <- exp(raster::predict(object=norm_stack_10, model=forage10[[3]], fun=pred_fun10))
pred11_forage <- exp(raster::predict(object=norm_stack_10, model=forage11[[3]], fun=pred_fun10))

##########################################################
########### survival clogit Predictions ##################

norm_stack_09_GLM <- stack(norm_stack_09@layers[[1]],
                           norm_stack_09@layers[[2]],
                           norm_stack_09@layers[[3]])
names(norm_stack_09_GLM) <- c("Green_Norm", "Wet_Norm", "Road_Dens_Norm")
norm_stack_10_GLM <- stack(norm_stack_10@layers[[1]],
                           norm_stack_10@layers[[2]],
                           norm_stack_10@layers[[3]])
names(norm_stack_10_GLM) <- c("Green_Norm", "Wet_Norm", "Road_Dens_Norm")

pred_09_df <- as.data.frame(norm_stack_09_GLM)
pred_09_df$fix <- 1
pred_10_df <- as.data.frame(norm_stack_10_GLM)
pred_10_df$fix <- 1

for (i in 1:length(name_list)) {
  all <- raster(Wet_09)
  all[] <- 0
  forage <- raster(Wet_09)
  forage[] <- 0
  
  if (i < 6) {
    pred <- predict(object=get(paste0("out",i))[[4]], newdata=pred_09_df, type='risk')
    pred_forage <- predict(object=get(paste0("forage",i))[[4]], newdata=pred_09_df, type='risk')
  } else {
    pred <- predict(object=get(paste0("out",i))[[4]], newdata=pred_10_df, type='risk')
    pred_forage <- predict(object=get(paste0("forage",i))[[4]], newdata=pred_10_df, type='risk')
  }
  all@data@values <- pred
  all <- all/(1 + all)
  forage@data@values <- pred_forage
  forage <- forage/(1 + forage)
  
  assign(paste0(name_list[i],"_All_Risk"), all)
  assign(paste0(name_list[i],"_Forage_Risk"), forage)
  
  print(i)
}

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Rerun_Results")
for (i in 1:length(name_list)) {
  writeRaster(get(paste0(name_list[i],"_All_Risk")), paste0(name_list[i],"_All_Risk.tif"), format='GTiff')
  writeRaster(get(paste0(name_list[i],"_Forage_Risk")), paste0(name_list[i],"_Forage_Risk.tif"), format='GTiff')
}

#####################################################################
##################### All Zebra Predictions #########################

read_all_final <- function(file_name) {
  output <- read_csv(paste0(file_name,"_All_Final.csv")) %>%
    dplyr::select("Green", "Wet", "Road_Dens",
                  "Green_Norm", "Wet_Norm", "Road_Dens_Norm", 
                  "ID", "binom", "fix")
  return(output)
}

read_forage_final <- function(file_name) {
  output <- read_csv(paste0(file_name,"_Foraging_Final.csv")) %>%
    dplyr::select("Green", "Wet", "Road_Dens",
                  "Green_Norm", "Wet_Norm", "Road_Dens_Norm", 
                  "ID", "binom", "fix")
  return(output)
}

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Rerun_Extraction_Results")
all_2009 <- read_all_final(name_list[1])
for (i in 2:5) {
  temp <- read_all_final(name_list[i])
  all_2009 <- rbind(all_2009, temp)
}
all_2010 <- read_all_final(name_list[6])
for (i in 7:11) {
  temp <- read_all_final(name_list[i])
  all_2010 <- rbind(all_2010, temp)
}

forage_2009 <- read_forage_final(name_list[1])
for (i in 2:5) {
  temp <- read_forage_final(name_list[i])
  forage_2009 <- rbind(forage_2009, temp)
}
forage_2010 <- read_forage_final(name_list[6])
for (i in 7:11) {
  temp <- read_forage_final(name_list[i])
  forage_2010 <- rbind(forage_2010, temp)
}

All_2009_clogit <- clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm + 
                             strata(fix) + strata(ID), data=all_2009)
All_2010_clogit <- clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm + 
                            strata(fix) + strata(ID), data=all_2010)
Foraging_2009_clogit <- clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm + 
                                  strata(fix) + strata(ID), data=forage_2009)
Foraging_2010_clogit <- clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm + 
                                 strata(fix) + strata(ID), data=forage_2010)

pred <- raster(Wet_09)
pred[] <- 0

pred_09_df <- as.data.frame(norm_stack_09_GLM)
pred_09_df$fix <- 1
pred_09_df$ID <- "AG068_2009"

pred_10_df <- as.data.frame(norm_stack_10_GLM)
pred_10_df$fix <- 1
pred_10_df$ID <- "AG068_2010"

all_2009_pred <- predict(object=All_2009_clogit, newdata=pred_09_df, type='risk')
all_2010_pred <- predict(object=All_2010_clogit, newdata=pred_10_df, type='risk')
forage_2009_pred <- predict(object=Foraging_2009_clogit, newdata=pred_09_df, type='risk')
forage_2010_pred <- predict(object=Foraging_2010_clogit, newdata=pred_10_df, type='risk')

pred@data@values <- all_2009_pred
pred <- pred/(1 + pred)
assign("All_2009_Risk", pred)

pred@data@values <- all_2010_pred
pred <- pred/(1 + pred)
assign("All_2010_Risk", pred)

pred@data@values <- forage_2009_pred
pred <- pred/(1 + pred)
assign("Forage_2009_Risk", pred)

pred@data@values <- forage_2010_pred
pred <- pred/(1 + pred)
assign("Forage_2010_Risk", pred)

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Rerun_Results")
writeRaster(All_2009_Risk, "Zebra_2009_All_Risk.tif", format='GTiff')
writeRaster(All_2010_Risk, "Zebra_2010_All_Risk.tif", format='GTiff')
writeRaster(Forage_2009_Risk, "Zebra_2009_Forage_Risk.tif", format='GTiff')
writeRaster(Forage_2010_Risk, "Zebra_2010_Forage_Risk.tif", format='GTiff')

###### Plot clogit Response (All and Foraging) with Points ########

plot(AG059_2009_All_Risk)
zebra09 %>% filter(ID == "AG059_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG059_2009_Forage_Risk)
zebra09 %>% filter(ID == "AG059_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG061_2009_All_Risk)
zebra09 %>% filter(ID == "AG061_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG061_2009_Forage_Risk)
zebra09 %>% filter(ID == "AG061_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG062_2009_All_Risk)
zebra09 %>% filter(ID == "AG062_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG062_2009_Forage_Risk)
zebra09 %>% filter(ID == "AG062_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG063_2009_All_Risk)
zebra09 %>% filter(ID == "AG063_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG063_2009_Forage_Risk)
zebra09 %>% filter(ID == "AG063_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG068_2009_All_Risk)
zebra09 %>% filter(ID == "AG068_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG068_2009_Forage_Risk)
zebra09 %>% filter(ID == "AG068_2009") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG063_2010_All_Risk)
zebra10 %>% filter(ID == "AG063_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG063_2010_Forage_Risk)
zebra10 %>% filter(ID == "AG063_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG068_2010_All_Risk)
zebra10 %>% filter(ID == "AG068_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG068_2010_Forage_Risk)
zebra10 %>% filter(ID == "AG068_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG252_2010_All_Risk)
zebra10 %>% filter(ID == "AG252_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG252_2010_Forage_Risk)
zebra10 %>% filter(ID == "AG252_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG253_2010_All_Risk)
zebra10 %>% filter(ID == "AG253_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG253_2010_Forage_Risk)
zebra10 %>% filter(ID == "AG253_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG255_2010_All_Risk)
zebra10 %>% filter(ID == "AG255_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG255_2010_Forage_Risk)
zebra10 %>% filter(ID == "AG255_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)

plot(AG256_2010_All_Risk)
zebra10 %>% filter(ID == "AG256_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)
plot(AG256_2010_Forage_Risk)
zebra10 %>% filter(ID == "AG256_2010") %>% as('Spatial') %>% points(cex=0.3, pch=19)

##########################################################
########### Compare Predictions by Cell ##################

diff.list <- list()
for (j in 1:length(name_list)){
  if (j < 6) {
    all <- raster::predict(object=norm_stack_09, model=get(paste0("out",j))[[4]], fun=predict, type='response')
    forage <- raster::predict(object=norm_stack_09, model=get(paste0("forage",j))[[4]], fun=predict, type='response')
  } else {
    all <- raster::predict(object=norm_stack_10, model=get(paste0("out",j))[[4]], fun=predict, type='response')
    forage <- raster::predict(object=norm_stack_10, model=get(paste0("forage",j))[[4]], fun=predict, type='response')
  }
  
  all_GLM <- (all - cellStats(all, 'mean')) / cellStats(all, 'sd')
  forage_GLM <- (forage - cellStats(forage, 'mean')) / cellStats(forage, 'sd')
  diff <- c()
  for (i in 1:length(norm1_GLM@data@values)) {
    if (all_GLM@data@values[i] > forage_GLM@data@values[i]) {
      diff[i] <- -1
    } else if (all_GLM@data@values[i] < forage_GLM@data@values[i]) {
      diff[i] <- 1
    } else {
      diff[i] <- 0
    }
  }
  diff.list[[j]] <- diff
}

test <- raster(dist_roads)
test[] <- 0
predictive.values <- c()
for (i in 1:length(name_list)) {
  test@data@values <- diff.list[[i]]
  if (i < 6) {
    pts <- zebra09 %>% filter(ID == name_list[i])
  } else {
    pts <- zebra10 %>% filter(ID == name_list[i])
  }
  vals <- raster::extract(test, pts)
  predictive.values[i] <- sum(vals)
}

