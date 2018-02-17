library(raster)
library(tidyverse)
library(lubridate)
library(biomod2)
library(sf)
library(dismo)

#########################################################

road_dens <- raster('ENP_Predictors/Road_Density.tif')

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

km_rast <- raster(extent(road_dens), res=c(1000,1000), crs=road_dens)

bio1 <- raster('Anthrax_GARP/wc2.0_bio_30s_01.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(bio1, "Anthrax_GARP/bio1_1km.tif")

bio7 <- raster('Anthrax_GARP/wc2.0_bio_30s_07.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(bio7, "Anthrax_GARP/bio7_1km.tif")

bio12 <- raster('Anthrax_GARP/wc2.0_bio_30s_12.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(bio12, "Anthrax_GARP/bio12_1km.tif")

bio13 <- raster('Anthrax_GARP/wc2.0_bio_30s_13.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(bio13, "Anthrax_GARP/bio13_1km.tif")

pH <- raster('Anthrax_GARP/pH_H20_250m.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(pH, "Anthrax_GARP/pH_H20_1km.tif")

OC <- raster('Anthrax_GARP/OrganicContent_250m.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(OC, "Anthrax_GARP/OrganicContent_1km.tif")

CEC <- raster('Anthrax_GARP/CationExchangeCapacity_250m.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(CEC, "Anthrax_GARP/CationExchangeCapacity_1km.tif")

Green_2009 <- raster('ENP_Predictors/Mean_Greenness_2009.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Green_2009, "Anthrax_GARP/Greenness_2009_1km.tif")

Green_2010 <- raster('ENP_Predictors/Mean_Greenness_2010.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Green_2010, "Anthrax_GARP/Greenness_2010_1km.tif")

Wet_2009 <- raster('ENP_Predictors/Mean_Wetness_2009.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Wet_2009, "Anthrax_GARP/Wetness_2009_1km.tif")

Wet_2010 <- raster('ENP_Predictors/Mean_Wetness_2010.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Wet_2010, "Anthrax_GARP/Wetness_2010_1km.tif")

Mean_NDVI_2009 <- raster('ENP_Predictors/Mean_NDVI_2009.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Mean_NDVI_2009, "Anthrax_GARP/Mean_NDVI_2009_1km.tif")

Max_NDVI_2009 <- raster('ENP_Predictors/Max_NDVI_2009.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Max_NDVI_2009, "Anthrax_GARP/Max_NDVI_2009_1km.tif")

Min_NDVI_2009 <- raster('ENP_Predictors/Min_NDVI_2009.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Min_NDVI_2009, "Anthrax_GARP/Min_NDVI_2009_1km.tif")

Range_NDVI_2009 <- raster('ENP_Predictors/Range_NDVI_2009.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Range_NDVI_2009, "Anthrax_GARP/Range_NDVI_2009_1km.tif")

Mean_NDVI_2010 <- raster('ENP_Predictors/Mean_NDVI_2010.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Mean_NDVI_2010, "Anthrax_GARP/Mean_NDVI_2010_1km.tif")

Max_NDVI_2010 <- raster('ENP_Predictors/Max_NDVI_2010.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Max_NDVI_2010, "Anthrax_GARP/Max_NDVI_2010_1km.tif")

Min_NDVI_2010 <- raster('ENP_Predictors/Min_NDVI_2010.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Min_NDVI_2010, "Anthrax_GARP/Min_NDVI_2010_1km.tif")

Range_NDVI_2010 <- raster('ENP_Predictors/Range_NDVI_2010.tif') %>%
  projectRaster(km_rast) %>% crop(zebra.ext) %>%
  resample(km_rast)
writeRaster(Range_NDVI_2010, "Anthrax_GARP/Range_NDVI_2010_1km.tif")

#########################################################

pH <- raster('Anthrax_GARP/pH_H20_1km.tif')
OC <- raster('Anthrax_GARP/OrganicContent_1km.tif')
CEC <- raster('Anthrax_GARP/CationExchangeCapacity_1km.tif')
bio1 <- raster("Anthrax_GARP/bio1_1km.tif")
bio7 <- raster("Anthrax_GARP/bio7_1km.tif")
bio12 <- raster("Anthrax_GARP/bio12_1km.tif")
bio13 <- raster("Anthrax_GARP/bio13_1km.tif")
mean_NDVI_2009 <- raster('Anthrax_GARP/Mean_NDVI_2009_1km.tif')
max_NDVI_2009 <- raster('Anthrax_GARP/Max_NDVI_2009_1km.tif')
min_NDVI_2009 <- raster('Anthrax_GARP/Min_NDVI_2009_1km.tif')
range_NDVI_2009 <- raster('Anthrax_GARP/Range_NDVI_2009_1km.tif')
wet_2009 <- raster('Anthrax_GARP/Wetness_2009_1km.tif')
green_2009 <- raster('Anthrax_GARP/Greenness_2009_1km.tif')
mean_NDVI_2010 <- raster('Anthrax_GARP/Mean_NDVI_2010_1km.tif')
max_NDVI_2010 <- raster('Anthrax_GARP/Max_NDVI_2010_1km.tif')
min_NDVI_2010 <- raster('Anthrax_GARP/Min_NDVI_2010_1km.tif')
range_NDVI_2010 <- raster('Anthrax_GARP/Range_NDVI_2010_1km.tif')
wet_2010 <- raster('Anthrax_GARP/Wetness_2010_1km.tif')
green_2010 <- raster('Anthrax_GARP/Greenness_2010_1km.tif')

carcass <- st_read('/Users/ericdougherty/Dropbox/EricDanaRSF/GIS/Etosha/AllAnthraxCarcasses09-10.shp')
#carcass <- st_read('/Anthrax_GARP/AllAnthraxCarcasses09-10.shp') %>%
carcass <- st_transform(carcass, '+init=epsg:32733') %>%
  mutate(Year = year(DATE))
carc.2009.sp <- carcass %>% filter(Year == 2009) %>% 
  mutate(resp.var = 1) %>% dplyr::select(resp.var) %>% as('Spatial') #%>% points(cex=0.3, pch=19)
carc.2010.sp <- carcass %>% filter(Year == 2010) %>% 
  mutate(resp.var = 1) %>% dplyr::select(resp.var) %>% as('Spatial') #%>% points(cex=0.3, pch=19)

expl_var_2009 <- stack (pH, OC, CEC, bio1, bio7, bio12,
                        bio13, mean_NDVI_2009, 
                        max_NDVI_2009, min_NDVI_2009,
                        range_NDVI_2009, wet_2009)
cov_2009 <- layerStats(expl_var_2009, 'pearson', na.rm=TRUE)
write.csv(cov_2009[[1]], "Anthrax_GARP/Covariance_Matrix_2009.csv")

expl_var_2010 <- stack (pH, OC, CEC, bio1, bio7, bio12,
                        bio13, mean_NDVI_2010, 
                        max_NDVI_2010, min_NDVI_2010,
                        range_NDVI_2010, wet_2010)
cov_2010 <- layerStats(expl_var_2010, 'pearson', na.rm=TRUE)
write.csv(cov_2010[[1]], "Anthrax_GARP/Covariance_Matrix_2010.csv")

random <- data.frame(randomPoints(expl_var_2009, 2*(nrow(carc.2009.sp))))
random$resp.var <- 0
true_2009 <- data.frame(carc.2009.sp@coords)
true_2009$resp.var <- 1
colnames(true_2009) <- c('x', 'y', 'resp.var')
bound_2009 <- rbind(true_2009, random)

random <- data.frame(randomPoints(expl_var_2010, 2*(nrow(carc.2010.sp))))
random$resp.var <- 0
true_2010 <- data.frame(carc.2010.sp@coords)
true_2010$resp.var <- 1
colnames(true_2010) <- c('x', 'y', 'resp.var')
bound_2010 <- rbind(true_2010, random)

###################################################################
###################################################################

biomod_2009 <- BIOMOD_FormatingData(resp.var = bound_2009$resp.var, 
                                    expl.var = expl_var_2009, 
                                    resp.xy = bound_2009[,c('x','y')],
                                    resp.name = 'B.anthracis_2009')

######## Run a single set of models to evaluate predictors ########
setwd('/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Anthrax_GARP/')
myBiomodModelOut <- BIOMOD_Modeling(data = biomod_2009,
                            models=c('GLM', 'GBM', 'GAM', 'CTA', 
                                     'ANN', 'SRE', 'FDA', 'MARS',
                                     'RF', 'MAXENT.Phillips'),
                            NbRunEval=1,
                            VarImport=1)
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                       chosen.models = 'all',
                                       em.by = 'all',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS','ROC'),
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional' )
get_evaluations(myBiomodEM)
var_import <- get_variables_importance(myBiomodModelOut)
var_import <- data.frame(var_import[,,1,1])
write.csv(var_import, "Variable_Importance_2009.csv")

# bio 13 is eliminated (covariance with bio12) - 0.187 vs. 0.132
# NDVI measures are eliminated (covariance with Wetness) - 0.416 vs. < 0.152
expl_var_2009_new <- stack(pH, OC, CEC, bio1, 
                           bio7, bio12, wet_2009)

######## Run all models five times to evaluate methods ########
myBiomodData <- BIOMOD_FormatingData(resp.var = bound_2009$resp.var, 
                                     expl.var = expl_var_2009_new, 
                                     resp.xy = bound_2009[,c('x','y')],
                                     resp.name = "B.anthracis_2009")
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,models = c('GLM', 'GBM', 'GAM', 'CTA', 
                                                             'ANN', 'SRE', 'FDA', 'MARS', 'RF'), 
                                     NbRunEval=5, 
                                     DataSplit = 80, 
                                     VarImport=5)
myBiomodModelOut

myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                       em.by = 'all',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS','ROC'),
                                       chosen.models = 'all', 
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
myBiomodEM
get_evaluations(myBiomodEM)

myBiomodProjection <- BIOMOD_Projection(myBiomodModelOut, expl_var_2009_new, 'ensemble')
myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection,
                                        EM.output = myBiomodEM)

myBiomodModelEval

myBiomodModelEval["ROC","Testing.data",,,]

get_variables_importance(myBiomodModelOut)

####### Run final predictors and methods 10 times #######

myBiomodData <- BIOMOD_FormatingData(resp.var = bound_2009$resp.var, 
                                     expl.var = expl_var_2009_new, 
                                     resp.xy = bound_2009[,c('x','y')],
                                     resp.name = "B.anthracis_2009")
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c('GLM', 'GBM', 'GAM', 'ANN',
                                               'FDA', 'MARS', 'RF'), 
                                    NbRunEval=10, 
                                    DataSplit = 80, 
                                    VarImport=10)
myBiomodModelOut
  
myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      em.by = 'all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.7),
                                      models.eval.meth = c('TSS','ROC'),
                                      chosen.models = c('B.anthracis.2009_AllData_Full_GLM', 
                                                        'B.anthracis.2009_AllData_Full_GBM', 
                                                        'B.anthracis.2009_AllData_Full_GAM', 
                                                        'B.anthracis.2009_AllData_Full_ANN', 
                                                        'B.anthracis.2009_AllData_Full_FDA', 
                                                        'B.anthracis.2009_AllData_Full_MARS', 
                                                        'B.anthracis.2009_AllData_Full_RF'), 
                                      prob.mean = TRUE,
                                      prob.cv = FALSE,
                                      prob.ci = FALSE,
                                      prob.ci.alpha = 0.05,
                                      prob.median = FALSE,
                                      committee.averaging = FALSE,
                                      prob.mean.weight = TRUE,
                                      prob.mean.weight.decay = 'proportional' )
myBiomodEM
get_evaluations(myBiomodEM)
  
myBiomodProjection <- BIOMOD_Projection(myBiomodModelOut, expl_var_2009_new, 'ensemble')
myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection,
                                         EM.output = myBiomodEM)

plot(myBiomodProjection, str.grep = "Full")
plot(myBiomodEF)

###################################################################
###################################################################

biomod_2010 <- BIOMOD_FormatingData(resp.var = bound_2010$resp.var, 
                                    expl.var = expl_var_2010, 
                                    resp.xy = bound_2010[,c('x','y')],
                                    resp.name = 'B.anthracis_2010')

######## Run a single set of models to evaluate predictors ########
setwd('/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Anthrax_GARP/')
myBiomodModelOut <- BIOMOD_Modeling(data = biomod_2010,
                                    models=c('GLM', 'GBM', 'GAM', 'CTA', 
                                             'ANN', 'SRE', 'FDA', 'MARS', 'RF'),
                                    NbRunEval=1,
                                    VarImport=1)
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                       chosen.models = 'all',
                                       em.by = 'all',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS','ROC'),
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional' )
get_evaluations(myBiomodEM)
var_import <- get_variables_importance(myBiomodModelOut)
var_import <- data.frame(var_import[,,1,1])
write.csv(var_import, "Variable_Importance_2010.csv")

# bio 13 is eliminated (covariance with bio12) - 0.352 vs. 0.237
# mean NDVI is eliminated (covariance with Wet) - 0.355 vs. < 0.038
# range NDVI is eliminated (covariance with min) - 0.033 vs. 0.015
expl_var_2010_new <- stack(pH, OC, CEC, bio1, 
                           bio7, bio12, max_NDVI_2010,
                           min_NDVI_2010, wet_2010)

######## Run all models five times to evaluate methods ########
myBiomodData <- BIOMOD_FormatingData(resp.var = bound_2010$resp.var, 
                                     expl.var = expl_var_2010_new, 
                                     resp.xy = bound_2010[,c('x','y')],
                                     resp.name = "B.anthracis_2010")
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,models = c('GLM', 'GBM', 'GAM', 'CTA', 
                                                             'ANN', 'SRE', 'FDA', 'MARS', 'RF'), 
                                     NbRunEval=5, 
                                     DataSplit = 80, 
                                     VarImport=5)
myBiomodModelOut

myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                       em.by = 'all',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       models.eval.meth = c('TSS','ROC'),
                                       chosen.models = 'all', 
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')
myBiomodEM
get_evaluations(myBiomodEM)

myBiomodProjection <- BIOMOD_Projection(myBiomodModelOut, expl_var_2010_new, 'ensemble')
myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection,
                                         EM.output = myBiomodEM)

myBiomodModelEval

myBiomodModelEval["ROC","Testing.data",,,]

get_variables_importance(myBiomodModelOut)

####### Run final predictors and methods 10 times #######

myBiomodData <- BIOMOD_FormatingData(resp.var = bound_2010$resp.var, 
                                     expl.var = expl_var_2010_new, 
                                     resp.xy = bound_2010[,c('x','y')],
                                     resp.name = "B.anthracis_2010")
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c('GLM', 'GBM', 'GAM', 'ANN',
                                               'FDA', 'MARS', 'RF'), 
                                    NbRunEval=10, 
                                    DataSplit = 80, 
                                    VarImport=10)
myBiomodModelOut

myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      em.by = 'all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.7),
                                      models.eval.meth = c('TSS','ROC'),
                                      chosen.models = c('B.anthracis.2010_AllData_Full_GLM', 
                                                        'B.anthracis.2010_AllData_Full_GBM', 
                                                        'B.anthracis.2010_AllData_Full_GAM', 
                                                        'B.anthracis.2010_AllData_Full_ANN', 
                                                        'B.anthracis.2010_AllData_Full_FDA', 
                                                        'B.anthracis.2010_AllData_Full_MARS', 
                                                        'B.anthracis.2010_AllData_Full_RF'), 
                                      prob.mean = TRUE,
                                      prob.cv = FALSE,
                                      prob.ci = FALSE,
                                      prob.ci.alpha = 0.05,
                                      prob.median = FALSE,
                                      committee.averaging = FALSE,
                                      prob.mean.weight = TRUE,
                                      prob.mean.weight.decay = 'proportional' )
myBiomodEM
get_evaluations(myBiomodEM)

myBiomodProjection <- BIOMOD_Projection(myBiomodModelOut, expl_var_2010_new, 'ensemble')
myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection,
                                         EM.output = myBiomodEM)

pdf('B.anthracis_2010_All_Models.pdf', width=12, height=10)
plot(myBiomodProjection, str.grep = "Full")
dev.off()

pdf('B.anthracis_2010_Ensemble.pdf', width=8, height=10)
plot(myBiomodEF)
dev.off()

BiomodEF <- raster("Anthrax_GARP/B.anthracis.2010/proj_ensemble/proj_ensemble_B.anthracis.2010_ensemble.grd")

writeRaster(BiomodEF, "B.anthracis_2010.tif", format = "GTiff")