library(raster)
library(tidyverse)
library(sf)
library(velox)

setwd("~/Box Sync/Dissertation/Behavioral_SSF")

########################################################

state.seq <- readRDS('HMM_State_Sequences.rds')

states09 <- data.frame(state.seq[[1]]$states)
names(states09) <- c('states')
for (i in 2:5) {
  states.temp <- data.frame(state.seq[[i]]$states)
  names(states.temp) <- c('states')
  states09 <- rbind(states09, states.temp)
}
states10 <- data.frame(state.seq[[6]]$states)
names(states10) <- c('states')
for (i in 7:11) {
  states.temp <- data.frame(state.seq[[i]]$states)
  names(states.temp) <- c('states')
  states10 <- rbind(states10, states.temp)
}

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

zebra09 <- read_csv("Zebra_Anthrax_2009_Cleaned.csv") %>%
  dplyr::select(x,y,date,ID) %>%
  cbind(states09) %>%
  dplyr::filter(!is.na(x)) %>%
  mutate(row =  seq(1,nrow(.),1)) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

zebra10 <- read_csv("Zebra_Anthrax_2010_Cleaned.csv") %>%
  dplyr::select(x,y,date,ID) %>%
  cbind(states10) %>%
  dplyr::filter(!is.na(x)) %>%
  mutate(row =  seq(1,nrow(.),1)) %>%
  st_as_sf(., coords = 1:2, crs = "+init=epsg:32733")

######################################################

road_dens <- raster('ENP_Predictors/Road_Density.tif')
Green_09 <- raster('ENP_Predictors/Mean_Greenness_2009.tif')
Green_10 <- raster('ENP_Predictors/Mean_Greenness_2010.tif')
Wet_09 <- raster('ENP_Predictors/Mean_Wetness_2009.tif')
Wet_10 <- raster('ENP_Predictors/Mean_Wetness_2010.tif')

pred_stack_09 <- stack(Green_09, Wet_09, road_dens)
pred_stack_10 <- stack(Green_10, Wet_10, road_dens)

rast <- raster(Green_09)
v_green_09 <- velox(setValues(rast, as.matrix(Green_09)))
green_09_extract <- v_green_09$extract_points(sp=zebra09)
v_wet_09 <- velox(setValues(rast, as.matrix(Wet_09)))
wet_09_extract <- v_wet_09$extract_points(sp=zebra09)
v_road <- velox(setValues(rast, as.matrix(road_dens)))
road_09_extract <- v_road$extract_points(sp=zebra09)
#zebra09$Green <- green_09_extract
#zebra09$Wet <- wet_09_extract
#zebra09$Road_Dens <- road_09_extract

v_green_10 <- velox(setValues(rast, as.matrix(Green_10)))
green_10_extract <- v_green_10$extract_points(sp=zebra10)
v_wet_10 <- velox(setValues(rast, as.matrix(Wet_10)))
wet_10_extract <- v_wet_10$extract_points(sp=zebra10)
v_road <- velox(setValues(rast, as.matrix(road_dens)))
road_10_extract <- v_road$extract_points(sp=zebra10)
#zebra10$Green <- green_10_extract
#zebra10$Wet <- wet_10_extract
#zebra10$Road_Dens <- road_10_extract

######################################################

means09 <- cellStats(pred_stack_09, stat='mean', na.rm=TRUE)
sd09 <- cellStats(pred_stack_09, stat='sd', na.rm=TRUE)
means10 <- cellStats(pred_stack_10, stat='mean', na.rm=TRUE)
sd10 <- cellStats(pred_stack_10, stat='sd', na.rm=TRUE)

# norm_stack_09 <- stack((Green_09 - means09[1]) / sd09[1],
#                        (Wet_09 - means09[2]) / sd09[2],
#                        (road_dens - means09[3])/ sd09[3])
# 
# norm_stack_10 <- stack((Green_10 - means10[1]) / sd10[1],
#                        (Wet_10 - means10[2]) / sd10[2],
#                        (road_dens - means10[3])/ sd10[3])


# Create csv of used points for both 'all' and 'foraging' (to be parsed later)
for (i in 1:length(name_list)) {
  if (i < 6) {
    temp <- zebra09 %>% 
      mutate(Green = as.vector(green_09_extract),
             Wet = as.vector(wet_09_extract),
             Road_Dens = as.vector(road_09_extract)) %>%
      filter(ID == name_list[i]) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 1)
    temp <- temp %>% mutate(Green_Norm = (Green - means09['Mean_Greenness_2009']) / sd09['Mean_Greenness_2009'],
             Wet_Norm = (Wet - means09['Mean_Wetness_2009']) / sd09['Mean_Wetness_2009'],
             Road_Dens_Norm = (Road_Dens - means09['Road_Density']) / sd09['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  } else {
    temp <- zebra10 %>% 
      mutate(Green = as.vector(green_10_extract),
             Wet = as.vector(wet_10_extract),
             Road_Dens = as.vector(road_10_extract)) %>%
      filter(ID == name_list[i]) %>%
      mutate(fix =  seq(1,nrow(.),1), binom = 1)
    temp <- temp %>% mutate(Green_Norm = (Green - means10['Mean_Greenness_2010']) / sd10['Mean_Greenness_2010'],
             Wet_Norm = (Wet - means10['Mean_Wetness_2010']) / sd10['Mean_Wetness_2010'],
             Road_Dens_Norm = (Road_Dens - means10['Road_Density']) / sd10['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  }
  write_csv(temp, paste0(name_list[i],"_Foraging_Used.csv"))
  
  #temp <- temp %>% mutate(fix = fix - 1)
  write_csv(temp, paste0(name_list[i],"_All_Used.csv"))
}

# Create csv of available points for 'foraging' (to be parsed later)
setwd("~/Box Sync/Dissertation/Behavioral_SSF/Rerun_Extraction_Results")
for (i in 1:length(name_list)) {
  if (i < 6) {
    state.temp <- read_csv(paste0(name_list[i],"_Foraging_Used.csv")) %>%
      dplyr::select(date, ID, states) 
    temp <- read.csv(paste0(name_list[i],"_Foraging_Available.csv")) %>%
      st_as_sf(., coords = 4:5, crs = "+init=epsg:32733") %>%
      mutate(Green = green_avail,
             Wet = wet_avail,
             Road_Dens = roads_avail) %>%
      bind_cols(., state.temp) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 0) %>% 
      mutate(Green_Norm = (Green - means09['Mean_Greenness_2009']) / sd09['Mean_Greenness_2009'],
             Wet_Norm = (Wet - means09['Mean_Wetness_2009']) / sd09['Mean_Wetness_2009'],
             Road_Dens_Norm = (Road_Dens - means09['Road_Density']) / sd09['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  } else {
    state.temp <- read_csv(paste0(name_list[i],"_Foraging_Used.csv")) %>%
      dplyr::select(date, ID, states) 
    temp <- read_csv(paste0(name_list[i],"_Foraging_Available.csv")) %>%
      st_as_sf(., coords = 4:5, crs = "+init=epsg:32733") %>%
      mutate(Green = green_avail,
             Wet = wet_avail,
             Road_Dens = roads_avail) %>%
      bind_cols(., state.temp) %>%
      mutate(fix =  seq(1,nrow(.),1), binom = 0) %>% 
      mutate(Green_Norm = (Green - means10['Mean_Greenness_2010']) / sd10['Mean_Greenness_2010'],
             Wet_Norm = (Wet - means10['Mean_Wetness_2010']) / sd10['Mean_Wetness_2010'],
             Road_Dens_Norm = (Road_Dens - means10['Road_Density']) / sd10['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  }
  write_csv(temp, paste0(name_list[i],"_Foraging_Available.csv"))
}

# Combine used and available, match available at time t with used at time t+1, only foraging points
for (k in 1:length(name_list)) {
  used <- read_csv(paste0(name_list[k],"_Foraging_Used.csv"))
  avail <- read_csv(paste0(name_list[k],"_Foraging_Available.csv")) %>% 
    dplyr::filter(states == 2)
  
  date_used <- as.numeric(used$date)
  date_avail <- as.numeric(avail$date)
  
  new_used <- c()
  new_avail <- c()
  w=1
  for (i in 1:nrow(avail)) {
    if ((date_avail[i] + 1200) %in% date_used) {
      which(date_used == (date_avail[i] + 1200))
      new_used[w] <- which(date_used == (date_avail[i] + 1200))
      new_avail[w] <- i
      w = w + 1
    }
  }
  
  used_temp <- data.frame(used[new_used[1],])
  avail_temp <- data.frame(avail[new_avail[1],])
  for (j in 2:length(new_used)) {
    used_temp[j,] <- used[new_used[j],]
    avail_temp[j,] <- avail[new_avail[j],]
  }
  
  used_temp <- mutate(used_temp, fix = seq(1,nrow(used_temp),1))
  avail_temp <- mutate(avail_temp, fix = seq(1,nrow(avail_temp),1))
  all_pts <- rbind(used_temp, avail_temp)
  write_csv(all_pts, paste0(name_list[k],"_Foraging_Final.csv"))
}

# Create csv of available points for 'all'
setwd("~/Box Sync/Dissertation/Behavioral_SSF/Rerun_Extraction_Results")
for (i in 1:length(name_list)) {
  if (i < 6) {
    state.temp <- read_csv(paste0(name_list[i],"_All_Used.csv")) %>%
      dplyr::select(date, ID, states) 
    temp <- read.csv(paste0(name_list[i],"_All_Available.csv")) %>%
      st_as_sf(., coords = 4:5, crs = "+init=epsg:32733") %>%
      mutate(Green = green_avail,
             Wet = wet_avail,
             Road_Dens = roads_avail) %>%
      bind_cols(., state.temp) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 0) %>% 
      mutate(Green_Norm = (Green - means09['Mean_Greenness_2009']) / sd09['Mean_Greenness_2009'],
             Wet_Norm = (Wet - means09['Mean_Wetness_2009']) / sd09['Mean_Wetness_2009'],
             Road_Dens_Norm = (Road_Dens - means09['Road_Density']) / sd09['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  } else {
    state.temp <- read_csv(paste0(name_list[i],"_All_Used.csv")) %>%
      dplyr::select(date, ID, states) 
    temp <- read_csv(paste0(name_list[i],"_All_Available.csv")) %>%
      st_as_sf(., coords = 4:5, crs = "+init=epsg:32733") %>%
      mutate(Green = green_avail,
             Wet = wet_avail,
             Road_Dens = roads_avail) %>%
      bind_cols(., state.temp) %>%
      mutate(fix =  seq(1,nrow(.),1), binom = 0) %>% 
      mutate(Green_Norm = (Green - means10['Mean_Greenness_2010']) / sd10['Mean_Greenness_2010'],
             Wet_Norm = (Wet - means10['Mean_Wetness_2010']) / sd10['Mean_Wetness_2010'],
             Road_Dens_Norm = (Road_Dens - means10['Road_Density']) / sd10['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  }
  write_csv(temp, paste0(name_list[i],"_All_Available.csv"))
}

# Combine used and available, match available at time t with used at time t+1
for (k in 1:length(name_list)) {
  used <- read_csv(paste0(name_list[k],"_All_Used.csv"))
  avail <- read_csv(paste0(name_list[k],"_All_Available.csv"))

  date_used <- as.numeric(used$date)
  date_avail <- as.numeric(avail$date)
  
  new_used <- c()
  new_avail <- c()
  w=1
  for (i in 1:nrow(avail)) {
    if ((date_avail[i] + 1200) %in% date_used) {
      which(date_used == (date_avail[i] + 1200))
      new_used[w] <- which(date_used == (date_avail[i] + 1200))
      new_avail[w] <- i
      w = w + 1
    }
  }
  
  used_temp <- data.frame(used[new_used[1],])
  avail_temp <- data.frame(avail[new_avail[1],])
  for (j in 2:length(new_used)) {
    used_temp[j,] <- used[new_used[j],]
    avail_temp[j,] <- avail[new_avail[j],]
  }
  
  used_temp <- mutate(used_temp, fix = seq(1,nrow(used_temp),1))
  avail_temp <- mutate(avail_temp, fix = seq(1,nrow(avail_temp),1))
  all_pts <- rbind(used_temp, avail_temp)
  write_csv(all_pts, paste0(name_list[k],"_All_Final.csv"))
}
