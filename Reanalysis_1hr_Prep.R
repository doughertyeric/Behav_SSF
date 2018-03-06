library(tidyverse)

##############################################

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

setwd("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/")
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

zebra09 <- read_csv("Zebra_Anthrax_2009_Cleaned.csv")[,-1]
zebra09 <- cbind(zebra09, states09)
zebra10 <- read_csv("Zebra_Anthrax_2010_Cleaned.csv")[,-1]
zebra10 <- cbind(zebra10, states10)
all_zebra <- rbind(zebra09, zebra10)

# 2009 zebra with states variable and means, then reduced to hourly fixes
temp <- all_zebra %>% filter(ID == name_list[1])
temp <- temp[1:((nrow(temp) %/% 3) * 3), ]
temp <- temp %>% mutate(hour = rep(1:(nrow(.)/3), each=3)) %>%
  mutate(mean = 0)
for (j in 1:(nrow(temp)-2)) {
  temp$mean[j] <- round(mean(c(temp$states[j], temp$states[(j+1)], temp$states[(j+2)]), na.rm=TRUE))
}
hour_09 <- temp[seq(1,nrow(temp), by=3), ]

for (i in 2:5) {
  temp <- all_zebra %>% filter(ID == name_list[i])
  temp <- temp[1:((nrow(temp) %/% 3) * 3), ]
  temp <- temp %>% mutate(hour = rep(1:(nrow(.)/3), each=3)) %>%
    mutate(mean = 0)
  
  for (j in 1:(nrow(temp)-2)) {
    temp$mean[j] <- round(mean(c(temp$states[j], temp$states[(j+1)], temp$states[(j+2)]), na.rm=TRUE))
  }
  
  temp <- temp[seq(1,nrow(temp), by=3), ]
  hour_09 <- rbind(hour_09, temp)
}

# 2010 zebra with states variable and means, then reduced to hourly fixes
temp <- all_zebra %>% filter(ID == name_list[6])
temp <- temp[1:((nrow(temp) %/% 3) * 3), ]
temp <- temp %>% mutate(hour = rep(1:(nrow(.)/3), each=3)) %>%
  mutate(mean = 0)
for (j in 1:(nrow(temp)-2)) {
  temp$mean[j] <- round(mean(c(temp$states[j], temp$states[(j+1)], temp$states[(j+2)]), na.rm=TRUE))
}
hour_10 <- temp[seq(1,nrow(temp), by=3), ]

for (i in 7:11) {
  temp <- all_zebra %>% filter(ID == name_list[i])
  temp <- temp[1:((nrow(temp) %/% 3) * 3), ]
  temp <- temp %>% mutate(hour = rep(1:(nrow(.)/3), each=3)) %>%
    mutate(mean = 0)
  
  for (j in 1:(nrow(temp)-2)) {
    temp$mean[j] <- round(mean(c(temp$states[j], temp$states[(j+1)], temp$states[(j+2)]), na.rm=TRUE))
  }
  
  temp <- temp[seq(1,nrow(temp), by=3), ]
  hour_10 <- rbind(hour_10, temp)
}

all_hours <- rbind(hour_09, hour_10)

##############################################

setwd("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Hour_SSF")
library(adehabitatLT)

hour_ltraj <- as.ltraj(xy = all_hours[,c('x','y')], date = all_hours$date,
                       id = all_hours$ID, proj4string = CRS('+init=epsg:32733'))

temp <- hour_ltraj[[1]] %>%
  mutate(ID = name_list[1])
all_hour <- temp
for (i in 2:4) {
  temp <- hour_ltraj[[i]] %>%
    mutate(ID = name_list[i])
  all_hour <- rbind(all_hour, temp)
}
temp <- hour_ltraj[[6]] %>%
  mutate(ID = name_list[5])
all_hour <- rbind(all_hour, temp)
all_hour <- cbind(all_hour, states = hour_09$states, hour = hour_09$hour, mean = hour_09$mean)
write.csv(all_hour, 'Zebra_2009_1hr.csv')

temp <- hour_ltraj[[5]] %>%
  mutate(ID = name_list[6])
all_hour <- temp
for (i in 7:11) {
  temp <- hour_ltraj[[i]] %>%
    mutate(ID = name_list[i])
  all_hour <- rbind(all_hour, temp)
}
all_hour <- cbind(all_hour, states = hour_10$states, hour = hour_10$hour, mean = hour_10$mean)
write.csv(all_hour, 'Zebra_2010_1hr.csv')

##############################################

library(fitdistrplus)

setwd("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Hour_SSF")
zeb_09_1hr <- read.csv('Zebra_2009_1hr.csv')[,-1]
zeb_10_1hr <- read.csv('Zebra_2010_1hr.csv')[,-1]
all_hour <- rbind(zeb_09_1hr, zeb_10_1hr)

gamma_params <- data.frame(matrix(0,11,3))
for (i in 1:length(name_list)) {
  temp <- all_hour %>% filter(ID == name_list[i]) %>%
    filter(!is.na(dist))
  step_lengths <- temp$dist
  fit_gamma <- fitdist(step_lengths, "gamma", method='mme')
  ests.gamma <- bootdist(fit_gamma, niter = 1000)
  gamma_params[i,1] <- i
  gamma_params[i,2] <- median(summary(ests.gamma)[1]$estim[,1])
  gamma_params[i,3] <- median(summary(ests.gamma)[1]$estim[,2])
}
colnames(gamma_params) <- c('ID', 'gamma.shape', 'gamma.rate')
write.csv(gamma_params, "Hour_Zebra_Step_Parameters.csv")

################################################

for (i in 1:length(name_list)) {
  temp <- all_hour %>% filter(ID == name_list[i]) %>%
    filter(!is.na(x)) %>%
    st_as_sf(coords=1:2) %>%
    st_set_crs('+init=epsg:32733')
  write_sf(temp, paste0(name_list[i],'_1hr.shp'))
}

################################################

foraging <- all_hour %>% filter(mean == 2)
gamma_params <- data.frame(matrix(0,11,3))
for (i in 1:length(name_list)) {
  temp <- foraging %>% filter(ID == name_list[i]) %>%
    filter(!is.na(dist))
  step_lengths <- temp$dist
  fit_gamma <- fitdist(step_lengths, "gamma", method='mme')
  ests.gamma <- bootdist(fit_gamma, niter = 1000)
  gamma_params[i,1] <- i
  gamma_params[i,2] <- median(summary(ests.gamma)[1]$estim[,1])
  gamma_params[i,3] <- median(summary(ests.gamma)[1]$estim[,2])
}
colnames(gamma_params) <- c('ID', 'gamma.shape', 'gamma.rate')
write.csv(gamma_params, "Hour_Foraging_Step_Parameters.csv")

###########################################

library(velox)

green09 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/ENP_Predictors/Mean_Greenness_2009.tif")
v_green09 <- velox(green09)
wet09 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/ENP_Predictors/Mean_Wetness_2009.tif")
v_wet09 <- velox(wet09)
road09 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/ENP_Predictors/Road_Density.tif")
v_road09 <- velox(road09)

green10 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/ENP_Predictors/Mean_Greenness_2010.tif")
v_green10 <- velox(green10)
wet10 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/ENP_Predictors/Mean_Wetness_2010.tif")
v_wet10 <- velox(wet10)
road10 <- raster("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/ENP_Predictors/Road_Density.tif")
v_road10 <- velox(road10)

for (i in 1:length(name_list)) {
  temp <- st_read(paste0('Hour_SSF/',name_list[i],'_1hr.shp'))
  if (i < 6) {
    used_green <- v_green09$extract_points(temp)
    used_wet <- v_wet09$extract_points(temp)
    used_road <- v_road09$extract_points(temp)
  } else {
    used_green <- v_green10$extract_points(temp)
    used_wet <- v_wet10$extract_points(temp)
    used_road <- v_road10$extract_points(temp) 
  }
  temp.sp <- temp %>% as('Spatial')
  dates <- data.frame(as.POSIXct(as.character(temp$date), format='%Y-%m-%d %H:%M:%S'))
  used_df <- data.frame(cbind(used_green, used_wet, used_road, temp.sp@coords, temp.sp@data$mean, 
                              date = dates[,1]))
  used_df$date <- as.POSIXct(used_df$date, format='%Y-%m-%d %H:%M:%S', origin = '1970-01-01 00:00:00')
  colnames(used_df) <- c('used_green', 'used_wet', 'used_road', 'x', 'y', 'mean', 'date')
  write.csv(used_df, paste0('Hour_SSF/',name_list[i],"_All_Used_1hr.csv"))
  
}

######################################################################
###### Use Python Scripts to obtain Available (All and Forage) #######
######################################################################

road_dens <- raster('ENP_Predictors/Road_Density.tif')
Green_09 <- raster('ENP_Predictors/Mean_Greenness_2009.tif')
Green_10 <- raster('ENP_Predictors/Mean_Greenness_2010.tif')
Wet_09 <- raster('ENP_Predictors/Mean_Wetness_2009.tif')
Wet_10 <- raster('ENP_Predictors/Mean_Wetness_2010.tif')

pred_stack_09 <- stack(Green_09, Wet_09, road_dens)
pred_stack_10 <- stack(Green_10, Wet_10, road_dens)

means09 <- cellStats(pred_stack_09, stat='mean', na.rm=TRUE)
sd09 <- cellStats(pred_stack_09, stat='sd', na.rm=TRUE)
means10 <- cellStats(pred_stack_10, stat='mean', na.rm=TRUE)
sd10 <- cellStats(pred_stack_10, stat='sd', na.rm=TRUE)

########################################################
############## FORMAT AND COMPILE `ALL` ################
########################################################

# Recompile ALL_USED points and normalize variables

for (i in 1:length(name_list)) {
  if (i < 6) {
    temp <- read_csv(paste0('Hour_SSF/',name_list[i],"_All_Used_1hr.csv")) %>% 
      mutate(Green = used_green,
             Wet = used_wet,
             Road_Dens = used_road,
             states = mean) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 1, ID = name_list[i])
    temp <- temp %>% mutate(Green_Norm = (Green - means09['Mean_Greenness_2009']) / sd09['Mean_Greenness_2009'],
                            Wet_Norm = (Wet - means09['Mean_Wetness_2009']) / sd09['Mean_Wetness_2009'],
                            Road_Dens_Norm = (Road_Dens - means09['Road_Density']) / sd09['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  } else {
    temp <- read_csv(paste0('Hour_SSF/',name_list[i],"_All_Used_1hr.csv")) %>%
      mutate(Green = used_green,
             Wet = used_wet,
             Road_Dens = used_road,
             states = mean) %>%
      mutate(fix =  seq(1,nrow(.),1), binom = 1, ID = name_list[i])
    temp <- temp %>% mutate(Green_Norm = (Green - means10['Mean_Greenness_2010']) / sd10['Mean_Greenness_2010'],
                            Wet_Norm = (Wet - means10['Mean_Wetness_2010']) / sd10['Mean_Wetness_2010'],
                            Road_Dens_Norm = (Road_Dens - means10['Road_Density']) / sd10['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  }
  used_df <- temp %>% mutate(fix = fix - 1)
  write_csv(used_df, paste0('Hour_SSF/Final/',name_list[i],"_All_Used_1hr.csv"))
}

# Recompile ALL_AVAILABLE points and normalize variables

for (i in 1:length(name_list)) {
  if (i < 6) {
    temp <- st_read(paste0('Hour_SSF/',name_list[i],'_1hr.shp'))
    temp.sp <- temp %>% as('Spatial')
    avail_df <- read_csv(paste0('Hour_SSF/',name_list[i],"_All_Available_1hr.csv")) %>% 
      mutate(Green = green_avail,
             Wet = wet_avail,
             Road_Dens = roads_avail) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 0, ID = name_list[i])
    dates <- data.frame(as.POSIXct(as.character(temp$date), format='%Y-%m-%d %H:%M:%S'))
    avail_df <- cbind(avail_df, mean=temp.sp@data$mean, date=dates[,1])
    avail_df$date <- as.POSIXct(avail_df$date, format='%Y-%m-%d %H:%M:%S', origin = '1970-01-01 00:00:00')
    avail_df <- avail_df %>% mutate(Green_Norm = (Green - means09['Mean_Greenness_2009']) / sd09['Mean_Greenness_2009'],
                            Wet_Norm = (Wet - means09['Mean_Wetness_2009']) / sd09['Mean_Wetness_2009'],
                            Road_Dens_Norm = (Road_Dens - means09['Road_Density']) / sd09['Road_Density'],
                            states = mean) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  } else {
    temp <- st_read(paste0('Hour_SSF/',name_list[i],'_1hr.shp'))
    temp.sp <- temp %>% as('Spatial')
    avail_df <- read_csv(paste0('Hour_SSF/',name_list[i],"_All_Available_1hr.csv")) %>% 
      mutate(Green = green_avail,
             Wet = wet_avail,
             Road_Dens = roads_avail) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 0, ID = name_list[i])
    dates <- data.frame(as.POSIXct(as.character(temp$date), format='%Y-%m-%d %H:%M:%S'))
    avail_df <- cbind(avail_df, mean=temp.sp@data$mean, date=dates[,1])
    avail_df$date <- as.POSIXct(avail_df$date, format='%Y-%m-%d %H:%M:%S', origin = '1970-01-01 00:00:00')
    avail_df <- avail_df %>% mutate(Green_Norm = (Green - means10['Mean_Greenness_2010']) / sd10['Mean_Greenness_2010'],
                            Wet_Norm = (Wet - means10['Mean_Wetness_2010']) / sd10['Mean_Wetness_2010'],
                            Road_Dens_Norm = (Road_Dens - means10['Road_Density']) / sd10['Road_Density'],
                            states = mean) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  }
  write_csv(avail_df, paste0('Hour_SSF/Final/',name_list[i],"_All_Available_1hr.csv"))
}

########################################

# Combine used and available, match available at time t with used at time t+1

for (k in 1:length(name_list)) {
  used <- read_csv(paste0('Hour_SSF/Final/',name_list[k],"_All_Used_1hr.csv")) %>%
    filter(!is.na(date))
  avail <- read_csv(paste0('Hour_SSF/Final/',name_list[k],"_All_Available_1hr.csv")) %>%
    filter(!is.na(date))
  
  date_used <- as.numeric(used$date)
  date_avail <- as.numeric(avail$date)
  
  new_used <- c()
  new_avail <- c()
  w=1
  for (i in 1:nrow(avail)) {
    if ((date_avail[i] + 3600) %in% date_used) {
      which(date_used == (date_avail[i] + 3600))
      new_used[w] <- which(date_used == (date_avail[i] + 3600))
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
  write_csv(all_pts, paste0('Hour_SSF/Final/',name_list[k],"_All_Final.csv"))
}


########################################################
########### FORMAT AND COMPILE `FORAGING` ##############
########################################################

# Recompile FORAGING_USED points and normalize variables
# Begin with All_Used file and extract points in state 2

for (i in 1:length(name_list)) {
  if (i < 6) {
    temp <- read_csv(paste0('Hour_SSF/',name_list[i],"_All_Used_1hr.csv")) %>% 
      mutate(Green = used_green,
             Wet = used_wet,
             Road_Dens = used_road,
             states = mean) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 1, ID = name_list[i]) %>%
      filter(states == 2)
    temp <- temp %>% mutate(Green_Norm = (Green - means09['Mean_Greenness_2009']) / sd09['Mean_Greenness_2009'],
                            Wet_Norm = (Wet - means09['Mean_Wetness_2009']) / sd09['Mean_Wetness_2009'],
                            Road_Dens_Norm = (Road_Dens - means09['Road_Density']) / sd09['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  } else {
    temp <- read_csv(paste0('Hour_SSF/',name_list[i],"_All_Used_1hr.csv")) %>%
      mutate(Green = used_green,
             Wet = used_wet,
             Road_Dens = used_road,
             states = mean) %>%
      mutate(fix =  seq(1,nrow(.),1), binom = 1, ID = name_list[i]) %>%
      filter(states == 2)
    temp <- temp %>% mutate(Green_Norm = (Green - means10['Mean_Greenness_2010']) / sd10['Mean_Greenness_2010'],
                            Wet_Norm = (Wet - means10['Mean_Wetness_2010']) / sd10['Mean_Wetness_2010'],
                            Road_Dens_Norm = (Road_Dens - means10['Road_Density']) / sd10['Road_Density']) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  }
  used_df <- temp %>% mutate(fix = fix - 1)
  write_csv(used_df, paste0('Hour_SSF/Final/',name_list[i],"_Foraging_Used_1hr.csv"))
}

# Recompile FORAGING_AVAILABLE points and normalize variables

for (i in 1:length(name_list)) {
  if (i < 6) {
    temp <- st_read(paste0('Hour_SSF/',name_list[i],'_1hr.shp'))
    temp.sp <- temp %>% as('Spatial')
    avail_df <- read_csv(paste0('Hour_SSF/',name_list[i],"_Foraging_Available_1hr.csv")) %>% 
      mutate(Green = green_avail,
             Wet = wet_avail,
             Road_Dens = roads_avail) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 0, ID = name_list[i])
    dates <- data.frame(as.POSIXct(as.character(temp$date), format='%Y-%m-%d %H:%M:%S'))
    avail_df <- cbind(avail_df, mean=temp.sp@data$mean, date=dates[,1])
    avail_df$date <- as.POSIXct(avail_df$date, format='%Y-%m-%d %H:%M:%S', origin = '1970-01-01 00:00:00')
    avail_df <- avail_df %>% mutate(Green_Norm = (Green - means09['Mean_Greenness_2009']) / sd09['Mean_Greenness_2009'],
                                    Wet_Norm = (Wet - means09['Mean_Wetness_2009']) / sd09['Mean_Wetness_2009'],
                                    Road_Dens_Norm = (Road_Dens - means09['Road_Density']) / sd09['Road_Density'],
                                    states = mean) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  } else {
    temp <- st_read(paste0('Hour_SSF/',name_list[i],'_1hr.shp'))
    temp.sp <- temp %>% as('Spatial')
    avail_df <- read_csv(paste0('Hour_SSF/',name_list[i],"_Foraging_Available_1hr.csv")) %>% 
      mutate(Green = green_avail,
             Wet = wet_avail,
             Road_Dens = roads_avail) %>%
      mutate(fix = seq(1,nrow(.),1), binom = 0, ID = name_list[i])
    dates <- data.frame(as.POSIXct(as.character(temp$date), format='%Y-%m-%d %H:%M:%S'))
    avail_df <- cbind(avail_df, mean=temp.sp@data$mean, date=dates[,1])
    avail_df$date <- as.POSIXct(avail_df$date, format='%Y-%m-%d %H:%M:%S', origin = '1970-01-01 00:00:00')
    avail_df <- avail_df %>% mutate(Green_Norm = (Green - means10['Mean_Greenness_2010']) / sd10['Mean_Greenness_2010'],
                                    Wet_Norm = (Wet - means10['Mean_Wetness_2010']) / sd10['Mean_Wetness_2010'],
                                    Road_Dens_Norm = (Road_Dens - means10['Road_Density']) / sd10['Road_Density'],
                                    states = mean) %>%
      dplyr::select(date, ID, Green, Wet, Road_Dens, Green_Norm, Wet_Norm, Road_Dens_Norm, fix, states, binom)
  }
  write_csv(avail_df, paste0('Hour_SSF/Final/',name_list[i],"_Foraging_Available_1hr.csv"))
}

########################################

# Combine used and available, match available at time t with used at time t+1

for (k in 1:length(name_list)) {
  used <- read_csv(paste0('Hour_SSF/Final/',name_list[k],"_Foraging_Used_1hr.csv")) %>%
    filter(!is.na(date))
  avail <- read_csv(paste0('Hour_SSF/Final/',name_list[k],"_Foraging_Available_1hr.csv")) %>%
    filter(!is.na(date))
  
  date_used <- as.numeric(used$date)
  date_avail <- as.numeric(avail$date)
  
  new_used <- c()
  new_avail <- c()
  w=1
  for (i in 1:nrow(avail)) {
    if ((date_avail[i] + 3600) %in% date_used) {
      which(date_used == (date_avail[i] + 3600))
      new_used[w] <- which(date_used == (date_avail[i] + 3600))
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
  write_csv(all_pts, paste0('Hour_SSF/Final/',name_list[k],"_Foraging_Final.csv"))
}

########################################################
############ Conduct SSF using 1 hr FINAL ##############
########################################################

norm_stack_09 <- stack((Green_09 - means09[1]) / sd09[1],
                       (Wet_09 - means09[2]) / sd09[2],
                       (road_dens - means09[3])/ sd09[3])

norm_stack_10 <- stack((Green_10 - means10[1]) / sd10[1],
                       (Wet_10 - means10[2]) / sd10[2],
                       (road_dens - means10[3])/ sd10[3])

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Hour_SSF/Final")

cond.logit.all <- function(i) {
  all_final <- read_csv(paste0(name_list[i],"_All_Final.csv")) %>%
    dplyr::select("Green", "Wet", "Road_Dens",
                  "Green_Norm", "Wet_Norm", "Road_Dens_Norm", 
                  "ID", "binom", "fix")
  assign(paste0(name_list[i],"_All"), all_final)
  assign(paste0(name_list[i],"_Corr"), cor(all_final[,1:4]))
  assign(paste0(name_list[i],"_Model"), clogistic(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm,
                                                  data=all_final, strata=fix))
  assign(paste0(name_list[i],"_clogit"), clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm +
                                                  strata(fix),
                                                data=all_final))
  assign(paste0(name_list[i],"_GLM"), glm(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm,
                                          data=all_final, family=binomial(link='logit')))
  
  out_list <- list(get(paste0(name_list[i],"_All")),
                   get(paste0(name_list[i],"_Corr")),
                   get(paste0(name_list[i],"_Model")),
                   get(paste0(name_list[i],"_clogit")),
                   get(paste0(name_list[i],"_GLM")))
  return(out_list)
}

cond.logit.forage <- function(i) {
  forage_final <- read_csv(paste0(name_list[i],"_Foraging_Final.csv")) %>%
    dplyr::select("Green", "Wet", "Road_Dens",
                  "Green_Norm", "Wet_Norm", "Road_Dens_Norm", 
                  "ID", "binom", "fix")
  assign(paste0(name_list[i],"_All"), forage_final)
  assign(paste0(name_list[i],"_Corr"), cor(forage_final[,1:4]))
  assign(paste0(name_list[i],"_Model"), clogistic(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm,
                                                  data=forage_final, strata=fix))
  assign(paste0(name_list[i],"_clogit"), clogit(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm + 
                                                  strata(fix), 
                                                data=forage_final))
  assign(paste0(name_list[i],"_GLM"), glm(binom ~ Wet_Norm + Green_Norm + Road_Dens_Norm,
                                          data=forage_final, family=binomial(link='logit')))
  
  out_list <- list(get(paste0(name_list[i],"_All")),
                   get(paste0(name_list[i],"_Corr")),
                   get(paste0(name_list[i],"_Model")),
                   get(paste0(name_list[i],"_clogit")),
                   get(paste0(name_list[i],"_GLM")))
  return(out_list)
}

############## Individual Predictions ###################

library(survival)
library(Epi)

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

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Hour_SSF/Results/")
for (i in 1:length(name_list)) {
  writeRaster(get(paste0(name_list[i],"_All_Risk")), paste0(name_list[i],"_All_Risk.tif"), format='GTiff')
  writeRaster(get(paste0(name_list[i],"_Forage_Risk")), paste0(name_list[i],"_Forage_Risk.tif"), format='GTiff')
}

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

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Hour_SSF/Final/")
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

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Hour_SSF/Results/")
writeRaster(All_2009_Risk, "Zebra_2009_All_Risk.tif", format='GTiff')
writeRaster(All_2010_Risk, "Zebra_2010_All_Risk.tif", format='GTiff')
writeRaster(Forage_2009_Risk, "Zebra_2009_Forage_Risk.tif", format='GTiff')
writeRaster(Forage_2010_Risk, "Zebra_2010_Forage_Risk.tif", format='GTiff')

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Figures")
pdf('Zebra_2009_All_Risk_1hr.pdf', height=6, width=8)
plot(All_2009_Risk)
dev.off()
pdf('Zebra_2010_All_Risk_1hr.pdf', height=6, width=8)
plot(All_2010_Risk)
dev.off()
pdf('Zebra_2009_Foraging_Risk_1hr.pdf', height=6, width=8)
plot(Forage_2009_Risk)
dev.off()
pdf('Zebra_2010_Foraging_Risk_1hr.pdf', height=6, width=8)
plot(Forage_2010_Risk)
dev.off()

################################################################
##################### Overlap Analysis #########################
################################################################


