library(tidyverse)

setwd("~/Box Sync/Dissertation/Behavioral_SSF/Original_Extraction_Results/")
# Files are now in 'Original_Extraction_Results' folder

########################################################

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

for (i in 1:length(name_list)) {
  used <- read_csv(paste0(name_list[i],"_All_Used.csv")) %>%
    mutate(ID = name_list[i], binom = 1, fix = seq(1,nrow(.),1))
  used_veg <- read_csv(paste0(name_list[i],"_All_Veg_Used.csv")) %>%
    mutate(ID = name_list[i], binom = 1, fix = seq(1,nrow(.),1))
  used_roads <- read_csv(paste0(name_list[i],"_All_Roads_Used.csv")) %>%
    mutate(ID = name_list[i], binom = 1, fix = seq(1,nrow(.),1))
  
  colnames(used) <- c("NDVI", "Green", "Dist_Roads", "Dist_Water", "Wet", "x", "y", "ID", "binom", "fix")
  colnames(used_roads) <- c("Road_Dens", "x", "y")
  all_used <- cbind(used[,c("x", "y", "NDVI", "Green", "Wet")], used_roads[,c('Road_Dens')],
                    used_veg[,c("prop_bare", "prop_steppe", "prop_grass", "prop_shrub", "prop_low", "prop_high")],
                    used[,c("binom", "fix", "ID")])
  
  avail <- read_csv(paste0(name_list[i],"_All_Available.csv")) %>%
    mutate(ID =name_list[i], binom = 0, fix = seq(1,nrow(.),1))
  avail_veg <- read_csv(paste0(name_list[i],"_All_Veg_Available.csv")) %>%
    mutate(ID = name_list[i], binom = 0, fix = seq(1,nrow(.),1))
  avail_roads <- read_csv(paste0(name_list[i],"_All_Roads_Available.csv")) %>%
    mutate(ID = name_list[i], binom = 1, fix = seq(1,nrow(.),1))
  
  colnames(avail) <- c("NDVI", "Green", "Dist_Roads", "Dist_Water", "Wet", "x", "y", "ID", "binom", "fix")
  colnames(avail_roads) <- c("Road_Dens", "x", "y")
  all_avail <- cbind(avail[,c("x", "y", "NDVI", "Green", "Wet")], avail_roads[,c('Road_Dens')],
                    avail_veg[,c("prop_bare", "prop_steppe", "prop_grass", "prop_shrub", "prop_low", "prop_high")],
                    avail[,c("binom", "fix", "ID")])
  
  all_pts <- rbind(all_used, all_avail)
  write.csv(all_pts, paste0(name_list[i], "_All_Final.csv"))
}

for (i in 1:length(name_list)) {
  used <- read_csv(paste0(name_list[i],"_Foraging_Used.csv")) %>%
    mutate(ID = name_list[i], binom = 1, fix = seq(1,nrow(.),1))
  used_veg <- read_csv(paste0(name_list[i],"_Foraging_Veg_Used.csv")) %>%
    mutate(ID = name_list[i], binom = 1, fix = seq(1,nrow(.),1))
  used_roads <- read_csv(paste0(name_list[i],"_Foraging_Roads_Used.csv")) %>%
    mutate(ID = name_list[i], binom = 1, fix = seq(1,nrow(.),1))
  
  colnames(used) <- c("NDVI", "Green", "Dist_Roads", "Dist_Water", "Wet", "x", "y", "ID", "binom", "fix")
  colnames(used_roads) <- c("Road_Dens", "x", "y")
  all_used <- cbind(used[,c("x", "y", "NDVI", "Green", "Wet")], used_roads[,c("Road_Dens")],
                    used_veg[,c("prop_bare", "prop_steppe", "prop_grass", "prop_shrub", "prop_low", "prop_high")],
                    used[,c("binom", "fix", "ID")])
  
  avail <- read_csv(paste0(name_list[i],"_Foraging_Available.csv")) %>%
    mutate(ID =name_list[i], binom = 0, fix = seq(1,nrow(.),1))
  avail_veg <- read_csv(paste0(name_list[i],"_Foraging_Veg_Available.csv")) %>%
    mutate(ID = name_list[i], binom = 0, fix = seq(1,nrow(.),1))
  avail_roads <- read_csv(paste0(name_list[i],"_Foraging_Roads_Available.csv")) %>%
    mutate(ID = name_list[i], binom = 0, fix = seq(1,nrow(.),1))
  
  colnames(avail) <- c("NDVI", "Green", "Dist_Roads", "Dist_Water", "Wet", "x", "y", "ID", "binom", "fix")
  colnames(avail_roads) <- c("Road_Dens", "x", "y")
  all_avail <- cbind(avail[,c("x", "y", "NDVI", "Green", "Wet")], avail_roads[,c("Road_Dens")],
                     avail_veg[,c("prop_bare", "prop_steppe", "prop_grass", "prop_shrub", "prop_low", "prop_high")],
                     avail[,c("binom", "fix", "ID")])
  
  all_pts <- rbind(all_used, all_avail)
  write.csv(all_pts, paste0(name_list[i], "_Foraging_Final.csv"))
}
