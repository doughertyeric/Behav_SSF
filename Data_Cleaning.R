library(tidyverse)
library(adehabitatLT)
library(lubridate)

setwd("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/")
name.list <- c("AG059", "AG061", "AG062", "AG063", "AG068", "AG252", "AG253", "AG255", "AG256")

setwd("/Users/ericdougherty/Box Sync/Dissertation/Behavioral_SSF/Uncleaned_Paths")
name.list <- c("AG059", "AG061", "AG062", "AG063", "AG068", "AG252", "AG253", "AG255", "AG256")

# Load in all zebra and compile a single dataset for ltraj object
i=1
data <- read.csv(paste(as.character(name.list[i]),"_Full.csv",sep=""), header=TRUE)
colnames(data) <- c("X", "Name", "Date", "Longitude", "Latitude", "Speed", "Direction", "Temp", "Altitude", "PDOP")
data$Date <- as.POSIXct(strptime(as.character(data$Date), tz='GMT', format='%m/%d/%y %H:%M:%S'))
data <- data[!duplicated(data$Date), ]
data$ID <- paste(name.list[i])

for (p in 2:length(name.list)) {
  data2 <- read.csv(paste(as.character(name.list[p]),"_Full.csv",sep=""))
  colnames(data2) <- c("X", "Name", "Date", "Longitude", "Latitude", "Speed", "Direction", "Temp",   "Altitude", "PDOP")
  data2$Date <- as.POSIXct(strptime(as.character(data2$Date), tz='GMT', format='%m/%d/%y %H:%M:%S'))
  data2 <- data2[!duplicated(data2$Date), ]
  data2$ID <- paste(name.list[p])
  data <- rbind(data, data2)
}

data <- data[!is.na(data$Longitude), ]
data <- data[data$Longitude < 17.2, ]
data <- data[data$Longitude > 14.4, ]
data <- data[data$Latitude > -19.5, ]
data <- data[data$Latitude < -18.5, ]

# Filter out anthrax season points from 2009 and 2010, new column with season as ID
# Apparent collar/satellite malfunction in early April 2009
anthrax.season.09 <- interval(ymd("2009-04-25"), ymd("2009-07-01")) 
anth.09 <- data %>% group_by(ID) %>% filter(Date %within% anthrax.season.09) %>% mutate(Season=paste0(ID,"_2009"))
#write.csv(anth.09, "Zebra_Anthrax_2009.csv")
anthrax.season.10 <- interval(ymd("2010-02-01"), ymd("2010-07-01"))
anth.10 <- data %>% group_by(ID) %>% filter(Date %within% anthrax.season.10) %>% mutate(Season=paste0(ID,"_2010"))
#write.csv(anth.10, "Zebra_Anthrax_2010.csv")
all.anth <- data.frame(rbind(anth.09, anth.10))

# Create ltraj object out of anthrax season points
latlong.pts <- SpatialPoints(all.anth[,c("Longitude","Latitude")], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
utm.pts <- spTransform(latlong.pts, CRS("+proj=utm +south +zone=33 +ellps=WGS84"))
anth.traj <- as.ltraj(xy = utm.pts@coords, date = all.anth$Date, id = all.anth$Season)
anth.traj <- anth.traj[summary(anth.traj)$nb.reloc > 10000]

refda <- strptime("00:00", "%H:%M", tz="GMT")
anth.NA <- setNA(anth.traj, refda, 20, units = "min")
traj.round <- sett0(anth.NA, refda, 20, units = "min")

# # Estimate gamma distribution parameters for each season/ID
# gamma.vars <- data.frame(matrix(0,length(traj.round),3))
# colnames(gamma.vars) <- c("ID", "shape", "rate")
# for (i in 1:length(traj.round)) {
#   temp <- traj.round[[i]]
#   step.lengths <- temp$dist[!is.na(temp$dist)]
#   fit_gamma <- fitdist(step.lengths, "gamma", method='mme')
#   ests.gamma <- bootdist(fit_gamma, niter = 1000)
#   gamma.vars[i,1] <- summary(traj.round)$id[i]
#   gamma.vars[i,2] <- ests.gamma$CI[1,1]
#   gamma.vars[i,3] <- ests.gamma$CI[2,1]
# }
# 
# # Estimate Pareto distribution parameters for each season/ID
# pareto.vars <- data.frame(matrix(0,length(traj.round),3))
# colnames(pareto.vars) <- c("ID", "scale", "shape")
# for (i in 1:length(traj.round)) {
#   temp <- traj.round[[i]]
#   step.lengths <- temp$dist[!is.na(temp$dist)]
#   fit_pareto <- fitdist(step.lengths, "pareto", method='mge')
#   ests.pareto <- bootdist(fit_pareto, niter = 1000)
#   pareto.vars[i,1] <- summary(traj.round)$id[i]
#   pareto.vars[i,2] <- ests.pareto$CI[1,1]
#   pareto.vars[i,3] <- ests.pareto$CI[2,1]
# }

dist.vars <- data.frame(matrix(0,length(traj.round),7))
colnames(dist.vars) <- c("ID", "gamma.shape", "gamma.rate", "gamma.GOF", "pareto.scale", "pareto.shape", "pareto.GOF")
for (i in 1:length(traj.round)) {
  temp <- traj.round[[i]]
  step.lengths <- temp$dist[!is.na(temp$dist)]
  fit_gamma <- fitdist(step.lengths, "gamma", method='mme')
  ests.gamma <- bootdist(fit_gamma, niter = 1000)
  fit_pareto <- fitdist(step.lengths, "pareto", method='mge')
  ests.pareto <- bootdist(fit_pareto, niter = 1000)
  gof <- gofstat(list(fit_gamma, fit_pareto), fitnames = c("gamma", "pareto"))
  dist.vars[i,1] <- summary(traj.round)$id[i]
  dist.vars[i,2] <- ests.gamma$CI[1,1]
  dist.vars[i,3] <- ests.gamma$CI[2,1]
  dist.vars[i,4] <- gof$ks[1]
  dist.vars[i,5] <- ests.pareto$CI[1,1]
  dist.vars[i,6] <- ests.pareto$CI[2,1]
  dist.vars[i,7] <- gof$ks[2]
}

#write.csv(dist.vars, "Zebra_StepDist_Parameters.csv")
