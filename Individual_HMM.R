library(tidyverse)
library(sf)
library(moveHMM)

#################################################################
################### Data Import and Prep ########################

zebra09 <- read_csv("Zebra_Anthrax_2009_Cleaned.csv") %>%
  dplyr::select(x,y,date,ID) 
zebra10 <- read_csv("Zebra_Anthrax_2010_Cleaned.csv") %>%
  dplyr::select(x,y,date,ID) 
all.data <- rbind(zebra09, zebra10)

#data09 <- prepData(trackData=data.frame(zebra09), type="UTM", 
#                 coordNames=c('x', 'y'))
#data10 <- prepData(trackData=data.frame(zebra10), type='UTM',
#                  coordNames=c('x', 'y'))

################################################################
###################### Two-State HMM ###########################

mu0 <- c(5,500) # step mean (two parameters: one for each state)
sigma0 <- c(2,250) # step SD
zeromass0 <- c(0,0) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)

angleMean0 <- c(0,0) # angle mean
kappa0 <- c(0.1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m <- fitHMM(data=data09, nbStates=2, stepPar0=stepPar0, anglePar0=anglePar0)

#################################################################
###################### Three-State HMM ##########################

mu0 <- c(5,50,500) # step mean (three parameters: one for each state)
sigma0 <- c(2,20,200) # step SD
zeromass0 <- c(0,0,0) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)

angleMean0 <- c(0,0,0) # angle mean
kappa0 <- c(1,1,0.1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m3 <- fitHMM(data=data09, nbStates=3, stepPar0=stepPar0, anglePar0=anglePar0)

#################################################################
#################### Loop over Individuals ######################

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

distributions <- list()
state.seq <- list()
for (i in 1:length(name_list)) {
  temp <- dplyr::filter(all.data, all.data$ID == name_list[i])
  temp.data <- prepData(trackData=data.frame(temp), type='UTM',
                        coordNames=c('x','y'))
  
  if (i %in% c(1,2,3,4,5,6,7,9,11)) {
    mu0 <- c(5,50,500) # step mean (three parameters: one for each state)
    sigma0 <- c(2,20,200) # step SD
    zeromass0 <- c(0,0,0) # step zero-mass
    stepPar0 <- c(mu0, sigma0, zeromass0)
    angleMean0 <- c(0,0,0) # angle mean
    kappa0 <- c(1,1,0.1) # angle concentration
    anglePar0 <- c(angleMean0,kappa0)
  } else {
    mu0 <- c(5,50,500) # step mean (three parameters: one for each state)
    sigma0 <- c(2,20,200) # step SD
    #zeromass0 <- c(0,0,0) # step zero-mass
    stepPar0 <- c(mu0, sigma0)#, zeromass0)
    angleMean0 <- c(0,0,0) # angle mean
    kappa0 <- c(1,1,0.1) # angle concentration
    anglePar0 <- c(angleMean0,kappa0)
  }
  
  m3 <- fitHMM(data=temp.data, nbStates=3, stepPar0=stepPar0, anglePar0=anglePar0)
  temp.list <- list(m3[[2]]$stepPar, m3[[2]]$anglePar)
  distributions[[i]] <- temp.list
  states <- viterbi(m3)
  state.seq[[i]] <- data.frame(cbind(m3[[1]], states))
}

for (i in 1:length(distributions)) {
  print(distributions[[i]][[1]][1,2])
}

for (i in 1:length(state.seq)) {
  print(table(state.seq[[i]][,3])[2]/nrow(state.seq[[i]]))
}

foraging.only <- list()
for (i in 1:length(state.seq)) {
  temp <- dplyr::filter(state.seq[[i]], state.seq[[i]]$states == 2) %>%
    dplyr::select(., -X.Intercept.)
  foraging.only[[i]] <- temp
}

#######################################################################
###################### Output Shapefiles ##############################

for (i in 1:length(foraging.only)) {
  forage <- dplyr::filter(foraging.only[[i]], !is.na(x)) %>%
    st_as_sf(coords=c('x','y'), crs ='+init=epsg:32733')
  st_write(forage, paste0(name_list[i],"_Foraging.shp"))
}

########################################################################
###################### Step Distributions ##############################

step.dist <- data.frame(matrix(0,11,3))
for (i in 1:length(distributions)) {
  step.mean = distributions[[i]][[1]][1,2]
  step.sd = distributions[[i]][[1]][2,2]
  shape = (step.mean^2) / (step.sd^2)
  rate = step.mean / (step.sd^2)
  step.dist[i,1] <- i
  step.dist[i,2] <- shape
  step.dist[i,3] <- rate
}
colnames(step.dist) <- c("ID", "gamma.shape", "gamma.rate")
write.csv(step.dist, "Foraging_StepDist_Parameters.csv")
