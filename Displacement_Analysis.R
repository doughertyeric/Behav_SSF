library(tidyverse)
library(lubridate)

#########################################

name_list = c('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009',
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010',
              'AG253_2010', 'AG255_2010', 'AG256_2010')

zebra09 <- read_csv("Zebra_Anthrax_2009_Cleaned.csv")[,-1]
zebra10 <- read_csv("Zebra_Anthrax_2010_Cleaned.csv")[,-1]
all_zebra <- rbind(zebra09, zebra10)

#########################################

times = c(1,2,3,4,5,6,7)
dist_list <- list()
for (i in 1:length(times)) {
  all_dists <- c()
  skip = ((times[i]*24*3) + 1)
  w = 1
  for (j in 1:length(name_list)) {
    temp <- all_zebra %>% filter(ID == name_list[j])
    for (k in 1:(nrow(temp) - skip)) {
      if (!is.na(temp$x[k]) && !is.na(temp$x[k + skip])) {
        x_diff <- temp$x[k + skip] - temp$x[k]
        y_diff <- temp$y[k + skip] - temp$y[k]
        all_dists[w] <- sqrt(x_diff^2 + y_diff^2)
        w = w + 1
      }
    }
  }
  dist_list[[i]] <- all_dists 
}

##########

library(fitdistrplus)

displacements <- data.frame(matrix(0,7,5))
for (i in 1:length(times)) {
  #plotdist(dist_list[[times[i]]], histo = TRUE, demp = TRUE)
  displacements[i,1] <- descdist(dist_list[[times[i]]], discrete=FALSE, boot=500)[[4]]

  #Fit Gamma distribution to step lengths
  fit_gamma <- fitdist(dist_list[[times[i]]], "gamma", method='mme')
  ests.gamma <- bootdist(fit_gamma, niter = 1000)
  displacements[i,2] <- median(summary(ests.gamma)[1]$estim[,1])
  displacements[i,3] <- median(summary(ests.gamma)[1]$estim[,2])
  displacements[i,4] <- qgamma(0.841, shape=displacements[i,2], rate=displacements[i,3])
  displacements[i,5] <- qgamma(0.975, shape=displacements[i,2], rate=displacements[i,3])
}
colnames(displacements) <- c('mean', 'shape', 'rate', '1SD', '2SD')
write.csv(displacements, "Zebra_Displacement_Parameters.csv")
