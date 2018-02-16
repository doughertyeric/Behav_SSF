library(tlocoh)
library(adehabitatLT)
library(move)

# Preparing first 2000 points of the toni path for testing
data(toni)
toni <- toni[!duplicated(toni$timestamp.utc),]
toni <- toni[1:2000,]

toni.sp.latlong <- SpatialPoints(toni[ , c("long","lat")], proj4string=CRS("+proj=longlat +ellps=WGS84"))
toni.sp.utm <- spTransform(toni.sp.latlong, CRS("+proj=utm +south +zone=36 +ellps=WGS84"))
toni.mat.utm <- coordinates(toni.sp.utm)
toni.gmt <- as.POSIXct(toni$timestamp.utc, tz="UTC")
local.tz <- "Africa/Johannesburg"
toni.localtime <- as.POSIXct(format(toni.gmt, tz=local.tz), tz=local.tz)

toni.move <- move(x=toni.mat.utm[,1],
                  y=toni.mat.utm[,2],
                  time=toni.localtime,
                  data=toni,
                  animal = as.factor(toni$id),
                  sensor = as.factor(toni$id),
                  proj=CRS("+proj=utm +south +zone=36"))

toni_ltraj  <- as(toni.move, 'ltraj')

#Create step.lengths vector for analysis
step.lengths <- toni_ltraj[[1]]$dist[!is.na(toni_ltraj[[1]]$dist)]

######################################

library(fitdistrplus)

plotdist(step.lengths, histo = TRUE, demp = TRUE)
descdist(step.lengths, discrete=FALSE, boot=500)

#Fit Gamma distribution to step lengths
fit_gamma <- fitdist(step.lengths, "gamma", method='mme')
ests.gamma <- bootdist(fit_gamma, niter = 1000)
summary(ests.gamma)

#Fit Pareto distribution to step lengths
fit_Pareto  <- fitdist(step.lengths, "pareto", method='mge')
ests.Pareto <- bootdist(fit_Pareto, niter = 1000)
summary(ests.Pareto)

#Plot results of both gamma and Pareto distribution fits
plot.legend <- c("gamma", "Pareto")
denscomp(list(fit_gamma, fit_Pareto), legendtext = plot.legend, ylim=c(0,0.002))
cdfcomp(list(fit_gamma, fit_Pareto), legendtext = plot.legend)
qqcomp(list(fit_gamma, fit_Pareto), legendtext = plot.legend)
ppcomp(list(fit_gamma, fit_Pareto), legendtext = plot.legend)

#Compare goodness of fit
gofstat(list(fit_gamma, fit_Pareto), fitnames = c("gamma", "Pareto"))



