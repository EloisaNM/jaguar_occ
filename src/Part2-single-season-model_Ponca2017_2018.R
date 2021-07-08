# Part2-Single-season Bayesian occupancy model
# Elildo Carvalho Jr @ ICMBio/CENAP 2020-08-17
# code based on templates from the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" by Marc Kéry & Michael Schaub (2012, Academic Press)
# and on templates from the JAGS translation available at:
# https://www.vogelwarte.ch/de/projekte/publikationen/bpa/code-for-running-bpa-using-jags


##Adapted by Eloisa REBIO Gurupi/ICMBio, 2021-03-10
##objetivo deste scr é analisar a probabilidade de ocupação da onca_pintada em relacao as covariaveis ambientais, antropicas e ecologicas.


#----- 1 - Load libraries-----
#library(dplyr)
#library(lubridate)
library(here)
library(R2jags)
library(rjags)
library(ggplot2)


#----- 2 - Source files-----
#source(here("bin", "ahumada_codes.R"))
#source(here("bin", "f-matrix-creator-experimental-probably-ok-but-need-check.R"))
#source(here("src", "fix_species_names.R")) # fix some names and remove "false" species


#----- 3 - Read and prepare data -----
Ponca_17_18 <- readRDS(here("data", "Ponca_17_18_team.rds")) 

y <- Ponca_17_18[,2:29]


SiteCovs <- Ponca_17_18[,30:42]
SiteCovs["log_biomass_pp_17"]<-c(log(SiteCovs$photoRates_biomass_pp_17))
SiteCovs["log_biomass_sp_17"]<-c(log(SiteCovs$photoRates_biomass_sp_17))

SiteCovs["log_biomass_pp_18"]<-c(log(SiteCovs$photoRates_biomass_pp_18))
SiteCovs["log_biomass_sp_18"]<-c(log(SiteCovs$photoRates_biomass_sp_18))
SiteCovs
#SiteCovs["log_distWater"]<- c(log(SiteCovs$dist_agua_sarivers))
#SiteCovs["log_distEdge"]<- c(log(SiteCovs$dist.to.edge))
SiteCovs[is.na(SiteCovs)] <- 0 
# check corr in SiteCovs
str(SiteCovs)
cor(SiteCovs)

# covariaveis com alta correlação: elevation, slope, block => indice de correlação >0.7
##
####
###covariaveis com média correlação elevation, slope e dist water => indice de correlação> 0.5
##

names(SiteCovs)

distWater <- SiteCovs[,1]
distEdge <- SiteCovs[,2]
elevation <- SiteCovs[,3]
slope <- SiteCovs[,4]
basalArea <- SiteCovs[,5]
treeDensity <- SiteCovs[,6]
biomass_pp_17 <- SiteCovs [,14]
biomass_sp_17 <- SiteCovs [,15]
fire<-SiteCovs [,11]
landcover<-SiteCovs[,12]
block <- SiteCovs [,13]
biomass_pp_18 <- SiteCovs [,16]
biomass_sp_18 <- SiteCovs [,17]
# Standardize covariates

mean.distWater <- mean(distWater, na.rm = TRUE)
sd.distWater <- sd(distWater[!is.na(distWater)])
distWater <- (distWater-mean.distWater)/sd.distWater     # Standardise distWater
distWater[is.na(distWater)] <- 0               # Impute zeroes (means)

mean.distEdge <- mean(distEdge, na.rm = TRUE)
sd.distEdge <- sd(distEdge[!is.na(distEdge)])
distEdge <- (distEdge-mean.distEdge)/sd.distEdge     # Standardise distEdge
distEdge[is.na(distEdge)] <- 0               # Impute zeroes (means)

mean.elevation <- mean(elevation, na.rm = TRUE)
sd.elevation <- sd(elevation[!is.na(elevation)])
elevation <- (elevation-mean.elevation)/sd.elevation     # Standardise elevation
elevation[is.na(elevation)] <- 0               # Impute zeroes (means)

mean.slope <- mean(slope, na.rm = TRUE)
sd.slope <- sd(slope[!is.na(slope)])
slope <- (slope-mean.slope)/sd.slope     # Standardise slope
slope[is.na(slope)] <- 0               # Impute zeroes (means)

mean.basalArea <- mean(basalArea, na.rm = TRUE)
sd.basalArea <- sd(basalArea[!is.na(basalArea)])
basalArea <- (basalArea-mean.basalArea)/sd.basalArea     # Standardise basalArea
basalArea[is.na(basalArea)] <- 0               # Impute zeroes (means)

mean.treeDensity <- mean(treeDensity, na.rm = TRUE)
sd.treeDensity <- sd(treeDensity[!is.na(treeDensity)])
treeDensity <- (treeDensity-mean.treeDensity)/sd.treeDensity     # Standardise treeDensity
treeDensity[is.na(treeDensity)] <- 0               # Impute zeroes (means)
treeDensity

mean.biomass_pp <- mean(biomass_pp, na.rm = TRUE)
sd.biomass_pp <- sd(biomass_pp[!is.na(biomass_pp)])
biomass_pp <- (biomass_pp-mean.biomass_pp)/sd.biomass_pp     # Standardise biomass_pp
biomass_pp[is.na(biomass_pp)] <- 0               # Impute zeroes (means)
biomass_pp

mean.biomass_sp <- mean(biomass_sp, na.rm = TRUE)
sd.biomass_sp <- sd(biomass_sp[!is.na(biomass_sp)])
biomass_sp <- (biomass_sp-mean.biomass_sp)/sd.biomass_sp     # Standardise biomass_sp
biomass_sp[is.na(biomass_sp)] <- 0               # Impute zeroes (means)
biomass_sp

mean.biomass_tt <- mean(biomass_tt, na.rm = TRUE)
sd.biomass_tt <- sd(biomass_sp[!is.na(biomass_tt)])
biomass_tt <- (biomass_sp-mean.biomass_tt)/sd.biomass_tt     # Standardise biomass_sp
biomass_tt[is.na(biomass_tt)] <- 0               # Impute zeroes (means)
biomass_tt


mean.fire <- mean(fire, na.rm = TRUE)
sd.fire <- sd(fire[!is.na(fire)])
fire <- (fire-mean.fire)/sd.fire     # Standardise biomass
fire[is.na(fire)] <- 0               # Impute zeroes (means)
fire

mean.landcover <- mean(landcover, na.rm = TRUE)
sd.landcover <- sd(landcover[!is.na(landcover)])
landcover <- (landcover-mean.landcover)/sd.landcover     # Standardise lancover
landcover[is.na(landcover)] <- 0               # Impute zeroes (means)
landcover



##dentre as variaveis com alta correlação (cor>0.7) elevation foi escolhida, pois varios trabalhos citam essa variavel relacionada (negativa) com a OP,

#----- 4 - Single-season occupancy model -----

# Specify model in JAGS language

# Specify model in JAGS language
sink(here("bin", "model.jags"))
cat("
model {

# Priors
alpha.psi ~ dnorm(0, 0.01)
#beta1.psi ~ dnorm(0, 0.01) # distwater sitecov 1
beta2.psi ~ dnorm(0, 0.01) # distEdge
beta3.psi ~ dnorm(0, 0.01) # elevation sitecov 3
#beta4.psi ~ dnorm(0, 0.01) # slope 4
beta5.psi ~ dnorm(0, 0.01) # basalArea 5
beta6.psi ~ dnorm(0, 0.01) # treeDensity 6
#beta7.psi ~ dnorm(0, 0.01) # biomass_pp (priorprey) 7
#beta8.psi ~ dnorm(0, 0.01) # biomass_sp (smallprey) 8
#beta9.psi ~ dnorm(0, 0.01) # fire 9
#beta10.psi ~ dnorm(0, 0.01) #landcover 10

alpha.p ~ dnorm(0, 0.01)
beta1.p ~ dnorm(0, 0.01) ##distwater sitecov # 
#beta5.p ~ dnorm(0, 0.01) ##basalArea
#beta3.p ~ dnorm(0, 0.01) ##elevation sitecov 3
#beta2.p ~ dnorm(0, 0.01)##distEdge
#beta8.p ~ dnorm(0, 0.01)
#beta7.p ~ dnorm(0, 0.01)

##random block effects
  for (i in 1:nblock){
         beta11.p[i] ~ dnorm(0, 2)
         beta11.psi[i] ~ dnorm(0, 2)
 } #i


# Likelihood
# Ecological model for the partially observed true state
for (i in 1:R) {
   z[i] ~ dbern(psi[i])                # True occurrence z at site i
   psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
   lpsi.lim[i] <- min(999, max(-999, lpsi[i]))
lpsi[i] <- alpha.psi + beta11.psi[block[i]]+beta2.psi*distEdge[i]+ beta3.psi*elevation[i] + beta5.psi*basalArea[i] + beta6.psi*treeDensity[i] #+ beta7.psi*biomass_pp[i] + beta8.psi*biomass_sp [i]


   # Observation model for the observations
   for (j in 1:T) {
      y[i,j] ~ dbern(mu.p[i,j])	# Detection-nondetection at i and j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
   lp[i,j] <- alpha.p + beta11.p[block[i]] + beta1.p*distWater[i]
   
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])                             # Number of occupied sites
mean.p <- exp(alpha.p) / (1 + exp(alpha.p))    # Sort of average detection
}
",fill = TRUE)
sink()


# Bundle data
#dataJAGS <- list(y = y, R = nrow(y), T = ncol(y), nblock=length(unique(block)), block=block, distWater=distWater, elevation=elevation, basalArea=basalArea, distEdge=distEdge, biomass_pp=biomass_pp, biomass_sp=biomass_sp)
dataJAGS <- list(y = y, R = nrow(y), T = ncol(y), nblock=length(unique(block)), block=block, distWater=distWater, elevation=elevation, basalArea=basalArea, distEdge=distEdge, treeDensity=treeDensity)

#

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
zst[zst == -Inf] <- 0
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parameters monitored
#params <- c("alpha.psi",  "beta11.p","beta11.psi", "beta1.p", "beta2.psi","beta5.psi", "beta3.psi","beta7.psi", "beta8.psi","mean.p", "occ.fs", "alpha.p", "z")
params <- c("alpha.psi",  "beta11.p","beta11.psi", "beta1.p", "beta2.psi","beta5.psi","beta6.psi", "beta3.psi","mean.p", "occ.fs", "alpha.p", "z")


# MCMC settings
ni <- 50000
nt <- 10
nb <- 20000
nc <- 3
#ni <- 200000
#nt <- 250
#nb <- 100000
#nc <- 3


# Call JAGS from R (BRT < 1 min)
out <- jags(dataJAGS, inits, params, here("bin", "model.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)

# Posterior distribution of the number of occupied sites in actual sample
hist(out$BUGSoutput$sims.list$occ.fs, nclass = 30, col = "gray", main = "", xlab = "Number of occupied sites", )
#abline(v = 10, lwd = 2) # The observed number

##model 1 (psi-constante, p~allcovariates): distWater was significtantly negative related to P.onca detection,
##model 2 (p~distwater, psi~all others covariates): distWater was significtive negative related to P.onca detection
##and biomass_pp was significative and positive related to Ponca occupancy 


# Distance to Water #distWater <- SiteCovs[,1] #beta1.psi ~ dnorm(0, 0.01) # dist.water
mcmc.sample <- out$BUGSoutput$n.sims
original.distWater <- SiteCovs[,1]
original.distWater.pred <- seq(min(original.distWater), max(original.distWater), length.out = 30)
distWater.pred <- (original.distWater.pred - mean.distWater)/sd.distWater
psi.pred.distWater <- rep(NA, length(distWater.pred))
for(i in 1:length(psi.pred.distWater)) {
  psi.pred.distWater[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta1.psi*distWater.pred[i])
}
array.psi.pred.distWater <- array(NA, dim = c(length(distWater.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.distWater[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta1.psi[i]*distWater.pred)
}

# Plot for a subsample of MCMC draws
# write as a function:
plot.distWater <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.distWater.pred, psi.pred.distWater, main = "", ylab = "Occupancy probability", xlab = "Distance to water (m)", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.distWater.pred, array.psi.pred.distWater[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.distWater.pred, psi.pred.distWater, type = "l", lwd = 3, col = "blue")
}
plot.distWater()

# save as jpeg
jpeg(here("results", "distWater_Ponca_11maio.jpg"), width = 800, height = 400) # Open jpeg file
plot.distWater()
dev.off()


# Distance to edge #distEdge <- SiteCovs[,2] #beta2.psi ~ dnorm(0, 0.01) # dist.edge
mcmc.sample <- out$BUGSoutput$n.sims
original.distEdge <- SiteCovs[,2]
original.distEdge.pred <- seq(min(original.distEdge), max(original.distEdge), length.out = 30)
distEdge.pred <- (original.distEdge.pred - mean.distEdge)/sd.distEdge
psi.pred.distEdge <- rep(NA, length(distEdge.pred))
for(i in 1:length(psi.pred.distEdge)) {
  psi.pred.distEdge[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta2.psi*distEdge.pred[i])
}
array.psi.pred.distEdge <- array(NA, dim = c(length(distEdge.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.distEdge[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta2.psi[i]*distEdge.pred)
}

# Plot for a subsample of MCMC draws
# write as a function:
plot.distEdge <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.distEdge.pred, psi.pred.distEdge, main = "", ylab = "Occupancy probability", xlab = "Distance to Edge (m)", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.distEdge.pred, array.psi.pred.distEdge[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.distEdge.pred, psi.pred.distEdge, type = "l", lwd = 3, col = "blue")
}

plot.distEdge()

# save as jpeg
jpeg(here("results", "dist_edge_Ponca_30mar.jpeg"), width = 800, height = 400) # Open jpeg file
plot.distEdge()
dev.off()

# #elevation <- SiteCovs[,3] #beta3
mcmc.sample <- out$BUGSoutput$n.sims
original.elevation <- SiteCovs[,3]
original.elevation.pred <- seq(min(original.elevation), max(original.elevation), length.out = 30)
elevation.pred <- (original.elevation.pred - mean.elevation)/sd.elevation
psi.pred.elevation <- rep(NA, length(elevation.pred))
for(i in 1:length(psi.pred.elevation)) {
  psi.pred.elevation[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta3.psi*elevation.pred[i])
}
array.psi.pred.elevation <- array(NA, dim = c(length(elevation.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.elevation[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta3.psi[i]*elevation.pred)
}

# Plot for a subsample of MCMC draws
# write as a function:
plot.elevation <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.elevation.pred, psi.pred.elevation, main = "", ylab = "Occupancy probability", xlab = "elevation (m)", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.elevation.pred, array.psi.pred.elevation[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.elevation.pred, psi.pred.elevation, type = "l", lwd = 3, col = "blue")
}

plot.elevation()

# save as jpeg
jpeg(here("results", "elevation_effect_Ponca12jan.jpeg"), width = 800, height = 400) # Open jpeg file
plot.elevation()
dev.off()


# SLOPE # slope <- SiteCovs[,4]###beta4

mcmc.sample <- out$BUGSoutput$n.sims
original.slope <- SiteCovs[,4]
original.slope.pred <- seq(min(original.slope), max(original.slope), length.out = 30)
slope.pred <- (original.slope.pred - mean.slope)/sd.slope
psi.pred.slope <- rep(NA, length(slope.pred))
for(i in 1:length(psi.pred.slope)) {
  psi.pred.slope[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta4.psi*slope.pred[i])
}
array.psi.pred.slope <- array(NA, dim = c(length(slope.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.slope[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta4.psi[i]*slope.pred)
}

# Plot for a subsample of MCMC draws
# write as a function:
plot.slope <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.slope.pred, psi.pred.slope, main = "", ylab = "Occupancy probability", xlab = "Slope (%)", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.slope.pred, array.psi.pred.slope[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.slope.pred, psi.pred.slope, type = "l", lwd = 3, col = "blue")
}
plot.slope()

## save as jpeg
jpeg(here("results", "slope_effect_Ponca14jan.jpg"), width = 800, height = 400) # Open jpeg file
plot.slope()
dev.off()


# treeBurned ###beta5 #treeBurned - SiteCovs[,5] 
mcmc.sample <- out$BUGSoutput$n.sims
original.treeBurned <- SiteCovs[,5]
original.treeBurned.pred <- seq(min(original.treeBurned), max(original.treeBurned), length.out = 30)
treeBurned.pred <- (original.treeBurned.pred - mean.treeBurned)/sd.treeBurned
psi.pred.treeBurned <- rep(NA, length(treeBurned.pred))
for(i in 1:length(psi.pred.treeBurned)) {
  psi.pred.treeBurned[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta5.psi*treeBurned.pred[i])
}
array.psi.pred.treeBurned <- array(NA, dim = c(length(treeBurned.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.treeBurned[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta5.psi[i]*treeBurned.pred)
}

# Plot for a subsample of MCMC draws
# write as a function:
plot.treeBurned <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.treeBurned.pred, psi.pred.treeBurned, main = "", ylab = "Occupancy probability", xlab = "Burned effect", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.treeBurned.pred, array.psi.pred.treeBurned[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.treeBurned.pred, psi.pred.treeBurned, type = "l", lwd = 3, col = "blue")
}

plot.treeBurned()

# save as jpeg
jpeg(here("results", "treeBurned_Ponca.jpg"), width = 800, height = 400) # Open jpeg file
plot.treeBurned()
dev.off()



# Basal area ###beta5 #basalArea - SiteCovs[,6] 
mcmc.sample <- out$BUGSoutput$n.sims
original.basalArea <- SiteCovs[,5]
original.basalArea.pred <- seq(min(original.basalArea), max(original.basalArea), length.out = 30)
basalArea.pred <- (original.basalArea.pred - mean.basalArea)/sd.basalArea
psi.pred.basalArea <- rep(NA, length(basalArea.pred))
for(i in 1:length(psi.pred.basalArea)) {
  psi.pred.basalArea[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta5.psi*basalArea.pred[i])
}
array.psi.pred.basalArea <- array(NA, dim = c(length(basalArea.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.basalArea[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta5.psi[i]*basalArea.pred)
}

# Plot for a subsample of MCMC draws
# write as a function:
plot.basalArea <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.basalArea.pred, psi.pred.basalArea, main = "", ylab = "Occupancy probability", xlab = "Basal area of trees (m2/ha)", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.basalArea.pred, array.psi.pred.basalArea[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.basalArea.pred, psi.pred.basalArea, type = "l", lwd = 3, col = "blue")
}

plot.basalArea()

# save as jpeg
jpeg(here("results", "basalArea_Ponca.jpg_11maio"), width = 800, height = 400) # Open jpeg file
plot.basalArea()
dev.off()


##beta7.psi ~ dnorm(0, 0.01) # treeDensity #treeDensity <- SiteCovs[,7]

mcmc.sample <- out$BUGSoutput$n.sims
original.treeDensity <- SiteCovs[,7]
original.treeDensity.pred <- seq(min(original.treeDensity), max(original.treeDensity), length.out = 30)
treeDensity.pred <- (original.treeDensity.pred - mean.treeDensity)/sd.treeDensity
psi.pred.treeDensity <- rep(NA, length(treeDensity.pred))
for(i in 1:length(psi.pred.treeDensity)) {
  psi.pred.treeDensity[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta7.psi*treeDensity.pred[i])
}
array.psi.pred.treeDensity <- array(NA, dim = c(length(treeDensity.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.treeDensity[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta7.psi[i]*treeDensity.pred)
}


# Plot for a subsample of MCMC draws
# write as a function:
plot.treeDensity <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.treeDensity.pred, psi.pred.treeDensity, main = "", ylab = "Occupancy probability", xlab = "treeDensity", cex.lab=1.2, cex.axis=1.2, ylim = c(0, 1), type = "l", lwd = 3, las=1,)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.treeDensity.pred, array.psi.pred.treeDensity[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.treeDensity.pred, psi.pred.treeDensity, type = "l", lwd = 3, col = "blue")
}

plot.treeDensity()

# save as jpeg
jpeg(here("results", "treeDensity_Ponca.jpg"), width = 800, height = 400) # Open jpeg file
plot.treeDensity()
dev.off()


# biomass_pp #beta7.psi #biomass <- SiteCovs [,7]

mcmc.sample <- out$BUGSoutput$n.sims
original.biomass_pp <- SiteCovs[,7]
original.biomass_pp.pred <- seq(min(original.biomass_pp), max(original.biomass_pp), length.out = 30)
biomass_pp.pred <- (original.biomass_pp.pred - mean.biomass_pp)/sd.biomass_pp
psi.pred.biomass_pp <- rep(NA, length(biomass_pp.pred))
for(i in 1:length(psi.pred.biomass_pp)) {
  psi.pred.biomass_pp[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta7.psi*biomass_pp.pred[i])
}
array.psi.pred.biomass_pp <- array(NA, dim = c(length(biomass_pp.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.biomass_pp[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta7.psi[i]*biomass_pp.pred)
}


# Plot for a subsample of MCMC draws
# write as a function:
plot.biomass_pp <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.biomass_pp.pred, psi.pred.biomass_pp, main = "", ylab = "Occupancy probability", xlab = "PreyBiomass(Kg/CT/dia)", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.biomass_pp.pred, array.psi.pred.biomass_pp[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.biomass_pp.pred, psi.pred.biomass_pp, type = "l", lwd = 3, col = "blue")
}

plot.biomass_pp()

# save as jpeg
jpeg(here("results", "preyprior_biomass_28abril.jpg"), width = 800, height = 400) # Open jpeg file
plot.biomass_pp()
dev.off()

# biomass_sp #beta9.psi #biomass <- SiteCovs [,9]

mcmc.sample <- out$BUGSoutput$n.sims
original.biomass_sp <- SiteCovs[,8]
original.biomass_sp.pred <- seq(min(original.biomass_sp), max(original.biomass_sp), length.out = 30)
biomass_sp.pred <- (original.biomass_sp.pred - mean.biomass_sp)/sd.biomass_sp
psi.pred.biomass_sp <- rep(NA, length(biomass_sp.pred))
for(i in 1:length(psi.pred.biomass_sp)) {
  psi.pred.biomass_sp[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta8.psi*biomass_sp.pred[i])
}
array.psi.pred.biomass_sp <- array(NA, dim = c(length(biomass_sp.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.biomass_sp[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta8.psi[i]*biomass_sp.pred)
}


# Plot for a subsample of MCMC draws
# write as a function:
plot.biomass_sp <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.biomass_sp.pred, psi.pred.biomass_sp, main = "", ylab = "Occupancy probability", xlab = "samllPreyBiomass(Kg)", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.biomass_sp.pred, array.psi.pred.biomass_sp[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.biomass_sp.pred, psi.pred.biomass_sp, type = "l", lwd = 3, col = "blue")
}

plot.biomass_sp()

# save as jpeg
jpeg(here("results", "smallprey_biomass_30mar.jpg"), width = 800, height = 400) # Open jpeg file
plot.biomass_sp()
dev.off()


# ndvi #beta10.psi # <- SiteCovs [,10]

mcmc.sample <- out$BUGSoutput$n.sims
original.ndvi <- SiteCovs[,10]
original.ndvi.pred <- seq(min(original.ndvi), max(original.ndvi), length.out = 30)
ndvi.pred <- (original.ndvi.pred - mean.ndvi)/sd.ndvi
psi.pred.ndvi <- rep(NA, length(ndvi.pred))
for(i in 1:length(psi.pred.ndvi)) {
  psi.pred.ndvi[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta10.psi*ndvi.pred[i])
}
array.psi.pred.ndvi <- array(NA, dim = c(length(ndvi.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.ndvi[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta10.psi[i]*ndvi.pred)
}


# Plot for a subsample of MCMC draws
# write as a function:
plot.ndvi <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.ndvi.pred, psi.pred.ndvi, main = "", ylab = "Occupancy probability", xlab = "NDVI", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.ndvi.pred, array.psi.pred.ndvi[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.ndvi.pred, psi.pred.ndvi, type = "l", lwd = 3, col = "blue")
}

plot.ndvi()

# save as jpeg
jpeg(here("results", "ndvi_Ponca_14jan.jpg"), width = 800, height = 400) # Open jpeg file
plot.ndvi()
dev.off()





##fire
mcmc.sample <- out$BUGSoutput$n.sims
original.fire <- SiteCovs[,11]
original.fire.pred <- seq(min(original.fire), max(original.fire), length.out = 30)
fire.pred <- (original.fire.pred - mean.fire)/sd.fire
psi.pred.fire <- rep(NA, length(fire.pred))
for(i in 1:length(psi.pred.fire)) {
  psi.pred.fire[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta11.psi*fire.pred[i])
}
array.psi.pred.fire <- array(NA, dim = c(length(fire.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.psi.pred.fire[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta11.psi[i]*fire.pred)
}

# Plot for a subsample of MCMC draws
# write as a function:
plot.fire <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.fire.pred, psi.pred.fire, main = "", ylab = "Occupancy probability", xlab = "fire frequency (m)", cex.lab=1.2, cex.axis=1.2,  ylim = c(0, 1), type = "l", lwd = 3, las=1)#frame.plot = FALSE)
  for (i in sub.set){
    lines(original.fire.pred, array.psi.pred.fire[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.fire.pred, psi.pred.fire, type = "l", lwd = 3, col = "blue")
}
plot.fire()

# save as jpeg
jpeg(here("results", "fire_Ponca_14jan.jpg"), width = 800, height = 400) # Open jpeg file
plot.fire()
dev.off()



