# Part2-Single-season Bayesian occupancy model
# Elildo Carvalho Jr @ ICMBio/CENAP 2020-08-17
# code based on templates from the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" by Marc Kéry & Michael Schaub (2012, Academic Press)
# and on templates from the JAGS translation available at:
# https://www.vogelwarte.ch/de/projekte/publikationen/bpa/code-for-running-bpa-using-jags

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
Mazama2017 <- readRDS(here("data", "Mazama_2017team.rds")) 
y <- Mazama2017[,2:15]
str(y)
y
SiteCovs <- Mazama2017[,16:23]
SiteCovs[is.na(SiteCovs)] <- 0 
# check corr in SiteCovs
str(SiteCovs)
cor(SiteCovs)

# elevation, slope and block cor > 0.7
### the elevation was choosed because brought better response in T.terrestris, T Pecari occupancy analisys, and
##it is a covariavel frequently associated in jaguar occupancy analisys, too. 


names(SiteCovs)
distWater <- SiteCovs[,1] ###sarivers adpt
distEdge <- SiteCovs[,2] ##area aberta maior que 10ha mabpiomas 2017
elevation <- SiteCovs[,3] #mean buffer 250m
slope <- SiteCovs[,4] ## mean buffer 250m
basalArea <- SiteCovs[,5]
treeDensity <- SiteCovs[,6]
fire<-SiteCovs [,7]
block<-SiteCovs [,8]

# Standardize covariates
mean.slope <- mean(slope, na.rm = TRUE)
sd.slope <- sd(slope[!is.na(slope)])
slope <- (slope-mean.slope)/sd.slope     # Standardise slope
slope[is.na(slope)] <- 0               # Impute zeroes (means)

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

mean.basalArea <- mean(basalArea, na.rm = TRUE)
sd.basalArea <- sd(basalArea[!is.na(basalArea)])
basalArea <- (basalArea-mean.basalArea)/sd.basalArea     # Standardise basalArea
basalArea[is.na(basalArea)] <- 0               # Impute zeroes (means)

mean.treeDensity <- mean(treeDensity, na.rm = TRUE)
sd.treeDensity <- sd(treeDensity[!is.na(treeDensity)])
treeDensity <- (treeDensity-mean.treeDensity)/sd.treeDensity     # Standardise treeDensity
treeDensity[is.na(treeDensity)] <- 0               # Impute zeroes (means)
treeDensity

mean.fire <- mean(fire, na.rm = TRUE)
sd.fire <- sd(fire[!is.na(fire)])
fire <- (fire-mean.fire)/sd.fire     # Standardise biomass
fire[is.na(fire)] <- 0               # Impute zeroes (means)
fire

mean.block <- mean(block, na.rm = TRUE)
sd.block <- sd(block[!is.na(block)])
block <- (block-mean.block)/sd.block     # Standardise biomass
block[is.na(block)] <- 0               # Impute zeroes (means)
block

#----- 4 - Single-season occupancy model -----

# Specify model in JAGS language
sink(here("bin", "model.jags"))
cat("
model {

# Priors
alpha.psi ~ dnorm(0, 0.01)
beta1.psi ~ dnorm(0, 0.01) # distwater sitecov 1
beta2.psi ~ dnorm(0, 0.01) # distEdge site cov 2
beta3.psi ~ dnorm(0, 0.01) # elevation sitecov 3
#beta4.psi ~ dnorm(0, 0.01) # slope 4
beta5.psi ~ dnorm(0, 0.01) # basalArea 5
#beta6.psi ~ dnorm(0, 0.01) # treeDensity 6
#beta7.psi ~ dnorm(0, 0.01) #fire 7

alpha.p ~ dnorm(0, 0.01)
#beta8.p ~ dnorm(0, 0.01) ##block
#beta1.p ~ dnorm(0, 0.01)
#beta3.p ~ dnorm(0, 0.01)
#beta4.p ~ dnorm(0, 0.01)
#beta5.p ~ dnorm(0, 0.01)

# Likelihood
# Ecological model for the partially observed true state
for (i in 1:R) {
   z[i] ~ dbern(psi[i])                # True occurrence z at site i
   psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
   lpsi.lim[i] <- min(999, max(-999, lpsi[i]))
lpsi[i] <- alpha.psi + beta1.psi*distWater[i] + beta2.psi*distEdge[i] + beta3.psi*elevation[i] +beta5.psi*basalArea[i]


   # Observation model for the observations
   for (j in 1:T) {
      y[i,j] ~ dbern(mu.p[i,j])	# Detection-nondetection at i and j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
   lp[i,j] <- alpha.p ##+ beta3.p*elevation[i]
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])                             # Number of occupied sites
mean.p <- exp(alpha.p) / (1 + exp(alpha.p))    # Sort of average detection
}
",fill = TRUE)
sink()


# Bundle data
dataJAGS <- list(y = y, R = nrow(y), T = ncol(y), distWater=distWater, distEdge=distEdge, elevation=elevation, basalArea=basalArea)

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
zst[zst == -Inf] <- 0
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parameters monitored

params <- c("alpha.psi", "beta1.psi", "beta2.psi", "beta3.psi","beta5.psi", "mean.p", "occ.fs", "alpha.p", "z")

# MCMC settings
ni <- 50000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(dataJAGS, inits, params, here("bin", "model.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)

# Posterior distribution of the number of occupied sites in actual sample
hist(out$BUGSoutput$sims.list$occ.fs, nclass = 30, col = "gray", main = "", xlab = "Number of occupied sites", )
#abline(v = 10, lwd = 2) # The observed number


