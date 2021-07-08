# Part2-Single-season Bayesian occupancy model
# Elildo Carvalho Jr @ ICMBio/CENAP 2020-08-17
# code based on templates from the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" by Marc Kéry & Michael Schaub (2012, Academic Press)
# and on templates from the JAGS translation available at:
# https://www.vogelwarte.ch/de/projekte/publikationen/bpa/code-for-running-bpa-using-jags


##Adapted by Eloisa REBIO Gurupi/ICMBio, 2021-03-10
##objetivo deste scr é analisar a probabilidade de deteccao da onca_pintada em relacao as covariaveis ambientais, antropicas e ecologicas.


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
Ponca2017 <- readRDS(here("data", "Ponca2017team.rds")) 
y <- Ponca2017[,2:15]
y

SiteCovs <- Ponca2017[,16:26]
SiteCovs
SiteCovs["log_biomass_pp"]<-c(log(SiteCovs$photoRates_biomass.x)) ##variavel com distribruiçao lognormal
SiteCovs["log_biomass_sp"]<-c(log(SiteCovs$photoRates_biomass.y)) ##variavel com distribruiçao lognormal

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
biomass_pp <- SiteCovs [,12] #log(kg/trap/day)
biomass_sp <- SiteCovs [,13] #log(kg/trap/day)
fire<-SiteCovs [,9]
landcover<-SiteCovs[,10]
block <- SiteCovs [,11]

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

mean.block <- mean(block, na.rm = TRUE)
sd.block <- sd(block[!is.na(block)])
block <- (block-mean.block)/sd.block     # Standardise block
block[is.na(block)] <- 0               # Impute zeroes (means)
block



##dentre as variaveis com alta correlação (elevation, slope, block, cor>0.7) elevation foi escolhida, 
##pois varios trabalhos citam essa variavel relacionada (negativa) com a OP,
## e tbem elevation teve respostas significativas nas analises de presas prioritarias. 
##como block foi adicionada para compensar possiveis diferenças na probabilidade 
##de detecção, ao considerar eletavion na deteccao o efeito do bloco estaria contemplado.
##(outra altertnativa seria adicionar bloco como variavel randomica)


##sobre as covariaveis ####
##após muitos modelos rodados com as presas prioritárias e com a OP selecionei apenas as variáveis 
##que tiveram alguma influencia sobre a ocupação das presas e que são as mesmas geralmente relacionadas ao uso de habitat da OP.
#são elas: cobertura (basal Area), distancia da agua, elevação,distancia de borda.

#----- 4 - Single-season occupancy model -----

# Specify model in JAGS language
sink(here("bin", "model.jags"))
cat("
model {

# Priors
alpha.psi ~ dnorm(0, 0.01)
#beta1.psi ~ dnorm(0, 0.01) # distwater sitecov 1
#beta2.psi ~ dnorm(0, 0.01) # distEdge
#beta3.psi ~ dnorm(0, 0.01) # elevation sitecov 3
#beta4.psi ~ dnorm(0, 0.01) # slope 4
#beta5.psi ~ dnorm(0, 0.01) # basalArea 5
#beta6.psi ~ dnorm(0, 0.01) # treeDensity 6
#beta7.psi ~ dnorm(0, 0.01) # biomass_pp (priorprey) 7
#beta8.psi ~ dnorm(0, 0.01) # biomass_sp (smallprey) 8
#beta9.psi ~ dnorm(0, 0.01) # fire 9
#beta10.psi ~ dnorm(0, 0.01) #landcover 10


alpha.p ~ dnorm(0, 0.01)
#beta11.p ~ dnorm(0, 0.01) ## block
beta1.p ~ dnorm(0, 0.01) ##distwater sitecov # 
beta5.p ~ dnorm(0, 0.01) ##basalArea
beta2.p ~ dnorm(0, 0.01)##distEdge
beta3.p ~ dnorm(0, 0.01) ##elevation
beta7.p ~ dnorm(0, 0.01)##biomass_pp (priorprey) 7 #glm analysis indicated biomass_pp related to dist_water
beta8.p ~ dnorm(0, 0.01)##biomass_sp (smallprey) 8 #glm analysis indicated biomass_sp related to basal area and elevation



# Likelihood
# Ecological model for the partially observed true state
for (i in 1:R) {
   z[i] ~ dbern(psi[i])                # True occurrence z at site i
   psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
   lpsi.lim[i] <- min(999, max(-999, lpsi[i]))
lpsi[i] <- alpha.psi 
 
   # Observation model for the observations
   for (j in 1:T) {
      y[i,j] ~ dbern(mu.p[i,j])	# Detection-nondetection at i and j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
    ## lp[i,j] <- alpha.p 
 lp[i,j] <- alpha.p + beta1.p*distWater[i] + beta3.p*elevation[i] + beta5.p*basalArea[i] + beta2.p*distEdge[i] + beta7.p*biomass_pp[i] + beta8.p*biomass_sp [i]
   
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])                             # Number of occupied sites
mean.p <- exp(alpha.p) / (1 + exp(alpha.p))    # Sort of average detection
}
",fill = TRUE)
sink()


# Bundle data
dataJAGS <- list(y = y, R = nrow(y), T = ncol(y), distWater=distWater, distEdge=distEdge, elevation=elevation, basalArea=basalArea, biomass_pp=biomass_pp, biomass_sp=biomass_sp)

#

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
zst[zst == -Inf] <- 0
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parameters monitored
params <- c("alpha.psi", "beta1.p", "beta2.p", "beta3.p", "beta5.p", "beta7.p", "beta8.p", "mean.p", "occ.fs", "alpha.p", "z")


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

