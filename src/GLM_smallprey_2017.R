
library(dplyr)
library(lubridate)
library(here)
library(fitdistrplus)

#----- 2 - Source files-----
here <- here::here # to avoid confusion with "here" function from lubridate
#source(here("bin", "camera trap analysis functions-10-06-18.R")) # using package here to build a path to the subdirectory "bin"
source(here("bin", "ahumada_codes.R"))
source(here("bin", "f-matrix-creator-experimental-probably-ok-but-need-check.R"))
source(here("src", "fix_species_names.R")) # fix some names and remove "false" species


## ----Load data-------

## Read covariate data


# Distance to water
dist.water <- read.csv(here("data", "dist_agua_conv_trsh6_13fev.csv"))
names(dist.water)[1] <- "Camera.Trap.Name"
dist.water$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", dist.water$Camera.Trap.Name)
dist.water$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", dist.water$Camera.Trap.Name)

#Sort by camera trap unit name for analysis
dist.water <- arrange(dist.water, Camera.Trap.Name)
head(dist.water)

# Distance to forest edge
dist.edge <- read.csv(here("data", "dist_pasto10ha.csv"))
names(dist.edge) <- c("Camera.Trap.Name", "dist.to.edge")
dist.edge$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", dist.edge$Camera.Trap.Name)
dist.edge$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", dist.edge$Camera.Trap.Name)


#Sort by camera trap unit name for analysis
dist.edge <- arrange(dist.edge, Camera.Trap.Name)
head(dist.edge)


# elevation
elev <- read.csv(here("data", "slope_elev.csv"))

names(elev)[1] <- "Camera.Trap.Name"
names(elev)[2] <- "slope_pto"
names(elev)[3] <- "elevation_pto"
names(elev)[4] <- "elevation_buf250"
names(elev)[5] <- "elevation_buf500"
elev$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", elev$Camera.Trap.Name)
elev$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", elev$Camera.Trap.Name)

#Sort by camera trap unit name for analysis
elev <- arrange(elev, Camera.Trap.Name)
head(elev)


# tree structure ## excluidas estações amostrais da area queimada
trees <- read.csv(here("data", "tree_structure.csv"))

#Sort by camera trap unit name for analysis
trees <- arrange(trees, Camera.Trap.Name)
head(trees)

##biomass prey prior >15kg /CT/dia##
biomass_sp <- readRDS(here("data", "cov_biomass_sp_2017.rds"))
#biomass_pp <- read.csv(here("data", "biomass_preyprior.csv"))
names(biomass_sp)[1] <- "Camera.Trap.Name"
biomass_sp ["log_biomass"] <- c(log(biomass_sp$photoRates_biomass))
biomass_sp <- arrange(biomass_sp, Camera.Trap.Name)
head(biomass_sp)




## create a single covariates dataframe 
covars <- merge(biomass_sp[,c(1,4)], dist.water[,c(1,4)],  by="Camera.Trap.Name", all.x = T) ##distagua sarivers
covars <- merge(covars, elev[,c(1,4)], by="Camera.Trap.Name", all.x = T) ##buf250
covars <- merge(covars, trees[,c(2,6)], by="Camera.Trap.Name", all.x = T)
covars <- merge(covars, dist.edge[,1:2], by="Camera.Trap.Name", all.x = T) ##ungulates
#covars <- merge(covars, fire[,c(1,2)], by="Camera.Trap.Name", all.x = T)
covars

names(covars)[2] <- "biomass"
names(covars)[3] <- "distWater"
names(covars)[4] <- "elevation"
names(covars)[5] <- "basalArea"
names(covars)[6] <- "distEdge"


saveRDS(covars, here("data","covars_glm_smallprey.rds" ))




###exploratory analysis ####
library(dplyr)
library(lubridate)
library(here)
library(fitdistrplus)

library(usdm) ## 

#load data
covars<- readRDS(here("data","covars_glm_smallprey.rds" ))


head(covars) ### funcao para olhar um conjunto reduzido dos seus dados
str(covars) ## funcao para olhar a estrutura dos seus dados, se eles são numéricos, caracteres, etc
covars[is.na(covars)] <- 0 

covsPonca <- covars


##what distribution probability 
#hist(covsPonca$biomass) #variavel resposta continua - parece lognormal##
hist(covsPonca$biomass) #lognormal##variável resposta log(tx biomass)##ok normalizou

hist(covsPonca$distWater) #continua
hist(covsPonca$elevation) #continua
hist(covsPonca$basalArea) #discreta
hist(covsPonca$distEdge) #continua

hist(basalArea)
hist(distEdge)

##praticamente nenhuma delas segue distribuição normal


## Verificar se existe outlier
boxplot(covsPonca$distWater)## sim 3
boxplot(covsPonca$elevation) ## 1
boxplot(covsPonca$basalArea) ##sim 3
boxplot(covsPonca$distEdge)## não
boxplot(covsPonca$biomass)##2

covsPonca <- data.frame(covsPonca[,c(2:6)])


cor(covsPonca)
?cor
test_normal <- fitdist(covsPonca$biomass, "norm")
plot(test_normal)
test_normal <- fitdist(covsPonca$distWater, "norm")
plot(test_normal)
test_normal <- fitdist(covsPonca$elevation, "norm")
plot(test_normal)
test_normal <- fitdist(covsPonca$basalArea, "nbinom")
plot(test_normal)

?fitdist

##fit better to normal distribution after log biomass transform 
?glm

#models
mdl1<-glm(biomass~distWater+elevation+basalArea+distEdge, data=covsPonca, family = "gaussian")
summary(mdl1)

mdl2<-glm(biomass~distWater, data=covsPonca, family = "gaussian")
summary(mdl2)

mdl3<-glm(biomass~elevation, data=covsPonca, family = "gaussian")

mdl4<-glm(biomass~basalArea, data=covsPonca, family = "gaussian")

mdl5<-glm(biomass~distWater+basalArea+distEdge, data=covsPonca, family = "gaussian")

mdl6<-glm(biomass~elevation+basalArea+distEdge, data=covsPonca, family = "gaussian")

mdl7<-glm(biomass~distWater+basalArea, data=covsPonca, family = "gaussian")

mdl8<-glm(biomass~elevation+basalArea, data=covsPonca, family = "gaussian")

mdl9<-glm(biomass~basalArea+distEdge, data=covsPonca, family = "gaussian")

mdl10<-glm(biomass~distWater+distEdge, data=covsPonca, family = "gaussian")

mdl11<-glm(biomass~elevation+distEdge, data=covsPonca, family = "gaussian")

mdl12<-glm(biomass~distEdge, data=covsPonca, family = "gaussian")

library(MuMIn)
model.sel(mdl1, mdl2, mdl3, mdl4, mdl5, mdl6, mdl7, mdl8, mdl9, mdl10, mdl11, mdl12, rank="AIC")


##best models => mdl2, mdl7, mdl10 (delta<2)

model.matrix(mdl2)
model.matrix((mdl7))

summary(mdl2)
summary(mdl7)
summary(mdl10)

call
lm(biomass~distWater)

