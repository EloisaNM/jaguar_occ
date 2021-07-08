
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

#load data
covars<- readRDS(here("data","covars_pp_glm.rds"))
covars_sp<- readRDS(here("data","covars_glm_smallprey.rds"))
names(covars_sp)[2] <- "biomass_sp" ##log_biomass_sp

covars<- merge(covars, covars_sp[,c(1,2)], by="Camera.Trap.Name", all.x = T, all.y = T)

head(covars) ### funcao para olhar um conjunto reduzido dos seus dados
str(covars) ## funcao para olhar a estrutura dos seus dados, se eles são numéricos, caracteres, etc
covars[is.na(covars)] <- 0 

covsPonca <- covars


##what distribution probability 
hist(covsPonca$biomass) #continua+##lognormal##variável resposta GLM analysis (log(tx biomass))
hist(covsPonca$log_biomass_pp)##ok normalizou +/-
hist(covsPonca$distWater) #continua +
hist(covsPonca$elevation) #continua +
hist(covsPonca$basalArea) #continua +
hist(covsPonca$distEdge) #continua +
hist(covsPonca$biomass_sp) ##normal #(log(tx biomass_sp))
##

## Verificar se existe outlier
boxplot(covsPonca$distWater)## 3
boxplot(covsPonca$elevation) ## 3
boxplot(covsPonca$basalArea) ##sim 3
boxplot(covsPonca$distEdge)## 1
boxplot(covsPonca$biomass)##2
boxplot(covsPonca$log_biomass_pp) #n
boxplot(covsPonca$biomass_sp) #n


covsPonca <- data.frame(covsPonca[,c(2:8)])

##Pearson's correlation 

cor(covsPonca)

##test normal variaveis respostas (biomass_pp e biomass_sp)
test_normal <- fitdist(covsPonca$biomass, "norm")
plot(test_normal)
test_normal <- fitdist(covsPonca$log_biomass_pp, "norm")
plot(test_normal)
test_normal <- fitdist(covsPonca$biomass_sp, "norm") 
plot(test_normal)


##fit better to normal distribution after log transformed 

