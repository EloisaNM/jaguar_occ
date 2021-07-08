# Part1-Prepare data

# Read in a TEAM data set and create and format so it is ready for wildlife community model
# Written by Jorge Ahumada @ Conservation International
# Adapted by Elildo Carvalho Jr @ ICMBio/CENAP, 2020-04-02

##Adapted by Eloisa REBIO Gurupi/ICMBio, 2021-03-10
##objetivo deste scr é gerar uma matriz de detecção com 12 ocasiões de 5 dias, da onca-pintada
#associada às covariáveis ambientais, antrópicas,e ecológicas (biomassa de presas grandes e de presas pequenas/CT/dia).
##considerando a disposição das estações amostrais em dois blocos, foi gerada também a covariavel (block) 
##relacionada a probabilidade de detecção.


#----- 1 - Load libraries-----
library(dplyr)
library(lubridate)
library(here)


#----- 2 - Source files-----
here <- here::here # to avoid confusion with "here" function from lubridate
#source(here("bin", "camera trap analysis functions-10-06-18.R")) # using package here to build a path to the subdirectory "bin"
source(here("bin", "ahumada_codes.R"))
source(here("bin", "f-matrix-creator-experimental-probably-ok-but-need-check.R"))
source(here("src", "fix_species_names.R")) # fix some names and remove "false" species


## ----Load data-------
dataRBG <- read.csv(here("data", "Wild_ID_RBG_2016to2019.csv"))


# some fixes
dataRBG$Sampling.Unit.Name <- as.factor(dataRBG$Camera.Trap.Name)
colnames(dataRBG)[9] <- "Photo.Time"
dataRBG$bin <- factor(dataRBG$bin)

# fix date formats (only needed if data was read with read.csv instead of f.readin.fix.data)
dataRBG$Photo.Date <- as.Date(dataRBG$Photo.Date)
dataRBG$Camera.Start.Date <- as.Date(dataRBG$Camera.Start.Date)
dataRBG$Camera.End.Date <- as.Date(dataRBG$Camera.End.Date)
dataRBG$Start.Date <- as.Date(dataRBG$Start.Date)
dataRBG$End.Date <- as.Date(dataRBG$End.Date)


# fix species names
f.fix.species.names(dataRBG)
dataRBG <- dataTemp # use new df created by function

#----- 4 - Extract binary presence/absence matrices for each species
species <- unique(dataRBG$bin)
cams <- unique(dataRBG$Camera.Trap.Name)
years <- unique(dataRBG$Sampling.Event)
secondPeriods <- 1:14 # we will use 11 occasions

# Separate different years - clunky?
dataRBG2016 <- dplyr::filter(dataRBG, Sampling.Event == 2016)
dataRBG2017 <- dplyr::filter(dataRBG, Sampling.Event == 2017)
dataRBG2018 <- dplyr::filter(dataRBG, Sampling.Event == 2018)
dataRBG2019 <- dplyr::filter(dataRBG, Sampling.Event == 2019)

# use only first 55 days of sampling
# since we are using only 1st 55 days, we must reset end dates using max photo date
f.update.end.data <- function(data, duration){
  new.end.date <- min(data$Start.Date)+duration
  df1 <- subset(data, Photo.Date <= new.end.date)
  for(i in 1:nrow(df1)){ 
    if (df1$End.Date[i] > new.end.date) {
      df1$End.Date[i] <- new.end.date
    }
  }
  df1$Camera.Start.Date <- df1$Start.Date
  df1$Camera.End.Date <- df1$End.Date
  assign("df1", df1, envir=.GlobalEnv)
} # End of function


f.update.end.data(dataRBG2017, 56)
dataRBG2017 <- df1
f.update.end.data(dataRBG2018, 56)
dataRBG2018 <- df1



# Create presence/absence matrices for each species each year
# matrix dimensions are all identical accross species and years

#paMats2016 <- f.matrix.creator3(dataRBG2016, cams, species)
#paMats2017 <- f.matrix.creator3(dataRBG2017, cams, species)
#paMats2018 <- f.matrix.creator3(dataRBG2018, cams, species)
#paMats2019 <- f.matrix.creator3(dataRBG2019, cams, species)

# before using f.matrix.creator check sampling duration
duration <- function(data) {
  sampling.days <- max(data$End.Date) - min(data$Start.Date) + 1
  return(sampling.days)
}

# 55 days ok
duration(dataRBG2017) 
duration(dataRBG2018) 

paMats2017 <- f.matrix.creator4(dataRBG2017, species, 14)
paMats2018 <- f.matrix.creator4(dataRBG2018, species, 14)

names(paMats2017)
names(paMats2018)

##"panthera onca 20


createSppData <- function(x) {
  for(i in 20:length(x)){
    df2 <- as.data.frame(paMats2017[x]) 
    colnames(df2) <- seq(20:length(colnames(df2))); colnames(df2) <- paste("X2017.", colnames(df2), sep="")
    df3 <- as.data.frame(paMats2018[x]) 
    colnames(df3) <- seq(20:length(colnames(df3))); colnames(df3) <- paste("X2018.", colnames(df3), sep="")
    bla <- cbind(df2, df3)
  }
  assign(paste("dataRBG_species", gsub(" ", "_", x), sep="_"), bla, envir = .GlobalEnv)
}

# check if it works
createSppData("Panthera onca")
dataRBG_species_Panthera_onca

dataRBG_species_Panthera_onca$Camera.Trap.Name <- rownames(dataRBG_species_Panthera_onca)






#----- 4 - Read covariate data

# Land cover Mapbiomas
cover <- read.csv(here("data", "mapbiom2017_cover500m.csv"))
names(cover)[1] <- "Camera.Trap.Name"
names(cover)[2] <- "landCover.500m.17"
cover$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", cover$Camera.Trap.Name)
cover$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", cover$Camera.Trap.Name)

#Sort by camera trap unit name for analysis
cover <- arrange(cover, Camera.Trap.Name)
head(cover)


##fire_2007_2016
fire <- read.csv(here("data", "fire_density_allyears_team2017.csv"))
names(fire)[1] <- "Camera.Trap.Name"
names(fire)[2] <- "fire_500m"
fire$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", fire$Camera.Trap.Name)
fire$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", fire$Camera.Trap.Name)

#Sort by camera trap unit name for analysis
fire <- arrange(fire, Camera.Trap.Name)
head(fire)


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
names(elev)[2] <- "slope"
names(elev)[3] <- "elevation_pto"
elev$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", elev$Camera.Trap.Name)
elev$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", elev$Camera.Trap.Name)

#Sort by camera trap unit name for analysis
elev <- arrange(elev, Camera.Trap.Name)
head(elev)

##slope
slope <- read.csv(here("data", "slope.csv"))
names(slope)[1] <- "Camera.Trap.Name"
names(slope)[2] <- "slope_buff250_mean"
names(slope)[3] <- "slope_buff250_stdev"
slope$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", slope$Camera.Trap.Name)
slope$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", slope$Camera.Trap.Name)

#Sort by camera trap unit name for analysis
slope <- arrange(slope, Camera.Trap.Name)
head(slope)

# tree structure ## excluidas estações amostrais da area queimada
#trees <- read.csv(here("data", "tree_structure.csv"))
trees <- read.csv(here("data", "trees.csv"))

#Sort by camera trap unit name for analysis
trees <- arrange(trees, Camera.Trap.Name)
head(trees)

##biomass prey prior >15kg /CT/dia##
#biomass_pp <- readRDS(here("data", "cov_rai_pp_2017.rds"))
biomass_pp <- readRDS(here("data", "cov_biomass_pp_2017.rds"))
names(biomass_pp)[1] <- "Camera.Trap.Name"
names(biomass_pp)[3] <- "photoRates_biomass_pp_17"
biomass_pp <- arrange(biomass_pp, Camera.Trap.Name)
head(biomass_pp)

# biomass small prey <15kg /CT/dia

biomass_sp <- readRDS(here("data", "cov_biomass_sp_2017.rds"))
names(biomass_sp)[1] <- "Camera.Trap.Name"
names(biomass_sp)[3] <- "photoRates_biomass_sp_17"

#Sort by camera trap unit name for analysis

biomass_sp <- arrange(biomass_sp, Camera.Trap.Name)
head(biomass_sp)

##blocos - sul e norte ##covariavel incluida devido a distancia e diferenças de relevo observada entre os
##grids de armadilhas fotograficas sul e norte
block <- read.csv(here("data", "blocos.csv"))
names(block)[1] <- "Camera.Trap.Name"
names(block)[2] <- "block" ## 
block$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", block$Camera.Trap.Name)
block$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", block$Camera.Trap.Name)

#Sort by camera trap unit name for analysis
block <- arrange(block, Camera.Trap.Name)
head(block)

##biomass prey prior >15kg /CT/dia##
#biomass_pp <- readRDS(here("data", "cov_rai_pp_2017.rds"))
biomass_pp_18 <- readRDS(here("data", "cov_biomass_pp_2018.rds"))
names(biomass_pp_18)[1] <- "Camera.Trap.Name"
names(biomass_pp_18)[3] <- "photoRates_biomass_pp_18"
biomass_pp_18 <- arrange(biomass_pp_18, Camera.Trap.Name)
head(biomass_pp_18)

# biomass small prey <15kg /CT/dia

biomass_sp_18 <- readRDS(here("data", "cov_biomass_sp_2018.rds"))
names(biomass_sp_18)[1] <- "Camera.Trap.Name"
names(biomass_sp_18)[3] <- "photoRates_biomass_sp_18"
#Sort by camera trap unit name for analysis

biomass_sp_18 <- arrange(biomass_sp_18, Camera.Trap.Name)
head(biomass_sp_18)


## create a single covariates dataframe 
covars <- merge(dist.water[,c(1,4)], dist.edge[,1:2], by="Camera.Trap.Name", all.x = T) ##distagua sarivers
covars <- merge(covars, elev[,c(1,4)], by="Camera.Trap.Name", all.x = T) ##sd_buf250
covars <- merge(covars, slope[,c(1,3)], by="Camera.Trap.Name", all.x = T)##st_dev buf250
covars <- merge(covars, trees[,c(1,4,3)], by="Camera.Trap.Name", all.x = T)
covars <- merge(covars, biomass_pp[,c(1,3)], by="Camera.Trap.Name", all.x = T) 
covars <- merge(covars, biomass_sp[,c(1,3)], by="Camera.Trap.Name", all.x = T)
covars <- merge(covars, biomass_pp_18[,c(1,3)], by="Camera.Trap.Name", all.x = T) 
covars <- merge(covars, biomass_sp_18[,c(1,3)], by="Camera.Trap.Name", all.x = T)
covars <- merge(covars, fire[,c(1,2)], by="Camera.Trap.Name", all.x = T)
covars <- merge(covars, cover[,c(1,2)], by="Camera.Trap.Name", all.x = T)
covars <- merge(covars, block[,c(1,2)], by="Camera.Trap.Name", all.x = T)
covars



## salvando covars para analise exploratoria
##write.csv( x = covars, file = here("results", "covars_Ponca.csv"), row.names = F)

# create a dataset for single-season model

Ponca17_18 <- dataRBG_species_Panthera_onca

##juntando dados de presença e ausencia com as covariaveis
Ponca_17_18 <- merge(Ponca17_18, covars, by="Camera.Trap.Name")
Ponca_17_18

# Save to disk

saveRDS(Ponca_17_18, here("data","Ponca_17_18_team.rds"))

