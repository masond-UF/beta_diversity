# 26 November 2020, David Mason—Final project ####
library(tidyverse)
library(vegan)
library(ade4)
library(packfor)
library(spdep)
library(betapart)
install.packages('ape')

spec <- read.csv("data/species.csv")
str(spec)
summary(spec)
spec <- spec[, colSums(spec != 0) > 0] # drop columns without observtions


env <- read.csv("data/env.csv")
str(env)
summary(env)
spatial <- env[,2:3] # extract spatial variables for later analysis
env <- env[,-c(1:4,16,17,18)] # drop site, date, spatial and texture proportions
env <- env[, colSums(env != 0) > 0] # drop columns without values
# Select species ####
source("code/biostats.r") # attach the code from biostats

occur <- foa.plots(spec)
rare <- which(occur[,2]<5)
common <- which(occur[,2]>95)

red_spec <- spec[,-c(rare,common)]

# Transform data ####
spe.hel <- decostand(spec, "hellinger") 
red_spe.hel <- decostand(red_spec, "hellinger")

# separate the categorical data
env_factors <- select(env, PERIOD, TEXTURE, OWNERSHIP, DISTURB, INTENSITY)
# separate the numerical data
env_numeric <- select(env, -PERIOD, -TEXTURE, -OWNERSHIP, -DISTURB, -INTENSITY)
# scale the data
env_numeric <- as.data.frame(scale(env_numeric))
# bring them back together 
env <- cbind(env_factors, env_numeric)
# rearrange the data so types (e.g., soil, landscape, overstory) are grouped
env <- select(env, PERIOD, OWNERSHIP, DISTURB, INTENSITY, TEXTURE, everything())
# Variable selection ####
# Full species matrix = 13% explained
spe.rda <- rda(spe.hel ~ ., env)
full_step.forward <- ordiR2step(rda(spe.hel ~ 1, data=env), 
scope=formula(spe.rda), R2scope = F, direction="forward", pstep=1000)

# Reduced species matrix = 16% explained
red_spe.rda <- rda(red_spe.hel ~ ., env)
red_step.forward <- ordiR2step(rda(red_spe.hel ~ 1, data=env), 
scope=formula(red_spe.rda), R2scope = F, direction="forward", pstep=1000)
# Conduct RDA on most parsimonious model ####
red_spe.rda <- rda(red_spe.hel ~ TEXTURE + OWNERSHIP + CANOPY + Quercus.alba + 
									 	DECID + DEVOP + Pinus.taeda + PERIOD + CROP + WATER.150 + 
									 	MG + PH + GRASS + ROCK, env)
red_spe.rda
# 32% of the variance?
# Create the PCNM variables ####
# normalize and center the coordinate data around (0,0)
spatial <- as.data.frame(scale(spatial, center=TRUE, scale=FALSE))
# create distance matrix
coord_dist <- dist(spatial) 
# create PCNM object
pcnm_obj <- pcnm(coord_dist)
# extract variables
pcnm_vec <- pcnm_obj$vectors
# visualize PCNM variables
plot(spatial, pch=15, cex=pcnm_vec[,3]*8+4) 
# extract threshold used to create PCNM. The
# distance matrix is truncated by this threshold,
# after which distances are filled with arbitrary high values
dmin <- pcnm_obj$thresh
# extract the number of positive EV
posEV <- length(pcnm_obj$values)
# extract PCNM variables with positive spatial autocorrelation
pcnm_obj

select <- which(pcnm_obj$Moran_I$Positive==TRUE) 
pcnm3=pcnm2[,select]

# Variance partitioning ####
red_spe.part <- varpart(red_spe.hel, ~TEXTURE+MG+PH+ROCK, 
												~CANOPY+Quercus.alba+Pinus.taeda,
												~DECID+DEVOP+CROP+WATER.150+GRASS+OWNERSHIP,
												~PERIOD+X+Y,
												data = env)
plot(red_spe.part,
		 digits = 2,
		 Xnames = c('Soil & ground', 'Overstory', 
		 					 'Landscape', 'Space & time'),
		 id.size = 1, bg = 2:5)

# Plot the variance paritioning

# Run anova on fractions ####
# Explore environmental changes between periods ####
