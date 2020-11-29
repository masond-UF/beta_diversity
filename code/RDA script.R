# 26 November 2020, David Mason—Final project ####
library(tidyverse)
library(vegan)
install.packages("betapart")
spec <- read.csv("data/species.csv")
str(spec)
summary(spec)
spec <- spec[, colSums(spec != 0) > 0] # drop columns without observtions


env <- read.csv("data/env.csv")
str(env)
summary(env)
env <- env[,3:96] # drop julian date and site
env <- env[,-c(12,13,14)] # drop soil texture proportions
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
									 	MG + PH + GRASS, env)
summary(red_spe.rda)
# 31% of the variance?

