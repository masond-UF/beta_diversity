# 26 November 2020, David Mason—Final project ####
library(tidyverse)
library(vegan)
spec <- read.csv("data/species.csv")
str(spec)
summary(spec)
spec <- spec[, colSums(spec != 0) > 0] # drop columns without observtions


env <- read.csv("data/env.csv")
str(env)
summary(env)
# Select species ####
source("code/biostats.r") # attach the code from biostats

occur <- foa.plots(spec)
rare <- which(occur[,2]<5)
common <- which(occur[,2]>95)

red_spec <- spec[,-c(rare,common)]

# Transform species ####
# Transform environmental ####
# Add sampling period
# Drop some variables
# Make sure scaled things are indeed scaled
