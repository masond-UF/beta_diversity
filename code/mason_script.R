# start 26 November 2020 / update 13 February 2021,
# David Mason—Beta diversity analysis ####
library(tidyverse)
library(vegan)
library(ade4)
library(packfor)
library(spdep)
library(betapart)
library(spdep)
library(corrplot)

options(scipen = 999)

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

# separate the categorical and ordinal data
env_factors <- select(env, PERIOD, TEXTURE, OWNERSHIP, DISTURB, INTENSITY)
# separate the numerical data
env_numeric <- select(env, -PERIOD, -TEXTURE, -OWNERSHIP, -DISTURB, -INTENSITY)
# scale the data
env_numeric <- as.data.frame(scale(env_numeric))
# bring them back together 
env <- cbind(env_factors, env_numeric)
# rearrange the data so types (e.g., soil, landscape, overstory) are grouped
env <- select(env, PERIOD, OWNERSHIP, DISTURB, INTENSITY, TEXTURE, everything())
# convert intensity into ordinal 
env$INTENSITY <- factor(env$INTENSITY, order = TRUE, 
                       levels = c("NONE", "LOW", "MED", "HIGH"))

# Check for collinearity ####
# Correlations among continuous variables
env.cont.corr <- env %>% 
	select() %>% 
	cor(method = c("kendall"))
corrplot(env.cont.corr, tl.cex = 0.5)

env.cont.corr.df <- as.data.frame(as.table(env.cont.corr))
env.cont.corr.filt <- env.cont.corr.df %>%  
											arrange(desc(Freq)) %>% 
											filter(Freq>0.7)

corr.env.cont.drop.list <- c("Morus.rubra", "Cercis.canadensis", "Ilex.decidua",
														 "S", "Asimina.triloba", "Carya.pallida")

# Correlations among categorical and continuous variables

# DIST
env.nominal.corr <- env %>% 
	select(-PERIOD, -TEXTURE, -OWNERSHIP, -INTENSITY) %>% 
	select(DISTURB, everything())

DIST <- env.nominal.corr[,1]
env.nominal.corr <- env.nominal.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = env.nominal.corr[,i] ~ DIST, 
					data = env.nominal.corr)
	rsq <- summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(env.nominal.corr)[[2]]) 
# No variables correlated with disturbance type (DIST) > 0.7

# PERIOD
env.nominal.corr <- env %>% 
	select(-DISTURB, -TEXTURE, -OWNERSHIP, -INTENSITY) %>% 
	select(PERIOD, everything())

PERIOD <- env.nominal.corr[,1]
env.nominal.corr <- env.nominal.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = env.nominal.corr[,i] ~ PERIOD, 
					data = env.nominal.corr)
	rsq <-  summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(env.nominal.corr)[[2]]) 
# soilNA correlated with sampling period >0.7

# TEXTURE
env.nominal.corr <- env %>% 
	select(-PERIOD, -DISTURB, -OWNERSHIP, -INTENSITY) %>% 
	select(TEXTURE, everything())

TEXTURE <- env.nominal.corr[,1]
env.nominal.corr <- env.nominal.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = env.nominal.corr[,i] ~ TEXTURE, 
					data = env.nominal.corr)
	rsq <-  summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(env.nominal.corr)[[2]]) 
# Quercus.pagoda correlated with soil texture >0.7

# OWNERSHIP
env.nominal.corr <- env %>% 
	select(-PERIOD, -DISTURB, -TEXTURE, -INTENSITY) %>% 
	select(OWNERSHIP, everything())

OWNERSHIP <- env.nominal.corr[,1]
env.nominal.corr <- env.nominal.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = env.nominal.corr[,i] ~ OWNERSHIP, 
					data = env.nominal.corr)
	rsq <-  summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(env.nominal.corr)[[2]]) 
# DEVHI correlated with land ownership  >0.7

# INTENSITY
env.ordinal.corr <- env %>% 
	select(-PERIOD, -DISTURB, -TEXTURE, -OWNERSHIP) %>% 
	select(INTENSITY, everything())

INTENSITY <- env.ordinal.corr[,1]
env.ordinal.corr <- env.ordinal.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = env.ordinal.corr[,i] ~ INTENSITY, 
					data = env.ordinal.corr)
	rsq <-  summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(env.nominal.corr)[[2]]) # nothing
# No variables correlated with disturbance intensity >0.7

# List of species correlated with categorical variables
corr.env.cat.drop.list <- c("soilNA", "Quercus.pagoda", "DEVHI")

# Combined list of correlated variables to drop
corr.env.comb.drop.list <- c(corr.env.cont.drop.list, corr.env.cat.drop.list)
# Drop the correlated variables
'%!in%' <- function(x,y){ 
	!('%in%'(x,y)) # make a function to do the opposite of %in%
}

env <- env[which(names(env) %!in% corr.env.comb.drop.list)]																

# Variable selection ####
# Full species matrix = 13% explained
spe.rda <- rda(spe.hel ~ ., env)
full_step.forward <- ordiR2step(rda(spe.hel ~ 1, data=env), 
scope=formula(spe.rda), R2scope = F, direction="forward", pstep=1000)

# Reduced species matrix = 16% explained
red_spe.rda <- rda(red_spe.hel ~ ., env)
red_step.forward <- ordiR2step(rda(red_spe.hel ~ 1, data=env), 
scope=formula(red_spe.rda), R2scope = F, direction="forward", pstep=1000)
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

# extract threshold used to create PCNM
dmin <- pcnm_obj$thresh

# which eigenvalues show positive spatial autoc
length(which(pcnm_obj$values>0.0000001))

# Add the spatial component onto the env matrix
env <- cbind(env, pcnm_vec)

# Variance partitioning ####
red_spe.part <- varpart(red_spe.hel, ~TEXTURE+MG+PH+ROCK, 
												~CANOPY+Quercus.alba+Pinus.taeda,
												~DECID+DEVOP+CROP+WATER.150+GRASS+OWNERSHIP,
												~PERIOD + PCNM1 + PCNM2 + PCNM3 + 
												 PCNM4 + PCNM5 + PCNM6 + PCNM7 +
												 PCNM8 + PCNM9 + PCNM10 + PCNM11 +
												 PCNM12 + PCNM13 + PCNM14 + 
												 PCNM15 + PCNM16 + PCNM17 + 
												 PCNM18 + PCNM19 + PCNM20 +
												 PCNM21 + PCNM22 + PCNM23 +
												 PCNM24 + PCNM25 + PCNM26 +
												 PCNM27 + PCNM28 + PCNM29 +
												 PCNM30 + PCNM31, data = env)
plot(red_spe.part, digits = 2, Xnames = c('Soil & ground', 'Overstory', 
		 'Landscape', 'Space & time'), id.size = 1, bg = 2:5)

# Plot the variance paritioning
# Run anova on fractions ####
rda.all <- rda(red_spe.hel ~ TEXTURE + MG + PH + ROCK +
					 CANOPY + Quercus.alba + Pinus.taeda + DECID + 
					 DEVOP + CROP + WATER.150 + GRASS + OWNERSHIP +
					 PERIOD + PCNM1 + PCNM2 + PCNM3 + 
					 PCNM4 + PCNM5 + PCNM6 + PCNM7 +
					 PCNM8 + PCNM9 + PCNM10 + PCNM11 +
					 PCNM12 + PCNM13 + PCNM14 + 
					 PCNM15 + PCNM16 + PCNM17 + 
					 PCNM18 + PCNM19 + PCNM20 +
					 PCNM21 + PCNM22 + PCNM23 +
					 PCNM24 + PCNM25 + PCNM26 +
					 PCNM27 + PCNM28 + PCNM29 +
					 PCNM30 + PCNM31, data = env)
anova(rda.all)

rda.soil.ground <- rda(red_spe.hel ~ TEXTURE + 
											 	MG + PH + ROCK, data = env)
anova(rda.soil.ground)

rda.overstory <- rda(red_spe.hel ~ CANOPY + Quercus.alba + 
					 					Pinus.taeda, data = env)
anova(rda.overstory)

rda.landscape <- rda(red_spe.hel ~ DECID + DEVOP + CROP + 
								WATER.150 + GRASS + OWNERSHIP, data = env)
anova(rda.landscape)

rda.space.time <- rda(red_spe.hel ~ 
					 PERIOD + PCNM1 + PCNM2 + PCNM3 + 
					 PCNM4 + PCNM5 + PCNM6 + PCNM7 +
					 PCNM8 + PCNM9 + PCNM10 + PCNM11 +
					 PCNM12 + PCNM13 + PCNM14 + 
					 PCNM15 + PCNM16 + PCNM17 + 
					 PCNM18 + PCNM19 + PCNM20 +
					 PCNM21 + PCNM22 + PCNM23 +
					 PCNM24 + PCNM25 + PCNM26 +
					 PCNM27 + PCNM28 + PCNM29 +
					 PCNM30 + PCNM31, data = env)
anova(rda.space.time)


