# original attempt at spatial

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
