

### GLCM

library(raster)
library(rgdal)
library(glcm)
library(doParallel)

# set working direction
home = "F:/Sp_Images/Site4/Raw"
setwd(home)

bandName = paste("MNF_", seq(1,10,1), sep="")

rasterNames <- rasterListNames(fileExtantion = ".tif", folder = "MNF")
rasterlist <- rasterList(fileExtantion = ".tif", folder = "MNF", rasterNames = bandName)

# MNF band to apply the transformation
band = 2

# create directory to stor GLCM images
dir.create("GLCM")

# initialize parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# run GLCM
for (i in 1:length(rasterlist)){
  rast <- rasterlist[[i]]
  rasterName = rasterNames[[i]]
  output <- GLCM(rast[[band]])
  out = paste( file.path(home, "GLCM"), "/", rasterName, "_GLCM.tif", sep="")
  writeRaster(output, filename=out, options="INTERLEAVE=BAND", overwrite=TRUE)
}
# stop parallel process
stopCluster(cl)
