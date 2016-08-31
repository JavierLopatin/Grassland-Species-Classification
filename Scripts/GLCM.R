

### GLCM

library(raster)
library(rgdal)
library(glcm)
library(doParallel)

# set working direction
#home = "/media/javier/JavierLopatin/Sp_Images/Site1/Raw/"
home = "F:/Sp_Images/Site4/Raw"
setwd(home)

processingFolder = "BN_MNF"
bandName = paste("MNF_", seq(1,10,1), sep="")

rasterNames <- rasterListNames(fileExtantion = ".tif", folder = processingFolder)
rasterlist <- rasterList(fileExtantion = ".tif", folder = processingFolder, rasterNames = bandName)

# MNF band to apply the transformation
band = 2

# create directory to stor GLCM images
dir.create(paste(processingFolder, "_GLCM", sep=""))

# initialize parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# run GLCM
for (i in 1:length(rasterlist)){
  rast <- rasterlist[[i]]
  rasterName = rasterNames[[i]]
  output <- GLCM(rast[[band]])
  out = paste( file.path(home, paste(processingFolder, "_GLCM", sep="")), "/", rasterName, "_GLCM.tif", sep="")
  writeRaster(output, filename=out, options="INTERLEAVE=BAND", overwrite=TRUE)
}
# stop parallel process
stopCluster(cl)
