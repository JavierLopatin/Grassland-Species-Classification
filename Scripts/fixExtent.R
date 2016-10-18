
home = "F:/Sp_Images/Site4/Raw"

setwd(home)

library(raster)
library(rgdal)

folder ="BN_MNF"

#### Source Functions from GitHub
source_github <- function(u) {
  # load package
  require(RCurl)
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script), envir=.GlobalEnv)
  detach("package:RCurl", unload=TRUE)
} 
source_github("https://raw.githubusercontent.com/JavierLopatin/Herbaceous-Species-Classification/master/Scripts/Functions.R")

# reference list
nameList_ref <- rasterListNames(".tif", ".")
raterRefList_ref <- rasterList(".tif", ".")

# target list
nameList_tar <- rasterListNames(".tif", folder)
raterRefList_tar <- rasterList(".tif", folder)

# check the names
data.frame(nameList_ref, nameList_tar)

# appy and save
for (i in 1:length(nameList_ref)){
  crs(raterRefList_tar[[i]]) <- crs(raterRefList_ref[[i]])  
  extent(raterRefList_tar[[i]]) <- extent(raterRefList_ref[[i]])  
  output = paste( nameList_tar[i], ".tif", sep="")
  writeRaster(raterRefList_tar[[i]], filename=output, options="INTERLEAVE=BAND", overwrite=TRUE)
}

