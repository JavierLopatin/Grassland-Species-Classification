## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

#### run Clasification!

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"

setwd(home)

load("results_canopy/fit_potVal.RData")
load("results_canopy/fit_rf.RData")
load("results_canopy/fit_potVal_BN.RData")
load("results_canopy/fit_potVal.RData")
load("results_canopy/fit_rf_BN.RData")
load("results_canopy/fit_potVal_MNF.RData")
load("results_canopy/fit_rf_MNF.RData")
load("results_canopy/fit_potVal_BN_MNF.RData")
load("results_canopy/fit_rf_BN_MNF.RData")
load("results_canopy/fit_potVal_GLCM.RData")

#####################
### load the data ###
#####################
# Site1
potVal <- read.table("data/potVal_all.csv", sep = ",", header = T)
rf     <- read.table("data/rf_all.csv", sep = ",", header = T)

potVal_BN <- read.table("data/potVal_BN_all.csv", sep = ",", header = T)
rf_BN     <- read.table("data/rf_BN_all.csv", sep = ",", header = T)

potVal_MNF <- read.table("data/potVal_MNF_all.csv", sep = ",", header = T)
rf_MNF     <- read.table("data/rf_MNF_all.csv", sep = ",", header = T)

potVal_BN_MNF <- read.table("data/potVal_BN_MNF_all.csv", sep = ",", header = T)
rf_BN_MNF     <- read.table("data/rf_BN_MNF_all.csv", sep = ",", header = T)

potVal_GLCM <- read.table("data/potVal_GLCM_all.csv", sep = ",", header = T)
rf_GLCM     <- read.table("data/rf_GLCM_all.csv", sep = ",", header = T)

potVal_BN_GLCM <- read.table("data/potVal_BN_GLCM_all.csv", sep = ",", header = T)
rf_BN_GLCM     <- read.table("data/rf_BN_GLCM_all.csv", sep = ",", header = T)

# wavelength
wl <- c( 398, 407, 415, 424, 432, 441, 450, 459, 468, 477, 486, 495, 504, 513, 522, 531, 540, 550, 558, 568,
         577, 587, 596, 605, 615, 624, 633, 643, 652, 661, 671, 680, 690, 699, 708, 717, 727, 736, 746, 755,
         765, 775, 784, 794, 803, 813, 822, 832, 842, 851, 861, 870, 880, 890, 899, 908, 918, 928, 937, 947, 957)
  
# load species cover dataset
#species <- read.table("Plots_Species.csv", header = T, sep=",")


#### Source Functions from GitHub
source_github <- function(u) {
  # load package
  require(RCurl)
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script), envir=.GlobalEnv)
  detach("package:RCurl", unload=TRUE)
} 
source_github("https://raw.githubusercontent.com/JavierLopatin/Herbaceous-Species-Classification/master/Scripts/0_Functions.R")

##########################
### Run Classification ###
##########################

#------------------------#
# Spectra
#------------------------#
# potVal
fit_potVal <- classificationEnsemble(potVal$Species, potVal[,3:length(potVal)], wl)
plot.classificationEnsemble( potVal[,3:length(potVal)]/10000, fit_potVal, xlab_tag=expression(lambda(nm)) )
save(fit_potVal, file="results_canopy/fit_potVal.RData")

# rip it off
fit_rf <- classificationEnsemble(rf$Species, rf[,3:length(rf)], wl)
plot.classificationEnsemble( rf[,3:length(rf)]/10000, fit_rf, xlab_tag=expression(lambda(nm)) )
save(fit_rf, file="results_canopy/fit_rf.RData")


#------------------------#
# Spectra BN
#------------------------#
# potVal
fit_potVal_BN <- classificationEnsemble(potVal_BN$Species, potVal_BN[,3:length(potVal_BN)], wl)
plot.classificationEnsemble( potVal_BN[,3:length(potVal_BN)], fit_potVal_BN, xlab_tag=expression(lambda(nm)) )
save(fit_potVal_BN, file="results_canopy/fit_potVal_BN.RData")

# rip it off
fit_rf_BN <- classificationEnsemble(rf_BN$Species, rf_BN[,3:length(rf_BN)], wl)
plot.classificationEnsemble( rf_BN[,3:length(rf_BN)], fit_rf_BN, xlab_tag=expression(lambda(nm)) )
save(fit_rf_BN, file="results_canopy/fit_rf_BN.RData")

#------------------------#
# MNF
#------------------------#
# potVal
fit_potVal_MNF <- classificationEnsemble(potVal_MNF$Species, potVal_MNF[,3:length(potVal_MNF)], seq(1,10,1))
plot.classificationEnsemble( potVal_MNF[,3:length(potVal_MNF)], fit_potVal_MNF, xlab_tag=expression(lambda(nm)) )
save(fit_potVal_MNF, file="results_canopy/fit_potVal_MNF.RData")

# rip it off
fit_rf_MNF <- classificationEnsemble(rf_MNF$Species, rf_MNF[,3:length(rf_MNF)], seq(1,10,1))
plot.classificationEnsemble( rf_MNF[,3:length(rf_MNF)], fit_rf_MNF, xlab_tag=expression(lambda(nm)) )
save(fit_rf_MNF, file="results_canopy/fit_rf_MNF.RData")

#------------------------#
# MNF BN
#------------------------#
# potVal
fit_potVal_BN_MNF <- classificationEnsemble(potVal_BN_MNF$Species, potVal_BN_MNF[,3:length(potVal_BN_MNF)], seq(1,10,1))
plot.classificationEnsemble( potVal_BN_MNF[,3:length(potVal_BN_MNF)], fit_potVal_BN_MNF, xlab_tag=expression(lambda(nm)) )
save(fit_potVal_BN_MNF, file="results_canopy/fit_potVal_BN_MNF.RData")

# rip it off
fit_rf_BN_MNF <- classificationEnsemble(rf_BN_MNF$Species, rf_BN_MNF[,3:length(rf_BN_MNF)], seq(1,10,1))
plot.classificationEnsemble( rf_BN_MNF[,3:length(rf_BN_MNF)], fit_rf_BN_MNF, xlab_tag=expression(lambda(nm)) )
save(fit_rf_BN_MNF, file="results_canopy/fit_rf_BN_MNF.RData")

#------------------------#
# GLCM
#------------------------#
# potVal
fit_potVal_GLCM <- classificationEnsemble(potVal_GLCM$Species, potVal_GLCM[,3:length(potVal_GLCM)], seq(1,6,1))
save(fit_potVal_GLCM, file="results_canopy/fit_potVal_GLCM.RData")

# rip it off
fit_rf_GLCM <- classificationEnsemble(rf_GLCM$Species, rf_GLCM[,3:length(rf_GLCM)], seq(1,6,1))
save(fit_rf_GLCM, file="results_canopy/fit_rf_GLCM.RData")

#------------------------#
# GLCM BN
#------------------------#
# potVal
fit_potVal_BN_GLCM <- classificationEnsemble(potVal_BN_GLCM$Species, potVal_BN_GLCM[,3:length(potVal_BN_GLCM)], seq(1,6,1))
save(fit_potVal_BN_GLCM, file="results_canopy/fit_potVal_BN_GLCM.RData")

# rip it off
fit_rf_BN_GLCM <- classificationEnsemble(rf_BN_GLCM$Species, rf_BN_GLCM[,3:length(rf_BN_GLCM)], seq(1,6,1))
save(fit_rf_BN_GLCM, file="results_canopy/fit_rf_BN_GLCM.RData")

# best model
fit_potVal$fit
fit_rf$fit
fit_potVal_BN$fit
fit_rf_BN$fit
fit_potVal_MNF$fit
fit_rf_MNF$fit
fit_potVal_BN_MNF$fit
fit_rf_BN_MNF$fit
fit_potVal_GLCM$fit
fit_rf_GLCM$fit
fit_potVal_BN_GLCM$fit
fit_rf_BN_GLCM$fit

# choose model
bestModel <- fit_rf_BN
bestData <- rf_BN

#####################################################
### Apply best model to the four sites separately ###
#####################################################

# obtain training data for the site
# Site1 
x = grep(1, bestData$Site) 
x1 <- bestData[x, ]
x1$Species <- factor(x1$Species)
# Site2 
x = grep(2, bestData$Site) 
x2 <- bestData[x, ]
x2$Species <- factor(x2$Species)
# Site3
x = grep(3, bestData$Site) 
x3 <- bestData[x, ]
x3$Species <- factor(x3$Species)
# Site4
x = grep(4, bestData$Site) 
x4 <- bestData[x, ]
x4$Species <- factor(x4$Species)


fit_site1 <- classificationEnsemble(x1$Species, x1[,3:length(x1)], wl)
plot.classificationEnsemble( x1[,3:length(x1)], fit_site1, xlab_tag=expression(lambda(nm)) )
fit_site1$fit
save(fit_site1, file="results_canopy/fit_site1.RData")

fit_site2 <- classificationEnsemble(x2$Species, x2[,3:length(x2)], wl)
fit_site2$fit
save(fit_site2, file="results_canopy/fit_site2.RData")

fit_site3 <- classificationEnsemble(x3$Species, x3[,3:length(x3)], wl)
fit_site3$fit
save(fit_site3, file="results_canopy/fit_site3.RData")

fit_site4 <- classificationEnsemble(x4$Species, x4[,3:length(x4)], wl)
fit_site4$fit
save(fit_site4, file="results_canopy/fit_site4.RData")


##################################
### predict model to the plots ###
##################################

library(raster)


# load images 
rasterDir = "F:/Sp_Images"
setwd(rasterDir)

# load rasters with "plot" pattern in the name
plots1 <-  rasterList(fileExtantion = ".tif", folder = "Site1/Raw/BN", dir=rasterDir, select="plot")
plots2 <-  rasterList(fileExtantion = ".tif", folder = "Site2/Raw/BN", dir=rasterDir, select="plot")
plots3 <-  rasterList(fileExtantion = ".tif", folder = "Site3/Raw/BN", dir=rasterDir, select="plot")
plots4 <-  rasterList(fileExtantion = ".tif", folder = "Site4/Raw/BN", dir=rasterDir, select="plot")



