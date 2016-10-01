## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

#### run Clasification!

library(raster)

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"

setwd(home)

load("results_canopy/fit_potVal_1.RData")
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


# load images 
rasterDir = "F:/Sp_Images"
setwd(rasterDir)

# load rasters with "plot" pattern in the name
plots1 <-  rasterList(fileExtantion = ".tif", folder = "Site1/Raw/BN", dir=rasterDir, select="plot")
plots2 <-  rasterList(fileExtantion = ".tif", folder = "Site2/Raw/BN", dir=rasterDir, select="plot")
plots3 <-  rasterList(fileExtantion = ".tif", folder = "Site3/Raw/BN", dir=rasterDir, select="plot")
plots4 <-  rasterList(fileExtantion = ".tif", folder = "Site4/Raw/BN", dir=rasterDir, select="plot")


##########################
### Run Classification ###
##########################

setwd(home)

outputDir = "D:/Sp_Images"

#------------------------#
# Spectra                #
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_1 <- tunningModels(data = potVal, Site = 1, wl)
save(fit_potVal_1, file="results_canopy/fit_potVal_1.RData")
# Bootstrap validation 
fit_boot_potVal_1 <- ApplyBootsClassification(data = potVal, Site = 1, rasterPlots = plots1,
                               en = fit_potVal_1, outDir = outputDir, modelTag = "potVal" )
save(fit_boot_potVal_1, file="results_canopy/fit_boot_potVal_1.RData")

#### Site 2
fit_potVal_2 <- tunningModels(data = potVal, Site = 2, wl)
save(fit_potVal_2, file="results_canopy/fit_potVal_2.RData")
# Bootstrap validation 
fit_boot_potVal_2 <- ApplyBootsClassification(data = potVal, Site = 2, rasterPlots = plots2, 
                                            en = fit_potVal_2, outDir = outputDir, modelTag = "potVal" )
save(fit_boot_potVal_2, file="results_canopy/fit_boot_potVal_2.RData")

#### Site 3
fit_potVal_3 <- tunningModels(data = potVal, Site = 3, wl)
save(fit_potVal_3, file="results_canopy/fit_potVal_3.RData")
# Bootstrap validation 
fit_boot_potVal_3 <- ApplyBootsClassification(data = potVal, Site = 3, rasterPlots = plots3, 
                                            en = fit_potVal_3, outDir = outputDir, modelTag = "potVal" )
save(fit_boot_potVal_3, file="results_canopy/fit_boot_potVal_3.RData")

#### Site 4
fit_potVal_4 <- tunningModels(data = potVal, Site = 4, wl)
save(fit_potVal_4, file="results_canopy/fit_potVal_4.RData")
# Bootstrap validation 
fit_boot_potVal_4 <- ApplyBootsClassification(data = potVal, Site = 4, rasterPlots = plots4, 
                                            en = fit_potVal_4, outDir = outputDir, modelTag = "potVal" )
save(fit_boot_potVal_4, file="results_canopy/fit_boot_potVal_4.RData")

##################
### rip it off ###
##################

#### Site 1
fit_rf_1 <- tunningModels(data = rf, Site = 1, wl)
save(fit_rf_1, file="results_canopy/fit_rf_1.RData")
# Bootstrap validation 
fit_boot_rf_1 <- ApplyBootsClassification(data = rf, Site = 1, rasterPlots = plots1, 
                                              en = fit_rf_1, outDir = outputDir, modelTag = "rf" )
save(fit_boot_rf_1, file="results_canopy/fit_boot_rf_1.RData")

#### Site 2
fit_rf_2 <- tunningModels(data = rf, Site = 2, wl)
save(fit_rf_2, file="results_canopy/fit_rf_2.RData")
# Bootstrap validation 
fit_boot_rf_2 <- ApplyBootsClassification(data = rf, Site = 2, rasterPlots = plots2, 
                                              en = fit_rf_2, outDir = outputDir, modelTag = "rf" )
save(fit_boot_rf_2, file="results_canopy/fit_boot_rf_2.RData")

#### Site 3
fit_rf_3 <- tunningModels(data = rf, Site = 3, wl)
save(fit_rf_3, file="results_canopy/fit_rf_3.RData")
# Bootstrap validation 
fit_boot_rf_3 <- ApplyBootsClassification(data = rf, Site = 3, rasterPlots = plots3, 
                                              en = fit_rf_3, outDir = outputDir, modelTag = "rf" )
save(fit_boot_rf_3, file="results_canopy/fit_boot_rf_3.RData")

#### Site 4
fit_rf_4 <- tunningModels(data = rf, Site = 4, wl)
save(fit_rf_4, file="results_canopy/fit_rf_4.RData")
# Bootstrap validation 
fit_boot_rf_4 <- ApplyBootsClassification(data = rf, Site = 4, rasterPlots = plots4, 
                                              en = fit_rf_4, outDir = outputDir, modelTag = "rf" )
save(fit_boot_rf_4, file="results_canopy/fit_boot_rf_4.RData")


#------------------------#
# Spectra BN
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_BN_1 <- tunningModels(data = potVal_BN, Site = 1, wl)
save(fit_potVal_BN_1, file="results_canopy/fit_potVal_BN_1.RData")
# Bootstrap validation 
fit_boot_potVal_BN_1 <- ApplyBootsClassification(data = potVal_BN, Site = 1, rasterPlots = plots1, 
                                              en = fit_potVal_BN_1, outDir = outputDir, modelTag = "potVal_BN" )
save(fit_boot_potVal_BN_1, file="results_canopy/fit_boot_potVal_BN_1.RData")

#### Site 2
fit_potVal_BN_2 <- tunningModels(data = potVal_BN, Site = 2, wl)
save(fit_potVal_BN_2, file="results_canopy/fit_potVal_BN_2.RData")
# Bootstrap validation 
fit_boot_potVal_BN_2 <- ApplyBootsClassification(data = potVal_BN, Site = 2, rasterPlots = plots2, 
                                              en = fit_potVal_BN_2, outDir = outputDir, modelTag = "potVal_BN" )
save(fit_boot_potVal_BN_2, file="results_canopy/fit_boot_potVal_BN_2.RData")

#### Site 3
fit_potVal_BN_3 <- tunningModels(data = potVal_BN, Site = 3, wl)
save(fit_potVal_BN_3, file="results_canopy/fit_potVal_BN_3.RData")
# Bootstrap validation 
fit_boot_potVal_BN_3 <- ApplyBootsClassification(data = potVal_BN, Site = 3, rasterPlots = plots3, 
                                              en = fit_potVal_BN_3, outDir = outputDir, modelTag = "potVal_BN" )
save(fit_boot_potVal_BN_3, file="results_canopy/fit_boot_potVal_BN_3.RData")

#### Site 4
fit_potVal_BN_4 <- tunningModels(data = potVal_BN, Site = 4, wl)
save(fit_potVal_BN_4, file="results_canopy/fit_potVal_BN_4.RData")
# Bootstrap validation 
fit_boot_potVal_BN_4 <- ApplyBootsClassification(data = potVal_BN, Site = 4, rasterPlots = plots4, 
                                              en = fit_potVal_BN_4, outDir = outputDir, modelTag = "potVal_BN" )
save(fit_boot_potVal_BN_4, file="results_canopy/fit_boot_potVal_BN_4.RData")

##################
### rip it off ###
##################

#### Site 1
fit_rf_BN_1 <- tunningModels(data = rf_BN, Site = 1, wl)
save(fit_rf_BN_1, file="results_canopy/fit_rf_BN_1.RData")
# Bootstrap validation 
fit_boot_rf_BN_1 <- ApplyBootsClassification(data = rf_BN, Site = 1, rasterPlots = plots1, 
                                          en = fit_rf_BN_1, outDir = outputDir, modelTag = "rf_BN" )
save(fit_boot_rf_BN_1, file="results_canopy/fit_boot_rf_BN_1.RData")

#### Site 2
fit_rf_BN_2 <- tunningModels(data = rf_BN, Site = 2, wl)
save(fit_rf_BN_2, file="results_canopy/fit_rf_BN_2.RData")
# Bootstrap validation 
fit_boot_rf_BN_2 <- ApplyBootsClassification(data = rf_BN, Site = 2, rasterPlots = plots2, 
                                          en = fit_rf_BN_2, outDir = outputDir, modelTag = "rf_BN" )
save(fit_boot_rf_BN_2, file="results_canopy/fit_boot_rf_BN_2.RData")

#### Site 3
fit_rf_BN_3 <- tunningModels(data = rf_BN, Site = 3, wl)
save(fit_rf_BN_3, file="results_canopy/fit_rf_BN_3.RData")
# Bootstrap validation 
fit_boot_rf_BN_3 <- ApplyBootsClassification(data = rf_BN, Site = 3, rasterPlots = plots3, 
                                          en = fit_rf_BN_3, outDir = outputDir, modelTag = "rf_BN" )
save(fit_boot_rf_BN_3, file="results_canopy/fit_boot_rf_BN_3.RData")

#### Site 4
fit_rf_BN_4 <- tunningModels(data = rf_BN, Site = 4, wl)
save(fit_rf_BN_4, file="results_canopy/fit_rf_BN_4.RData")
# Bootstrap validation 
fit_boot_rf_BN_4 <- ApplyBootsClassification(data = rf_BN, Site = 4, rasterPlots = plots4, 
                                          en = fit_rf_BN_4, outDir = outputDir, modelTag = "rf_BN" )
save(fit_boot_rf_BN_4, file="results_canopy/fit_boot_rf_BN_4.RData")

#------------------------#
# MNF
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_MNF_1 <- tunningModels(data = potVal_MNF, Site = 1, wl)
save(fit_potVal_MNF_1, file="results_canopy/fit_potVal_MNF_1.RData")
# Bootstrap validation 
fit_boot_potVal_MNF_1 <- ApplyBootsClassification(data = potVal_MNF, Site = 1, rasterPlots = plots1, 
                                                 en = fit_potVal_MNF_1, outDir = outputDir, modelTag = "potVal_MNF" )
save(fit_boot_potVal_MNF_1, file="results_canopy/fit_boot_potVal_MNF_1.RData")

#### Site 2
fit_potVal_MNF_2 <- tunningModels(data = potVal_MNF, Site = 2, wl)
save(fit_potVal_MNF_2, file="results_canopy/fit_potVal_MNF_2.RData")
# Bootstrap validation 
fit_boot_potVal_MNF_2 <- ApplyBootsClassification(data = potVal_MNF, Site = 2, rasterPlots = plots2, 
                                                 en = fit_potVal_MNF_2, outDir = outputDir, modelTag = "potVal_MNF" )
save(fit_boot_potVal_MNF_2, file="results_canopy/fit_boot_potVal_MNF_2.RData")

#### Site 3
fit_potVal_MNF_3 <- tunningModels(data = potVal_MNF, Site = 3, wl)
save(fit_potVal_MNF_3, file="results_canopy/fit_potVal_MNF_3.RData")
# Bootstrap validation 
fit_boot_potVal_MNF_3 <- ApplyBootsClassification(data = potVal_MNF, Site = 3, rasterPlots = plots3, 
                                                 en = fit_potVal_MNF_3, outDir = outputDir, modelTag = "potVal_MNF" )
save(fit_boot_potVal_MNF_3, file="results_canopy/fit_boot_potVal_MNF_3.RData")

#### Site 4
fit_potVal_MNF_4 <- tunningModels(data = potVal_MNF, Site = 4, wl)
save(fit_potVal_MNF_4, file="results_canopy/fit_potVal_MNF_4.RData")
# Bootstrap validation 
fit_boot_potVal_MNF_4 <- ApplyBootsClassification(data = potVal_MNF, Site = 4, rasterPlots = plots4, 
                                                 en = fit_potVal_MNF_4, outDir = outputDir, modelTag = "potVal_MNF" )
save(fit_boot_potVal_MNF_4, file="results_canopy/fit_boot_potVal_MNF_4.RData")

##################
### rip it off ###
##################

#### Site 1
fit_rf_MNF_1 <- tunningModels(data = rf_MNF, Site = 1, wl)
save(fit_rf_MNF_1, file="results_canopy/fit_rf_MNF_1.RData")
# Bootstrap validation 
fit_boot_rf_MNF_1 <- ApplyBootsClassification(data = rf_MNF, Site = 1, rasterPlots = plots1, 
                                             en = fit_rf_MNF_1, outDir = outputDir, modelTag = "rf_MNF" )
save(fit_boot_rf_MNF_1, file="results_canopy/fit_boot_rf_MNF_1.RData")

#### Site 2
fit_rf_MNF_2 <- tunningModels(data = rf_MNF, Site = 2, wl)
save(fit_rf_MNF_2, file="results_canopy/fit_rf_MNF_2.RData")
# Bootstrap validation 
fit_boot_rf_MNF_2 <- ApplyBootsClassification(data = rf_MNF, Site = 2, rasterPlots = plots2, 
                                             en = fit_rf_MNF_2, outDir = outputDir, modelTag = "rf_MNF" )
save(fit_boot_rf_MNF_2, file="results_canopy/fit_boot_rf_MNF_2.RData")

#### Site 3
fit_rf_MNF_3 <- tunningModels(data = rf_MNF, Site = 3, wl)
save(fit_rf_MNF_3, file="results_canopy/fit_rf_MNF_3.RData")
# Bootstrap validation 
fit_boot_rf_MNF_3 <- ApplyBootsClassification(data = rf_MNF, Site = 3, rasterPlots = plots3, 
                                             en = fit_rf_MNF_3, outDir = outputDir, modelTag = "rf_MNF" )
save(fit_boot_rf_MNF_3, file="results_canopy/fit_boot_rf_MNF_3.RData")

#### Site 4
fit_rf_MNF_4 <- tunningModels(data = rf_MNF, Site = 4, wl)
save(fit_rf_MNF_4, file="results_canopy/fit_rf_MNF_4.RData")
# Bootstrap validation 
fit_boot_rf_MNF_4 <- ApplyBootsClassification(data = rf_MNF, Site = 4, rasterPlots = plots4, 
                                             en = fit_rf_MNF_4, outDir = outputDir, modelTag = "rf_MNF" )
save(fit_boot_rf_MNF_4, file="results_canopy/fit_boot_rf_MNF_4.RData")


#------------------------#
# MNF BN
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_BN_MNF_1 <- tunningModels(data = potVal_BN_MNF, Site = 1, wl)
save(fit_potVal_BN_MNF_1, file="results_canopy/fit_potVal_BN_MNF_1.RData")
# Bootstrap validation 
fit_boot_potVal_BN_MNF_1 <- ApplyBootsClassification(data = potVal_BN_MNF, Site = 1, rasterPlots = plots1, 
                                                 en = fit_potVal_BN_MNF_1, outDir = outputDir, modelTag = "potVal_BN_MNF" )
save(fit_boot_potVal_BN_MNF_1, file="results_canopy/fit_boot_potVal_BN_MNF_1.RData")

#### Site 2
fit_potVal_BN_MNF_2 <- tunningModels(data = potVal_BN_MNF, Site = 2, wl)
save(fit_potVal_BN_MNF_2, file="results_canopy/fit_potVal_BN_MNF_2.RData")
# Bootstrap validation 
fit_boot_potVal_BN_MNF_2 <- ApplyBootsClassification(data = potVal_BN_MNF, Site = 2, rasterPlots = plots2, 
                                                 en = fit_potVal_BN_MNF_2, outDir = outputDir, modelTag = "potVal_BN_MNF" )
save(fit_boot_potVal_BN_MNF_2, file="results_canopy/fit_boot_potVal_BN_MNF_2.RData")

#### Site 3
fit_potVal_BN_MNF_3 <- tunningModels(data = potVal_BN_MNF, Site = 3, wl)
save(fit_potVal_BN_MNF_3, file="results_canopy/fit_potVal_BN_MNF_3.RData")
# Bootstrap validation 
fit_boot_potVal_BN_MNF_3 <- ApplyBootsClassification(data = potVal_BN_MNF, Site = 3, rasterPlots = plots3, 
                                                 en = fit_potVal_BN_MNF_3, outDir = outputDir, modelTag = "potVal_BN_MNF" )
save(fit_boot_potVal_BN_MNF_3, file="results_canopy/fit_boot_potVal_BN_MNF_3.RData")

#### Site 4
fit_potVal_BN_MNF_4 <- tunningModels(data = potVal_BN_MNF, Site = 4, wl)
save(fit_potVal_BN_MNF_4, file="results_canopy/fit_potVal_BN_MNF_4.RData")
# Bootstrap validation 
fit_boot_potVal_BN_MNF_4 <- ApplyBootsClassification(data = potVal_BN_MNF, Site = 4, rasterPlots = plots4, 
                                                 en = fit_potVal_BN_MNF_4, outDir = outputDir, modelTag = "potVal_BN_MNF" )
save(fit_boot_potVal_BN_MNF_4, file="results_canopy/fit_boot_potVal_BN_MNF_4.RData")

##################
### rip it off ###
##################

#### Site 1
fit_rf_BN_MNF_1 <- tunningModels(data = rf_BN_MNF, Site = 1, wl)
save(fit_rf_BN_MNF_1, file="results_canopy/fit_rf_BN_MNF_1.RData")
# Bootstrap validation 
fit_boot_rf_BN_MNF_1 <- ApplyBootsClassification(data = rf_BN_MNF, Site = 1, rasterPlots = plots1, 
                                             en = fit_rf_BN_MNF_1, outDir = outputDir, modelTag = "rf_BN_MNF" )
save(fit_boot_rf_BN_MNF_1, file="results_canopy/fit_boot_rf_BN_MNF_1.RData")

#### Site 2
fit_rf_BN_MNF_2 <- tunningModels(data = rf_BN_MNF, Site = 2, wl)
save(fit_rf_BN_MNF_2, file="results_canopy/fit_rf_BN_MNF_2.RData")
# Bootstrap validation 
fit_boot_rf_BN_MNF_2 <- ApplyBootsClassification(data = rf_BN_MNF, Site = 2, rasterPlots = plots2, 
                                             en = fit_rf_BN_MNF_2, outDir = outputDir, modelTag = "rf_BN_MNF" )
save(fit_boot_rf_BN_MNF_2, file="results_canopy/fit_boot_rf_BN_MNF_2.RData")

#### Site 3
fit_rf_BN_MNF_3 <- tunningModels(data = rf_BN_MNF, Site = 3, wl)
save(fit_rf_BN_MNF_3, file="results_canopy/fit_rf_BN_MNF_3.RData")
# Bootstrap validation 
fit_boot_rf_BN_MNF_3 <- ApplyBootsClassification(data = rf_BN_MNF, Site = 3, rasterPlots = plots3, 
                                             en = fit_rf_BN_MNF_3, outDir = outputDir, modelTag = "rf_BN_MNF" )
save(fit_boot_rf_BN_MNF_3, file="results_canopy/fit_boot_rf_BN_MNF_3.RData")

#### Site 4
fit_rf_BN_MNF_4 <- tunningModels(data = rf_BN_MNF, Site = 4, wl)
save(fit_rf_BN_MNF_4, file="results_canopy/fit_rf_BN_MNF_4.RData")
# Bootstrap validation 
fit_boot_rf_BN_MNF_4 <- ApplyBootsClassification(data = rf_BN_MNF, Site = 4, rasterPlots = plots4, 
                                             en = fit_rf_BN_MNF_4, outDir = outputDir, modelTag = "rf_BN_MNF" )
save(fit_boot_rf_BN_MNF_4, file="results_canopy/fit_boot_rf_BN_MNF_4.RData")

#------------------------#
# GLCM
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_GLCM_1 <- tunningModels(data = potVal_GLCM, Site = 1, wl)
save(fit_potVal_GLCM_1, file="results_canopy/fit_potVal_GLCM_1.RData")
# Bootstrap validation 
fit_boot_potVal_GLCM_1 <- ApplyBootsClassification(data = potVal_GLCM, Site = 1, rasterPlots = plots1, 
                                                 en = fit_potVal_GLCM_1, outDir = outputDir, modelTag = "potVal_GLCM" )
save(fit_boot_potVal_GLCM_1, file="results_canopy/fit_boot_potVal_GLCM_1.RData")

#### Site 2
fit_potVal_GLCM_2 <- tunningModels(data = potVal_GLCM, Site = 2, wl)
save(fit_potVal_GLCM_2, file="results_canopy/fit_potVal_GLCM_2.RData")
# Bootstrap validation 
fit_boot_potVal_GLCM_2 <- ApplyBootsClassification(data = potVal_GLCM, Site = 2, rasterPlots = plots2, 
                                                 en = fit_potVal_GLCM_2, outDir = outputDir, modelTag = "potVal_GLCM" )
save(fit_boot_potVal_GLCM_2, file="results_canopy/fit_boot_potVal_GLCM_2.RData")

#### Site 3
fit_potVal_GLCM_3 <- tunningModels(data = potVal_GLCM, Site = 3, wl)
save(fit_potVal_GLCM_3, file="results_canopy/fit_potVal_GLCM_3.RData")
# Bootstrap validation 
fit_boot_potVal_GLCM_3 <- ApplyBootsClassification(data = potVal_GLCM, Site = 3, rasterPlots = plots3, 
                                                 en = fit_potVal_GLCM_3, outDir = outputDir, modelTag = "potVal_GLCM" )
save(fit_boot_potVal_GLCM_3, file="results_canopy/fit_boot_potVal_GLCM_3.RData")

#### Site 4
fit_potVal_GLCM_4 <- tunningModels(data = potVal_GLCM, Site = 4, wl)
save(fit_potVal_GLCM_4, file="results_canopy/fit_potVal_GLCM_4.RData")
# Bootstrap validation 
fit_boot_potVal_GLCM_4 <- ApplyBootsClassification(data = potVal_GLCM, Site = 4, rasterPlots = plots4, 
                                                 en = fit_potVal_GLCM_4, outDir = outputDir, modelTag = "potVal_GLCM" )
save(fit_boot_potVal_GLCM_4, file="results_canopy/fit_boot_potVal_GLCM_4.RData")

##################
### rip it off ###
##################

#### Site 1
fit_rf_GLCM_1 <- tunningModels(data = rf_GLCM, Site = 1, wl)
save(fit_rf_GLCM_1, file="results_canopy/fit_rf_GLCM_1.RData")
# Bootstrap validation 
fit_boot_rf_GLCM_1 <- ApplyBootsClassification(data = rf_GLCM, Site = 1, rasterPlots = plots1,
                                             en = fit_rf_GLCM_1, outDir = outputDir, modelTag = "rf_GLCM" )
save(fit_boot_rf_GLCM_1, file="results_canopy/fit_boot_rf_GLCM_1.RData")

#### Site 2
fit_rf_GLCM_2 <- tunningModels(data = rf_GLCM, Site = 2, wl)
save(fit_rf_GLCM_2, file="results_canopy/fit_rf_GLCM_2.RData")
# Bootstrap validation 
fit_boot_rf_GLCM_2 <- ApplyBootsClassification(data = rf_GLCM, Site = 2, rasterPlots = plots2, 
                                             en = fit_rf_GLCM_2, outDir = outputDir, modelTag = "rf_GLCM" )
save(fit_boot_rf_GLCM_2, file="results_canopy/fit_boot_rf_GLCM_2.RData")

#### Site 3
fit_rf_GLCM_3 <- tunningModels(data = rf_GLCM, Site = 3, wl)
save(fit_rf_GLCM_3, file="results_canopy/fit_rf_GLCM_3.RData")
# Bootstrap validation 
fit_boot_rf_GLCM_3 <- ApplyBootsClassification(data = rf_GLCM, Site = 3, rasterPlots = plots3, 
                                             en = fit_rf_GLCM_3, outDir = outputDir, modelTag = "rf_GLCM" )
save(fit_boot_rf_GLCM_3, file="results_canopy/fit_boot_rf_GLCM_3.RData")

#### Site 4
fit_rf_GLCM_4 <- tunningModels(data = rf_GLCM, Site = 4, wl)
save(fit_rf_GLCM_4, file="results_canopy/fit_rf_GLCM_4.RData")
# Bootstrap validation 
fit_boot_rf_GLCM_4 <- ApplyBootsClassification(data = rf_GLCM, Site = 4, rasterPlots = plots4, 
                                             en = fit_rf_GLCM_4, outDir = outputDir, modelTag = "rf_GLCM" )
save(fit_boot_rf_GLCM_4, file="results_canopy/fit_boot_rf_GLCM_4.RData")

#------------------------#
# GLCM BN
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_BN_GLCM_1 <- tunningModels(data = potVal_BN_GLCM, Site = 1, wl)
save(fit_potVal_BN_GLCM_1, file="results_canopy/fit_potVal_BN_GLCM_1.RData")
# Bootstrap validation 
fit_boot_potVal_BN_GLCM_1 <- ApplyBootsClassification(data = potVal_BN_GLCM, Site = 1, rasterPlots = plots1, 
                                                 en = fit_potVal_BN_GLCM_1, outDir = outputDir, modelTag = "potVal_BN_GLCM" )
save(fit_boot_potVal_BN_GLCM_1, file="results_canopy/fit_boot_potVal_BN_GLCM_1.RData")

#### Site 2
fit_potVal_BN_GLCM_2 <- tunningModels(data = potVal_BN_GLCM, Site = 2, wl)
save(fit_potVal_BN_GLCM_2, file="results_canopy/fit_potVal_BN_GLCM_2.RData")
# Bootstrap validation 
fit_boot_potVal_BN_GLCM_2 <- ApplyBootsClassification(data = potVal_BN_GLCM, Site = 2, rasterPlots = plots2, 
                                                 en = fit_potVal_BN_GLCM_2, outDir = outputDir, modelTag = "potVal_BN_GLCM" )
save(fit_boot_potVal_BN_GLCM_2, file="results_canopy/fit_boot_potVal_BN_GLCM_2.RData")

#### Site 3
fit_potVal_BN_GLCM_3 <- tunningModels(data = potVal_BN_GLCM, Site = 3, wl)
save(fit_potVal_BN_GLCM_3, file="results_canopy/fit_potVal_BN_GLCM_3.RData")
# Bootstrap validation 
fit_boot_potVal_BN_GLCM_3 <- ApplyBootsClassification(data = potVal_BN_GLCM, Site = 3, rasterPlots = plots3, 
                                                 en = fit_potVal_BN_GLCM_3, outDir = outputDir, modelTag = "potVal_BN_GLCM" )
save(fit_boot_potVal_BN_GLCM_3, file="results_canopy/fit_boot_potVal_BN_GLCM_3.RData")

#### Site 4
fit_potVal_BN_GLCM_4 <- tunningModels(data = potVal_BN_GLCM, Site = 4, wl)
save(fit_potVal_BN_GLCM_4, file="results_canopy/fit_potVal_BN_GLCM_4.RData")
# Bootstrap validation 
fit_boot_potVal_BN_GLCM_4 <- ApplyBootsClassification(data = potVal_BN_GLCM, Site = 4, rasterPlots = plots4, 
                                                 en = fit_potVal_BN_GLCM_4, outDir = outputDir, modelTag = "potVal_BN_GLCM" )
save(fit_boot_potVal_BN_GLCM_4, file="results_canopy/fit_boot_potVal_BN_GLCM_4.RData")

##################
### rip it off ###
##################

#### Site 1
fit_rf_BN_GLCM_1 <- tunningModels(data = rf_BN_GLCM, Site = 1, wl)
save(fit_rf_BN_GLCM_1, file="results_canopy/fit_rf_BN_GLCM_1.RData")
# Bootstrap validation 
fit_boot_rf_BN_GLCM_1 <- ApplyBootsClassification(data = rf_BN_GLCM, Site = 1, rasterPlots = plots1, 
                                             en = fit_rf_BN_GLCM_1, outDir = outputDir, modelTag = "rf_BN_GLCM" )
save(fit_boot_rf_BN_GLCM_1, file="results_canopy/fit_boot_rf_BN_GLCM_1.RData")

#### Site 2
fit_rf_BN_GLCM_2 <- tunningModels(data = rf_BN_GLCM, Site = 2, wl)
save(fit_rf_BN_GLCM_2, file="results_canopy/fit_rf_BN_GLCM_2.RData")
# Bootstrap validation 
fit_boot_rf_BN_GLCM_2 <- ApplyBootsClassification(data = rf_BN_GLCM, Site = 2, rasterPlots = plots2, 
                                             en = fit_rf_BN_GLCM_2, outDir = outputDir, modelTag = "rf_BN_GLCM" )
save(fit_boot_rf_BN_GLCM_2, file="results_canopy/fit_boot_rf_BN_GLCM_2.RData")

#### Site 3
fit_rf_BN_GLCM_3 <- tunningModels(data = rf_BN_GLCM, Site = 3, wl)
save(fit_rf_BN_GLCM_3, file="results_canopy/fit_rf_BN_GLCM_3.RData")
# Bootstrap validation 
fit_boot_rf_BN_GLCM_3 <- ApplyBootsClassification(data = rf_BN_GLCM, Site = 3, rasterPlots = plots3, 
                                             en = fit_rf_BN_GLCM_3, outDir = outputDir, modelTag = "rf_BN_GLCM" )
save(fit_boot_rf_BN_GLCM_3, file="results_canopy/fit_boot_rf_BN_GLCM_3.RData")

#### Site 4
fit_rf_BN_GLCM_4 <- tunningModels(data = rf_BN_GLCM, Site = 4, wl)
save(fit_rf_BN_GLCM_4, file="results_canopy/fit_rf_BN_GLCM_4.RData")
# Bootstrap validation 
fit_boot_rf_BN_GLCM_4 <- ApplyBootsClassification(data = rf_BN_GLCM, Site = 4, rasterPlots = plots4, 
                                             en = fit_rf_BN_GLCM_4, outDir = outputDir, modelTag = "rf_BN_GLCM" )
save(fit_boot_rf_BN_GLCM_4, file="results_canopy/fit_boot_rf_BN_GLCM_4.RData")













