################################################################################
## R-Script - 3_Classification_CanopyLevel.R                                  ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## Manuscript: Hyperspectral classification of grassland species: towards an  ##
## application for semi-automatic field surveys                               ##
##                                                                            ##
## description: Classification procidure for canopy-level spectras            ## 
##                                                                            ##
################################################################################

library(raster)

home = "D:/Sp_Images"

setwd(home)

#load("Class_Canopy.RData")

## load the data
potVal_spec    <- read.table("Data/potVal_spec.csv", sep = ",", header = T)
rf_spec        <- read.table("Data/rf_spec.csv", sep = ",", header = T)
potVal_spec_BN <- read.table("Data/potVal_spec_BN.csv", sep = ",", header = T)
rf_spec_BN     <- read.table("Data/rf_spec_BN.csv", sep = ",", header = T)
potVal_MNF     <- read.table("Data/potVal_MNF.csv", sep = ",", header = T)
rf_MNF         <- read.table("Data/rf_MNF.csv", sep = ",", header = T)
potVal_MNF_BN  <- read.table("Data/potVal_MNF_BN.csv", sep = ",", header = T)
rf_MNF_BN      <- read.table("Data/rf_MNF_BN.csv", sep = ",", header = T)
potVal_GLCM    <- read.table("Data/potVal_GLCM.csv", sep = ",", header = T)
rf_GLCM        <- read.table("Data/rf_GLCM.csv", sep = ",", header = T)
potVal_GLCM_BN <- read.table("Data/potVal_GLCM_BN.csv", sep = ",", header = T)
rf_GLCM_BN     <- read.table("Data/rf_GLCM_BN.csv", sep = ",", header = T)


# wavelength
wl <- c( 398, 407, 415, 424, 432, 441, 450, 459, 468, 477, 486, 495, 504, 513, 522, 531, 540, 550, 558, 568,
         577, 587, 596, 605, 615, 624, 633, 643, 652, 661, 671, 680, 690, 699, 708, 717, 727, 736, 746, 755,
         765, 775, 784, 794, 803, 813, 822, 832, 842, 851, 861, 870, 880, 890, 899, 908, 918, 928, 937, 947, 957)
  
# load species cover dataset
species <- read.table("Data/Plots_Species.csv", header = T, sep=",")


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


## load plot images 
raster_spec <-  rasterList(fileExtantion = ".tif", folder = "Plots/Plots_spec", dir=home)
r <- stack( paste0(home, "/Plots/Plots_spec/plot_17.dat") )
raster_spec[[11]] <- r; names(raster_spec) <- c("plot_10", "plot_11", "plot_12", "plot_13", "plot_14", "plot_15", "plot_16", "plot_18", "plot_19", "plot_9", "plot_17")

raster_spec_BN <-  rasterList(fileExtantion = ".tif", folder = "Plots/Plots_spec_BN", dir=home)
raster_MNF     <-  rasterList(fileExtantion = ".tif", folder = "Plots/Plots_MNF", dir=home)
raster_MNF_BN  <-  rasterList(fileExtantion = ".tif", folder = "Plots/Plots_MNF_BN", dir=home)
raster_GLCM    <-  rasterList(fileExtantion = ".tif", folder = "Plots/Plots_GLCM", dir=home)
raster_GLCM_BN <-  rasterList(fileExtantion = ".tif", folder = "Plots/Plots_GLCM_BN", dir=home)


##########################
### Run Classification ###
##########################

setwd(home)

save.image("Class_Canopy.RData")

subplotDir = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp/subplots"

#------------------------#
# Spectra                #
#------------------------#

ApplyModels(valData = species, 
            potVal = potVal_spec, 
            rf = rf_spec, 
            raster_List = raster_spec, 
            wl = wl, 
            modelTag = "spect",
            boots = 10)

#------------------------#
# Spectra BN             #
#------------------------#

ApplyModels(valData = species, 
            potVal = potVal_spec_BN, 
            rf = rf_spec_BN,
            raster_List = raster_spec_BN, 
            wl = wl, 
            modelTag = "spect_BN",
            boots = 10)

#------------------------#
# MNF                    #
#------------------------#

ApplyModels(valData = species, 
            potVal = potVal_MNF, 
            rf = rf_MNF, 
            raster_List = raster_MNF, 
            wl = seq(1,10,1), 
            modelTag = "MNF",
            boots = 10)

#------------------------#
# MNF BN                 #
#------------------------#

ApplyModels(valData = species, 
            potVal = potVal_MNF_BN, 
            rf = rf_MNF_BN, 
            raster_List = raster_MNF_BN, 
            wl =  seq(1,10,1), 
            modelTag = "MNF_BN",
            boots = 10)
