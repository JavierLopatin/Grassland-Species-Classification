## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

#### run Clasification!

library(raster)

home = "F:/Sp_Images/temp"

setwd(home)

## load the data
potVal_spec    <- read.table("Data/potVal_spec.csv", sep = ",", header = T)
rf_spec        <- read.table("Data/rf_spec.csv", sep = ",", header = T)
potVal_spec_BN <- read.table("Data/potVal_spec_BN.csv", sep = ",", header = T)
rf_spec_BN     <- read.table("Data/rf_spec_BN.csv", sep = ",", header = T)
potVal_MNF     <- read.table("Data/potVal_MNF.csv", sep = ",", header = T)
rf_MNF         <- read.table("Data/rf_MNF.csv", sep = ",", header = T)
potVal_MNF_BN  <- read.table("Data/potVal_MNF_BN.csv", sep = ",", header = T)
rf_MNF_BN      <- read.table("Data/rf_MNF_BN.csv", sep = ",", header = T)
potVal_GLCM    <- read.table("Data/potVal_GLCM_all.csv", sep = ",", header = T)
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


#------------------------#
# Spectra                #
#------------------------#
valData = species
potVal = potVal_spec
rf = rf_spec
rasterList = raster_spec

ApplyModels <- function(valData, potVal, rf, rasterList, wl){
  
  ## Apply functions tunningModels and BootsClassification per plot
  
  for (i in 1:length(rasterList)){ # four sites
    
    # obtain the validation data per plot
    raster = rasterList[[i]]
    plot = unique(na.omit(as.numeric(unlist(strsplit( names( rasterList )[[i]], "[^0-9]+")))))
    plot_name = paste0("plot_", plot)  
    
    x = grep( plot, valData$Plot )  
    data = valData[x, ] 
    data$Species <- factor( data$Species ) # reset species Levels
    # get species to classify in the plot
    classes = unique( data$Species )
    
    # obtain subset of data of potVal and rf that include these species
    x = grep( paste(classes, collapse = "|") , potVal$Species )
    data_potVal = potVal[x, ]
    data_potVal$Species <- factor(data_potVal$Species) 
    
    x = grep( paste(classes, collapse = "|") , rf$Species )
    data_rf = rf[x, ]
    data_rf$Species <- factor(data_rf$Species) 
    
    ### apply tunningModels function
    fit_potVal <-  tunningModels(classes = data_potVal$Species, 
                                 spectra = data_potVal[, 3:length( data_potVal )],
                                 wl = wl)
    
    fit_rf     <-  tunningModels(classes = data_rfl$Species, 
                                 spectra = data_rf[, 3:length( data_rf )],
                                 wl = wl)
    
    
    }
  
  
  
}





#------------------------#
# Spectra BN
#------------------------#

##############
### potVal ###
##############

#### Site 1
#fit_potVal_BN_1 <- tunningModels(data = potVal_BN, Site = 1, wl)
#save(fit_potVal_BN_1, file="results_canopy/fit_potVal_BN_1.RData")
# Bootstrap validation 
fit_boot_potVal_BN_1 <- ApplyBootsClassification(data = potVal_BN, Site = 1, rasterPlots = plots1_BN, boots=10, 
                                              en = fit_potVal_BN_1, outDir = outputDir, modelTag = "potVal_BN" )

#### Site 2
#fit_potVal_BN_2 <- tunningModels(data = potVal_BN, Site = 2, wl)
#save(fit_potVal_BN_2, file="results_canopy/fit_potVal_BN_2.RData")
# Bootstrap validation 
fit_boot_potVal_BN_2 <- ApplyBootsClassification(data = potVal_BN, Site = 2, rasterPlots = plots2_BN, boots=10,
                                              en = fit_potVal_BN_2, outDir = outputDir, modelTag = "potVal_BN" )

#### Site 3
#fit_potVal_BN_3 <- tunningModels(data = potVal_BN, Site = 3, wl)
#save(fit_potVal_BN_3, file="results_canopy/fit_potVal_BN_3.RData")
# Bootstrap validation 
fit_boot_potVal_BN_3 <- ApplyBootsClassification(data = potVal_BN, Site = 3, rasterPlots = plots3_BN, boots=10,
                                              en = fit_potVal_BN_3, outDir = outputDir, modelTag = "potVal_BN" )

#### Site 4
#fit_potVal_BN_4 <- tunningModels(data = potVal_BN, Site = 4, wl)
#save(fit_potVal_BN_4, file="results_canopy/fit_potVal_BN_4.RData")
# Bootstrap validation 
fit_boot_potVal_BN_4 <- ApplyBootsClassification(data = potVal_BN, Site = 4, rasterPlots = plots4_BN, boots=10,
                                              en = fit_potVal_BN_4, outDir = outputDir, modelTag = "potVal_BN" )

##################
### rip it off ###
##################

#### Site 1
#fit_rf_BN_1 <- tunningModels(data = rf_BN, Site = 1, wl)
#save(fit_rf_BN_1, file="results_canopy/fit_rf_BN_1.RData")
# Bootstrap validation 
fit_boot_rf_BN_1 <- ApplyBootsClassification(data = rf_BN, Site = 1, rasterPlots = plots1_BN, boots=10, 
                                          en = fit_rf_BN_1, outDir = outputDir, modelTag = "rf_BN" )

#### Site 2
#fit_rf_BN_2 <- tunningModels(data = rf_BN, Site = 2, wl)
#save(fit_rf_BN_2, file="results_canopy/fit_rf_BN_2.RData")
# Bootstrap validation 
fit_boot_rf_BN_2 <- ApplyBootsClassification(data = rf_BN, Site = 2, rasterPlots = plots2_BN,  boots=10,
                                          en = fit_rf_BN_2, outDir = outputDir, modelTag = "rf_BN" )

#### Site 3
fit_rf_BN_3 <- tunningModels(data = rf_BN, Site = 3, wl)
save(fit_rf_BN_3, file="results_canopy/fit_rf_BN_3.RData")
# Bootstrap validation 
fit_boot_rf_BN_3 <- ApplyBootsClassification(data = rf_BN, Site = 3, rasterPlots = plots3_BN,  boots=10,
                                          en = fit_rf_BN_3, outDir = outputDir, modelTag = "rf_BN" )

#### Site 4
fit_rf_BN_4 <- tunningModels(data = rf_BN, Site = 4, wl)
save(fit_rf_BN_4, file="results_canopy/fit_rf_BN_4.RData")
# Bootstrap validation 
fit_boot_rf_BN_4 <- ApplyBootsClassification(data = rf_BN, Site = 4, rasterPlots = plots4_BN,  boots=10,
                                          en = fit_rf_BN_4, outDir = outputDir, modelTag = "rf_BN" )

#------------------------#
# MNF
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_MNF_1 <- tunningModels(data = potVal_MNF, Site = 1, wl= seq(1,10,1))
save(fit_potVal_MNF_1, file="results_canopy/fit_potVal_MNF_1.RData")
# Bootstrap validation 
fit_boot_potVal_MNF_1 <- ApplyBootsClassification(data = potVal_MNF, Site = 1, rasterPlots = plots1_MNF, boots=10, 
                                                 en = fit_potVal_MNF_1, outDir = outputDir, modelTag = "potVal_MNF" )

#### Site 2
fit_potVal_MNF_2 <- tunningModels(data = potVal_MNF, Site = 2, wl= seq(1,10,1))
save(fit_potVal_MNF_2, file="results_canopy/fit_potVal_MNF_2.RData")
# Bootstrap validation 
fit_boot_potVal_MNF_2 <- ApplyBootsClassification(data = potVal_MNF, Site = 2, rasterPlots = plots2_MNF,  boots=10,
                                                 en = fit_potVal_MNF_2, outDir = outputDir, modelTag = "potVal_MNF" )

#### Site 3
fit_potVal_MNF_3 <- tunningModels(data = potVal_MNF, Site = 3, wl= seq(1,10,1))
save(fit_potVal_MNF_3, file="results_canopy/fit_potVal_MNF_3.RData")
# Bootstrap validation 
fit_boot_potVal_MNF_3 <- ApplyBootsClassification(data = potVal_MNF, Site = 3, rasterPlots = plots3_MNF,  boots=10,
                                                 en = fit_potVal_MNF_3, outDir = outputDir, modelTag = "potVal_MNF" )

#### Site 4
fit_potVal_MNF_4 <- tunningModels(data = potVal_MNF, Site = 4, wl= seq(1,10,1))
save(fit_potVal_MNF_4, file="results_canopy/fit_potVal_MNF_4.RData")
# Bootstrap validation 
fit_boot_potVal_MNF_4 <- ApplyBootsClassification(data = potVal_MNF, Site = 4, rasterPlots = plots4_MNF, boots=10, 
                                                 en = fit_potVal_MNF_4, outDir = outputDir, modelTag = "potVal_MNF" )

##################
### rip it off ###
##################

#### Site 1
fit_rf_MNF_1 <- tunningModels(data = rf_MNF, Site = 1, wl = seq(1,10,1))
save(fit_rf_MNF_1, file="results_canopy/fit_rf_MNF_1.RData")
# Bootstrap validation 
fit_boot_rf_MNF_1 <- ApplyBootsClassification(data = rf_MNF, Site = 1, rasterPlots = plots1_MNF,  boots=10,
                                             en = fit_rf_MNF_1, outDir = outputDir, modelTag = "rf_MNF" )

#### Site 2
fit_rf_MNF_2 <- tunningModels(data = rf_MNF, Site = 2, wl= seq(1,10,1))
save(fit_rf_MNF_2, file="results_canopy/fit_rf_MNF_2.RData")
# Bootstrap validation 
fit_boot_rf_MNF_2 <- ApplyBootsClassification(data = rf_MNF, Site = 2, rasterPlots = plots2_MNF, boots=10, 
                                             en = fit_rf_MNF_2, outDir = outputDir, modelTag = "rf_MNF" )

#### Site 3
fit_rf_MNF_3 <- tunningModels(data = rf_MNF, Site = 3, wl= seq(1,10,1))
save(fit_rf_MNF_3, file="results_canopy/fit_rf_MNF_3.RData")
# Bootstrap validation 
fit_boot_rf_MNF_3 <- ApplyBootsClassification(data = rf_MNF, Site = 3, rasterPlots = plots3_MNF, boots=10,
                                             en = fit_rf_MNF_3, outDir = outputDir, modelTag = "rf_MNF" )

#### Site 4
fit_rf_MNF_4 <- tunningModels(data = rf_MNF, Site = 4, wl= seq(1,10,1))
save(fit_rf_MNF_4, file="results_canopy/fit_rf_MNF_4.RData")
# Bootstrap validation 
fit_boot_rf_MNF_4 <- ApplyBootsClassification(data = rf_MNF, Site = 4, rasterPlots = plots4_MNF, boots=10, 
                                             en = fit_rf_MNF_4, outDir = outputDir, modelTag = "rf_MNF" )


#------------------------#
# MNF BN
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_BN_MNF_1 <- tunningModels(data = potVal_BN_MNF, Site = 1, wl= seq(1,10,1))
save(fit_potVal_BN_MNF_1, file="results_canopy/fit_potVal_BN_MNF_1.RData")
# Bootstrap validation 
fit_boot_potVal_BN_MNF_1 <- ApplyBootsClassification(data = potVal_BN_MNF, Site = 1, rasterPlots = plots1_BN_MNF, boots=10, 
                                                 en = fit_potVal_BN_MNF_1, outDir = outputDir, modelTag = "potVal_BN_MNF" )

#### Site 2
fit_potVal_BN_MNF_2 <- tunningModels(data = potVal_BN_MNF, Site = 2, wl= seq(1,10,1))
save(fit_potVal_BN_MNF_2, file="results_canopy/fit_potVal_BN_MNF_2.RData")
# Bootstrap validation 
fit_boot_potVal_BN_MNF_2 <- ApplyBootsClassification(data = potVal_BN_MNF, Site = 2, rasterPlots = plots2_BN_MNF, boots=10, 
                                                 en = fit_potVal_BN_MNF_2, outDir = outputDir, modelTag = "potVal_BN_MNF" )

#### Site 3
fit_potVal_BN_MNF_3 <- tunningModels(data = potVal_BN_MNF, Site = 3, wl= seq(1,10,1))
save(fit_potVal_BN_MNF_3, file="results_canopy/fit_potVal_BN_MNF_3.RData")
# Bootstrap validation 
fit_boot_potVal_BN_MNF_3 <- ApplyBootsClassification(data = potVal_BN_MNF, Site = 3, rasterPlots = plots3_BN_MNF,  boots=10,
                                                 en = fit_potVal_BN_MNF_3, outDir = outputDir, modelTag = "potVal_BN_MNF" )

#### Site 4
fit_potVal_BN_MNF_4 <- tunningModels(data = potVal_BN_MNF, Site = 4, wl= seq(1,10,1))
save(fit_potVal_BN_MNF_4, file="results_canopy/fit_potVal_BN_MNF_4.RData")
# Bootstrap validation 
fit_boot_potVal_BN_MNF_4 <- ApplyBootsClassification(data = potVal_BN_MNF, Site = 4, rasterPlots = plots4_BN_MNF,  boots=10,
                                                 en = fit_potVal_BN_MNF_4, outDir = outputDir, modelTag = "potVal_BN_MNF" )

##################
### rip it off ###
##################

#### Site 1
fit_rf_BN_MNF_1 <- tunningModels(data = rf_BN_MNF, Site = 1, wl= seq(1,10,1))
save(fit_rf_BN_MNF_1, file="results_canopy/fit_rf_BN_MNF_1.RData")
# Bootstrap validation 
fit_boot_rf_BN_MNF_1 <- ApplyBootsClassification(data = rf_BN_MNF, Site = 1, rasterPlots = plots1_BN_MNF, boots=10, 
                                             en = fit_rf_BN_MNF_1, outDir = outputDir, modelTag = "rf_BN_MNF" )

#### Site 2
fit_rf_BN_MNF_2 <- tunningModels(data = rf_BN_MNF, Site = 2, wl= seq(1,10,1))
save(fit_rf_BN_MNF_2, file="results_canopy/fit_rf_BN_MNF_2.RData")
# Bootstrap validation 
fit_boot_rf_BN_MNF_2 <- ApplyBootsClassification(data = rf_BN_MNF, Site = 2, rasterPlots = plots2_BN_MNF, boots=10, 
                                             en = fit_rf_BN_MNF_2, outDir = outputDir, modelTag = "rf_BN_MNF" )

#### Site 3
fit_rf_BN_MNF_3 <- tunningModels(data = rf_BN_MNF, Site = 3, wl= seq(1,10,1))
save(fit_rf_BN_MNF_3, file="results_canopy/fit_rf_BN_MNF_3.RData")
# Bootstrap validation 
fit_boot_rf_BN_MNF_3 <- ApplyBootsClassification(data = rf_BN_MNF, Site = 3, rasterPlots = plots3_BN_MNF,  boots=10,
                                             en = fit_rf_BN_MNF_3, outDir = outputDir, modelTag = "rf_BN_MNF" )

#### Site 4
fit_rf_BN_MNF_4 <- tunningModels(data = rf_BN_MNF, Site = 4, wl= seq(1,10,1))
save(fit_rf_BN_MNF_4, file="results_canopy/fit_rf_BN_MNF_4.RData")
# Bootstrap validation 
fit_boot_rf_BN_MNF_4 <- ApplyBootsClassification(data = rf_BN_MNF, Site = 4, rasterPlots = plots4_BN_MNF,  boots=10,
                                             en = fit_rf_BN_MNF_4, outDir = outputDir, modelTag = "rf_BN_MNF" )

#------------------------#
# GLCM
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_GLCM_1 <- tunningModels(data = potVal_GLCM, Site = 1, wl= seq(1,6,1))
save(fit_potVal_GLCM_1, file="results_canopy/fit_potVal_GLCM_1.RData")
# Bootstrap validation 
fit_boot_potVal_GLCM_1 <- ApplyBootsClassification(data = potVal_GLCM, Site = 1, rasterPlots = plots1_GLCM,  boots=10,
                                                 en = fit_potVal_GLCM_1, outDir = outputDir, modelTag = "potVal_GLCM" )
save(fit_boot_potVal_GLCM_1, file="results_canopy/fit_boot_potVal_GLCM_1.RData")

#### Site 2
fit_potVal_GLCM_2 <- tunningModels(data = potVal_GLCM, Site = 2, wl= seq(1,6,1))
save(fit_potVal_GLCM_2, file="results_canopy/fit_potVal_GLCM_2.RData")
# Bootstrap validation 
fit_boot_potVal_GLCM_2 <- ApplyBootsClassification(data = potVal_GLCM, Site = 2, rasterPlots = plots2_GLCM, boots=10, 
                                                 en = fit_potVal_GLCM_2, outDir = outputDir, modelTag = "potVal_GLCM" )

#### Site 3
fit_potVal_GLCM_3 <- tunningModels(data = potVal_GLCM, Site = 3, wl= seq(1,6,1))
save(fit_potVal_GLCM_3, file="results_canopy/fit_potVal_GLCM_3.RData")
# Bootstrap validation 
fit_boot_potVal_GLCM_3 <- ApplyBootsClassification(data = potVal_GLCM, Site = 3, rasterPlots = plots3_GLCM, boots=10,
                                                 en = fit_potVal_GLCM_3, outDir = outputDir, modelTag = "potVal_GLCM" )

#### Site 4
fit_potVal_GLCM_4 <- tunningModels(data = potVal_GLCM, Site = 4, wl= seq(1,6,1))
save(fit_potVal_GLCM_4, file="results_canopy/fit_potVal_GLCM_4.RData")
# Bootstrap validation 
fit_boot_potVal_GLCM_4 <- ApplyBootsClassification(data = potVal_GLCM, Site = 4, rasterPlots = plots4_GLCM, boots=10, 
                                                 en = fit_potVal_GLCM_4, outDir = outputDir, modelTag = "potVal_GLCM" )

##################
### rip it off ###
##################

#### Site 1
fit_rf_GLCM_1 <- tunningModels(data = rf_GLCM, Site = 1, wl= seq(1,6,1))
save(fit_rf_GLCM_1, file="results_canopy/fit_rf_GLCM_1.RData")
# Bootstrap validation 
fit_boot_rf_GLCM_1 <- ApplyBootsClassification(data = rf_GLCM, Site = 1, rasterPlots = plots1_GLCM, boots=10,
                                             en = fit_rf_GLCM_1, outDir = outputDir, modelTag = "rf_GLCM" )

#### Site 2
fit_rf_GLCM_2 <- tunningModels(data = rf_GLCM, Site = 2, wl= seq(1,6,1))
save(fit_rf_GLCM_2, file="results_canopy/fit_rf_GLCM_2.RData")
# Bootstrap validation 
fit_boot_rf_GLCM_2 <- ApplyBootsClassification(data = rf_GLCM, Site = 2, rasterPlots = plots2_GLCM,  boots=10,
                                             en = fit_rf_GLCM_2, outDir = outputDir, modelTag = "rf_GLCM" )

#### Site 3
fit_rf_GLCM_3 <- tunningModels(data = rf_GLCM, Site = 3, wl= seq(1,6,1))
save(fit_rf_GLCM_3, file="results_canopy/fit_rf_GLCM_3.RData")
# Bootstrap validation 
fit_boot_rf_GLCM_3 <- ApplyBootsClassification(data = rf_GLCM, Site = 3, rasterPlots = plots3_GLCM,  boots=10,
                                             en = fit_rf_GLCM_3, outDir = outputDir, modelTag = "rf_GLCM" )

#### Site 4
fit_rf_GLCM_4 <- tunningModels(data = rf_GLCM, Site = 4, wl= seq(1,6,1))
save(fit_rf_GLCM_4, file="results_canopy/fit_rf_GLCM_4.RData")
# Bootstrap validation 
fit_boot_rf_GLCM_4 <- ApplyBootsClassification(data = rf_GLCM, Site = 4, rasterPlots = plots4_GLCM,  boots=10,
                                             en = fit_rf_GLCM_4, outDir = outputDir, modelTag = "rf_GLCM" )

#------------------------#
# GLCM BN
#------------------------#

##############
### potVal ###
##############

#### Site 1
fit_potVal_BN_GLCM_1 <- tunningModels(data = potVal_BN_GLCM, Site = 1, wl= seq(1,6,1))
save(fit_potVal_BN_GLCM_1, file="results_canopy/fit_potVal_BN_GLCM_1.RData")
# Bootstrap validation 
fit_boot_potVal_BN_GLCM_1 <- ApplyBootsClassification(data = potVal_BN_GLCM, Site = 1, rasterPlots = plots1_BN_GLCM, boots=10, 
                                                 en = fit_potVal_BN_GLCM_1, outDir = outputDir, modelTag = "potVal_BN_GLCM" )

#### Site 2
fit_potVal_BN_GLCM_2 <- tunningModels(data = potVal_BN_GLCM, Site = 2, wl= seq(1,6,1))
save(fit_potVal_BN_GLCM_2, file="results_canopy/fit_potVal_BN_GLCM_2.RData")
# Bootstrap validation 
fit_boot_potVal_BN_GLCM_2 <- ApplyBootsClassification(data = potVal_BN_GLCM, Site = 2, rasterPlots = plots2_BN_GLCM,  boots=10,
                                                 en = fit_potVal_BN_GLCM_2, outDir = outputDir, modelTag = "potVal_BN_GLCM" )

#### Site 3
fit_potVal_BN_GLCM_3 <- tunningModels(data = potVal_BN_GLCM, Site = 3, wl= seq(1,6,1))
save(fit_potVal_BN_GLCM_3, file="results_canopy/fit_potVal_BN_GLCM_3.RData")
# Bootstrap validation 
fit_boot_potVal_BN_GLCM_3 <- ApplyBootsClassification(data = potVal_BN_GLCM, Site = 3, rasterPlots = plots3_BN_GLCM,  boots=10, 
                                                 en = fit_potVal_BN_GLCM_3, outDir = outputDir, modelTag = "potVal_BN_GLCM" )

#### Site 4
fit_potVal_BN_GLCM_4 <- tunningModels(data = potVal_BN_GLCM, Site = 4, wl= seq(1,6,1))
save(fit_potVal_BN_GLCM_4, file="results_canopy/fit_potVal_BN_GLCM_4.RData")
# Bootstrap validation 
fit_boot_potVal_BN_GLCM_4 <- ApplyBootsClassification(data = potVal_BN_GLCM, Site = 4, rasterPlots = plots4_BN_GLCM,  boots=10,
                                                 en = fit_potVal_BN_GLCM_4, outDir = outputDir, modelTag = "potVal_BN_GLCM" )

##################
### rip it off ###
##################

#### Site 1
fit_rf_BN_GLCM_1 <- tunningModels(data = rf_BN_GLCM, Site = 1, wl= seq(1,6,1))
save(fit_rf_BN_GLCM_1, file="results_canopy/fit_rf_BN_GLCM_1.RData")
# Bootstrap validation 
fit_boot_rf_BN_GLCM_1 <- ApplyBootsClassification(data = rf_BN_GLCM, Site = 1, rasterPlots = plots1_BN_GLCM,  boots=10,
                                             en = fit_rf_BN_GLCM_1, outDir = outputDir, modelTag = "rf_BN_GLCM" )

#### Site 2
fit_rf_BN_GLCM_2 <- tunningModels(data = rf_BN_GLCM, Site = 2, wl= seq(1,6,1))
save(fit_rf_BN_GLCM_2, file="results_canopy/fit_rf_BN_GLCM_2.RData")
# Bootstrap validation 
fit_boot_rf_BN_GLCM_2 <- ApplyBootsClassification(data = rf_BN_GLCM, Site = 2, rasterPlots = plots2_BN_GLCM,  boots=10,
                                             en = fit_rf_BN_GLCM_2, outDir = outputDir, modelTag = "rf_BN_GLCM" )

#### Site 3
fit_rf_BN_GLCM_3 <- tunningModels(data = rf_BN_GLCM, Site = 3, wl= seq(1,6,1))
save(fit_rf_BN_GLCM_3, file="results_canopy/fit_rf_BN_GLCM_3.RData")
# Bootstrap validation 
fit_boot_rf_BN_GLCM_3 <- ApplyBootsClassification(data = rf_BN_GLCM, Site = 3, rasterPlots = plots3_BN_GLCM,  boots=10,
                                             en = fit_rf_BN_GLCM_3, outDir = outputDir, modelTag = "rf_BN_GLCM" )


#### Site 4
fit_rf_BN_GLCM_4 <- tunningModels(data = rf_BN_GLCM, Site = 4, wl= seq(1,6,1))
save(fit_rf_BN_GLCM_4, file="results_canopy/fit_rf_BN_GLCM_4.RData")
# Bootstrap validation 
fit_boot_rf_BN_GLCM_4 <- ApplyBootsClassification(data = rf_BN_GLCM, Site = 4, rasterPlots = plots4_BN_GLCM, boots=10, 
                                             en = fit_rf_BN_GLCM_4, outDir = outputDir, modelTag = "rf_BN_GLCM" )






