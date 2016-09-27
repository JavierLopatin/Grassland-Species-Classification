## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

#### run Clasification!

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
# home = "~/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"

setwd(home)

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
plot.classificationEnsemble( potVal_GLCM[,3:length(potVal_GLCM)], fit_potVal_GLCM, xlab_tag=expression(lambda(nm)) )
save(fit_potVal_GLCM, file="results_canopy/fit_potVal_GLCM.RData")

# rip it off
fit_rf_GLCM <- classificationEnsemble(rf_GLCM$Species, rf_GLCM[,3:length(rf_GLCM)], seq(1,6,1))
plot.classificationEnsemble( rf_GLCM[,3:length(rf_GLCM)], fit_rf_GLCM, xlab_tag=expression(lambda(nm)) )
save(fit_rf_GLCM, file="results_canopy/fit_rf_GLCM.RData")

#------------------------#
# GLCM BN
#------------------------#
# potVal
fit_potVal_BN_GLCM <- classificationEnsemble(potVal_BN_GLCM$Species, potVal_BN_GLCM[,3:length(potVal_BN_GLCM)], seq(1,6,1))
plot.classificationEnsemble( potVal_BN_GLCM[,3:length(potVal_BN_GLCM)], fit_potVal_BN_GLCM, xlab_tag=expression(lambda(nm)) )
save(fit_potVal_BN_GLCM, file="results_canopy/fit_potVal_BN_GLCM.RData")

# rip it off
fit_rf_BN_GLCM <- classificationEnsemble(rf_BN_GLCM$Species, rf_BN_GLCM[,3:length(rf_BN_GLCM)], seq(1,6,1))
plot.classificationEnsemble( rf_BN_GLCM[,3:length(rf_BN_GLCM)], fit_rf_BN_GLCM, xlab_tag=expression(lambda(nm)) )
save(fit_rf_BN_GLCM, file="results_canopy/fit_rf_BN_GLCM.RData")
