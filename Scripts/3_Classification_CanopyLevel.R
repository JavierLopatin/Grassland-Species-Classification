## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

#### run Clasification!

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
# home = "~/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"

setwd(home)

load("fit_potVal_1_spectra.RData")
load("fit_rf_1_spectra.RData")
load("fit_potVal_1_spectraBN.RData")
load("fit_rf_1_spectraBN.RData")
load("fit_potVal_1.RData")

### load the data
# Site1
potVal_1 <- read.table("data/potVal_trainningAreas_site1.csv", sep = ",", header = T)
rf_1     <- read.table("data/rf_trainningAreas_site1.csv", sep = ",", header = T)

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

##########################
## Site 1
##########################

#------------------------#
# Spectra
#------------------------#

# potVal
fit_potVal_1_spectra <- classificationEnsemble(potVal_1$Species, potVal_1[,2:62], wl)
plot.classificationEnsemble( potVal_1[,2:62]/10000, fit_potVal_1_spectra)
save(fit_potVal_1_spectra, file="fit_potVal_1_spectra1.Rdata")

# rip it off
fit_rf_1_spectra <-  classificationEnsemble(rf_1$Species, rf_1[,2:62], wl)
plot.classificationEnsemble( rf_1[,2:62]/10000, fit_rf_1_spectra)
save(fit_rf_1_spectra, file="fit_rf_1_spectra.Rdata")

#------------------------#
# Spectra BN
#------------------------#

# Normalization
BN <- function(x){ x / sqrt (rowSums (x ^ 2))}

# potVal
spect_bn <- BN(potVal_1[,2:62])
fit_potVal_1_spectraBN <- classificationEnsemble(potVal_1$Species, spect_bn, wl)
plot.classificationEnsemble( spect_bn, fit_potVal_1_spectraBN)
save(fit_potVal_1_spectraBN, file="fit_potVal_1_spectraBN.Rdata")

# rip it off
spect_bn <- BN(rf_1[,2:62])
fit_rf_1_spectraBN <-  classificationEnsemble(rf_1$Species, spect_bn, wl)
plot.classificationEnsemble( spect_bn, fit_rf_1_spectraBN)
save(fit_rf_1_spectraBN, file="fit_rf_1_spectraBN.Rdata")

#------------------------#
# MNF
#------------------------#

# potVal

# rip it off

#------------------------#
# MNF BN
#------------------------#

# potVal

# rip it off


#------------------------#
# GLCM
#------------------------#

# potVal

# rip it off


#------------------------#
# GLCM BN
#------------------------#

