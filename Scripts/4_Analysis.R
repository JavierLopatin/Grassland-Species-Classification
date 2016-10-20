## R-Script - Analysis
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
rasterDir = "D:/Sp_Images"

setwd(home)

load("outputGOF.RData")
load("potVal_cover.RData")
load("rf_cover.RData")

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

## load the data
fit_potVal <-  read.table("Data/Fits_potVal.csv", sep = ",", header = T)
fit_rf     <-  read.table("Data/Fits_rf.csv", sep = ",", header = T)

## Obtain cover fits
setwd(rasterDir)

potVal_cover <- coverSummary("potVal")
rf_cover     <- coverSummary("rf")

potVal_cover$Species <- factor(potVal_cover$Species)
rf_cover$Species     <- factor(rf_cover$Species)

setwd(home)

save(potVal_cover, file="potVal_cover.RData")
save(rf_cover,     file="rf_cover.RData")

# Obtain R2, RMSE and Bias per species
gof_pv <- GOF(potVal_cover)
gof_rf <- GOF(rf_cover)

# erase empty "Species"
# gof_pv <- gof_pv[- seq(1,24,1),]
# gof_rf <- gof_rf[- seq(1,24,1),]

save(gof_pv, file="gof_pv.RData")
save(gof_rf, file="gof_rf.RData")

# get GOF per model/tag
#x = grep("potVal", outputGOF$Validation)
#gof_pv <- outputGOF[x, ]

#x = grep("rf", outputGOF$Validation)
#gof_rf <- outputGOF[x, ]

############################
### Best model rf MNF_BN ###
############################

x = grep("RF", rf_cover$Model)
y <- rf_cover[x, ] 

x = grep("MNF_BN", y$Normalization)
bestModel <- y[x, ]

x = grep("rf", outputGOF$Validation)
y <- outputGOF[x, ] 

x = grep("RF", y$Models)
y <- y[x, ] 

x = grep("MNF_BN", y$Normalization)
gof_best <- y[x, ]

# count for well, miss and over classifications
bestModel <- ClassPresence(bestModel)

### Analysis of architectural complexities
x = grep(paste0(c(18,19), collapse="|"), bestModel$Plot)
complex1 <- bestModel[x, ]

x = grep(paste0(c(13,14), collapse="|"), bestModel$Plot)
complex2 <- bestModel[x, ]

x = grep(paste0(c(9,10,11,12), collapse="|"), bestModel$Plot)
complex3 <- bestModel[x, ]

x = grep(paste0(c(15,16,17), collapse="|"), bestModel$Plot)
complex4 <- bestModel[x, ]

plot(complex1$Observed, complex1$Predicted, xlim=c(0,100), ylim=c(0,100))
plot(complex2$Observed, complex2$Predicted, xlim=c(0,100), ylim=c(0,100))
plot(complex3$Observed, complex3$Predicted, xlim=c(0,100), ylim=c(0,100))
plot(complex4$Observed, complex4$Predicted, xlim=c(0,100), ylim=c(0,100))

gof_1 <- GOFbest(complex1)
gof_1$complex <- "complex1"
gof_2 <- GOFbest(complex2)
gof_2$complex <- "complex2"
gof_3 <- GOFbest(complex3)
gof_3$complex <- "complex3"
gof_4 <- GOFbest(complex4)
gof_4$complex <- "complex4"

gof_complex <- rbind(gof_1, gof_2, gof_3, gof_4)

### analysis per cover percentage
cov1 <- bestModel[bestModel$Observed < 20, ]
cov2 <- bestModel[bestModel$Observed >= 20 & bestModel$Observed < 40, ]
cov3 <- bestModel[bestModel$Observed >= 40 & bestModel$Observed < 60, ]
cov4 <- bestModel[bestModel$Observed >= 60 & bestModel$Observed < 80, ]
cov5 <- bestModel[bestModel$Observed >= 80 & bestModel$Observed < 100, ]

cov1$CoverRange <- "cov1"
cov2$CoverRange <- "cov2"
cov3$CoverRange <- "cov3"
cov4$CoverRange <- "cov4"
cov5$CoverRange <- "cov5"

bestModel <- rbind(cov1, cov2, cov3, cov4, cov5)

### Analysis per PFT


### Analysis of resolutions
setwd(rasterDir)

## load the data
rf_MNF_BN <- read.table("Data/rf_MNF_BN.csv", sep = ",", header = T)

# load species cover dataset
species <- read.table("Data/Plots_Species.csv", header = T, sep=",")

## load plot images 
# resample done with gdal_translate, as: 
# Windows --> for %i in (*.tif) do gdal_translate -tr 2 2 %i 2/%i 
# Linux   --> FILES=inPath/*tif
#             for f in $FILES
#             do
#               gdal_translate gdal_translate -tr 2 2 ${f} 2/${f}
#             done

raster2 <- rasterList(fileExtantion = ".tif", folder = "Plots/resolution_analysis/2", dir=rasterDir)
raster4 <- rasterList(fileExtantion = ".tif", folder = "Plots/resolution_analysis/4", dir=rasterDir)
raster6 <- rasterList(fileExtantion = ".tif", folder = "Plots/resolution_analysis/6", dir=rasterDir)
raster8 <- rasterList(fileExtantion = ".tif", folder = "Plots/resolution_analysis/8", dir=rasterDir)
raster10 <- rasterList(fileExtantion = ".tif", folder = "Plots/resolution_analysis/10", dir=rasterDir)
raster12 <- rasterList(fileExtantion = ".tif", folder = "Plots/resolution_analysis/12", dir=rasterDir)

res2 <- tunningModels(classes = rf_MNF_BN$Species, 
                      spectra = rf_MNF_BN[, 3:length( rf_MNF_BN )], 
                      wl=seq(1,10,1))

BootsClassificationBest(classes = rf_MNF_BN$Species, 
                        spectra = rf_MNF_BN[, 3:length( rf_MNF_BN )],
                        en = res2, 
                        raster = raster2, 
                        boots = 10, 
                        outDir = file.path(rasterDir, "2"), 
                        modelTag = paste0("2_", modelTag),
                        plotName = plot_name)





