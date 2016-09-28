## R-Script - 3_Classification_LeafLevel.R
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
setwd(home)

# load required libraries
library(hyperSpec)

# load the data
data <- read.table("/data/LeafHerbaceous.txt", sep = "", header = T)

# load species names
SpNames <- read.table("/data/SpNamesLeafClip.csv", sep = "", header = T)
# add to data
data$Species <- SpNames$Species
# erase bad data
data <- subset(data, Species != "Erase")
data$Species <- factor(data$Species)
# spectral bands
spectra <- data[, 2:(length(data)-1)]

# create an hyperSpec object
new("hyperSpec")
hyperASD <- new("hyperSpec", spc=spectra, data=data, wavelength = seq(350, 2500, 1),
                label=list(spc="Reflectance", .wavelength =  expression(lambda(nm))))

plot(hyperASD)
plot(sample(hyperASD, 3))
plot(hyperASD[2,, 390~1000])
plot(hyperASD, "spcprctl5")

# eliminate the noisy 350-390 bands
hyperASD <- hyperASD[,, c(390~max)]
hyperASD <- spc.bin (hyperASD, 10)
plot(sample(hyperASD, 10))
plot(hyperASD, "spcprctl5")
nwl(hyperASD)

# set a data with the AISA+ spectral characteristics
hyperAISA <- hyperASD[,, c(390~990)]
nwl(hyperAISA)
plot(hyperAISA)

par(mfrow=c(1,2),lend = 1, mai = c(1.2, 1.2, 0.5, 0.5))
plot(hyperASD, "spcprctl5")
plot(hyperAISA, "spcprctl5")

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

# With the AISA+ band seting
fitAISA <- classificationEnsemble(hyperAISA@data$Species, hyperAISA$spc, hyperAISA@wavelength)
plot.classificationEnsemble(hyperAISA$spc, fitAISA)
save(fitAISA, file="fitAISALeaf.Rdata")

# ASD full range
fitASD <-  classificationEnsemble(hyperASD@data$Species, hyperASD$spc, hyperASD@wavelength)
plot.classificationEnsemble(hyperASD$spc, fitASD)
save(fitASD, file="fitASDLeaf.Rdata")

##################################
## Bootstrap significance test ###
##################################

boot_test <- significanceTest_LeafLevel(data, fitASD, fitAISA)
save(boot_test, file="boot_testLeaf.Rdata")

# Hist do not have to overlap with 0, otherwise is not significant
# The black "zero-line" needs to be left of the blue "alpha-line". The green line is just the upper quantile.
par(mfrow=c(1,2), mar=c(2,3,3,1))
main <- c("OA Leaf level", "Kappa Leaf level")
for(i in 1:2){
  hist(unlist(boot_test$boot_test[i]), main=main[i], col="grey", border="white", xlab="", ylab="")
  abline(v=quantile(unlist(boot_test$boot_test[i]), probs=c(0.05, 0.95)), col=c("blue", "green"))
  abline(v=0, col=c("black"))
  box()
}
