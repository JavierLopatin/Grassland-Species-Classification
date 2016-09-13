## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 


home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
# home = "~/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"

setwd(home)

load("fitAISALeaf.RData")
load("fitASDLeaf.RData")
load("boot_testLeaf.Rdata")

library(hyperSpec)

# load the data
data <- read.table("LeafHerbaceous.txt", sep = "", header = T)

# add Species
SpNames <- read.table("SpNamesLeafClip.csv", sep = "", header = T)
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

####################
### Plot results ###
####################

# feature selection
pdf(file = "Figures/plot_hyper.pdf", width=10, height=10)
par(mfrow = c(2,1), mai = c(0,0,0,0))
plot.classificationEnsemble(hyperASD$spc, fitASD)
mtext("A", side=3, line=0.5, adj=0, cex=1.3)

plot.classificationEnsemble(hyperAISA$spc, fitAISA)
mtext("B", side=3, line=0.5, adj=0, cex=1.3)
dev.off()

# Confusion matrix
library(beanplot)

# OA and Kappa
beanplot( boot_test$fit$OA.ASD, boot_test$fit$OA.AISA, boot_test$fit$kappa.ASD, boot_test$fit$kappa.AISA, 
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", 
          ll = 0, side = "b", log="", main = "Squared Pearson's correlation coefficient", 
          names=c("Total", "Tree", "Shrub", "Herb"), ylab = expression(r^2), ylim = c(0,1), 
          yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1)


