## R-Script - Classification
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 


home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"
# home = "~/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"

setwd(home)

pkgs<-c("caret", "raster", "rgdal", "e1071", "gtools", "doParallel", "autopls", "maptools", "hyperSpec")
lapply(pkgs, require, character.only=T)

# load the data
data <- read.table("LeafHerbaceous.txt", sep = "", header = T)
spectra <- data[, 2:length(data)]

# add Species
SpNames <- read.table("SpNamesLeafClip.csv", sep = "", header = T)
# add to data
data$Species <- SpNames$Species

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
plot(sample(hyperASD, 10))
plot(hyperASD, "spcprctl5")

# set a data with the AISA+ spectral characteristics
hyperAISA <- hyperASD[,, c(390~990)]
hyperAISA <- spc.bin (hyperAISA, 10)
nwl(hyperAISA)
plot(hyperAISA, "spcprctl5")

par(mfrow=c(1,2),lend = 1, mai = c(1.2, 1.2, 0.5, 0.5))
plot(hyperASD, "spcprctl5")
plot(hyperAISA, "spcprctl5")
