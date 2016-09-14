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

# Hist do not have to overlap with 0, otherwise is not significant
# The black "zero-line" needs to be left of the blue "alpha-line". The green line is just the upper quantile.
par(mfrow=c(1,2), mar=c(2,3,3,1))
main <- c("OA Leaf level", "Kappa LEaf level")
for(i in 1:2){
  hist(unlist(boot_test$boot_test[i]), main=main[i], col="grey", border="white", xlab="", ylab="")
  abline(v=quantile(unlist(boot_test$boot_test[i]), probs=c(0.05, 0.95)), col=c("blue", "green"))
  abline(v=0, col=c("black"))
  box()
}

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

# PA and UA
PA.ASD  <- do.call("rbind", boot_test$fit$PA.ASD); PA.ASD <- na.omit(PA.ASD)
PA.AISA <- do.call("rbind", boot_test$fit$PA.AISA); PA.AISA <- na.omit(PA.AISA)
UA.ASD  <- do.call("rbind", boot_test$fit$UA.ASD); UA.ASD <- na.omit(UA.ASD)
UA.AISA <- do.call("rbind", boot_test$fit$UA.AISA); UA.AISA <- na.omit(UA.AISA)

## plot
pdf(file = "Figures/confMatrix_Leaf.pdf", width=10, height=10)
m <- matrix(c(1,2,3,4,5,1,2,3,4,6,7,7,7,7,7), nrow=5, ncol=3)
layout(m)
par(mar = c(3,6,2,0))
# PA 
i=1
beanplot( PA.ASD[,i], PA.AISA[,i],  PA.ASD[,i+1], PA.AISA[,i+1], PA.ASD[,i+2], PA.AISA[,i+2], PA.ASD[,i+3], PA.AISA[,i+3], PA.ASD[,i+4], PA.AISA[,i+4],
          PA.ASD[,i+5], PA.AISA[,i+5],  PA.ASD[,i+6], PA.AISA[,i+6], PA.ASD[,i+7], PA.AISA[,i+7], PA.ASD[,i+8], PA.AISA[,i+8], PA.ASD[,i+9], PA.AISA[,i+9],
          PA.ASD[,i+10], PA.AISA[,i+10],  PA.ASD[,i+11], PA.AISA[,i+11], PA.ASD[,i+12], PA.AISA[,i+12], PA.ASD[,i+13], PA.AISA[,i+13], PA.ASD[,i+14], PA.AISA[,i+14],
          PA.ASD[,i+15], PA.AISA[,i+15],  PA.ASD[,i+16], PA.AISA[,i+16], PA.ASD[,i+17], PA.AISA[,i+17], PA.ASD[,i+18], PA.AISA[,i+18], PA.ASD[,i+19], PA.AISA[,i+19],
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", names=seq(i,i+19,1),
          ll = 0, side = "b", log="", main = "PA", ylab = "", yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1)

i=21
beanplot( PA.ASD[,i], PA.AISA[,i],  PA.ASD[,i+1], PA.AISA[,i+1], PA.ASD[,i+2], PA.AISA[,i+2], PA.ASD[,i+3], PA.AISA[,i+3], PA.ASD[,i+4], PA.AISA[,i+4],
          PA.ASD[,i+5], PA.AISA[,i+5],  PA.ASD[,i+6], PA.AISA[,i+6], PA.ASD[,i+7], PA.AISA[,i+7], PA.ASD[,i+8], PA.AISA[,i+8], PA.ASD[,i+9], PA.AISA[,i+9],
          PA.ASD[,i+10], PA.AISA[,i+10],  PA.ASD[,i+11], PA.AISA[,i+11], PA.ASD[,i+12], PA.AISA[,i+12], PA.ASD[,i+13], PA.AISA[,i+13], PA.ASD[,i+14], PA.AISA[,i+14],
          PA.ASD[,i+15], PA.AISA[,i+15],  PA.ASD[,i+16], PA.AISA[,i+16], PA.ASD[,i+17], PA.AISA[,i+17], PA.ASD[,i+18], PA.AISA[,i+18], PA.ASD[,i+19], PA.AISA[,i+19],PA.ASD[,i+20], PA.AISA[,i+20],
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", names=seq(i,i+20,1),
          ll = 0, side = "b", log="", main = "PA", ylab = "", yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1)

# UA 
i=1
beanplot( UA.ASD[,i], UA.AISA[,i],  UA.ASD[,i+1], UA.AISA[,i+1], UA.ASD[,i+2], UA.AISA[,i+2], UA.ASD[,i+3], UA.AISA[,i+3], UA.ASD[,i+4], UA.AISA[,i+4],
          UA.ASD[,i+5], UA.AISA[,i+5],  UA.ASD[,i+6], UA.AISA[,i+6], UA.ASD[,i+7], UA.AISA[,i+7], UA.ASD[,i+8], UA.AISA[,i+8], UA.ASD[,i+9], UA.AISA[,i+9],
          UA.ASD[,i+10], UA.AISA[,i+10],  UA.ASD[,i+11], UA.AISA[,i+11], UA.ASD[,i+12], UA.AISA[,i+12], UA.ASD[,i+13], UA.AISA[,i+13], UA.ASD[,i+14], UA.AISA[,i+14],
          UA.ASD[,i+15], UA.AISA[,i+15],  UA.ASD[,i+16], UA.AISA[,i+16], UA.ASD[,i+17], UA.AISA[,i+17], UA.ASD[,i+18], UA.AISA[,i+18], UA.ASD[,i+19], UA.AISA[,i+19],
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", names=seq(i,i+19,1),
          ll = 0, side = "b", log="", main = "UA", ylab = "Accurracy [0-1]", yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1)

i=21
beanplot( UA.ASD[,i], UA.AISA[,i],  UA.ASD[,i+1], UA.AISA[,i+1], UA.ASD[,i+2], UA.AISA[,i+2], UA.ASD[,i+3], UA.AISA[,i+3], UA.ASD[,i+4], UA.AISA[,i+4],
          UA.ASD[,i+5], UA.AISA[,i+5],  UA.ASD[,i+6], UA.AISA[,i+6], UA.ASD[,i+7], UA.AISA[,i+7], UA.ASD[,i+8], UA.AISA[,i+8], UA.ASD[,i+9], UA.AISA[,i+9],
          UA.ASD[,i+10], UA.AISA[,i+10],  UA.ASD[,i+11], UA.AISA[,i+11], UA.ASD[,i+12], UA.AISA[,i+12], UA.ASD[,i+13], UA.AISA[,i+13], UA.ASD[,i+14], UA.AISA[,i+14],
          UA.ASD[,i+15], UA.AISA[,i+15],  UA.ASD[,i+16], UA.AISA[,i+16], UA.ASD[,i+17], UA.AISA[,i+17], UA.ASD[,i+18], UA.AISA[,i+18], UA.ASD[,i+19], UA.AISA[,i+19],UA.ASD[,i+20], UA.AISA[,i+20],
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median",  names=seq(i,i+20,1),
          ll = 0, side = "b", log="", main = "UA", ylab = "", yaxs = "i",cex.lab=1.3, cex.axis=1.3, las=1)


# OA and Kappa
beanplot( unlist(boot_test$fit$OA.ASD), unlist(boot_test$fit$OA.AISA), 
          unlist(boot_test$fit$kappa.ASD), unlist(boot_test$fit$kappa.AISA), 
          col = list("black", "gray"), border = NA, innerboerder=NA, beanlines="median", 
          ll = 0, side = "b", log="", main = "", names=c("OA", "Kappa"), ylab = "", ylim=c(0.4,0.9),
          yaxs = "i", cex.lab=1.3, cex.axis=1.3, las=1)
text(c(1,2),y=0.85, labels="*", cex=2.5)

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("center", legend=c("Full range", "Until red edge"), fill=c("black", "gray"), bty="n", cex=1.5)

# legend
names = gsub("Class: ", "",  colnames(PA.ASD))
names = gsub("_", " ", names)
names = paste(seq(1,41,1), names, sep=".- ")
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("left", title="Species", legend=names, bty = "n", ncol=1, cex=1.5, text.font=3)

dev.off()
