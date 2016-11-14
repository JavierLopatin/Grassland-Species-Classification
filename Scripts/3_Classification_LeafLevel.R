##############################################################################
## R-Script - 3_Classification_LeafLevel.R                                  ##
## author: Javier Lopatin                                                   ##
## mail: javierlopatin@gmail.com                                            ##  
##                                                                          ##
## description: 
##
## Manuscript: 
##
## last changes: 
##                                                                          ##
##############################################################################

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
lines(svrcf)

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
fitAISA <- tunningModels(hyperAISA@data$Species, hyperAISA$spc, hyperAISA@wavelength)
save(fitAISA, file="fitAISALeaf.Rdata")

# ASD full range
fitASD <-  tunningModels(hyperASD@data$Species, hyperASD$spc, hyperASD@wavelength)
save(fitASD, file="fitASDLeaf.Rdata")

############################################
### Growth form analysis with full range ###
############################################

ASD <- data.frame( classes=hyperASD@data$classes, hyperASD$spc)

graminoids_varImport_ASD <- subset(ASD, classes == "Grass_sp9" | classes == "Nardus_stricta" 
                               | classes == "Grass_Sp_23"   | classes == "Setaria_pumila" 
                               | classes == "Elymus_repens" | classes == "Echinochloa_crus-galli"
                               | classes == "Panicum_capillare")
graminoids_varImport_ASD$classes <- factor(graminoids_varImport_ASD$classes)

fobs_varImport_ASD <- subset(ASD, classes == "Prunella_vulgaris" | classes == "Sp_2" 
                         | classes == "Hypochaeris_radicata" | classes == "Trifolium_pratense" 
                         | classes == "Trifolium_repens" | classes == "Conyza_canadensis" 
                         | classes == "Potentilla_reptans" | classes == "Taraxacum_officinale" 
                         | classes == "Galium_sp" | classes == "Bellis perennis" 
                         | classes == "Glechoma_hederacea" | classes == "Medicago_lupulina" 
                         | classes == "Minuartia_hybrida" | classes == "Plantago_lancelota" 
                         | classes == "Geranium_pusillum" | classes == "Plantago_major" 
                         | classes == "Potentilla_2" | classes == "Achillea_millefolium"
                         | classes == "Oxalis_stricta" | classes == "Medicago_arabica" 
                         | classes == "Echium_vulgare" | classes == "Erigoron_annuus" 
                         | classes == "Senecio_vulgaris" | classes == "Filago_arvensis" 
                         | classes == "Anagallis_arvensis" | classes == "Daucum_carota" 
                         | classes == "Medicago_sativa " | classes == "Rumex_obtusifolius" 
                         | classes == "Convolvulus_sepium" | classes == "Verbena_officinalis" 
                         | classes == "Urtica_dioica" | classes == "Cichorium_intybus" 
                         | classes == "Solidago_gigantea" | classes == "Polygonum_persicaria" 
                         | classes == "Oenothera_biennis" | classes == "Arthemisia_vulgaris" 
                         | classes == "Anthemis_arvensis")
fobs_varImport_ASD$classes <- factor(fobs_varImport_ASD$classes)

### Variable importance

library(vegan)

normalize <- function(x) { (x-min(x))/(max(x)-min(x)) }

# all species 


leaf_mrpp = matrix(NA, ncol=length( hyperASD$spc[1,] ), nrow=2)
for(i in 1:length( hyperASD$spc[1,] )){
  obj_mrpp = mrpp(dat =  hyperASD$spc[,i], grouping = ASD$classes, parallel = 16, 
                  distance = "mahalanobis", permutations = 500)
  leaf_mrpp[1,i] = obj_mrpp$A
  leaf_mrpp[2,i] = obj_mrpp$Pvalue
}

save(leaf_mrpp, file = "bestImp.RData")

# graminoids
leaf_mrpp_gram = matrix(NA, ncol=length( hyperASD$spc[1,] ), nrow=2)
for(i in 1:length( hyperASD$spc[1,] )){
  obj_mrpp = mrpp(dat =  graminoids_varImport_ASD[,2:length(ASD)][i], grouping = graminoids_varImport_ASD$classes, 
                  parallel = 16, distance = "mahalanobis", permutations = 500)
  leaf_mrpp_gram[1,i] = obj_mrpp$A
  leaf_mrpp_gram[2,i] = obj_mrpp$Pvalue
}

save(leaf_mrpp_gram, file = "Gramm_Imp_ASD.RData")

# Forbs
leaf_mrpp_forbs = matrix(NA, ncol=length( hyperASD$spc[1,] ), nrow=2)
for(i in 1:length( hyperASD$spc[1,] )){
  obj_mrpp = mrpp(dat =  leaf_mrpp_forbs[,2:length(ASD)][i], grouping = leaf_mrpp_forbs$classes, 
                  parallel = 16, distance = "mahalanobis", permutations = 500)
  leaf_mrpp_forbs[1,i] = obj_mrpp$A
  leaf_mrpp_forbs[2,i] = obj_mrpp$Pvalue
}

save(leaf_mrpp_forbs, file = "Fobs_Imp_ASD.RData")

plot (hyperASD@wavelength, leaf_mrpp[1,], type="l", main="mrpp")
lines(hyperASD@wavelength, leaf_mrpp_gram[1,], type="l", lty=2, col="blue")
lines(hyperASD@wavelength, leaf_mrpp_forbs[1,], type="l", lty=3, col="red")


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

save.image("ClassLeafLevel.RData")

#########################
### PA and UA results ### 
#########################

pa <- unlist(boot_test$fit$PA.ASD)
ua <- unlist(boot_test$fit$UA.ASD)

PA_UA <- data.frame(Species=names(pa), PA=pa, UA=ua )

# get the median value per specie
dummy_matrix <- matrix( ncol = ncol(PA_UA), nrow = length( unique(PA_UA$Species) ) )
colnames(dummy_matrix) <- colnames(PA_UA)
dummy_matrix[,1] <- as.character( unique(PA_UA$Species) )

for (i in 1:length(unique(PA_UA$Species))){#(levelsNumber-1)
  x = grep( unique(PA_UA$Species)[i], PA_UA$Species )
  sp = PA_UA[x, ]
  pa =  sp$PA 
  ua = sp$UA 
  pa_med = median(pa) 
  ua_med = median(ua)
  dummy_matrix[i,2] <- as.numeric(pa_med)
  dummy_matrix[i,3] <- as.numeric(ua_med)
}

PA_UA <- as.data.frame(dummy_matrix)
PA_UA$PA <- as.numeric( as.character(PA_UA$PA) )
PA_UA$UA <- as.numeric( as.character(PA_UA$UA) )

write.table(PA_UA, file = "Leaf_PA_UA.csv", sep = ",", col.names = T, row.names = F)
