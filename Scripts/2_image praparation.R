## R-Script - Preparing images to process
## author: Javier Lopatin 
## mail: javierlopatin@gmail.com
## Manuscript: 
## last changes: 

pkgs<-c("rgdal", "raster", "doParallel")
lapply(pkgs, require, character.only=T)


# set working direction
home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp/"

setwd(home)

# make a list of all .tif files
rast_list = list.files("TIFF", pattern = ".tif") 
# import rasters
rasterlist <- list()
for (i in 1:length(rast_list)){
  ras <- stack(paste(home, "/TIFF/", rast_list[i], sep=""))
  rasterlist[[i]] <- ras
  }
# assign names to the list
a = gsub("0605_samples_", x=rast_list, replacement ="")
a = gsub(".tif", x=a, replacement ="")
a = gsub("-", x=a, replacement ="_")
names(rasterlist) <- a

# assign band names to rasters
# 61 band spectra
spectra <- c("398nm", "407nm", "415nm", "424nm", "432nm", "441nm", "450nm", "459nm",
             "468nm", "477nm", "486nm", "495nm", "504nm", "513nm", "522nm", "531nm",
             "540nm", "550nm", "558nm", "568nm", "577nm", "587nm", "596nm", "605nm",
             "615nm", "624nm", "633nm", "643nm", "652nm", "661nm", "671nm", "680nm",
             "690nm", "699nm", "708nm", "717nm", "727nm", "736nm", "746nm", "755nm",
             "765nm", "775nm", "784nm", "794nm", "803nm", "813nm", "822nm", "832nm",
             "842nm", "851nm", "861nm", "870nm", "880nm", "890nm", "899nm", "908nm", 
             "918nm", "928nm", "937nm", "947nm", "957nm")
# first 10 MNF compomponents
MNF <- list()
for (i in 1:10){
  r <- paste("MNF", i, sep="")
  MNF[i] <- r 
}
# gray level co-occurence matrix of the MNF components
var <- list()
homo <- list()
contrast <- list()
diss <- list()
entropy <- list()
sec_mom <- list()
for (i in 1:10){
  var[i] <- paste("Variance", i, sep="")
  homo[i] <- paste("Homogeneity", i, sep="")
  contrast[i] <- paste("Contrast", i, sep="")
  diss[i] <- paste("Dissimilarity", i, sep="")
  entropy[i] <- paste("Entropy", i, sep="")
  sec_mom[i] <- paste("Second_Moment", i, sep="")
}
GLCM <- c(unlist(var), unlist(homo), unlist(contrast), unlist(diss),
          unlist(entropy), unlist(sec_mom)) 

bandNames <- c(spectra, unlist(MNF), GLCM)

#set names
for (i in 1:length(rast_list)){
  names(rasterlist[[i]]) <- bandNames
}

## Extract species training areas
# make a list of the shepfiles with the trainning areas
setwd(home)
shp_list = list.files("shp/points", pattern = ".shp")
shp_list = gsub('.{4}$', '', shp_list)
# import shepfiles
setwd(shpDir)
shapefiles <- list()
for (i in 1:length(shp_list)){
  shp <- readOGR(dsn=".", layer= shp_list[i])
  shapefiles[[i]] <- shp
}

# assign names to the list
names(shapefiles) <- shp_list

## extract values per specie
###
# initialize parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

x1 <- data.frame(extract(rasterlist$plot2, shapefiles$Background_p2))
x2 <- data.frame(extract(rasterlist$plot3, shapefiles$Background_p3))
Background <- data.frame(rbind(x1, x2))
Background[["spp"]] <- rep("Background", length(Background[,1]))
Background[["class"]] <- rep("Back", length(Background[,1]))
### spp 1
x1 <- extract(rasterlist$spp_1_8_1, shapefiles$spp1_1)
x2 <- extract(rasterlist$spp_1_8_2, shapefiles$spp1_2)
Potentilla_neumanniana <- data.frame(rbind(x1, x2))
Potentilla_neumanniana[["spp"]] <- rep("Potentilla_neumanniana", length(Potentilla_neumanniana[,1]))
Potentilla_neumanniana[["class"]] <- rep("spp01", length(Potentilla_neumanniana[,1]))
### spp 2
x1 <- extract(rasterlist$spp_1_8_1, shapefiles$spp2_1)
x2 <- extract(rasterlist$spp_1_8_2, shapefiles$spp2_2)
Glechoma_hederacea <- data.frame(rbind(x1, x2))
Glechoma_hederacea[["spp"]] <- rep("Glechoma_hederacea", length(Glechoma_hederacea[,1]))
Glechoma_hederacea[["class"]] <- rep("spp02", length(Glechoma_hederacea[,1]))
### spp 3-5
x1 <- extract(rasterlist$spp_1_8_1, shapefiles$spp3_1)
x2 <- extract(rasterlist$spp_1_8_2, shapefiles$spp3_2)
x3 <- extract(rasterlist$plot3, shapefiles$spp3_p3)
Festuca_pratensis <- data.frame(rbind(x1, x2, x3))
Festuca_pratensis[["spp"]] <- rep("Festuca_pratensis", length(Festuca_pratensis[,1]))
Festuca_pratensis[["class"]] <- rep("spp03", length(Festuca_pratensis[,1]))
### spp 6
x1 <- extract(rasterlist$spp_1_8_1, shapefiles$spp6_1)
x2 <- extract(rasterlist$spp_1_8_2, shapefiles$spp6_2)
x2 <- extract(rasterlist$plot3, shapefiles$spp6_p3)
Plantago_lancelota <- data.frame(rbind(x1, x2, x3))
Plantago_lancelota[["spp"]] <- rep("Plantago_lancelota", length(Plantago_lancelota[,1]))
Plantago_lancelota[["class"]] <- rep("spp06", length(Plantago_lancelota[,1]))
### spp 7
x1 <- extract(rasterlist$spp_1_8_1, shapefiles$spp7_1)
x2 <- extract(rasterlist$spp_1_8_2, shapefiles$spp7_2)
Bellis_perennis <- data.frame(rbind(x1, x2))
Bellis_perennis[["spp"]] <- rep("Bellis_perennis", length(Bellis_perennis[,1]))
Bellis_perennis[["class"]] <- rep("spp07", length(Bellis_perennis[,1]))
### spp 7 flower
x1 <- extract(rasterlist$spp_1_8_1, shapefiles$spp7_1_flower)
x2 <- extract(rasterlist$spp_1_8_2, shapefiles$spp7_2_flower)
Bellis_perennis_flower <- data.frame(rbind(x1, x2))
Bellis_perennis_flower[["spp"]] <- rep("Bellis_perennis_flower", length(Bellis_perennis_flower[,1]))
Bellis_perennis_flower[["class"]] <- rep("spp07_f", length(Bellis_perennis_flower[,1]))
### spp 8
x1 <- extract(rasterlist$spp_1_8_1, shapefiles$spp8_1)
x2 <- extract(rasterlist$spp_1_8_2, shapefiles$spp8_2)
x3 <- extract(rasterlist$plot3, shapefiles$spp8_p3)
Taraxacum_officinale <- data.frame(rbind(x1, x2, x3))
Taraxacum_officinale[["spp"]] <- rep("Taraxacum_officinale", length(Taraxacum_officinale[,1]))
Taraxacum_officinale[["class"]] <- rep("spp08", length(Taraxacum_officinale[,1]))
### spp 8 flower
x1 <- extract(rasterlist$spp_1_8_1, shapefiles$spp8_1_flower)
x2 <- extract(rasterlist$spp_1_8_2, shapefiles$spp8_2_flower)
x2 <- extract(rasterlist$plot3, shapefiles$spp8_flower_p3)
Taraxacum_officinale_flower <- data.frame(rbind(x1, x2, x3))
Taraxacum_officinale_flower[["spp"]] <- rep("Taraxacum_officinale_flower", length(Taraxacum_officinale_flower[,1]))
Taraxacum_officinale_flower[["class"]] <- rep("spp08_f", length(Taraxacum_officinale_flower[,1]))
### spp 9
Klee <- data.frame(extract(rasterlist$spp_9.12_2, shapefiles$spp9))
Klee[["spp"]] <- rep("Klee", length(Klee[,1]))
Klee[["class"]] <- rep("spp9", length(Klee[,1]))
### spp 12
Moss <- data.frame(extract(rasterlist$spp_9.12_2, shapefiles$spp12))
Moss[["spp"]] <- rep("Moss", length(Moss[,1]))
Moss[["class"]] <- rep("spp12", length(Moss[,1]))
### spp 13
Gallium_sp <- data.frame(extract(rasterlist$spp_13.15.16_1, shapefiles$spp13))
Gallium_sp[["spp"]] <- rep("Gallium_sp", length(Gallium_sp[,1]))
Gallium_sp[["class"]] <- rep("spp13", length(Gallium_sp[,1]))
### spp 14
Rumex_sp <- data.frame(extract(rasterlist$plot2, shapefiles$Spp14_p2))
Rumex_sp[["spp"]] <- rep("Rumex_sp", length(Rumex_sp[,1]))
Rumex_sp[["class"]] <- rep("spp14", length(Rumex_sp[,1]))
### spp 15
Sauergras <- data.frame(extract(rasterlist$spp_13.15.16_1, shapefiles$spp15))
Sauergras[["spp"]] <- rep("Sauergras", length(Sauergras[,1]))
Sauergras[["class"]] <- rep("spp15", length(Sauergras[,1]))
### spp 16
Trifolium_pratense <- data.frame(extract(rasterlist$spp_13.15.16_1, shapefiles$spp16))
Trifolium_pratense[["spp"]] <- rep("Trifolium_pratense", length(Trifolium_pratense[,1]))
Trifolium_pratense[["class"]] <- rep("spp16", length(Trifolium_pratense[,1]))
### spp 17 
Asteraceae <- data.frame(extract(rasterlist$plot2, shapefiles$spp17_p2))
Asteraceae[["spp"]] <- rep("Asteraceae", length(Asteraceae[,1]))
Asteraceae[["class"]] <- rep("spp17", length(Asteraceae[,1]))
### spp 19
Fabaceae_sp1 <-  data.frame(extract(rasterlist$plot3, shapefiles$spp19_p3))
Fabaceae_sp1[["spp"]] <- rep("Fabaceae_sp1", length(Fabaceae_sp1[,1]))
Fabaceae_sp1[["class"]] <- rep("spp19", length(Fabaceae_sp1[,1]))
### spp 20
Kleiner_Odermennig <- data.frame(extract(rasterlist$plot3, shapefiles$spp20_p3))
Kleiner_Odermennig[["spp"]] <- rep("Kleiner_Odermennig", length(Kleiner_Odermennig[,1]))
Kleiner_Odermennig[["class"]] <- rep("spp21", length(Kleiner_Odermennig[,1]))
### spp 22
Ranunculus_acris <- data.frame(extract(rasterlist$plot3, shapefiles$spp22_p3))
Ranunculus_acris[["spp"]] <- rep("Ranunculus_acris", length(Ranunculus_acris[,1]))
Ranunculus_acris[["class"]] <- rep("spp21", length(Ranunculus_acris[,1]))
### spp 23
Achillea_millefolium <- data.frame(extract(rasterlist$plot4, shapefiles$spp_23_p4))
Achillea_millefolium[["spp"]] <- rep("Achillea_millefolium", length(Achillea_millefolium[,1]))
Achillea_millefolium[["class"]] <- rep("spp21", length(Achillea_millefolium[,1]))
### spp 28
Centaurea_jacea <- data.frame(extract(rasterlist$plot5, shapefiles$spp28_p5))
Centaurea_jacea[["spp"]] <- rep("Centaurea_jacea", length(Centaurea_jacea[,1]))
Centaurea_jacea[["class"]] <- rep("spp28", length(Centaurea_jacea[,1]))
### spp 30
x1 <- data.frame(extract(rasterlist$plot5, shapefiles$spp30_p5))
x2 <- data.frame(extract(rasterlist$plot7, shapefiles$spp30_p7))
Fabaceae_sp2 <- rbind(x1, x2) 
Fabaceae_sp2[["spp"]] <- rep("Fabaceae_sp2", length(Fabaceae_sp2[,1]))
Fabaceae_sp2[["class"]] <- rep("spp30", length(Fabaceae_sp2[,1]))
### spp 31
Tragopogon_pratensis <- data.frame(extract(rasterlist$plot5, shapefiles$spp31_p5))
Tragopogon_pratensis[["spp"]] <- rep("Tragopogon_pratensis", length(Tragopogon_pratensis[,1]))
Tragopogon_pratensis[["class"]] <- rep("spp31", length(Tragopogon_pratensis[,1]))
### spp 32
x1 <- data.frame(extract(rasterlist$spp_32_1, shapefiles$spp32_1))
x2 <- data.frame(extract(rasterlist$spp_32_2, shapefiles$spp32_2))
Strawberry <- rbind(x1, x2) 
Strawberry[["spp"]] <- rep("Strawberry", length(Strawberry[,1]))
Strawberry[["class"]] <- rep("spp32", length(Strawberry[,1]))

# stop parallel process
stopCluster(cl)

## list of species
species.all <- list(Background, Potentilla_neumanniana, Glechoma_hederacea, Festuca_pratensis, Plantago_lancelota, Bellis_perennis,
                    Bellis_perennis_flower, Taraxacum_officinale, Taraxacum_officinale_flower, Klee, Moss, Gallium_sp, Rumex_sp, Sauergras, 
                    Trifolium_pratense, Asteraceae, Fabaceae_sp1, Kleiner_Odermennig, Ranunculus_acris,Achillea_millefolium,  Centaurea_jacea, 
                    Fabaceae_sp2, Tragopogon_pratensis, Strawberry)

species.all <- do.call("rbind", species.all)

### Assing species per plot
plot1 <- na.omit(rbind(Background, Potentilla_neumanniana, Glechoma_hederacea, Festuca_pratensis, Plantago_lancelota, Bellis_perennis,
               Bellis_perennis_flower, Taraxacum_officinale, Taraxacum_officinale_flower, Klee, Moss))
plot2 <- na.omit(rbind(Background, Potentilla_neumanniana, Glechoma_hederacea, Festuca_pratensis, Plantago_lancelota, Bellis_perennis,
               Bellis_perennis_flower, Taraxacum_officinale, Taraxacum_officinale_flower, Klee, Moss, Gallium_sp, Rumex_sp, 
               Sauergras, Trifolium_pratense, Asteraceae))
plot3 <- na.omit(rbind(Background, Festuca_pratensis, Plantago_lancelota, Bellis_perennis, Bellis_perennis_flower, 
               Taraxacum_officinale, Taraxacum_officinale_flower, Moss, Fabaceae_sp1, Ranunculus_acris, Tragopogon_pratensis, 
               Achillea_millefolium, Kleiner_Odermennig))
plot4 <- na.omit(rbind(Background, Glechoma_hederacea, Festuca_pratensis, Plantago_lancelota, Taraxacum_officinale, Taraxacum_officinale_flower,
               Klee, Moss, Sauergras, Trifolium_pratense, Asteraceae, Fabaceae_sp1, Achillea_millefolium))
plot5 <- na.omit(rbind(Background, Festuca_pratensis, Plantago_lancelota, Taraxacum_officinale, Taraxacum_officinale_flower, Klee, Sauergras, 
               Trifolium_pratense, Centaurea_jacea, Fabaceae_sp2, Tragopogon_pratensis, Achillea_millefolium))
plot6 <- na.omit(rbind(Background, Festuca_pratensis, Plantago_lancelota, Taraxacum_officinale, Moss, Sauergras,
                       Trifolium_pratense, Strawberry, Achillea_millefolium, Kleiner_Odermennig))
plot7 <- na.omit(rbind(Background, Festuca_pratensis, Plantago_lancelota, Taraxacum_officinale, Taraxacum_officinale_flower, Sauergras, Trifolium_pratense, 
               Asteraceae, Fabaceae_sp1, Fabaceae_sp2, Achillea_millefolium, Kleiner_Odermennig))

plot.all <- list(plot1, plot2, plot3, plot4, plot5, plot6, plot7)
names(plot.all) <- c("plot1", "plot2", "plot3", "plot4", "plot5", "plot6", "plot7")


# export plots
setwd(home)
write.table(plot1, file = "plot1_train.txt", sep = " ")
write.table(plot2, file = "plot2_train.txt", sep = " ")
write.table(plot3, file = "plot3_train.txt", sep = " ")
write.table(plot4, file = "plot4_train.txt", sep = " ")
write.table(plot5, file = "plot5_train.txt", sep = " ")
write.table(plot6, file = "plot6_train.txt", sep = " ")
write.table(plot7, file = "plot7_train.txt", sep = " ")

save.image("Class.RData")

