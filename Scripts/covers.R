library(raster)

home = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp"


setwd(home)

Iter=10
ObservedSpecies = species
shpDir = "C:/Users/Lopatin/Dropbox/PhD/Grass_single_spp_segmentation/Single_spp/plot_extent"

shpMaskName = "plot_18"
plotNumber = 18
rasterDir = "D:/Sp_Images/plot_18_MNF_PLS_potVal_MNF"

obstainCovers(rasterDir = rasterDir, shpDir = shpDir, ObservedSpecies = species, plotNumber = plotNumber, shpMaskName = shpMaskName, Iter=Iter)

rasterDir = "D:/Sp_Images/plot_18_MNF_PLS_rf_MNF"

obstainCovers(rasterDir = rasterDir, shpDir = shpDir, ObservedSpecies = species, plotNumber = plotNumber, shpMaskName = shpMaskName, Iter=Iter)

rasterDir = "D:/Sp_Images/plot_18_MNF_RF_potVal_MNF"

obstainCovers(rasterDir = rasterDir, shpDir = shpDir, ObservedSpecies = species, plotNumber = plotNumber, shpMaskName = shpMaskName, Iter=Iter)

rasterDir = "D:/Sp_Images/plot_18_MNF_RF_rf_MNF"

obstainCovers(rasterDir = rasterDir, shpDir = shpDir, ObservedSpecies = species, plotNumber = plotNumber, shpMaskName = shpMaskName, Iter=Iter)
rasterDir = "D:/Sp_Images/plot_18_MNF_SVM_potVal_MNF"

obstainCovers(rasterDir = rasterDir, shpDir = shpDir, ObservedSpecies = species, plotNumber = plotNumber, shpMaskName = shpMaskName, Iter=Iter)

rasterDir = "D:/Sp_Images/plot_18_MNF_SVM_rf_MNF"

obstainCovers(rasterDir = rasterDir, shpDir = shpDir, ObservedSpecies = species, plotNumber = plotNumber, shpMaskName = shpMaskName, Iter=Iter)

