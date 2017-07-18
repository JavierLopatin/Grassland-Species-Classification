
################################################################################
## R-Script - 1_Preprocessing.R                                               ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##  
##                                                                            ##
## Manuscript: Mapping plant species in mixed grassland communities using     ##
##             close range imaging spectroscopy                               ##
##                                                                            ##
## description: Preprocessing AISA+ images                                    ## 
##                                                                            ##
################################################################################

library(raster)
library(rgdal)

# set working direction
home = "F/Sp_Images/Site1/Raw"
setwd(home)

rasterNames <- rasterListNames(fileExtantion = ".dat", folder = "Raw")
rasterlist <- rasterList(fileExtantion = ".dat", folder = "Raw")

# extract reference reflectance
spectralon <- list()
for(i in 1:length(rasterNames)){
  # select the with/gray reference area
  rast <- rasterlist[[i]]
  x11()
  plotRGB(rast, 31,18,5, stretch="lin")
  ext <- drawExtent()
  a <- extract(rast, ext, fun=mean)
  b <- as.numeric(a)
  spectralon[[i]] <- b
  dev.off()
}

# set the band name to the spectra
AISAbandnames <- c("398nm", "407nm", "415nm", "424nm", "432nm", "441nm", "450nm", "459nm",
               "468nm", "477nm", "486nm", "495nm", "504nm", "513nm", "522nm", "531nm",
               "540nm", "550nm", "558nm", "568nm", "577nm", "587nm", "596nm", "605nm",
               "615nm", "624nm", "633nm", "643nm", "652nm", "661nm", "671nm", "680nm",
               "690nm", "699nm", "708nm", "717nm", "727nm", "736nm", "746nm", "755nm",
               "765nm", "775nm", "784nm", "794nm", "803nm", "813nm", "822nm", "832nm",
               "842nm", "851nm", "861nm", "870nm", "880nm", "890nm", "899nm", "908nm", 
               "918nm", "928nm", "937nm", "947nm", "957nm")


# Spectralon reflectance values
ref = c(2000, 3100, 4300, 5800, 6800, 7800, 8900, 9500, 9600, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800
      , 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800
      , 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9700, 9300, 9300, 9200)


# Plot spectralon spectra
plot(seq(395,955,9.3), unlist(spectralon[1]), type="l", xlab="Wavelength [nm]", ylab="Refectance", ylim=c(0,10000))
for(i in 1:length(spectralon)){  lines(seq(395,955,9.3), unlist(spectralon[i+1]))  }
lines(seq(395,955,9.3), ref, col="red")
legend("topright", c("Observed reference spectra", "Base reference spectra"), lty =1, col=c("black", "red"))

# calibrate relativ deviances with spectralon (between scans)
offset <- matrix(NA, nrow=length(AISAbandnames), ncol=length(spectralon))
for(i in 1:61){
  offset[, i] =  ref / as.numeric(spectralon[[i]])  
}
offset

# create directory to stor BN images
dir.create("Processed")

# Appy correction
setwd(file.path(home, "Processed"))

for(i in 1:length(rasterlist)){
  rast <- rasterlist[[i]]
  # list of the raster bands
  band_list <- stack()
  # loop thru each band of each scan
  for(j in 1:61){
    a <- rast[[j]] * offset[j,i]
    band_list <- addLayer(band_list, a) 
  }  
  # Export the corrected images
  output = paste(rasterNames[i], ".tif", sep="")
  writeRaster(band_list, filename=output, options="INTERLEAVE=BAND", overwrite=TRUE)
}



