
## Preprocessing AISA+ images

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


# spectralon values
ref_spectralon = c(2000, 3100, 4300, 5800, 6800, 7800, 8900, 9500, 9600, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800
                   , 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800
                   , 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9700, 9300, 9300, 9200)

ref_GrayPannel = c(1368.071, 1637.420, 1719.419, 1730.734, 1732.168, 1731.766, 1732.232, 1732.717, 1732.277, 1729.486, 1724.009, 1716.584,
                   1706.202, 1693.486, 1684.633, 1681.018, 1674.633, 1666.777, 1665.688, 1668.840, 1672.203, 1669.678, 1662.034, 1655.917,
                   1657.575, 1663.587, 1669.863, 1671.691, 1667.264, 1658.895, 1650.662, 1647.841, 1654.812, 1667.135, 1676.851, 1685.421,
                   1693.352, 1701.121, 1696.755, 1696.767, 1696.670, 1695.870, 1696.382, 1695.727, 1694.088, 1692.767, 1692.079, 1691.755,
                   1691.591, 1692.392, 1692.879, 1692.952, 1693.440, 1693.822, 1696.213, 1695.322, 1693.605, 1692.638, 1693.317, 1693.062,
                   1692.026)

# reference pannel
ref = ref_GrayPannel

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



