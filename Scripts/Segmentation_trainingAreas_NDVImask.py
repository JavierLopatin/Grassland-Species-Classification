#!/usr/bin/env python

# Import python modules
import os, rasterio, rsgislib
from numpy import multiply
from glob import glob
from rsgislib.segmentation import segutils
from rsgislib import vectorutils, imagecalc

def ApplyMaskNDVI(NDVI, imgMask, outMask):
    """
    Function to mask out NDVI areas under 0.2
    """
    # open NDVI
    r = rasterio.open(NDVI)
    r2 = r.read()
    # open image to appy the mask
    img = rasterio.open(imgMask)
    img2 = img.read()
    # prepare mask for NDVI under 0.2
    r2[r2 <= 0.2] = 0
    r2[r2 > 0.2] = 1
    # apply mask
    mask = multiply(img2, r2)
    mask[mask == 0] = "nan"
    # save data
    new_dataset = rasterio.open(outMask, 'w', driver='GTiff',
               height=img.shape[0], width=img.shape[1],
               count=img2.shape[0], dtype=str(mask.dtype),
               crs=img.crs, transform=img.transform)
    new_dataset.write(mask)
    new_dataset.close()


# Create a list of rasters
rasterList = glob("*.tif")
# Delete all images of the plots and leav only the ones from the calibration samples
rasterList = [x for x in rasterList if "plot" not in x]

# Create folders to store results if thay do no exist
if not os.path.exists("clumps"):
    os.makedirs("clumps")

if not os.path.exists("shp"):
    os.makedirs("shp")

##################### Perform Segmentation #####################

for i in range(len(rasterList)):

    # The input image for the segmentation
    inputImage = rasterList[i]
    # The output segments (clumps) image
    clumpsFile = "clumps/"+rasterList[i][:-4]+"_Clump.kea"
    # The output clump means image (for visualsation)
    meanImage = "clumps/"+rasterList[i][:-4]+"_Segments.kea" 
    # The output shapefile
    shapeOut = "shp/"+rasterList[i][:-4]+".shp"

    # run segmentation
    segutils.runShepherdSegmentation(inputImage, clumpsFile,
                    meanImage, numClusters=100, minPxls=1)

    # NDVI
    NDVI = "clumps/"+rasterList[i][:-4]+"_NDVI.kea" 
    imagecalc.calcNDVI(meanImage, 31, 43, NDVI)
    # mask out everything below NDVI 0.2
    maskNDVI = "clumps/"+rasterList[i][:-4]+"_maskNDVI.tif" 
    ApplyMaskNDVI(NDVI, clumpsFile, maskNDVI)

    # run polygonization
    vectorutils.polygoniseRaster(clumpsFile, shapeOut)