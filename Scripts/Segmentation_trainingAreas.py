#!/usr/bin/env python

"""
Python scritp to perform k-mean base segmentation to 
automatic delineation of homogeneous training areas 
for classification.

The segmentation algorithm is based on:
     
"""

# Import python modules
import os
from glob import glob
from rsgislib.segmentation import segutils
from rsgislib import vectorutils

# Create a list of rasters
rasterList = glob("*.tif")
# Delete all images of the plots and leav only the ones from the calibration samples
rasterList = [x for x in rasterList if "plot" not in x]

# Create folders to store results if thay do no exist
if not os.path.exists("temp"):
    os.makedirs("temp")

if not os.path.exists("shp"):
    os.makedirs("shp")

##################### Perform Segmentation #####################

for i in range(len(rasterList)):

    # The input image for the segmentation
    inputImage = rasterList[i]
    # The output segments (clumps) image
    clumpsFile = "temp/"+rasterList[i][:-4]+"_Clump.kea"
    # The output clump means image (for visualsation)
    meanImage = "temp/"+rasterList[i][:-4]+"_Segments.kea" 
    # The output shapefile
    shapeOut = "shp/"+rasterList[i][:-4]+".shp"

    # run segmentation
    segutils.runShepherdSegmentation(inputImage, clumpsFile,
                    meanImage, numClusters=100, minPxls=1)

    # run polygonization
    vectorutils.polygoniseRaster(clumpsFile, shapeOut)