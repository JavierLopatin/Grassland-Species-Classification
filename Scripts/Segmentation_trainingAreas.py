#!/usr/bin/env python

"""
<<<<<<< HEAD
=======
########################################################################################################
>>>>>>> 96cce89f0ae51338f8a2442e52f6222b2df7bf22
Python scritp to perform k-mean base segmentation to 
automatic delineation of homogeneous training areas 
for classification.

<<<<<<< HEAD
The segmentation algorithm is based on:
     
=======
Author: Javier Lopatin
Email: javierlopatin@gmail.com
Date: 09/08/2016
Version: 1.0

The segmentation algorithm is based on:
   Clewley, D.; Bunting, P.; Shepherd, J.; Gillingham, S.; Flood, N.; Dymond, J.; Lucas, R.; 
   Armston, J.; Moghaddam, M. A Python-Based Open Source System for Geographic Object-Based Image 
   Analysis (GEOBIA) Utilizing Raster Attribute Tables. Remote Sensing 2014, 6, 6111-6135.
   http://www.mdpi.com/2072-4292/6/7/6111   
########################################################################################################
>>>>>>> 96cce89f0ae51338f8a2442e52f6222b2df7bf22
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
<<<<<<< HEAD
    vectorutils.polygoniseRaster(clumpsFile, shapeOut)
=======
    vectorutils.polygoniseRaster(clumpsFile, shapeOut)
>>>>>>> 96cce89f0ae51338f8a2442e52f6222b2df7bf22
