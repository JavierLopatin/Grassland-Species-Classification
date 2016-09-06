#!/usr/bin/env python

# Import python modules
import os
from glob import glob
from rsgislib.segmentation import segutils

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
# A temporary path for layers generated during the
# segmentation process. The directory will be created
# and deleted during processing.
tmpPath = "./tmp/"
# The number of clusters (k) in the KMeans.
numClusters = 30
# The minimum object size in pixels.
minObjectSize = 50
# The distance threshold to prevent merging.
# this has been set to an arbitrarily large 
# number to disable this function.  
distThres = 1000000
# The sampling of the input image for input
# to the KMeans
imgSampling = 100
# Maximum number of iterations within KMeans
maxKMeanIter = 200

# RSGISLib function call to execute the segmentation
for i in range(len(rasterList)):
    # The input image for the segmentation
    inputImage = rasterList[i]
    # The output segments (clumps) image
    segmentClumps = "/clumps/"+rasterList[i][:-4]+"clumps.tif"
    # The output clump means image (for visualsation) 
    outputMeanSegments = "/clumps/"+rasterList[i][:-4]+"segments.tif"
    segutils.runShepherdSegmentation(inputImage, segmentClumps, outputMeanSegments, 
                                     tmpPath, "GTiff", False, False, False, numClusters, 
                                     minObjectSize, distThres, None, imgSampling, maxKMeanIter)
