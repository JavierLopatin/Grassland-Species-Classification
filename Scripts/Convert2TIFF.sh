# Created by Javier Lopatin (javierlopatin@gmail.com)
#
# A simple script to use gdal_translate to convert files 
# to the GTiff format.
#

# Inputs:
# $1 is the input directory 
# $2 is the output directory
# $3 is the input files extension

FILES=$1/*.$3
for f in $FILES
do
  echo "Processing $f file..."
  filename=`basename ${f} .${3}`
  echo "Output: ${2}/${filename}.tif"
  gdal_translate -of GTiff ${f} ${2}/${filename}.tif
  gdalcalcstats ${2}/${filename}.tif -ignore 0
done
