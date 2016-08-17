# Created by Javier Lopatin (javierlopatin@mac.com)
#
# A simple script to use gdal_edit.py change images heather to set files 
# to the real resolution.
# This example is set to work with all TIF imagery of a folder
#

# Inputs:
# $1 is the input directory 
# $2 is the input files extension
# $3 is the output resolution

# e.g.: sh SetResolution.sh ~/Documentos/Spp tif 0.0025

FILES=$1/*.$2
for f in $FILES
do
  echo "Processing $f file..."
  filename="basename ${f}.$2"
  python gdal_edit.py -tr $3 $3 ${f}
done
