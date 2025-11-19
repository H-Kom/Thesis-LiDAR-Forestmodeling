#!/bin/bash

# Define file
FILENAME="$1"

# Define scenario       
SCENARIONAME="$2"    

# Base directory
BASE_DIR="/media/hanna/EXTERNAL_USB/Master/iLand/"

# Source directory
SRC_DIR="${BASE_DIR}${SCENARIONAME}/output/"

# Base destination directory
BASE_DEST_DIR="/media/hanna/EXTERNAL_USB/Master/data/iland/output/"

# Destination directory
DEST_DIR="${BASE_DEST_DIR}${FILENAME}/"


for file in "${SRC_DIR}"*.sqlite*; do
    [ -e "$file" ] || continue  
    basename=$(basename "$file")

    if [[ "$basename" == RCP* ]]; then
        mv "$file" "${DEST_DIR}${SCENARIONAME}_${basename}"
    else
        mv "$file" "${DEST_DIR}${basename}"
    fi
done



echo "Move Output to $DEST_DIR"
