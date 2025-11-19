#!/bin/bash

# Define file
FILENAME="$1"
# Define scenario        
SCENARIONAME="$2"    

# Base directory
BASE_DIR="/media/hanna/EXTERNAL_USB/Master/data/iland/input/"

# Source directory
SRC_DIR="${BASE_DIR}${FILENAME}/"

# Base destination directory
BASE_DEST_DIR="/media/hanna/EXTERNAL_USB/Master/iLand/"

# Destination directory
DEST_DIR="${BASE_DEST_DIR}${SCENARIONAME}/"


# copy files
cp "${SRC_DIR}"*.asc* "${DEST_DIR}gis/"
cp "${SRC_DIR}environment.txt" "${DEST_DIR}gis/"
cp "${SRC_DIR}tree_init.txt" "${DEST_DIR}init/"
cp "${SRC_DIR}"*_climate.sqlite* "${DEST_DIR}database/"


echo "Input copied to $DEST_DIR"
