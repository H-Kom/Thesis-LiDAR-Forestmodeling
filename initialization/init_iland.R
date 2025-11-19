# ==============================================================================
# Script Name:        preparing input files for iLand model
# Description:        Automated preparation of input files for iLand simulations.
#                     This script processes LiDAR-derived tree crown polygons and
#                     dominant tree species raster data, computes DBH
#                     using allometric models, and generates all necessary
#                     input files for iLand
#
# Author:             Hanna Komischke
# Date:               [2025-11-18]
#
# Input Data:         - LiDAR-derived tree crown polygons (.gpkg)
#                     - Dominant tree species raster (.tif)
#                     - Allometric model coefficients (CSV)
# Output Data:        - width_height.csv
#                     - objectid_crs.tif and objectid.asc
#                     - soid.asc
#                     - environment.txt
#                     - tree_init.txt
#
# ==============================================================================

# directories and arguments  ---------------------------------------------------

setwd("/media/hanna/EXTERNAL_USB/Master/data/")

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

cat("Processing file:", filename, "\n")

result_dir <- paste0("iland/input/",filename,"/")
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}


# packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(sf)
})

# functions ---------------------------------------------------------------

# height to dbh (cm)
height_to_dbh <- function(h,b0,b1,b2){
  d <- b0 + b1 * (h)^b2
  return(d)
}

# ==============================================================================
# data =========================================================================
# ==============================================================================

# allometric coefficients
Parker_coef <- read.csv("TreeSpecies/Parker_dbh_coef.csv")

# aoi ---------------------------------------------------------------------
aoi <- vect(paste0("itcs_ams3d/sachsenwald_ams3d/",filename,".gpkg"))

print("AOI loaded.")

# tree species map --------------------------------------------------------

# load data
blick <- rast("TreeSpecies/Dominant_Species_Class.tif")

# create df with tree species parameter
blick_code <- read.csv("TreeSpecies/tree_species_code.csv")

 
iland_species <- c("Birch","Beech","Douglas fir","Oak","Alder","Spruce","Pine","Larch","Fir","Maple","Ash","Hornbeam","Willow","Poplar")
iland_species_code <- c("bepe","fasy","psme","quro","algl","piab","pisy","lade","abal","acps","frex","cabe","saca","potr")

iland_code <- data.frame("species" = iland_species,"code"= iland_species_code)

allo_species <- iland_code %>% 
  left_join(Parker_coef, by = "species") %>%
  left_join(blick_code, by = "species") %>%
  mutate(ID = case_when(
    is.na(ID) & species %in% c("Maple","Ash","Hornbeam")   ~ 16,
    is.na(ID) & species %in% c("Willow","Poplar")      ~ 17,
    TRUE                                               ~ ID  
  )) %>% rename(Dominant_Species_Class = ID)

print("tree species map loaded.")

# ==============================================================================
# processing ===================================================================
# ==============================================================================

# get width and height  ---------------------------------------------------

# get Bounding Box of aoi
bbox <- round(ext(aoi)+100, digits=-2)

# get width and height
width <- bbox[2]-bbox[1]        ### xml file <width>
height <- bbox[4]-bbox[3]       ### xml file <height>

width_height <- cbind("width"=width,"height"=height)

write.csv(width_height,paste0(result_dir,"width_height.csv"))

print("Widht and Height extracted and saved.")

# environment grid --------------------------------------------------------

## object id ### xml file <environmentGrid>

# create raster
bbox_ras <- rast(extent=bbox, resolution=100)
crs(bbox_ras) <- crs(aoi)

# set ids
values(bbox_ras) <- cells(bbox_ras)
names(bbox_ras) <- "id"

bbox_ras_crs <- bbox_ras

# save raster as tif (with crs)
writeRaster(bbox_ras_crs, paste0(result_dir,"objectid_crs.tif"), overwrite=TRUE)

# remove crs informations
crs(bbox_ras) <-""
ext(bbox_ras) <- c(0, width, 0, height)  # xmin, xmax, ymin, ymax

# save raster
writeRaster(bbox_ras, paste0(result_dir,"objectid.asc"), overwrite=TRUE)

print("Environment Grid saved.")

## soid ### xml file <standGrid>

bbox_ras_10 <- rast(extent=bbox, resolution=10)

# set all values 1
values(bbox_ras_10) <-1

# rasterize aoi
aoi_rast_10 <- rasterize(aoi,bbox_ras_10, touches=T)

# set -2 as "forested outside pixel"
aoi_rast_10[is.na(aoi_rast_10)] <- -2

# remove crs informations
crs(aoi_rast_10) <-""
ext(aoi_rast_10) <- c(0, width, 0, height)  # xmin, xmax, ymin, ymax

# save raster
writeRaster(aoi_rast_10, paste0(result_dir,"soid.asc"), overwrite=TRUE)

print("Stand Grid saved.")

# environment file ---------------------------------------------------------

### xml file <environmentFile>

# extract coords
coords_ras <- as.data.frame(bbox_ras,xy=T) 

# combine data
env_aoi <- coords_ras[,c("id","x","y")]
env_aoi$model.initialization.file <- "init/tree_init.txt"

# set values
env_aoi$model.site.availableNitrogen <- 70
env_aoi$model.site.soilDepth <- 200
env_aoi$model.site.pctSand <- 82
env_aoi$model.site.pctSilt <- 10
env_aoi$model.site.pctClay <- 8

# save
write.table(env_aoi,paste0(result_dir,"environment.txt"), sep = " ", row.names = F)

print("Environment file created and saved.")

# single tree input file --------------------------------------------------

### xml file <initialization><file>init/real_tree_init.txt

# filter
aoi_f <- aoi[aoi$area > 4, ]

#--- Position

# Load the data into an sf object 
crowns_sf <- st_as_sf(aoi_f, wkt = "geometry")  

# Calculate the centroid for each geometry
crowns_sf$centroid <- st_centroid(crowns_sf$geometry)

# Extract the coordinates of the centroid
crowns_sf$centroid_x <- st_coordinates(crowns_sf$centroid)[, 1]
crowns_sf$centroid_y <- st_coordinates(crowns_sf$centroid)[, 2]

print("Crown parameter extracted.")

#-- Tree Spieces

# crop and project
bbox_wgs84 <- project(bbox_ras_crs, crs(blick))
blick_c <- crop(blick, ext(bbox_wgs84)+20)
blick_proj <- project(blick_c,crs(aoi))

# extract values and set code
species_vals_raw <- terra::extract(blick_proj, vect(crowns_sf$centroid)) 
species_vals_1 <- species_vals_raw %>% filter(!Dominant_Species_Class %in% c(16,17)) %>% left_join(allo_species,by = join_by(Dominant_Species_Class))
species_vals_2 <- species_vals_raw %>% filter(Dominant_Species_Class %in% c(16,17)) %>% left_join(allo_species,by = join_by(Dominant_Species_Class), relationship="many-to-many") %>% group_by(ID) %>%
  slice_sample(n = 1) %>%  
  ungroup() 

species_vals_allo <- bind_rows(species_vals_1, species_vals_2)

# join
species_vals <- species_vals_raw %>%
  left_join(species_vals_allo, by = "ID")

print("Species extracted.")


#--- DBH
crowns_sf$species <- species_vals$code
crowns_sf$b0 <- species_vals$b0
crowns_sf$b1 <- species_vals$b1
crowns_sf$b2 <- species_vals$b2

# Convert the sf object to a data frame
crowns_df <- as.data.frame(crowns_sf)

# remove trees without species
crowns_df_true <- crowns_df %>% filter(!is.na(species))

# select cols of interest
trees <- crowns_df_true %>% dplyr::select(centroid_x,centroid_y,Z_tree_top, species,b0,b1,b2)

# set correct names
names(trees) <- c("x","y","height","species","b0","b1","b2")

# correct coords
xmin <- ext(bbox)[1]
ymin <- ext(bbox)[3]
trees$x <- trees$x-xmin
trees$y <- trees$y-ymin

# get dbh
trees$dbh <- pmax(height_to_dbh(trees$height,trees$b0,trees$b1,trees$b2),7)

print("DBH computed.")

# select relevant parameter
trees_clean <- trees %>% select(x,y,dbh,height,species)

# save
write.table(trees_clean,paste0(result_dir,"tree_init.txt"), sep = ";", row.names = F)

print("Inventory file saved.")
