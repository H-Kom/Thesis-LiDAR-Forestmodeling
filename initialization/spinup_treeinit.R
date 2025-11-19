# ==============================================================================
# Script Name:        Prepare initial tree and sapling data for iLand spin-up
# Description:        Generates tree initialization data using a 500-year spin-up
#                     routine. Processes historic output, computes canopy LAI and
#                     height rasters, filters trees and saplings, and creates
#                     the combined tree_init_spinup.txt input file for iLand.
#
# Author:             Hanna Komischke
# Date:               2025-11-18
#
# Input Data:         - Historic iLand output SQLite files (historic_output.sqlite)
#                     - AOI raster (objectid.asc)
# Output Data:        - tree_init_spinup.txt (tree initialization input for iLand)
# ==============================================================================


# directories and arguments -----------------------------------------------

setwd("/media/hanna/EXTERNAL_USB/Master/")

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

# packages ----------------------------------------------------------------

library(RSQLite)
library(dplyr)
library(tidyr)
library(terra)


# functions ---------------------------------------------------------------


get_output <- function(filename) {
  output_list <- list()
  db <- dbConnect(SQLite(), paste0("data/iland/output/",filename,"/historic_output.sqlite"))
  on.exit(dbDisconnect(db), add = TRUE)
  output_names <- dbListTables(db)
  
  for (tbl in output_names) {
    
    cols <- dbListFields(db, tbl)
    
    if ("year" %in% cols) {
      if (grepl("landscape", tbl, ignore.case = TRUE)) {
        query <- paste0("SELECT * FROM ", tbl)
      } else if (grepl("tree|saplingsdetail", tbl, ignore.case = TRUE)) {
        query <- paste0("SELECT * FROM ", tbl, " WHERE year IN (0,500)")
      } else {
        query <- paste0("SELECT * FROM ", tbl)
      }
    } else {
      query <- paste0("SELECT * FROM ", tbl)
    }
    
    result <- dbGetQuery(db, query)
    output_list[[tbl]] <- result
  }
  
  return(output_list)
}

expand_saplings <- function(df) {
  df %>%
    mutate(
      x_index = position %% 50,
      y_index = floor(position / 50)
    ) %>%
    uncount(weights = round(n_represented), .remove = FALSE) %>%
    mutate(
      x_offset = x_index * 2 + runif(n(), 0, 2),
      y_offset = y_index * 2 + runif(n(), 0, 2)
    )
}

# ==============================================================================
# trees ========================================================================
# ==============================================================================

# data --------------------------------------------------------------------
# load trees at start andafter 500 years
hist_output <- get_output(filename)

tree <- hist_output$tree

tree_0 <- tree %>% filter(year==0) 
tree_500 <- tree %>% filter(year==500) 


# LAI raster ------------------------------------------------------------------

# Parameter
cellsize <- 10
radius <- 10  # Moving Window Radius in m

# get range
x_range <- range(tree_0$x)
y_range <- range(tree_0$y)

# create raster
r <- rast(
  xmin = min(tree_0$x), xmax = max(tree_0$x),
  ymin = min(tree_0$y), ymax = max(tree_0$y),
  resolution = cellsize
)

# start trees as vector
trees_vect <- vect(tree_0, geom = c("x","y"))

# LeafArea-Raster
leaf_r <- rasterize(trees_vect, r, field = "leafArea_m2", fun = sum, background = 0)

# Height-Raster
height_r <- rasterize(trees_vect, r, field = "height", fun = median, background = NA)

# circular moving window
window <- focalMat(leaf_r, radius, type = "circle")

# focal
leaf_lai <- focal(leaf_r, w = window, fun = sum, na.policy = "omit") / (cellsize^2)

window <- focalMat(height_r, radius, type = "circle")
window[window > 0] <- 1

suppressWarnings(height_canopy <- focal(height_r, w = window, fun=function(x, ...) max(x, na.rm=TRUE)))
height_canopy[!is.finite(height_canopy)] <- NA

# as dataframe
tree_cover_mw <- c(leaf_lai, height_canopy) %>% 
  as.data.frame(xy = TRUE)

names(tree_cover_mw) <- c("x","y","leaf_area","height_canopy")

# correct coordinates
tree_cover_mw <- tree_cover_mw %>%
  mutate(
    cell_x = floor(x / cellsize),
    cell_y = floor(y / cellsize)
  ) %>%
  select(cell_x, cell_y, leaf_area, height_canopy)

tree_with_cells <- tree_500 %>%  
  mutate(
    cell_x = floor(x / cellsize),
    cell_y = floor(y / cellsize)
  )

# filter leaf_area > 2
valid_cells <- tree_cover_mw %>%
  filter(leaf_area >= 2) %>%
  select(cell_x, cell_y, height_canopy)


# combine and filter for height
trees_new <- tree_with_cells %>%
  inner_join(valid_cells, by = c("cell_x", "cell_y"))  %>% filter(height <= (height_canopy)) 

# filter
trees_filtered <- trees_new %>% select(x,y,dbh,species,height)

# combine 
tree_init_new <- tree_0 %>% select(names(trees_filtered)) %>% rbind(trees_filtered)


# ==============================================================================
# saplings =====================================================================
# ==============================================================================

# get saplings
sap_all <- hist_output$saplingdetail

# filter
sap <- sap_all %>% filter( height >= 1,year ==500)

# get objectid information
objectid <-rast(paste0("data/iland/input/",filename,"/objectid.asc"))
ru_df <- as.data.frame(terra::xyFromCell(objectid, 1:ncell(objectid)))
ru_df$id <- terra::values(objectid)
colnames(ru_df) <- c("ru_x", "ru_y", "ru")

# correct coordinates
sap_expanded <- expand_saplings(sap)  %>% left_join(ru_df, by = "ru")%>% na.omit() %>% mutate(
  x = ru_x - 50 + x_offset,
  y = ru_y - 50 + y_offset
) %>% select(x,y,dbh,species,height)


# save --------------------------------------------------------------------

# combine and save
trees_spa <- rbind(tree_init_new,sap_expanded)


write.table(trees_spa,paste0("data/iland/input/",filename,"/tree_init_spinup.txt"), sep = ";", row.names = F)

