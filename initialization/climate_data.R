# ==============================================================================
# Script Name:  create climate input data for Iland model
# Description:  Automated preparation of climate input data for iLand
#               simulations. This script processes both historic and future
#               climate datasets, crops and resamples them to the area of
#               interest (AOI), calculates vapor pressure deficit (VPD), and
#               saves the resulting data into SQLite databases for use in iLand.
#               Finally, it updates the environment.txt file with the generated
#               climate tables.
#
# Author:       Hanna Komischke
# Date:         2025-11-18
#
# Input Data:   - LiDAR-derived AOI raster (objectid_crs.tif)
#               - Historic climate NetCDF files (tasmin, tasmax, tas_, pr,
#                 hurs, rsds)
#               - Future climate NetCDF files for RCP scenarios (tasmin, tasmax,
#                 tas_, pr, huss, ps, rsds)
#
# Output Data:  - SQLite climate databases per scenario:
#                 historic_climate.sqlite, RCP26_climate.sqlite,
#                 RCP45_climate.sqlite, RCP85_climate.sqlite
#               - Updated environment.txt file with model.climate.tableName
#                 entries
# ==============================================================================

# directories and arguments  ---------------------------------------------------

setwd("/media/hanna/EXTERNAL_USB/Master/data/")

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

cat("Processing file:", filename, "\n")

result_dir <- paste0("iland/input/",filename,"/")


# packages ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(RSQLite)
  library(dplyr)
})
# functions ---------------------------------------------------------------

load_var_for_aoi_historic <- function(RCP,aoi_path, var_name) {
  # List files
  nc_files <-
    list.files(paste0("climate/",RCP,"/DWD/"),pattern = var_name,
               recursive = TRUE,
               full.names = TRUE)
  
  r_list <- list()
  
  for (i in seq_along(nc_files)) {
    nc_path <- nc_files[i]
    r <- rast(nc_path)
    
    # dates as names
    names(r) <- format(time(r), "%Y-%m-%d")
    
    if (i == 1) {
      # load AOI
      rast_objectid <- rast(aoi_path)
      
      # project
      rast_objectid_proj <- project(rast_objectid, crs(r))
      bbox_proj <- ext(rast_objectid_proj)
      bbox <- ext(rast_objectid)
    }
    
    #crop
    r_c <- crop(r, ext(bbox_proj)+1000, snap = "out")
    r_list[[i]] <- resample(project(r_c, crs(rast_objectid)),rast_objectid, method="bilinear")
  }
  
  r_stack <- c(rast_objectid, rast(r_list))
  
  df <- as.data.frame(values(r_stack))
  
  
  return(df)
}

load_var_for_aoi <- function(RCP,aoi_path, var_name) {
  # List files
  nc_files <-
    list.files(paste0("climate/",RCP,"/germany_25832/"),pattern = var_name,
               recursive = TRUE,
               full.names = TRUE)
  
  r_list <- list()
  
  for (i in seq_along(nc_files)) {
    nc_path <- nc_files[i]
    suppressWarnings(r <- rast(nc_path))
    
    # dates as names
    names(r) <- format(time(r), "%Y-%m-%d")
    
    if (i == 1) {
      # load AOI
      rast_objectid <- rast(aoi_path)
      crs(rast_objectid) <- crs(r)
    }
    
    #crop
    r_c <-   suppressWarnings(crop(r, ext(rast_objectid)+ 10000, snap="out"))
    
    #resample
    r_list[[i]] <- resample(r_c, rast_objectid, method="bilinear")
  }
  
  r_stack <- c(rast_objectid, rast(r_list))
  
  df <- as.data.frame(values(r_stack))
  
  
  return(df)
}

#### functions for VPD ####

# function(Tetens - formula)
calc_vpd <- function(T, RH) {
  #  if(T >= 0){
  es <- 0.6108 * exp((17.27 * T) / (T + 237.3))  # in kPa
  #  }else{
  #  es <- 0.6108 * exp((21.875 * T) / (T + 265.5))  # in kPa
  #  }
  
  ea <- es * (RH / 100)
  vpd <- es - ea  
  
  # Set negative values to zero
  vpd[vpd < 0] <- 0
  return(vpd)
}

calc_vpd_SH <- function(T, SH, P) {
  # T in Kelvin
  # SH = specific humidity in kg/kg
  # P in Pa (air pressure)
  
  # 1. Convert temperature from Kelvin to Celsius
  T_C <- T - 273.15
  
  # 2. Convert pressure from Pa to kPa
  P_kPa <- P / 1000
  
  # 3. Saturation vapor pressure (es) in kPa
  #  if(T_C >= 0){
  es <- 0.6108 * exp((17.27 * T_C) / (T_C + 237.3))  # in kPa
  #  }else{
  #    es <- 0.6108 * exp((21.875 * T) / (T + 265.5))  # in kPa
  #  }
  # 4. Actual vapor pressure (ea) in kPa
  ea <- (SH * P_kPa) / (0.622 + 0.378 * SH)
  
  # 5. Vapor Pressure Deficit (VPD) = es - ea
  vpd <- es - ea
  
  # Set negative values to zero
  vpd[vpd < 0] <- 0
  return(vpd)
}

# ==============================================================================
# future climate ===============================================================
# ==============================================================================

# source : https://cds.climate.copernicus.eu/datasets/projections-cordex-domains-single-levels?tab=download

RCP_names <- c("RCP26","RCP45","RCP85")

for (RCP in RCP_names) {
  
  aoi_path <- paste0(result_dir,"objectid_crs.tif")
  
  df_tasmin <- load_var_for_aoi(RCP,aoi_path, "tasmin")  # 2m temperature K
  df_tasmax <-
    load_var_for_aoi(RCP,aoi_path, "tasmax") # Maximum 2m temperature in the last 24 hours K
  df_tas <-
    load_var_for_aoi(RCP,aoi_path, "tas_") # Minimum 2m temperature in the last 24 hours
  print("Temperature data loaded.")
  df_pr <-
    load_var_for_aoi(RCP,aoi_path, "pr") # Mean precipitation flux kg.m-2.s-1
  print("Precipitation data loaded.")
  df_huss <-
    load_var_for_aoi(RCP,aoi_path, "huss") # 2m surface specific humidity Dimensionsless
  print("Humidity data loaded.")
  df_ps <- load_var_for_aoi(RCP,aoi_path, "ps") # Surface pressure Pa
  print("Pressure data loaded.")
  df_rsds <-
    load_var_for_aoi(RCP,aoi_path, "rsds") # Surface solar radiation downwards W.m-2
  print("Radiation data loaded.")
  print(paste(RCP,"data loaded."))
  
  
  # prepare data --------------------------------------------------------------
  
  df_climates <- df_tas %>%
    mutate(climate = id) %>%
    dplyr::select(id,climate)
  
  dates <- as.Date(names(df_tas)[-1])  # time steps
  
  df_list <- vector("list", nrow(df_climates))
  
  for (i in seq(nrow(df_climates))) {
    df_list[[i]] <- data.frame(
      climate = i,
      year   = as.numeric(format(dates, "%Y")),
      month  = as.numeric(format(dates, "%m")),
      day    = as.numeric(format(dates, "%d")),
      #tas = as.numeric(df_tas[i, -1]),
      min_temp = as.numeric(df_tasmin[i,-1]) - 273.15,
      max_temp = as.numeric(df_tasmax[i,-1]) - 273.15,
      prec = as.numeric(df_pr[i,-1]) * 86400,
      #rsds = as.numeric(df_rsds[i, -1]),
      rad = as.numeric(df_rsds[i,-1]) * 86400 / 1e6,
      # rsds in W/m² → MJ/m²/Tag
      #huss = as.numeric(df_huss[i, -1]),
      vpd = calc_vpd_SH(
        as.numeric(df_tas[i,-1]),
        as.numeric(df_huss[i,-1]),
        as.numeric(df_ps[i,-1])
      )
    )
    
  }
  print(paste(RCP,"data prepared."))
  
  # save as sql database ----------------------------------------------------
  
  # connect with sql
  db <- dbConnect(SQLite(), paste0(result_dir,RCP,"_climate.sqlite"))
  
  # df per id
  for (i in seq_along(df_list)) {
    df <- df_list[[i]]
    
    # names
    table_name <- paste0("climate_", df$climate[1])
    
    # write
    dbWriteTable(db, table_name, df, overwrite = TRUE)
  }
  
  # disconnect
  dbDisconnect(db)
  
  print(paste(RCP,"data saved."))
  
}

# ==============================================================================
# historic climate ===============================================================
# ==============================================================================

RCP <- "historic" 

aoi_path <- paste0(result_dir,"objectid_crs.tif")

df_tasmin <- load_var_for_aoi_historic(RCP,aoi_path, "tasmin")  # 2m temperature K
df_tasmax <-
  load_var_for_aoi_historic(RCP,aoi_path, "tasmax") # Maximum 2m temperature in the last 24 hours K
df_tas <-
  load_var_for_aoi_historic(RCP,aoi_path, "tas_") # Minimum 2m temperature in the last 24 hours

df_pr <-
  load_var_for_aoi_historic(RCP,aoi_path, "pr") # Mean precipitation flux kg.m-2.s-1
df_hurs <-
  load_var_for_aoi_historic(RCP,aoi_path, "hurs") # 2m surface relative humidity 
df_rsds <-
  load_var_for_aoi_historic(RCP,aoi_path, "rsds") # global radiation W.m-2

print(paste(RCP,"data loaded."))

# prepare data --------------------------------------------------------------
df_climates <- df_tas %>%
  mutate(climate = id) %>%
  dplyr::select(id,climate)

dates <- as.Date(names(df_tas)[-1])  # time steps


df_list <- vector("list", nrow(df_climates))

for (i in seq(nrow(df_climates))) {
  df_list[[i]] <- data.frame(
    climate = i, 
    year   = as.numeric(format(dates, "%Y")),
    month  = as.numeric(format(dates, "%m")),
    day    = as.numeric(format(dates, "%d")),
    #tas = as.numeric(df_tas[i, -1]),
    min_temp = as.numeric(df_tasmin[i, -1]),
    max_temp = as.numeric(df_tasmax[i, -1]),
    prec = as.numeric(df_pr[i, -1]),
    #rsds = as.numeric(df_rsds[i, -1]),
    rad = as.numeric(df_rsds[i, -1]) * 86400 / 1e6,  # rsds in W/m² → MJ/m²/Tag
    #hurs = as.numeric(df_hurs[i, -1]),
    vpd = calc_vpd(as.numeric(df_tas[i, -1]),as.numeric(df_hurs[i, -1]))
  )
}
print(paste(RCP,"data prepared."))


# save as sql database ----------------------------------------------------

# connect with sql
db <- dbConnect(SQLite(), paste0(result_dir,RCP,"_climate.sqlite"))

# df per id
for (i in seq_along(df_list)) {
  df <- df_list[[i]]
  
  # names
  table_name <- paste0("climate_", df$climate[1])
  
  # write
  dbWriteTable(db, table_name, df, overwrite = TRUE)
}

# disconnect
dbDisconnect(db)

print(paste(RCP,"data saved."))

# ==============================================================================
# set climate in environment file ==============================================
# ==============================================================================

df_climates$model.climate.tableName <- paste0("climate_",df_climates$climate)
df_climates <- df_climates %>% dplyr::select(-climate)
env_aoi <- read.table(paste0(result_dir,"environment.txt"), 
                      header = TRUE, sep = " ")

env_aoi_climate <- env_aoi %>% left_join(df_climates,by = join_by(id)) 

# save
write.table(env_aoi_climate,paste0(result_dir,"environment.txt"), 
            sep = " ", row.names = FALSE, quote = FALSE)


print("Environment saved.")

