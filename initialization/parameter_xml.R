# ==============================================================================
# Script Name:        Set/Update parameters in iLand project files
# Description:        Updates width, height, random sampling, and external seed
#                     settings in iLand XML project files for different RCP 
#                     scenarios. 
#
# Author:             Hanna Komischke
# Date:               2025-11-18
#
# Input Data:         - width_height.csv 
#                     - iLand project XML files per scenario
# Output Data:        - Updated iLand project XML files
#
# ==============================================================================

# directories and arguments ----------------------------------------------------

setwd("/media/hanna/EXTERNAL_USB/Master/")

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
scenario <- args[2] 


# packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(xml2)
  library(dplyr)
})


# set values --------------------------------------------------------------

# read values
dat <- read.csv(paste0("data/iland/input/",filename,"/width_height.csv"))

# loop over different RCP files
RCP_names <- c("historic","RCP26","RCP45","RCP85")

for (RCP in RCP_names) {
  
  # read XML file
  doc <- read_xml(paste0("iLand/",scenario,"/projectFile_example_",RCP,"_baseline.xml"))
  
  # change width
  xml_find_all(doc, ".//width") %>%
    xml_set_text(as.character(dat$width))
  
  # change height
  xml_find_all(doc, ".//height") %>%
    xml_set_text(as.character(dat$height))
  
  # change time 
  set.seed(161)
  random_sample <- rep(sample(60:79), 8)
  xml_find_all(doc, ".//randomSamplingList") %>%
    xml_set_text(as.character(paste(random_sample, collapse = ",")))
  
  # external seeds
  xml_find_all(doc, ".//externalSeedEnabled") %>%
    xml_set_text("true")
  
  xml_find_all(doc, ".//externalSeedSpecies") %>%
    xml_set_text(paste(c("bepe","fasy","psme","quro","algl","piab","pisy","lade","abal","acps","frex","cabe","saca","potr"), collapse = ", "))
  
  
  # save
  write_xml(doc, paste0("iLand/",scenario,"/projectFile_example_",RCP,"_baseline.xml"))
}

