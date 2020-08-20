# Script for estimating monthly SNOTEL snowfall fractions
# Uses optimized rain-snow temperature from Jennings et al. (2018)

# Load packages
library(raster)
library(sp)
library(tidyverse)

################################################################################
################################# Import data ##################################
################################################################################

# Daily snotel data
snotel <- readRDS("../data/snotel/processed/snotel_daily_formatted.RDS")

# PRISM tair and tdew
# Only works on KSJ local machine with external drive -- data v big
prism_tair_path <- "/Volumes/files/climate_data/prism/daily/tair/"
prism_tdew_path <- "/Volumes/files/climate_data/prism/daily/tdew/"

# list all the met files
met_files <- data.frame(tair_files = list.files(pattern = "\\.bil$"), 
                        stringsAsFactors = FALSE) %>% 
  mutate(year = as.numeric(substr(ppt_files, start = 24, stop = 27)),
         month = as.numeric(substr(ppt_files, start = 28, stop = 29)),
         wyear = ifelse(month %in% 10:12, year + 1, year))

# Make a vector of the years to be analyzed
years <- sort(unique(met_files$year)) # 