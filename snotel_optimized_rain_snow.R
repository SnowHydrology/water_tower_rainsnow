# Script for estimating monthly SNOTEL snowfall fractions
# Uses optimized rain-snow temperature from Jennings et al. (2018)

# Load packages
library(raster)
library(sp)
library(tidyverse)

################################################################################
################################# Import data ##################################
################################################################################

# Rain-snow threshold produced by binary logistic regression
t50_map_logi <- raster("../data/geospatial/jennings_et_al_2018_file4_temp50_raster.tif")

# Rain-snow threshold produced by linear regression
t50_map_regr <- raster("../data/geospatial/jennings_et_al_2018_file5_temp50_linregr_raster.tif")

# SNOTEL locations
snotel_meta <- read.csv("../data/snotel/SNOTEL_summary.csv") %>% 
  rename(site_id = SNOTEL.ID,
         state = SNOTEL.State,
         site_name = SNOTEL.Name,
         elev_ft = SNOTEL.Elev..ft.,
         lat = SNOTEL.Lat,
         lon = SNOTEL.Lon,
         start = Data.Starting.Date,
         end = Data.Ending.Date) %>% 
  filter(state != "AK")

# List the snotel data files
snotel_path = "../data/snotel/station_data/"
snotel_files <- data.frame(files = list.files(path = snotel_path,
                                              pattern = "*.txt")) %>% 
                             mutate(lat = as.numeric(str_sub(files, 6, 10)),
                                    lon = as.numeric(str_sub(files, 15, 21)))

# Join to metadata by lat-lon
snotel_meta <- left_join(snotel_meta,
                         snotel_files,
                         by = c("lat", "lon"))

# Identify snotel data column names (from PNNL site)
snotel_cols <- c("year", "month", "day", "ppt_in", "tmax_f", "tmin_f", "tavg_f", "swe_in")

# Loop through station metadata, import data, and identify site
snotel <- data.frame()
for(i in 1:length(snotel_meta$site_id)){
  # Import data to temporary files
  tmp <- read.table(paste0(snotel_path, snotel_meta[i, "files"]), 
                    header = FALSE, col.names = snotel_cols) %>% 
    mutate(site_id = snotel_meta[i, "site_id"])
  
  # Bind to existing data
  snotel <- bind_rows(snotel, tmp)
}

################################################################################
#################### Extract rain-snow temperature threshold ###################
################################################################################

# Convert snotel data to spatialpointsdataframe
coords <- select(snotel_meta, lon, lat)
id <- select(snotel_meta, site_id)
crs_snotel <- CRS("+init=epsg:4326")

# Convert obs to a SpatialPointsDataFrame
snotel_sp <- SpatialPointsDataFrame(coords = coords, 
                                 data = id, 
                                 proj4string = crs_snotel)

# Extract from threshold maps
snotel_meta$t50_logi <- raster::extract(t50_map_logi, snotel_sp)
snotel_meta$t50_regr <- raster::extract(t50_map_regr, snotel_sp)

# Merge into single column (only 3 sites needed regression map)
snotel_meta <- snotel_meta %>% 
  mutate(t50 = case_when(is.na(t50_logi) ~ t50_regr,
                         is.na(t50_regr) ~ t50_logi)) %>% 
  select(-t50_regr, -t50_logi) # remove old columns

################################################################################
####################### Compute monthly snowfall fractions #####################
################################################################################

# Conversion factor for degrees F to C
# and inches to mm
f_to_c = function(x){(x - 32) * (5/9)}
in_to_mm = function(x){x * 25.4}

# Join snotel t50 to data
snotel <- left_join(snotel, 
                    select(snotel_meta, site_id, t50), 
                    by = "site_id")

# Convert air temperature and compute snowfall
snotel <- snotel %>% 
  mutate(tavg_c = f_to_c(tavg_f),
         snowfall_in = case_when(tavg_c <= t50 ~ ppt_in,
                              TRUE ~ 0),
         del_swe_in = c(0, diff(swe_in)))

# Compute monthly snowfall fraction per site
snotel_snowfall_fraction <- snotel %>% 
  filter(!is.na(ppt_in) & !is.na(snowfall_in)) %>% 
  group_by(site_id, month, year) %>% 
  summarise(ppt_mm = sum(ppt_in) %>% in_to_mm(), 
            snowfall_mm = sum(snowfall_in) %>% in_to_mm(),
            snowfall_frac = snowfall_mm / ppt_mm, 
            n_obs = n())

# Compute monthly melt per site
snotel_melt <- snotel %>% 
  filter(del_swe_in < 0) %>% 
  group_by(site_id, month, year) %>% 
  summarise(melt_mm = sum(del_swe_in) %>% in_to_mm(), 
            melt_days = n())

################################################################################
################################# Export data ##################################
################################################################################

saveRDS(object = snotel_snowfall_fraction,
        file = "../data/snotel/processed/snotel_monthly_snowfall_fraction_t50_optimized.RDS")

saveRDS(object = snotel_melt,
        file = "../data/snotel/processed/snotel_monthly_melt.RDS")



