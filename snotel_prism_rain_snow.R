# Script for estimating monthly SNOTEL snowfall fractions
# Uses optimized rain-snow temperature from Jennings et al. (2018)

# Load packages
library(raster)
library(sp)
library(tidyverse)
library(velox)
library(foreach)
library(doMC); registerDoMC(cores = 4)
library(reshape2)
library(humidity) # devtools::install_github("https://github.com/SnowHydrology/humidity)

################################################################################
################################# Import data ##################################
################################################################################

# Daily snotel data
snotel <- readRDS("../data/snotel/processed/snotel_daily_formatted.RDS")

# PRISM tair and tdew
# Only works on KSJ local machine with external drive -- data v big
prism_tair_path <- "/Volumes/files/climate_data/prism/daily/tair/"
prism_tdew_path <- "~/Downloads/tdew/"

# list all the met files
met_files <- data.frame(tair_files = list.files(pattern = "\\.bil$",
                                                path = prism_tair_path), 
                        stringsAsFactors = F) %>% 
  mutate(date = as.Date(str_sub(tair_files, -16, -9), format = "%Y%m%d")) %>% 
  filter(date <= as.Date("2018-12-31")) %>% 
  left_join(., 
            data.frame(tdew_files = list.files(pattern = "\\.bil$",
                                               path = prism_tdew_path), 
                       stringsAsFactors = F) %>% 
              mutate(date = as.Date(str_sub(tdew_files, -16, -9), format = "%Y%m%d")),
            by = "date") %>% 
  arrange(date) %>% 
  mutate(year = year(date))

# Make a vector of the years to be analyzed
years <- sort(unique(met_files$year)) 

# Import SNOTEL stations
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

# Convert snotel data to spatialpointsdataframe
coords <- dplyr::select(snotel_meta, lon, lat)
id <- dplyr::select(snotel_meta, site_id)
crs_snotel <- CRS("+init=epsg:4326")

# Convert obs to a SpatialPointsDataFrame
snotel_sp <- SpatialPointsDataFrame(coords = coords, 
                                    data = id, 
                                    proj4string = crs_snotel)

###############################################################################
##########################  Extract PRISM Data  ###############################
###############################################################################

# Run parallel loop to extract tair data by year

###########start timer###############
start.time = Sys.time()

# Parallel loop
tair_list <- 
  foreach(i = 1:length(years), .errorhandling = "pass", .verbose = T) %dopar% {
    
    #Filter to just 1 year
    tmp.files <- filter(met_files, year == years[i])
    
    #Put prism data into raster stack velox object
    prism.v <- velox(stack(paste0(prism_tair_path, tmp.files$tair_files)))
    
    #Extract the data using the locs of interest
    tair_extract <- prism.v$extract_points(sp = snotel_sp)
    
    #Output the extracted data to the ppt list
    tair_extract
  }

###########end timer###############
end.time = Sys.time()
end.time - start.time
###########end###############

# Run parallel loop to extract tair data by year

###########start timer###############
start.time = Sys.time()

# Parallel loop
tdew_list <- 
  foreach(i = 1:length(years), .errorhandling = "pass", .verbose = T) %dopar% {
    
    #Filter to just 1 year
    tmp.files <- filter(met_files, year == years[i])
    
    #Put prism data into raster stack velox object
    prism.v <- velox(stack(paste0(prism_tdew_path, tmp.files$tdew_files)))
    
    #Extract the data using the locs of interest
    tdew_extract <- prism.v$extract_points(sp = snotel_sp)
    
    #Output the extracted data to the ppt list
    tdew_extract
  }

###########end timer###############
end.time = Sys.time()
end.time - start.time
###########end###############


###############################################################################
##########################  Format and Join Data  #############################
###############################################################################

# Put the list into a vertical data frame format

#Make dummy data frame
met.df <- data.frame()

#Loop through each list element
for(i in 1:length(tair_list)){
  
  #Make data frame of tair values from list element i and add the lat-lon
  tmp.df <- as.data.frame(tair_list[[i]]) %>%
    mutate(site_id = snotel_meta$site_id)
  
  #melt the data frame
  tmp.df.melt <- melt(tmp.df, id.vars = "site_id",
                      variable.name = "column", 
                      value.name = "tair_c") %>% 
    dplyr::select(., -column)
  
  #Add the date
  tmp.date.seq <- rep(seq.Date(from = as.Date(paste0(years[i], "-01-01")),
                               to = as.Date(paste0(years[i], "-12-31")),
                               by = "1 day"), each = length(tmp.df$site_id)) #repeats each date n times
  tmp.df.melt$date <- tmp.date.seq
  
  #Make data frame of tdew values from list element i and add the lat-lon
  tmp.df <- as.data.frame(tdew_list[[i]]) %>%
    mutate(site_id = snotel_meta$site_id)
  
  #melt the data frame
  tmp.df.melt2 <- melt(tmp.df, id.vars = "site_id",
                      variable.name = "column", 
                      value.name = "tdew_c") %>% 
    dplyr::select(., -column)
  
  #Add the date
  tmp.df.melt2$date <- tmp.date.seq
  
  # Join the data
  tmp.df.melt <- left_join(tmp.df.melt, tmp.df.melt2,
                           by = c("site_id", "date"))
  
  #Bind to data frame
  met.df <- bind_rows(met.df, tmp.df.melt)
}


###############################################################################
####################  Add PRISM to SNOTEL and Analyze  ########################
###############################################################################

# Add date to SNOTEL
snotel <- snotel %>% 
  mutate(date = as.Date(paste(year, month, day, sep = "-")))

# Join data
snotel <- left_join(snotel, met.df,
                    by = c("site_id", "date"))

# Compute relative humidity
snotel <- snotel %>% 
  mutate(rh = relhum(tair_c, tdew_c), 
         rh = case_when(rh > 100 ~ 100,
                        TRUE ~ rh))

# Function for computing snowfall probability
# This is optimized binary logistic regression from Jennings et al. (2018)
p_snow <- function(TAIR, RH){
  (1/(1 + exp(-10.04 + 1.41 * TAIR + 0.09 * RH)))
}

# Compute snowfall probability
# >= 0.5 = snow
# And snowfall
snotel <- snotel %>% 
  mutate(snow_prob = p_snow(tavg_c, rh),
         snowfall_in_binlog = case_when(snow_prob >= 0.5 ~ ppt_in,
                                        TRUE ~ 0))

# Inch to mm function
in_to_mm = function(x){x * 25.4}

# Compute monthly snowfall fraction per site
snotel_snowfall_fraction <- snotel %>% 
  filter(!is.na(ppt_in) & !is.na(snowfall_in_binlog)) %>% 
  group_by(site_id, month, year) %>% 
  summarise(ppt_mm = sum(ppt_in) %>% in_to_mm(), 
            snowfall_mm = sum(snowfall_in_binlog) %>% in_to_mm(),
            snowfall_frac = snowfall_mm / ppt_mm, 
            n_obs = n())

################################################################################
################################# Export data ##################################
################################################################################

saveRDS(object = snotel,
        file = "../data/snotel/processed/snotel_daily_formatted_withRH.RDS")

saveRDS(object = snotel_snowfall_fraction,
        file = "../data/snotel/processed/snotel_monthly_snowfall_fraction_prism_binlog.RDS")


