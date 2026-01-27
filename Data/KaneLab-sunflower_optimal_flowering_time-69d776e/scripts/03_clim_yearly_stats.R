#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
############    Purpose: Create climate variables 
############     for each location in each year
############- 
############            By: Eliza Clark
############       Last modified: 09/23/2024
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


# Packages ----------------------------------------------------------------

library(tidyverse) 
library(ggpubr)
library(GGally)
library(factoextra)


deriv_data_filepath <- "derived_data"
path2climatedata <- "data/"

# Data --------------------------------------------------------------------

## Climate Data ----------------

# # Daymet daily data
# 
#Pogoda server path
if(dir.exists("/home/elizaic/sunflower_geospatial/Sunflower_Geospatial/data/")){
  path2climatedata <- "/home/elizaic/sunflower_geospatial/Sunflower_Geospatial/data/"
}
# 
# #Laptop local path
# if(dir.exists("C:/Users/elcl6271/Local Documents")){
#   path2climatedata <- "C:/Users/elcl6271/Local Documents/CU SUNFLOWER/"
# }

clim_data <- read.csv(paste0(path2climatedata, "daymet_timeseries_cleaned.csv"))
# fresh copy of the climate data, because it takes forever to read in
clim_data_fresh <- clim_data


# uncomment if you want to run the SPEI drought index
# will need to download the SPEI dataset from online
# # SPEI drought index monthly data
# spei_data <- read.csv(paste0(path2climatedata, "spei_6mo.csv"))
# spei_data_fresh <- spei_data
# # spei_data <- spei_data_fresh
# 
# 
# # clean up dataset
# spei_data <- spei_data %>% group_by(location, latitude, longitude, Year, Month) %>%
#   summarise(
#     n = n(),
#     SPEI = mean(SPEI, na.rm = T),
#     sd_spei = sd(SPEI, na.rm = T)
#   ) %>% 
#   dplyr::select(-sd_spei, -n) %>%
#   ungroup()


## Phenology Data
deriv_data <- read.csv("data/sunflower_data_simple_v1.csv")

### fix up columns
deriv_data$Irrigated <- as.factor(deriv_data$Irrigated)
deriv_data$Oil_Confection <- as.factor(deriv_data$Oil_Confection)
deriv_data$Location <- as.factor(deriv_data$Location)
deriv_data$Unif_Name <- as.factor(deriv_data$Unif_Name)
deriv_data$flower_50_doy <- as.integer(deriv_data$flower_50_doy)
deriv_data$Year <- as.factor(deriv_data$Year)
deriv_data <- deriv_data %>% 
  mutate(county_state = paste(garden_county, State, sep = '_'))

# ## Geocoordinates Data - Location Data
# geocoordinates_mapped <- read.csv('lookup_tables/geocoordinates_mapped.csv')

## List of locations that are in the downstream analyses
R2_indata <- read.csv(paste0(deriv_data_filepath,"/R2_indata.csv"))


# Check that average summer temps for all years are similar to summer temps
# for 1980-2000, which is similar to the 'historical' temperatures that 
# the climate models compare to.
# For both max and min, the temperatures are similar for both sets of years.
# This means that it's not a problem to add the future anomaly to all 
# years of data, even though some of it is not 'historical'.

# Look at historical average summer temperatures
clim_data %>% 
  filter(year < 2001 & yday > 152 & yday < 243) %>%
  distinct() %>%
  summarise(tamxdegc_Ave = mean(tmaxdegc),
            tmindegc_ave = mean(tmindegc))

# All years summer average temperature
clim_data %>% 
  filter(yday > 152 & yday < 243) %>%
  distinct() %>%
  summarise(tamxdegc_Ave = mean(tmaxdegc),
            tmindegc_ave = mean(tmindegc))


# Prepare data ----------------------------------------------------------------



clim_yearly <- clim_data %>% filter(!is.na(year)) %>% 
  group_by(location, year, latitude, longitude) %>%
  summarize(
    n = n()
  ) %>% 
  mutate(year = as.factor(year))

#average planting, harvest, flower days by location
loc_averages_doys <- deriv_data %>% #distinct(Year, planting_doy, harvest_doy, flower_50_doy, Location) %>%
  rename(year = Year, location = Location) %>%
  mutate(location = as.character(location)) %>%
  group_by(location) %>%
  summarize(
    planting_doy_ave = mean(planting_doy, na.rm = T),
    harvest_doy_ave = mean(harvest_doy, na.rm = T),
    flower_50_doy_ave = mean(flower_50_doy, na.rm = T)
  ) 
  
loc_yr_ave_doys <- deriv_data %>% #distinct(Year, planting_doy, harvest_doy, flower_50_doy, Location) %>%
  rename(year = Year, location = Location) %>%
  mutate(location = as.character(location)) %>%
  group_by(location, year) %>%
  summarize(
    planting_doy = mean(planting_doy, na.rm = T),
    harvest_doy = mean(harvest_doy, na.rm = T),
    flower_50_doy = mean(flower_50_doy, na.rm = T)
  )


# join yearly climate names to site phenology averages 
# and join to actual phenology values for years that have them
# and create new columns that take the actual value, if available, if not take average site value
clim_yearly <- clim_yearly %>% 
  left_join(., loc_averages_doys) %>%
  left_join(., loc_yr_ave_doys) %>%
    mutate(planting_doy_combined = ifelse(is.na(planting_doy), planting_doy_ave, planting_doy),
           harvest_doy_combined = ifelse(is.na(harvest_doy), harvest_doy_ave, harvest_doy),
           flower_doy_combined = ifelse(is.na(flower_50_doy), flower_50_doy_ave, flower_50_doy)
           ) %>%
  ungroup()            


#location, named_location, year list of values
named_loc_lookup <- deriv_data %>% distinct(Location, named_location, Year) %>%
  rename(year = Year, location = Location) %>%
  arrange(location, year)

# add named_location to clim_yearly
clim_yearly <- clim_yearly %>% left_join(., named_loc_lookup) 

# add county info
geocoord_lookup <- deriv_data %>% 
  distinct(Location, lat, lon, county_state, State) %>% arrange(Location)


# geocoordinates_mapped <- geocoordinates_mapped %>%
#   mutate(county_state = paste(garden_county, State, sep = "_"))
# geocoord_lookup <- geocoordinates_mapped %>%
#   distinct(Location, lat, lon, county_state, State) %>% arrange(Location)

clim_yearly <- clim_yearly %>% 
  left_join(., geocoord_lookup, relationship = 'many-to-one', 
            by = c("location" = 'Location', 'latitude' = 'lat', 'longitude' = 'lon'))


# and make year into an integer
# and replace NaN with NA
# !! not this (and remove sites that don't have any phenology data)
# remove sites that are not in the downstream analyses
# remove unncessary columns - keep phenology_combined columns
## NOTE: Do not change the order of columns in dplyr::select!!
## The columns are referenced by position, not name in the functions!
clim_yearly <- clim_yearly %>% 
  mutate(year = as.character(year)) %>% 
  mutate_all(function(x) ifelse(is.nan(x), NA, x)) %>% 
  filter(location %in% R2_indata$Location) %>%
  # filter(!is.na(planting_doy_combined) & !is.na(flower_doy_combined)) %>%
  dplyr::select(location, year, latitude, longitude, 
                n, planting_doy_combined, harvest_doy_combined, flower_doy_combined,
                county_state, State) 



# check there is only one planting doy per location per year
clim_yearly %>% add_count(location, latitude, longitude, year) %>%
  filter(nn > 1)




# find average window around planting
# Doy 125 to 250 will be development time for basically all locations
# including planting and flowering time
clim_yearly %>% 
  summarise(
    ave_planting = mean(planting_doy_combined, na.rm = T), 
    sd_planting = sd(planting_doy_combined, na.rm = T),
    ave_flower = mean(flower_doy_combined, na.rm = T), 
    sd_flower = sd(flower_doy_combined, na.rm = T))
clim_yearly %>% ggplot() +
  geom_histogram(aes(x = flower_doy_combined, fill = 'flowering'), alpha = 0.5) +
  geom_histogram(aes(x = planting_doy_combined, fill = 'planting'), alpha = 0.5) +
  geom_histogram(aes(x = harvest_doy_combined, fill = 'harvest'), alpha = 0.5)



# make into a list for faster calculations
clim_yearly_list <- clim_yearly %>% group_split(location, latitude, longitude)

# clim_yearly_list[[2]]

# FUNCTIONS --------------------------------------------------------------

{
## Growing Degree Days ----------------------------------------------------

# site = clim_yearly_list[[48]]
# climate.dat <- clim_data

gdd_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # calculate growing degree days
  tbase = 6.7
  clim_site <- clim_site %>% mutate(GDD = pmax(0, (tmaxdegc + tmindegc)/2 - tbase))
  
  # generate cumulative degree days for the 70 days after planting
  gdd_70 <- NULL
  
  # start = planting doy; end = planting doy + 70
  for (cur_yr in site$year) { 
    start_cur = site[site$year == cur_yr, 6] %>% as.numeric()
    end_cur = start_cur + 70
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    if (nrow(dat_cur) == 0) {
      gdd_cur <- NA } 
    else {
      gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    }
    # gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    gdd_70 <- c(gdd_70, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(CUM_GDD_70 = gdd_70)

  site
}

## Growing Degree Days - upper threshold -------------------------------------

gdd_modave_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # calculate growing degree days - modified average formula
  tbase = 6.7
  tmax = 26
  clim_site <- clim_site %>% mutate(GDD_modave = pmax(0, (
    ifelse(tmaxdegc > tmax, tmax, tmaxdegc) + ifelse(tmindegc < tbase, tbase, tmindegc)
  ) / 2 - tbase))
  
  # generate cumulative degree days for the 70 days after planting
  clim_out <- NULL
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 6] %>% as.numeric()
    end_cur = start_cur + 70
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    if (nrow(dat_cur) == 0) {
      gdd_cur <- NA } 
    else {
      gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    }
    # gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(GDD_modave = clim_out)

  site
}

## Future GDD - upper threshold ---------------------------------------------
# Can choose between an average of RCP 4.5 and 8.5 scenarios, or each scenario by itself.
# Could update code to have all of these run. For now, uncomment the appropriate line below
# and leave the other two commented out.

gdd_mod_future_calc <- function(site, climate.dat){
  # projected increases in min and max temperatures (average of RCP 4.5 and RCP 8.5)
  # temp_increases <- matrix(c(2.918611, 3.3522), nrow = 2, ncol = 1, dimnames = list(c("min", "max"), c("increase_C")))
  # projected increases in min and max temperatures (RCP 8.5)
  # temp_increases <- matrix(c(3.4389, 3.9339), nrow = 2, ncol = 1, dimnames = list(c("min", "max"), c("increase_C")))
  # projected increases in min and max temperatures (RCP 4.5)
  temp_increases <- matrix(c(2.4411, 2.8861), nrow = 2, ncol = 1, dimnames = list(c("min", "max"), c("increase_C")))
  
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # calculate growing degree days - modified average formula
  tbase = 6.7
  tmax = 26
  clim_site <- clim_site %>% mutate(GDD_modave = pmax(0, (
    ifelse(tmaxdegc + temp_increases[2] > tmax, tmax, tmaxdegc + temp_increases[2]) + 
      ifelse(tmindegc + temp_increases[1] < tbase, tbase, tmindegc + temp_increases[1])
  ) / 2 - tbase))
  
  # generate cumulative degree days for the 70 days after planting
  clim_out <- NULL
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 6] %>% as.numeric()
    end_cur = start_cur + 70
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    if (nrow(dat_cur) == 0) {
      gdd_cur <- NA } 
    else {
      gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    }
    # gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(GDD_future = clim_out)
  
  site
}

## Fixed period GDD - upper threshold ---------------------------------------------

gdd_mod_fixed_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # calculate growing degree days - modified average formula
  tbase = 6.7
  tmax = 26
  clim_site <- clim_site %>% mutate(GDD_modave = pmax(0, (
    ifelse(tmaxdegc > tmax, tmax, tmaxdegc) + ifelse(tmindegc < tbase, tbase, tmindegc)
  ) / 2 - tbase))
  
  # generate cumulative degree days for the 70 days after planting
  clim_out <- NULL
  
  start = 125
  end = 250
  for (cur_yr in site$year) {
    dat_cur = clim_site %>% filter(between(yday, start, end) &
                                     year == cur_yr)
    
    gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(GDD_fixed = clim_out)
  
  site
}

## Vegetative GDD-mod ------------------------------------------------------------

gdd_veg_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # calculate growing degree days - modified average formula
  tbase = 6.7
  tmax = 26
  clim_site <- clim_site %>% mutate(GDD_modave = pmax(0, (
    ifelse(tmaxdegc > tmax, tmax, tmaxdegc) + ifelse(tmindegc < tbase, tbase, tmindegc)
  ) / 2 - tbase))
  
  # generate cumulative degree days for the 70 days after planting
  clim_out <- NULL
  
  # start = mean(site$planting_doy, na.rm = T)
  # end = mean(site$flower_50_doy, na.rm = T) - 21

  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 6] %>% as.numeric()
    end_cur = site[site$year == cur_yr, 8]  %>% as.numeric() - 21
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(GDD_veg = clim_out)
  
  site
}

## Reproductive GDD-mod ------------------------------------------------------------

gdd_repro_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # calculate growing degree days - modified average formula
  tbase = 6.7
  tmax = 26
  clim_site <- clim_site %>% mutate(GDD_modave = pmax(0, (
    ifelse(tmaxdegc > tmax, tmax, tmaxdegc) + ifelse(tmindegc < tbase, tbase, tmindegc)
  ) / 2 - tbase))
  
  # generate cumulative degree days for the 70 days after planting
  clim_out <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) - 21
  # end = mean(site$flower_50_doy, na.rm = T) + 30
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8]  %>% as.numeric() - 21
    end_cur = site[site$year == cur_yr, 8]  %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(GDD_repro = clim_out)
  
  site
}

## Long-season GDD-mod ------------------------------------------------------------

gdd_long_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # calculate growing degree days - modified average formula
  tbase = 6.7
  tmax = 26
  clim_site <- clim_site %>% mutate(GDD_modave = pmax(0, (
    ifelse(tmaxdegc > tmax, tmax, tmaxdegc) + ifelse(tmindegc < tbase, tbase, tmindegc)
  ) / 2 - tbase))
  
  # generate cumulative degree days for the 70 days after planting
  clim_out <- NULL
  
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 6] %>% as.numeric()
    end_cur = site[site$year == cur_yr, 8]  %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    gdd_cur <- sum(dat_cur$GDD, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(GDD_long = clim_out)
  
  site
}

## Heat Stress Degree Days ----------------------------------------------------

gddheat_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # calculate heat stress degree days
  tbase = 26
  clim_site <- clim_site %>% mutate(GDD_heat = pmax(0, (tmaxdegc + tmindegc)/2 - tbase))
  
  # generate heat degree days 21 days before and 30 days after average flowering
  gdd_heat <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) - 21
  # end = mean(site$flower_50_doy, na.rm = T) + 30
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8] %>% as.numeric() - 21
    end_cur = site[site$year == cur_yr, 8] %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    gdd_cur <- sum(dat_cur$GDD_heat, na.rm = T)
    gdd_heat <- c(gdd_heat, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(GDD_heat = gdd_heat)

  site
}

## Total Precipitation ------------------------------------------------------

totprecip_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude) %>%
    distinct(date, .keep_all = T)
  
  # generate total precip 21 days before and 30 days after average flowering
  clim_out <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) - 21
  # end = mean(site$flower_50_doy, na.rm = T) + 30
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8] %>% as.numeric() - 21
    end_cur = site[site$year == cur_yr, 8] %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    gdd_cur <- sum(dat_cur$prcpmmday, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(tot_precip = clim_out)

  site
}

## Precipitation Frequency --------------------------------------------------

precipfreq_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude)%>%
    distinct(date, .keep_all = T)
  
  # generate number of days of precip 
  # betweem 21 days before and 30 days after average flowering
  
  clim_out <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) - 21
  # end = mean(site$flower_50_doy, na.rm = T) + 30
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8] %>% as.numeric() - 21
    end_cur = site[site$year == cur_yr, 8] %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    gdd_cur <- nrow(dat_cur %>% filter(prcpmmday > 0))
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(precip_freq = clim_out)

  site
}

## Precipitation Days without --------------------------------------------------

precipdayswithout_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude)%>%
    distinct(date, .keep_all = T)
  
  # generate longest number of consecutive days without precipitation 
  # between 21 days before and 30 days after average flowering
  clim_out <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) - 21
  # end = mean(site$flower_50_doy, na.rm = T) + 30
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8] %>% as.numeric() - 21
    end_cur = site[site$year == cur_yr, 8] %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    x <- rle(dat_cur$prcpmmday == 0)
    gdd_cur <- max(x$lengths[x$values == TRUE])
    
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(precip_dayswithout = clim_out)

  site
}

## First Frost Day ----------------------------------------------------------

frostday_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude)%>%
    distinct(date, .keep_all = T)
  
  # generate first day with temp below -4.44C after 21 days before average flowering
  clim_out <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) - 21

  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8] %>% as.numeric() - 21
    # end_cur = site[site$year == cur_yr, 8] %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, 365) &
                                     year == cur_yr)
    
    frost_days <- dat_cur %>% filter(tmindegc < -3.9) %>% dplyr::select(yday)
    gdd_cur <- ifelse(nrow(frost_days) > 0, min(frost_days), 366)

    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(frost_day = clim_out)
  site
}

## SPEI Drought Index ----------------------------------------------------------

spei_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude)
  
  # get spei index for the average month of flowering 
  clim_out <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) %>% as.Date(., origin = "2001-01-01") %>% month()
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8] %>% as.numeric() %>% 
      as.Date(., origin = "2001-01-01") %>% month()
    dat_cur = clim_site %>% filter(Month == start_cur &
                                     Year == cur_yr)
    
    gdd_cur <- dat_cur %>% dplyr::select(SPEI) %>% as.numeric()
    # gdd_cur <- pmin(gdd_cur, 366)
    
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(SPEI = clim_out)
  site
}

## Solar Radiation ------------------------------------------------------------

solrad_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude)%>%
    distinct(date, .keep_all = T)
  
  # generate average solar radiation  
  # between 21 days before and 30 days after average flowering
  clim_out <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) - 21
  # end = mean(site$flower_50_doy, na.rm = T) + 30
  
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8] %>% as.numeric() - 21
    end_cur = site[site$year == cur_yr, 8] %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    gdd_cur <- mean(dat_cur$sradWm2, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(sol_rad = clim_out)
  
  site
}

## Daylength ----------------------------------------------------------------

daylen_calc <- function(site, climate.dat){
  # climate data for site
  clim_site <- climate.dat %>% filter(location %in% site$location &
                                        latitude %in% site$latitude &
                                        longitude %in% site$longitude)%>%
    distinct(date, .keep_all = T)
  
  # generate mean daylength 
  # between 21 days before and 30 days after average flowering
  clim_out <- NULL
  
  # start = mean(site$flower_50_doy, na.rm = T) - 21
  # end = mean(site$flower_50_doy, na.rm = T) + 30
  # 
  for (cur_yr in site$year) {
    start_cur = site[site$year == cur_yr, 8] %>% as.numeric() - 21
    end_cur = site[site$year == cur_yr, 8] %>% as.numeric() + 30
    dat_cur = clim_site %>% filter(between(yday, start_cur, end_cur) &
                                     year == cur_yr)
    
    gdd_cur <- mean(dat_cur$dayls, na.rm = T)
    clim_out <- c(clim_out, gdd_cur)
  }
  
  # output
  site <- site %>% 
    mutate(day_length = clim_out)
  
  site
}

}

# Run Functions ---------------------------------------------------------
# Each functions will add that climate variables to the list
# All of these lines work, but the take a while to run, so only uncomment
# the ones that you want to run.
clim_yearly_list <- lapply(clim_yearly_list, gdd_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, gddheat_calc, climate.dat = clim_data)

clim_yearly_list <- lapply(clim_yearly_list, gdd_modave_calc, climate.dat = clim_data)

clim_yearly_list <- lapply(clim_yearly_list, gdd_mod_future_calc, climate.dat = clim_data)

clim_yearly_list <- lapply(clim_yearly_list, gdd_mod_fixed_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, gdd_veg_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, gdd_repro_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, gdd_long_calc, climate.dat = clim_data)

clim_yearly_list <- lapply(clim_yearly_list, totprecip_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, precipfreq_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, precipdayswithout_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, spei_calc, climate.dat = spei_data)

clim_yearly_list <- lapply(clim_yearly_list, frostday_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, solrad_calc, climate.dat = clim_data)

# clim_yearly_list <- lapply(clim_yearly_list, daylen_calc, climate.dat = clim_data)


# Collapse the list into a dataframe
clim_yearly_out <- do.call(rbind.data.frame, clim_yearly_list)

# # add state column
# clim_yearly_out <- clim_yearly_out %>% left_join(., deriv_data %>% dplyr::select(Location, State) %>%
#                                                    rename(location = Location) %>% 
#                                                    distinct() %>% filter(!is.na(State))) %>%
#   relocate(State, .after = location) %>%
#   mutate(EW = ifelse(longitude > median(clim_yearly_out$longitude), "E", "W"), .after = State,
#          State_EW = paste(State, EW, sep = "_"))

# Write out yearly climate data
write.csv(clim_yearly_out, paste0(deriv_data_filepath,"/yearly_climate_data.csv"))





# Correlations between climate variables  ------------------------------------

#extract column number of the first and last columns that contain climate data
col1 <- which(colnames(clim_yearly_out) == "CUM_GDD_70")
col2 <-which(colnames(clim_yearly_out) == "frost_day")
# col2 <-which(colnames(clim_yearly_out) == "day_length")


# simple correlation plots
clim_yearly_out %>% dplyr::select(col1:col2) %>% pairs()

clim_yearly_out %>% add_count(location) %>%
  filter(nn > 5) %>%
  GGally::ggpairs(mapping = aes(color = as.factor(location)), columns = col1:col2)



# Correlation plots by location

# Example plot
# by site 
ggplot(clim_yearly_out, aes(x = planting_doy_combined, y = GDD_modave, color = location)) +
  geom_point() +
  geom_smooth(se = F, method = 'lm') +
  guides(color = 'none') +
  theme_classic()


# Custom function for correlation plots with groups
my_cor_plot <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(alpha = 0.2) + 
    geom_smooth(se = F, method = 'lm', ...)
  p
}

# custom function to plot name of variable in the diagonal
my_diag <- function(data, mapping, ...){
  ggally_text(rlang::as_label(mapping$x), col="black", size=4) +
    theme_void()}


# Correlation plots with lines for each location
clim_yearly_out %>% add_count(location) %>%
  filter(nn > 5) %>%
  ggpairs(columns = col1:col2,
         mapping = aes(color = as.factor(location)),
         lower = list(continuous = my_cor_plot),
         diag = list("continuous" = my_diag))
#export to pdf 28x20in

# Correlations by State
clim_yearly_out %>% add_count(State) %>%
  filter(nn > 3) %>%
  ggpairs(columns = col1:col2,
          mapping = aes(color = as.factor(State)),
          lower = list(continuous = my_cor_plot),
          diag = list("continuous" = my_diag))


# PCA of climate variables -------------------------------------------------
#!! Most of these need updated subsetting, if not all climate variables
# were run above

# check for missing values
colSums(is.na(clim_yearly_out))

# PCA with fewer variables
clim_scaled_light <- clim_yearly_out %>% 
  select(CUM_GDD_70:frost_day, -GDD_future) %>% scale() 
  # select(CUM_GDD_70:day_length, -GDD_future, -GDD_repro, - GDD_veg, -GDD_long) %>% scale() # if all climate variables were made above
data_pca_light <- clim_scaled_light %>%
  princomp()
summary(data_pca_light)
data_pca_light$loadings
# PCA biplot with state
fviz_pca(data_pca_light, col.var = 'black', col.ind = 'grey80', label = 'var', repel = T,
         # habillage = clim_yearly_out$State, 
         addEllipses = F, ellipse.level = 0.95
)

#!! This requires a different dataset from a later file
# PCA with variables in analysis
pca_mod <- optim_clim_sc %>% select(GDD_mod_sitemean_sc:frost_anomaly_sc) 
data_pca_mod <- pca_mod %>% princomp()
summary(data_pca_mod)
data_pca_mod$loadings
fviz_pca(data_pca_mod, col.var = 'black', col.ind = 'grey80', label = 'var', repel = T,
         # habillage = optim_clim_sc$State,
         addEllipses = F, ellipse.level = 0.95
)


# PCAs with all variables
clim_scaled <- clim_yearly_out %>% select(CUM_GDD_70:frost_day) %>% scale() 
data_pca <- clim_scaled %>%
  princomp()
summary(data_pca)
data_pca$loadings


fviz_eig(data_pca, addlabels = T)

# PCA biplot with state
fviz_pca(data_pca, col.var = 'black', col.ind = 'grey80', label = 'var', repel = T,
         # habillage = clim_yearly_out$State, 
         addEllipses = F, ellipse.level = 0.95
         )

# biplot with location (messy)
fviz_pca(data_pca, col.var = 'black', label = 'var',
         habillage = clim_yearly_out$location)

# # biplot with E/W
# fviz_pca(data_pca, col.var = 'black', label = 'var',
#          habillage = clim_yearly_out$EW,
#          addEllipses = TRUE, ellipse.level = 0.95
# )
# 
# # biplot with E/W & State
# fviz_pca(data_pca, col.var = 'black', label = 'var',
#          habillage = clim_yearly_out$State_EW,
#          addEllipses = TRUE, ellipse.level = 0.95
# )

# biplot with Year
fviz_pca(data_pca, col.var = 'black', label = 'var',
         habillage = clim_yearly_out$year,
)




# PCA of Temp-related variables
clim_scaled_temp <- clim_yearly_out %>% select(CUM_GDD_70:GDD_repro, frost_day:day_length) %>% scale() 
data_pca_temp <- clim_scaled_temp %>%
  princomp()
summary(data_pca_temp)

# PCA biplot with state
fviz_pca(data_pca_temp, col.var = 'black', label = 'var',
         habillage = clim_yearly_out$State, 
         addEllipses = TRUE, ellipse.level = 0.95
)


# PCA of Precip-related variables
clim_scaled_precip <- clim_yearly_out %>% select(tot_precip:SPEI) %>% scale() 
data_pca_prec <- clim_scaled_precip %>%
  princomp()
summary(data_pca_prec)

# PCA biplot with state
fviz_pca(data_pca_prec, col.var = 'black', label = 'var',
         habillage = clim_yearly_out$State, 
         addEllipses = TRUE, ellipse.level = 0.95
)


# PCA of for each state
clim_scaled_ND <- clim_yearly_out %>% 
  filter(State == 'ND') %>%
  select(CUM_GDD_70:day_length) %>% scale() 
data_pca_ND <- clim_scaled_ND %>%
  princomp()
summary(data_pca_ND)

clim_scaled_SD <- clim_yearly_out %>% 
  filter(State == 'SD') %>%
  select(CUM_GDD_70:day_length) %>% scale() 
data_pca_SD <- clim_scaled_SD %>%
  princomp()
summary(data_pca_SD)

clim_scaled_KS <- clim_yearly_out %>% 
  filter(State == 'KS') %>%
  select(CUM_GDD_70:day_length) %>% scale() 
data_pca_KS <- clim_scaled_KS %>%
  princomp()
summary(data_pca_KS)

clim_scaled_TX <- clim_yearly_out %>% 
  filter(State == 'TX') %>%
  select(CUM_GDD_70:day_length) %>% scale() 
data_pca_TX <- clim_scaled_TX %>%
  princomp()
summary(data_pca_TX)

clim_scaled_NE <- clim_yearly_out %>% 
  filter(State == 'NE') %>%
  select(CUM_GDD_70:day_length) %>% scale() 
data_pca_NE <- clim_scaled_NE %>%
  princomp()
summary(data_pca_NE)

# PCA biplot with state
ggarrange(
fviz_pca(data_pca_ND, col.var = 'black', label = 'var', title = 'ND',
         habillage = clim_yearly_out %>% filter(State == 'ND') %>% .$location,
         addEllipses = TRUE, ellipse.level = 0.95 ) +
  theme(legend.position = 'none'),
fviz_pca(data_pca_SD, col.var = 'black', label = 'var', title = 'SD',
         habillage = clim_yearly_out %>% filter(State == 'SD') %>% .$location,
         addEllipses = TRUE, ellipse.level = 0.95 ) +
  theme(legend.position = 'none'),
fviz_pca(data_pca_KS, col.var = 'black', label = 'var', title = 'KS',
         habillage = clim_yearly_out %>% filter(State == 'KS') %>% .$location,
         addEllipses = TRUE, ellipse.level = 0.95 ) +
  theme(legend.position = 'none'),
fviz_pca(data_pca_NE, col.var = 'black', label = 'var', title = 'NE',
         habillage = clim_yearly_out %>% filter(State == 'NE') %>% .$location,
         addEllipses = TRUE, ellipse.level = 0.95 ) +
  theme(legend.position = 'none'),
fviz_pca(data_pca_TX, col.var = 'black', label = 'var', title = 'TX',
         habillage = clim_yearly_out %>% filter(State == 'TX') %>% .$location,
         addEllipses = TRUE, ellipse.level = 0.95) +
  theme(legend.position = 'none'),
fviz_pca(data_pca, col.var = 'black', label = 'var', title = 'ALL',
         habillage = clim_yearly_out$State, 
         addEllipses = TRUE, ellipse.level = 0.95)
)



