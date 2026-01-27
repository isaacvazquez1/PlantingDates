
# Summary Statistics for planting date and flowering time
# EIC March 2025



# Packages ----------------------------------------------------------------
library(tidyverse)

deriv_data_filepath <- "derived_data"

# Data --------------------------------------------------------------------


## R2 data subset ---------------------------------------------------------
# This is the largest subset of data that was included in analyses

R2_indata <- read.csv(paste0(deriv_data_filepath,"/R2_indata.csv"))


# how many named locations, hybrids, years
R2_indata %>% count(named_location) %>% nrow() # 65 unique locations
R2_indata %>% count(county_state) %>% nrow() # 37 location groups (counties)
R2_indata %>% group_by(county_state, Year) %>% n_groups() # 341 site years

R2_indata %>% distinct(county_state, Year) %>% group_by(county_state) %>% 
  summarise(n = n()) %>%
  summarise(mean = mean(n), median = median(n), sd = sd(n),
            min = min(n), max = max(n)) # 9 years +/- 10.5 years per county

R2_indata %>% count(county_state, Year) %>%
  summarise(mean = mean(n), median = median(n), sd = sd(n), min = min(n), max = max(n)) # 43 +/-25(sd) hybrids per site-year 
R2_indata %>% count(Unif_Name) %>% nrow() # 2158 hybrids
R2_indata %>% count(Year) %>% nrow() # 46 years



##
## 06a subset - climate ---------------------------------------------------
## set of data after excluding trials not used in climate analyses

optim_clim <- read.csv(paste0(deriv_data_filepath,"/optim_clim_sc_county.csv"))


# number of unique sites
optim_clim %>% group_by(named_location) %>% n_groups() 
optim_clim %>% group_by(county_state) %>% n_groups() # 26 counties
optim_clim %>% group_by(Year) %>% n_groups() # 44 years

### number of site-years
optim_clim %>% group_by(named_location, Year) %>% n_groups() # 255 site years






