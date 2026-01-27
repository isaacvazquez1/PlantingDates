#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
############    Purpose: Useful Functions for 
############          Sunflower Analyses
############            By: Eliza Clark
############       Last modified: 09/09/2024
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



# Optimum Timing from Gaussian model ----------------------------
# From Edwards & Crone, 2021
# Function to calculate optimum timing, etc. from model coefficients
pheno_calc=function(coefs #vector or matrix of coefficients, 
                    # ordered from b0 to b2 (if matrix, columns ordered that way)
){ 
  #convert vector to matrix if needed 
  if(!is.matrix(coefs)){
    coefs=matrix(coefs, nrow=1)
  } 
  res=data.frame(mu=numeric(nrow(coefs)), 
                 sd=numeric(nrow(coefs)), 
                 fp=numeric(nrow(coefs)), 
                 fst=numeric(nrow(coefs)), 
                 N=numeric(nrow(coefs)), 
                 h=numeric(nrow(coefs)), 
                 badfit = (coefs[,2] < 0) | (coefs[,3]>0)
  ) 
  rownames(res)=rownames(coefs) 
  res$mu=coefs[,2]/(-2*coefs[,3]) 
  res$sd=sqrt(1/(-2*coefs[,3])) 
  res$fp=qnorm(.9, mean=res$mu, sd=res$sd)-qnorm(.1, mean=res$mu, sd=res$sd) 
  res$fst=qnorm(.1, mean=res$mu, sd=res$sd) 
  res$N=sqrt(pi)/sqrt(-coefs[,3])*exp(coefs[,1]+coefs[,2]^2/(-4*coefs[,3])) 
  res$h=exp(coefs[,1] + coefs[,2]*res$mu + coefs[,3]*res$mu^2) 
  res$ofs=-sqrt(coefs[,2]^2 - 4 * coefs[,1] * coefs[,3])/coefs[,3] 
  res$ofs[is.na(res$ofs)]=0 
  return(res)
}


# Unscale Gaussian Coefficients ----------------------------
# From Edwards & Crone, 2021
# Function for reversing scaling of coefficient estimates 
coefs_unscaled = function(fixefs, mean, sd){ 
  #fixefs: the 3 coefficients of the regression fitted to scaled predictor 
  # in either matrix/data frame form or vector form #mean: mean of original unscaled predictor
  #sd: sd of original unscaled predictor 
  ## the equations below can be obtained by doing the algebra 
  ## based on the transformation of scaling 
  
  #if a vector, turn into a 1-row matrix
  if(is.null(nrow(fixefs))){
    fixefs=matrix(fixefs, nrow=1)
  }
  fixefs[,1] <- fixefs[,1] - fixefs[,2] * mean / sd + fixefs[,3] * mean^2 / sd^2
  fixefs[,2] <- fixefs[,2] / sd - 2 * fixefs[,3] * mean / sd^2
  fixefs[,3] <- fixefs[,3] / sd^2
  colnames(fixefs)=c("beta0", "beta1", "beta2") 
  return(fixefs)
}

# Unscale Linear Coefficients ----------------------------
# Function to unscale coefficients of linear model
lin_coef_unscale <- function(fixefs, mean, sd) {
  fixefs[1] <- fixefs[1] - fixefs[2] * mean / sd
  fixefs[2] <- fixefs[2] / sd
  names(fixefs)=c("lm_beta0", "lm_beta1") 
  return(fixefs)
}


# Unscale Multiple Linear Coefficients ----------------------------
# Function to unscale coefficients of linear model, when multiple years are 
# run as one model. This function is used in gaussfit.
lin_coef_unscale_mat <- function(fixefs, mean, sd) {
  fixefs[,1] <- fixefs[,1] - fixefs[,2] * mean / sd
  fixefs[,2] <- fixefs[,2] / sd
  colnames(fixefs)=c("lm_beta0", "lm_beta1") 
  return(fixefs)
}

# Lists of sites ----------------------------------------------------
# This function takes the full dataset and makes a list of sites 
# that can be used for downstream analyses.
# Currently, there are four accepted types of lists of sites that 
# can be generated and can be specified in the type argument.
# Accepted types: 'common_exact', 'common_city', 'all_exact', 'all_city'
# The cutoff for common sites can be adjusted manually, but default is 
# 5, thus locations with 5 or more site-yrs will be included in the common type lists

# IMPORTANT NOTE: All lists generated will only have the sites that 
# contain flower_50_doy data and do not 
# contain irrigated, NO, late, recrop, or short in the Location name !

# data = deriv_data

list_sites <- function(data, type, cutoff = 5, exclude_double = FALSE){

  # types: common_exact, all_exact, all_county
  
  if(type == 'common_exact'){
    common_sites_exactloc <- data %>% 
      filter(!is.na(flower_50_doy)) %>% 
      filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
      distinct(lat, lon, named_location, Year) %>%
      count(lat, lon, named_location) %>% 
      arrange(desc(n)) %>% 
      filter(n >= cutoff) %>%
      left_join(., data %>% select(Location, lat, lon, named_location, Year) %>%   
                  filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
                  distinct())
    
    exactloc_list <- common_sites_exactloc %>% group_split(named_location)
    names(exactloc_list) <- sapply(exactloc_list, \(x) x$named_location[1])
    return(exactloc_list)
  } 
  
  # if(type == 'common_city'){
  #   city_excl <- data %>% 
  #     filter(!is.na(flower_50_doy)) %>% 
  #     filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
  #     select(Year, Location, named_location, geolocation_precision, lat, lon) %>%
  #     distinct() %>% 
  #     mutate(named_loc_trunc = gsub('[[:punct:]].*', "", named_location)) %>%
  #     filter(geolocation_precision != 'county') %>%
  #     count(named_loc_trunc) %>% 
  #     arrange(desc(n)) %>% 
  #     filter(n < cutoff) 
  #   
  #   common_sites_city <- data %>% 
  #     filter(!is.na(flower_50_doy)) %>% 
  #     filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
  #     select(Year, Location, named_location, geolocation_precision, lat, lon) %>%
  #     distinct() %>%  
  #     mutate(named_loc_trunc = gsub('[[:punct:]].*', "", named_location)) %>%
  #     filter(geolocation_precision != 'county') %>%
  #     filter(!named_loc_trunc %in% city_excl$named_loc_trunc)
  #   
  #   city_list <- common_sites_city %>% group_split(named_loc_trunc)
  #   names(city_list) <- sapply(city_list, \(x) x$named_loc_trunc[1])
  #   return(city_list)
  # } 
  
  if(type == 'all_exact'){
    all_dryoil_sites <- data %>% 
      filter(!is.na(flower_50_doy)) %>% 
      filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
      distinct(Year, Location, named_location, lat, lon)
    
    all_exactloc_list <- all_dryoil_sites %>% group_split(named_location)
    names(all_exactloc_list) <- sapply(all_exactloc_list, \(x) x$named_location[1])
    return(all_exactloc_list)
  }
  
  # if(type == 'all_city'){
  #   # city group - all dryland oil sites with flowering doy data, grouped by city grouping
  #   all_sites_city <- data %>% 
  #     filter(!is.na(flower_50_doy)) %>% 
  #     filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
  #     distinct(Year, Location, named_location, garden_city, lat, lon) 
  #   
  #   all_city_list <- all_sites_city %>% group_split(garden_city)
  #   names(all_city_list) <- sapply(all_city_list, \(x) x$garden_city[1])
  #   return(all_city_list)
  # }
  
  if(type == 'all_county'){
    # city group - all dryland oil sites with flowering doy data, grouped by city grouping
    all_sites_county <- data %>% 
      filter(!is.na(flower_50_doy)) %>% 
      filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
      mutate(county_state = paste(garden_county, State, sep = '_')) %>%
      distinct(Year, Location, named_location, garden_city, county_state, lat, lon) 

      if(exclude_double){
      # exclude trials that are in the same county in the same year, but have 
        # planting dates that are different by more than 2 days
      all_sites_county <- data %>% 
        filter(!is.na(flower_50_doy)) %>% 
        filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
        # get the sites that have more than one trial in the same county in the same year
        inner_join(data %>% 
                     filter(!is.na(flower_50_doy)) %>% 
                     filter(!grepl('irrigated|NO|late|recrop|short|additional', Location)) %>%
                     select(garden_county, Year, State, Location) %>%
                     distinct() %>%
                     group_by(garden_county, Year, State) %>%
                     summarise(
                       n = n()
                     ) %>% filter(n > 1), 
                   by = c("garden_county", "Year", "State")) %>%
        # add how many hybrids are in each trial
        add_count(garden_county, Location, Year, name = "n_hybrids") %>%
        distinct(garden_county, Location, Year, planting_doy, planting_date, n_hybrids) %>%
        group_by(garden_county, Year) %>%
        # filter to keep only trials in which planting dates are different by more than 2 days
        filter(max(planting_doy) - min(planting_doy) > 2) %>%
        arrange(n_hybrids) %>%
        # keep the trial with fewer hybrids (to exclude from the list of sites)
        slice_head(n = -1) %>%
        # join back to list of sites, to remove the smaller duplicate trial from the list
        anti_join(all_sites_county, ., by = c('Location', 'Year'))
    }
    
    all_sites_county <- all_sites_county %>% group_split(county_state)
    names(all_sites_county) <- sapply(all_sites_county, \(x) x$county_state[1])
    return(all_sites_county)
  }
  
  stop('Check that argument type is accepted')
  
}
# #Example usage
# sites_out = list_sites(deriv_data, type = 'common_exact')
# sites_out[[3]]
