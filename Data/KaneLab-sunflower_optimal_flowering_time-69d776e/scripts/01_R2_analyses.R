
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
############   Purpose: This file extracts R2 and regression  
############        parameters from the gaussian model
############          of optimal flowering time 
############       (flowering doy predicting yield)
############           for downstream analyses
############        and runs meta-regressions on R2.
############.
############.        
############               By: Eliza Clark
############             Created: 07/26/2024
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


# Packages ----------------------------------------------------------------

library(MASS)
library(tidyverse) 
library(ggpubr)

select <- dplyr::select

# Data ------------------------------------------------------------------
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


deriv_data_filepath <- "derived_data"
figure_filepath <- "figures"

# Functions -------------------------------------------------------------

# Custom functions
source("scripts/00_analysis_functions.R")



# Function to generate R2 data for both gaussian and linear models for 
# a list of site-years. Works on both exactlocation list and citylocation list.

# # Development
# site_list = list_sites(deriv_data, type = 'all_county', exclude_double = TRUE)
# i=7
# year_i = 1989

R2_data <- function(site_list){
  
  # make empty data frame of the output 
  R2_out <- data.frame(
    Year = NA,
    # Location = NA,
    named_location = NA,
    group = NA,
    State = NA,
    lat = NA,
    lon = NA,
    n_site_yr = NA,
    spread_flower_site_yr = NA,
    pseudoR2_gau = NA,
    badfit_gau = NA,
    beta0_gau = NA,
    beta1_gau = NA,
    beta2_gau = NA,
    AIC_gau = NA,
    R2_lm = NA,
    beta0_lm = NA,
    beta1_lm = NA,
    beta1_lm_sig = NA,
    AIC_lm = NA,
    preferred_mod = NA,
    mod_for_cens = NA,
    higher_R2 = NA
  )
  
  for (i in 1:length(site_list)) {
    
    # subset data into each of the common sites - 
    # this is done by both lat/long and location name.
    # Add a column to scale flowering doy by the mean/sd of the entire location
    
    site_dat_i <- deriv_data %>% 
      semi_join(., site_list[[i]], by = c('Year', 'Location', 'named_location')) %>%
      filter(
        # lat == site_list$lat[i] &
        #   lon == site_list$lon[i] &
        #   Location == site_list$Location[i] &
        # named_location %in% unique(site_list[[i]]$named_location) &
        #   Location %in% unique(site_list[[i]]$Location) &
          !is.na(flower_50_doy) &
          !is.na(yield_lb_acre)) %>% 
      mutate(flower_50_doy_sc = scale(flower_50_doy, scale = T))

    # # Exclude sites that have a flowering day spread of less than 5 days
    # exclude_low_spread <- site_dat_i %>%
    #   group_by(Year, named_location) %>%
    #   summarize(spr = max(flower_50_doy, na.rm = T) - min(flower_50_doy, na.rm = T) + 1) %>%
    #   arrange(spr) %>%
    #   filter(spr < 5)
    # site_dat_i <- site_dat_i %>% anti_join(., exclude_low_spread, by = join_by(Year, named_location))

    # # Exclude sites that only have one year
    # exclude_singletons <- site_dat_i %>%
    #   # add_count(Year, named_location)
    #   group_by(Year) %>%
    #   summarise(n = n())
    # if(nrow(exclude_singletons) <= 1) return(NA)

      group <- names(site_list[i])
      State <- unique(site_dat_i$State)
    
    for (year_i in unique(site_dat_i$Year)) {
      
      # subset location data into one year, to fit the model on each year separately
      site_dat_i_yr <- site_dat_i %>% 
        filter(Year == year_i)
      
      # site_i_nm <- unique(site_dat_i$Location) %>% as.character()
      named_loc_i <- ifelse(length(unique(site_dat_i_yr$named_location)) == 1,
                            unique(site_dat_i_yr$named_location),
                            paste(unique(site_dat_i_yr$named_location)[1],
                                  unique(site_dat_i_yr$named_location)[2], sep = ' & ')) %>% 
        as.character()
      
      
      # fit the gaussian model
      glm_tmp <- glm(yield_lb_acre ~  
                       flower_50_doy_sc + 
                       I(flower_50_doy_sc^2),
                     family = gaussian(link = "log"),
                     data = site_dat_i_yr)
      
      # grab coefficients and name them
      coefs_i <- coefficients(glm_tmp)
      names(coefs_i) <- c('beta0sc', 'beta1sc', 'beta2sc')
      
      # unscale the coefficients
      coef.unsc <- coefs_unscaled(coefs_i,
                                  mean = mean(site_dat_i$flower_50_doy),
                                  sd = sd(site_dat_i$flower_50_doy))  
      
      # see if the model fit a curve that's not gaussian
      badfit_i = (coef.unsc[2] < 0) | (coef.unsc[3] > 0)
      
      #calculate McFadden's psuedo R2
      pseudoR2_i <- round(with(summary(glm_tmp), 1 - deviance/null.deviance), digits = 4)
      
      
      
      # fit the linear model
      lm_tmp <- glm(yield_lb_acre ~  
                     flower_50_doy_sc,
                   family = gaussian(link = "log"),
                   data = site_dat_i_yr)
      
      # grab coefficients and name them
      coefs_lm_i <- coefficients(lm_tmp)
      names(coefs_lm_i) <- c('beta0sc', 'beta1sc')
      
      # unscale the coefficients
      coef.lm.unsc <- lin_coef_unscale(coefs_lm_i,
                                       mean = mean(site_dat_i$flower_50_doy),
                                       sd = sd(site_dat_i$flower_50_doy))  
      
      # is beta1 (slope) significant?
      beta1_lm_sig = summary(lm_tmp)$coefficients[2,4]
      
      # extract R2 from linear model
      # R2_lm_i <- round(summary(lm_tmp)$r.squared, digits = 4)
      R2_lm_i <- round(with(summary(lm_tmp), 1 - deviance/null.deviance), digits = 4)
      
      # find model preferred by AIC
      preferred_mod <- ifelse(AIC(glm_tmp) < AIC(lm_tmp), 'gaussian', 'linear')
      
      # which model has higher r2
      higher_R2 <- ifelse(pseudoR2_i> R2_lm_i, 'gaussian', 'linear')
      
      n_site_yr <- nrow(site_dat_i_yr)
      
      spread_flower_site_yr <- as.numeric(max(site_dat_i_yr$flower_50_doy, na.rm = T) - 
        min(site_dat_i_yr$flower_50_doy, na.rm = T) + 1)
      
      # Which model should optimal flowering time be based on?
      # if badfit = FALSE --> Gaussian
      # if badfit = TRUE --> check delta AIC, check p val
      
      mod_for_cens <- ifelse(badfit_i == FALSE, 'gauss',
                             ifelse(AIC(glm_tmp) - AIC(lm_tmp) < -2, 'exclude',
                                    ifelse(beta1_lm_sig > 0.05, 'exclude',
                                           'linear')))
      
      # generate the output for this site-year
      output <- c(year_i, 
                  # as.character(site_i_nm), 
                  named_loc_i,
                  group,
                  State,
                  unique(site_dat_i_yr$lat), 
                  unique(site_dat_i_yr$lon), 
                  n_site_yr,
                  spread_flower_site_yr,
                  pseudoR2_i,
                  badfit_i,
                  round(coef.unsc[1], digits = 8),
                  round(coef.unsc[2], digits = 8),
                  round(coef.unsc[3], digits = 8),
                  round(AIC(glm_tmp),4),
                  R2_lm_i,
                  round(coef.lm.unsc[1], digits = 4),
                  round(coef.lm.unsc[2], digits = 4),
                  round(beta1_lm_sig, digits = 4),
                  round(AIC(lm_tmp),4),
                  preferred_mod,
                  mod_for_cens,
                  higher_R2
                  
      )
      
      # bind the site-yr data to the big data set
      R2_out <- rbind(R2_out, output)
      
    }
    
  }
  
  # clean up and assign appropriate classes to columns
  R2_out <- R2_out %>% filter(!is.na(Year))
  R2_out$Year <- as.numeric(R2_out$Year)
  R2_out$lat <- as.numeric(R2_out$lat)
  R2_out$lon <- as.numeric(R2_out$lon)
  R2_out$pseudoR2_gau <- as.numeric(R2_out$pseudoR2_gau)
  R2_out$beta0_gau <- as.numeric(R2_out$beta0_gau)
  R2_out$beta1_gau <- as.numeric(R2_out$beta1_gau)
  R2_out$beta2_gau <- as.numeric(R2_out$beta2_gau)
  R2_out$R2_lm <- as.numeric(R2_out$R2_lm)
  R2_out$beta0_lm <- as.numeric(R2_out$beta0_lm)
  R2_out$beta1_lm <- as.numeric(R2_out$beta1_lm)
  R2_out$beta1_lm_sig <- as.numeric(R2_out$beta1_lm_sig)
  R2_out$AIC_gau <- as.numeric(R2_out$AIC_gau)
  R2_out$AIC_lm <- as.numeric(R2_out$AIC_lm)
  R2_out$n_site_yr <- as.numeric(R2_out$n_site_yr)
  R2_out$spread_flower_site_yr <- as.numeric(R2_out$spread_flower_site_yr)
  
  
  R2_out <- R2_out %>% mutate(AIC_delta = AIC_gau - AIC_lm,
                              n_predictors = 2,
                              site_yr = paste(Year, named_location, sep = "_"),
                              decade = ifelse(Year >= 1970 & Year < 1980, "1",
                                              ifelse(Year >= 1980 & Year < 1990, "2",
                                                     ifelse(Year >= 1990 & Year < 2000, "3",
                                                            ifelse(Year >= 2000 & Year < 2010, "4",
                                                                   ifelse(Year >= 2010 & Year < 2020, '5',
                                                                          ifelse(Year >= 2020 & Year < 2030, '6',
                                                                                 'timeless'))))))
                              )
  
  return(R2_out)
}

# Write out R2 and AIC data ------------------------------------------------------

# # common sites, exact location
# R2_exactloc_data <- list_sites(deriv_data, type = 'common_exact') %>% 
#   R2_data()
# R2_exactloc_data <- R2_exactloc_data %>% mutate(subset = "common_exactloc")



# all sites, no grouping - exact location
R2_allexact_data <- list_sites(deriv_data, type = 'all_exact') %>%
  R2_data()
R2_allexact_data <- R2_allexact_data %>% mutate(subset = "all_exactloc")

# Exclude trials with low spread of flowering
# low spread
# R2_allexact_data %>% arrange(spread_flower_site_yr)
R2_allexact_data <- R2_allexact_data %>%
  filter(spread_flower_site_yr >= 5)



# all sites, grouped by county
R2_allcounty_data <- list_sites(deriv_data, type = 'all_county', exclude_double = TRUE) %>%
  R2_data()
R2_allcounty_data <- R2_allcounty_data %>% mutate(subset = "all_countyloc")


# Exclude trials with low spread of flowering - for county data
# low spread
# R2_allcounty_data %>% arrange(spread_flower_site_yr)
R2_allcounty_data <- R2_allcounty_data %>%
  filter(spread_flower_site_yr >= 5)


# check to see what is different between optim_clim and R2_allcounty_data
# this can be uncommented and run, after you run 02_gaus_opt_flow_time.R or 
# load the optim_county data from the resulting csv

# anti_join(optim_county, R2_allcounty_data, by = join_by('county_state' == 'group', Year))
# #there should be 10 sites here - the 10 singletons excluded when run in 01_gaus_opt_flow_time.R
# anti_join(R2_allcounty_data, optim_county, by = join_by('group' == 'county_state', Year))


R2_data_out <- rbind(
  # R2_exactloc_data, 
  R2_allexact_data, 
  R2_allcounty_data)

write.csv(R2_data_out, file.path(deriv_data_filepath, "R2_AIC_data.csv"),
          row.names = FALSE)

# write out the full dataset that goes into making R2 for counties
county_trials <- do.call(rbind, list_sites(deriv_data, type = 'all_county', exclude_double = TRUE) )

R2_indata <- deriv_data %>% 
  semi_join(., county_trials, by = c('Year', 'Location', 'named_location')) %>%
  filter(
    !is.na(flower_50_doy) &
      !is.na(yield_lb_acre)) 
# Exclude sites that have a flowering day spread of less than 5 days
exclude_low_spread <- R2_indata %>%
  group_by(Year, named_location) %>%
  summarize(spr = max(flower_50_doy, na.rm = T) - min(flower_50_doy, na.rm = T) + 1) %>%
  arrange(spr) %>%
  filter(spr < 5)
R2_indata <- R2_indata %>% anti_join(., exclude_low_spread, by = join_by(Year, named_location))

write.csv(R2_indata, file.path(deriv_data_filepath, "R2_indata.csv"),
          row.names = FALSE)


# write out the full dataset that goes into making R2 for exact locations
exact_trials <- do.call(rbind, list_sites(deriv_data, type = 'all_exact') )

R2_indata_exact <- deriv_data %>% 
  semi_join(., exact_trials, by = c('Year', 'Location', 'named_location')) %>%
  filter(
    !is.na(flower_50_doy) &
      !is.na(yield_lb_acre)) 
# Exclude sites that have a flowering day spread of less than 5 days
exclude_low_spread_exact <- R2_indata_exact %>%
  group_by(Year, named_location) %>%
  summarize(spr = max(flower_50_doy, na.rm = T) - min(flower_50_doy, na.rm = T) + 1) %>%
  arrange(spr) %>%
  filter(spr < 5)
R2_indata_exact <- R2_indata_exact %>% anti_join(., exclude_low_spread_exact, by = join_by(Year, named_location))

write.csv(R2_indata_exact, file.path(deriv_data_filepath, "R2_indata_exact.csv"),
          row.names = FALSE)


# Run a meta-analysis ---------------------------------------------------

# Packages --------------------------------------------------------------
library(tidyverse)
library(metafor)
library(devtools)
# library(brms)
# library(tidybayes)
# library(broom)
# library(marginaleffects)

devtools::install_github("daniel1noble/orchaRd", ref = "main", force = TRUE)
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, orchaRd, emmeans,
               ape, phytools, flextable)

# Load Data ------------------------------------------------------------
# This is the R2 output csv of the four different ways to group the data
R2_full_data <- read.csv(paste0(deriv_data_filepath,"/R2_AIC_data.csv"))


# Calculate Adjusted R2 ---------------------------------------------------

R2_full_data <- R2_full_data %>% 
  mutate(
  adjR2_gau = 1 - (1 - pseudoR2_gau) * (n_site_yr - 1) / (n_site_yr - n_predictors - 1), 
  adjR2_gau2 = ifelse(adjR2_gau <= 0, 0.00001, adjR2_gau),
  adj_R2_gau_beta = ((adjR2_gau + 1)/2),
  .after = pseudoR2_gau
)

quantile(R2_full_data$adjR2_gau, probs = c(0.025, 0.5,  0.975))

# Calculate effect sizes and variance of effect ---------------------------
# plain R2
R2_full_data <- escalc(
  measure = "R2",
  r2i = pseudoR2_gau,
  mi = n_predictors,
  ni = n_site_yr,
  slab = site_yr,
  data = R2_full_data,
  append = T
)

# plain R2 on adjusted R2
R2_full_data <- escalc(
  measure = "R2",
  r2i = adjR2_gau2,
  mi = n_predictors,
  ni = n_site_yr,
  slab = site_yr,
  data = R2_full_data,
  append = T,
  var.names = c("yi.adj", "vi.adj")
)


# Z transformed R2
R2_full_data <- escalc(
  measure = "ZR2",
  r2i = pseudoR2_gau,
  mi = n_predictors,
  ni = n_site_yr,
  slab = site_yr,
  data = R2_full_data,
  append = T,
  var.names = c("yi.z", "vi.z")
)

# Z transformed on adjusted R2
R2_full_data <- escalc(
  measure = "ZR2",
  r2i = adjR2_gau2,
  mi = n_predictors,
  ni = n_site_yr,
  slab = site_yr,
  data = R2_full_data,
  append = T,
  var.names = c("yi.adj.z", "vi.adj.z")
)

# add this column of 1s, if you don't want to weight by sample size
R2_full_data <- R2_full_data %>%
  mutate(W = rep(1, nrow(R2_full_data)))


# Meta-Regressions -----------------------------------------------------
#### MUST SUBSET DATA HERE!!!!!!!
# R2_full_data MUST BE SUBSET BY COLUMN 'subset' BEFORE RUNNING A META-ANALYSIS!!!



###
## All sites, grouped by EXACT LOCATION ---------------------------------------
#subset
R2_allexact_data <- R2_full_data %>% filter(subset == 'all_exactloc') %>%
  # mutate(county = group, .after = group) %>%
  # add scaled parameters
  mutate(Year_sc = scale(Year, scale = T),
         n_site_yr_sc = scale(n_site_yr),
         n_site_yr2_sc = scale(I(n_site_yr^2)))

# global average 
ave_r2_allexact <- rma.mv(yi ~ n_site_yr_sc + n_site_yr2_sc, vi, 
                         random = ~ 1 | group,
                         data = R2_allexact_data)
summary(ave_r2_allexact)
orchaRd::orchard_plot(ave_r2_allexact, mod = "1", group = 'group', xlab = "R2")
orchaRd::caterpillars(ave_r2_allexact, mod = '1', group = 'group', xlab = "R2")

# mods
## Year
year_r2_allexact <- rma.mv(yi ~ Year_sc + n_site_yr_sc + n_site_yr2_sc, V = vi, 
                          random = ~ 1 | group, 
                          data = R2_allexact_data)
summary(year_r2_allexact)
orchaRd::mod_results(year_r2_allexact, mod = 'Year_sc', group = 'group') %>%
  orchaRd::bubble_plot(mod = 'Year_sc', group = 'group', xlab = 'Year', ylab = 'R2') +
  geom_abline(slope = 0, intercept = 0, col = 'black', lty = 2) +
  # theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

# Customized plot of results
orchaRd::mod_results(year_r2_allexact, mod = 'Year_sc', group = 'group') %>%
  ggplot() +
  geom_jitter(aes(x = data$moderator, y = data$yi, size = 1/sqrt(data$vi)), 
             shape = 21, color = 'grey50', fill = 'grey90', alpha = 0.95,
             position = position_jitter(height = 0, width = 0.2)) +
  geom_line(aes(x = mod_table$moderator, y = mod_table$estimate), 
            linewidth = 1) +
  geom_ribbon(aes(x = mod_table$moderator, 
                  ymin = mod_table$lowerCL, ymax = mod_table$upperCL),
              alpha = 0.2) +
  labs(y = bquote(R^2), x = "Year", size = "Precision (1/SE)") +
  # guides(size = guide_legend(position = 'bottom')) +
  geom_abline(slope = 0, intercept = 0, col = 'black', lty = 2) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = 'horizontal',
        # legend.position.inside = c(0.7,0.88),
        # legend.justification.inside = c('left','top'),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        # legend.background = element_rect(color = 'grey90')
        ) 
# ggsave("temp_plots/opt_manuscript_figures/year_r2_allexact.png", 
#        width = 6, height = 4)
  

## State
R2_allexact_data$State <- factor(R2_allexact_data$State, 
                                  levels = c("TX", "CO", "WY", "KS", "NE", "SD", "MN", "ND"))
state_r2_allexact <- rma.mv(yi ~ State + n_site_yr_sc + n_site_yr2_sc, V = vi, 
                           random = ~ 1 | group, 
                           data = R2_allexact_data)
summary(state_r2_allexact)
orchaRd::mod_results(state_r2_allexact, mod = "State", group = "group")
orchaRd::caterpillars(state_r2_allexact, mod = 'State', group = 'group', xlab = "R2")

orchaRd::orchard_plot(state_r2_allexact, mod = 'State', group = 'group', 
                      xlab = bquote(R^2), angle = 0) +
  labs(ylab = "State") +
  theme_bw(base_size = 16) +
  guides(size = guide_legend(position = 'inside')) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = 'horizontal',
        legend.position.inside = c(0.75,0.035),
        # legend.justification.inside = c('right','bottom'),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        legend.background = element_blank()
        )
# ggsave("temp_plots/opt_manuscript_figures/state_r2_allexact.png", 
#        width = 6.5, height = 5)




## All Sites, COUNTY grouping -------------------------------------------------
#subset
R2_allcounty_data <- R2_full_data %>% filter(subset == 'all_countyloc') %>%
  mutate(county = group, .after = group) %>%
  # add scaled parameters
  mutate(Year_sc = scale(Year, scale = T),
         n_site_yr_sc = scale(n_site_yr),
         n_site_yr2_sc = scale(n_site_yr^2))


# sample size bias
ggplot(R2_allcounty_data) +
  geom_point(aes(x = n_site_yr, y = pseudoR2_gau, color = "R2")) +
  geom_point(aes(x = n_site_yr, y = adjR2_gau, color = "adjR2")) +
  geom_smooth(aes(x = n_site_yr, y = pseudoR2_gau, color = "R2"), 
              formula = y ~ log(x), method = 'glm',  se = T) +
  geom_smooth(aes(x = n_site_yr, y = adjR2_gau, color = "adjR2"), 
              formula = y ~ log(x), method = 'glm', se = T) +
  geom_smooth(aes(x = n_site_yr, y = adjR2_gau2, color = "adjR2_0"), 
              formula = y ~ log(x), method = 'glm', se = T) +
  geom_point(aes(x = n_site_yr, y = adjR2_gau2, color = "adjR2_0")) 



### REGULAR R2 ----------------------------------------------------------
# non-adjusted
# graph of sample size bias
ggplot(R2_allcounty_data) +
  geom_point(aes(x = n_site_yr, y = pseudoR2_gau), alpha = 0.3) +
  geom_smooth(aes(x = n_site_yr, y = pseudoR2_gau), 
              formula = y ~ log(x), method = 'glm',  se = T, color = 'black') +
  labs(x = "Sample size", y = bquote(R^2)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
ggsave(paste0(figure_filepath,"/county/R2_sample_size.png"), 
       width = 6, height = 4)

# global average 
ave_r2_allcounty2 <- rma.mv(yi ~ n_site_yr_sc + n_site_yr2_sc, vi, 
                           random = ~ 1 | group,
                           data = R2_allcounty_data)
summary(ave_r2_allcounty2)
# plot(ave_r2_allcounty2)
orchaRd::orchard_plot(ave_r2_allcounty2, mod = "1", group = 'group', xlab = "R2")
orchaRd::caterpillars(ave_r2_allcounty2, mod = '1', group = 'group', xlab = "R2")
# 
# new_mods <- R2_allcounty_data %>% 
#   select(n_site_yr_sc, n_site_yr2_sc) %>% as.matrix()
# predict(ave_r2_allcounty2, newmods = new_mods, addx = TRUE) %>% as.data.frame() %>%
#   cbind(R2_allcounty_data %>% select(yi)) %>%
#   ggplot(aes(y = pred, x = yi)) +
#   geom_point( alpha = 0.3)

# global average, excluding Western Colorado & Coastal Texas to check for high leverage
ave_r2_allcounty2_leverage <- rma.mv(yi ~ n_site_yr_sc + n_site_yr2_sc, vi, 
                            random = ~ 1 | group,
                            data = R2_allcounty_data %>% 
                              filter(group != "Montezuma_CO") %>%
                              filter(group != "Nueces_TX"))
summary(ave_r2_allcounty2_leverage)
# # dot plot histogram
# ggplot() +
#   geom_dotplot(data = R2_allcounty_data, aes(x = yi, fill = group)) +
#   geom_dotplot(data = R2_allcounty_data %>% 
#                  filter(group == "Montezuma_CO") %>%
#                  filter(group == "Nueces_TX"), 
#                aes(x = yi), color = 'red') 


## Combine Year and State Effects
R2_allcounty_data$State <- factor(R2_allcounty_data$State, 
                                  levels = c("TX", "CO", "WY", "KS", "NE", "SD", "MN", "ND"))

state_Year_r2_allcounty2 <- rma.mv(yi ~ Year_sc + State + n_site_yr_sc + n_site_yr2_sc, V = vi, 
                                   random = ~ 1 | group/State, 
                                   data = R2_allcounty_data)
summary(state_Year_r2_allcounty2)

# constants to unscale Year
yr_sd_unsc <- sd(R2_allcounty_data$Year)
yr_mean_unsc <- mean(R2_allcounty_data$Year)


# quick plot of state effect
orchaRd::orchard_plot(state_Year_r2_allcounty2,
                      mod = 'State', group = 'group', 
                      xlab = bquote(R^2), angle = 0) +
  labs(ylab = "State") +
  theme_bw(base_size = 16) +
  guides(size = guide_legend(position = 'inside')) +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.direction = 'horizontal',
        legend.position.inside = c(0.25,0.035),
        # legend.justification.inside = c('left','bottom'),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        legend.background = element_blank()
  )

# how much does R2 change per decade?
coef(state_Year_r2_allcounty2)[2] * yr_sd_unsc * 10

# quick plot of year effect
orchaRd::mod_results(state_Year_r2_allcounty2, mod = 'Year_sc', 
                     # at = list("scale(n_site_yr)" = 0, "scale(I(n_site_yr^2))" = 0), 
                     group = 'group') %>%
  orchaRd::bubble_plot(mod = 'Year_sc', group = 'group', xlab = 'Year_sc', ylab = 'R2') +
  geom_abline(slope = 0, intercept = 0, col = 'black', lty = 2) +
  # theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

## State figure
r2_state_plot <- orchaRd::orchard_plot(state_Year_r2_allcounty2,
                      mod = 'State', group = 'group', 
                      xlab = bquote(R^2), angle = 0) +
  labs(x = "State") +
  theme_bw(base_size = 16) +
  guides(size = guide_legend(position = 'inside', title = "Precision (1/SE)")) +
  theme(panel.grid = element_blank(),
        legend.position = "none", # change this to "inside" to replace the precision scale legend
        legend.direction = 'horizontal',
        legend.position.inside = c(0.25,0.035),
        legend.justification.inside = c('left','bottom'),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        legend.background = element_blank()
  ) 
# remove vertical line
r2_state_plot <- gginnards::delete_layers(r2_state_plot, match_type = "GeomHline")
# remove the sample size annotations from plot
r2_state_plot <- gginnards::delete_layers(r2_state_plot, match_type = "GeomText")
  
# add vertical line and CI at mean
r2_state_plot <- r2_state_plot +
  geom_hline(yintercept = 0.089, linetype = 2, linewidth = 0.7) +
  annotate("rect", ymin = 0.0514, ymax = 0.1264, xmin = -Inf, xmax = Inf,
           fill = "grey50", alpha = 0.4) 
r2_state_plot

ggsave(paste0(figure_filepath,"/county/state_r2_county_unadj.png"),
       width = 6.5, height = 5)

### Year figure
orchaRd::mod_results(state_Year_r2_allcounty2, mod = 'Year_sc', group = 'group') %>% 
  ggplot() +
  geom_jitter(aes(x = data$moderator * yr_sd_unsc + yr_mean_unsc, 
                  y = data$yi, size = 1/sqrt(data$vi)), 
              shape = 21, color = 'grey50', fill = 'grey90', alpha = 0.95,
              position = position_jitter(height = 0, width = 0.2)) +
  geom_line(aes(x = mod_table$moderator * yr_sd_unsc + yr_mean_unsc, 
                y = mod_table$estimate), 
            linewidth = 1) +
  geom_ribbon(aes(x = mod_table$moderator * yr_sd_unsc + yr_mean_unsc, 
                  ymin = mod_table$lowerCL, ymax = mod_table$upperCL),
              alpha = 0.2) +
  labs(y = bquote(R^2), x = "Year", size = "Precision (1/SE)") +
  # annotate("text", x = 1978, y = 0.95, label = "Year trend: -0.01 (-0.02, -0.006)", hjust = 0) +
  # guides(size = guide_legend(position = 'bottom')) +
  geom_abline(slope = 0, intercept = 0, col = 'black', lty = 2) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = 'horizontal',
        # legend.position.inside = c(0.7,0.88),
        # legend.justification.inside = c('left','top'),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        # legend.background = element_rect(color = 'grey90')
  ) 
ggsave(paste0(figure_filepath,"/county/year_r2_county_unadj.png"),
       width = 6, height = 4)




# Other versions of analyses - not used
## Year by itself
year_r2_allcounty2 <- rma.mv(yi ~ Year_sc + n_site_yr_sc + n_site_yr2_sc, V = vi,
                            random = ~ 1 | group,
                            data = R2_allcounty_data)
summary(year_r2_allcounty2)
plot(residuals(year_r2_allcounty2))
orchaRd::mod_results(year_r2_allcounty2, mod = 'Year_sc',
                     # at = list("scale(n_site_yr)" = 0, "scale(I(n_site_yr^2))" = 0),
                     group = 'group') %>%
  orchaRd::bubble_plot(mod = 'Year_sc', group = 'group', xlab = 'Year_sc', ylab = 'R2') +
  geom_abline(slope = 0, intercept = 0, col = 'black', lty = 2) +
  # theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

# how much does R2 change per decade?
coef(year_r2_allcounty2)[2] * yr_sd_unsc * 10
orchaRd::mod_results(year_r2_allcounty2,
                     at = list(Year_sc = c(0, 10/yr_sd_unsc, 20/yr_sd_unsc)), by = "Year_sc", group = 'group')$mod_table %>%
  mutate(Year = condition * yr_sd_unsc + yr_mean_unsc)




# State
R2_allcounty_data$State <- factor(R2_allcounty_data$State,
                                  levels = c("TX", "CO", "WY", "KS", "NE", "SD", "MN", "ND"))
state_r2_allcounty2 <- rma.mv(yi ~ State + n_site_yr_sc + n_site_yr2_sc, V = vi,
                             random = ~ 1 | group,
                             data = R2_allcounty_data)
summary(state_r2_allcounty2)
orchaRd::mod_results(state_r2_allcounty2, mod = "State", group = "group")
orchaRd::caterpillars(state_r2_allcounty2, mod = 'State', group = 'group', xlab = "R2")
orchaRd::orchard_plot(state_r2_allcounty2, mod = "State", group = 'group', xlab = "R2")





## Latitude
lat_r2_allcounty2 <- rma.mv(yi ~ lat + n_site_yr_sc + n_site_yr2_sc, V = vi, 
                              random = ~ 1 | group, 
                              data = R2_allcounty_data %>% filter(lat > 35))
summary(lat_r2_allcounty2)
orchaRd::mod_results(lat_r2_allcounty2, mod = "lat", group = "group")
orchaRd::mod_results(lat_r2_allcounty2, mod = 'lat', group = 'group') %>%
  orchaRd::bubble_plot(mod = 'lat', group = 'group', xlab = 'Latitude', ylab = 'R2') +
  geom_abline(slope = 0, intercept = 0, col = 'black', lty = 2) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())





### ADJUSTED R2 --------------------------------------------------------------
# global average 
ave_r2_allcounty <- rma.mv(yi.adj, vi.adj, 
                           random = ~ 1 | group,
                           data = R2_allcounty_data)
summary(ave_r2_allcounty)
orchaRd::orchard_plot(ave_r2_allcounty, mod = "1", group = 'group', xlab = "R2")
orchaRd::caterpillars(ave_r2_allcounty, mod = '1', group = 'group', xlab = "R2")

## Year
year_r2_allcounty <- rma.mv(yi.adj ~ Year, V = vi.adj, 
                            random = ~ 1 | group, 
                            data = R2_allcounty_data)
summary(year_r2_allcounty)
orchaRd::mod_results(year_r2_allcounty, mod = 'Year', group = 'group') %>%
  orchaRd::bubble_plot(mod = 'Year', group = 'group', xlab = 'Year', ylab = 'R2') +
  geom_abline(slope = 0, intercept = 0, col = 'black', lty = 2) +
  # theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())


## State
state_r2_allcounty <- rma.mv(yi.adj ~ State, V = vi.adj, 
                             random = ~ 1 | group, 
                             data = R2_allcounty_data)
summary(state_r2_allcounty)
orchaRd::mod_results(state_r2_allcounty, mod = "State", group = "group")
orchaRd::caterpillars(state_r2_allcounty, mod = 'State', group = 'group', xlab = "R2")
orchaRd::orchard_plot(state_r2_allcounty, mod = "State", group = 'group', xlab = "R2")





# ## R2 and OPT_RELATION optimal (not used) -----------------------------------
# # Not used below here
# R2_allcounty_data %>% left_join(optim_county %>% select(county_state, Year, opt_relation_to_data),
#                                 by = c("county" = "county_state", "Year")) %>% 
#   # lme4::lmer(yi.adj ~ opt_relation_to_data + n_site_yr_sc + n_site_yr2_sc + (1|county), data = .) %>%
#   # emmeans::emmeans(pairwise~opt_relation_to_data)
#   group_by(opt_relation_to_data) %>%
#   summarise(R2 = mean(yi.adj),
#             n = n(),
#             ave_n_site_yr = mean(n_site_yr))
# 





