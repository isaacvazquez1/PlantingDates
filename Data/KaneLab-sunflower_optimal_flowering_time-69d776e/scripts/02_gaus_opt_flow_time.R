#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
############    Purpose: Predict and plot the optimal
############     flowering time of sunflower
############ Partially adapted from: Edwards & Crone, 2021
############.            
############.
############            By: Eliza Clark
############       Last modified: 09/10/2024
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


# Packages ----------------------------------------------------------------

library(MASS)
library(tidyverse) 
library(ggpubr)
library(lme4)
library(lmerTest)
select <- dplyr::select

# Custom Functions 
source("scripts/00_analysis_functions.R")

deriv_data_filepath <- "derived_data"
figure_filepath <- "figures"

# Data --------------------------------------------------------------------

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
  
# deriv_data %>% filter(grepl('Hitchcock', named_location)) %>%
#   select(named_location, Location, garden_city, county_state, Year) %>% distinct() %>% View()


## AIC/R2 Data
AIC_R2_data <- read.csv(paste0(deriv_data_filepath,"/R2_AIC_data.csv"))

### must subset AIC/R2 data to the correct subset before doing analyses!
### Subsets: 'all_exactloc', 'all_countyloc'



# Functions ----------------------------------------------------------------


## FUN: Fit Gaussian Model -------------------------------------------------
# the function gaussfit takes the name of a site, filters the data to include
# only datapoints with flow_50_doy and yield_lb_acre and runs a gaussian model
# on each year of the data. Then it calculates standard errors/CI with a 
# bootstrap method. It also calculates predicted yield across the range of 
# the flowering doy data for that site. It returns a list of three dataframes
# The first is deriv_site, which is the data for one site that was used to 
# do everything. This is useful for plotting. The second list element is 
# pheno_metrics, which contains all of the calculated metrics. Most importantly
# in this dataframe are mu (optimum), mu_## (standard errors and confidence
# intervals from the bootstrap), and y2, y3, cens, and y_na are variables 
# created that deal with whether the optimum is censored (above or below the
# range of data). These are used to fit brms models.

# # Development/Troubleshooting
# data = deriv_data
# # Create list of sites:
# site = sites_county[[25]]
# # AIC subset
# AIC = R2_county_data
# loc.type = 'named_location'

# # # #
# # site = "Carrington"
# # site = sites_all_exact[[20]]
# # AIC = R2_allexac_data
# # loc.type = "Location"
# # # #
# # site = sites_common_exact[[4]]
# loc.type = 'named_location'
# # AIC <- R2_commonexac_data



gaussfit <- function(data, site, AIC, loc.type) {
  
  
  ### Smaller dataset for convenience
  if(loc.type == "Location") {
    deriv_site <- data %>% 
      filter(Location == site &
               !is.na(flower_50_doy) &
               !is.na(yield_lb_acre) &
               !grepl('irrigated|NO|late|recrop|additional', Location))
  } else if(loc.type == "named_location"){
    deriv_site <- data %>%
      semi_join(., site, by = c('Year', 'Location', 'named_location')) %>%
      filter(
        # named_location %in% site$named_location &
        #        Location %in% site$Location &
        #        lat %in% site$lat &
        #        lon %in% site$lon &
               !is.na(flower_50_doy) &
               !is.na(yield_lb_acre))
  } else print("Check that loc.type is an accepted value")
  
  # Exclude sites that have a flowering day spread of less than 5 days
  exclude_low_spread <- deriv_site %>%
    group_by(Year, named_location) %>% 
    summarize(spr = max(flower_50_doy) - min(flower_50_doy) + 1) %>% 
    arrange(spr) %>% 
    filter(spr < 5)
  deriv_site <- deriv_site %>% anti_join(., exclude_low_spread, by = join_by(Year, named_location))
  
  # Exclude sites that only have one year
  exclude_singletons <- deriv_site %>%
    # add_count(Year, named_location)
    group_by(Year) %>%
    summarise(n = n())
  if(nrow(exclude_singletons) <= 1) return(NA)

  # create scaled flowering day
  deriv_site$flower_50_doy_sc <- scale(deriv_site$flower_50_doy, scale = T)
  
  # remove this site-yr that doesn't run
  if("Cheyenne_HPAL" %in% unique(deriv_site$named_location)) {
    deriv_site <- deriv_site %>% filter(Year != 2013)
  }
  
  ### Smaller AIC dataset
  AIC_site <- AIC %>% filter(named_location %in% site$named_location)
  
  
  ####
  # Fit Gauss Model ----------
  
  
  mod.lin <- glm(yield_lb_acre ~ 0 + Year + 
                   flower_50_doy_sc:Year + 
                   I(flower_50_doy_sc^2):Year,
                 family = gaussian(link = "log"),
                 data = deriv_site)
  
  # grabbing coefficients
  coef.nam <- names(coef(mod.lin))
  
  # sanity check to make sure there are 3 coeffs per year
  if(any(table(gsub(":.*", "", coef.nam))!=3)){
    print('Check your coefficients! Seems wrong!')
  }  
  
  coef.uniq <- unique(gsub(":.*", "", coef.nam))
  
  reg.store = NULL
  for(cur.coef in coef.uniq){
    ind = grep(cur.coef, coef.nam)
    reg.store = rbind(reg.store,
                      coefficients(mod.lin)[ind])
  }
  
  rownames(reg.store) <- coef.uniq
  colnames(reg.store) <- c('beta0sc', 'beta1sc', 'beta2sc')
  
  # unscaled coefficients
  coef.unsc <- coefs_unscaled(reg.store,
                              mean = mean(deriv_site$flower_50_doy),
                              sd = sd(deriv_site$flower_50_doy))
  
  
  ####
  # Fit Linear Model ---------
  mod.lm <- glm(yield_lb_acre ~ 0 + Year + 
                 flower_50_doy_sc:Year,
               family = gaussian(link = "log"),
               data = deriv_site)
  
  # grabbing coefficients
  coef.nam.lm <- names(coef(mod.lm))
  
  # sanity check to make sure there are 3 coeffs per year
  if(any(table(gsub(":.*", "", coef.nam.lm))!=2)){
    print('Check your coefficients! Seems wrong!')
  }  
  
  coef.uniq.lm <- unique(gsub(":.*", "", coef.nam.lm))
  
  reg.store.lm = NULL
  for(cur.coef in coef.uniq.lm){
    ind = grep(cur.coef, coef.nam.lm)
    reg.store.lm = rbind(reg.store.lm,
                         coefficients(mod.lm)[ind])
  }
  
  rownames(reg.store.lm) <- coef.uniq.lm
  colnames(reg.store.lm) <- c('beta0sc', 'beta1sc')
  
  # unscaled coefficients
  coef.unsc.lm <- lin_coef_unscale_mat(reg.store.lm,
                                       mean = mean(deriv_site$flower_50_doy),
                                       sd = sd(deriv_site$flower_50_doy))
  
  # Check that two model methods are doing the same thing -------------
  equal_coefs_check <- AIC_site %>% select(Year, named_location, beta0_gau, beta1_gau, beta2_gau, beta0_lm, beta1_lm) %>%
    left_join(., coef.unsc.lm %>% as.data.frame() %>%
                mutate(Year = row.names(coef.unsc.lm) %>%
                         gsub(pattern = "Year", replacement = "")) %>%
                mutate(Year = as.numeric(Year))) %>%
    left_join(., coef.unsc %>% as.data.frame() %>%
                mutate(Year = row.names(coef.unsc) %>%
                         gsub(pattern = "Year", replacement = "")) %>%
                mutate(Year = as.numeric(Year))) %>%
    mutate(equal_coefs = ifelse(beta0_gau - beta0 < 0.01  &
                                  beta1_gau - beta1 < 0.01 &
                                  beta2_gau - beta2 < 0.01 &
                                  beta0_lm - lm_beta0 < 0.01 &
                                  beta1_lm - lm_beta1< 0.01,
                                "GOOD", "FALSE")) %>%
    filter(equal_coefs == FALSE)
  
  ####


  # variance covariance matrix for error estimation
  
  vcov.list = list()
  for (cur.coef in coef.uniq) {
    cur.cov = vcov(mod.lin)[grep(cur.coef, rownames(vcov(mod.lin))),
                            grep(cur.coef, colnames(vcov(mod.lin)))]
    vcov.list[[cur.coef]] <- cur.cov
  }
  
  # check to make sure columns are Year, doy, doy^2
  vcov.list[[1]]
  
  
  ## Predict model points --------------
  
  range_flower_doy = range(deriv_site$flower_50_doy)
  
  x.pred = seq(range_flower_doy[1], range_flower_doy[2], by = 0.1)
  
  pred_out = NA
  
  for(year.cur in unique(droplevels(deriv_site$Year))){
    
    pred_dat <- deriv_site[deriv_site$Year == year.cur,]
    #gaussian
    coef.cur = coef.unsc[rownames(coef.unsc) == paste0("Year", year.cur),]
    y.pred = exp(coef.cur[[1]] + coef.cur[[2]] * x.pred + coef.cur[[3]] * x.pred^2)
    #linear
    coef.cur.lm = coef.unsc.lm[rownames(coef.unsc.lm) == paste0("Year", year.cur),]
    y.pred.lm = exp(coef.cur.lm[[1]] + coef.cur.lm[[2]] * x.pred)
    
    Year = rep(year.cur, length(x.pred))
    pred_curyear <- cbind(Year, x.pred, y.pred, y.pred.lm)
    pred_out <- rbind(pred_out, pred_curyear)
  }
  
  pred_out <- na.omit(pred_out)
  pred_out <- data.frame(pred_out)
  pred_out$x.pred <- as.numeric(pred_out$x.pred)
  pred_out$y.pred <- as.numeric(pred_out$y.pred)
  pred_out$y.pred.lm <- as.numeric(pred_out$y.pred.lm)
  
  # str(pred_out)
  
  
  ## Calculate phenological metrics ---------
  
  pheno_metrics <- pheno_calc(coef.unsc)
  pheno_metrics$Year <- rownames(pheno_metrics)
  # str(pheno_metrics)
  pheno_metrics$Year <- gsub(pattern = "Year", replacement = "", x = as.character(pheno_metrics$Year))
  
  # figure out if estimate is within range of data
  range_flower_doy <- deriv_site %>% group_by (Year) %>%
    summarise(
      flower_min = min(flower_50_doy),
      flower_max = max(flower_50_doy),
      flower_median = median(flower_50_doy)
    )
  pheno_metrics <- pheno_metrics %>%
    mutate(outside_data = ifelse(between(mu, range_flower_doy$flower_min, range_flower_doy$flower_max),
                                 'in_range',
                                 ifelse(mu >= range_flower_doy$flower_median, "greater_than",
                                        "less_than")))
  
  ## Parametric bootstrapping ---------
  
  # place to store simulated data
  dat.sim = list()
  
  # place to store results
  boot.res = data.frame(Year = rownames(reg.store),
                        mu_025 = numeric(nrow(reg.store)),
                        mu_975 = numeric(nrow(reg.store)),
                        mu_16 = numeric(nrow(reg.store)),
                        mu_84 = numeric(nrow(reg.store)),
                        sigma_025 = numeric(nrow(reg.store)),
                        sigma_975 = numeric(nrow(reg.store)),
                        sigma_16 = numeric(nrow(reg.store)),
                        sigma_84 = numeric(nrow(reg.store)),
                        N_025 = numeric(nrow(reg.store)),
                        N_975 = numeric(nrow(reg.store)),
                        N_16 = numeric(nrow(reg.store)),
                        N_84 = numeric(nrow(reg.store)),
                        h_025 = numeric(nrow(reg.store)),
                        h_975 = numeric(nrow(reg.store)),
                        h_16 = numeric(nrow(reg.store)),
                        h_84 = numeric(nrow(reg.store)),
                        bad = numeric(nrow(reg.store))
  )
  
  nsample = 10000
  
  for (i in 1:nrow(reg.store)) {
    cur.year = rownames(reg.store)[i]
    cur.mean = reg.store[i,]
    
    cur.est = coefs_unscaled(cur.mean, mean = mean(deriv_site$flower_50_doy),
                             sd = sd(deriv_site$flower_50_doy))
    
    #make simulated data
    sim = mvrnorm(n = round(nsample*2.5),
                  mu = cur.mean,
                  Sigma = vcov.list[[cur.year]]
    )
    sim.unsc = coefs_unscaled(sim, mean = mean(deriv_site$flower_50_doy),
                              sd = sd(deriv_site$flower_50_doy)
    )
    boot.bad = sim.unsc[,2]<0 | sim.unsc[,3]>0
    
    if(mean(boot.bad)>0.6){  # break loop if more than 60% are bad
      boot.res[i,2:17] = NA      # and add NA to the rest of the columns
      boot.res$bad[i] = mean(boot.bad)
      
      next
    }
    
    boot.res$bad[i] = mean(boot.bad)
    
    sim.unsc = sim.unsc[!boot.bad,][1:nsample,]
    sim = sim[!boot.bad,][1:nsample,]
    
    sim.mu= -sim.unsc[,2]/(2*sim.unsc[,3])
    sim.sigma=sqrt(-1/(2*sim.unsc[,3]))
    sim.N=sqrt(2*pi/(-2*sim.unsc[,3]))*exp(sim.unsc[,1]+(sim.unsc[,2]^2)/(-4*sim.unsc[,3]))
    sim.h=exp(sim.unsc[,1] + sim.unsc[,2]*sim.mu + sim.unsc[,3]*sim.mu^2)
    
    boot.res$mu_025[i]=quantile(sim.mu, probs=0.025)
    boot.res$mu_975[i]=quantile(sim.mu, probs=0.975)
    boot.res$mu_16[i]=quantile(sim.mu, probs=pnorm(-1))
    boot.res$mu_84[i]=quantile(sim.mu, probs=pnorm(1))
    boot.res$sigma_025[i]=quantile(sim.sigma, probs=0.025)
    boot.res$sigma_975[i]=quantile(sim.sigma, probs=0.975)
    boot.res$sigma_16[i]=quantile(sim.sigma, probs=pnorm(-1))
    boot.res$sigma_84[i]=quantile(sim.sigma, probs=pnorm(1))
    boot.res$N_025[i]=quantile(sim.N, probs=0.025)
    boot.res$N_975[i]=quantile(sim.N, probs=0.975)
    boot.res$N_16[i]=quantile(sim.N, probs=pnorm(-1))
    boot.res$N_84[i]=quantile(sim.N, probs=pnorm(1))
    boot.res$h_025[i]=quantile(sim.h, probs=0.025)
    boot.res$h_975[i]=quantile(sim.h, probs=0.975)
    boot.res$h_16[i]=quantile(sim.h, probs=pnorm(-1))
    boot.res$h_84[i]=quantile(sim.h, probs=pnorm(1))
    
    dat.sim[[paste0(cur.year)]] = cbind(sim, sim.unsc, sim.mu, sim.sigma, sim.N, sim.h)
  }
  
  boot.res$Year <- gsub(pattern = "Year", replacement = "", x = as.character(boot.res$Year))
  
  ## Data Censoring -------------
  
  # making new columns for all estimate mu (optimal flowering doy)
  # mu = y = estimated optimal flowering doy
  # y2 = lower bound of censored interval
  # y3 = upper bound of censored interval
  # cens = type of censoring = "interval"
  
  # combine bootstrap standard error columns to pheno_metric dataframe (convenience)
  pheno_metrics_comb <- pheno_metrics %>% 
    left_join(., boot.res %>% 
                dplyr::select(Year, mu_025, mu_975, mu_16, mu_84)) %>%
    mutate(Year = as.integer(Year)) %>%
    left_join(., AIC_site %>% dplyr::select(Year, AIC_gau, AIC_lm, mod_for_cens, beta1_lm)) %>%
    mutate(Year = as.factor(Year))
    
  
  # combine planting and harvest dates with pheno_metric dataframe for lower
  # and upper bound calculations
  mean_site_harvest = mean(deriv_site$harvest_doy, na.rm = T)
  mean_site_planting = mean(deriv_site$planting_doy, na.rm = T)
  
  pheno_metrics_comb <- pheno_metrics_comb %>% 
    left_join(., deriv_site %>% 
                distinct(Year, planting_doy, harvest_doy, flower_50_doy) %>% arrange(Year) %>%
                group_by(Year) %>%
                summarise(
                  planting_doy = mean(planting_doy, na.rm = T),
                  harvest_doy = mean(harvest_doy, na.rm = T),
                  min_flow_50_doy = min(flower_50_doy, na.rm = T),
                  max_flow_50_doy = max(flower_50_doy, na.rm = T),
                ) %>%
                mutate(harvest_doy = ifelse(is.na(harvest_doy), mean_site_harvest, harvest_doy),
                       planting_doy = ifelse(is.na(planting_doy), mean_site_planting, planting_doy))
              )
  
  pheno_metrics_comb <- pheno_metrics_comb %>%
    mutate(
      # calculate lower bound of censored variable (y2)
      y2 = ifelse(
        # gauss doesn't fit and linear isn't sig
        test = mod_for_cens == 'exclude',
        yes = NA,
        no = ifelse(
          # linear and opt is before data
          test = mod_for_cens == 'linear' &
            beta1_lm < 0,
          yes = planting_doy + 40,
          no = ifelse(
            # linear and opt is after data
            test = mod_for_cens == 'linear' &
              beta1_lm > 0,
            yes = max_flow_50_doy,
            no = ifelse(
              # gaussian and opt is in range
              test = mod_for_cens == 'gauss' &
                outside_data == 'in_range',
              yes = mu_16,
              no = ifelse(
                # gaussian and opt is greater than range
                test = mod_for_cens == 'gauss' &
                  outside_data == 'greater_than',
                yes = mu_16,
                no = ifelse(
                  # gauss and opt is less than range
                  test = mod_for_cens == 'gauss' &
                    outside_data == 'less_than',
                  yes = planting_doy + 40,
                  no = 'HELP_ME'
                )
              )
            )
          )
        )
      ),
      # upper bound of censoring variable - y3
      y3 = ifelse(
        # gauss doesn't fit and linear isn't sig
        test = mod_for_cens == 'exclude',
        yes = NA,
        no = ifelse(
          # linear and opt is before data
          test = mod_for_cens == 'linear' &
            beta1_lm < 0,
          yes = min_flow_50_doy,
          no = ifelse(
            # linear and opt is after data
            test = mod_for_cens == 'linear' &
              beta1_lm > 0,
            yes = harvest_doy - 35,
            no = ifelse(
              # gaussian and opt is in range
              test = mod_for_cens == 'gauss' &
                outside_data == 'in_range',
              yes = mu_84,
              no = ifelse(
                # gaussian and opt is greater than range
                test = mod_for_cens == 'gauss' &
                  outside_data == 'greater_than',
                yes = harvest_doy - 35,
                # gauss and opt is less than range
                no = ifelse(
                  # gauss and opt is less than range
                  test = mod_for_cens == 'gauss' &
                    outside_data == 'less_than',
                  yes = mu_84,
                  no = 'HELP_ME'
                )
              )
            )
          )
        )
      ),
      # optimum/censoring intervals is before, after, or within data range
      opt_relation_to_data = ifelse(
        # gauss doesn't fit and linear isn't sig
        test = mod_for_cens == 'exclude',
        yes = 'unclear',
        no = ifelse(
          # linear and opt is before data
          test = mod_for_cens == 'linear' &
            beta1_lm < 0,
          yes = 'opt_before',
          no = ifelse(
            # linear and opt is after data
            test = mod_for_cens == 'linear' &
              beta1_lm > 0,
            yes = 'opt_after',
            no = ifelse(
              # gaussian and opt is in range
              test = mod_for_cens == 'gauss' &
                outside_data == 'in_range',
              yes = 'opt_within',
              no = ifelse(
                # gaussian and opt is greater than range
                test = mod_for_cens == 'gauss' &
                  outside_data == 'greater_than',
                yes = 'opt_after',
                # gauss and opt is less than range
                no = ifelse(
                  # gauss and opt is less than range
                  test = mod_for_cens == 'gauss' &
                    outside_data == 'less_than',
                  yes = 'opt_before',
                  no = 'HELP_ME'
                )
              )
            )
          )
        )
      ),
      cens = 'interval',
      y_na = ifelse(badfit == FALSE, mu, NA)
      
    ) 
  
  # old censoring logic
    #                    ifelse(outside_data == 'in_range', mu_16, 
    #                           ifelse(outside_data == 'greater_than', mu_16,
    #                                  yday(ymd(planting_date))+21))),
    #        y3 = ifelse(badfit == TRUE, NA,
    #                    ifelse(outside_data == 'in_range', mu_84, 
    #                           ifelse(outside_data == 'less_than', mu_84,
    #                                  yday(ymd(harvest_date))-21))),
    #        cens = 'interval',
    #        y_na = ifelse(badfit == FALSE, mu, NA)
    # )
  
  pheno_metrics_comb <- pheno_metrics_comb %>%
    mutate(y2 = as.numeric(y2),
           y3 = as.numeric(y3))
  
  # sanity check that y3 is greater than y2 - should be all TRUE or NA
  unique(pheno_metrics_comb$y3 > pheno_metrics_comb$y2)
  
  # Function Output -----------------
  out = list(deriv_site,
             pheno_metrics_comb,
             pred_out,
             equal_coefs_check,
             loc.type)
  names(out)[1] <- 'deriv_site'
  names(out)[2] <- 'pheno_metrics'
  names(out)[3] <- 'pred_out'
  names(out)[4] <- 'coefs_check'
  names(out)[5] <- 'loc.type'
  
  
  return(out)
}


## FUN: Plot Optimum & CI ----------------------------------------------------------
# The function plotopt takes the output from the function gaussfit and 
# plots the data, fitted curve, standard errors, and confidence intervals
# and saves a png of the plot. The plot will save in the working directory 
# (I think, this functionality not tested), or to the path you specify.

# # troubleshooting
# model = gauss_exactloc_out[[3]]
# path = path2plots
# 
# gauss_exactloc_out[[2]]$deriv_site %>% View()
# gauss_exactloc_out[[2]]$pheno_metrics %>% View()

plotopt <- function(model){
  
  # phenological/optimum flowering time curves
  if(model$loc.type == "Location") {
    site = unique(model$deriv_site$Location)
    site_lab = site
  } else if(model$loc.type == "named_location") {
    site = paste(unique(model$deriv_site$named_location), collapse = "")
    site_lab = unique(model$deriv_site$named_location)[1]
  }
  
  year_name <- model$deriv_site %>% select(Year, named_location) %>% distinct() %>% arrange(Year)
  year_labs <- paste0(year_name$Year, " at " , year_name$named_location)
  names(year_labs) <- year_name$Year %>% as.character()
  
  plot_opt_ci <- ggplot() +
    geom_point(data = model$deriv_site, 
               aes(x = flower_50_doy, y = yield_lb_acre), alpha = 0.5, size = 0.5) +
    geom_line(data = model$pred_out %>% 
                mutate(y.pred = ifelse(y.pred > 4000, NA, y.pred)), 
              aes(x = x.pred, y = y.pred), color = 'royalblue') +
    geom_rect(data = model$pheno_metrics, 
              aes(xmin = mu_025, xmax = mu_975, ymin=-Inf, ymax=Inf),
              fill = 'goldenrod1', alpha = 0.2) + # confidence interval
    geom_rect(data = model$pheno_metrics, 
              aes(xmin = mu_16, xmax = mu_84, ymin=-Inf, ymax=Inf),
              fill = 'tomato1', alpha = 0.2) + # SE
    geom_vline(data = model$pheno_metrics %>% filter(badfit == FALSE), aes( xintercept = mu),
               color = 'tomato3') +
    facet_wrap(~ Year, scales = 'free', labeller = labeller(Year = year_labs)) +
    labs(title = site) +
    theme_bw(base_size = 16)

  plot_opt_ci
    
  
}

## FUN: Plot Gaus & Linear curves ----------------------------------------------------------
# The function plotopt takes the output from the function gaussfit and 
# plots the data, fitted gaussian curve, and fitted line,
# and saves a png of the plot. The plot will save in the working directory 
# (I think, this functionality not tested), or to the path you specify.

# # Development
# model <- gauss_exactloc_out[[3]]
# AIC <- R2_commonexac_data

plotlingauss <- function(model, AIC){
  
  # phenological/optimum flowering time curves
  if(model$loc.type == "Location") {
    site = unique(model$deriv_site$Location)
    site_lab = site
  } else if(model$loc.type == "named_location") {
    site = paste(unique(model$deriv_site$named_location), collapse = "")
    site_lab = unique(model$deriv_site$named_location)[1]
  }
  
  year_name <- model$deriv_site %>% dplyr::select(Year, named_location) %>% distinct() %>% arrange(Year)
  year_labs <- paste0(year_name$Year, " at " , year_name$named_location)
  names(year_labs) <- year_name$Year %>% as.character()
  AIC_annot <- year_name %>% left_join(., AIC %>% dplyr::select(Year, named_location, AIC_gau, AIC_lm, mod_for_cens) %>%
                                         mutate(Year = as.factor(Year))) %>%
    mutate(label_g = paste0("AIC_gau = ",round(AIC_gau,2)),
           label_l = paste0("AIC_lm = ",round(AIC_lm,2)),
           # label_mod = paste0("chosen mod: ", mod_for_cens),
           color_g = ifelse(AIC_gau < AIC_lm, '1',
                            '2'),
           color_l = ifelse(AIC_gau > AIC_lm, '1',
                            '2')
           )
  yr_spread <- model$deriv_site %>% group_by(Year) %>% 
    summarize(spr = max(flower_50_doy) - min(flower_50_doy)) %>% 
    arrange(Year)
  
  
  plot_lin_gauss <- ggplot() +
    geom_point(data = model$deriv_site, 
               aes(x = flower_50_doy, y = yield_lb_acre), alpha = 0.5, size = 0.5) +
    geom_line(data = model$pred_out %>%
                mutate(y.pred = ifelse(y.pred > 4000, NA, y.pred)),
              aes(x = x.pred, y = y.pred), color = 'royalblue') +
    geom_line(data = model$pred_out %>% 
                mutate(y.pred.lm = ifelse(y.pred.lm > 4000, NA, y.pred.lm)) %>% 
                mutate(y.pred.lm = ifelse(y.pred.lm < 0, NA, y.pred.lm)), 
              aes(x = x.pred, y = y.pred.lm), color = 'mediumpurple3') +
    geom_text(data = AIC_annot, aes(x=200, y = -Inf, label = label_g, color = color_g),
              hjust = 0, vjust = -1) +
    geom_text(data = AIC_annot, aes(x=220, y = -Inf, label = label_l, color = color_l),
              hjust = 0, vjust = -1) +
    geom_text(data = AIC_annot, aes(x=200, y = Inf, label = mod_for_cens), color = 'red',
              hjust = 0, vjust = 1) +
    geom_text(data = yr_spread, aes(x=220, y = Inf-1000, label = spr), color = 'navy',
              hjust = -1, vjust = 1) +
    # annotate('text', x = Inf, y = 0.9, hjust = 1, vjust = 1,
    #          label = AIC_label(mod)) +
    geom_vline(data = model$pheno_metrics %>% filter(badfit == FALSE), aes( xintercept = mu),
               color = 'tomato3') +
    geom_vline(data = model$pheno_metrics, aes( xintercept = y2),
               color = 'darkorange2') +
    geom_vline(data = model$pheno_metrics, aes( xintercept = y3),
               color = 'darkorange2') +
    facet_wrap(~ Year, scales = 'free', labeller = labeller(Year = year_labs)) +
    scale_color_manual(values = c('forestgreen', 'black')) +
    guides(color = 'none') +
    labs(title = site) +
    theme_bw(base_size = 16)
  
  plot_lin_gauss
  
  
}


## FUN: Write out FILES ------------------------------------------------------------
# object = gauss_county_out
# file_name = "optimum_flowering_time.csv"
write_optima_exact <- function(object, file_name){
  optims <- lapply(object, FUN = function(x) x$pheno_metrics %>% cbind() %>% as.data.frame)
  optims <- do.call(rbind, optims)
  optims <- optims %>% mutate(named_location = rownames(.), .before = mu)
  optims <- optims %>% mutate(named_location = gsub(pattern = "\\..|\\...", 
                                                    replacement = "", x = named_location))
  rownames(optims) <- seq(1,nrow(optims),1)
  
  # original data for each site
  gaus_in_data <-lapply(object, FUN = function(x) x$deriv_site %>% cbind() %>% as.data.frame)
  gaus_in_data <- do.call(rbind, gaus_in_data)
  rownames(gaus_in_data) <- seq(1,nrow(gaus_in_data),1)
  gaus_in_data_uni <- gaus_in_data %>% select(named_location, Location, Year, State, lat, lon, county_state) %>% 
    distinct()
  
  # add in other important columns to optima dataset from original data
  optims <- optims %>% left_join(., gaus_in_data_uni,
                                 join_by(named_location, Year), 
                                 multiple = 'any'
  ) %>%
    relocate(Location, State, lat, lon, Year, planting_doy, harvest_doy, min_flow_50_doy, max_flow_50_doy,  .after = named_location) %>%
    rename(location = Location,
           latitude = lat,
           longitude = lon)
  
  # write optima data to a file
  write.csv(optims, file = paste0(deriv_data_filepath,"/", file_name))
}

# object = gauss_county_out
# file_name = "optimum_flowering_time.csv"
write_optima_county <- function(object, file_name){
  optims <- lapply(object, FUN = function(x) x$pheno_metrics %>% cbind() %>% as.data.frame)
  optims <- do.call(rbind, optims)
  optims <- optims %>% mutate(county_state = rownames(.), .before = mu)
  optims <- optims %>% mutate(county_state = gsub(pattern = "\\..|\\...", 
                                                    replacement = "", x = county_state))
  rownames(optims) <- seq(1,nrow(optims),1)
  
  # original data for each site
  gaus_in_data <-lapply(object, FUN = function(x) x$deriv_site %>% cbind() %>% as.data.frame)
  gaus_in_data <- do.call(rbind, gaus_in_data)
  rownames(gaus_in_data) <- seq(1,nrow(gaus_in_data),1)
  gaus_in_data_uni <- gaus_in_data %>% select(named_location, Location, Year, State, lat, lon, county_state) %>% 
    distinct()
  
  # add in other important columns to optima dataset from original data
  optims <- optims %>% left_join(., gaus_in_data_uni,
                                 join_by(county_state, Year), 
                                 multiple = 'any'
  ) %>%
    relocate(Location, State, lat, lon, Year, planting_doy, harvest_doy, min_flow_50_doy, max_flow_50_doy,  .after = named_location) %>%
    rename(location = Location,
           latitude = lat,
           longitude = lon)
  
  # write optima data to a file
  write.csv(optims, file = paste0(deriv_data_filepath,"/", file_name))
}

write_loc_lookup <- function(object, file_name){
  # original data for each site
  gaus_in_data <-lapply(object, FUN = function(x) x$deriv_site %>% cbind() %>% as.data.frame)
  gaus_in_data <- do.call(rbind, gaus_in_data)
  rownames(gaus_in_data) <- seq(1,nrow(gaus_in_data),1)
  gaus_in_data_uni <- gaus_in_data %>% select(named_location, county_state, Location, Year, State, lat, lon) %>% 
    distinct()
  
  # write original data to a file
  write.csv(gaus_in_data_uni, file = paste0(deriv_data_filepath,"/", file_name))
}





# Run Analysis  -------------------------------------------------------------
# What are reasonable ways to censor the optimal flowering time?

# Minimum time to flower ever is 41.3 days, so y2 can be planting day + 40
range(deriv_data$flower_50pct_days_past_planting, na.rm = T)

# Flowering can occur as little as 14 days before harvest, but most are 
# between 40-100 days before harvest. y3 can be censored harvest doy - 35
range(deriv_data$harvest_doy - deriv_data$flower_50_doy, na.rm = T)
hist(deriv_data$harvest_doy - deriv_data$flower_50_doy, breaks = 50)




# ## Common sites, exact location ----------------------------------------------
# 
# # Create list of sites:
# sites_common_exact <- list_sites(deriv_data, type = 'common_exact')
# 
# # Subset AIC/R2 data
# R2_commonexac_data <- AIC_R2_data %>% filter(subset == 'common_exactloc')
# 
# ### Gaussfit
# gauss_exactloc_out <- lapply(sites_common_exact, 
#                              gaussfit, 
#                              data = deriv_data, 
#                              AIC = R2_commonexac_data,
#                              loc.type = 'named_location')
# 
# # check models where coefs are not the same
# coef_check <- lapply(gauss_exactloc_out, FUN = function(x) x$coefs_check %>% cbind() %>% as.data.frame)
# coef_check <- do.call(rbind, coef_check)
# 
# ### Plotopt
# plot_opt_plots <- lapply(gauss_exactloc_out, plotopt)
# plot_opt_plots[[1]]
# 
# pdf("temp_plots/ec/optimal_flower_time.pdf", width = 25, height = 16, onefile = TRUE)
# plot_opt_plots
# dev.off()
# 
# ### Plotlingauss
# 
# plot_lingauss_plots <- lapply(gauss_exactloc_out, plotlingauss, AIC = R2_commonexac_data)
# plot_lingauss_plots[[5]]
# 
# pdf("temp_plots/ec/linear_gauss_curves.pdf", width = 16, height = 12, onefile = TRUE)
# plot_lingauss_plots
# dev.off()
# 
# 
# # # Write common-exact models to R object 
# # saveRDS(gauss_exactloc_out, file = "data_derived/gauss_exact_mods.RData")
# # 


## ALL sites, EXACT location ------------------------------------------------

# Create list of sites:
sites_all_exact <- list_sites(deriv_data, type = 'all_exact')
# AIC subset
R2_allexac_data <- AIC_R2_data %>% filter(subset == 'all_exactloc')


### Gaussfit
gauss_all_out <- lapply(sites_all_exact, 
                             gaussfit, 
                             data = deriv_data, 
                             AIC = R2_allexac_data,
                             loc.type = 'named_location')

# Remove elements of list that are NA (only have one site-yr)
gauss_all_out <- Filter(function(a) any(!is.na(a)), gauss_all_out)

# Write all-exact models to an r-file
saveRDS(gauss_all_out, file = paste0(deriv_data_filepath,"/gauss_allexact_mods.RData"))

## Plotlingauss
lingauss_all_plots <- lapply(gauss_all_out, plotlingauss, AIC = R2_allexac_data)
lingauss_all_plots[[5]]

pdf(paste0(figure_filepath,"/linear_gauss_curves_all.pdf"), width = 16, height = 12, onefile = TRUE)
lingauss_all_plots
dev.off()

# Exact Loc write out optim and named loc data
write_optima_exact(gauss_all_out, "optimum_flowering_time.csv")
write_loc_lookup(gauss_all_out, "exact_loc_lookup.csv")


### SUM STATS on EXACT data ------------------------------------

# run this line if you haven't run the rest of the script
gauss_all_out <- readRDS(paste0(deriv_data_filepath,"/gauss_allexact_mods.RData"))

# make data into big dataframe
exact_in_data <- purrr::map(names(gauss_all_out), ~gauss_all_out[[.x]]$deriv_site %>% 
                           cbind() %>% as.data.frame %>% 
                           mutate(optim_group = .x)
) %>% bind_rows()

# how many named locations, hybrids, years
exact_in_data %>% count(named_location) %>% nrow() # 20 locations
exact_in_data %>% group_by(named_location, Year) %>% n_groups() # 296 site years

exact_in_data %>% distinct(named_location, Year) %>% group_by(named_location) %>% 
  summarise(n = n()) %>%
  summarise(mean = mean(n), sd = sd(n)) # 15 years +/- 11.5 years per location

exact_in_data %>% count(named_location, Year) %>%
  summarise(mean = mean(n), sd = sd(n)) # 43 +/-25(sd) hybrids per site-year 
exact_in_data %>% count(Unif_Name) %>% nrow() # 2039 hybrids
exact_in_data %>% count(Year) %>% nrow() # 46 years




## ALL Sites, COUNTY location -----------------------------------------------
# Create list of sites:
# exclude Minot-NuSun 2000 & Brookings 1998 because multipe trials were planted
# in these sites and years and these are the smaller trials
sites_county <- list_sites(deriv_data, type = 'all_county', exclude_double = TRUE)
  
# AIC subset
R2_county_data <- AIC_R2_data %>% filter(subset == 'all_countyloc')

### Gaussfit
gauss_county_out <- lapply(sites_county, 
                        gaussfit, 
                        data = deriv_data, 
                        AIC = R2_county_data,
                        loc.type = 'named_location')

# Remove elements of list that are NA (only have one site-yr)
gauss_county_out <- Filter(function(a) any(!is.na(a)), gauss_county_out)

# Write all-exact models to an r-file
saveRDS(gauss_county_out, file = paste0(deriv_data_filepath,"/gauss_county_mods.RData"))

## Plotlingauss
lingauss_county_plots <- lapply(gauss_county_out, plotlingauss, AIC = R2_county_data)
lingauss_county_plots[[3]]

pdf(paste0(figure_filepath,"/linear_gauss_curves_county.pdf"), width = 16, height = 12, onefile = TRUE)
lingauss_county_plots
dev.off()

# COUNTY write out optim and named loc data
write_optima_county(gauss_county_out, "optimum_flowering_time_county.csv")
write_loc_lookup(gauss_county_out, "county_lookup.csv")

### SUM STATS on COUNTY data ------------------------------------

# run this line if you haven't run the rest of the script
gauss_county_out <- readRDS(paste0(deriv_data_filepath,"/gauss_county_mods.RData"))

# make data into big dataframe
county_in_data <- purrr::map(names(gauss_county_out), ~gauss_county_out[[.x]]$deriv_site %>% 
                              cbind() %>% as.data.frame %>% 
                              mutate(optim_group = .x)
) %>% bind_rows()

# how many named locations, hybrids, years
county_in_data %>% count(named_location) %>% nrow() # 54 unique locations
county_in_data %>% count(county_state) %>% nrow() # 27 location groups (counties)
county_in_data %>% group_by(county_state, Year) %>% n_groups() # 330 site years

county_in_data %>% distinct(county_state, Year) %>% group_by(county_state) %>% 
  summarise(n = n()) %>%
  summarise(mean = mean(n), median = median(n), sd = sd(n)) # 12.3 years +/- 10.9 years per county

county_in_data %>% count(county_state, Year) %>%
  summarise(mean = mean(n), median = median(n), sd = sd(n), min = min(n), max = max(n)) # 43 +/-25(sd) hybrids per site-year 
county_in_data %>% count(Unif_Name) %>% nrow() # 2131 hybrids
county_in_data %>% count(Year) %>% nrow() # 46 years

# how variable is each hybrid in DPP?
county_in_data %>% group_by(Unif_Name) %>% 
  summarise(n = n(),
            mean = mean(flower_50pct_days_past_planting), 
            range = max(flower_50pct_days_past_planting) - min(flower_50pct_days_past_planting),
            sd = sd(flower_50pct_days_past_planting)) %>% 
  filter(n > 1) %>% 
  reframe(mean_ra = mean(range, na.rm = T), 
          sd_ra = sd(range, na.rm = T),
          CV_ra = sd_ra/mean_ra,
          min_ra = min(range, na.rm = T),
          max_ra = max(range, na.rm = T))



# Single site-year graphs - for presentations ---------------------------------


# Site that has optimum
car2000 <- gauss_county_out$Foster_ND
car2000_plant_date <- car2000$deriv_site$planting_date %>% median() %>% ymd() %>% yday()

car_00 <- ggplot() +
  geom_point(data = car2000$deriv_site %>% 
               filter(Year == 2000),
             aes(x = flower_50_doy - car2000_plant_date,
                 y = yield_lb_acre), alpha = 0.3, size = 2) +
  geom_line(data = car2000$pred_out %>% 
              filter(Year == 2000) %>%
              mutate(y.pred = ifelse(y.pred > 4000 | y.pred < 1500, NA, y.pred),
                     x.pred = ifelse(y.pred > 4000 | y.pred < 1500, NA, x.pred)),
            aes(x = x.pred - car2000_plant_date, 
                y = y.pred), color = "#4EA72E", linewidth = 1.5) +
  labs(title = "Suitable", subtitle = "Carrington ND, 2000", x = "Days to Flowering", y = "Yield") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(color = "#4EA72E", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))
car_00
# ggsave("temp_plots/ec/car_opt.png", width = 3, height = 2.5, units = 'in', dpi = 300)



# Site with early flowering
hays1995 <- gauss_county_out$Ellis_KS
hays1995_plant_date <- hays1995$deriv_site$planting_date %>% median() %>% ymd() %>% yday()

hay_95 <- ggplot() +
  geom_point(data = hays1995$deriv_site %>% 
               filter(Year == 1995),
             aes(x = flower_50_doy - car2000_plant_date,
                 y = yield_lb_acre), alpha = 0.5, size = 2) +
  geom_line(data = hays1995$pred_out %>% 
              filter(Year == 1995) %>%
              mutate(y.pred.lm = ifelse(y.pred.lm > 1300 | y.pred.lm < 600, NA, y.pred.lm),
                     x.pred = ifelse(y.pred > 1300 | y.pred < 600, NA, x.pred)),
            aes(x = x.pred - car2000_plant_date, 
                y = y.pred.lm), color = '#79CDCD', linewidth = 1.5) +
  labs(title = "Mismatch: Early", subtitle = "Hays KS, 1995", x = "Days to Flowering", y = "") +
  lims(x = c(83, 98)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(color = "#79CDCD", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))
hay_95
# ggsave("temp_plots/ec/hays_early.png", width = 3, height = 2.5, units = 'in', dpi = 300)



# Site with late flowering
lang2014 <- gauss_county_out$Cavalier_ND
lang2014_plant_date <- lang2014$deriv_site$planting_date %>% median() %>% ymd() %>% yday()

lan_00 <- ggplot() +
  geom_point(data = lang2014$deriv_site %>% 
               filter(Year == 2014),
             aes(x = flower_50_doy - car2000_plant_date,
                 y = yield_lb_acre), alpha = 0.5, size = 2) +
  geom_line(data = lang2014$pred_out %>% 
              filter(Year == 2014) %>%
              mutate(y.pred.lm = ifelse(y.pred.lm > 3300 | y.pred.lm < 1900, NA, y.pred.lm),
                     x.pred = ifelse(y.pred > 3300 | y.pred < 1900, NA, x.pred)),
            aes(x = x.pred - car2000_plant_date, 
                y = y.pred.lm), color = '#E97132', linewidth = 1.5) +
  labs(title = "Mismatch: Late", subtitle = "Langdon ND, 2014", x = "Days to Flowering", y = "") +
  lims(x = c(63,77)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(color = "#E97132", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))
lan_00
# ggsave("temp_plots/ec/lang_late.png", width = 3, height = 2.5, units = 'in', dpi = 300)



#2019 Williston - no optimum or trend
will2019 <- gauss_county_out$Williams_ND
will2019_plant_date <- will2019$deriv_site$planting_date %>% median() %>% ymd() %>% yday()

wil_19 <- ggplot() +
  geom_point(data = will2019$deriv_site %>% 
               filter(Year == 2019),
             aes(x = flower_50_doy - car2000_plant_date,
                 y = yield_lb_acre), alpha = 0.5, size = 2) +
  labs(title = "No optimum", subtitle = "Williston ND, 2019", x = "Days to Flowering", y = "") +
  lims(x = c(60,80), y = c(700, 1500)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(color = "grey30", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))
wil_19
# ggsave("temp_plots/ec/Will_no_opt.png", width = 3, height = 2.5, units = 'in', dpi = 300)

ggarrange(car_00, hay_95, lan_00, wil_19, nrow = 1, ncol = 4)
ggsave(paste0(figure_filepath,"/county/single_site_yr_examples.png"),
       width = 16, height = 4.5)




# ANALYSES ON OPTIMA ------------------------------------------------------
#-------------------------------------------------------------------------#

## PERFORMANCE-OPT CATEGORY ----------------------------------------------------
## How does performance change when a trial achieves or does not achieve optimum?

# performance data
## set of data used to calculate all the optima
gauss_county_out <- readRDS(paste0(deriv_data_filepath,"/gauss_county_mods.RData"))
subset_data <- purrr::map(names(gauss_county_out), ~gauss_county_out[[.x]]$deriv_site %>% 
                           cbind() %>% as.data.frame %>% 
                           mutate(optim_group = .x)
) %>% bind_rows()

# average performance of each trial
trial_performance_opt <- subset_data %>% 
  group_by(county_state, Year) %>% 
  summarise(trial_yield = mean(yield_lb_acre, na.rm = T)) %>% 
  left_join(subset_data %>%
              distinct(county_state, Year, lat))

# optimum data
optim_county <- read.csv(paste0(deriv_data_filepath,"/optimum_flowering_time_county.csv"))

# join them
trial_performance_opt <- trial_performance_opt %>% 
  left_join(., optim_county %>% 
              select(county_state, Year, mu, mu_025, mu_975, opt_relation_to_data) %>%
              mutate(Year = as.factor(Year))) %>%
  mutate(Year = as.numeric(as.character(Year))) %>%
  mutate(opt_relation_to_data = factor(opt_relation_to_data, 
                                       levels = c('opt_within', 'opt_after', 'opt_before', 'unclear'),
                                       labels = c('Suitable', 'Early', 'Late', 'Unclear')))

sum(is.na(trial_performance_opt$opt_relation_to_data))

# filter to only repeated trials - there are no counties that have only one trial
trial_performance_opt %>% count(county_state) %>% filter(n == 1)

# summary of performance by opt_relation category
trial_performance_opt %>% group_by(opt_relation_to_data) %>% 
  summarise(mean_yield = mean(trial_yield, na.rm = T),
            sd_yield = sd(trial_yield, na.rm = T),
            n = n()) 

# plots by grouping variables
trial_performance_opt %>% 
  group_by(opt_relation_to_data, Year) %>% 
  summarise(mean_yield = mean(trial_yield, na.rm = T)) %>%
  ggplot() +
  geom_point(aes(x = opt_relation_to_data, y = mean_yield)) +
  facet_wrap(~Year, scales = 'fixed')

trial_performance_opt %>% 
  group_by(opt_relation_to_data, county_state) %>% 
  summarise(mean_yield = mean(trial_yield, na.rm = T)) %>%
  ggplot() +
  geom_point(aes(x = opt_relation_to_data, y = mean_yield)) +
  facet_wrap(~county_state, scales = 'fixed')

# model performance
# performance ~ site + opt_relation
mod_perf <- lmer(trial_yield ~ opt_relation_to_data + (1|Year) + (1|county_state) , data = trial_performance_opt)
summary(mod_perf)
emmeans::emmeans(mod_perf, pairwise~opt_relation_to_data, type = 'response', adjust = 'none')

# percent reductions in yield compared to suitable
yield_suitable <- emmeans::emmeans(mod_perf, ~opt_relation_to_data, type = 'response') %>% as.data.frame() %>% 
  filter(opt_relation_to_data == 'Suitable') %>% pull(emmean)

# nicer plot of emmeans result
emmeans::emmeans(mod_perf, ~opt_relation_to_data, type = 'response') %>% as.data.frame() %>% 
  mutate(percent_diff = (emmean - yield_suitable)/yield_suitable) %>% 
  mutate(percent_diff2 = paste0(round(percent_diff*100, 0), "%")) %>%
  mutate(percent_diff2 = ifelse(opt_relation_to_data == "Suitable", NA, percent_diff2)) %>%
  ggplot() +
  geom_jitter(data = trial_performance_opt, aes(x = opt_relation_to_data, y = trial_yield), 
              alpha = 0.2, color = 'grey50', shape = 16, position = position_jitter(width = 0.1)) +
  geom_pointrange(aes(x = opt_relation_to_data, y = emmean, ymin = lower.CL, ymax = upper.CL), 
                  size = .8) +
  geom_hline(aes(yintercept = emmean[1]), linetype = 2, size = .5) +
  # add confidence intervals in text
  # geom_text(aes(x = opt_relation_to_data, y = 100, label = paste0("[",round(lower.CL,0),", ",round(upper.CL,0),"]"))) +
  geom_text(aes(x = opt_relation_to_data, y = 300, label = percent_diff2)) +
  labs(y = "Mean Yield (lb/acre)", x = "Suitability of Flowering", color = 'Latitude') +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

ggsave(paste0(figure_filepath,"/county/yield_optimum_relationship.png"), width = 6, height = 4, dpi = 300)

ggsave(paste0(figure_filepath,"/county/yield_optimum_relationship.pdf"), 
       width = 6, height = 4, dpi = 600)



## ABSOLUTE YIELD CHANGE -----------------------------------------------------------
# How much is yield predicted to decrease when flowering is +/- 7 days from optimum

R2_full_data <- read.csv(paste0(deriv_data_filepath,"/R2_AIC_data.csv"))

R2_allcounty_data <- R2_full_data %>% filter(subset == 'all_countyloc') %>%
  mutate(county = group, .after = group)

# average yield and flowering window of each trial
trial_yld_flow_summary <- subset_data %>% 
  group_by(county_state, Year) %>% 
  summarise(trial_yield = mean(yield_lb_acre, na.rm = T),
            first_flow = min(flower_50_doy),
            last_flow = max(flower_50_doy)) %>% 
  left_join(subset_data %>%
              distinct(county_state, Year, lat))


# optimum data
optim_county <- read.csv(paste0(deriv_data_filepath,"/optimum_flowering_time_county.csv"))
abs_yield_change <- full_join(
  # coefficients
  R2_allcounty_data %>% select(county, Year, mod_for_cens, starts_with("beta"), -beta1_lm_sig),
  optim_county %>% 
    select(county_state, Year, mu, mod_for_cens, opt_relation_to_data)
    # mutate(Year = as.numeric(as.character(Year))) %>%
    # rename(county = county_state),
  ,
  by = c('county' = 'county_state', "Year", 'mod_for_cens')
  ) %>%
  full_join(., trial_yld_flow_summary %>% 
              mutate(Year = as.integer(as.character(Year))), 
            by = c('county' = 'county_state', "Year")
) %>%
  filter(!is.na(mu))

# visualizing the curves, for sanity checks
tst_data <- data.frame(flow_day = seq((214-60),(214+60), by = 1))
tst_beta0 <-abs_yield_change$beta0_gau[1]
tst_beta1 <-abs_yield_change$beta1_gau[1]
tst_beta2 <-abs_yield_change$beta2_gau[1]
tst_lin0 <-abs_yield_change$beta0_lm[1]
tst_lin1 <-abs_yield_change$beta1_lm[1]
tst_data <- tst_data %>%
  mutate(yield = tst_beta0 + tst_beta1 * (flow_day) + tst_beta2 * (flow_day)^2,
         yield_unlog = exp(yield),
         linear = tst_lin0 + tst_lin1 * flow_day,
         linear_unlog = exp(linear))
ggplot(tst_data, aes(x = flow_day, y = yield_unlog)) +
  geom_line()
ggplot(tst_data, aes(x = flow_day, y = linear_unlog)) +
  geom_line()

# calculate predicted change in yield with 7 day shift in flowering
abs_yield_change <- abs_yield_change %>%
  mutate(
    # for quadratic models 
    gaus_delta0 = exp(beta0_gau + beta1_gau * (mu) + beta2_gau * (mu)^2),
    gaus_delta7 = exp(beta0_gau + beta1_gau * (mu + 7) + beta2_gau * (mu + 7)^2),
    gaus_delta_pct = (gaus_delta0 - gaus_delta7) / gaus_delta0,
    # for linear models - positive slope
    lin_delta_pos0 = exp(beta0_lm + beta1_lm * last_flow),
    lin_delta_pos7 = exp(beta0_lm + beta1_lm * (last_flow-7)),
    lin_delta_pos_pct = (lin_delta_pos0 - lin_delta_pos7) / lin_delta_pos0,
    # for linear models - negative slope
    lin_delta_neg0 = exp(beta0_lm + beta1_lm * first_flow),
    lin_delta_neg7 = exp(beta0_lm + beta1_lm * (first_flow+7)),
    lin_delta_neg_pct = (lin_delta_neg0 - lin_delta_neg7) / lin_delta_neg0,
    # new column that combines the previous columns based on whether that 
    # trial has an optimum before, after, within, or not at all
    delta_yield_corrected = as.numeric(ifelse(opt_relation_to_data == 'opt_within', gaus_delta0 - gaus_delta7,
                                   ifelse(opt_relation_to_data == 'opt_after', lin_delta_pos0 - lin_delta_pos7,
                                          ifelse(opt_relation_to_data == 'opt_before', lin_delta_neg0 - lin_delta_neg7,
                                                 ifelse(opt_relation_to_data == 'unclear', 0, "HELP ME"))))),
    delta_yield_pct = as.numeric(ifelse(opt_relation_to_data == 'opt_within', gaus_delta_pct,
                                              ifelse(opt_relation_to_data == 'opt_after', lin_delta_pos_pct,
                                                     ifelse(opt_relation_to_data == 'opt_before', lin_delta_neg_pct,
                                                            ifelse(opt_relation_to_data == 'unclear', 0, "HELP ME"))))),
    delta_yield_pct = delta_yield_pct*100
    )

# bootstrap confidence intervals around the mean
library(boot)

# Bootstrap function to calculate median
median_boot <- function(data, indices) {
  sample_data <- data[indices]
  return(median(sample_data, na.rm = TRUE))
}

# Run bootstrap
set.seed(123)
boot_results <- boot(
  data = abs_yield_change %>% filter(opt_relation_to_data != 'unclear') %>% .$delta_yield_pct,
  statistic = median_boot,
  R = 10000
)

# Get 95% CI using percentile method
boot.ci(boot_results, type = "perc")
boot_ci <- boot.ci(boot_results, type = "perc")$perc[4:5]


# daysoff = 5
# dec_yield_fn <- function(daysoff){
# out <- abs_yield_change %>%
#   mutate(
#     # for quadratic models 
#     gaus_delta7 = exp(beta0_gau + beta1_gau * (mu) + beta2_gau * (mu)^2) -
#       exp(beta0_gau + beta1_gau * (mu + daysoff) + beta2_gau * (mu + daysoff)^2),
#     # for linear models - positive slope
#     lin_delta_pos = exp(beta0_lm + beta1_lm * last_flow) - 
#       exp(beta0_lm + beta1_lm * (last_flow-daysoff)),
#     # for linear models - negative slope
#     lin_delta_neg = exp(beta0_lm + beta1_lm * first_flow) - 
#       exp(beta0_lm + beta1_lm * (first_flow+daysoff)),
#     # new column that combines the previous columns based on whether that 
#     # trial has an optimum before, after, within, or not at all
#     delta_yield_corrected = as.numeric(ifelse(opt_relation_to_data == 'opt_within', gaus_delta7,
#                                               ifelse(opt_relation_to_data == 'opt_after', lin_delta_pos,
#                                                      ifelse(opt_relation_to_data == 'opt_before', lin_delta_neg,
#                                                             ifelse(opt_relation_to_data == 'unclear', 0, "HELP ME"))))),
#     delta_yield_pct = delta_yield_corrected/trial_yield*100
#   )
#     out$delta_yield_pct
# }
# dec_yield_fn(1)

# summary stats - average decrease in yield 
abs_yield_change %>% filter(opt_relation_to_data != 'unclear') %>%
  summarise(mean_abs = mean(delta_yield_corrected, na.rm = T),
            mean_pct = mean(delta_yield_pct, na.rm = T),
            median_pct = median(delta_yield_pct))
median_yield_change <- as.numeric(abs_yield_change %>% filter(opt_relation_to_data != 'unclear') %>%
  summarise(median_pct = median(delta_yield_pct)))
abs_yield_change %>%
  group_by(opt_relation_to_data) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(proportion = n / sum(n))

# abs_yield_change_3day %>%
#   summarise(mean_abs = mean(delta_yield_corrected, na.rm = T),
#             mean_pct = mean(delta_yield_pct, na.rm = T))


yld_dec_plt <- ggplot(abs_yield_change %>% filter(opt_relation_to_data != 'unclear'), 
       aes(x = delta_yield_pct)) +
  geom_histogram(bins = 32, fill = 'royalblue', alpha = 0.6) +
  geom_vline(xintercept = median_yield_change, color = 'black', linetype = 5, size = 0.7) +
  # annotate('text', x = median_yield_change + 4, y = 25, hjust = 0, label = paste0(round(median_yield_change), "% median decrease in yield"), 
  #          color = 'black', size = 4) +
  annotate("rect", xmin = boot_ci[1], xmax = boot_ci[2], ymin = -Inf, ymax = Inf,
           fill = "grey50", alpha = 0.4) +
  labs(x = "Percent decrease in yield with\n7 day shift from optimal flowering", y = "Number of trials") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
yld_dec_plt

ggsave(paste0(figure_filepath,"/county/yld_dec_7_days.png"),
       width = 6.5, height = 5)

# Run code from 01_R2_analyses to get the r2_state_plot 
# in the global environment - this is just for publicaiton figure
ggarrange(r2_state_plot, yld_dec_plt, labels = "AUTO")

ggsave(paste0(figure_filepath,"/county/importance_2panel.pdf"),
       width = 10, height = 4, dpi = 600)

# Checking that change in yield does not depend on group
ggplot(abs_yield_change, aes(y = delta_yield_pct, x = opt_relation_to_data)) +
  geom_boxplot()

lm(delta_yield_pct~ opt_relation_to_data, 
   data = abs_yield_change %>% filter(opt_relation_to_data != 'unclear')) %>%
  emmeans::emmeans(pairwise~opt_relation_to_data)

 

###################################################--
## Questions about optimums ------------------------------------------------
# optims <- read.csv("data_derived/optimum_flowering_time.csv")
# 
#   
# gauss_opts <- optims %>% filter(mod_for_cens == 'gauss' & outside_data == 'in_range')
# gauss_opts$named_location <- as.factor(gauss_opts$named_location)
# gauss_opts$State <- as.factor(gauss_opts$State)
# str(gauss_opts)
# 
# # Does the optimum change over time or lat/lon?
# opt_1 <- lmer(mu ~ planting_doy + latitude + longitude + Year + (1|named_location) , data = gauss_opts)
# summary(opt_1)
# anova(opt_1, type = 3)
# 
# ggarrange(
# ggplot(gauss_opts, aes(x = planting_doy, y = mu, color = named_location)) +
#   geom_point() +
#   geom_smooth(method = 'lm', se = F) +
#   labs(y = "Optimal Flowering DOY") +
#   theme_bw(base_size = 16) +
#   theme(panel.grid = element_blank()),
# 
# ggplot(gauss_opts, aes(x = Year, y = mu, color = named_location)) +
#   geom_point() +
#   geom_smooth(method = 'lm', se = F) +
#   labs(y = "Optimal Flowering DOY") +
#   theme_bw(base_size = 16) +
#   theme(panel.grid = element_blank()),
# 
# ggplot(gauss_opts, aes(x = latitude, y = mu)) +
#   geom_point(aes(color = named_location)) +
#   geom_smooth(method = 'lm', se = F)+
#   labs(y = "Optimal Flowering DOY") +  
#   theme_bw(base_size = 16) +
#   theme(panel.grid = element_blank()),
# 
# ggplot(gauss_opts, aes(x = longitude, y = mu)) +
#   geom_point(aes(color = named_location)) +
#   geom_smooth(method = 'lm', se = F) +
#   labs(y = "Optimal Flowering DOY") +  
#   theme_bw(base_size = 16) +
#   theme(panel.grid = element_blank()),
# nrow = 2, ncol = 2, common.legend = T, legend = 'right'
# )
# 
# # Is optimal flowering time different at different sites?
# opt_2 <- lm(mu ~ planting_doy * named_location , data = gauss_opts)
# summary(opt_2)
# Anova(opt_2, type = 3)
# 
# emmeans(opt_2, ~named_location|planting_doy) %>%
#   as.data.frame() %>% 
#   mutate(opt_date = as.Date(emmean, origin = "2001-01-01"),
#          upper.CL = as.Date(upper.CL, origin = "2001-01-01"),
#          lower.CL = as.Date(lower.CL, origin = "2001-01-01")) %>%
#   left_join(., gauss_opts %>% select(named_location, State) %>% distinct()) %>%
#   ggplot() +
#   geom_point(data = gauss_opts, aes(x = named_location, y = as.Date(mu, origin = "2001-01-01")),
#              alpha = 0.3, position = position_jitter(width = 0.25)) +
#   geom_pointrange(aes(x = named_location, y = opt_date, 
#                       ymin = lower.CL, ymax = upper.CL, color = State), 
#              position = position_dodge(width = 0.5), fatten = 4) +
#   facet_grid(~State, scales = 'free_x', space = 'free') +
#   scale_y_date(date_breaks = "1 week", date_labels = "%b%d") +
#   labs(y = "Optimal Flowering Day") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# 
# emmeans(opt_2, ~named_location|planting_doy, at = list(planting_doy = c(135, 150, 165))) %>%
#   as.data.frame() %>% 
#   mutate(planting_doy = as.factor(planting_doy)) %>%
#   left_join(., gauss_opts %>% select(named_location, State) %>% distinct()) %>%
#   ggplot(aes(x = planting_doy, y = emmean)) +
#   geom_line(aes(group = named_location)) +
#   geom_point(aes(color = State), size = 2) +
#   facet_wrap(~State)+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# # How has optimal flowering time changed over time?
# opt_3 <- lmer(mu ~ planting_doy + Year + (1|named_location) , data = gauss_opts)
# summary(opt_3)
# ggplot(gauss_opts, aes(x = Year, y = mu, color = named_location)) +
#   geom_point() +
#   geom_smooth(method = 'lm', se = F) +
#   labs(y = "Optimal Flowering DOY")
# 
# # What proportion of sites have an optima inside data across years?
# proportion_optima <- optims %>% group_by(Year) %>%
#   summarize(
#     total_sites = n(),
#     # sites_with_opt = sum(outside_data == 'in_range'),
#     prop_optimal = sum(outside_data == 'in_range')/total_sites,
#     prop_early = sum(outside_data == 'greater_than')/total_sites,
#     prop_late = sum(outside_data == 'less_than')/total_sites
#   )
# lm(prop_has_opt ~ Year, data = proportion_optima) %>%
#   summary()
# proportion_optima %>%
#   ggplot(aes(x = Year, y = prop_has_opt)) +
#   geom_point() +
#   geom_smooth(method = lm)
# proportion_optima %>%
#   pivot_longer(cols = starts_with("prop"), 
#                names_to = "flowering_was", 
#                names_prefix = "prop_",
#                values_to = "proportion") %>%
#   ggplot(aes(x = Year, y = proportion, fill = flowering_was)) +
#   geom_col()
#   
#   
# 
# # What proportion of sites have an optima inside data across locations?
# proportion_optima_loc <- optims %>% group_by(location) %>%
#   summarize(
#     total_sites = n(),
#     sites_with_opt = sum(outside_data == 'in_range'),
#     prop_has_opt = sites_with_opt/total_sites
#   ) %>% filter(total_sites > 1) %>%
#   mutate(location = as.factor(location))
# lm(prop_has_opt ~ location, data = proportion_optima_loc) %>%
#   summary()
#   Anova(type = 3)
# proportion_optima_loc %>%
#   ggplot(aes(x = location, y = prop_has_opt)) +
#   geom_point() +
#   geom_smooth(method = lm)


