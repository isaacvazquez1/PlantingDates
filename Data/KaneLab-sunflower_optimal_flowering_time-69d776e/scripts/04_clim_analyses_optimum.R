# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
############    Purpose: What climate parameters
############     explain variation in optimal flowering time?
############ -
############            By: Eliza Clark & Sarah E.
############       Last modified : 10/15/2024
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(brms)
library(tidybayes)
library(rstan)
library(GGally)
library(ggmcmc)
library(lisa)
library(ggpmisc)
library(tinytable)


select <- dplyr::select

deriv_data_filepath <- "derived_data"
figure_filepath <- "figures"


# Should the analysis be done with locations grouped by county or by exact location?
by_county = TRUE

if(by_county == FALSE){
  path2plots = paste0(figure_filepath,"/exact_location/")
}
if(by_county == TRUE){
  path2plots = paste0(figure_filepath,"/county/")
}


# Data ---------------------------------------------------------------------
# Combine optima dataset with the climate dataset

# By county Data
if(by_county == TRUE){
optims <- read.csv(paste0(deriv_data_filepath,"/optimum_flowering_time_county.csv"))
}

# By Exact location data
if(by_county == FALSE){
optims <- read.csv(paste0(deriv_data_filepath,"/optimum_flowering_time.csv"))
}


clim_yearly_out <- read.csv(paste0(deriv_data_filepath,"/yearly_climate_data.csv"))

# get clim_yearly dataset ready
# Fix up Year column
# add group column = either county_state or location
clim_yearly_out_tmp <- clim_yearly_out %>%
  mutate(year = year(as.Date(as.character(year), format = "%Y"))) %>%
  rename(Year = year) %>%
  mutate(group = if(by_county) county_state, .after = State)



# calculate site means and anomalies of climate variables
# group by exact or county location
# if by exact, then group by location, latitude, and longitude
# if by county, group only by county_state
if(by_county == FALSE){
  clim_yearly_out_tmp <- clim_yearly_out_tmp %>%
    mutate(GDD_mod_sitemean = ave(GDD_modave, location, latitude, longitude, FUN = mean),
           GDD_anomaly = GDD_modave - GDD_mod_sitemean,
           GDD_fixed_sitemean = ave(GDD_fixed, location, latitude, longitude, FUN = mean),
           GDD_fixed_anomaly = GDD_fixed - GDD_fixed_sitemean,
           tot_precip_sitemean = ave(tot_precip, location, latitude, longitude, FUN = mean),
           precip_anomaly = tot_precip - tot_precip_sitemean,
           frost_sitemean = ave(frost_day, location, latitude, longitude, FUN = mean),
           frost_anomaly = frost_day - frost_sitemean,
           GDD_future_sitemean = ave(GDD_future, location, latitude, longitude, FUN = mean),
           GDD_future_anomaly = GDD_future - GDD_mod_sitemean)
}

if(by_county == TRUE){
  # clim_yearly data, but only one set for each unique lat lon within each county,
  # so the sitemean is not biased by how many trials were in each part of the county
  county_lat_lons <- clim_yearly_out_tmp %>%
    select(group, latitude, longitude, Year, 
           GDD_modave, GDD_fixed, tot_precip, frost_day, GDD_future) %>%
    distinct(group, latitude, longitude, Year, .keep_all = T)
  
  # sanity check that each lat/long only has one set of climate data
  county_lat_lons %>%
    group_by(group, latitude, longitude) %>%
    summarise(n = n()) %>%
    filter(n > 44)
  
  # Calculate site means, which is average climate for all lat lons within each county.
  # Calculate anomalies, which is the difference between county site mean and 
  # the site value for lat lon each year.
  # GDD_future_anomaly is the difference between the future GDD and 
  # the current site mean GDD
  clim_yearly_out_tmp <- county_lat_lons %>%
    mutate(GDD_mod_sitemean = ave(GDD_modave, group, FUN = mean),
           GDD_anomaly = GDD_modave - GDD_mod_sitemean,
           GDD_fixed_sitemean = ave(GDD_fixed, group, FUN = mean),
           GDD_fixed_anomaly = GDD_fixed - GDD_fixed_sitemean,
           tot_precip_sitemean = ave(tot_precip, group, FUN = mean),
           precip_anomaly = tot_precip - tot_precip_sitemean,
           frost_sitemean = ave(frost_day, group, FUN = mean),
           frost_anomaly = frost_day - frost_sitemean,
           GDD_future_sitemean = ave(GDD_future, group, FUN = mean),
           GDD_future_anomaly = GDD_future - GDD_mod_sitemean)
}



# combine climate and optima datasets

if(by_county == TRUE){
  optim_clim <- left_join(optims %>% mutate(group = county_state), 
                          clim_yearly_out_tmp, join_by(group, Year, latitude, longitude))
}

if(by_county == FALSE){
  optim_clim <- left_join(optims, 
                          clim_yearly_out_tmp, join_by(county_state, location, latitude, longitude, State, Year)) %>%
    select(- X.y, -n, -X.x) %>%
    relocate(county_state, .after = State) %>%
    mutate(group = named_location, .after = county_state)
}

# check for na's - all NAs should be before 1980 (before climate dataset starts)
colSums(is.na(optim_clim))
optim_clim %>% filter(is.na(GDD_mod_sitemean)) %>% select(location, named_location, Year)
optim_clim %>% filter(is.na(y3) & mod_for_cens != 'exclude') %>% 
  select(location, named_location, Year, y2, y3)

# check for sites where y2 >= y2
optim_clim %>% filter(y2 >= y3) %>% select(location, county_state, named_location, Year, y2, y3)


# scale parameters
# Exclude some site-years
# !! Need to go back and check out why this point is y2=y3
optim_clim_sc <- optim_clim %>% 
  filter(mod_for_cens != 'exclude' & # remove site-years without any optimum
           !is.na(y3) & # one site-year that doesn't have y3 due to no harvest date recorded
           # y2 < y3,  # one site-yr where y2=y3
         planting_doy > 75, # remove super early texas outlier in planting
         !(is.na(GDD_mod_sitemean))) %>% # remove site-years without climate data
  select(
    group, named_location, location, county_state, Year, latitude, longitude,
    planting_doy, mu, sd, y2, y3, cens,
    y_na, outside_data, mod_for_cens, 
    opt_relation_to_data, GDD_mod_sitemean:GDD_future_anomaly
  ) %>% 
  mutate(across(GDD_mod_sitemean:GDD_future_anomaly, ~ as.numeric(scale(.x)), .names = "{col}_sc")) %>%
  filter(!is.na(GDD_anomaly_sc))


# write out optim_clim_sc (filtered version with scaling) 
# and optim clim (no filtering) as a csv
if(by_county == FALSE){
  write.csv(optim_clim_sc, paste0(deriv_data_filepath,"/optim_clim_sc.csv"), row.names = FALSE)
  write.csv(optim_clim, paste0(deriv_data_filepath,"/optim_clim.csv"), row.names = FALSE)
}

if(by_county == TRUE){
  write.csv(optim_clim_sc, paste0(deriv_data_filepath,"/optim_clim_sc_county.csv"), row.names = FALSE)
  write.csv(optim_clim, paste0(deriv_data_filepath,"/optim_clim_county.csv"), row.names = FALSE)
}




## Planting DOY & GDD site mean CORRELATION -------------------------------
# Is planting later in hotter years?
fixed_vs_relative_gdd <- optim_clim_sc %>%
  ggplot() +
  geom_point(aes(x = GDD_fixed_sitemean, y = planting_doy, color = latitude)) +
  geom_smooth(aes(x = GDD_fixed_sitemean, y = planting_doy), method = 'lm', 
              color = 'dodgerblue3', alpha = 0.3) +
  geom_point(aes(x = GDD_mod_sitemean, y = planting_doy, color = latitude)) +
  geom_smooth(aes(x = GDD_mod_sitemean, y = planting_doy), color = 'brown3', 
              linetype = 1, method = 'lm', alpha = 0.3) +
  annotate(geom = 'text', x = 1400, y = 190, label = "Fixed window\nGDD = doy 125 to 250", color = 'dodgerblue4') +
  annotate(geom = 'text', x = 800, y = 190, label = "Relative to planting\nGDD = planting doy + 70", color = 'brown3') +
  stat_regline_equation(aes(x = GDD_fixed_sitemean, y = planting_doy),
                        label.x = 1325, label.y = 185, color = 'dodgerblue4', hjust = 0) +
  stat_regline_equation(aes(x = GDD_mod_sitemean, y = planting_doy),
                        label.x = 725, label.y = 185, color = 'brown3', hjust = 0) +
  # stat_poly_eq(aes(x = GDD_mod_sitemean, y = planting_doy), 
  #              use_label(c("eq", 'p')), inherit.aes = F) +
  labs(x = "Growing degree days", y = 'Planting day of year') +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
fixed_vs_relative_gdd

## Distribution of optimal dates -----------------------------------------------
  optim_clim_sc %>% 
  mutate(combined_mu = ifelse(opt_relation_to_data == 'opt_within', mu,
                              ifelse(opt_relation_to_data == 'opt_after', y2, 
                                     ifelse(opt_relation_to_data == 'opt_before', y3, NA))),
         opt_relation_to_data = factor(opt_relation_to_data, 
                                       levels = c("opt_before", "opt_after", "opt_within"),
                                       labels = c("Optimal timing\n   less than value", "Optimal timing\n   greater than value", "Optimal timing")),
         combined_dpp = combined_mu - planting_doy,
         opt_censored = ifelse(opt_relation_to_data == 'opt_within', 'Uncensored', "Censored")) %>%
  ggplot(aes(x = combined_dpp)) +
  # geom_freqpoly(aes(x = combined_dpp, color = opt_relation_to_data), bins = 30) +
  # geom_density(aes(x = combined_dpp, y = after_stat(count), fill = opt_relation_to_data), 
  #              alpha = 0.5, bw = 5) +
  # geom_area(aes(y = after_stat(count), fill = opt_relation_to_data), alpha = .75, stat = "bin") +
  geom_histogram(aes(x = combined_dpp, fill = opt_relation_to_data), position = 'stack', 
                 alpha = 0.75, color = 'white', linewidth = 0) +
  # geom_dotplot(aes(x = combined_dpp, fill = opt_relation_to_data), method = "histodot", 
  #              binwidth = 2, stackgroups = T) +
  # scale_fill_manual(values = c('grey50', 'forestgreen')) +
  scale_fill_manual(values = c('orange','darkslategray3','forestgreen',   'grey50')) +
  # scale_fill_brewer(palette = "Dark2", direction = -1) +
  guides(fill = guide_legend(keywidth = 1.5, label.theme = element_text(size = 10))) +
  labs(x = "Days past planting", y = "Count of trials", fill = "") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.75, 0.9),
        legend.background = element_blank())


if(by_county == TRUE){
  ggsave(paste0(path2plots, "dpp_distribution.png"), width = 5, height = 5)
}

# range of intervals
optim_clim_sc %>% 
  mutate(combined_mu = ifelse(opt_relation_to_data == 'opt_within', mu,
                              ifelse(opt_relation_to_data == 'opt_after', y2, 
                                     ifelse(opt_relation_to_data == 'opt_before', y3, NA))),
         combined_dpp = combined_mu - planting_doy,
         opt_censored = ifelse(opt_relation_to_data == 'opt_within', 'Uncensored', "Censored")) %>%
  summarise(min(combined_dpp), max(combined_dpp),
            min(y2-planting_doy), max (y3-planting_doy))

# CV of optimal DPP
optim_clim_sc %>% 
  mutate(combined_mu = ifelse(opt_relation_to_data == 'opt_within', mu,
                              ifelse(opt_relation_to_data == 'opt_after', y2, 
                                     ifelse(opt_relation_to_data == 'opt_before', y3, NA))),
         combined_dpp = combined_mu - planting_doy,
         opt_censored = ifelse(opt_relation_to_data == 'opt_within', 'Uncensored', "Censored")) %>%
  summarise(
    mean = mean(combined_dpp),
    sd = sd(combined_dpp),
    CV = sd(combined_dpp)/mean(combined_dpp))
  
# change in optimal flowering time over time for each county
optim_clim_sc %>% 
  mutate(combined_mu = ifelse(opt_relation_to_data == 'opt_within', mu,
                              ifelse(opt_relation_to_data == 'opt_after', y2, 
                                     ifelse(opt_relation_to_data == 'opt_before', y3, NA))),
         opt_relation_to_data = factor(opt_relation_to_data, 
                                       levels = c("opt_before", "opt_after", "opt_within"),
                                       labels = c("Optimal timing\n   less than value", "Optimal timing\n   greater than value", "Optimal timing")),
         combined_dpp = combined_mu - planting_doy,
         opt_censored = ifelse(opt_relation_to_data == 'opt_within', 'Uncensored', "Censored")) %>%
  ggplot(aes(x = Year, y = combined_dpp)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = county_state), method = 'lm', , alpha = 0.5) +
  stat_poly_eq(aes(color = county_state, label = ..p.value.label..)) +
  facet_wrap(~county_state) 
  


optim_clim_sc %>%
  ggplot() +
  geom_point(aes(x = GDD_fixed_anomaly, y = planting_doy, color = as.factor(Year))) +
  geom_smooth(aes(x = GDD_fixed_anomaly, y = planting_doy), method = 'lm') +
  facet_wrap(~group, scales = 'free') 


  
# custom trace plot geom
fk <- lisa_palette("FridaKahlo", n = 31, type = "continuous")

geom_trace <- function(subtitle = NULL,
                       xlab = "iteration",
                       xbreaks = 0:4 * 500) {
  list(
    annotate(
      geom = "rect",
      xmin = 0, xmax = 1000, ymin = -Inf, ymax = Inf,
      fill = fk[16], alpha = 1 / 2, linewidth = 0
    ),
    geom_line(linewidth = 1 / 3),
    scale_color_manual(values = fk[c(3, 8, 27, 31)]),
    scale_x_continuous(xlab, breaks = xbreaks, expand = c(0, 0)),
    labs(subtitle = subtitle),
    theme(panel.grid = element_blank())
  )
}

# DOY MODEL ---------------------------------------------------------------------


## Create stan data ---------------------------------------------------------
# original model running everything as censored
# predict censored optimal flowering day with degree days

# make some columns into factors
optim_clim_sc <- optim_clim_sc %>%
  mutate(group = as.factor(group), 
         cens = as.factor(cens),
         outside_data = as.factor(outside_data))

# use brms to format the data as per a censored model
sdata1 <- make_standata(
  y2 | cens(cens, y3) ~ GDD_mod_sitemean_sc + GDD_anomaly_sc +
    tot_precip_sitemean_sc + precip_anomaly_sc +
    # frost_sitemean_sc +
    frost_anomaly_sc +
    planting_doy + (1 | group),
  data = optim_clim_sc,
  family = gaussian
)

# add in the estimated peak (mu) and se on peak (sd) from
# the gaussian models

#!! we shoudl check this but we think this is right
# that mu and sd are the right params to pick out of the flowering curves


sdata1$sigma_meas <- optim_clim_sc$sd
sdata1$measured_y <- optim_clim_sc$mu


# pass stan a new variable that is sigma_meas but if that doesn't exist
# just fill with a very big number. If you pass the data in with NAs I think
# they will be dropped becauase Stan hates NAs
# these never enter the likelihood so it's fine
# !! but we should triple check at some point these rows don't have
# crazy fitted values, which they would if these values were entering the
# likelihood.
# or check by making a different crazy value and ensuring you get the
# same results
sdata1$se <- ifelse(is.na(sdata1$sigma_meas), 9999, sdata1$sigma_meas)
sdata1$measured_y <- ifelse(is.na(sdata1$sigma_meas), -9999, sdata1$measured_y)

# create a new indicator value,
# 2 if interval censored, 3 if measurement error
# these could be any value, I kept 2 for interval censored in keeping
# with how brms codes it

## !! this should be made off 'in range' instead for 3
# if data is "in range" and it has an estimate of sd, then cens_2 should be 3 (measurement error),
# if out of range or does not have a sd measured, then 2 (interval censored)
sdata1$cens_2 <- ifelse(optim_clim_sc$outside_data == "in_range" & 
                          !is.na(optim_clim_sc$sd), 3, 2)


# remove the original censored variable that the make_standata function
# put in. We won't use this anymore.
sdata1$cens <- NULL

sdata1$sigma_meas <- NULL # EC: this contains NAs and isn't used later

## STAN Model --------------------------------------------------------------
# define the bake-your-own-frankenlikelihood model
# this is a copy of the brms model with the likelihood modified
# I added SCE on any lines that I added/modified from the censored model
# so you can see what I did
stanmod_2 <- "// generated with brms 2.19.0 then modified
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable for interval censoring (SCE: actually this is left)
  int<lower=2,upper=3> cens_2[N];  // SCE: indicates censoring vs meas error, 2 is cens, 3 is meas error
  vector[N] measured_y;  // SCE estimate for non-interval censored
  vector<lower=0>[N] se;  // SCE known error on optimal estimate for non-interval censored
  vector[N] rcens;  // right censor points for interval censoring
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  vector<lower=0>[N] se2 = square(se);
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
}
// priors block from interval censored model
// this was based on the mean_y and sd_y per brms specs
// We can also just use the default student t. I am
// not entirely sure why the defaults differ between the two
// !! in either case would be good to plot posteriors on prior
// to make sure they are sufficiently uniformative for the data range
//
//transformed parameters {
  //  vector[N_1] r_1_1;  // actual group-level effects
  //  real lprior = 0;  // prior contributions to the log posterior
  //  r_1_1 = (sd_1[1] * (z_1[1]));
  //  lprior += normal_lpdf(Intercept | 221.3135, 30.02485 * 100);
  //  lprior += normal_lpdf(sigma | 0, 30.02485)
  //    - 1 * normal_lccdf(0 | 0, 30.02485);
  //  lprior += student_t_lpdf(sd_1 | 3, 0, 9.9)
  //    - 1 * student_t_lccdf(0 | 3, 0, 9.9);
  //}

transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  r_1_1 = (sd_1[1] * (z_1[1]));
  lprior += student_t_lpdf(Intercept | 3, 219.3, 9.2);
  lprior += student_t_lpdf(sigma | 3, 0, 9.2)
  - 1 * student_t_lccdf(0 | 3, 0, 9.2);
  lprior += student_t_lpdf(sd_1 | 3, 0, 9.2)
  - 1 * student_t_lccdf(0 | 3, 0, 9.2);
}


model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b;
    for (n in 1:N) {
      // special treatment of censored data
      if (cens_2[n] == 2) {
        target += log_diff_exp(
          normal_lcdf(rcens[n] | mu[n], sigma),
          normal_lcdf(Y[n] | mu[n], sigma)
        );
        // SCE added measurement error lines and removed options for l/r censored
      } else if (cens_2[n] == 3) {
      // error is going to be the sum of the error on the measurement and the error on the estimate
      // this is taken from a brms se model
        target += normal_lpdf(Y[n] | mu[n], sqrt(square(sigma) + se2[n]));
      }
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // real mu_estimate = Intercept + Xc * b;
}
"

## Run model -----------------------------------------------------------------
mean_y <- mean(optim_clim$y_na, na.rm = T)

inits <- list(Intercept = mean_y)
inits_list <- list(inits, inits, inits, inits)

# actual model needs to be run with longer burnin/warmup
fit1 <- stan(
  model_code = stanmod_2, # Stan program
  data = sdata1, # named list of data
  chains = 4, # number of Markov chains
  warmup = 500, # number of warmup iterations per chain
  iter = 4500, # total number of iterations per chain
  cores = 4, # number of cores (could use one per chain)
  # refresh = 0,             # no progress shown
  init = inits_list, # init for intercept needs to be specified
  control = list(max_treedepth = 10)
)

# names of parameters
ggs(fit1, burnin = TRUE) %>%
  distinct(Parameter) %>%
  View()

# trace plot !! doesn't really work right now
ggs(fit1, burnin = TRUE) %>%
  filter(Parameter == "b_Intercept") %>%
  mutate(
    chain = factor(Chain),
    intercept = value
  ) %>%
  ggplot(aes(x = Iteration, y = intercept, color = chain)) +
  geom_trace(subtitle = "fit1 (default settings)") +
  scale_y_continuous(breaks = c(0, 650, 1300), limits = c(NA, 1430))

## Prior only model  --------------------------------------------------------
# it's overkill to run this many iterations
# but this is a good shortcut to look at your priors
sdata1_prior <- sdata1
sdata1_prior$prior_only <- 1
fit1_prior <- stan(
  model_code = stanmod_2, # Stan program
  data = sdata1_prior, # named list of data
  chains = 4, # number of Markov chains
  warmup = 50, # number of warmup iterations per chain
  iter = 1000, # total number of iterations per chain
  cores = 4, # number of cores (could use one per chain)
  # refresh = 0,             # no progress shown
  # init = inits_list, #seemed ok without specifying inits
  control = list(max_treedepth = 10),
  init = inits_list
)

# you can plot posteriors on priors like this
# if you want to make sure whatever is in there is not crazy.
# this is for fixefs only, the blue is so tiny you can't see it with
# the prior turned on which indicates the prior is way uninformative
fit1 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_mod_sitemean_sc",
        param_order == 2 ~ "GDD_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  ggplot(aes(x = param_name, y = b)) +
  stat_halfeye(
    color = "blue", fill = "blue", alpha = 0.2,
    side = "top"
  ) +
  stat_halfeye(
    color = "red", fill = "red", alpha = 0.2,
    side = "bottom",
    data = fit1_prior %>%
      spread_draws(b[param_order]) %>%
      mutate(
        param_name =
          case_when(
            param_order == 1 ~ "GDD_mod_sitemean_sc",
            param_order == 2 ~ "GDD_anomaly_sc",
            param_order == 3 ~ "tot_precip_sitemean_sc",
            param_order == 4 ~ "precip_anomaly_sc",
            param_order == 5 ~ "frost_anomaly_sc",
            param_order == 6 ~ "planting_doy",
            TRUE ~ NA_character_
          )
      )
  )

## Results -----------------------------------------------------------------
# sanity check the x matrix is column order by variable name
# as provided in the model formula
# head (optim_clim_sc$tot_precip_sitemean_sc)
# head (sdata1$X[,4])
print(summary(fit1)$summary)

# can look at summary of fixed effects only
summary(fit1, pars = c("Intercept", "b"), probs = c(0.05, 0.95))$summary

# tidybayes way to look at more or less the same thing
# I like that these color the ones where the CI is different from 0
# I did spot check that the data has the covariates entered in the
# order they are entered in the model formula statement.
fit1 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_mod_sitemean_sc",
        param_order == 2 ~ "GDD_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  group_by(param_name) %>%
  median_qi() %>%
  print()


summarise_draws(tidy_draws(fit1))

## Make predictions -------------------------------------------------------


# get min, mean, and max of all variables, to see the interval on which to predict
ranges <- data.frame(GDD_site = c(range(optim_clim_sc$GDD_mod_sitemean_sc), mean(optim_clim_sc$GDD_mod_sitemean_sc)),
                     GDD_anom = c(range(optim_clim_sc$GDD_anomaly_sc), mean(optim_clim_sc$GDD_anomaly_sc)),
                     precip_site = c(range(optim_clim_sc$tot_precip_sitemean_sc), mean(optim_clim_sc$tot_precip_sitemean_sc)),
                     precip_anom = c(range(optim_clim_sc$precip_anomaly_sc), mean(optim_clim_sc$precip_anomaly_sc)),
                     frost_anom = c(range(optim_clim_sc$frost_anomaly_sc), mean(optim_clim_sc$frost_anomaly_sc)),
                     planting_doy = c(range(optim_clim_sc$planting_doy), mean(optim_clim_sc$planting_doy))
)
rownames(ranges) <- c("min", "max", 'mean')
ranges <- t(ranges)

# how many predictions per factor to make
reps <- 25

# functions to create sequences and reps for new data to predict on
myseq <- function(row){seq(ranges[row,1], ranges[row,2], length.out = reps)}
myrep <- function(row){rep(ranges[row,3], reps)}

# create new data on which to predict
# this is marginal predictions - one variable change at a time, the rest are at their mean
new_data <- data.frame(GDD_mod_sitemean_sc = c(myseq(1), rep(myrep(1), 5)),
           GDD_anomaly_sc = c(myrep(2), myseq(2), rep(myrep(2), 4)),
           tot_precip_sitemean_sc = c(rep(myrep(3), 2), myseq(3), rep(myrep(3), 3)),
           precip_anomaly_sc = c(rep(myrep(4), 3), myseq(4), rep(myrep(4), 2)),
           frost_anomaly_sc = c(rep(myrep(5), 4), myseq(5), rep(myrep(5), 1)),
           planting_doy = c(rep(myrep(6), 5), myseq(6)),
           group = rep(c('GDD_mod_sitemean_sc', 'GDD_anomaly_sc', 'tot_precip_sitemean_sc', 'precip_anomaly_sc', 'frost_anomaly_sc', 'planting_doy'), each = 25)
           )

# Extract draws of betas and intercept
betas <- fit1 %>%
  tidy_draws() %>%
  select(starts_with("b")) 
colnames(betas) <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b_Intercept')

# For each row of new data, multiply it out by the betas, 
# then take the median_qi and add that summary to the table
pred_y = NULL
for(i in 1:nrow(new_data)){
  pred_i <- median_qi(betas$b_Intercept + 
                        betas$b1 * new_data$GDD_mod_sitemean_sc[i] +
                        betas$b2 * new_data$GDD_anomaly_sc[i] +
                        betas$b3 * new_data$tot_precip_sitemean_sc[i] +
                        betas$b4 * new_data$precip_anomaly_sc[i] +
                        betas$b5 * new_data$frost_anomaly_sc[i] +
                        betas$b6 * new_data$planting_doy[i])
  pred_y <- rbind(pred_y, pred_i)
}

# combine prediction summaries to the new data
new_data1 <- cbind(new_data, pred_y) 

# for plotting, we need the predictor that varies (x) and the prediction (y)
# and we don't need the other variables that are set to their means
# put these x's and y's into elements of a list 
make_pred_list <- function(predictor, new.dat){
  new.dat %>% 
    filter(group == predictor) %>% 
    select({{predictor}}, 7:10) %>%
    rename(x = {{predictor}}) %>%
    relocate(group, .after = ymax)
}
predictions_list <- as.list(unique(new_data1$group))
predictions_list <- lapply(predictions_list, make_pred_list, new.dat = new_data1)

# get means and standard deviations of the unscaled climate 
# predictors for unscaling
ranges_unsc <- data.frame(GDD_mod_sitemean_sc = c(sd(optim_clim_sc$GDD_mod_sitemean), mean(optim_clim_sc$GDD_mod_sitemean)),
                          GDD_anomaly_sc = c(sd(optim_clim_sc$GDD_anomaly), mean(optim_clim_sc$GDD_anomaly)),
                          tot_precip_sitemean_sc = c(sd(optim_clim_sc$tot_precip_sitemean), mean(optim_clim_sc$tot_precip_sitemean)),
                          precip_anomaly_sc = c(sd(optim_clim_sc$precip_anomaly), mean(optim_clim_sc$precip_anomaly)),
                          frost_anomaly_sc = c(sd(optim_clim_sc$frost_anomaly), mean(optim_clim_sc$frost_anomaly)),
                          planting_doy = c(1, 0)
)
rownames(ranges_unsc) <- c("sd", 'mean')
ranges_unsc <- t(ranges_unsc)

# little function to add a column of the unscaled predictors 
# to the predictions list
unscale_predictions <- function(predictions){
  group_i <- unique(predictions$group)
  predictions <- predictions %>% mutate(x_unsc = 
                           predictions$x * 
                             ranges_unsc[which(rownames(ranges_unsc)==group_i),1] +
                             ranges_unsc[which(rownames(ranges_unsc)==group_i),2],
                           .after = x)
}
predictions_list_unsc <- lapply(predictions_list, unscale_predictions)




## Plots of trends ---------------------------------------------------------

# get a table of the coefficents to put on the plot
coefs <- fit1 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_mod_sitemean_sc",
        param_order == 2 ~ "GDD_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  group_by(param_name) %>%
  median_qi() %>% arrange(param_order)

coefs_unsc <- coefs %>% mutate(b_unsc = b / ranges_unsc[,1],
                 b.lower_unsc = b.lower / ranges_unsc[,1],
                 b.upper_unsc = b.upper / ranges_unsc[,1])

# make a table of x-axis labels
labels <- data.frame(predictor = gsub('_sc', '', coefs_unsc$param_name),
                     label = c("Degree days site mean (°C)",
                               "Degree days anomaly (°C)",
                               "Precipitation site mean (mm)", 
                               "Precipitation anomaly (mm)",
                               "First frost anomaly",
                               "Planting day of year")
)
# function that will make each plot, with trend line, 
# points and intervals of data, and coefficients
# dat = predictions_list_unsc[[3]]
plot_trends2 <- function(dat){
  predictor <- gsub('_sc', '', as.character(unique(dat$group)))
  coef_i <- coefs_unsc %>%
    filter(grepl(predictor, param_name)) %>%
    select(b_unsc, b.lower_unsc, b.upper_unsc) %>%
    round(digits = 3)
  
  ggplot() +
    # vertical lines for censoring intervals
    geom_linerange(
      data = optim_clim_sc %>% filter(outside_data != 'in_range' | is.na(sd)),
      aes(x = .data[[predictor]], 
          ymax = as.Date("2019-12-31") + y3, 
          ymin = as.Date("2019-12-31") + y2, 
          color = latitude),
      alpha = 0.5
    ) +
    # points and lines for measurement error points
    geom_pointrange(
      data = optim_clim_sc %>% filter(outside_data == 'in_range' & !is.na(sd)),
      aes(x = .data[[predictor]], 
          y = as.Date("2019-12-31") + mu, 
          ymax = as.Date("2019-12-31") + mu+sd, 
          ymin = as.Date("2019-12-31") + mu-sd, 
          color = latitude),
      alpha = 0.5
    ) +
    # error of prediction from model
    geom_ribbon(
      data = dat,
      aes(x = x_unsc, ymin = as.Date("2019-12-31") + ymin, ymax = as.Date("2019-12-31") + ymax),
      fill = "grey75", alpha = 0.75
    ) +
    # trend line of prediction from model
    geom_line(
      data = dat,
      aes(x = x_unsc, y = as.Date("2019-12-31") + y), 
      linetype = as.numeric(ifelse(coef_i[2] * coef_i[3] > 0, 1, 2)), 
      linewidth = 0.8
    ) +
    # coefficient annotation
    # annotate("text",
    #          x = -Inf, y = as.Date("2019-12-31") + Inf,
    #          label = paste0("Slope = ", coef_i[1], "\nCI: ", coef_i[2], ", ", coef_i[3]),
    #          color = ifelse(coef_i[2] * coef_i[3] > 0, "red", "black"), size = 5, hjust = 0, vjust = 1
    # ) +
    # formatting
    scale_y_date(date_breaks = "2 weeks", date_labels = "%b %d") +
    labs(
      x = labels[which(labels$predictor == predictor),2],
      y = "Optimal flowering date",
      color = "Latitude"
    ) +
    # guides(color = "none") +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank())
}

# run the function
regression_plots <- lapply(predictions_list_unsc, plot_trends2)
# reorder regression_plots to be 1, 2, 6, 3, 4, 5
regression_plots <- regression_plots[c(1, 2, 6, 3, 4, 5)]

# arrange resulting plots
ggarrange(plotlist = regression_plots, common.legend = T, legend = 'right',
          labels = "AUTO", ncol = 3, nrow = 2)

ggsave(paste0(path2plots,"DOY_opt_clim.png"),
       width = 13, height = 8)



## check out why top left panel looks weird
regression_plots[[1]]

# The planting doy is later at sites with high GDD, 
# and there is a positive trend with opt flow time and planting doy,
# which is positive. So, that's mostly the trend that we see here, but
# the prediction line in at an average planting day.
ggplot() +
  # vertical lines for censoring intervals
  geom_linerange(
    data = optim_clim_sc %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = GDD_mod_sitemean_sc, ymax = y3, ymin = y2, color = planting_doy),
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = optim_clim_sc %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = GDD_mod_sitemean_sc, y = mu, ymax = mu+sd, ymin = mu-sd, color = planting_doy),
    alpha = 0.5
  ) 

# We don't have this same relationship with 
# precipitation site means, for example
ggplot() +
  # vertical lines for censoring intervals
  geom_linerange(
    data = optim_clim_sc %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = tot_precip_sitemean_sc, ymax = y3, ymin = y2, color = planting_doy),
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = optim_clim_sc %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = tot_precip_sitemean_sc, y = mu, ymax = mu+sd, ymin = mu-sd, color = planting_doy),
    alpha = 0.5
  )


## Heatmap of temp x planting doy ----------------------------------------------

# Quantiles of GDD anomaly (colder year, middle year, hotter year)
gdd_anom_quants <- quantile(optim_clim_sc$GDD_anomaly_sc, c(.05, .50, .95)) %>% as.numeric()
# quantile(optim_clim_sc$GDD_future_anomaly_sc, c(.25, .50, .75)) %>% as.numeric()

# function that will generate sequence with a buffer around the min and max
# row = 1
myseq2 <- function(row){seq(ranges[row,1] - 2, ranges[row,2] + 3, length.out = reps)}


# effect of GDD site mean and planting day at three levels of GDD_anomaly
new_data2 <- rbind(
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants[1], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants[2], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants[3], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  )
)

# predict y points
pred_y2 = NULL
for(i in 1:nrow(new_data2)){
  pred_i <- median_qi(betas$b_Intercept + 
                        betas$b1 * new_data2$GDD_mod_sitemean_sc[i] +
                        betas$b2 * new_data2$GDD_anomaly_sc[i] +
                        betas$b3 * new_data2$tot_precip_sitemean_sc[i] +
                        betas$b4 * new_data2$precip_anomaly_sc[i] +
                        betas$b5 * new_data2$frost_anomaly_sc[i] +
                        betas$b6 * new_data2$planting_doy[i])
  pred_y2 <- rbind(pred_y2, pred_i)
}

# combine prediction summaries to the new data
new_data2 <- cbind(new_data2, pred_y2) 

new_data2 <- new_data2 %>% mutate(GDD_mod_sitemean = 
                                        GDD_mod_sitemean_sc * 
                                        ranges_unsc[1,1] +
                                        ranges_unsc[1,2],
                                  GDD_anomaly = 
                                    GDD_anomaly_sc * 
                                    ranges_unsc[2,1] +
                                    ranges_unsc[2,2]
                                  )

new_data2$GDD_anomaly_sc <- factor(new_data2$GDD_anomaly_sc, 
                  labels = c("Cold year\n(5th percentile)", 
                             "Average year\n(50th percentile)", 
                             "Warm year\n(95th percentile)"))

# add columns that say if the planting day/site GDD is actually realistic
realistic_ranges <- optim_clim %>% filter(planting_doy > 90, !is.na(GDD_modave)) %>% 
  select(Year, GDD_modave, planting_doy) %>%
  mutate(GDD_modave_1 = GDD_modave - 35,
         GDD_modave_2 = GDD_modave + 35,
         planting_doy1 = planting_doy - 4,
         planting_doy2 = planting_doy + 4
  )
# range_GDD <- range(optim_clim %>% filter(planting_doy > 90, Year > 1980) %>% select(GDD_modave), na.rm = T)
# range_planting <- range(optim_clim %>% filter(planting_doy > 90, Year > 1980) %>% select(planting_doy), na.rm = T)

new_data2 <- new_data2 %>% mutate(realistic = 
                                    ifelse(apply(cbind(new_data2$GDD_mod_sitemean,new_data2$planting_doy), 1, function (p){
                                      any(realistic_ranges$GDD_modave_1 <= p[1] & 
                                            realistic_ranges$GDD_modave_2 >= p[1] &
                                            realistic_ranges$planting_doy1 <= p[2] & 
                                            realistic_ranges$planting_doy2 >= p[2])
                                    }
                                    ) ,
                                    'yes', 'no'),
                                  y_realistic = ifelse(realistic == 'yes', y, NA))

new_data2 <- new_data2 %>%
  filter(GDD_mod_sitemean > 480 & GDD_mod_sitemean < 1280 & planting_doy > 117 & planting_doy < 198)

ggplot() +
  # error of prediction from model
  geom_ribbon(
    data = new_data2,
    aes(x = GDD_mod_sitemean, ymin = ymin, ymax = ymax, group = as.factor(planting_doy)),
    fill = "grey75", alpha = 0.75
  ) +
  # trend line of prediction from model
  geom_line(
    data = new_data2,
    aes(x = GDD_mod_sitemean, y = y, color = as.factor(planting_doy)), linewidth = 0.8
  ) +
  facet_wrap(~GDD_anomaly_sc) +
  # formatting
  labs(
    y = "Optimal Flowering Day of Year",
    color = "Planting day"
  ) +
  # guides(color = "none") +
  theme_bw(base_size = 16)


# heat map of optimal flowering day of year
ggplot(data = new_data2, aes(x = GDD_mod_sitemean, 
                             y = as.Date("2019-12-31") + planting_doy, 
                             fill = as.Date("2019-12-31") + y_realistic)) +
  geom_tile(linewidth = 0.1) +
  facet_wrap(~GDD_anomaly_sc) +
  # geom_hex(data = optim_clim %>% filter(planting_doy > 90, Year > 1980),
  #          aes(x = GDD_modave, 
  #              y = as.Date("2019-12-31") + planting_doy),
  #          bins = 20, color = 'black', fill = 'white', alpha = 0,
  #          linewidth = 0.9) +
  annotate(geom = 'rect', xmin = 879, xmax = 916, 
           ymin = as.Date("2019-12-31") + 152.5, ymax = as.Date("2019-12-31") + 156,
           alpha = 0, color = 'black', linewidth = 1) +
  labs(fill = "Optimal flowering \nday of year",
       y = "Planting day of year",
       x = "Average site degree days") +
  scale_y_date(date_breaks = '2 weeks', date_labels = "%b %d", expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  # viridis::scale_fill_viridis(trans = "date", na.value = 'grey75') +
  # scale_fill_distiller(palette = "RdPu", trans = "date", direction = -1, na.value = 'grey75') +
  scico::scale_fill_scico(palette = "vik", trans = "date", direction = -1, na.value = 'grey75') +
  guides(fill = guide_colorbar(direction = 'vertical')) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'right',
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'white', color = 'white'),
        strip.text = element_text(size = 11, face = 'bold'))

# Figure out which combinations of parameters are realistic# Figure out which combinations of parameters are realistic# Figure out which combinations of parameters are realistic
# new_data2 %>% mutate(y_realistic = ifelse())
ggplot() +
  geom_hex(data = optim_clim %>% filter(planting_doy > 90, Year >1980),
           aes(x = GDD_modave, 
               y = as.Date("2019-12-31") + planting_doy),
           # bins = 20,
           binwidth = c(35,4),
           color = 'black', fill = 'black', alpha = 0.5,
           linewidth = 0.9) 






## Climate change heat map ----------------------------------
# the following works

# what is the anomaly from current climate
fut_anomaly <- optim_clim_sc %>%
  summarise(
    curr_GDD_anom = mean(GDD_anomaly),
    fut_GDD_anom = mean(GDD_future_anomaly),
    fut_GDD_anom_sc = (fut_GDD_anom - ranges_unsc[2,2])/ranges_unsc[2,1]
  ) %>% as.numeric()


# Quantiles of GDD anomaly (colder year, middle year, hotter year)
gdd_anom_quants_fut <- quantile(optim_clim_sc$GDD_anomaly_sc, c(.05, .50, .95)) %>% as.numeric()

# function that will generate sequence with a buffer around the min and max
# row = 1
# myseq2 <- function(row){seq(ranges[row,1] - 2, ranges[row,2] + 3, length.out = reps)}

# effect of GDD site mean and planting day at three levels of GDD_anomaly
new_data2_fut <- rbind(
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants_fut[1], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants_fut[2], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants_fut[3], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(fut_anomaly[3], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  )
)

# predict y points
pred_y2 = NULL
for(i in 1:nrow(new_data2_fut)){
  pred_i <- median_qi(betas$b_Intercept + 
                        betas$b1 * new_data2_fut$GDD_mod_sitemean_sc[i] +
                        betas$b2 * new_data2_fut$GDD_anomaly_sc[i] +
                        betas$b3 * new_data2_fut$tot_precip_sitemean_sc[i] +
                        betas$b4 * new_data2_fut$precip_anomaly_sc[i] +
                        betas$b5 * new_data2_fut$frost_anomaly_sc[i] +
                        betas$b6 * new_data2_fut$planting_doy[i])
  pred_y2 <- rbind(pred_y2, pred_i)
}

# combine prediction summaries to the new data
new_data2_fut <- cbind(new_data2_fut, pred_y2) 

new_data2_fut <- new_data2_fut %>% mutate(GDD_mod_sitemean = 
                                    GDD_mod_sitemean_sc * ranges_unsc[1,1] + ranges_unsc[1,2]#,
)

# rename levels of anomaly, for facet labels
new_data2_fut$GDD_anomaly_sc <- factor(new_data2_fut$GDD_anomaly_sc, 
                                   labels = c("Cold year\n(5th percentile)", 
                                              "Average year\n(50th percentile)", 
                                              "Warm year\n(95th percentile)",
                                              "Projected mid-century\naverage year"))

# grey out unrealistic points
new_data2_fut <- new_data2_fut %>% mutate(realistic = 
                                    ifelse(apply(cbind(new_data2_fut$GDD_mod_sitemean,new_data2_fut$planting_doy), 1, function (p){
                                      any(realistic_ranges$GDD_modave_1 <= p[1] & 
                                            realistic_ranges$GDD_modave_2 >= p[1] &
                                            realistic_ranges$planting_doy1 <= p[2] & 
                                            realistic_ranges$planting_doy2 >= p[2])
                                    } ) ,'yes', 'no'),
                                  y_realistic = ifelse(realistic == 'yes', y, NA))

# remove extra points
new_data2_fut <- new_data2_fut %>%
  filter(GDD_mod_sitemean > 480 & GDD_mod_sitemean < 1280 & planting_doy > 117 & planting_doy < 198)

# add lables with flowering dates to each panel
date_labs <- new_data2_fut %>% filter(between(GDD_mod_sitemean, 879, 916) & 
                           between(planting_doy, 152.5, 156) ) %>% 
  select(GDD_anomaly_sc, y) %>%
  mutate(label = (format(as.Date("2019-12-31") + y, "%b-%d") ),
         ABC = LETTERS[1:4])
# 
# date_labs <- data.frame(GDD_anomaly_sc = c("Cold year\n(5th percentile)", 
#                       "Average year\n(50th percentile)", 
#                       "Warm year\n(95th percentile)",
#                       "Projected mid-century\naverage year"),
#            label = c("Aug. 16", "Aug. 10", "Aug. 4", "Aug. 2"))
date_labs$GDD_anomaly_sc <- factor(date_labs$GDD_anomaly_sc, levels = c("Cold year\n(5th percentile)", 
                                                                        "Average year\n(50th percentile)", 
                                                                        "Warm year\n(95th percentile)",
                                                                        "Projected mid-century\naverage year"))

# heat map 
ggplot(data = new_data2_fut, aes(x = GDD_mod_sitemean, 
                             y = as.Date("2019-12-31") + planting_doy, 
                             fill = as.Date("2019-12-31") + y_realistic)) +
  geom_tile(linewidth = 0.1) +
  facet_wrap(~GDD_anomaly_sc, ncol = 4) +
  # geom_hex(data = optim_clim %>% filter(planting_doy > 90, Year > 1980),
  #          aes(x = GDD_modave, 
  #              y = as.Date("2019-12-31") + planting_doy),
  #          bins = 20, color = 'black', fill = 'white', alpha = 0,
  #          linewidth = 0.9) +
  annotate(geom = 'rect', xmin = 879, xmax = 916, 
           ymin = as.Date("2019-12-31") + 152.5, ymax = as.Date("2019-12-31") + 156,
           alpha = 0, color = 'black', linewidth = 1) +
  geom_text(data = date_labs, aes(x = 650, y = as.Date("2019-12-31") + 190,
                                  label = label, group = GDD_anomaly_sc), inherit.aes = F,
            # label.size = 1, fill = 'grey75',
            size = 4.5) +
  annotate(geom = 'rect', xmin = 555, xmax = 555+185,
           ymin = as.Date("2019-12-31") + 184, ymax = as.Date("2019-12-31") + 195,
           alpha = 0, color = 'black', linewidth = 1) +
  labs(fill = "Optimal flowering \nday of year",
       y = "Planting day of year",
       x = "Current average site degree days") +
  geom_text(data = date_labs, aes(x = 510, y = as.Date("2019-12-31") + 196.5,
                                  label = ABC, group = GDD_anomaly_sc), inherit.aes = F,
            size = 5, fontface = 'bold') +
  scale_y_date(date_breaks = '2 weeks', date_labels = "%b %d", expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  # viridis::scale_fill_viridis(trans = "date", direction = -1, na.value = 'grey75') +
  # scale_fill_gradientn(colors = terrain.colors(3), trans = "date", na.value = 'grey75') +
  # scale_fill_gradient(low = 'blue', high = 'orange', trans = 'date', na.value = 'grey75') +
  # scale_fill_gradientn(colors = c('palegreen3', 'white', 'sienna3'), scales::rescale(c(-.5, 0, 0.5)), trans = "date", na.value = 'grey75')+
  # scico::scale_fill_scico(palette = "vik", trans = "date", direction = -1, na.value = 'grey75') +
  scico::scale_fill_scico(palette = "cork", trans = "date", direction = -1, na.value = 'grey75') +
  guides(fill = guide_colorbar(direction = 'vertical')) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'right',
        # panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'white', color = 'white'),
        strip.text = element_text(size = 11, face = 'bold'))

ggsave(paste0(path2plots,"DOY_future.png"),
       width = 16, height = 6)


## Predictor Residuals -----------------------------------------------------

# predict GDD_site_mean_sc with planting_doy

resid_mod_1 <- brm(GDD_mod_sitemean_sc ~ planting_doy, 
                    data = optim_clim_sc,
                    family = gaussian,
                    prior = c(prior(normal(0, 100), class = Intercept),
                              prior(normal(0, 1), class = b),
                              prior(normal(0, 10), class = sigma)),
                    chains = 4, cores = 4,
                    iter = 2000, warmup = 500
)

print(resid_mod_1)

f_1 <- fitted(resid_mod_1) %>% as_tibble() %>% bind_cols(optim_clim_sc)

res_p_1 <- f_1 %>%
  ggplot(aes(x = as.Date("2019-12-31") + planting_doy, y = GDD_mod_sitemean_sc)) +
  geom_point(aes(color = latitude)) +
  geom_line(aes(y = Estimate)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(x = "Planting DOY", y = "GDD site mean (scaled)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

r_1 <- residuals(resid_mod_1) %>% as_tibble() %>% bind_cols(optim_clim_sc)

#!! alternative with residuals with all predictors
# Only run the following two lines to make the following plots with residuals
# from model with all predictors. Object r_1 will overwrite previous object!!
# resid_mod_1_all <- brm(GDD_mod_sitemean_sc ~ planting_doy + GDD_anomaly_sc + 
#                          tot_precip_sitemean_sc + precip_anomaly_sc + frost_anomaly_sc, 
#                    data = optim_clim_sc,
#                    family = gaussian,
#                    prior = c(prior(normal(0, 100), class = Intercept),
#                              prior(normal(0, 1), class = b),
#                              prior(normal(0, 10), class = sigma)),
#                    chains = 4, cores = 4,
#                    iter = 2000, warmup = 500
# )
# r_1 <- residuals(resid_mod_1_all) %>% as_tibble() %>% bind_cols(optim_clim_sc)

# as.Date("2019-12-31") +
res_p_2 <- ggplot(data = r_1 %>% filter(outside_data == 'in_range' & !is.na(sd)),
         aes(x = Estimate, y = as.Date("2019-12-31") + mu)
         ) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_linerange(
    data = r_1 %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = Estimate,
        ymax =  as.Date("2019-12-31") + y3,
        ymin =  as.Date("2019-12-31") + y2,
        color = latitude), inherit.aes = F,
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = r_1 %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = Estimate,
        y = as.Date("2019-12-31") + mu,
        ymax = as.Date("2019-12-31") + mu+sd,
        ymin = as.Date("2019-12-31") + mu-sd,
        color = latitude),
    alpha = 0.5
  ) +
  # error of prediction from model
  geom_ribbon(
    data = predictions_list_unsc[[1]],
    aes(x = x, ymin = as.Date("2019-12-31") + ymin, ymax = as.Date("2019-12-31") + ymax),
    fill = "grey75", alpha = 0.75, inherit.aes = F
  ) +
  geom_line(
    data = predictions_list_unsc[[1]],
    aes(x = x, y = as.Date("2019-12-31") + y), linewidth = 0.8, linetype = 1
  ) +
  annotate("text",
           x = -Inf, y = as.Date("2019-12-31") + Inf,
           label = paste0("Slope = ", format(round(coefs_unsc[1,11], 3), nsmall = 3),
                          "\nCI: ", round(coefs_unsc[1,12],3), ", ",
                          round(coefs_unsc[1,13],3)),
           color = 'red', size = 5, hjust = 0, vjust = 1
  ) +
  scale_y_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  # stat_smooth(method = 'lm', fullrange = T, color = 'black') +
  labs(x = "GDD site mean residuals", y = "Optimal flowering date") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
res_p_2

# fit model the other way for planting doy residuals
resid_mod_2 <- brm(planting_doy ~ GDD_mod_sitemean_sc, 
                   data = optim_clim_sc,
                   family = gaussian,
                   prior = c(prior(normal(0, 100), class = Intercept),
                             prior(normal(0, 1), class = b),
                             prior(normal(0, 10), class = sigma)),
                   chains = 4, cores = 4,
                   iter = 2000, warmup = 500
)

print(resid_mod_2)

r_2 <- residuals(resid_mod_2) %>% as_tibble() %>% bind_cols(optim_clim_sc)

#!! alternative with residuals with all predictors
# Only run the following two lines to make the following plots with residuals
# from model with all predictors. Object r_2 will overwrite previous object!!
# resid_mod_2_all <- brm(planting_doy ~ GDD_mod_sitemean_sc + GDD_anomaly_sc + 
#                          tot_precip_sitemean_sc + precip_anomaly_sc + frost_anomaly_sc, 
#                        data = optim_clim_sc,
#                        family = gaussian,
#                        prior = c(prior(normal(0, 100), class = Intercept),
#                                  prior(normal(0, 1), class = b),
#                                  prior(normal(0, 10), class = sigma)),
#                        chains = 4, cores = 4,
#                        iter = 2000, warmup = 500
# )
# r_2 <- residuals(resid_mod_2_all) %>% as_tibble() %>% bind_cols(optim_clim_sc)


# c(sd(optim_clim_sc$frost_anomaly), mean(optim_clim_sc$frost_anomaly))
planting_doy_rescaled <- predictions_list_unsc[[6]] %>%
  mutate(x_resc = (x - mean(optim_clim_sc$planting_doy))) #/ sd(optim_clim_sc$planting_doy))

# as.Date("2019-12-31") +
res_p_3 <- ggplot(data = r_2 %>% filter(outside_data == 'in_range' & !is.na(sd)),
       aes(x = Estimate, y = as.Date("2019-12-31") + mu)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_linerange(
    data = r_2 %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = Estimate,
        ymax =  as.Date("2019-12-31") + y3,
        ymin =  as.Date("2019-12-31") + y2,
        color = latitude), inherit.aes = F,
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = r_2 %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = Estimate,
        y = as.Date("2019-12-31") + mu,
        ymax = as.Date("2019-12-31") + mu+sd,
        ymin = as.Date("2019-12-31") + mu-sd,
        color = latitude),
    alpha = 0.5
  ) +
  geom_ribbon(
    data = planting_doy_rescaled,
    aes(x = x_resc, ymin = as.Date("2019-12-31") + ymin, ymax = as.Date("2019-12-31") + ymax),
    fill = "grey75", alpha = 0.75, inherit.aes = F
  ) +
  geom_line(
    data = planting_doy_rescaled,
    aes(x = x_resc, y = as.Date("2019-12-31") + y), linewidth = 0.8, linetype = 1
  ) +
  annotate("text",
           x = -Inf, y = as.Date("2019-12-31") + Inf,
           label = paste0("Slope = ", format(round(coefs_unsc[6,11], 3), nsmall = 3),
                          "\nCI: ", round(coefs_unsc[6,12],3), ", ",
                          format(round(coefs_unsc[6,13],3)), nsmall = 3),
           color = 'red', size = 5, hjust = 0, vjust = 1
  ) +
  scale_y_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  # stat_smooth(method = 'lm', fullrange = T, color = 'black') +
  labs(x = "Planting day of year residuals", y = "Optimal flowering date") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
res_p_3

ggarrange(res_p_1, res_p_2, res_p_3, nrow = 1, widths = c(0.85,1,1), common.legend = T, legend = 'right')


# DOY MODEL - FIXED GDD ---------------------------------------------------

## Stan data -------------------------------------------------------------
# use brms to format the data as per a censored model
sdata4 <- make_standata(
  y2 | cens(cens, y3) ~ GDD_fixed_sitemean_sc + GDD_fixed_anomaly_sc +
    tot_precip_sitemean_sc + precip_anomaly_sc +
    # frost_sitemean_sc +
    frost_anomaly_sc +
    planting_doy + (1 | group),
  data = optim_clim_sc,
  family = gaussian
)

# add in the estimated peak (mu) and se on peak (sd) from
# the gaussian models
sdata4$sigma_meas <- optim_clim_sc$sd
sdata4$measured_y <- optim_clim_sc$mu


# pass stan a new variable that is sigma_meas but if that doesn't exist
# just fill with a very big number. If you pass the data in with NAs I think
# they will be dropped becauase Stan hates NAs
# these never enter the likelihood so it's fine
# !! but we should triple check at some point these rows don't have
# crazy fitted values, which they would if these values were entering the
# likelihood.
# or check by making a different crazy value and ensuring you get the
# same results
sdata4$se <- ifelse(is.na(sdata4$sigma_meas), 9999, sdata4$sigma_meas)
sdata4$measured_y <- ifelse(is.na(sdata4$sigma_meas), -9999, sdata4$measured_y)

# create a new indicator value,
# 2 if interval censored, 3 if measurement error
# these could be any value, I kept 2 for interval censored in keeping
# with how brms codes it

## !! this should be made off 'in range' instead for 3
# if data is "in range" and it has an estimate of sd, then cens_2 should be 3 (measurement error),
# if out of range or does not have a sd measured, then 2 (interval censored)
sdata4$cens_2 <- ifelse(optim_clim_sc$outside_data == "in_range" & 
                          !is.na(optim_clim_sc$sd), 3, 2)


# remove the original censored variable that the make_standata function
# put in. We won't use this anymore.
sdata4$cens <- NULL

sdata4$sigma_meas <- NULL # EC: this contains NAs and isn't used later

## Run model -----------------------------------------------------------------

# actual model needs to be run with longer burnin/warmup
fit4 <- stan(
  model_code = stanmod_2, # Stan program
  data = sdata4, # named list of data
  chains = 4, # number of Markov chains
  warmup = 500, # number of warmup iterations per chain
  iter = 4500, # total number of iterations per chain
  cores = 4, # number of cores (could use one per chain)
  # refresh = 0,             # no progress shown
  init = inits_list, # init for intercept needs to be specified
  control = list(max_treedepth = 10)
)

# names of parameters
ggs(fit4, burnin = TRUE) %>%
  distinct(Parameter)

# trace plot !! doesn't really work right now
ggs(fit4, burnin = TRUE) %>%
  filter(Parameter == "b_Intercept") %>%
  mutate(
    chain = factor(Chain),
    intercept = value
  ) %>%
  ggplot(aes(x = Iteration, y = intercept, color = chain)) +
  geom_trace(subtitle = "fit4 (default settings)") +
  scale_y_continuous(breaks = c(0, 650, 1300), limits = c(NA, 1430))

## Results -----------------------------------------------------------------
# sanity check the x matrix is column order by variable name
# as provided in the model formula
# head (optim_clim_sc$tot_precip_sitemean_sc)
# head (sdata1$X[,4])
print(summary(fit4)$summary)

# can look at summary of fixed effects only
summary(fit4, pars = c("Intercept", "b"), probs = c(0.05, 0.95))$summary

# tidybayes way to look at more or less the same thing
# I like that these color the ones where the CI is different from 0
# I did spot check that the data has the covariates entered in the
# order they are entered in the model formula statement.
fit4 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_fixed_sitemean_sc",
        param_order == 2 ~ "GDD_fixed_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  group_by(param_name) %>%
  median_qi() %>%
  print()

summarise_draws(tidy_draws(fit4))



## Make predictions -------------------------------------------------------


# get min, mean, and max of all variables, to see the interval on which to predict
ranges_fixed <- data.frame(GDD_site = c(range(optim_clim_sc$GDD_fixed_sitemean_sc), mean(optim_clim_sc$GDD_fixed_sitemean_sc)),
                     GDD_anom = c(range(optim_clim_sc$GDD_fixed_anomaly_sc), mean(optim_clim_sc$GDD_fixed_anomaly_sc)),
                     precip_site = c(range(optim_clim_sc$tot_precip_sitemean_sc), mean(optim_clim_sc$tot_precip_sitemean_sc)),
                     precip_anom = c(range(optim_clim_sc$precip_anomaly_sc), mean(optim_clim_sc$precip_anomaly_sc)),
                     frost_anom = c(range(optim_clim_sc$frost_anomaly_sc), mean(optim_clim_sc$frost_anomaly_sc)),
                     planting_doy = c(range(optim_clim_sc$planting_doy), mean(optim_clim_sc$planting_doy))
)
rownames(ranges_fixed) <- c("min", "max", 'mean')
ranges_fixed <- t(ranges_fixed)


# functions to create sequences and reps for new data to predict on
myseq_fixed <- function(row){seq(ranges_fixed[row,1], ranges_fixed[row,2], length.out = reps)}
myrep_fixed <- function(row){rep(ranges_fixed[row,3], reps)}

# create new data on which to predict
# this is marginal predictions - one variable change at a time, the rest are at their mean
new_data_fixed <- data.frame(GDD_fixed_sitemean_sc = c(myseq_fixed(1), rep(myrep_fixed(1), 5)),
                       GDD_fixed_anomaly_sc = c(myrep_fixed(2), myseq_fixed(2), rep(myrep_fixed(2), 4)),
                       tot_precip_sitemean_sc = c(rep(myrep_fixed(3), 2), myseq_fixed(3), rep(myrep_fixed(3), 3)),
                       precip_anomaly_sc = c(rep(myrep_fixed(4), 3), myseq_fixed(4), rep(myrep_fixed(4), 2)),
                       frost_anomaly_sc = c(rep(myrep_fixed(5), 4), myseq_fixed(5), rep(myrep_fixed(5), 1)),
                       planting_doy = c(rep(myrep_fixed(6), 5), myseq_fixed(6)),
                       group = rep(c('GDD_fixed_sitemean_sc', 'GDD_fixed_anomaly_sc', 
                                     'tot_precip_sitemean_sc', 'precip_anomaly_sc', 
                                     'frost_anomaly_sc', 'planting_doy'), each = reps)
)

# Extract draws of betas and intercept
betas_fixed <- fit4 %>%
  tidy_draws() %>%
  select(starts_with("b")) 
colnames(betas_fixed) <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b_Intercept')

# For each row of new data, multiply it out by the betas, 
# then take the median_qi and add that summary to the table
pred_y = NULL
for(i in 1:nrow(new_data_fixed)){
  pred_i <- median_qi(betas_fixed$b_Intercept + 
                        betas_fixed$b1 * new_data_fixed$GDD_fixed_sitemean_sc[i] +
                        betas_fixed$b2 * new_data_fixed$GDD_fixed_anomaly_sc[i] +
                        betas_fixed$b3 * new_data_fixed$tot_precip_sitemean_sc[i] +
                        betas_fixed$b4 * new_data_fixed$precip_anomaly_sc[i] +
                        betas_fixed$b5 * new_data_fixed$frost_anomaly_sc[i] +
                        betas_fixed$b6 * new_data_fixed$planting_doy[i])
  pred_y <- rbind(pred_y, pred_i)
}

# combine prediction summaries to the new data
new_data_fixed1 <- cbind(new_data_fixed, pred_y) 

# # for plotting, we need the predictor that varies (x) and the prediction (y)
# # and we don't need the other variables that are set to their means
# # put these x's and y's into elements of a list 
# make_pred_list <- function(predictor, new.dat){
#   new.dat %>% 
#     filter(group == predictor) %>% 
#     select({{predictor}}, 7:10) %>%
#     rename(x = {{predictor}}) %>%
#     relocate(group, .after = ymax)
# }
predictions_list_fixed <- as.list(unique(new_data_fixed1$group))
predictions_list_fixed <- lapply(predictions_list_fixed, make_pred_list, new.dat = new_data_fixed1)

# unscale the predictors
ranges_unsc_fixed <- data.frame(GDD_fixed_sitemean_sc = c(sd(optim_clim_sc$GDD_fixed_sitemean), mean(optim_clim_sc$GDD_fixed_sitemean)),
                                GDD_fixed_anomaly_sc = c(sd(optim_clim_sc$GDD_fixed_anomaly), mean(optim_clim_sc$GDD_fixed_anomaly)),
                                tot_precip_sitemean_sc = c(sd(optim_clim_sc$tot_precip_sitemean), mean(optim_clim_sc$tot_precip_sitemean)),
                                precip_anomaly_sc = c(sd(optim_clim_sc$precip_anomaly), mean(optim_clim_sc$precip_anomaly)),
                                frost_anomaly_sc = c(sd(optim_clim_sc$frost_anomaly), mean(optim_clim_sc$frost_anomaly)),
                                planting_doy = c(1, 0)
)
rownames(ranges_unsc_fixed) <- c("sd", 'mean')
ranges_unsc_fixed <- t(ranges_unsc_fixed)

# little function to add a column of the unscaled predictors 
# to the predictions list
unscale_predictions_fixed <- function(predictions){
  group_i <- unique(predictions$group)
  predictions <- predictions %>% mutate(x_unsc = 
                                          predictions$x * 
                                          ranges_unsc_fixed[which(rownames(ranges_unsc_fixed)==group_i),1] +
                                          ranges_unsc_fixed[which(rownames(ranges_unsc_fixed)==group_i),2],
                                        .after = x)
}

predictions_list_fixed_unsc <- lapply(predictions_list_fixed, unscale_predictions_fixed)

# Plot of trends -----------------------------------------------------------

# get a table of the coefficents to put on the plot
coefs_fixed <- fit4 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_fixed_sitemean_sc",
        param_order == 2 ~ "GDD_fixed_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  group_by(param_name) %>%
  median_qi() %>% arrange(param_order)

coefs_fixed_unsc <- coefs_fixed %>% mutate(b_unsc = b / ranges_unsc_fixed[,1],
                                                   b.lower_unsc = b.lower / ranges_unsc_fixed[,1],
                                                   b.upper_unsc = b.upper / ranges_unsc_fixed[,1])

# make a table of x-axis labels
labels_fixed <- data.frame(predictor = gsub('_sc', '', coefs_fixed$param_name),
                           label = c("Fixed period GDD site mean (°C)",
                                     "Fixed period GDD anomaly (°C)",
                                     "Precipitation site mean (mm)",
                                     "Precipitation anomaly (mm)",
                                     "First frost anomaly",
                                     "Planting day of year")
)

# function that will make each plot, with trend line, 
# points and intervals of data, and coefficients
# dat = predictions_list_fixed_unsc[[1]]
plot_trends2_fixed <- function(dat){
  predictor <- gsub('_sc', '', as.character(unique(dat$group)))
  coef_i <- coefs_fixed_unsc %>%
    filter(grepl(predictor, param_name)) %>%
    select(b_unsc, b.lower_unsc, b.upper_unsc) %>%
    round(digits = 3)
  
  ggplot() +
    # vertical lines for censoring intervals
    geom_linerange(
      data = optim_clim_sc %>% filter(outside_data != 'in_range' | is.na(sd)),
      aes(x = .data[[predictor]], 
          ymax = as.Date("2019-12-31") + y3, 
          ymin = as.Date("2019-12-31") + y2, color = latitude),
      alpha = 0.5
    ) +
    # points and lines for measurement error points
    geom_pointrange(
      data = optim_clim_sc %>% filter(outside_data == 'in_range' & !is.na(sd)),
      aes(x = .data[[predictor]], 
          y = as.Date("2019-12-31") + mu, 
          ymax = as.Date("2019-12-31") + mu+sd, 
          ymin = as.Date("2019-12-31") + mu-sd, color = latitude),
      alpha = 0.5
    ) +
    # error of prediction from model
    geom_ribbon(
      data = dat,
      aes(x = x_unsc, 
          ymin = as.Date("2019-12-31") + ymin, 
          ymax = as.Date("2019-12-31") + ymax),
      fill = "grey75", alpha = 0.75
    ) +
    # trend line of prediction from model
    geom_line(
      data = dat,
      aes(x = x_unsc, y = as.Date("2019-12-31") + y), 
      linetype = as.numeric(ifelse(coef_i[2] * coef_i[3] > 0, 1, 2)), 
      linewidth = 0.8
    ) +
    # coefficient annotation
    # annotate("text",
    #          x = -Inf, y = as.Date("2019-12-31") + Inf,
    #          label = paste0("Slope = ", coef_i[1], "\nCI: ", coef_i[2], ", ", coef_i[3]),
    #          color = ifelse(coef_i[2] * coef_i[3] > 0, "red", "black"), size = 5, hjust = 0, vjust = 1
    # ) +
    # formatting
    scale_y_date(date_breaks = "2 weeks", date_labels = "%b %d") +
    labs(
      x = labels_fixed[which(labels_fixed$predictor == predictor),2],
      y = "Optimal flowering date",
      color = "Latitude"
    ) +
    # guides(color = "none") +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank())
}

# run the function
regression_plots_fixed <- lapply(predictions_list_fixed_unsc, plot_trends2_fixed)
regression_plots_fixed <- regression_plots_fixed[c(1, 2, 6, 3, 4, 5)]

# arrange resulting plots and add a title
annotate_figure(
  ggarrange(plotlist = regression_plots_fixed, common.legend = T, 
            legend = 'right', labels = 'AUTO', nrow = 2, ncol = 3),
  top = text_grob("Fixed GDD period & flowering date model", face = "bold", size = 16))
ggsave(paste0(path2plots, "fixed_DOY_opt_clim.png"), 
       width = 13, height = 8, units = "in")

## Predictor Residuals -----------------------------------------------------

# predict GDD_fixed_sitemean_sc with planting_doy

resid_mod_1_fixed <- brm(GDD_fixed_sitemean_sc ~ planting_doy,
                         data = optim_clim_sc,
                         family = gaussian,
                         prior = c(prior(normal(0, 100), class = Intercept),
                                   prior(normal(0, 1), class = b),
                                   prior(normal(0, 10), class = sigma)),
                         chains = 4, cores = 4,
                         iter = 2000, warmup = 500
)

print(resid_mod_1_fixed)

f_1_fixed <- fitted(resid_mod_1_fixed) %>% as_tibble() %>% bind_cols(optim_clim_sc)

res_p_1_fixed_doy <- f_1_fixed %>%
  ggplot(aes(x = as.Date("2019-12-31") + planting_doy, y = GDD_fixed_sitemean_sc)) +
  geom_point(aes(color = latitude)) +
  geom_line(aes(y = Estimate)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(x = "Planting DOY", y = "Fixed GDD site mean (scaled)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

r_1_fixed <- residuals(resid_mod_1_fixed) %>% as_tibble() %>% bind_cols(optim_clim_sc)

res_p_2_fixed <- ggplot(data = r_1_fixed %>% filter(outside_data == 'in_range' & !is.na(sd)),
                            aes(x = Estimate, y =  mu)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_linerange(
    data = r_1_fixed %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = Estimate,
        ymax =   y3,
        ymin =   y2,
        color = latitude), inherit.aes = F,
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = r_1_fixed %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = Estimate,
        y = mu,
        ymax = mu+sd,
        ymin =  mu-sd,
        color = latitude),
    alpha = 0.5
  ) +
  stat_smooth(method = 'lm', fullrange = T, color = 'black') +
  labs(x = "Fixed GDD site mean residuals", y = "Optimal flowering day of year") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

resid_mod_2_fixed <- brm(planting_doy ~ GDD_fixed_sitemean_sc,
                         data = optim_clim_sc,
                         family = gaussian,
                         prior = c(prior(normal(0, 100), class = Intercept),
                                   prior(normal(0, 1), class = b),
                                   prior(normal(0, 10), class = sigma)),
                         chains = 4, cores = 4,
                         iter = 2000, warmup = 500
)

print(resid_mod_2_fixed)

r_2_fixed <- residuals(resid_mod_2_fixed) %>% as_tibble() %>% bind_cols(optim_clim_sc)

# as.Date("2019-12-31") +
res_p_3_fixed <- ggplot(data = r_2_fixed %>% filter(outside_data == 'in_range' & !is.na(sd)),
                            aes(x = Estimate, y = mu)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_linerange(
    data = r_2_fixed %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = Estimate,
        ymax =  y3,
        ymin =  y2,
        color = latitude), inherit.aes = F,
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = r_2_fixed %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = Estimate,
        y =  mu,
        ymax =  mu+sd,
        ymin =  mu-sd,
        color = latitude),
    alpha = 0.5
  ) +
  stat_smooth(method = 'lm', fullrange = T, color = 'black') +
  labs(x = "Planting day of year residuals", y = "Optimal flowering day of year") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

ggarrange(res_p_1_fixed_doy, res_p_2_fixed, res_p_3_fixed, nrow = 1, common.legend = T, legend = 'right')

# compare correlations of GDD and fixed GDD side-by-side
ggarrange(res_p_1, res_p_1_fixed_doy, nrow = 1, common.legend = T, legend = 'right')













#############################################################################.
#############################################################################.
#############################################################################.
# DPP MODEL ---------------------------------------------------
# Option b - model as days after planting.
# to do this I think you just have to subtract the planting_doy (scalar)
# from the optimal flowering code. Alternatively you could go back and remake
# the optimals directly on the days after planting, but as it's just subtracting
# a constant I don't think it would differ

# for simplicity so as not to have to rename everything in the model we'll
# just make a new dataset with the same variables y2, y3, and mu
# representing days after planting
optim_clim_sc_2 <- optim_clim_sc %>%
  mutate(
    y2 = y2 - planting_doy,
    y3 = y3 - planting_doy,
    mu = mu - planting_doy
  )

# use brms to format the data as per a censored model
sdata2 <- make_standata(
  y2 | cens(cens, y3) ~ GDD_mod_sitemean_sc + GDD_anomaly_sc +
    tot_precip_sitemean_sc + precip_anomaly_sc +
    frost_anomaly_sc +
    planting_doy + (1 | group),
  data = optim_clim_sc_2,
  family = gaussian
)


sdata2$sigma_meas <- optim_clim_sc$sd
sdata2$measured_y <- optim_clim_sc$mu

# bogus values bc stan hates NAs
# these never enter the likelihood so it's fine
sdata2$se <- ifelse(is.na(sdata2$sigma_meas), 9999, sdata2$sigma_meas)
sdata2$measured_y <- ifelse(is.na(sdata2$sigma_meas), -9999, sdata2$measured_y)

# create a new indicator value,
# 2 if interval censored, 3 if measurement error
sdata2$cens_2 <- ifelse(!is.na(optim_clim_sc$sd), 3, 2)
sdata2$cens <- NULL

# define the stanvars
mean_y <- mean(optim_clim_sc_2$mu, na.rm = T)

# get some reasonable initial values
# !! ideally you want to vary these over your chains though
# seems like they are mixing ok?
inits_2 <- list(Intercept = mean(optim_clim_sc_2$mu, na.rm = T))
inits_list_2 <- list(inits_2, inits_2, inits_2, inits_2)

# refit as days past planting is the response variable
fit2 <- stan(
  model_code = stanmod_2, # Stan program
  data = sdata2, # named list of data
  chains = 4, # number of Markov chains
  warmup = 500, # number of warmup iterations per chain
  iter = 4500, # total number of iterations per chain
  cores = 4, # number of cores (could use one per chain)
  # refresh = 0,             # no progress shown
  # init = inits_list,
  control = list(max_treedepth = 10),
  init = inits_list_2
)

## Results -----------------------------------------------------------------
summarise_draws(tidy_draws(fit2))

print(summary(fit2)$summary)

summary(fit2, pars = c("Intercept", "b"), probs = c(0.05, 0.95))$summary

fit2 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_mod_sitemean_sc",
        param_order == 2 ~ "GDD_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  group_by(param_name) %>%
  median_qi() %>%
  print()





## Make predictions -------------------------------------------------------

# Extract draws of betas and intercept from days past planting
betas_dpp <- fit2 %>%
  tidy_draws() %>%
  select(starts_with("b")) 
colnames(betas_dpp) <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b_Intercept')

# For each row of new data, multiply it out by the betas, 
# then take the median_qi and add that summary to the table
pred_y = NULL
for(i in 1:nrow(new_data)){
  pred_i <- median_qi(betas_dpp$b_Intercept + 
                        betas_dpp$b1 * new_data$GDD_mod_sitemean_sc[i] +
                        betas_dpp$b2 * new_data$GDD_anomaly_sc[i] +
                        betas_dpp$b3 * new_data$tot_precip_sitemean_sc[i] +
                        betas_dpp$b4 * new_data$precip_anomaly_sc[i] +
                        betas_dpp$b5 * new_data$frost_anomaly_sc[i] +
                        betas_dpp$b6 * new_data$planting_doy[i])
  pred_y <- rbind(pred_y, pred_i)
}

# combine prediction summaries to the new data
new_data_dpp <- cbind(new_data, pred_y) 

# make the list of predictions for each predictor 
predictions_list_dpp <- as.list(unique(new_data_dpp$group))
predictions_list_dpp <- lapply(predictions_list_dpp, make_pred_list, new.dat = new_data_dpp)

predictions_list_dpp_unsc <- lapply(predictions_list_dpp, unscale_predictions)




## Plots of trends ---------------------------------------------------------

# # make a table of x-axis labels
# labels <- data.frame(predictor = gsub('_sc', '', coefs_unsc$param_name),
#                      label = c("Degree days site mean (°C)",
#                                "Degree days yearly anomaly (°C)",
#                                "Precipitation site mean (mm)", 
#                                "Precipitation yearly anomaly (mm)",
#                                "Frost day of year yearly anomaly",
#                                "Planting day of year")
# )

# get a table of the coefficents to put on the plot
coefs_dpp <- fit2 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_mod_sitemean_sc",
        param_order == 2 ~ "GDD_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  group_by(param_name) %>%
  median_qi() %>% arrange(param_order)

coefs_dpp_unsc <- coefs_dpp %>% mutate(b_unsc = b / ranges_unsc[,1],
                               b.lower_unsc = b.lower / ranges_unsc[,1],
                               b.upper_unsc = b.upper / ranges_unsc[,1])


# function that will make each plot, with trend line, 
# points and intervals of data, and coefficients
# dat = predictions_list_dpp_unsc[[6]]
plot_trends2_dpp <- function(dat){
  predictor <- gsub('_sc', '', as.character(unique(dat$group)))
  coef_i <- coefs_dpp_unsc %>%
    filter(grepl(predictor, param_name)) %>%
    select(b_unsc, b.lower_unsc, b.upper_unsc) %>%
    round(digits = 3)
  
  ggplot() +
    # vertical lines for censoring intervals
    geom_linerange(
      data = optim_clim_sc_2 %>% filter(outside_data != 'in_range' | is.na(sd)),
      aes(x = .data[[predictor]], ymax = y3, ymin = y2, color = latitude),
      alpha = 0.5
    ) +
    # points and lines for measurement error points
    geom_pointrange(
      data = optim_clim_sc_2 %>% filter(outside_data == 'in_range' & !is.na(sd)),
      aes(x = .data[[predictor]], y = mu, ymax = mu+sd, ymin = mu-sd, color = latitude),
      alpha = 0.5
    ) +
    # error of prediction from model
    geom_ribbon(
      data = dat,
      aes(x = x_unsc, ymin = ymin, ymax = ymax),
      fill = "grey75", alpha = 0.75
    ) +
    # trend line of prediction from model
    geom_line(
      data = dat,
      aes(x = x_unsc, y = y), 
      linetype = as.numeric(ifelse(coef_i[2] * coef_i[3] > 0, 1, 2)), 
      linewidth = 0.8
    ) +
    # coefficient annotation
    # annotate("text",
    #          x = -Inf, y = Inf,
    #          label = paste0("Slope = ", coef_i[1], "\nCI: ", coef_i[2], ", ", coef_i[3]),
    #          color = ifelse(coef_i[2] * coef_i[3] > 0, "red", "black"), size = 5, hjust = 0, vjust = 1
    # ) +
    # formatting
    labs(
      x = labels[which(labels$predictor == predictor),2],
      y = "Optimal days to flowering",
      color = "Latitude"
    ) +
    # guides(color = "none") +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank())
}

# run the function
regression_plots_dpp <- lapply(predictions_list_dpp_unsc, plot_trends2_dpp)
regression_plots_dpp <- regression_plots_dpp[c(1, 2, 6, 3, 4, 5)]

# arrange resulting plots
ggarrange(plotlist = regression_plots_dpp, common.legend = T, legend = 'right',
          labels = "AUTO", ncol = 3, nrow = 2)
ggsave(paste0(path2plots, "DPP_opt_clim.png"), 
       width = 13, height = 8)
ggsave(paste0(path2plots, "DPP_opt_clim.pdf"), 
       width = 13, height = 8, dpi = 600)

## Predictor Residuals -----------------------------------------------------

# predict GDD_site_mean_sc with planting_doy

# resid_mod_1 <- brm(GDD_mod_sitemean_sc ~ planting_doy, 
#                    data = optim_clim_sc,
#                    family = gaussian,
#                    prior = c(prior(normal(0, 100), class = Intercept),
#                              prior(normal(0, 1), class = b),
#                              prior(normal(0, 10), class = sigma)),
#                    chains = 4, cores = 4,
#                    iter = 2000, warmup = 500
# )
# 
# print(resid_mod_1)
# 
f_1_dpp <- fitted(resid_mod_1) %>% as_tibble() %>% bind_cols(optim_clim_sc_2)

res_p_1 <- f_1_dpp %>%
  ggplot(aes(x = as.Date("2019-12-31") + planting_doy, y = GDD_mod_sitemean_sc)) +
  geom_point(aes(color = latitude)) +
  geom_line(aes(y = Estimate)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(x = "Planting date", y = "GDD site mean (scaled)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

r_1_dpp <- residuals(resid_mod_1) %>% as_tibble() %>% bind_cols(optim_clim_sc_2)

res_p_2_dpp <- ggplot(data = r_1_dpp %>% filter(outside_data == 'in_range' & !is.na(sd)),
                  aes(x = Estimate, y =  mu)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_linerange(
    data = r_1_dpp %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = Estimate,
        ymax =   y3,
        ymin =   y2,
        color = latitude), inherit.aes = F,
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = r_1_dpp %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = Estimate,
        y = mu,
        ymax = mu+sd,
        ymin =  mu-sd,
        color = latitude),
    alpha = 0.5
  ) +
  annotate("text",
           x = -Inf, y = as.Date("2019-12-31") + Inf,
           label = paste0("Slope = ", format(round(coefs_dpp_unsc[1,11], 3), nsmall = 3),
                          "\nCI: ", round(coefs_dpp_unsc[1,12],3), ", ",
                          format(round(coefs_dpp_unsc[1,13],3)), nsmall = 3),
           color = 'red', size = 5, hjust = 0, vjust = 1
  ) +
  stat_smooth(method = 'lm', fullrange = T, color = 'black') +
  labs(x = "GDD site mean residuals", y = "Optimal days to flowering") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
res_p_2_dpp

# resid_mod_2 <- brm(planting_doy ~ GDD_mod_sitemean_sc, 
#                    data = optim_clim_sc,
#                    family = gaussian,
#                    prior = c(prior(normal(0, 100), class = Intercept),
#                              prior(normal(0, 1), class = b),
#                              prior(normal(0, 10), class = sigma)),
#                    chains = 4, cores = 4,
#                    iter = 2000, warmup = 500
# )
# 
# print(resid_mod_2)

r_2_dpp <- residuals(resid_mod_2) %>% as_tibble() %>% bind_cols(optim_clim_sc_2)

# as.Date("2019-12-31") +
res_p_3_dpp <- ggplot(data = r_2_dpp %>% filter(outside_data == 'in_range' & !is.na(sd)),
                  aes(x = Estimate, y = mu)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_linerange(
    data = r_2_dpp %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = Estimate,
        ymax =  y3,
        ymin =  y2,
        color = latitude), inherit.aes = F,
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = r_2_dpp %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = Estimate,
        y =  mu,
        ymax =  mu+sd,
        ymin =  mu-sd,
        color = latitude),
    alpha = 0.5
  ) +
  annotate("text",
           x = -Inf, y = as.Date("2019-12-31") + Inf,
           label = paste0("Slope = ", format(round(coefs_dpp_unsc[6,11], 3), nsmall = 3),
                          "\nCI: ", round(coefs_dpp_unsc[6,12],3), ", ",
                          format(round(coefs_dpp_unsc[6,13],3)), nsmall = 3),
           color = 'black', size = 5, hjust = 0, vjust = 1
  ) +
  stat_smooth(method = 'lm', fullrange = T, color = 'black') +
  labs(x = "Planting day of year residuals", y = "Optimal days to flowering") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
res_p_3_dpp

ggarrange(res_p_1, res_p_2_dpp, res_p_3_dpp, nrow = 1, common.legend = T, 
          legend = 'right', labels = "AUTO")

##!!!!!
# ggsave("temp_plots/opt_manuscript_figures/DPP_resids.png", 
#        width = 12, height = 4)

## Plot of temp x planting doy ----------------------------------------------

# effect of GDD site mean and planting day at three levels of GDD_anomaly
# make all combinations of x predictors
new_data2_dpp <- rbind(
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants[1], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-2, ranges[6,2]+1, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants[2], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-2, ranges[6,2]+1, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants[3], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-2, ranges[6,2]+1, length.out = reps), each = reps)
  )
)

# predict y points
pred_y2 = NULL
for(i in 1:nrow(new_data2_dpp)){
  pred_i <- median_qi(betas_dpp$b_Intercept + 
                        betas_dpp$b1 * new_data2_dpp$GDD_mod_sitemean_sc[i] +
                        betas_dpp$b2 * new_data2_dpp$GDD_anomaly_sc[i] +
                        betas_dpp$b3 * new_data2_dpp$tot_precip_sitemean_sc[i] +
                        betas_dpp$b4 * new_data2_dpp$precip_anomaly_sc[i] +
                        betas_dpp$b5 * new_data2_dpp$frost_anomaly_sc[i] +
                        betas_dpp$b6 * new_data2_dpp$planting_doy[i])
  pred_y2 <- rbind(pred_y2, pred_i)
}

# combine prediction summaries to the new data
new_data2_dpp <- cbind(new_data2_dpp, pred_y2) 

# unscale GDD_sitemean (GDD anomaly doesn't need to be unscaled, because 
# it's a factor in this graph (only 3 levels))
new_data2_dpp <- new_data2_dpp %>% mutate(GDD_mod_sitemean = 
                                    GDD_mod_sitemean_sc * 
                                    ranges_unsc[1,1] +
                                    ranges_unsc[1,2])

new_data2_dpp$GDD_anomaly_sc <- factor(new_data2_dpp$GDD_anomaly_sc, 
                                   labels = c("Colder year\n(25th percentile)", 
                                              "Average year\n(50th percentile)", 
                                              "Warmer year\n(75th percentile)"))

# add columns that say if the planting day/site GDD is actually realistic
# uses objects defined above
new_data2_dpp <- new_data2_dpp %>% mutate(realistic = 
                                    ifelse(apply(cbind(new_data2_dpp$GDD_mod_sitemean, new_data2_dpp$planting_doy), 1, function(p){
                                      any(realistic_ranges$GDD_modave_1 <= p[1] & 
                                            realistic_ranges$GDD_modave_2 >= p[1] &
                                            realistic_ranges$planting_doy1 <= p[2] & 
                                            realistic_ranges$planting_doy2 >= p[2])
                                    }
                                    ) ,
                                    'yes', 'no'),
                                  y_realistic = ifelse(realistic == 'yes', y, NA))

ggplot() +
  # error of prediction from model
  geom_ribbon(
    data = new_data2_dpp,
    aes(x = GDD_mod_sitemean, ymin = ymin, ymax = ymax, group = as.factor(planting_doy)),
    fill = "grey75", alpha = 0.75
  ) +
  # trend line of prediction from model
  geom_line(
    data = new_data2_dpp,
    aes(x = GDD_mod_sitemean, y = y, color = as.factor(planting_doy)), linewidth = 0.8
  ) +
  facet_wrap(~GDD_anomaly_sc) +
  # formatting
  labs(
    y = "Optimal days past planting",
    color = "Planting day"
  ) +
  # guides(color = "none") +
  theme_bw(base_size = 16)


# heat map of optimal flowering day of year
ggplot(data = new_data2_dpp, aes(x = GDD_mod_sitemean, 
                             y = as.Date("2019-12-31") + planting_doy, 
                             fill = y_realistic)) +
  geom_tile(linewidth = 0.1) +
  facet_wrap(~GDD_anomaly_sc) +
  # geom_hex(data = optim_clim %>% filter(planting_doy > 92, Year > 1980),
  #          aes(x = GDD_modave,
  #              y = as.Date("2019-12-31") + planting_doy),
  #          bins = 20, color = 'black', fill = 'white', alpha = 0,
  #          linewidth = 0.9) +
  labs(fill = "Optimal flowering \ndays past planting",
       y = "Planting day of year",
       x = "Average site degree days") +
  scale_y_date(date_breaks = '2 weeks', date_labels = "%b %d", expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  # viridis::scale_fill_viridis(trans = "date") +
  scico::scale_fill_scico(palette = "vik", na.value = 'grey75') +
  guides(fill = guide_colorbar(direction = 'vertical')) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'right',
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'white', color = 'white'),
        strip.text = element_text(size = 11, face = 'bold'))

# Figure out which combinations of parameters are realistic# Figure out which combinations of parameters are realistic# Figure out which combinations of parameters are realistic
# new_data2 %>% mutate(y_realistic = ifelse())
ggplot() +
  geom_hex(data = optim_clim %>% filter(planting_doy > 92, Year > 1980),
           aes(x = GDD_modave, 
               y = as.Date("2019-12-31") + planting_doy),
           # bins = 20,
           binwidth = c(35,4),
           color = 'black', fill = 'black', alpha = 0.5,
           linewidth = 0.9) 


## Climate change heat map DPP -------------------------------------------
# effect of GDD site mean and planting day at three levels of GDD_anomaly
new_data2_fut_dpp <- rbind(
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants_fut[1], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants_fut[2], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(gdd_anom_quants_fut[3], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  ),
  data.frame(GDD_mod_sitemean_sc = rep(myseq2(1), reps),
             GDD_anomaly_sc = rep(fut_anomaly[3], reps*reps),
             tot_precip_sitemean_sc = rep(myrep(3), reps),
             precip_anomaly_sc = rep(myrep(4), reps),
             frost_anomaly_sc = rep(myrep(5), reps),
             planting_doy = rep(seq(ranges[6,1]-10, ranges[6,2]+10, length.out = reps), each = reps)
  )
)

# predict y points
pred_y2 = NULL
for(i in 1:nrow(new_data2_fut_dpp)){
  pred_i <- median_qi(betas_dpp$b_Intercept +
                        betas_dpp$b1 * new_data2_fut_dpp$GDD_mod_sitemean_sc[i] +
                        betas_dpp$b2 * new_data2_fut_dpp$GDD_anomaly_sc[i] +
                        betas_dpp$b3 * new_data2_fut_dpp$tot_precip_sitemean_sc[i] +
                        betas_dpp$b4 * new_data2_fut_dpp$precip_anomaly_sc[i] +
                        betas_dpp$b5 * new_data2_fut_dpp$frost_anomaly_sc[i] +
                        betas_dpp$b6 * new_data2_fut_dpp$planting_doy[i])
  pred_y2 <- rbind(pred_y2, pred_i)
}

# combine prediction summaries to the new data
new_data2_fut_dpp <- cbind(new_data2_fut_dpp, pred_y2)

new_data2_fut_dpp <- new_data2_fut_dpp %>% mutate(GDD_mod_sitemean =
                                            GDD_mod_sitemean_sc * ranges_unsc[1,1] + ranges_unsc[1,2]#,
)

# rename levels of anomaly, for facet labels
new_data2_fut_dpp$GDD_anomaly_sc <- factor(new_data2_fut_dpp$GDD_anomaly_sc,
                                       labels = c("Cold year\n(5th percentile)",
                                                  "Average year\n(50th percentile)",
                                                  "Warm year\n(95th percentile)",
                                                  "Projected mid-century\naverage year"))

# grey out unrealistic points
new_data2_fut_dpp <- new_data2_fut_dpp %>% mutate(realistic =
                                            ifelse(apply(cbind(new_data2_fut_dpp$GDD_mod_sitemean,new_data2_fut_dpp$planting_doy), 1, function (p){
                                              any(realistic_ranges$GDD_modave_1 <= p[1] &
                                                    realistic_ranges$GDD_modave_2 >= p[1] &
                                                    realistic_ranges$planting_doy1 <= p[2] &
                                                    realistic_ranges$planting_doy2 >= p[2])
                                            } ) ,'yes', 'no'),
                                          y_realistic = ifelse(realistic == 'yes', y, NA))

# remove extra points
new_data2_fut_dpp <- new_data2_fut_dpp %>%
  filter(GDD_mod_sitemean > 480 & GDD_mod_sitemean < 1280 & planting_doy > 117 & planting_doy < 198)

# add lables with flowering dates to each panel
date_labs_dpp <- new_data2_fut_dpp %>% filter(between(GDD_mod_sitemean, 879, 916) &
                                        between(planting_doy, 152.5, 156) ) %>%
  select(GDD_anomaly_sc, y) %>%
  mutate(label = paste0(round(y, 0), '\ndays' ),
         ABC = LETTERS[1:4])

date_labs_dpp$GDD_anomaly_sc <- factor(date_labs_dpp$GDD_anomaly_sc, levels = c("Cold year\n(5th percentile)",
                                                                        "Average year\n(50th percentile)",
                                                                        "Warm year\n(95th percentile)",
                                                                        "Projected mid-century\naverage year"))

# heat map
ggplot(data = new_data2_fut_dpp, aes(x = GDD_mod_sitemean,
                                 y = as.Date("2019-12-31") + planting_doy,
                                 fill = y_realistic)) +
  geom_tile(linewidth = 0.1) +
  facet_wrap(~GDD_anomaly_sc, ncol = 4) +
  # geom_hex(data = optim_clim %>% filter(planting_doy > 90, Year > 1980),
  #          aes(x = GDD_modave,
  #              y = as.Date("2019-12-31") + planting_doy),
  #          bins = 20, color = 'black', fill = 'white', alpha = 0,
  #          linewidth = 0.9) +
  annotate(geom = 'rect', xmin = 879, xmax = 916,
           ymin = as.Date("2019-12-31") + 152.5, ymax = as.Date("2019-12-31") + 156,
           alpha = 0, color = 'black', linewidth = 1) +
  geom_text(data = date_labs_dpp, aes(x = 650, y = as.Date("2019-12-31") + 190,
                                  label = label, group = GDD_anomaly_sc), inherit.aes = F,
            # label.size = 1, fill = 'grey75',
            size = 4.5) +
  annotate(geom = 'rect', xmin = 555, xmax = 555+185,
           ymin = as.Date("2019-12-31") + 184, ymax = as.Date("2019-12-31") + 195,
           alpha = 0, color = 'black', linewidth = 1) +
  labs(fill = "Optimal days\nto flowering",
       y = "Planting day of year",
       x = "Current average site degree days") +
  geom_text(data = date_labs_dpp, aes(x = 510, y = as.Date("2019-12-31") + 196.5,
                                  label = ABC, group = GDD_anomaly_sc), inherit.aes = F,
            size = 5, fontface = 'bold') +
  scale_y_date(date_breaks = '2 weeks', date_labels = "%b %d", expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  # viridis::scale_fill_viridis(trans = "date", direction = -1, na.value = 'grey75') +
  # scale_fill_gradientn(colors = terrain.colors(3), trans = "date", na.value = 'grey75') +
  # scale_fill_gradient(low = 'blue', high = 'orange', trans = 'date', na.value = 'grey75') +
  # scale_fill_gradientn(colors = c('palegreen3', 'white', 'sienna3'), scales::rescale(c(-.5, 0, 0.5)), trans = "date", na.value = 'grey75')+
  scico::scale_fill_scico(palette = "cork", direction = -1, na.value = 'grey75') +
  guides(fill = guide_colorbar(direction = 'vertical')) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'right',
        # panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'white', color = 'white'),
        strip.text = element_text(size = 11, face = 'bold'))

ggsave(paste0(path2plots, "DPP_future.png"),
       width = 16, height = 6)
ggsave(paste0(path2plots, "DPP_future.pdf"),
       width = 16, height = 6, dpi = 600)


# quantiles of projected days to flowering in the future
# and other stuff that can be uncommented
new_data2_fut_dpp %>% 
  filter(!is.na(y_realistic),
         GDD_anomaly_sc == "Projected mid-century\naverage year") %>%
  # what is the y_realistic at 50th percentile
  summarise(median = median(y_realistic),
            lower = quantile(y_realistic, 0.05),
            upper = quantile(y_realistic, 0.95))
  # ggplot() +
  # geom_density(aes(x = y_realistic, fill = "Future Climates"),
  #              alpha = 0.3) 
  # proportion of realistic points that are below 55 days
  # summarise(prop_70 = sum(y_realistic < 70) / n(),
  #           prop_60 = sum(y_realistic < 60) / n(),
  #           prop_58 = sum(y_realistic < 58) / n())

# quantiles of current average optimal days to flowering
new_data2_fut_dpp %>% 
  filter(!is.na(y_realistic),
         GDD_anomaly_sc == "Average year\n(50th percentile)") %>%
  # what is the y_realistic at 50th percentile for current average
  summarise(median = median(y_realistic),
            lower = quantile(y_realistic, 0.05),
            upper = quantile(y_realistic, 0.95))


############################################################################.
#############################################################################.
#############################################################################.
#############################################################################.




# DPP MODEL - FIXED GDD ------------------------------------------------------

# use the same optim_clim_sc_2 dataset as for DPP above
optim_clim_sc_2 <- optim_clim_sc %>%
  mutate(
    y2 = y2 - planting_doy,
    y3 = y3 - planting_doy,
    mu = mu - planting_doy
  )

# use brms to format the data as per a censored model
# instead of GDD_mod_sitemean_sc and GDD_anomaly_sc
# we are using GDD_fixed_sitemean_sc and GDD_fixed_anomaly_sc
sdata3 <- make_standata(
  y2 | cens(cens, y3) ~ GDD_fixed_sitemean_sc + GDD_fixed_anomaly_sc +
    tot_precip_sitemean_sc + precip_anomaly_sc +
    frost_anomaly_sc +
    planting_doy + (1 | group),
  data = optim_clim_sc_2,
  family = gaussian
)


sdata3$sigma_meas <- optim_clim_sc$sd
sdata3$measured_y <- optim_clim_sc$mu

# bogus values bc stan hates NAs
# these never enter the likelihood so it's fine
sdata3$se <- ifelse(is.na(sdata3$sigma_meas), 9999, sdata3$sigma_meas)
sdata3$measured_y <- ifelse(is.na(sdata3$sigma_meas), -9999, sdata3$measured_y)

# create a new indicator value,
# 2 if interval censored, 3 if measurement error
sdata3$cens_2 <- ifelse(!is.na(optim_clim_sc$sd), 3, 2)
sdata3$cens <- NULL

# define the stanvars
mean_y <- mean(optim_clim_sc_2$mu, na.rm = T)

# get some reasonable initial values
# !! ideally you want to vary these over your chains though
# seems like they are mixing ok?
inits_2 <- list(Intercept = mean(optim_clim_sc_2$mu, na.rm = T))
inits_list_2 <- list(inits_2, inits_2, inits_2, inits_2)

# refit as days past planting is the response variable
fit3 <- stan(
  model_code = stanmod_2, # Stan program
  data = sdata3, # named list of data
  chains = 4, # number of Markov chains
  warmup = 500, # number of warmup iterations per chain
  iter = 4500, # total number of iterations per chain
  cores = 4, # number of cores (could use one per chain)
  # refresh = 0,             # no progress shown
  # init = inits_list,
  control = list(max_treedepth = 10),
  init = inits_list_2
)

## Results -----------------------------------------------------------------
summarise_draws(tidy_draws(fit3))

print(summary(fit3)$summary)

summary(fit3, pars = c("Intercept", "b"), probs = c(0.05, 0.95))$summary

fit3 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_fixed_sitemean_sc",
        param_order == 2 ~ "GDD_fixed_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  group_by(param_name) %>%
  median_qi() %>%
  print()


## Make predictions -------------------------------------------------------

# Extract draws of betas and intercept from days past planting
betas_dpp_fixed <- fit3 %>%
  tidy_draws() %>%
  select(starts_with("b")) 
colnames(betas_dpp_fixed) <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b_Intercept')

# For each row of new data, multiply it out by the betas, 
# then take the median_qi and add that summary to the table
pred_y = NULL
for(i in 1:nrow(new_data_fixed)){
  pred_i <- median_qi(betas_dpp_fixed$b_Intercept + 
                        betas_dpp_fixed$b1 * new_data_fixed$GDD_fixed_sitemean_sc[i] +
                        betas_dpp_fixed$b2 * new_data_fixed$GDD_fixed_anomaly_sc[i] +
                        betas_dpp_fixed$b3 * new_data_fixed$tot_precip_sitemean_sc[i] +
                        betas_dpp_fixed$b4 * new_data_fixed$precip_anomaly_sc[i] +
                        betas_dpp_fixed$b5 * new_data_fixed$frost_anomaly_sc[i] +
                        betas_dpp_fixed$b6 * new_data_fixed$planting_doy[i])
  pred_y <- rbind(pred_y, pred_i)
}

# combine prediction summaries to the new data
new_data_dpp_fixed <- cbind(new_data_fixed, pred_y) 

# make the list of predictions for each predictor 
predictions_list_dpp_fixed <- as.list(unique(new_data_dpp_fixed$group))
predictions_list_dpp_fixed <- lapply(predictions_list_dpp_fixed, make_pred_list, new.dat = new_data_dpp_fixed)

# Unscale predictions
predictions_list_dpp_fixed_unsc <- lapply(predictions_list_dpp_fixed, unscale_predictions_fixed)



## Plots of trends ---------------------------------------------------------



# get a table of the coefficents to put on the plot
coefs_dpp_fixed <- fit3 %>%
  spread_draws(b[param_order]) %>%
  mutate(
    param_name =
      case_when(
        param_order == 1 ~ "GDD_fixed_sitemean_sc",
        param_order == 2 ~ "GDD_fixed_anomaly_sc",
        param_order == 3 ~ "tot_precip_sitemean_sc",
        param_order == 4 ~ "precip_anomaly_sc",
        param_order == 5 ~ "frost_anomaly_sc",
        param_order == 6 ~ "planting_doy",
        TRUE ~ NA_character_
      )
  ) %>%
  group_by(param_name) %>%
  median_qi() %>% arrange(param_order)

coefs_dpp_fixed_unsc <- coefs_dpp_fixed %>% mutate(b_unsc = b / ranges_unsc_fixed[,1],
                                       b.lower_unsc = b.lower / ranges_unsc_fixed[,1],
                                       b.upper_unsc = b.upper / ranges_unsc_fixed[,1])

# make a table of x-axis labels
labels_fixed <- data.frame(predictor = gsub('_sc', '', coefs_dpp_fixed$param_name),
                     label = c("Fixed period GDD site mean (°C)",
                               "Fixed period GDD anomaly (°C)",
                               "Precipitation site mean (mm)",
                               "Precipitation anomaly (mm)",
                               "First frost anomaly",
                               "Planting day of year")
)

# function that will make each plot, with trend line, 
# points and intervals of data, and coefficients
# dat = predictions_list_dpp_fixed_unsc[[1]]
plot_trends2_dpp_fixed <- function(dat){
  predictor <- gsub('_sc', '', as.character(unique(dat$group)))
  coef_i <- coefs_dpp_fixed_unsc %>%
    filter(grepl(predictor, param_name)) %>%
    select(b_unsc, b.lower_unsc, b.upper_unsc) %>%
    round(digits = 3)
  
  ggplot() +
    # vertical lines for censoring intervals
    geom_linerange(
      data = optim_clim_sc_2 %>% filter(outside_data != 'in_range' | is.na(sd)),
      aes(x = .data[[predictor]], ymax = y3, ymin = y2, color = latitude),
      alpha = 0.5
    ) +
    # points and lines for measurement error points
    geom_pointrange(
      data = optim_clim_sc_2 %>% filter(outside_data == 'in_range' & !is.na(sd)),
      aes(x = .data[[predictor]], y = mu, ymax = mu+sd, ymin = mu-sd, color = latitude),
      alpha = 0.5
    ) +
    # error of prediction from model
    geom_ribbon(
      data = dat,
      aes(x = x_unsc, ymin = ymin, ymax = ymax),
      fill = "grey75", alpha = 0.75
    ) +
    # trend line of prediction from model
    geom_line(
      data = dat,
      aes(x = x_unsc, y = y), 
      linetype = as.numeric(ifelse(coef_i[2] * coef_i[3] > 0, 1, 2)), 
      linewidth = 0.8
    ) +
    # coefficient annotation
    # annotate("text",
    #          x = -Inf, y = Inf,
    #          label = paste0("Slope = ", coef_i[1], "\nCI: ", coef_i[2], ", ", coef_i[3]),
    #          color = ifelse(coef_i[2] * coef_i[3] > 0, "red", "black"), size = 5, hjust = 0, vjust = 1
    # ) +
    # formatting
    labs(
      x = labels_fixed[which(labels_fixed$predictor == predictor),2],
      y = "Optimal days to flowering",
      color = "Latitude"
    ) +
    # guides(color = "none") +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank())
}

# run the function
regression_plots_dpp_fixed <- lapply(predictions_list_dpp_fixed_unsc, plot_trends2_dpp_fixed)
regression_plots_dpp_fixed <- regression_plots_dpp_fixed[c(1, 2, 6, 3, 4, 5)]

# # arrange resulting plots
# ggarrange(plotlist = regression_plots_dpp_fixed, common.legend = T, legend = 'right')

# arrange resulting plots and add a title
annotate_figure(
  ggarrange(plotlist = regression_plots_dpp_fixed, common.legend = T, 
            legend = 'right', labels = 'AUTO', nrow = 2, ncol = 3),
  top = text_grob("Fixed GDD period & days to flowering model", face = "bold", size = 16))
ggsave(paste0(path2plots, "fixed_DPP_opt_clim.png"), 
       width = 13, height = 8, units = "in")

## Predictor Residuals -----------------------------------------------------

# predict GDD_fixed_sitemean_sc with planting_doy

resid_mod_1_fixed <- brm(GDD_fixed_sitemean_sc ~ planting_doy,
                   data = optim_clim_sc,
                   family = gaussian,
                   prior = c(prior(normal(0, 100), class = Intercept),
                             prior(normal(0, 1), class = b),
                             prior(normal(0, 10), class = sigma)),
                   chains = 4, cores = 4,
                   iter = 2000, warmup = 500
)

print(resid_mod_1_fixed)

f_1_dpp_fixed <- fitted(resid_mod_1_fixed) %>% as_tibble() %>% bind_cols(optim_clim_sc_2)

res_p_1_fixed <- f_1_dpp_fixed %>%
  ggplot(aes(x = as.Date("2019-12-31") + planting_doy, y = GDD_fixed_sitemean_sc)) +
  geom_point(aes(color = latitude)) +
  geom_line(aes(y = Estimate)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(x = "Planting DOY", y = "Fixed GDD site mean (scaled)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

r_1_dpp_fixed <- residuals(resid_mod_1_fixed) %>% as_tibble() %>% bind_cols(optim_clim_sc_2)

res_p_2_dpp_fixed <- ggplot(data = r_1_dpp_fixed %>% filter(outside_data == 'in_range' & !is.na(sd)),
                      aes(x = Estimate, y =  mu)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_linerange(
    data = r_1_dpp_fixed %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = Estimate,
        ymax =   y3,
        ymin =   y2,
        color = latitude), inherit.aes = F,
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = r_1_dpp_fixed %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = Estimate,
        y = mu,
        ymax = mu+sd,
        ymin =  mu-sd,
        color = latitude),
    alpha = 0.5
  ) +
  stat_smooth(method = 'lm', fullrange = T, color = 'black') +
  labs(x = "Fixed GDD site mean residuals", y = "Optimal days past planting") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

resid_mod_2_fixed <- brm(planting_doy ~ GDD_fixed_sitemean_sc,
                   data = optim_clim_sc,
                   family = gaussian,
                   prior = c(prior(normal(0, 100), class = Intercept),
                             prior(normal(0, 1), class = b),
                             prior(normal(0, 10), class = sigma)),
                   chains = 4, cores = 4,
                   iter = 2000, warmup = 500
)

print(resid_mod_2_fixed)

r_2_dpp_fixed <- residuals(resid_mod_2_fixed) %>% as_tibble() %>% bind_cols(optim_clim_sc_2)

# as.Date("2019-12-31") +
res_p_3_dpp_fixed <- ggplot(data = r_2_dpp_fixed %>% filter(outside_data == 'in_range' & !is.na(sd)),
                      aes(x = Estimate, y = mu)) +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
  geom_linerange(
    data = r_2_dpp_fixed %>% filter(outside_data != 'in_range' | is.na(sd)),
    aes(x = Estimate,
        ymax =  y3,
        ymin =  y2,
        color = latitude), inherit.aes = F,
    alpha = 0.5
  ) +
  # points and lines for measurement error points
  geom_pointrange(
    data = r_2_dpp_fixed %>% filter(outside_data == 'in_range' & !is.na(sd)),
    aes(x = Estimate,
        y =  mu,
        ymax =  mu+sd,
        ymin =  mu-sd,
        color = latitude),
    alpha = 0.5
  ) +
  stat_smooth(method = 'lm', fullrange = T, color = 'black') +
  labs(x = "Planting day of year residuals", y = "Optimal days past planting") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

ggarrange(res_p_1_fixed, res_p_2_dpp_fixed, res_p_3_dpp_fixed, nrow = 1, common.legend = T, legend = 'right')

# compare correlations of GDD and fixed GDD side-by-side
ggarrange(res_p_1, res_p_1_fixed, common.legend = T, legend = 'right')


# ALL RESIDUAL PLOTS --------------------------------------------------------
ggarrange(
  fixed_vs_relative_gdd,
  ggarrange(res_p_2_dpp, res_p_3_dpp, res_p_2, res_p_3, 
            nrow = 2, ncol = 2, common.legend = T, legend = 'right', 
            labels = c("B", "C", "D", "E")),
  nrow = 2, ncol = 1, common.legend = T, 
  heights = c(0.7,1),
  legend = 'none', labels = (c("A", ""))
  )

ggsave(paste0(path2plots,"ALL_resids.png"),
       width = 12, height = 16)

# TABLES of All Results ---------------------------------------------------------

# data frame mapping labels for predictors
table_labels <- data.frame(predictor = c(coefs_unsc$param_name, coefs_fixed_unsc$param_name[1:2]),
                           label = c("GDD site mean",
                                     "GDD anomaly",
                                     "Precip site mean",
                                     "Precip anomaly",
                                     "First frost anomaly",
                                     "Planting date",
                                     "site-GDD site mean",
                                     "site-GDD anomaly")
)

# print separate tables for DPP and the other models, when running by county
if(by_county == TRUE){
  


# table for DPP results
tt_dpp <- coefs_dpp_unsc %>% 
  mutate(ci = paste0(round(b.lower_unsc, 3), ", ", round(b.upper_unsc, 3)),
         b_unsc = round(b_unsc, 3),
         Model = "DPP") %>%
  left_join(., table_labels, by = c("param_name" = "predictor")) %>%
  select(Model, label, b_unsc, ci) %>% 
  as.data.frame() %>%
  rename("Predictor" = 'label',
         "Slope" = "b_unsc",
         "95% CI" = "ci") %>%
  tt(theme = 'void') %>%
  style_tt(i = which(coefs_dpp_unsc$b.lower_unsc * coefs_dpp_unsc$b.upper_unsc > 0) %>% as.numeric,
           j = c(1:4),
           bold = TRUE) %>%
  style_tt(i = 1, j = 1,
           rowspan = 6) %>%
  style_tt(i = c(1), line= "t") %>%
  style_tt(i = c(6), line= "b", line_width = 0.15)
tt_dpp

tt_dpp %>%
  save_tt(paste0(path2plots,"table_dpp.png"), overwrite = T)


## table for all other results
all_table_coefs <- rbind(
# coefs_dpp_unsc %>% mutate(Model = "DPP"),
coefs_unsc %>% mutate(Model = "DOY"),
coefs_dpp_fixed_unsc %>% mutate(Model = "DPP (site-GDD)"),
coefs_fixed_unsc %>% mutate(Model = "DOY (site-GDD)")
) 
tt_all <- all_table_coefs %>%
  mutate(ci = paste0(round(b.lower_unsc, 3), ", ", round(b.upper_unsc, 3)),
         b_unsc = round(b_unsc, 3)) %>%
  left_join(., table_labels, by = c("param_name" = "predictor")) %>%
  select(Model, label, b_unsc, ci) %>% 
  as.data.frame() %>%
  rename("Predictor" = 'label',
         "Slope" = "b_unsc",
         "95% CI" = "ci") %>%
  tt(theme = 'void') %>%
  style_tt(i = which(all_table_coefs$b.lower_unsc * all_table_coefs$b.upper_unsc > 0) %>% as.numeric,
           j = c(1:4),
           bold = TRUE) %>%
  style_tt(i = c(1, 7, 13, 19), j = 1,
           rowspan = 6, align ="r") %>%
  style_tt(i = c(7, 13), line= "t") %>%
  style_tt(i = c(1), line= "t", line_width = 0.15) %>%
  style_tt(i = c(18), line= "b", line_width = 0.15)
tt_all

tt_all %>%
  save_tt(paste0(path2plots,"table_all.png"), overwrite = T)
}


# print all results as a single table, for supplement
if(by_county == FALSE){
  ## table for all results
  all_table_coefs <- rbind(
    coefs_dpp_unsc %>% mutate(Model = "DPP"),
    coefs_unsc %>% mutate(Model = "DOY"),
    coefs_dpp_fixed_unsc %>% mutate(Model = "DPP (site-GDD)"),
    coefs_fixed_unsc %>% mutate(Model = "DOY (site-GDD)")
  ) 
  tt_all <- all_table_coefs %>%
    mutate(ci = paste0(round(b.lower_unsc, 3), ", ", round(b.upper_unsc, 3)),
           b_unsc = round(b_unsc, 3)) %>%
    left_join(., table_labels, by = c("param_name" = "predictor")) %>%
    select(Model, label, b_unsc, ci) %>% 
    as.data.frame() %>%
    rename("Predictor" = 'label',
           "Slope" = "b_unsc",
           "95% CI" = "ci") %>%
    tt(theme = 'void') %>%
    style_tt(i = which(all_table_coefs$b.lower_unsc * all_table_coefs$b.upper_unsc > 0) %>% as.numeric,
             j = c(1:4),
             bold = TRUE) %>%
    style_tt(i = c(1, 7, 13, 19), j = 1,
             rowspan = 6, align ="r") %>%
    style_tt(i = c(7, 13, 19), line= "t") %>%
    style_tt(i = c(1), line= "t", line_width = 0.15) %>%
    style_tt(i = c(24), line= "b", line_width = 0.15)
  tt_all
  
  tt_all %>%
    save_tt(paste0(path2plots,"table_all.png"), overwrite = T)
}


#############################################################################.
#############################################################################.
#############################################################################.
#############################################################################.
#############################################################################.
############################################################################.
############################################################################.
############################################################################.
# how the sausage was made. -------------------------------------------------
# used the brms make_stancode function to get the stan code
# for the censored model. Wrote to a text file and inspected
# make stancode for the censored model does require user specificed
# prior on intercept and sigma

# Take the stancode you want and modify the likelihood
# mod<-brms::make_stancode(formula = brmsformula(y2 | cens(cens, y3) ~ GDD_mod_sitemean + GDD_anomaly +
#       tot_precip_sitemean + precip_anomaly +
#       frost_sitemean + frost_anomaly +
#       planting_doy + (1|named_location)),
#     data = optim_clim,
#     family = gaussian,
#     prior = c(prior(normal(221.3135, 30.02485 * 100), class = Intercept),
#               prior(normal(0,30.02485), class = sigma)))
# writeLines(mod, "clipboard")


# repeat for a measurement error model
# inspect what is different (fortunately most of it is identical so the
# data set up and parameter estimates largely the same! Just need to
# pass a few more variables in the data and make the likelihood differ
# based on whether it's a measurement error site-yr or a censored site-yr
# tst_meas<-brms::make_stancode(formula = brmsformula(mu | se(sd, sigma = TRUE) ~ GDD_mod_sitemean + GDD_anomaly +
#                                                              tot_precip_sitemean + precip_anomaly +
#                                                              frost_sitemean + frost_anomaly +
#                                                              planting_doy + (1|named_location)),
#                         data = optim_clim,
#                         family = 'gaussian')
# writeLines(tst_meas, "clipboard")

# stan itself doesn't have all these, may have to code up separately or
# switch to rstanarm
# bayes_R2(fit1)
# loo_cv/waic things would also be good to do but we need to add the loglikelihood
# to the stan model first
# loo::loo(fit1)

