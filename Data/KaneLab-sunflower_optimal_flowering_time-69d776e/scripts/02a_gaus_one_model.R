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



# Fit Gaussian model to all sites in a single model ----------------------


# Sites
sites_county <- list_sites(deriv_data, type = 'all_county', exclude_double = TRUE)

# collapse site_county into a dataframe
sites_county_df <- do.call(rbind.data.frame, sites_county)


# AIC subset
R2_county_data <- AIC_R2_data %>% filter(subset == 'all_countyloc')


# Filter deriv_data to only include sites in sites_county_df and 
# exclude no yield and no phenology
deriv_site_1mod <- deriv_data %>%
  semi_join(., sites_county_df, by = c('Year', 'Location', 'named_location')) %>%
  filter(
    # named_location %in% site$named_location &
    #        Location %in% site$Location &
    #        lat %in% site$lat &
    #        lon %in% site$lon &
    !is.na(flower_50_doy) &
      !is.na(yield_lb_acre))

# Exclude sites that have a flowering day spread of less than 5 days
exclude_low_spread_1mod <- deriv_site_1mod %>%
  group_by(Year, named_location) %>% 
  summarize(spr = max(flower_50_doy) - min(flower_50_doy) + 1) %>% 
  arrange(spr) %>% 
  filter(spr < 5)
deriv_site_1mod <- deriv_site_1mod %>% anti_join(., exclude_low_spread_1mod, by = join_by(Year, named_location))

# Don't need to exclude singletons, so leave as FALSE
exclude_singleton_sites = FALSE

if(exclude_singleton_sites) {
# Exclude sites that only have one year
exclude_singletons_1mod <- deriv_site_1mod %>%
  distinct(county_state, Year) %>%
  # add_count(Year, named_location)
  group_by(county_state) %>%
  summarise(n = n()) %>%
  filter(n == 1)
deriv_site_1mod <- deriv_site_1mod %>% anti_join(., exclude_singletons_1mod, by = join_by(Year, named_location))
}

# create scaled flowering day
deriv_site_1mod$flower_50_doy_sc <- scale(deriv_site_1mod$flower_50_doy, scale = T)

# remove this site-yr that doesn't run
# if("Cheyenne_HPAL" %in% unique(deriv_site$named_location)) {
#   deriv_site <- deriv_site %>% filter(Year != 2013)
# }

# Recode Year to be Site-Year
deriv_site_1mod <- deriv_site_1mod %>%
  mutate(Year_actual = Year,
         Year = as.factor(as.numeric(factor(paste(named_location, Year_actual, sep = "_")))+1000))

# Fit Gaussian model
mod.lin_1mod <- glm(yield_lb_acre ~ 0 + Year + 
                 flower_50_doy_sc:Year + 
                 I(flower_50_doy_sc^2):Year,
               family = gaussian(link = "log"),
               data = deriv_site_1mod)
summary(mod.lin_1mod)

# grabbing coefficients
coef.nam <- names(coef(mod.lin_1mod))

# sanity check to make sure there are 3 coeffs per year
if(any(table(gsub(":.*", "", coef.nam))!=3)){
  print('Check your coefficients! Seems wrong!')
}  

coef.uniq <- unique(gsub(":.*", "", coef.nam))

reg.store = NULL
for(cur.coef in coef.uniq){
  ind = grep(cur.coef, coef.nam)
  reg.store = rbind(reg.store,
                    coefficients(mod.lin_1mod)[ind])
}

rownames(reg.store) <- coef.uniq
colnames(reg.store) <- c('beta0sc', 'beta1sc', 'beta2sc')

# unscaled coefficients
coef.unsc <- coefs_unscaled(reg.store,
                            mean = mean(deriv_site_1mod$flower_50_doy),
                            sd = sd(deriv_site_1mod$flower_50_doy))

# Check that all together model and individual models are doing the same thing -------------
coef.unsc %>% as.data.frame() %>%
  mutate(Year = row.names(coef.unsc),
         Year = gsub(Year, pattern = "Year", replacement = "")) %>%
  left_join(deriv_site_1mod %>% select(Year, named_location, Trial_ID, Year_actual), by = "Year")

equal_coefs_check <- R2_county_data %>% select(Year, named_location, beta0_gau, beta1_gau, beta2_gau, beta0_lm, beta1_lm) %>%
  mutate(Year = as.factor(Year)) %>%
  # left_join(., coef.unsc.lm %>% as.data.frame() %>%
  #             mutate(Year = row.names(coef.unsc.lm) %>%
  #                      gsub(pattern = "Year", replacement = "")) %>%
  #             mutate(Year = as.numeric(Year))) %>%
  full_join(., coef.unsc %>% as.data.frame() %>%
              mutate(Year = row.names(coef.unsc),
                     Year = gsub(Year, pattern = "Year", replacement = "")) %>%
              left_join(deriv_site_1mod %>% select(Year, named_location, Trial_ID, Year_actual) %>% distinct(), 
                        by = "Year") %>%
              rename(Year_trialcode = Year), 
            by = c("Year" = "Year_actual", "named_location")) %>% 
  mutate(equal_coefs = ifelse(abs(beta0_gau - beta0 < 0.01)  &
                                abs(beta1_gau - beta1 < 0.01) &
                                abs(beta2_gau - beta2 < 0.01),
                              "GOOD", "FALSE")) %>%
  filter(equal_coefs == FALSE)
  # head(15)
View(equal_coefs_check)

anti_join(coef.unsc %>% as.data.frame() %>%
            mutate(Year = row.names(coef.unsc),
                   Year = gsub(Year, pattern = "Year", replacement = "")) %>%
            left_join(deriv_site_1mod %>% select(Year, named_location, Trial_ID, Year_actual) %>% distinct(), 
                      by = "Year") %>%
            rename(Year_trialcode = Year),
          R2_county_data %>% select(Year, named_location, beta0_gau, beta1_gau, beta2_gau, beta0_lm, beta1_lm) %>%
            mutate(Year = as.factor(Year)),
          by = c("Year_actual" = "Year", "named_location")
          )

opt_from_multiplemodels <- read.csv("derived_data/optimum_flowering_time_county.csv")


anti_join(coef.unsc %>% as.data.frame() %>%
            mutate(Year = row.names(coef.unsc),
                   Year = gsub(Year, pattern = "Year", replacement = "")) %>%
            left_join(deriv_site_1mod %>% select(Year, named_location, Trial_ID, Year_actual) %>% distinct(), 
                      by = "Year") %>%
            rename(Year_trialcode = Year),
          opt_from_multiplemodels %>%
            mutate(Year = as.factor(Year)),
          by = c("Year_actual" = "Year", "named_location")
)

# Calculate phenological metrics from coefficients -------------------
pheno_metrics <- pheno_calc(coef.unsc)
pheno_metrics$Year <- rownames(pheno_metrics)
# str(pheno_metrics)
pheno_metrics$Year <- gsub(pattern = "Year", replacement = "", x = as.character(pheno_metrics$Year))

# figure out if estimate is within range of data
range_flower_doy <- deriv_site_1mod %>% group_by (Year) %>%
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

