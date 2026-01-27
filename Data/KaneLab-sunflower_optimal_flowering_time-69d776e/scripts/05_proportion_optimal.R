# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
############    Purpose: How does the proportion of sites
############     with early/late/optimal timing change across
############      time and location?
############ -
############            By: Eliza Clark
############       Last modified : 1/22/2025
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggpmisc)
library(VGAM)
library(ggpubr)
library(brms)
library(tidybayes)

deriv_data_filepath <- "derived_data"
figure_filepath <- "figures"

# Data --------------------------------------------------------------------
by_county = TRUE

if(by_county == TRUE){
  optim_clim_county <- read.csv(paste0(deriv_data_filepath,"/optim_clim_county.csv"))
  path2plots = paste0(figure_filepath,"/county/")
}

if(by_county == FALSE){
  optim_clim_county <- read.csv(paste0(deriv_data_filepath,"/optim_clim.csv"))
  path2plots = paste0(figure_filepath,"/exact_location/")
}




# Multinomial regressions --------------------------------------------------
str(optim_clim_county)

# What is average latitude of each county?
# This will be used to join with other dataframes.
county_latitudes <- optim_clim_county %>% 
  group_by(group, State) %>%
  summarise(county_latitude = mean(latitude),
            county_longitude = mean(longitude))


## By Year -------------------------------------------------------------------
# Fit multinomial regression model to see how the proportion of 
# optimal, early, late, and unclear
# trials changes over time.

# data
optim_clim_county_multi <- optim_clim_county %>% 
  # select(Year, opt_relation_to_data) %>%
  filter(Year > 1979) %>%
  mutate(opt_relation_to_data = factor(opt_relation_to_data, levels = c("opt_within", "opt_after", "opt_before", "unclear"))) %>%
  add_count(group)

# Fit model
year_multi_mod <- vglm(opt_relation_to_data ~ Year, family = multinomial(refLevel = 1, parallel = FALSE), 
                       data = optim_clim_county_multi)
summary(year_multi_mod)
anova(year_multi_mod)
# Likelihood ratio test
null_multi_mod <- vglm(opt_relation_to_data ~ 1, multinomial(refLevel = 1), 
                       data = optim_clim_county_multi)
lrtest(year_multi_mod, null_multi_mod)

# Predicted probabilities
year_multi_predict <- data.frame(Year = 1980:2023, 
                                 predict(year_multi_mod, newdata = data.frame(Year = 1980:2023),
                                         type = 'response'))


#########################
year_multi_predict %>% 
  mutate(across(-Year, ~ c(NA, diff(.) / diff(Year)[1]), .names = "slope_{.col}")) %>% 
  select(Year, starts_with("slope_")) %>%
  pivot_longer(-Year, names_to = "Category", values_to = "Slope")

###
# Extract coefficients and standard errors from summary
coef_summary_vg <- summary(year_multi_mod)@coef3  # matrix of Estimate, Std. Error, z value

# Convert to data frame
ci_df_vg <- as.data.frame(coef_summary_vg)

# Keep only the Year slopes (not intercepts)
ci_df_vg <- ci_df_vg %>%
  rownames_to_column(var = "Term") %>%
  filter(grepl("^Year", Term)) %>%   # keep only Year coefficients
  rename(
    Estimate = Estimate,
    StdError = `Std. Error`
  )

# Calculate 95% confidence intervals (on log-odds scale)
ci_df_vg <- ci_df_vg %>%
  mutate(
    CI_lower = Estimate - 1.96 * StdError,
    CI_upper = Estimate + 1.96 * StdError,
    
    # Also compute odds ratio and its CI
    OR = exp(Estimate),
    OR_lower = exp(CI_lower),
    OR_upper = exp(CI_upper)
  )
##################


# create dataframe that totals the number of sites that are 
# optimal, early, late, or unclear by year.
# Then pivot longer.
proportion_optima <- optim_clim_county %>% filter(Year > 1979) %>%
  group_by(Year) %>%
  summarize(
    total_sites = n(),
    # sites_with_opt = sum(outside_data == 'in_range'),
    count_suitable = sum(opt_relation_to_data == 'opt_within', na.rm = T),
    count_early = sum(opt_relation_to_data == 'opt_after', na.rm = T),
    count_late = sum(opt_relation_to_data == 'opt_before', na.rm = T),
    count_unclear = sum(opt_relation_to_data == 'unclear', na.rm = T),
    prop_suitable = count_suitable/total_sites,
    prop_early = count_early/total_sites,
    prop_late = count_late/total_sites,
    prop_unclear = count_unclear/total_sites
  ) %>%
  pivot_longer(cols = starts_with(c("prop", "count")),
               names_to = c(".value", "flowering_was" ),
               names_sep = '_')
proportion_optima$flowering_was <- factor(proportion_optima$flowering_was, 
                                          levels = c("suitable", "early", "late", "unclear"))

# Plot
year_mult <- ggplot(data = year_multi_predict %>%
        pivot_longer(cols = -Year,
                     names_to = 'flowering_was',
                     values_to = 'prob') %>%
         mutate(flowering_was = factor(flowering_was, 
                                      levels = c("opt_within", "opt_after", "opt_before", "unclear"),
                                      labels = c("suitable", "early", "late", "unclear"))), 
       aes(x = Year, y = prob, color = flowering_was)) +
  geom_line(linewidth = 1.4, linetype = 2) +
  # add points
  geom_jitter(data = proportion_optima, aes(x = Year, y = prop, color = flowering_was), 
              alpha = 0.5, size = 2) +
  scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "Predicted probability", color = "Flowering was:", x = "Year") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'right',
        panel.grid = element_blank())
year_mult


## By Year BRMS ----------------------------------------------
yr_mod_brm <- brm(
  bf(opt_relation_to_data ~ Year + (1 | county_state)),
  data = optim_clim_county_multi,
  family = categorical(refcat = "unclear"),
  prior = c(
    prior(normal(0, 5), class = b, dpar = muoptafter),
    prior(normal(0, 5), class = b, dpar = muoptbefore),
    prior(normal(0, 5), class = b, dpar = muoptwithin),
    prior(exponential(1), class = sd, dpar = muoptafter),
    prior(exponential(1), class = sd, dpar = muoptbefore),
    prior(exponential(1), class = sd, dpar = muoptwithin)
  ),
  adapt_delta = 0.9,
  chains = 4, cores = 4, iter = 2000, seed = 5,
  file = "derived_data/yr_mod_brm" 
)
summary(yr_mod_brm)

# extract conditional effects
ce_yr <- brms::conditional_effects(yr_mod_brm, categorical = TRUE,
                                      re_formula = NA)

ce_yr$`Year:cats__` |>
  ggplot() +
  geom_jitter(data = proportion_optima %>%
                mutate(flowering_was = factor(flowering_was,
                                              labels = c("opt_within", "opt_after", "opt_before", "unclear"),
                                              levels = c("suitable", "early", "late", "unclear"))), 
              aes(x = Year, y = prop, color = flowering_was), 
              alpha = 0.2, size = 2) +
  geom_ribbon(
    aes(Year, y = estimate__,
        ymin = lower__, ymax = upper__,
        group = cats__),
    # color = cats__, fill = cats__),
    alpha = 0.1
  ) +
  geom_line(aes(Year, estimate__, color = cats__),
            linetype = 2, linewidth = 1.2) +
  scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50'),
                     labels = c("suitable", "early", "late", "unclear")) +
  # scale_fill_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "Predicted probability", color = "Flowering was:", x = "GDD site mean") +
  # guides(color = 'none') +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 12)
  ) +
  coord_cartesian(ylim = c(0, 1))
ggsave(paste0(path2plots, "year_mult.png"),
       width = 8, height = 7)

## By Location ---------------------------------------------------------------
# Fit multinomial regression model to see how the proportion of 
# optimal, early, late, and unclear
# trials changes over Location.

# Fit model
county_multi_mod <- vglm(opt_relation_to_data ~ group, multinomial(refLevel = 1), 
                       data = optim_clim_county_multi)
summary(county_multi_mod)
anova(county_multi_mod, type = 3)
# Likelihood ratio test
# same null model as above
lrtest(county_multi_mod, null_multi_mod)
# Predicted probabilities
county_multi_predict <- data.frame(group = unique(optim_clim_county_multi$group), 
                                 predict(county_multi_mod, 
                                         newdata = data.frame(group = unique(optim_clim_county_multi$group)),
                                         type = 'response'))
# sample sizes
optim_clim_county_multi %>%
  group_by(group) %>%
  summarize(n = n())
# Graph
county_mult <- ggplot(data = county_multi_predict %>%
                        # get data in correct format
         pivot_longer(cols = -group,
                      names_to = 'flowering_was',
                      values_to = 'prob') %>%
         full_join(., county_latitudes, by = 'group') %>%
         mutate(flowering_was = factor(flowering_was, 
                                       levels = c("opt_within", "opt_after", "opt_before", "unclear"),
                                       labels = c("optimal", "early", "late", "unclear")),
                group = factor(group, levels = county_latitudes %>%
                                 arrange(county_latitude) %>% pull(group) ))
       , 
       aes(x = group, y = prob, fill = flowering_was)) +
  # add columns
  geom_col(position = position_stack(reverse = TRUE)) +
  # add sample size annotation
  geom_text(data = optim_clim_county_multi %>%
              group_by(group) %>%
              summarize(n = n()),
            aes(x = group, y = 0.1, label = n), inherit.aes = FALSE,
            size = 4) +
  scale_fill_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "Predicted probability", fill = "Flowering was:", x = "County (S to N)") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
        )
county_mult

ggarrange(year_mult, county_mult, ncol = 2, labels = 'AUTO', widths = c(1, 1.8),
          common.legend = T, legend = 'top')
# ggsave(paste0(path2plots, "year_county_mult.png"), 
#        width = 12, height = 7)

# look at whether there is a trend with sample size
ggplot(data = county_multi_predict %>%
         # get data in correct format
         pivot_longer(cols = -group,
                      names_to = 'flowering_was',
                      values_to = 'prob') %>%
         full_join(., county_latitudes, by = 'group') %>%
         full_join(., optim_clim_county_multi %>% group_by(group) %>% summarize(n = n()),
                   by = 'group') %>%
         mutate(flowering_was = factor(flowering_was, 
                                       levels = c("opt_within", "opt_after", "opt_before", "unclear"),
                                       labels = c("optimal", "early", "late", "unclear")),
                group = factor(group, levels = optim_clim_county_multi %>% group_by(group) %>% summarize(n = n()) %>%
                                 arrange(n) %>% pull(group) ))
       , 
       aes(x = group, y = prob, fill = flowering_was)) +
  # add columns
  geom_col(position = position_stack(reverse = TRUE)) +
  # add sample size annotation
  geom_text(data = optim_clim_county_multi %>%
              group_by(group) %>%
              summarize(n = n()),
            aes(x = group, y = 0.1, label = n), inherit.aes = FALSE,
            size = 4) +
  scale_fill_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "Predicted probability", fill = "Flowering was:", x = "County (S to N)") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )


## By Climate ---------------------------------------------------------------
# Fit multinomial regression model to see how the proportion of 
# optimal, early, late, and unclear
# trials changes over climate variables


# Fit model
levels(optim_clim_county_multi$opt_relation_to_data)
GDD_multi_mod <- vglm(opt_relation_to_data ~ GDD_mod_sitemean + GDD_anomaly + 
                        tot_precip_sitemean + precip_anomaly + frost_anomaly, 
                      multinomial(refLevel = 4), 
                      data = optim_clim_county_multi %>% filter(!is.na(GDD_anomaly)))
# GDD_multi_mod_simple <- vglm(opt_relation_to_data ~ GDD_mod_sitemean + precip_anomaly, 
#                       multinomial(refLevel = 4), 
#                       data = optim_clim_county_multi %>% filter(!is.na(GDD_anomaly)))
# AIC(GDD_multi_mod_simple)

summary(GDD_multi_mod)
anova(GDD_multi_mod, type = 3)
AIC(GDD_multi_mod)
# Likelihood ratio test
null_multi_mod2 <- vglm(opt_relation_to_data ~ 1, multinomial(refLevel = 4), 
                       data = optim_clim_county_multi %>% filter(!is.na(GDD_anomaly)))
AIC(null_multi_mod2)
lrtest(GDD_multi_mod, null_multi_mod2)

null_multi_mod_g <- vglm(opt_relation_to_data ~ GDD_anomaly + tot_precip_sitemean + precip_anomaly + frost_anomaly, 
                      multinomial(refLevel = 4), 
                      data = optim_clim_county_multi %>% filter(!is.na(GDD_anomaly)))
lrtest(GDD_multi_mod, null_multi_mod_g)

null_multi_mod_f <- vglm(opt_relation_to_data ~ GDD_mod_sitemean + GDD_anomaly + tot_precip_sitemean + precip_anomaly, 
                         multinomial(refLevel = 4), 
                         data = optim_clim_county_multi %>% filter(!is.na(GDD_anomaly)))
lrtest(GDD_multi_mod, null_multi_mod_f)



# Predicted probabilities
## find min/max of climate params to predict across
clim_summary <- optim_clim_county_multi %>%
  summarise(
    min_GDD = min(GDD_mod_sitemean, na.rm = T),
    mean_GDD = mean(GDD_mod_sitemean, na.rm = T),
    max_GDD = max(GDD_mod_sitemean, na.rm = T),
    min_precip_anom = min(precip_anomaly, na.rm = T),
    mean_precip_anom = mean(precip_anomaly, na.rm = T),
    max_precip_anom = max(precip_anomaly, na.rm = T),
    precip_anom_025 = quantile(precip_anomaly, 0.05, na.rm = T),
    precip_anom_975 = quantile(precip_anomaly, 0.95, na.rm = T),
    min_frost = min(frost_anomaly, na.rm = T),
    mean_frost = mean(frost_anomaly, na.rm = T),
    max_frost = max(frost_anomaly, na.rm = T),
    frost_05 = quantile(frost_anomaly, 0.05, na.rm = T),
    frost_95 = quantile(frost_anomaly, 0.95, na.rm = T),
    mean_GDD_anom = mean(GDD_anomaly, na.rm = T),
    mean_precip = mean(tot_precip_sitemean, na.rm = T)
  )

## new data to predict across, with GDD_sitemean ranging from min to max
# and all other climate variables at their mean
GDD_multi_predict_newdat <- data.frame(GDD_mod_sitemean = c(seq(clim_summary$min_GDD, clim_summary$max_GDD, length.out = 30),
                                                       rep(clim_summary$mean_GDD, 30*2)),
                                       precip_anomaly = c(rep(clim_summary$mean_precip_anom, 30),
                                                          seq(clim_summary$min_precip_anom, clim_summary$max_precip_anom, length.out = 30),
                                                          rep(clim_summary$mean_precip_anom, 30)),
                                       frost_anomaly = c(rep(clim_summary$mean_frost, 30*2),
                                                         seq(clim_summary$min_frost, clim_summary$max_frost, length.out = 30)), 
                                       GDD_anomaly = rep(clim_summary$mean_GDD_anom, 30*3),
                                       tot_precip_sitemean = rep(clim_summary$mean_precip, 30*3),
                                       var = rep(c("GDD_mod_sitemean", "precip_anomaly", "frost_anomaly"), each = 30))

## predict
GDD_multi_predict <- data.frame(GDD_multi_predict_newdat, 
                                predict(GDD_multi_mod, newdata = GDD_multi_predict_newdat,
                                        type = 'response')) %>%
  pivot_longer(cols = -c(GDD_mod_sitemean, precip_anomaly, frost_anomaly, GDD_anomaly, tot_precip_sitemean, var),
               names_to = 'flowering_was',
               values_to = 'prob') %>%
  mutate(flowering_was = factor(flowering_was, 
                                levels = c("opt_within", "opt_after", "opt_before", "unclear"),
                                labels = c("suitable", "early", "late", "unclear")))

ggarrange(
  GDD_multi_predict %>% filter(var == 'GDD_mod_sitemean') %>%
    ggplot() +
    geom_line(aes(x = GDD_mod_sitemean, y = prob, color = flowering_was),
              linewidth = 1.4, linetype = 2) +
    scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
    labs(y = "Predicted probability", color = "Flowering was:", x = "GDD site mean") +
    theme_bw(base_size = 16) +
    theme(legend.position = 'right',
          panel.grid = element_blank()),
  # GDD_multi_predict %>% filter(var == 'precip_anomaly') %>%
  #   ggplot() +
  #   geom_line(aes(x = precip_anomaly, y = prob, color = flowering_was),
  #             linewidth = 1.4, linetype = 2) +
  #   scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  #   labs(y = "Predicted probability", color = "Flowering was:", x = "Precip. anomaly") +
  #   theme_bw(base_size = 16) +
  #   theme(legend.position = 'right',
  #         panel.grid = element_blank()),
  GDD_multi_predict %>% filter(var == 'frost_anomaly') %>%
    ggplot() +
    geom_line(aes(x = frost_anomaly, y = prob, color = flowering_was),
              linewidth = 1.4, linetype = 2) +
    scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
    labs(y = "Predicted probability", color = "Flowering was:", x = "First frost anomaly") +
    theme_bw(base_size = 16) +
    theme(legend.position = 'right',
          panel.grid = element_blank()),
  ncol = 2, labels = 'AUTO', common.legend = T)




## different package for climate multinomial regression to get error bars on plot
library(nnet)

test <- multinom(opt_relation_to_data ~ GDD_mod_sitemean + GDD_anomaly + 
                   tot_precip_sitemean + precip_anomaly + frost_anomaly, 
                 data = optim_clim_county_multi %>% filter(!is.na(GDD_anomaly)))

summary(test)
predict(test, newdata = GDD_multi_predict_newdat, "probs")
fit.eff <- effects::Effect(c("GDD_mod_sitemean"), test, xlevels = 50)
plot(fit.eff)
gddmulti_plot <- data.frame(fit.eff$model.matrix, fit.eff$prob, fit.eff$lower.prob, 
           fit.eff$upper.prob) %>% 
  pivot_longer(cols = -c(X.Intercept., GDD_mod_sitemean, precip_anomaly, frost_anomaly, GDD_anomaly, tot_precip_sitemean), 
               names_to = c(".value", "flowering_was"),
             names_pattern = "(.*).(opt_within|opt_after|opt_before|unclear)") %>%
  mutate(flowering_was = factor(flowering_was, 
                                levels = c("opt_within", "opt_after", "opt_before", "unclear"),
                                labels = c("suitable", "early", "late", "unclear"))) %>%
  ggplot() +
  geom_ribbon(aes(x = GDD_mod_sitemean, ymin = L.prob, ymax = U.prob,
                  group = flowering_was), alpha = 0.1) +
  geom_line(aes(x = GDD_mod_sitemean, y = prob, color = flowering_was),
            linewidth = 1.4, linetype = 2) +
  scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "Predicted probability", color = "Flowering was:", x = "GDD site mean") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'right',
        panel.grid = element_blank())

fit.eff2 <- effects::Effect(c("frost_anomaly"), test, xlevels = 50)
plot(fit.eff2)
frostmulti_plot <- data.frame(fit.eff2$model.matrix, fit.eff2$prob, fit.eff2$lower.prob, 
           fit.eff2$upper.prob) %>% 
  pivot_longer(cols = -c(X.Intercept., GDD_mod_sitemean, precip_anomaly, frost_anomaly, GDD_anomaly, tot_precip_sitemean), 
               names_to = c(".value", "flowering_was"),
               names_pattern = "(.*).(opt_within|opt_after|opt_before|unclear)") %>%
  mutate(flowering_was = factor(flowering_was, 
                                levels = c("opt_within", "opt_after", "opt_before", "unclear"),
                                labels = c("suitable", "early", "late", "unclear"))) %>%
  ggplot() +
  geom_ribbon(aes(x = frost_anomaly, ymin = L.prob, ymax = U.prob,
                  group = flowering_was), alpha = 0.1) +
  geom_line(aes(x = frost_anomaly, y = prob, color = flowering_was),
            linewidth = 1.4, linetype = 2) +
  scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "Predicted probability", color = "Flowering was:", x = "First frost anomaly") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'right',
        panel.grid = element_blank())

ggarrange(gddmulti_plot, frostmulti_plot,
  ncol = 2, labels = 'AUTO', common.legend = T)







## Plot across GDD with 3 levels of frost anomaly

GDD_multi_predict_newdat_frostxgdd <- data.frame(GDD_mod_sitemean = c(rep(seq(clim_summary$min_GDD, clim_summary$max_GDD, length.out = 30),3)),
                                       frost_anomaly = c(rep(clim_summary$mean_frost, 30),
                                                         rep(clim_summary$frost_05, 30),
                                                         rep(clim_summary$frost_95, 30)), 
                                       precip_anomaly = c(rep(clim_summary$mean_precip_anom, 30*3)),
                                       GDD_anomaly = rep(clim_summary$mean_GDD_anom, 30*3),
                                       tot_precip_sitemean = rep(clim_summary$mean_precip, 30*3),
                                       var = rep(c("mean_frost", "early_frost", "late_frost"), each = 30))
## predict
predict_frostxgdd <- data.frame(GDD_multi_predict_newdat_frostxgdd, 
                                predict(GDD_multi_mod, newdata = GDD_multi_predict_newdat_frostxgdd,
                                        type = 'response')) %>%
  pivot_longer(cols = -c(GDD_mod_sitemean, precip_anomaly, frost_anomaly, GDD_anomaly, tot_precip_sitemean, var),
               names_to = 'flowering_was',
               values_to = 'prob') %>%
  mutate(flowering_was = factor(flowering_was, 
                                levels = c("opt_within", "opt_after", "opt_before", "unclear"),
                                labels = c("suitable", "early", "late", "unclear")),
         frost_anom_factor = factor(frost_anomaly, 
                                    levels = c(clim_summary$frost_05, clim_summary$mean_frost, clim_summary$frost_95),
                                    labels = c("Early frost", "Average frost", "Late frost")))


ggplot(data = predict_frostxgdd,
       aes(x = GDD_mod_sitemean, y = prob, color = flowering_was)) +
  geom_line(linewidth = 1.4, linetype = 2) +
  facet_wrap(~frost_anom_factor) +
  scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "Predicted probability", color = "Flowering was:", x = "GDD site mean") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'top',
        panel.grid = element_blank())


## By Climate BRMS ----------------------------------------------------------


multi_mod_brm <- brm(
  bf(opt_relation_to_data ~ GDD_mod_sitemean + GDD_anomaly + 
       tot_precip_sitemean + precip_anomaly + frost_anomaly + (1 | county_state)),
  data = optim_clim_county_multi %>% filter(!is.na(GDD_anomaly)),
  family = categorical(refcat = "opt_within"),
  prior = c(
    prior(normal(0, 5), class = b, dpar = muoptafter),
    prior(normal(0, 5), class = b, dpar = muoptbefore),
    prior(normal(0, 5), class = b, dpar = muunclear),
    prior(exponential(1), class = sd, dpar = muoptafter),
    prior(exponential(1), class = sd, dpar = muoptbefore),
    prior(exponential(1), class = sd, dpar = muunclear)
  ),
  chains = 4, cores = 4, iter = 2000, seed = 5,
  file = "derived_data/multi_mod_brm",
  file_refit = 'on_change'
)

# library(bayesplot)
# mcmc_trace(multi_mod_brm)
# library(shinystan)
# launch_shinystan(multi_mod_brm)

prior_summary(multi_mod_brm)
summary(multi_mod_brm)


# NOTE: Here's some code for plotting the conditional/marginal effects of the 
# multinomial you fit. re_formula = NA marginalizes over the group effects.
# This is probably the best way to start thinking about what this models says
# within it's set of assumptions. Plotting bayesian multinomials is a tad tricky!
# 

# Plotting 

# extract conditional effects
ce_multi <- brms::conditional_effects(multi_mod_brm, categorical = TRUE,
                                      re_formula = NA)

### ---- GDD site mean ----
gddmean_cond <- ce_multi$`GDD_mod_sitemean:cats__` |>
  ggplot() +
  geom_ribbon(
    aes(GDD_mod_sitemean, y = estimate__,
        ymin = lower__, ymax = upper__,
        group = cats__),
        # color = cats__, fill = cats__),
    alpha = 0.1
  ) +
  geom_line(aes(GDD_mod_sitemean, estimate__, color = cats__), 
            linetype = 2, linewidth = 1.2) +
  scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50'),
                     labels = c("suitable", "early", "late", "unclear")) +
  # scale_fill_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "Predicted probability", color = "Flowering was:   ", x = "GDD site mean") + # need spaces for saving fig
  # guides(color = 'none') +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1))
gddmean_cond

### ---- frost anomaly ----
frostanom_cond <- ce_multi$`frost_anomaly:cats__` |>
  ggplot() +
  geom_ribbon(
    aes(frost_anomaly, y = estimate__,
        ymin = lower__, ymax = upper__,
        group = cats__),
    # color = cats__, fill = cats__),
    alpha = 0.1
  ) +
  geom_line(aes(frost_anomaly, estimate__, color = cats__), 
            linetype = 2, linewidth = 1.2) +
  scale_color_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50'),
                     labels = c("suitable", "early", "late", "unclear")) +
  # scale_fill_manual(values = c('forestgreen', 'darkslategray3', 'orange', 'grey50')) +
  labs(y = "", color = "Flowering was:   ", x = "First frost anomaly") + # need spaces for saving fig
  guides(fill = 'none') +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 1))
frostanom_cond

# ---- arrange them
ggpubr::ggarrange(
  gddmean_cond,
  frostanom_cond, ncol = 2, nrow = 1,
  common.legend = TRUE,
  labels = "AUTO"
)

ggsave(paste0(path2plots, "climate_mult_full.png"), 
       width = 12*.75, height = 7*.75, dpi = 600)
ggsave(paste0(path2plots, "climate_mult_full.pdf"), 
       width = 12*.75, height = 7*.75, dpi = 600)


# Pie chart -----------------------------------------------------------------
optim_clim_county %>%
  # filter(Year > 1979) %>%
  summarize(
    total_sites = n(),
    count_suitable = sum(opt_relation_to_data == 'opt_within', na.rm = T),
    count_early = sum(opt_relation_to_data == 'opt_after', na.rm = T),
    count_late = sum(opt_relation_to_data == 'opt_before', na.rm = T),
    count_unclear = sum(opt_relation_to_data == 'unclear', na.rm = T)
  ) %>%
  pivot_longer(cols = starts_with('count'),
               names_to = 'flowering_was',
               names_prefix = 'count_') %>%
  arrange(desc(flowering_was)) %>%
  mutate(prop = value/total_sites * 100,
         ypos = cumsum(prop) - 0.5 * prop, 
         # flowering_was = factor(flowering_was, levels = c("optimal", "early", "late", "unclear"))
  ) %>%
  ggplot(aes(x = "", y = prop, fill = flowering_was)) +
  geom_bar(width = 1, stat = 'identity') +
  geom_text(aes(y = ypos, label = paste0(round(prop, digits = 0),"%")), color = "white", size = 5) +
  coord_polar("y") +
  scale_fill_manual(values = c( 'darkslategray3', 'orange', 'forestgreen','grey50')) +
  labs(fill = "Flowering was:") +
  theme_void() +
  theme(legend.position = 'bottom')


