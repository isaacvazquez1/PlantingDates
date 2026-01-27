# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
############    Purpose: Trends in Climate Data
############ -    
############ -
############            By: Eliza Clark
############       Last modified : 11/20/2024
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggpubr)

figure_filepath <- "figures"
deriv_data_filepath <- "derived_data"


# Data --------------------------------------------------------------------

clim_yearly_out <- read.csv(paste0(deriv_data_filepath,"/yearly_climate_data.csv"))

optims <- read.csv(paste0(deriv_data_filepath,"/optimum_flowering_time_county.csv"))




# CLIMATE TRENDS ------------------------------------------------------------


# filter clim_yearly_out to only include locations in optims
clim_yearly_filtered_locs <- clim_yearly_out %>% 
  filter(location %in% unique(optims$location)) %>%
  distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(county_state = gsub("_", ", ", county_state))
unique(clim_yearly_filtered_locs$location)
unique(clim_yearly_filtered_locs$county_state)



# GDD trends ---------------------------------------------------------------
# Climate trends by location
# GDD anomaly trend over years
clim_loess_trend <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(GDD_mod_sitemean = ave(GDD_fixed, location, latitude, longitude, FUN = mean),
         GDD_anomaly = GDD_fixed - GDD_mod_sitemean,
         location = factor(location, levels = unique(.$location[order(.$latitude)])) ,
         county_state = factor(county_state, levels = unique(.$county_state[order(.$latitude)])), 
         ave5yr = slider::slide_dbl(GDD_anomaly, mean, .before = 2, .after = 2, .complete = F)) %>%
  ggplot(aes(x = year, y = GDD_anomaly, color = county_state)) +
  geom_jitter(alpha = 0.5, width = 0.15, height = 0) +
  geom_line(aes(y = ave5yr), linewidth = 1.5, alpha = 0.5) + 
  # geom_smooth(se = F) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(28)) +
  labs(title = "Growing Degree Days",
       y = 'GDD anomaly', x = 'Year', color = 'County (S to N)') +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        panel.grid = element_blank())
clim_loess_trend

# GDD_modave slope across years for each location
clim_slopes <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(GDD_mod_sitemean = ave(GDD_fixed, location, latitude, longitude, FUN = mean),
         GDD_anomaly = GDD_fixed - GDD_mod_sitemean) %>%
  group_by(county_state) %>%
  do(broom::tidy(lm(data = ., formula = GDD_fixed ~ year))) %>% 
  filter(term == 'year') %>% 
  mutate(positive = ifelse(estimate > 0, 'Positive', 'Negative'),
         sig = ifelse(p.value <= 0.05, 'Sig.', 'Non-sig.')) %>%
  left_join(., clim_yearly_filtered_locs %>% distinct(county_state, latitude, longitude), by = 'county_state') %>%
  mutate(county_state = factor(county_state, levels = unique(.$county_state[order(.$latitude)])) ) %>%
  ggplot(aes(x = county_state, y = estimate, color = positive)) +
  # geom_smooth(method = 'lm', se = F) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error,
                    linetype = sig), linewidth = .9) +
  geom_point(aes(shape = sig), size = 3, fill = 'white') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_linetype_manual(values = c('Sig.' = 'solid', 'Non-sig.' = 'dashed')) +
  scale_shape_manual(values = c('Sig.' = 16, 'Non-sig.' = 21)) +
  scale_color_manual(values = c('Positive' = 'palegreen4', 'Negative' = 'darkorange2')) +
  labs(y = 'GDD trend', x = 'County (S to N)', 
       color = 'Trend direction', linetype = "Trend significance") +
  guides(shape = 'none', fill = 'none') +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank())
clim_slopes

ggarrange(clim_loess_trend, clim_slopes, ncol = 1, nrow = 2, legend = 'right',
          labels = "AUTO", heights = c(1.25, 1))


# # Breakpoint Segmented Analysis
# library(segmented)
# 
# lm_car1 <- lm(GDD_fixed ~ year, data = clim_yearly_filtered_locs %>% filter(location == 'Carrington'))
# seg_car1 <- segmented(lm_car1, seg.Z = ~ year, psi = 2000)
# seg_car1$psi
# slope(seg_car1)
# fitted(seg_car1)
# 
# 
# # write function to do breakpoint segmented analysis for each location in data
# seg_slope <- function(loc, data = clim_yearly_filtered_locs){
#   lm_loc <- lm(GDD_fixed ~ year, data = data %>% filter(location == loc))
#   seg_loc <- segmented(lm_loc, seg.Z = ~ year, psi = 2000)
#   return(slope(seg_loc))
# }
# 
# seg_breakpoint <- function(loc, data = clim_yearly_filtered_locs){
#   lm_loc <- lm(GDD_fixed ~ year, data = data %>% filter(location == loc))
#   seg_loc <- segmented(lm_loc, seg.Z = ~ year, psi = 2000)
#   return(as.data.frame(seg_loc$psi))
# }
# 
# seg_slope(loc = unique(clim_yearly_filtered_locs$location)[2])
# seg_breakpoint(loc = unique(clim_yearly_filtered_locs$location)[2])
# 
# # run get slope of segments over all locations
# sapply(as.character(unique(clim_yearly_filtered_locs$location)), seg_slope)
# sapply(as.character(unique(clim_yearly_filtered_locs$location)), seg_breakpoint) %>% t()
# 

# Precip trends ---------------------------------------------------------------

# Climate trends by location
# GDD anomaly trend over years
prec_5yr_trend <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(precip_sitemean = ave(tot_precip, location, latitude, longitude, FUN = mean),
         precip_anomaly = tot_precip - precip_sitemean) %>%
  mutate(county_state = factor(county_state, levels = unique(.$county_state[order(.$latitude)])) ) %>%
  mutate(ave5yr = slider::slide_dbl(precip_anomaly, mean, .before = 2, .after = 2, .complete = F)) %>%
  ggplot(aes(x = year, y = precip_anomaly, color = county_state)) +
  geom_jitter(alpha = 0.5, width = 0.15, height = 0) +
  geom_line(aes(y = ave5yr), size = 1.5, alpha = 0.5) + 
  # geom_smooth(se = F) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(length(unique(clim_yearly_filtered_locs$county_state)))) +
  labs(title = 'Precipitation',
       y = 'Precipitation anomaly', x = 'Year', color = 'County (S to N)') +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        panel.grid = element_blank())
prec_5yr_trend

# GDD_modave slope across years for each location
prec_slopes <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(precip_sitemean = ave(tot_precip, location, latitude, longitude, FUN = mean),
         precip_anomaly = tot_precip - precip_sitemean) %>%
  group_by(county_state) %>%
  do(broom::tidy(lm(data = ., formula = tot_precip ~ year))) %>% 
  filter(term == 'year') %>% 
  mutate(positive = ifelse(estimate > 0, 'Positive', 'Negative'),
         sig = ifelse(p.value <= 0.05, 'Sig.', 'Non-sig.')) %>%
  left_join(., clim_yearly_filtered_locs %>% distinct(county_state, latitude, longitude), by = 'county_state') %>%
  mutate(county_state = factor(county_state, levels = unique(.$county_state[order(.$latitude)])) ) %>%
  ggplot(aes(x = county_state, y = estimate, color = positive)) +
  # geom_smooth(method = 'lm', se = F) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error,
                    linetype = sig), size = .9) +
  geom_point(aes(shape = sig), size = 3, fill = 'white') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_linetype_manual(values = c('Sig.' = 'solid', 'Non-sig.' = 'dashed')) +
  scale_shape_manual(values = c('Sig.' = 16, 'Non-sig.' = 21)) +
  scale_color_manual(values = c('Positive' = 'palegreen4', 'Negative' = 'darkorange2')) +
  labs(y = 'Precipitation trend', x = 'County (S to N)', 
       color = 'Trend direction', linetype = "Trend significance") +
  guides(shape = 'none', fill = 'none') +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        legend.box.margin = margin(20, 55, 20, 55))
prec_slopes

ggarrange(prec_5yr_trend, prec_slopes, ncol = 1, nrow = 2, legend = 'right',
          labels = c("C", "D"), heights = c(1.25, 1))


# Frost day trends ---------------------------------------------------------------

# Climate trends by location
# GDD anomaly trend over years
frost_5yr_trend <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(frost_sitemean = ave(frost_day, location, latitude, longitude, FUN = mean),
         frost_anomaly = frost_day - frost_sitemean) %>%
  mutate(county_state = factor(county_state, levels = unique(.$county_state[order(.$latitude)])) ) %>%
  mutate(ave5yr = slider::slide_dbl(frost_anomaly, mean, .before = 2, .after = 2, .complete = F)) %>%
  ggplot(aes(x = year, y = frost_anomaly, color = county_state)) +
  geom_jitter(alpha = 0.5, width = 0.15, height = 0) +
  geom_line(aes(y = ave5yr), size = 1.5, alpha = 0.5) + 
  # geom_smooth(se = F) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(length(unique(clim_yearly_filtered_locs$county_state)))) +
  labs(title = 'Frost day',
    y = 'Frost day anomaly', x = 'Year', color = 'County (S to N)') +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        panel.grid = element_blank())
frost_5yr_trend

# GDD_modave slope across years for each location
frost_slopes <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(frost_sitemean = ave(frost_day, location, latitude, longitude, FUN = mean),
         frost_anomaly = frost_day - frost_sitemean) %>%
  group_by(county_state) %>%
  do(broom::tidy(lm(data = ., formula = frost_day ~ year))) %>% 
  filter(term == 'year') %>% 
  mutate(positive = ifelse(estimate > 0, 'Positive', 'Negative'),
         sig = ifelse(p.value <= 0.05, 'Sig.', 'Non-sig.')) %>%
  left_join(., clim_yearly_filtered_locs %>% distinct(county_state, latitude, longitude), by = 'county_state') %>%
  mutate(county_state = factor(county_state, levels = unique(.$county_state[order(.$latitude)])) ) %>%
  ggplot(aes(x = county_state, y = estimate, color = positive)) +
  # geom_smooth(method = 'lm', se = F) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error,
                    linetype = sig), size = .9) +
  geom_point(aes(shape = sig), size = 3, fill = 'white') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_linetype_manual(values = c('Sig.' = 'solid', 'Non-sig.' = 'dashed')) +
  scale_shape_manual(values = c('Sig.' = 16, 'Non-sig.' = 21)) +
  scale_color_manual(values = c('Positive' = 'palegreen4', 'Negative' = 'darkorange2')) +
  labs(y = 'Frost day trend', x = 'County (S to N)', 
       color = 'Trend direction', linetype = "Trend significance") +
  guides(shape = 'none', fill = 'none') +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        legend.box.margin = margin(20, 55, 20, 55))
frost_slopes

ggarrange(frost_5yr_trend, frost_slopes, ncol = 1, nrow = 2, legend = 'right',
          labels = c("E", "F"), heights = c(1.25, 1))


# Combine all plots ----------------------------------------------------------
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

ggarrange(
ggarrange(clim_loess_trend, frost_5yr_trend, prec_5yr_trend,
          ncol = 3, nrow = 1, legend = 'right', legend.grob = get_legend(clim_loess_trend),
          labels = c("A", "C", "E"), #heights = c(1.25, 1),
          common.legend = T),
ggarrange(clim_slopes, frost_slopes, prec_slopes, 
          ncol = 3, nrow = 1, legend = 'right', legend.grob = get_legend(frost_slopes),
                 labels = c("B", "D", "F"), #heights = c(1.25, 1),
          common.legend = T),
nrow = 2, ncol = 1
)

ggsave(paste0(figure_filepath,"/climate_trends.png"), width = 20, height = 11, dpi = 600)

# # This makes essentially the same thing, but with a little more room around the legend
# ggarrange(
# ggarrange(clim_loess_trend, clim_slopes, 
#           ncol = 1, legend = 'none', labels = 'AUTO'),
# ggarrange(frost_5yr_trend, frost_slopes, 
#           ncol = 1, legend = 'none', labels = c("C", "D")),
# ggarrange(prec_5yr_trend, prec_slopes, 
#           ncol = 1, legend = 'none', labels = c("E", "F")),
# ggarrange(get_legend(clim_loess_trend), get_legend(clim_slopes),
#           ncol = 1),
# ncol = 4)






#. --------------------------------------------------------------------------
# Climate Maps ----------------------------------------------------------
world <- map_data("world")
states <- map_data("state")

# GDD
clim_GDD <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(GDD_mod_sitemean = ave(GDD_modave, location, latitude, longitude, FUN = mean),
         GDD_anomaly = GDD_modave - GDD_mod_sitemean) %>%
  mutate(location = factor(location, levels = unique(.$location[order(.$latitude)])) ) %>%
  distinct(location, latitude, longitude, GDD_mod_sitemean) %>%
  distinct(location, .keep_all = T) %>%
  filter(location != "Ellis") %>%
ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey95') +
  geom_point(aes(x = longitude, y = latitude, color = GDD_mod_sitemean), size = 4) +
  scale_color_gradient(low="dodgerblue3", high="red2") +
  coord_sf(xlim = c(-110, -92), ylim = c(25.5, 50), expand = F) +
  labs(x = '', y = '', color = 'County mean GDD') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())
clim_GDD


# Precip
clim_prec <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(prec_mod_sitemean = ave(tot_precip, location, latitude, longitude, FUN = mean),
         prec_anomaly = tot_precip - prec_mod_sitemean) %>%
  mutate(location = factor(location, levels = unique(.$location[order(.$latitude)])) ) %>%
  distinct(location, latitude, longitude, prec_mod_sitemean) %>%
  distinct(location, .keep_all = T) %>%
  filter(location != "Ellis") %>%
ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey95') +
  geom_point(aes(x = longitude, y = latitude, color = prec_mod_sitemean), size = 4) +
  scale_color_gradient(low="darkorange2", high="springgreen3") +
  coord_sf(xlim = c(-110, -92), ylim = c(25.5, 50), expand = F) +
  labs(x = '', y = '', color = 'County mean Precipitation') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())
clim_prec

# Frost Day
clim_frost <- clim_yearly_filtered_locs %>% 
  # distinct(latitude, longitude, year, .keep_all = T) %>%
  mutate(frost_sitemean = ave(frost_day, location, latitude, longitude, FUN = mean),
         frost_anomaly = frost_day - frost_sitemean) %>%
  mutate(location = factor(location, levels = unique(.$location[order(.$latitude)])) ) %>%
  distinct(location, latitude, longitude, frost_sitemean) %>%
  distinct(location, .keep_all = T) %>%
  filter(location != "Ellis") %>%
ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = 'black', fill = 'grey95') +
  geom_point(aes(x = longitude, y = latitude, color = frost_sitemean), size = 4) +
  scale_color_gradient(low="dodgerblue3", high="red2") +
  coord_sf(xlim = c(-110, -92), ylim = c(25.5, 50), expand = F) +
  labs(x = '', y = '', color = 'County mean \nFrost Day') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())
clim_frost

ggarrange(clim_GDD, clim_prec, clim_frost, ncol = 3, nrow = 1, legend = 'right')
# ggsave("temp_plots/ec/site_mean_clim.png", 
#        width = 12, height = 5)
