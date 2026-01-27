# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
####
#### Simple map of points used in optimal flowering time analyses
#### By : Eliza Clark
#### March 2025
####
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


# Packages --------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)

figure_filepath <- "figures"
deriv_data_filepath <- "derived_data"

# Data ------------------------------------------------------------------------


# All exact locations used in optimal flowering time analyses
R2_indata_exact <- read.csv(paste0(deriv_data_filepath,"/R2_indata_exact.csv"))


# All data with county grouping 
R2_indata <- read.csv(paste0(deriv_data_filepath,"/R2_indata.csv"))


# climate data
clim_yearly_out <- read.csv(paste0(deriv_data_filepath,"/yearly_climate_data.csv"))


# Map data 
world <- map_data("world")
states <- map_data("state")
counties <- map_data('county')


# Map for manuscript ---------------------------------------------------------
# grouped by county

# find number of years/trials in each county
data_map_c <- R2_indata %>%
  distinct(named_location, Location, county_state, Year, State, lat, lon) %>%
  group_by(county_state) %>%
  mutate(n_years = n()) %>% 
  summarise(lat = mean(lat), 
            lon = mean(lon),
            n_years = mean(n_years))


# climate data
#  GDD averages
GDD_county_aves <- clim_yearly_out %>% 
  filter(county_state %in% unique(data_map_c$county_state)) %>%
  distinct(latitude, longitude, year, .keep_all = T) %>% 
  mutate(GDD_fix_sitemean = ave(GDD_fixed, county_state, FUN = mean),
         GDD_fix_anomaly = GDD_fixed - GDD_fix_sitemean) %>%
  distinct(county_state, GDD_fix_sitemean) %>% 
  arrange(county_state)

data_map_c <- data_map_c %>% 
  left_join(., GDD_county_aves, by = c("county_state" = "county_state"))
data_map_c


main_map <- ggplot() +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group), color = 'grey80', fill = 'grey95') +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = 'black', alpha = 0) +
  geom_point(data = data_map_c, aes(x = lon, y = lat, size = n_years, color = GDD_fix_sitemean)) +
  scico::scale_color_scico(palette = "roma", direction = -1, na.value = 'grey75') +
  coord_sf(xlim = c(-110, -92), ylim = c(25.5, 50), expand = F) +
  labs(color = 'GDD', size = 'Years\nof data') +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank())
main_map
inset_map <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), color = 'black', fill = 'white', alpha = 1) +
  geom_rect(aes(xmin = -110, xmax = -92, ymin = 25.5, ymax = 50), color = 'black', fill = NA, linewidth = 1) +
  coord_sf(xlim = c(-128, -65), ylim = c(18, 57), expand = F) +
  theme_void() +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.5),
    plot.background = element_rect(fill = "grey85")
  )
inset_map

ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map,
            height = 0.2,
            x = -0.23,
            y = 0.08
  )

ggsave(paste0(figure_filepath,"/county/ugly_map.png"), width = 6.2, height = 6, dpi = 500)



# Maps to compare which sites are included with exact and county grouping
# With Exact grouping -------------------------------------------------------


optim_data_map <- R2_indata_exact %>% 
  distinct(named_location, Year, lat, lon) %>%
  group_by(lat, lon) %>%
  mutate(n_years = n()) %>% 

    summarise(lat = mean(lat), 
              lon = mean(lon),
              n_years = mean(n_years)) %>%
  left_join(., R2_indata_exact %>% distinct(Location, named_location, lat, lon))
optim_data_map


exact_loc_map <- ggplot() +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group), color = 'grey80', fill = 'grey95') +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = 'black', alpha = 0) +
  # geom_label(data = optim_data_map, aes(x = lon, y = lat, label = named_location), size = 2, color = 'black', fontface = 'bold', nudge_y = 0.2) +
  geom_point(data = optim_data_map, aes(x = lon, y = lat, size = n_years), color = "black") +
  coord_sf(xlim = c(-110, -92), ylim = c(25.5, 50), expand = F) +
  labs(x = '', y = '', color = '', size = 'Years\nof data',  title = "Exact Location Grouping") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())
exact_loc_map


# With County grouping -------------------------------------------------------

county_map <- ggplot() +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group), color = 'grey80', fill = 'grey95') +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = 'black', alpha = 0) +
  # geom_label(data = data_map_c, aes(x = lon, y = lat, label = county_state), size = 2, color = 'black', fontface = 'bold', nudge_y = 0.2) +
  geom_point(data = data_map_c, aes(x = lon, y = lat, size = n_years), color = "black") +
  coord_sf(xlim = c(-110, -92), ylim = c(25.5, 50), expand = F) +
  labs(x = '', y = '', color = '', size = 'Years\nof data', title = "County Grouping") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())
county_map

# comparison of exact and county grouping locations included
ggarrange(exact_loc_map, county_map, ncol = 2, nrow = 1)
# ggsave("temp_plots/opt_manuscript_figures/ugly_map_county_vs_exact.png", width = 18, height = 12, dpi = 500)





