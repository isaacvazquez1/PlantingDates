# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
############    Purpose: GRIN Accession data summary stats
############-     
############ -
############            By: Eliza Clark 
############       Last modified : 11/20/2024
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggpubr)


figure_filepath <- "figures"


# Data ----------------------------------------------------------------

# Download GRIN data from
# https://npgsweb.ars-grin.gov/gringlobal/descriptordetail?id=79038
grin_data <- read.csv("data/grin_dayflow.csv")

deriv_data <- read.csv("data/sunflower_data_simple_v1.csv")
deriv_data$Irrigated <- as.factor(deriv_data$Irrigated)
deriv_data$Oil_Confection <- as.factor(deriv_data$Oil_Confection)
deriv_data$Location <- as.factor(deriv_data$Location)
deriv_data$Unif_Name <- as.factor(deriv_data$Unif_Name)
deriv_data$flower_50_doy <- as.integer(deriv_data$flower_50_doy)
deriv_data$Year <- as.factor(deriv_data$Year)



# Histograms of genetic variation in days to flowering ----------------------


# GRIN DATA

# filter out non-Helianthus annuus species
grin_data <- grin_data %>% 
  filter(taxon == "Helianthus annuus")

# histogram of days to flowering
ggplot(grin_data, aes(x = observation_value)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Days to Flowering in Accessions\nof Helianthus annuus",
       x = "Days to Flowering",
       y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank())
  # theme_minimal()

grin_data %>% summarise(prop_55 = sum(observation_value < 55) / n())

# summarise grin data for average value for each accession
grin_data_summarised <- grin_data %>% group_by(accession_id) %>%
  summarise(mean_obs_val = mean(observation_value, na.rm = T),
            sd_obs_val = sd(observation_value, na.rm = T))
grin_data_summarised %>%
  ggplot(aes(x = mean_obs_val)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Days to Flowering in Accessions\nof Helianthus annuus",
       x = "Days to Flowering",
       y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank())
grin_data_summarised %>% summarise(prop_58 = sum(mean_obs_val < 58) / n())


# DERIV_DATA
# average number of days to flowering for each genotype
deriv_data_summarised <- deriv_data %>%
  filter(!is.na(flower_50_doy) &
           !is.na(yield_lb_acre) &
           !grepl('irrigated|NO|late|recrop|additional', Location)) %>%
  group_by(Unif_Name) %>%
  summarise(mean_flower_dpp = mean(flower_50pct_days_past_planting, na.rm = T),
            sd_flower_dpp = sd(flower_50pct_days_past_planting, na.rm = T))

deriv_data_summarised %>%
  ggplot(aes(x = mean_flower_dpp)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Average Days to Flowering in Commerical Trials",
       x = "Days to Flowering",
       y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank())
deriv_data_summarised %>% summarise(prop_58 = sum(mean_flower_dpp < 58, na.rm = T) / n())
# only 3% of genotypes had an average flowering date below 58 days (which is projected median
# optimal flowering time in the future)



deriv_data %>%
  filter(!is.na(flower_50_doy) &
           !is.na(yield_lb_acre) &
           !grepl('irrigated|NO|late|recrop|additional', Location)) %>%
  ggplot(aes(x = flower_50pct_days_past_planting)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Days to Flowering in Commerical Trials",
       x = "Days to Flowering",
       y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank())
deriv_data %>% summarise(prop_58 = sum(flower_50pct_days_past_planting < 58, na.rm = T) / n())
# Only 8% of records had a flowering date below 58 days (which is projected 
# median optimal flowering time in the future)
deriv_data %>% filter(flower_50pct_days_past_planting <= 58) %>% 
  distinct(Unif_Name) %>% nrow() / (deriv_data %>% nrow())

# Both datasets
# average number of days to flowering for each genotype
deriv_data %>%
  filter(!is.na(flower_50_doy) &
           !is.na(yield_lb_acre) &
           !grepl('irrigated|NO|late|recrop', Location)) %>%
  group_by(Unif_Name) %>%
  summarise(mean_flower_dpp = mean(flower_50pct_days_past_planting, na.rm = T),
            sd_flower_dpp = sd(flower_50pct_days_past_planting, na.rm = T)) %>%
  ggplot() +
  geom_histogram(aes(x = mean_flower_dpp, fill = "Commercial varieties"),
               alpha = 0.3, binwidth = 2) +
  geom_histogram(data = grin_data_summarised, aes(x = mean_obs_val, fill = "USDA accessions"),
               alpha = 0.3, binwidth = 2) +
  # geom_density(aes(x = mean_flower_dpp, fill = "Commercial Varieties"),
  #                 alpha = 0.3) +
  # geom_density(data = grin_data, aes(x = observation_value, fill = "USDA Accessions"),
  #              alpha = 0.3) +
  geom_vline(xintercept = c(64), linetype = 1, linewidth = 1) + # current median optimal days to flowering
  geom_vline(xintercept = c(58), linetype = 2, linewidth = 1) + # future median opt 
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  annotate("text", x = 66, y = 300, label = "Current median\noptimum", hjust = 0, size = 5) +
  annotate("text", x = 56, y = 300, label = "Future median\noptimum", hjust = 1, size = 5) +
  # annotate('segment', x = 76, y = 300, yend = 300, xend = 65.5, linewidth = 0.75,
    # arrow = arrow(length = unit(0.4, 'cm'))) +
  labs(x = "Days to flowering",
       y = "Number of genotypes") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.75, 0.4),
        legend.title = element_blank()
  )

# ggsave("temp_plots/ec/accessions_2.png", width = 8, height = 6, dpi = 300)

ggsave(paste0(figure_filepath,"/variation_DPP.png"), 
       width = 8, height = 6, dpi = 300)
ggsave(paste0(figure_filepath,"/variation_DPP.pdf"), 
       width = 8, height = 6, dpi = 600)

# Both datasets
# no averaging of genotypes
deriv_data %>%
  filter(!is.na(flower_50_doy) &
           !is.na(yield_lb_acre) &
           !grepl('irrigated|NO|late|recrop', Location)) %>%
  ggplot() +
  geom_histogram(aes(x = flower_50pct_days_past_planting, fill = "Commercial varieties"),
                 alpha = 0.3, binwidth = 2) +
  geom_histogram(data = grin_data, aes(x = observation_value, fill = "USDA accessions"),
                 alpha = 0.3, binwidth = 2) +
  # geom_density(aes(x = mean_flower_dpp, fill = "Commercial Varieties"),
  #                 alpha = 0.3) +
  # geom_density(data = grin_data, aes(x = observation_value, fill = "USDA Accessions"),
  #              alpha = 0.3) +
  geom_vline(xintercept = c(65), linetype = 1, linewidth = 1) + # current median optimal days to flowering
  geom_vline(xintercept = c(55), linetype = 2, linewidth = 1) + # future median opt 
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  annotate("text", x = 66, y = 1000, label = "Current median\noptimum", hjust = 0, size = 5) +
  annotate("text", x = 54, y = 1000, label = "Future median\noptimum", hjust = 1, size = 5) +
  # annotate('segment', x = 76, y = 300, yend = 300, xend = 65.5, linewidth = 0.75,
  # arrow = arrow(length = unit(0.4, 'cm'))) +
  labs(x = "Days to flowering",
       y = "Number of datapoints\n(genotypes x trials)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.75, 0.4),
        legend.title = element_blank()
  )


deriv_data %>%
  filter(!is.na(flower_50_doy) &
           !is.na(yield_lb_acre) &
           !grepl('irrigated|NO|late|recrop', Location)) %>%
  # group_by(Unif_Name) %>%
  # summarise(mean_flower_dpp = mean(flower_50pct_days_past_planting, na.rm = T),
  #           sd_flower_dpp = sd(flower_50pct_days_past_planting, na.rm = T)) %>%
  summarise(median = median(flower_50pct_days_past_planting, na.rm = T),
          lower = quantile(flower_50pct_days_past_planting, 0.05, na.rm = T),
          upper = quantile(flower_50pct_days_past_planting, 0.95, na.rm = T))
