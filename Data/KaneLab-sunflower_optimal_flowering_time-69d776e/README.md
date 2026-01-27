# Scripts for "Climate drives variation in optimal phenology: 46 years of multi-environment trials in sunflower" 

Eliza Clark

Created: 17 April 2025  
Last updated: 28 April 2025


## How to run

### Download input data files

Put all input files in a directory that you make called `data/` in the root of the repository. The scripts will look for the data files in this directory.

1. *sunflower_data_simple_v1.csv* Download this file from Ag Data Commons. [Citation with doi.]

2. *daymet_timeseries_cleaned.csv* This is the climate data from Daymet. To get this file, run the make file in the scripts/climate directory. 

3. *grin_dayflow.csv* Download this data from  GRIN-Global at https://npgsweb.ars-grin.gov/gringlobal/descriptordetail?id=79038 

### Run Scripts

After downloading the above data, you can run the scripts in order. After running all scripts once, all intermediate files will be saved and each script may be run without running the others first. The data derived from these files and resulting figures will be saved in the corresponding folders. All output is already in the repository, so you can skip this step if you just want to look at the results.

### Descriptions of scripts
The description of each script, files it uses, and files it creates is in the word document "Script_Descriptions" in the scripts directory.

### List of figures and tables in manuscript
Table 1. county/table_dpp.png  
Figure 1. -- [map not made in R]  
Figure 2. county/state_r2_county_unadj.png  
Figure 3. county/climate_mult_full.png  
Figure 4. county/DPP_opt_clim.png  
Figure 5. county/DPP_future.png  
Figure 6. variation_DPP.png  

**Supplement**  
Table S1. county/table_all.png  
Table S2.  exact_location/table_all.png  
Figure S1. climate_trends.png  
Figure S2. county/R2_sample_size.png  
Figure S3. county/single_site_yr_examples.png  
Figure S4. county/year_r2_county_unadj.png  
Figure S5. county/year_county_mult.png  
Figure S6. county/yield_optimum_relationship.png  
Figure S7. county/dpp_distribution.png  
Figure S8. county/ALL_resids.png  
Figure S9. county/DOY_opt_clim.png  
Figure S10.  exact location/DPP_opt_clim.png  
Figure S11. county/fixed_DPP_opt_clim.png  
Figure S12. county/fixed_DOY_opt_clim.png  
Figure S13. county/DOY_future.png  
Figure S14. nass_planting_dates.png (not in repository)  


