################################################################################
# Extract and compile daymetr data.

# Miles Alan Moore
# July 18, 2024
################################################################################

rm(list=ls())

packages <- c("dplyr", "lubridate", "stringr", "daymetr")

# Install any missing packages
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load all packages
lapply(packages, library, character.only = TRUE)

# Pathing ---------------------------
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is passed
if (length(args) != 1) {
  stop("Please provide exactly one argument: the DATA_DIR variable in the makefile.")
}

data_dir <- args[1]

daymet_intermediate_dir <- file.path(data_dir, "daymetr")
unique_locs_out_filename <- file.path(daymet_intermediate_dir, "unique_locs_for_daymet.csv")
out_filename <- file.path(data_dir, "daymet_timeseries_cleaned.csv")

if (!dir.exists(daymet_intermediate_dir)) dir.create(daymet_intermediate_dir)

# User Options
extract_daymet <- TRUE # FALSE for debugging purposes really, TRUE 99.99% of the time.

################################################################################
# Extract daymet
################################################################################

# Batch retrieve daymet ...
if(extract_daymet){

  # Remove existing daymet intermediate files
  unlink(list.files(file.path(data_dir, "daymetr"), full.names = TRUE))

  # Read in geocords, filter to distinct, drop NA row & reformat for daymetr
  geocoords <- read.csv(file.path(data_dir, "sunflower_data_simple_v1.csv")) |> 
    dplyr::select(Location, lat, lon) |> 
    dplyr::rename(site = Location, latitude = lat, longitude = lon) |> 
    dplyr::distinct() |> 
    dplyr::filter(!is.na(latitude)) |> 
    dplyr::mutate(site_yr_id = paste0(site,"_",as.character(row_number()))) |> 
    dplyr::select(-site) |> 
    dplyr::rename(site = site_yr_id) |> 
    dplyr::relocate(site, .before = latitude)

  geocoords |> write.csv(unique_locs_out_filename, row.names = FALSE)

  dm_batch <- daymetr::download_daymet_batch(file_location = unique_locs_out_filename,
                                             internal = FALSE, path = daymet_intermediate_dir,
                                             start = 1980, end = 2023)
} else {

  geocoords <- read.csv(unique_locs_out_filename)

}   


# Munge results
flist <- list.files(daymet_intermediate_dir) 

read_in_daymet <- function(f){
  
  # Get metadata (sitename, lat, lon) from fname and header row
  cat('Parsing', f, '...')
  sitename <- gsub("_1980_2023.csv", "", f)
  hdr <- readLines(file.path(daymet_intermediate_dir, f),n = 1)
  lat <- stringr::str_extract(hdr, "Latitude: (\\d+\\.\\d+) ", group = 1) |> 
    as.numeric()
  lon <- stringr::str_extract(hdr, " Longitude: (\\-\\d+\\.\\d+)", group = 1) |> 
    as.numeric()
  
  # Read in and make metadata columns
  df <- read.csv(file.path(daymet_intermediate_dir, f), skip = 6)
  df$site <- sitename
  df$latitude <- lat
  df$longitude <- lon 
  colnames(df) <- gsub("\\.\\.*?", "", colnames(df))

  cat('Done.\n')
  return(df)
}

clean_daymet <- function(df){
  cat('Reformatting daymet data...\n')
  df <- df |> 
  dplyr::full_join(geocoords, by = c('site', 'latitude', 'longitude')) |> 
    dplyr::rename(location = site) |> 
    dplyr::mutate(
      date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"),
      location = gsub("_\\d+$", "", location)
      ) |>
    dplyr::relocate(location, .before = year) |> 
    dplyr::relocate(latitude, .after = location) |> 
    dplyr::relocate(longitude, .after = latitude) |> 
    dplyr::relocate(date, .after = longitude)

  cat('Finished.\n')
  return(df)
}

daymet_dat <- lapply(flist, function(f){
  if (f != 'unique_locs_for_daymet.csv'){  
    df <- read_in_daymet(f)
    return(df)
  } else {
    return(data.frame())
  }
  } ) |> dplyr::bind_rows() |> 
  clean_daymet()

# Write it out! ## -------------------------------------------------------------
cat('Finished extracting DayMet! File written out to', out_filename, '\n')
daymet_dat |> 
  write.csv(out_filename, row.names = FALSE)

# End of Script #
