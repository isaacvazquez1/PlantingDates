# Climate data compilation
The script *0_01_extract_daymet_climate_data.R* will install any missing R packages (see **Requirements** below), download the daymet data into an intermediate data one site at a time, and finally compile all the data into a single file: *daymet_timeseries_cleaned.csv* which will be output into `path/to/output/dir`

## Usage

`Rscript 0_01_extract_daymet_climate_data.R path/to/output/dir`


## Requirements

-   R \~ 4.3.2

### R Packages:

-   `dplyr` \>=1.1.4
-   `lubridate` \>=1.9.3
-   `stringr` \>= 1.5.1
-   `daymetr` 1.7.1
