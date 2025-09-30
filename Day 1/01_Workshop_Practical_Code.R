# ==============================================================================
# PRACTICAL SESSION: Linking Climate Data with Health Data for Time-series Analysis
# ==============================================================================
# Author: Shivang (Karolinska Institutet)
# Workshop: Population Health Impacts from Climate Extremes and Climatic Factors
# Purpose: Learn how to process spatial climate data and link it with health data
# ==============================================================================

# LEARNING OBJECTIVES:
# 1. Understand spatial intersection of population and administrative boundaries
# 2. Calculate population-weighted temperature exposures by district
# 3. Create time series datasets for distributed lag non-linear models (DLNM)
# 4. Link climate exposures with health outcomes data

# ==============================================================================
# PART A: SPATIAL DATA PREPARATION AND INTERSECTION
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PACKAGE LOADING AND SETUP
# ------------------------------------------------------------------------------

# Load required packages for spatial analysis and data manipulation
# install.packages(c("terra","sf","data.table","exactextractr")) # Run if needed
library(terra)        # For raster and vector spatial data handling
library(sf)           # For simple features (vector) spatial operations
library(data.table)   # For fast data manipulation
library(exactextractr) # For precise raster extraction from polygons
library(lubridate)    # For date/time manipulation
library(gnm)          # For generalized non-linear models
library(dlnm)         # For distributed lag non-linear models
library(ggplot2)      # For data visualization

# ------------------------------------------------------------------------------
# 1. PROJECT SETUP AND FILE PATHS
# ------------------------------------------------------------------------------

# EXERCISE 1: Update the base directory path to your working directory
# Set your base folder path (MODIFY THIS PATH!)
base_dir <- "REPLACE/ME/WITH/YOUR/PATH"  # Uncomment and modify this line

# Define file paths for input data
pop_path <- file.path(base_dir, "totalbefolkning_stockholm_lan.gpkg")      # Population grid data
admin_path <- file.path(base_dir, "stockholm_lan_district.gpkg")           # Administrative districts
out_ix_gpkg <- file.path(base_dir, "district_population.gpkg")             # Output intersection file
tas_nc_path <- file.path(base_dir, "tas_NORDIC-3_SMHI-UERRA-Harmonie_RegRean_v1_Gridpp_v1.0.1_day_20180701-20180731.nc") # Temperature data

# VARIABLES EXPLAINED:
# - pop_path: Contains population counts in grid cells (Ruta = grid cell, POP = population)
# - admin_path: Contains administrative district boundaries (distriktskod = district code, distriktsnamn = district name)
# - tas_nc_path: NetCDF file with daily temperature data (tas = temperature at surface)

# TO EXPLORE CLIMATE DATA: GO TO THIS LINK - https://opendata-download-metanalys.smhi.se/gridclim/
# ------------------------------------------------------------------------------
# 2. LOAD SPATIAL DATA
# ------------------------------------------------------------------------------

# Load population grid data as vector (polygon) format
popu <- terra::vect(pop_path)
cat("Population data loaded. Number of grid cells:", nrow(popu), "\n")

# Load administrative district boundaries
admin <- terra::vect(admin_path)
cat("Administrative data loaded. Number of districts:", nrow(admin), "\n")

# EXERCISE 2: Examine the data structure
# Check the column names and first few rows of each dataset
print("Population data columns:")
print(names(popu))
print("Administrative data columns:")
print(names(admin))

# ------------------------------------------------------------------------------
# 3. SPATIAL INTERSECTION AND AREA CALCULATIONS
# ------------------------------------------------------------------------------

# CONCEPT: Spatial intersection creates new polygons where population grid cells
# overlap with administrative districts. This allows us to calculate how much of
# each population grid cell falls within each district.

# Perform spatial intersection between population grid and administrative boundaries
cat("Performing spatial intersection...\n")
ix <- terra::intersect(popu, admin)
cat("Intersection complete. Number of intersected polygons:", nrow(ix), "\n")

# Calculate area of each intersected polygon in square meters
ix$area <- terra::expanse(ix)

# Calculate area fraction in square kilometers
# This represents the area of each intersection polygon
ix$frac <- ix$area / 1e+06  # Convert from m² to km²

# EXERCISE 3: Understanding the intersection
# The 'frac' variable represents the area of intersection pieces
# Why is this important for population weighting?

# Save the intersection results for use in Part B
writeVector(ix, out_ix_gpkg, filetype = "GPKG", overwrite = TRUE)
cat("Intersection data saved to:", out_ix_gpkg, "\n")

# ------------------------------------------------------------------------------
# 4. VISUALIZATION AND QUALITY CHECK
# ------------------------------------------------------------------------------

# Create side-by-side plots to visualize the spatial intersection
par(mfrow = c(1,2))

# Plot 1: Original administrative districts with population grid overlay
plot(admin, main = "Districts with Population Grid", 
     col = "lightblue", border = "blue", lwd = 2)
plot(popu, add = TRUE, border = "grey80", lwd = 0.5)

# Plot 2: Intersected polygons (result of spatial intersection)
plot(ix, main = "Intersected Polygons (Grid × District)", 
     col = topo.colors(10), border = "grey80", lwd = 0.3, values = ix$POP)
max(ix$POP)
# Reset plotting parameters
par(mfrow = c(1,1))

# EXERCISE 4: Interpretation
# Look at the plots above. What do you observe about the intersection results?
# How many intersection polygons were created compared to the original datasets?

# ==============================================================================
# PART B: POPULATION-WEIGHTED TEMPERATURE CALCULATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 5. HELPER FUNCTIONS FOR POPULATION WEIGHTING
# ------------------------------------------------------------------------------

# FUNCTION 1: Calculate population-weighted averages
# This function computes population-weighted temperature for each district
weighted_population <- function(exposure, shap, days) {
  # Convert spatial data and exposure values to data.table format for fast processing
  shap <- as.data.table(as.data.frame(shap))
  exposure <- as.data.table(exposure)
  
  # Attach extracted temperature values to shapefile data
  shap[, var := exposure[[1]]]
  
  # Remove polygons with missing temperature data (e.g., over water bodies)
  cat('Removing NA values...\n')
  shap <- shap[complete.cases(shap$var)]
  exposure <- exposure[complete.cases(exposure[[1]])]
  
  # POPULATION WEIGHTING CALCULATION:
  # Step 1: Calculate population contribution of each intersection polygon
  shap[, pop_inc := POP * frac]  # Population × area fraction
  
  # Step 2: Calculate weights (proportion of total district population)
  shap[, weights := pop_inc / sum(pop_inc), by = distriktskod]
  
  # EXERCISE 5: Understanding weights
  # The 'weights' variable represents what proportion of the district's total
  # population lives in each intersection polygon. Why do we need this?
  
  # Assign formatted date-time column names to exposure data
  names(exposure) <- format(days, "%Y-%m-%dT%H:%M:%SZ")
  
  # Apply population weights to temperature values
  cat('Calculating population-weighted temperature...\n')
  exposure <- exposure * shap$weights
  
  # Add district identifier for grouping
  exposure[, location := shap$distriktskod]
  
  # Aggregate weighted values by district (sum of weighted temperatures)
  exposure <- exposure[, lapply(.SD, sum), by = location]
  
  # Reshape from wide to long format: one row per district-date combination
  cat('Reshaping data to long format...\n')
  exposure <- melt(exposure, id.vars = "location", 
                   variable.name = "date", value.name = 'tas')
  
  # Sort by location for consistent output
  exposure <- exposure[order(location)]
  
  cat('Population weighting complete!\n')
  return(exposure)
}

# FUNCTION 2: Process NetCDF climate data and calculate district averages
create_dte <- function(file, v) {
  # Load temperature raster data from NetCDF file
  nc <- terra::rast(file)
  cat("Loaded NetCDF with", nlyr(nc), "time layers\n")
  
  # Extract date information from raster time dimension
  sday <- format(min(terra::time(nc)), "%Y-%m-%d")
  eday <- format(max(terra::time(nc)), "%Y-%m-%d")
  year <- as.numeric(format(min(terra::time(nc)), "%Y"))
  
  cat("Processing data from", sday, "to", eday, "\n")
  
  # Generate sequence of daily timestamps
  days <- seq(from = as.POSIXct(sday, tz = 'UTC'), 
              to = as.POSIXct(eday, tz = 'UTC') + 1*60*60*23, 
              by = '1 day')
  
  # Assign date names to raster layers
  names(nc) <- days
  
  # Reproject raster to match vector coordinate system
  # Coordinate Reference System (CRS) is a standard that 
  # tells maps how to convert real-world locations into 
  # coordinates (and back) so different layers line up correctly.
  nc <- terra::project(nc, terra::crs(v))
  
  # Extract mean temperature for each intersection polygon
  cat("Extracting temperature values for each polygon...\n")
  e <- exact_extract(nc, st_as_sf(v), fun = 'mean', max_cells_in_memory = 1e+09)
  
  # Calculate population-weighted district averages
  dte <- weighted_population(e, v, days)
  
  return(dte)
}

# ------------------------------------------------------------------------------
# 6. CALCULATE POPULATION-WEIGHTED TEMPERATURE BY DISTRICT
# ------------------------------------------------------------------------------

# Load the intersection data created in Part A
district_pop <- st_read(out_ix_gpkg)
cat("Loaded intersection data with", nrow(district_pop), "polygons\n")

# Process temperature data and calculate population-weighted district averages
cat("Processing temperature data...\n")
dte <- create_dte(tas_nc_path, district_pop)

# EXERCISE 6: Temperature unit conversion
# Convert temperature from Kelvin to Celsius
# QUESTION: What value should we subtract from Kelvin to get Celsius?
value #  <- VALUE  Kelvin to Celsius conversion factor
dte[, tas := tas - value]  # tas [K] → tas [°C]

# Display summary of the processed data
cat("Temperature data summary:\n")
print(summary(dte$tas))
cat("Date range:", range(dte$date), "\n")
cat("Number of districts:", length(unique(dte$location)), "\n")

# Save the processed temperature data
csv_out <- file.path(base_dir, "district_tas.csv")
fwrite(dte, csv_out)
cat("Temperature data saved to:", csv_out, "\n")

# ==============================================================================
# PART C: LINKING CLIMATE DATA WITH HEALTH DATA
# ==============================================================================

# ------------------------------------------------------------------------------
# 7. CREATE COMPREHENSIVE TIME SERIES FOR DLNM ANALYSIS
# ------------------------------------------------------------------------------

# CONCEPT: DLNM (Distributed Lag Non-linear Models) require a complete time series
# with daily observations for each individual, including days without health events.

# Load Asthmaexacerbations data (cases when health outcomes occurred)
# including primary care, hospitalizations, and emergency visits
cases <- fread(paste0(base_dir, "health_events.csv"))
cat("Loaded health events data with", nrow(cases), "observations\n")

# Sort cases by individual ID and date for consistent processing
setorder(cases, id, date)

# EXERCISE 7: Understanding the health data structure
# Examine the cases dataset
print("Health events data structure:")
print(str(cases))
print("First few rows:")
print(head(cases))

# Create complete time series for DLNM analysis
# Date range: January 01, 2018 to December 31, 2018 (includes lag period from Dec 21, 2017 (10 days prior))
date_range <- seq(as.Date("2017-12-21"), as.Date("2018-12-31"), by = "day")
all_ids <- unique(cases$id)

cat("Creating time series for", length(all_ids), "individuals over", 
    length(date_range), "days\n")

# Generate all possible combinations of individual ID and date
time_series <- data.table(expand.grid(id = all_ids, date = date_range))
setorder(time_series, id, date)

# ------------------------------------------------------------------------------
# 8. ADD HEALTH OUTCOME INDICATORS
# ------------------------------------------------------------------------------

# Mark health events with outcome = 1
cases[, outcome := 1]

# Merge health events with complete time series
# All unmatched dates get outcome = 0 (no health event)
time_series <- merge(time_series, cases, by = c("id", "date"), all.x = TRUE)
time_series[is.na(outcome), outcome := 0]

# EXERCISE 8: Understanding the outcome variable
# What does outcome = 1 vs outcome = 0 represent?
# How many days have outcome = 1 vs outcome = 0?
cat("Health outcome summary:\n")
print(table(time_series$outcome))

# ------------------------------------------------------------------------------
# 9. ADD TEMPERATURE EXPOSURE DATA
# ------------------------------------------------------------------------------

# Load daily temperature data for the study location
temp_data <- fread(paste0(base_dir, "tas_212098_daily.csv"))
temp_data$date <- as.Date(temp_data$date)
names(temp_data)[2] <- "tmean"  # Standardize temperature column name

# Convert temperature from Kelvin to Celsius
temp_data$tmean <- temp_data$tmean - 273.15

# EXERCISE 9: Temperature data exploration
cat("Temperature data summary:\n")
print(summary(temp_data$tmean))

# Merge temperature data with time series
time_series <- merge(time_series, temp_data, by = "date", all.x = TRUE)
setorder(time_series, id, date)

# ------------------------------------------------------------------------------
# 10. ADD INDIVIDUAL CHARACTERISTICS (CONFOUNDERS)
# ------------------------------------------------------------------------------

# Load cohort information containing individual characteristics
cohort_info <- fread(paste0(base_dir, "cohort_info.csv"))

# Add smoking status as a potential confounder
time_series <- merge(time_series, cohort_info[, .(id, smoking)], by = "id", all.x = TRUE)

# Final sorting by individual ID and date
setorder(time_series, id, date)

# EXERCISE 10: Final dataset examination
cat("Final time series dataset summary:\n")
print(str(time_series))
cat("Number of observations:", nrow(time_series), "\n")
cat("Number of individuals:", length(unique(time_series$id)), "\n")
cat("Date range:", range(time_series$date), "\n")

# Save the complete time series dataset for DLNM analysis
final_output <- paste0(base_dir, "time_series.csv")
fwrite(time_series, file = final_output)
cat("Complete time series dataset saved to:", final_output, "\n")

# ==============================================================================
# END OF PRACTICAL SESSION
# ==============================================================================
