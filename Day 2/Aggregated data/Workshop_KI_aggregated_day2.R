
##########################################################################################################################
##########################################################################################################################
#######  
#######              Workshop in Population Health Impacts from Climate Extremes and Climatic Factors
#######
#######     Modelling Associations Between Climatic Factors & Health Outcomes Using Aggregated-Level Data - Day 2
#######
#######                              Blanca Paniello-Castillo (blanca.paniello@isglobal.org)
#######                                             ISGlobal, Barcelona, Spain
#######                                 Universitat Pompeu Fabra (UPF), Barcelona, Spain
#######  
##########################################################################################################################
##########################################################################################################################

# Cleaning environment
rm(list = ls())

####################### Loading libraries

# ---- Core Modeling ----
library(dlnm)        # Crossbasis, crosspred, crossreduce
library(splines)     # Natural splines (ns) for temporal trend and seasonality
library(tsModel)     # Lagged time series structures for DLNM models

# ---- Epidemiology & Meta-analysis ----
library(mixmeta)     # Meta-analysis of reduced estimates (second stage)
library(mvmeta)      # Alternative older package, only if explicitly needed

# ---- Data Handling ----
library(dplyr)       # Data wrangling and joins (filter, left_join, etc.)
library(tibble)      # rownames_to_column() and tidy tibble manipulation
library(lubridate)   # Date handling functions like wday(), year(), month()
library(ISOweek)     # ISOweek2date() for epidemiological week handling

# ---- Visualisation ----
library(ggplot2)     # For advanced visualizations (optional but useful)
library(scales)      # Better axis formatting and scales for ggplot2
# Base R plotting functions (plot, lines, polygon) also work without extra libraries

# ---- Statistical Utilities ----
library(MASS)        # For GLM-related functions (quasipoisson, etc.)


###################################################################################
########### Data Preparation for DLNM Analysis
###################################################################################

# --- Set working directory (adjust path to your system) ---
setwd("Z:/proj/bpaniello/Conferences_Workshops/KI_Oct2025") # Change to project folder


##############################################################
########### 1. Load and Prepare Temperature Data
##############################################################

# Load weekly temperature data
load("datatable_temp_es_se_fr.RData")

# Ensure columns are numeric
datatable_temp$year <- as.numeric(datatable_temp$year)  # Year of observation
datatable_temp$woy  <- as.numeric(datatable_temp$woy)   # Week of year
datatable_temp$temp <- as.numeric(datatable_temp$temp)  # Temperature (Kelvin)

# Convert temperature from Kelvin to Celsius
datatable_temp$temp <- datatable_temp$temp - 273.15


##############################################################
########### 2. Load and Prepare Mortality Data
##############################################################

# Load weekly mortality data
load("datatable_mort_es_se_fr.RData")

# Ensure numeric format
datatable_mort$year <- as.numeric(datatable_mort$year)
datatable_mort$woy  <- as.numeric(datatable_mort$woy)   # Week of year
datatable_mort$mort <- as.numeric(datatable_mort$mort)  # Mortality count

# Merge mortality with combined temperature table
datatable_data <- left_join(datatable_mort, datatable_temp,
                            by = c("location", "year", "woy"))

# Clean up memory by removing temporary datasets
rm(datatable_mort, datatable_temp)

# Add actual calendar date (Thursday of each ISO week) 
# ISOweek format: YYYY-Www-d (d = weekday, 4 = Thursday)
datatable_data$date <- ISOweek2date(
  paste0(datatable_data$year, "-W", 
         sprintf("%02d", datatable_data$woy), "-", 4)
)

# Print total mortality count for a quick sanity check
print(paste0("Total Mortality: ", sum(datatable_data$mort, na.rm = TRUE)))


##############################################################
########### 4. Load Metadata Table
##############################################################

# Load METATABLE: contains additional descriptive information
# (e.g., locations, region codes, classification variables)

load("METATABLE_es_se_fr.RData")

##############################################################
########### 5. Load Functions
##############################################################

setwd("Z:/proj/bpaniello/p20231006_BP_mortality_changes_COVID/code/R functions")

# Required Codes
source("FWALD.R") # Wald Test

###################################################################################
############## Data Preparation and Parameter Definition
###################################################################################

# Here we are going to define the parameters that we will use later on in the first and second stage

##############################################################
########### 6. Define Calibration and Prediction Periods
##############################################################

# Here we define two key time windows:
#   1) Calibration Period:
#        - Used to fit the model and estimate the temperature–mortality relationship.
#   2) Prediction Period:
#        - The period where the calibrated model is applied to estimate
#          temperature-related mortality.
#
# These periods can be the same (in-sample prediction) or different 
# (out-of-sample prediction). 
#
# Example:
#   - SAME periods: analyze and validate historical impacts in one timeframe.
#   - DIFFERENT periods: apply past risk estimates to later years to explore
#     how mortality would look under previous conditions (counterfactual analysis).
#
# In this script, both periods are set to the same dates for simplicity.

##############################
# Calibration Period
##############################

# Calibration period: 2000-2019
day1_cali <- 6   ; mon1_cali <- 1   ; yea1_cali <- 2000  # Start: 06-Jan-2000
day2_cali <- 26  ; mon2_cali <- 12  ; yea2_cali <- 2019  # End:   26-Dec-2019

date1_cali <- as.Date(paste(yea1_cali, mon1_cali, day1_cali, sep = "-"))
date2_cali <- as.Date(paste(yea2_cali, mon2_cali, day2_cali, sep = "-"))

# Ensure both dates are Thursdays (ISO week reference day)
if(lubridate::wday(date1_cali, week_start = 1) != 4 | 
   lubridate::wday(date2_cali, week_start = 1) != 4){
  stop("ERROR: Invalid Calibration Dates !!!")
}

##############################
# Prediction Period
##############################

# Prediction period: 2000-2019
day1_pred <- 6   ; mon1_pred <- 1   ; yea1_pred <- 2000
day2_pred <- 26  ; mon2_pred <- 12  ; yea2_pred <- 2019

date1_pred <- as.Date(paste(yea1_pred, mon1_pred, day1_pred, sep = "-"))
date2_pred <- as.Date(paste(yea2_pred, mon2_pred, day2_pred, sep = "-"))

# Ensure both dates are Thursdays (ISO week reference day)
if(lubridate::wday(date1_pred, week_start = 1) != 4 | 
   lubridate::wday(date2_pred, week_start = 1) != 4){
  stop("ERROR: Invalid Prediction Dates !!!")
}

##############################################################
########### 7. Create Vectors and Parameters
##############################################################

# ---- Percentiles for cumulative exposure-response predictions ----
PRED_PRC <- sort(unique(c(
  seq( 0.0,   1.0, 0.1),  # very low range
  seq( 1.5,   5.0, 0.5),  # low range
  seq( 6.0,  94.0, 1.0),  # mid range
  seq(95.0,  98.5, 0.5),  # high range
  seq(99.0, 100.0, 0.1)   # very high range
) / 100))

# Validate percentile vector 
if(any(0 > PRED_PRC | PRED_PRC > 1)){ 
  stop("ERROR: Invalid Percentile Vector for Cumulative Exposure-Response !!!")
}

# ---- Minimum Mortality Temperature (MMT) Calculation ----
# Fixed lower/upper percentiles
MIN_PMMT <-  5  / 100  # 5th percentile
MAX_PMMT <- 100 / 100  # 100th percentile

# ---- Exposure-response spline settings ----
VAR_FUN <- "ns"           # Spline type ("ns" = natural spline)
VAR_DEG <- NA             # Spline degree (NA = default)

# Knot locations (percentiles of temperature distribution)
VAR_PRC <- c(10, 50, 90) / 100

# ---- Lag-response settings ----
MIN_LAG <- 0  # Minimum lag in weeks
if(MIN_LAG < 0) stop("ERROR: Invalid MIN_LAG !!!")

MAX_LAG <- 3  # Maximum lag in weeks
if(MAX_LAG <= MIN_LAG) stop("ERROR: Invalid MAX_LAG !!!")

# Degrees of Freedom per Year for the Seasonal and Long-Term Trends
DF_SEAS = 8
if( DF_SEAS <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ) }

##############################################################
########### 8. Region and Country Vectors
##############################################################

# Region codes (e.g., NUTS regions)
vREG <- METATABLE$location # Getting the code of the regions
nREG <- length(vREG) # Getting the number of regions

# English region names
vREG_NAME <- METATABLE$name_eng

# Country code of each region
vCOU_of_REG <- METATABLE$nuts0

# Unique countries present
vCOU <- unique(METATABLE$nuts0) # Getting the code of the countries
nCOU <- length(vCOU) # Getting the number of countries

# Assign readable country names
# Example: Rename "SE" to "Sweden"
vCOU_NAME <- vCOU
vCOU_NAME[which(vCOU_NAME == "SE")] <- "Sweden"
vCOU_NAME[which(vCOU_NAME == "ES")] <- "Spain"
vCOU_NAME[which(vCOU_NAME == "FR")] <- "France"

##############################################################
########### 9. Restrict Data to Selected Periods
##############################################################

# It's important to include additional weeks (MAX_LAG - 3 weeks) at the start for lagged calculations
datatable_cali <- datatable_data[
  which(date1_cali <= datatable_data$date &
          datatable_data$date <= date2_cali + max(7 * MAX_LAG, 0)), ]

datatable_pred <- datatable_data[
  which(date1_pred <= datatable_data$date &
          datatable_data$date <= date2_pred + max(7 * MAX_LAG, 0)), ]

##############################################################
########### 10. Create Data Lists by Region
##############################################################

# Split data into lists, one element per region
datalist_cali <- lapply(vREG, function(x) datatable_cali[datatable_cali$location == x, ])
names(datalist_cali) <- vREG

datalist_pred <- lapply(vREG, function(x) datatable_pred[datatable_pred$location == x, ])
names(datalist_pred) <- vREG

##############################################################
########### 11. Calculate Total Demographics
##############################################################

print(paste0("Number of Regions: ", nREG))
print(paste0("Number of Countries: ", nCOU))

##############################################################
########### 12. Add Week-of-Period Variable
##############################################################

# Sequential week number within each region
# Used as a time spline covariate in DLNM models
# We will use it adjust for seasonality in the first stage
for(iREG in 1:nREG){
  datalist_cali[[iREG]]$wop <- 1:length(datalist_cali[[iREG]]$woy)
}
rm(iREG)

##############################################################
########### 13. Model Formula and Random Effects
##############################################################

# Meta-analysis formula for Best Linear Unbiased Predictions (BLUP)
# TEMP_AVG = mean temperature
# TEMP_IQR = interquartile range temperature
# COEF_MODEL = Coefficients for each of the regions

# Aside from TEMP_AVG and TEMP_IQR, we could also add other meta-predictors
FORMULA_META <- COEF_MODEL ~ TEMP_AVG + TEMP_IQR 


###################################################################################
########### 1st Stage - Location Specific DLNM
###################################################################################

print("")
print("= Calculation of the Location-Specific Associations =")
print("")

##############################################################
########### 14. Initialize Storage Structures
##############################################################

# Reduced Coefficients matrix
# Stores regression coefficients for each region
COEF_MODEL <- matrix(
  data = NA, 
  nrow = nREG, 
  ncol = length(VAR_PRC) + 1, 
  dimnames = list(vREG)
)

# Reduced covariance matrices (one per region)
VCOV_MODEL <- vector("list", nREG)
names(VCOV_MODEL) <- vREG

# Cumulative exposure-response predictions before meta-analysis
CROSS_PRED_REG_NOMETA <- vector("list", nREG)
names(CROSS_PRED_REG_NOMETA) <- vREG


##############################################################
########### 15. Loop over Regions
##############################################################

for(iREG in 1:nREG){
  
  print(paste0("Region ", iREG, ": ", vREG_NAME[iREG], " (", vREG[iREG], ")"))
  
  # ------------------------------
  # Define Formulas
  # ------------------------------
  # FORMULA_SEA: seasonality-only model
  # FORMULA_CRB: seasonality + temperature cross-basis
  # wop = week-of-period variable for spline
  
  FORMULA_SEA <- mort ~ ns(wop, df = round(DF_SEAS*length(wop)*7/365.25))
  FORMULA_CRB <- mort ~ ns(wop, df = round(DF_SEAS*length(wop)*7/365.25)) + CROSS_BASIS
  
  # ------------------------------
  # Fit Seasonality Model
  # ------------------------------
  GLM_MODEL_SEA <- glm(
    formula = FORMULA_SEA,
    data = datalist_cali[[iREG]],
    family = quasipoisson, 
    na.action = "na.exclude")
  
  # Predict seasonality
  datalist_cali[[iREG]]$mort_pred_seas8 <- predict(GLM_MODEL_SEA, type = "response")
  rm(GLM_MODEL_SEA)
  
  # ------------------------------
  # Create Cross-Basis for Temperature
  # ------------------------------
  CROSS_BASIS <- crossbasis(
    datalist_cali[[iREG]]$temp,
    lag = c(MIN_LAG, MAX_LAG),
    argvar = list(
      fun = VAR_FUN,
      knots = quantile(datalist_cali[[iREG]]$temp, VAR_PRC, na.rm = TRUE),
      Boundary.knots = range(datalist_cali[[iREG]]$temp, na.rm = TRUE)),
    arglag = list(fun = "integer")
  )
  
  # ------------------------------
  # Fit Cross-Basis Model
  # ------------------------------
  GLM_MODEL_CRB <- glm(
    formula = FORMULA_CRB, 
    data = datalist_cali[[iREG]], 
    family = quasipoisson, 
    na.action = "na.exclude"
  )
  
  # ------------------------------
  # Cumulative Exposure-Response (Before Centring)
  # ------------------------------
  CROSS_PRED_REG_NOMETA[[iREG]] <- crosspred(
    CROSS_BASIS, 
    GLM_MODEL_CRB, 
    at = quantile(datalist_cali[[iREG]]$temp, PRED_PRC, na.rm = TRUE))
  
  # ------------------------------
  # Determine Minimum Mortality Temperature (MMT)
  # ------------------------------
  # The goal here is to find the temperature(s) associated with the lowest predicted mortality.
  # CROSS_PRED_REG_NOMETA[[iREG]]$allRRfit contains the predicted relative risks (RRs) 
  # for a range of temperatures defined in predvar. 
  
  #   - Select the global minimum within a predefined percentile range of temperatures (MIN_PMMT to MAX_PMMT).
  # This ensures we choose a biologically and statistically plausible MMT even if no clear local minima exist.
  
  # Use global minimum within predefined percentile range
  MMT <- CROSS_PRED_REG_NOMETA[[iREG]]$predvar[
    which(PRED_PRC == MIN_PMMT) - 1 + which.min(
      CROSS_PRED_REG_NOMETA[[iREG]]$allRRfit[
        which(PRED_PRC == MIN_PMMT):which(PRED_PRC == MAX_PMMT)
      ]
    )
  ]
  
  # ------------------------------
  # Cumulative Exposure-Response (With Centring at MMT)
  # ------------------------------
  CROSS_PRED_REG_NOMETA[[iREG]] <- crosspred(
    CROSS_BASIS, 
    GLM_MODEL_CRB,
    at = quantile(datalist_cali[[iREG]]$temp, PRED_PRC, na.rm = TRUE),
    bylag = 1, 
    cen = MMT
  )
  
  # ------------------------------
  # Extract Reduced Coefficients and Covariance Matrices
  # ------------------------------
  REDUCED <- crossreduce(CROSS_BASIS, GLM_MODEL_CRB, cen = MMT)
  COEF_MODEL[iREG,] <- coef(REDUCED)
  VCOV_MODEL[[iREG]] <- vcov(REDUCED)
  rm(REDUCED)
  
  # Clean up temporary objects for current iteration
  rm(FORMULA_SEA, FORMULA_CRB, CROSS_BASIS, GLM_MODEL_CRB, MMT)
  
} # End of region loop

rm(iREG)


##################################################################################
########### 2nd Stage - Meta-analysis and BLUPs
##################################################################################

print("")
print("= Calculation of the Best Linear Unbiased Predictions =")

##############################################################
########### 16. Define Temperature Meta-Predictors
##############################################################

# TEMP_AVG: Mean temperature in each region
# TEMP_IQR: Interquartile range (IQR) of temperature in each region
TEMP_AVG = sapply(datalist_cali, function(x) mean(x$temp, na.rm = TRUE))
TEMP_IQR = sapply(datalist_cali, function(x) IQR(x$temp, na.rm = TRUE))

# Placeholder for Minimum Mortality Temperature (MMT) per region
MMT_REG = array(NA, dim = c(nREG), dimnames = list(vREG))

# Placeholder for cumulative exposure-response predictions (BLUPs) per region
CROSS_PRED_REG_META = vector("list", nREG)
names(CROSS_PRED_REG_META) = vREG


##############################################################
########### 17. Multivariate Meta-Analysis of Reduced Coefficients
##############################################################

# Combines first-stage coefficients across regions using random or fixed-effects models
# mixmeta() uses the variance-covariance matrices (VCOV_MODEL) for weighting, also 
# including country-level random effects in meta-analysis
# The algorithm iterates until convergence, which is reached when the log-likelihood changes very little between iterations.

fCOU = factor(vCOU_of_REG)
fREG = factor(vREG)

MULTIVAR = mixmeta(FORMULA_META,
                   VCOV_MODEL,
                   data = data.frame(vREG = vREG),
                   control = list(showiter = TRUE, igls.inititer = 10),
                   method = "reml",
                   random =~ 1 | fCOU/fREG)

rm(fCOU,fREG)

print(summary(MULTIVAR)) # Print meta-analysis summary


##############################################################
########### 18. Wald Test for Meta-Predictors
##############################################################

# Evaluates the significance of meta-predictors (e.g., TEMP_AVG, TEMP_IQR)
# FWALD() calculates the Wald statistic and p-value

if(length(summary(MULTIVAR)$lab$p) > 1){
  for(iMETAPRED in 1:length(summary(MULTIVAR)$lab$p)){
    print(paste0("Wald Test of ", summary(MULTIVAR)$lab$p[iMETAPRED], 
                 ": p = ", 
                 sprintf("%.10f", FWALD(MULTIVAR, 
                                        summary(MULTIVAR)$lab$p[iMETAPRED]))))
  }
  rm(iMETAPRED)
}

##############################################################
########### 19. Calculate Best Linear Unbiased Predictions (BLUPs)
##############################################################

# BLUPs combine first-stage estimates with overall meta-estimates
# Shrinks extreme region-specific estimates toward the overall mean
BLUP = blup(MULTIVAR, vcov = TRUE)  # Returns both predicted coefficients and variance-covariance

##############################################################
########### 20. Onebasis Function
##############################################################

# Creates a spline basis for temperature to model non-linear relationships
# Input parameters:
#   - fun: type of spline ("ns" = natural spline, "bs" = B-spline)
#   - knots: internal knots
#   - Boundary.knots: min/max values (for ns)
#   - degree: polynomial degree (for bs)
# Output: basis matrix for prediction

for(iREG in 1:nREG){
  
  # Create basis for temperature
  BASIS_VAR = onebasis(quantile(datalist_cali[[iREG]]$temp, PRED_PRC, na.rm = TRUE),
                       fun = VAR_FUN,
                       knots = quantile(datalist_cali[[iREG]]$temp, VAR_PRC, na.rm = TRUE),
                       Boundary.knots = range(datalist_cali[[iREG]]$temp, na.rm = TRUE))
  
  # Percentiles of temperature for prediction
  PRC_VAR = quantile(datalist_cali[[iREG]]$temp, PRED_PRC, na.rm = TRUE)
  
  ##############################################################
  ########### 21. Cumulative Exposure-Response (Before Centering)
  ##############################################################
  
  # crosspred() predicts relative risks (RR) using basis matrix and BLUP coefficients
  PRED_MORT <- suppressMessages(crosspred(BASIS_VAR, 
                                          coef = BLUP[[iREG]]$blup,
                                          vcov = BLUP[[iREG]]$vcov,
                                          model.link = "log",
                                          at = PRC_VAR))
  
  ##############################################################
  ########### 22. Identify Local Minima (Potential MMT)
  ##############################################################
  # MMT (Minimum Mortality Temperature) is chosen as the local minimum with lowest RR
  
  MMT_REG[iREG] = PRC_VAR[ which(PRED_PRC == MIN_PMMT) - 1 +
                             which.min(PRED_MORT$allRRfit[which(
                               PRED_PRC == MIN_PMMT) : which(
                                 PRED_PRC == MAX_PMMT)])]
  

  ##############################################################
  ########### 23. Cumulative Exposure-Response (Centered at MMT)
  ##############################################################
  
  
  CROSS_PRED_REG_META[[iREG]] = crosspred(BASIS_VAR,
                                          coef = BLUP[[iREG]]$blup,
                                          vcov = BLUP[[iREG]]$vcov,
                                          model.link = "log",
                                          at = PRC_VAR,
                                          cen = MMT_REG[iREG])
  rm(BASIS_VAR, PRC_VAR)
}
rm(iREG)

##############################################################
########### 24. Calculation of Pooled Coefficients Across Regions
##############################################################

# Combine region-specific predictions into a pooled exposure-response for plotting or overall estimates

print("")
print("= Calculation of the Pooled Coefficients =")

# Placeholder for pooled cumulative exposure-response
CROSS_PRED_TOT_META = vector("list", 1)
names(CROSS_PRED_TOT_META) = "SP, SE, and FR"


##############################################################
########### 25. Prepare New Data for Meta-Predictor Based Predictions
##############################################################

# Create a data frame of the average temperature and IQR across all regions
# This will be used to predict the pooled exposure-response from the meta-analysis model
NEW_DATA = data.frame(
  TEMP_AVG  = mean(tapply(TEMP_AVG, vREG, mean)),   # Avg. of regional mean temperatures
  TEMP_IQR  = mean(tapply(TEMP_IQR, vREG, mean))    # Avg. of regional temperature IQRs
)


##############################################################
########### 26. Predict Pooled Coefficients Using Meta-Analysis Results
##############################################################

# predict() generates pooled coefficients for the new data
# vcov = TRUE returns the variance-covariance matrix for uncertainty estimates
MULTIVAR_PRED = predict(MULTIVAR, NEW_DATA, vcov = TRUE, format = "list")
rm(MULTIVAR, NEW_DATA)


##############################################################
########### 27. Calculate Average Multi-Location Temperature Percentiles
##############################################################

# rowMeans() combines the temperature percentiles across all regions
# jitter() adds a small random noise to avoid exact ties
set.seed(13041975)
POOLED_TEMP_AVG = rowMeans(sapply(datalist_cali, function(x) quantile(jitter(x$temp), PRED_PRC, na.rm = TRUE)))


##############################################################
########### 28. Create Spline Basis for Pooled Temperature
##############################################################

# Generate a spline basis matrix for pooled temperature using same spline type as for individual regions
BASIS_VAR = onebasis(POOLED_TEMP_AVG,
                     fun = VAR_FUN,
                     knots = POOLED_TEMP_AVG[paste0(100 * VAR_PRC, ".0%")],
                     Boundary.knots = range(POOLED_TEMP_AVG, na.rm = TRUE))


##############################################################
########### 29. Predict Pooled Cumulative Exposure-Response
##############################################################

# Multiply the basis matrix by predicted coefficients to get relative risks (RR)
PRED_MORT = BASIS_VAR %*% MULTIVAR_PRED$fit


##############################################################
########### 30.  Identify Local Minima for Pooled MMT
##############################################################

# Use global minimum within predefined percentile range
POOLED_CENT_IND = which(PRED_PRC == MIN_PMMT) - 1 + which.min(
  PRED_MORT[which(PRED_PRC == MIN_PMMT) : which(PRED_PRC == MAX_PMMT)])

rm(PRED_MORT)


##############################################################
########### 31. Determine Percentile and Value of Pooled MMT
##############################################################

# Convert the percentile index into a value between 0 and 100
POOLED_CENT_PRC = pmin(pmax(100 * PRED_PRC[POOLED_CENT_IND], 0), 100)
rm(POOLED_CENT_IND)

# Extract the corresponding temperature value for centering the exposure-response
TEMP_AVG_CEN = POOLED_TEMP_AVG[paste0(POOLED_CENT_PRC, ".0%")]
rm(POOLED_CENT_PRC)


##############################################################
########### 32. Calculate Pooled Cumulative Exposure-Response Centered at MMT
##############################################################

# crosspred() computes the cumulative exposure-response using pooled coefficients
# The 'cen' argument centers the RR curve at the pooled MMT
CROSS_PRED_TOT_META[[1]] = crosspred(BASIS_VAR, 
                                     coef = MULTIVAR_PRED$fit, 
                                     vcov = MULTIVAR_PRED$vcov,
                                     model.link = "log", 
                                     at = POOLED_TEMP_AVG, 
                                     cen = TEMP_AVG_CEN)


##############################################################
########### 33. Identify Local Minima in Pooled Curve and Select
###########     Minimum Mortality Temperature (MMT) for Pooled Prediction
##############################################################

# Use global minimum within defined percentile range if no local minima
MMT_TOT = CROSS_PRED_TOT_META[[1]]$predvar[which(
  PRED_PRC == MIN_PMMT) - 1 + which.min(
    CROSS_PRED_TOT_META[[1]]$allRRfit[which(
      PRED_PRC == MIN_PMMT) : which(
        PRED_PRC == MAX_PMMT)])]


##############################################################
########### 34. Recalculate Pooled Exposure-Response Centered at Final MMT
##############################################################

CROSS_PRED_TOT_META[[1]] = crosspred(BASIS_VAR, 
                                     coef = MULTIVAR_PRED$fit, 
                                     vcov = MULTIVAR_PRED$vcov,
                                     model.link = "log", 
                                     at = POOLED_TEMP_AVG, 
                                     cen = MMT_TOT)

# Remove intermediate variables to free memory
rm(MULTIVAR_PRED, POOLED_TEMP_AVG, BASIS_VAR, TEMP_AVG_CEN, COEF_MODEL, VCOV_MODEL)



##################################################################################
########### Plots and figures
##################################################################################

# Set working directory to where outputs and plots will be saved
setwd("Z:/proj/bpaniello/Conferences_Workshops/KI_Oct2025")

##################################################################################
########### Figure: Mortality and Seasonality Time Series 
##################################################################################

# This figure shows observed mortality and the estimated seasonal mortality component for each region over time.

print("")
print("= Figure: Mortality and Seasonality Time Series =")
print("")

# Create output folder for plots if it doesn't exist
foldout = paste0("./dataout/Seasonality_Plot/")
if (!file_test("-d", foldout)) { dir.create(foldout, recursive = TRUE) }

# Create a PDF for saving plots (10 x 14 inches)
pdf(paste0(foldout, "TIME_SERIES_MORT_SEA_REG.pdf"), width = 10, height = 14)

# Layout: 7 plots per page, arranged in 1 column
layout(matrix(seq(1 * 7), nrow = 7, byrow = TRUE))

# Configure graphical parameters for better readability
par(mar = c(3.5, 4, 3, 1),   # margins: bottom, left, top, right
    oma = c(1, 1, 1, 1),      # outer margins
    mgp = c(2, 0.7, 0),       # axis title, labels, line distance
    las = 1,                   # axis labels horizontal
    cex.axis = 0.8,            # smaller axis labels
    cex.main = 0.9,            # smaller plot titles
    cex.lab = 0.9)             # smaller axis titles

# Loop over all regions
for (iREG in 1:nREG) {
  # Plot observed mortality
  plot(datalist_cali[[iREG]]$date, datalist_cali[[iREG]]$mort, 
       col = rgb(0.00, 0.00, 0.00), lwd = 1.5, lty = 1, type = "l",
       ylim = range(c(datalist_cali[[iREG]]$mort, datalist_cali[[iREG]]$mort_pred_seas8), na.rm = TRUE), 
       main = paste0(vREG_NAME[iREG], " (", vREG[iREG], ")" ),
       xlab = "Time (year)", ylab = "Mortality (deaths)")
  
  # Add red line showing seasonal component of mortality
  lines(datalist_cali[[iREG]]$date, datalist_cali[[iREG]]$mort_pred_seas8, 
        col = rgb(1.00, 0.00, 0.00), lwd = 1.5, lty = 1)
  
  # Highlight calibration period as shaded area
  polygon(c(date1_cali, date2_cali, date2_cali, date1_cali), 
          10^6 * c(-1, -1, +1, +1), col = rgb(1, 1, 1, 0.1), border = FALSE)
  
  # Horizontal line at zero
  abline(h = 0, col = rgb(0, 0, 0), lwd = 1, lty = 1)
  
  # Add legend only for the first plot
  if (iREG == 1) {
    legend("top", c("Mortality", "Seasonality"), 
           col = c("black", "red"), lwd = 2, lty = 1, box.lty = 0, horiz = TRUE, bg = "transparent")
  }
  
  # Customize x-axis with yearly ticks
  axis(1, at = seq(min(datalist_cali[[iREG]]$date), max(datalist_cali[[iREG]]$date), by = "years"),
       labels = format(seq(min(datalist_cali[[iREG]]$date), max(datalist_cali[[iREG]]$date), by = "years"), "%Y"))
}

dev.off()

##################################################################################
########### Figure: Cumulative Exposure-Response of Individual Regions 
########### Before and After the Meta-Analysis
##################################################################################

# This section plots the cumulative exposure-response curves for temperature vs mortality
# Both before and after pooling data with a meta-analysis.

print("")
print("= Cumulative Exposure-Response of Individual Regions =")
print("")

# Create output folder for cumulative exposure-response plots
foldout = paste0("./dataout/CumulativeExpResp/")
if(!file_test("-d", foldout)){ dir.create(foldout, recursive = TRUE) }

# -------------------------------
# BEFORE META-ANALYSIS (Individual Regions)
# -------------------------------
pdf(paste0(foldout, "PLOT_CUMU_ERF_REG_NOMETA.pdf"), width = 15, height = 21)
layout(matrix(seq(5*7), nrow = 7, byrow = TRUE))  # 7 rows, 5 columns
par(mex = 0.8, mgp = c(2.5,1,0), las = 0)

for(iREG in 1:nREG){
  plot(CROSS_PRED_REG_NOMETA[[iREG]]$predvar,
       CROSS_PRED_REG_NOMETA[[iREG]]$allRRfit, 
       col = rgb(0.00,0.00,0.00), lwd = 4, lty = 1, type = "l", 
       xlab = expression(paste("Temperature (", degree, "C)")),
       ylab = "Relative Risk", 
       main = paste0(vREG_NAME[iREG], " (", vREG[iREG], ")"), 
       ylim = c(1, max(CROSS_PRED_REG_NOMETA[[iREG]]$allRRfit)), axes = T)
  
  # Add shaded area for 95% CI
  polygon(c(CROSS_PRED_REG_NOMETA[[iREG]]$predvar, 
            rev(CROSS_PRED_REG_NOMETA[[iREG]]$predvar)),
          c(CROSS_PRED_REG_NOMETA[[iREG]]$allRRlow,
            rev(CROSS_PRED_REG_NOMETA[[iREG]]$allRRhigh)),
          col = rgb(0.00,0.00,0.00,1/3), border = FALSE)
  
  # Draw main relative risk line
  lines(CROSS_PRED_REG_NOMETA[[iREG]]$predvar,
        CROSS_PRED_REG_NOMETA[[iREG]]$allRRfit, 
        col = rgb(0.00,0.00,0.00), lwd = 4, lty = 1)
  
  # Reference line at RR = 1
  abline(h = 1, lwd = 1, lty = 1, col = rgb(0.00,0.00,0.00))
}
rm(iREG)
dev.off()


# -------------------------------
# AFTER META-ANALYSIS (Individual Regions)
# -------------------------------
pdf(paste0(foldout, "PLOT_CUMU_ERF_REG_META.pdf"), width = 15, height = 21)
layout(matrix(seq(5*7), nrow = 7, byrow = TRUE))
par(mex = 0.8, mgp = c(2.5,1,0), las = 0)

for(iREG in 1:nREG){
  plot(CROSS_PRED_REG_META[[iREG]]$predvar, 
       CROSS_PRED_REG_META[[iREG]]$allRRfit, 
       col = rgb(0.00,0.00,0.00), lwd = 4, lty = 1, type = "l",
       xlab = expression(paste("Temperature (", degree, "C)")),
       ylab = "Relative Risk", 
       main = paste0(vREG_NAME[iREG], " (", vREG[iREG], ")"),
       ylim = c(1, max(CROSS_PRED_REG_META[[iREG]]$allRRfit)), axes = T)
  
  # 95% CI shaded area
  polygon(c(CROSS_PRED_REG_META[[iREG]]$predvar, 
            rev(CROSS_PRED_REG_META[[iREG]]$predvar)), 
          c(CROSS_PRED_REG_META[[iREG]]$allRRlow, 
            rev(CROSS_PRED_REG_META[[iREG]]$allRRhigh)), 
          col = rgb(0.00,0.00,0.00,1/3), border = FALSE)
  
  # Main RR line
  lines(CROSS_PRED_REG_META[[iREG]]$predvar,
        CROSS_PRED_REG_META[[iREG]]$allRRfit, 
        col = rgb(0.00,0.00,0.00), lwd = 4, lty = 1)
  
  # Reference line at RR = 1
  abline(h = 1, lwd = 1, lty = 1, col = rgb(0.00,0.00,0.00))
}
rm(iREG)
dev.off()


# -------------------------------
# OVERALL CUMULATIVE EXPOSURE-RESPONSE (All Regions Combined)
# -------------------------------
pdf(paste0(foldout, "PLOT_CUMU_ERF_TOT_META.pdf"), width = 4, height = 4)
layout(matrix(seq(1*1), nrow = 1, byrow = TRUE))
par(mex = 0.8, mgp = c(2.5,1,0), las = 0)

plot(CROSS_PRED_TOT_META[[1]]$predvar, CROSS_PRED_TOT_META[[1]]$allRRfit,
     col = rgb(0.00,0.00,0.00), lwd = 4, lty = 1, type = "l",
     xlab = expression(paste("Temperature (", degree, "C)")), 
     ylab = "Relative Risk",
     main = "SP, SE, FR",
     ylim = c(1, max(CROSS_PRED_TOT_META[[1]]$allRRfit)), axes = T)

# Add shaded confidence interval
polygon(c(CROSS_PRED_TOT_META[[1]]$predvar, rev(CROSS_PRED_TOT_META[[1]]$predvar)), 
        c(CROSS_PRED_TOT_META[[1]]$allRRlow, rev(CROSS_PRED_TOT_META[[1]]$allRRhigh)),
        col = rgb(0.00,0.00,0.00,1/3), border = FALSE)

# Main RR line
lines(CROSS_PRED_TOT_META[[1]]$predvar, CROSS_PRED_TOT_META[[1]]$allRRfit, 
      col = rgb(0.00,0.00,0.00), lwd = 4, lty = 1)

# Mark Minimum Mortality Temperature (MMT)
abline(v = MMT_TOT, col = "black", lwd = 1, lty = 1)
mtext(paste("MMT =", round(MMT_TOT, 2)), side = 3, adj = 1, line = 1, col = "black", cex = 0.75)

dev.off()

# Clean up
rm(foldout)

##################################################################################
########### END OF PRACTICAL SESSION
##################################################################################

