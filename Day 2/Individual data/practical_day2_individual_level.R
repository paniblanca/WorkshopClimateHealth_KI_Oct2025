# ==============================================================================
# PRACTICAL SESSION: Modelling Associations Between Climatic Factors & Health 
# Outcomes Using Individual-Level Data
# ==============================================================================
# Author: Isabel Walter (isabel.walter@ki.se)
# Workshop: Population Health Impacts from Climate Extremes and Climatic Factors
# Date: 02-10-2025
# ==============================================================================
#
# LEARNING OBJECTIVES:
# 1. Explore an individual-level case time series dataset of asthma exacerbations
# 2. Learn how to model the association between daily mean temperature 
#    and asthma exacerbations using cross-basis functions
# 3. Study subgroup differences (e.g., smokers vs non-smokers)
# 4. Predict and interpret exposure–response curves
#
# ==============================================================================

# ==============================================================================
# SETUP AND DATA LOADING
# ==============================================================================

# ------------------------------------------------------------------------------
# CLEAN ENVIRONMENT AND LOAD PACKAGES
# ------------------------------------------------------------------------------
rm(list = ls())


# Install required packages (only run once if not installed)
packages <- c("dplyr", "gnm", "dlnm", "ggplot2", "lubridate")

installed <- packages %in% installed.packages()
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load required packages
library(dplyr)
library(gnm)
library(dlnm)
library(ggplot2)
library(lubridate)

# ------------------------------------------------------------------------------
# SET WORKING DIRECTORY AND LOAD DATA
# ------------------------------------------------------------------------------
# NOTE: Adjust the path to your own system
setwd("/Users/isabelwalter/Desktop/PhD/Workshop/Practicals/Day 2 - individual level")

df <- read.csv("data/time_series.csv")

# ============================================================================== 
# EXPLORING THE DATASET
# ==============================================================================

# -------------------------------
# Q1. How many rows and columns are in the dataset?
# -------------------------------

# -------------------------------
# Q2. How many unique individuals are there?
# -------------------------------

# -------------------------------
# Q3. What does the distribution of temperature (temp) look like?
# -------------------------------
ggplot(df, aes(x = temp)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Distribution of daily mean temperature")

summary(df$temp)

# -------------------------------
# Q4. How many days of follow-up are there on average per individual?
# -------------------------------

# -------------------------------
# Q5. On average, how many asthma exacerbations per person?
# -------------------------------

# -------------------------------
# Q6. How many individuals had no outcomes at all?
# Also, remove them from the dataset, as we want cases only.
# -------------------------------
outcomes_per_person <- df %>%
  group_by(id) %>%
  summarise(total_outcomes = sum(outcome))

# Count individuals with zero outcomes
sum(outcomes_per_person$total_outcomes == 0)

# Get their IDs
zero_outcome_ids <- outcomes_per_person %>%
  filter(total_outcomes == 0) %>%
  pull(id)

# Remove those individuals from the dataset
df_clean <- df %>%
  filter(!id %in% zero_outcome_ids)

# Check: number of rows before vs after
nrow(df)
nrow(df_clean)

df <- df_clean

# -------------------------------
# Q7. How many individuals are smokers (ever smoking == 1)?
# -------------------------------
smokers <- df %>%
  group_by(id) %>%
  summarise(is_smoker = any(smoking == 1))

sum(smokers$is_smoker == TRUE)

# -------------------------------
# Q8. Do smokers have more exacerbations on average?
# -------------------------------

# -------------------------------
# Q9. Explore the data for one individual:
#   (a) Plot the time series of asthma exacerbations.
#   (b) Plot the time series of mean daily air temperature.
#   (c) Is the number of asthma exacerbations realistic?
#   (d) Does there seem to be a correlation between 
#       temperature and asthma exacerbations?
# -------------------------------
person_id <- 5
df_person <- df[df$id == person_id, ]
df_person$date <- as.Date(df_person$date)

# (a) Exacerbations: vertical segments ending at the point
ggplot(subset(df_person, outcome == 1), aes(x = date, y = 1)) +
  geom_segment(aes(xend = date, yend = 0), linetype = "dashed", color = "black") +
  geom_point(color = "darkgreen", size = 2) +
  scale_y_continuous(limits = c(0, 1.2), breaks = 1, labels = "Day with Event") +
  scale_x_date(date_labels = "%b") +
  labs(title = paste("Time series of Exacerbations (Person", person_id, ")"),
       x = "", y = "") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank())

# (b) Temperature: solid axes, same x-axis style
ggplot(df_person, aes(x = date, y = temp)) +
  geom_line(color = "darkgreen") +
  scale_x_date(date_labels = "%b") +
  labs(title = paste("Mean Daily Temperature (Person", person_id, " area)"),
       x = "", y = "Temperature (°C)") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank())

# Plot both in one figure
ggplot(df_person, aes(x = date)) +
  # Temperature line
  geom_line(aes(y = temp), color = "darkgreen") +
  # Exacerbation events as points
  geom_point(data = subset(df_person, outcome == 1),
             aes(y = min(temp) - 2),
             color = "firebrick", size = 2) +
  geom_segment(data = subset(df_person, outcome == 1),
               aes(xend = date, y = min(temp) - 2, yend = min(temp) - 1.5),
               color = "firebrick", linetype = "dashed") +
  scale_x_date(date_labels = "%b") +
  labs(title = paste("Temperature and Exacerbations (Person", person_id, ")"),
       x = "", y = "Temperature (°C)") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.minor = element_blank())

# ==============================================================================
# DATA MODELLING
# ==============================================================================

# -------------------------------
# Q10. What are knots, and why do we need them when modelling temperature?
# -------------------------------

# -------------------------------
# Q11. How do we define the knots and centering values for temperature?
# -------------------------------

# In this dataset we define:
# - Internal knots at the 10th, 75th, and 90th percentiles
# - Boundary knots at the 1st and 99th percentiles
# - A centering value at the 50th percentile (median)

temp_knots <- c(-2.510472, 16.020317, 19.799403)  # 10%, 75%, 90%
temp_b_knots <- c(-8.547488, 25.588451)           # 1%, 99%
cen <- 6.914142                                   # 50%

# -------------------------------
# Q12. What is the stratum variable, and why do we need it?
# -------------------------------
df$stratum <- factor(paste(df$id, month(df$date), sep = "-"))

# -------------------------------
# Q13. Which two dimensions are combined in a cross-basis?
# -------------------------------

# -------------------------------
# Q14. Create the crossbasis using the code in the script. What arguments do 
# we hand crossbasis()? 
# -------------------------------
cbtemp <- crossbasis(df$temp, lag = 21,
                      argvar = list(fun = "ns", 
                                    knots = temp_knots, 
                                    Boundary.knots = temp_b_knots),
                      arglag = list(fun = "ns", 
                                    knots = logknots(21, nk = 3)),
                      group = df$id)

# -------------------------------
# Q15. Fit a poisson regression model with the crossbasis as exposure, 
# day of week as a covariate, and the stratum variable as fixed effects. 
# -------------------------------
day_of_week <- factor(wday(df$date))
#mod <- gnm(..)

# ==============================================================================
# PART D: CROSSPREDICTION
# ==============================================================================

# -------------------------------
# Q16. Predict and plot the exposure–response curve  
# for temperature from our model using the provided code.
# -------------------------------

# Create cross-prediction object
cptemp <- crosspred(cbtemp, mod, cen = cen, by = 1.5)

# Convert to a tidy dataframe
plot_data <- data.frame(
  exposure = cptemp$predvar,
  fit      = cptemp$allRRfit,
  lower    = cptemp$allRRlow,
  upper    = cptemp$allRRhigh
)

# Plot using ggplot
ggplot(plot_data, aes(x = exposure, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.2, fill = "#1b9e77", color = NA) +
  geom_line(linewidth = 1, color = "#1b9e77") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = cen, linewidth = 1.2, alpha = 0.5, color = "grey") +
  annotate("text", x = cen + 2, y = max(plot_data$fit, na.rm = TRUE) * 0.9,
           label = "50th percentile", color = "grey20", size = 3) +
  labs(title = "Exposure-response curve: Asthma exacerbations",
       x = "Temperature (°C)",
       y = "Incidence rate ratio (IRR)") +
  scale_y_continuous(breaks = seq(0, 5, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-10, 25, by = 5), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-10, 25), ylim = c(0, 5)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(color = "white", fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

# -------------------------------
# Q17. How do we interpret the exposure–response curve?
# -------------------------------

# -------------------------------
# Q18. Compare the temperature–asthma exacerbation association 
# between smokers and non-smokers by running the provided code. Interpret 
# the plot.
# -------------------------------

# Step 1. Subset the dataset
# - Smokers: any day with smoking == 1 for that person
# - Non-smokers: all other individuals
df_smok <- df %>% 
  group_by(id) %>% 
  filter(any(smoking == 1)) %>% 
  ungroup()

df_non <- df %>% 
  group_by(id) %>% 
  filter(all(smoking == 0)) %>% 
  ungroup()

# Step 2. Create strata variables (subject-month) for each subset
df_smok$stratum <- factor(paste(df_smok$id, month(df_smok$date), sep = "-"))
df_non$stratum  <- factor(paste(df_non$id,  month(df_non$date),  sep = "-"))

# Step 3. Build cross-basis for temperature
cbtemp_smok <- crossbasis(df_smok$temp, lag = 21,
                           argvar = list(fun = "ns", knots = temp_knots, Boundary.knots = temp_b_knots),
                           arglag = list(fun = "ns", knots = logknots(21, nk = 3)),
                           group = df_smok$id)

cbtemp_non <- crossbasis(df_non$temp, lag = 21,
                          argvar = list(fun = "ns", knots = temp_knots, Boundary.knots = temp_b_knots),
                          arglag = list(fun = "ns", knots = logknots(21, nk = 3)),
                          group = df_non$id)

# Step 4. Fit conditional Poisson regression models
mod_smok <- gnm(outcome ~ cbtemp_smok + factor(wday(df_smok$date)),
                eliminate = df_smok$stratum, data = df_smok,
                family = poisson)

mod_non <- gnm(outcome ~ cbtemp_non + factor(wday(df_non$date)),
               eliminate = df_non$stratum, data = df_non,
               family = poisson)

# Step 5. Cross-predictions
cpt_smok <- crosspred(cbtemp_smok, mod_smok, cen = cen, by = 1.5)
cpt_non  <- crosspred(cbtemp_non,  mod_non,  cen = cen, by = 1.5)

# Step 6. Tidy for plotting
plot_smok <- data.frame(
  exposure = cpt_smok$predvar,
  fit      = cpt_smok$allRRfit,
  lower    = cpt_smok$allRRlow,
  upper    = cpt_smok$allRRhigh,
  group    = "Smoker"
)

plot_non <- data.frame(
  exposure = cpt_non$predvar,
  fit      = cpt_non$allRRfit,
  lower    = cpt_non$allRRlow,
  upper    = cpt_non$allRRhigh,
  group    = "Non-smoker"
)

plot_data <- rbind(plot_smok, plot_non)

# Step 7. Plot both curves together
ggplot(plot_data, aes(x = exposure, y = fit, 
                      color = group, fill = group, linetype = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = cen, linewidth = 1.2, alpha = 0.5, color = "grey") +
  annotate("text", x = cen + 2.5, y = 3.5, label = "50th percentile",
           color = "grey20", size = 3) +
  labs(title = "Exposure–response curves: Asthma exacerbations",
       x = "Temperature (°C)",
       y = "Incidence rate ratio (IRR)") +
  scale_color_manual(values = c("Non-smoker" = "#1b9e77", "Smoker" = "#d95f02")) +
  scale_fill_manual(values = c("Non-smoker" = "#1b9e77", "Smoker" = "#d95f02")) +
  scale_linetype_manual(values = c("Non-smoker" = "solid", "Smoker" = "dashed")) +
  scale_y_continuous(breaks = seq(0, 5, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-10, 25, by = 5), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-10, 25), ylim = c(0, 5)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.position = "right",
    legend.title = element_blank()
  )

# ==============================================================================
# END OF PRACTICAL SESSION
# ==============================================================================