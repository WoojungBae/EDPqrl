# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# setwd("/Users/woojung/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA")
# load(paste0("~/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNImcmc.RData"))
# load("~/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNImcmc+post.RData")

# setwd("/Users/woojung/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNI")
# install.packages("ADNIMERGE2.tar.gz", repos = NULL, type = "source")
library(ADNIMERGE2)
library(Rcpp)
library(RcppArmadillo)

library(survival)

library(dplyr)
library(tidyverse)

library(viridis)
library(akima)
library(scales)
library(survminer)
library(grid)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
setwd("/Users/WooJung/Documents/Rproject/EDPqrl/source")

# Load R code
source("EDPqrl_r.R")

# Load cpp code
sourceCpp("EDPqrl_cpp.cpp")

# seed
set.seed(1)

# ------------------------------------------------------------------------
# Load Data --------------------------------------------------------------
# ------------------------------------------------------------------------
setwd("/Users/woojung/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA")

# ?ADAE
# ?ADQS
# ?ADRS
# ?ADSL
# colnames(ADAE)
# colnames(ADQS)
# colnames(ADRS)
# colnames(ADSL)

# Time to dimentia (survival)
# ADSL; ENRLDT Date Date of Enrollment Date: YYY-MM-DD
# ADSL; EOSDT Date End of Study Date Date: YYY-MM-DD
# ADSL; EOSDTF character End of Study Date Imputation Flag , with character options: M
# ADRS; ADY numeric Analysis Relative Day Range value: -1849, ..., 7208
# ADRS; ADT Date Analysis Date Date: YYY-MM-DD
survival_data <- ADRS %>%
  filter(PARAMCD == "DX") %>%
  # filter(COLPROT != "ADNI4") %>%
  arrange(USUBJID, ADT) %>%
  group_by(USUBJID) %>%
  summarise(
    Baseline_DX = first(AVALC),
    ENRLDT = first(ENRLDT),
    # If all diagnoses are missing, return NA. Otherwise, check for "DEM".
    Status = if(all(is.na(AVALC))) {
      NA_real_
    } else {
      as.numeric(any(AVALC == "DEM", na.rm = TRUE))
    },
    # Safely calculate DXDT (Date)
    DXDT = if(first(AVALC) %in% "DEM") {
      first(ADT)
    } else if(any(AVALC == "DEM", na.rm = TRUE)) {
      dem_dates <- ADT[which(AVALC == "DEM")]
      if(any(!is.na(dem_dates))) {
        min(dem_dates, na.rm = TRUE)
      } else if(any(!is.na(ADT))) {
        max(ADT, na.rm = TRUE)
      } else {
        as.Date(NA) # Returns proper NA for Date columns
      }
    } else {
      if(any(!is.na(ADT))) max(ADT, na.rm = TRUE) else as.Date(NA)
    },
    
    # # Safely calculate T_days (Numeric)
    # T_days = pmax(
    #   if(first(AVALC) %in% "DEM") {
    #     0
    #   } else if(any(AVALC == "DEM", na.rm = TRUE)) {
    #     dem_days <- ADY[which(AVALC == "DEM")]
    #     if(any(!is.na(dem_days))) {
    #       min(dem_days, na.rm = TRUE)
    #     } else if(any(!is.na(ADY))) {
    #       max(ADY, na.rm = TRUE)
    #     } else {
    #       NA_real_ 
    #     }
    #   } else {
    #     if(any(!is.na(ADY))) max(ADY, na.rm = TRUE) else NA_real_
    #   }, 
    #   0
    # ),
    # Safely calculate T_days (Leaves NAs intact)
    T_days = if(first(AVALC) %in% "DEM") {
      0
    } else {
      pmax(as.numeric(difftime(DXDT, ENRLDT, units = "days")), 0)
    },
    .groups = "drop"
  )

# AGE numeric Age (in Years) Range value: 50.5, ..., 91.4
# SEX character Sex , with character options: Male, Female
# EDUC numeric Education Range value: 4, ..., 20
# APOE character APOE Genotype , with character options: ε4/ε4, ε3/ε3, ε3/ε4, ε2/ε3, ε2/ε2, ε2/ε4
# AMYSTAT character Baseline Amyloid Status , with character options: Elevated, Non Elevated
adsl_analysis <- ADSL %>%
  dplyr::select(USUBJID, AGE, SEX, EDUC, APOE, AMYSTAT) %>%
  dplyr::mutate(
    Amyloid_Positive = dplyr::case_when(
      AMYSTAT == "Elevated" ~ 1,
      AMYSTAT == "Non Elevated" ~ 0,
      TRUE ~ NA_real_
    ),
    APOE4_Carrier = dplyr::case_when(
      is.na(APOE) ~ NA_real_,        # If missing, keep it as NA
      grepl("4", APOE) ~ 1,          # If it contains a 4, they are a carrier
      TRUE ~ 0                       # Otherwise, they are not
    ),
    APOE4_Alleles = stringr::str_count(APOE, "4")
  )

# Re-run the merge, sort, and reorder columns
merge_data <- survival_data %>%
  dplyr::left_join(adsl_analysis, by = "USUBJID") %>%
  dplyr::mutate(
    T_years = T_days/365.25
  ) %>%
  # Sort alphabetically/numerically by Subject ID
  dplyr::arrange(USUBJID) %>%
  # Explicitly define the exact order of columns you want
  dplyr::select(
    USUBJID, 
    ENRLDT, DXDT, T_years, T_days, 
    Status, Baseline_DX, Amyloid_Positive, 
    APOE, APOE4_Carrier, APOE4_Alleles,
    AGE, SEX, EDUC
  )

exclusion_summary <- merge_data %>%
  dplyr::summarise(
    # --- 1. Individual Exclusion Counts ---
    Missing_Amyloid = sum(is.na(Amyloid_Positive)),
    Missing_APOE = sum(is.na(APOE)),
    Baseline_DEM = sum(Baseline_DX == "DEM", na.rm = TRUE),
    Time_is_Zero_or_Less = sum(T_years <= 0, na.rm = TRUE),
    Time_is_NA = sum(is.na(T_years)),
    # --- 2. Specific Intersections ---
    # Missing Exposure AND Baseline DEM
    Missing_Amyloid_and_Base_DEM = sum(is.na(Amyloid_Positive) & Baseline_DX == "DEM", na.rm = TRUE),
    # Missing Confounder AND Baseline DEM
    Missing_APOE_and_Base_DEM = sum(is.na(APOE) & Baseline_DX == "DEM", na.rm = TRUE),
    # Missing BOTH Amyloid and APOE
    Missing_Both_Amyloid_APOE = sum(is.na(Amyloid_Positive) & is.na(APOE), na.rm = TRUE),
    # Baseline DEM, BUT they have valid follow-up time (T_years > 0)
    Base_DEM_with_Valid_Time = sum(Baseline_DX == "DEM" & T_years > 0, na.rm = TRUE),
    # --- 3. Total Unique Subjects Dropped ---
    # Anyone missing the exposure, missing the confounder, starting with dementia, 
    # or having invalid follow-up time gets dropped from the final CN/MCI analytic cohort
    Total_Excluded_Subjects = sum(
      is.na(Amyloid_Positive) | is.na(APOE) | Baseline_DX == "DEM" | T_years <= 0 | is.na(T_years), 
      na.rm = TRUE
    )
  )
t(exclusion_summary)

final_data <- merge_data %>%
  dplyr::filter(
    # Exclude missing time
    !is.na(T_years),
    # Exclude T_years = 0 (and negative screening times)
    T_years > 0,
    # Exclude Baseline Dementia (Keep only CN and MCI, and implicitly drop NAs in diagnosis)
    Baseline_DX != "DEM",
    # Baseline_DX == "MCI",
    # EXCLUDE missing Amyloid (because it's the primary exposure Z)
    !is.na(Amyloid_Positive)
  ) %>%
  dplyr::mutate(
    SEX = ifelse(SEX == "Female", 1, 0)
  ) %>%
  # Explicitly define the exact order of columns you want
  dplyr::select(
    T_years, Status, 
    Amyloid_Positive, 
    APOE4_Carrier, Baseline_DX, SEX, AGE, EDUC
    # USUBJID, 
    # ENRLDT, DXDT, T_years, T_days, 
    # Status, Baseline_DX, Amyloid_Positive, 
    # APOE, APOE4_Carrier, APOE4_Alleles,
    # AGE, SEX, EDUC
  )

# Baseline summary
# 1) Final N
N_final <- nrow(final_data)
# 2) Event count and percent (assumes Status is 1=event, 0=censored)
n_event <- sum(final_data$Status == 1, na.rm = TRUE)
pct_event <- 100 * n_event / N_final
# 3) NEW PRIMARY EXPOSURE: Baseline Amyloid Positivity ---
n_amyloid <- sum(final_data$Amyloid_Positive == 1, na.rm = TRUE)
pct_amyloid <- 100 * n_amyloid / N_final
# 4) KEY CONFOUNDER: APOE ε4 carrier count and percent ---
n_carrier <- sum(final_data$APOE4_Carrier == 1, na.rm = TRUE)
pct_carrier <- 100 * n_carrier / N_final
# 5) Baseline covariate summaries
age_mean <- mean(final_data$AGE, na.rm = TRUE)
age_sd   <- sd(final_data$AGE, na.rm = TRUE)
# ADNI ADSL dataset uses "Female" and "Male" strings
pct_female <- 100 * mean(final_data$SEX, na.rm = TRUE)
educ_mean <- mean(final_data$EDUC, na.rm = TRUE)
educ_sd   <- sd(final_data$EDUC, na.rm = TRUE)
# 6) Missingness (Double-checking your final analytic cohort)
miss_tbl <- final_data %>%
  dplyr::summarise(
    Missing_Baseline_DX = sum(is.na(Baseline_DX)),
    Missing_Amyloid     = sum(is.na(Amyloid_Positive)),
    Missing_APOE        = sum(is.na(APOE4_Carrier)),
    Missing_AGE         = sum(is.na(AGE)),
    Missing_SEX         = sum(is.na(SEX)),
    Missing_EDUC        = sum(is.na(EDUC))
  )
list(
  N_final = N_final,
  events = c(n_event = n_event, pct_event = pct_event),
  amyloid_positive = c(n_amyloid = n_amyloid, pct_amyloid = pct_amyloid),
  apoe4_carrier = c(n_carrier = n_carrier, pct_carrier = pct_carrier),
  age = c(mean = age_mean, sd = age_sd),
  female = c(pct_female = pct_female),
  educ = c(mean = educ_mean, sd = educ_sd),
  missing = miss_tbl
)

dfADNI <- final_data %>%
  dplyr::transmute(
    Tobs = T_years,
    Dobs = Status,
    Zobs = Amyloid_Positive,
    Xobs = cbind(
      Baseline_DX = ifelse(Baseline_DX == "MCI", 1, 0),
      APOE4_Carrier,
      SEX = SEX,
      AGE = AGE,
      EDUC = EDUC
    )
  )

km_fit <- survfit(Surv(Tobs, Dobs) ~ Zobs, data = dfADNI)
km_plot <- ggsurvplot(
  km_fit, 
  data = dfADNI,
  pval = TRUE,                
  conf.int = TRUE, 
  conf.int.alpha = 0.15,          
  # --- 1. Risk Table Updates ---
  risk.table = "nrisk_cumevents", 
  risk.table.height = 0.25,  
  risk.table.col = "black",       
  # --- 2. Color & Linetype Updates ---
  linetype = c("solid", "dashed"), 
  palette = c("red", "blue"),      # <--- Updated to Red (Z=0) and Blue (Z=1)
  linewidth = 1.2,                 
  # --- 3. Fix the Censoring Marks ---
  censor.shape = 124,              
  censor.size = 3,                 
  # --- Center the Title Here ---
  ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5)),      
  # --- Labels ---
  title = "Kaplan-Meier Survival Curve by Baseline Amyloid Status",
  xlab = "Time (Years)",
  ylab = "Probability of Remaining Dementia-Free",
  legend.title = "Amyloid Status",
  legend.labs = c("Non-Elevated (Z=0)", "Elevated (Z=1)")
)
# print(km_plot)

# quartz(type = "pdf", file = "EDPqrl_KM.pdf", width = 8, height = 6)
# print(km_plot)
# dev.off()

# 1. Create the subgroup data using Zobs and Dobs
df_subgroups <- dfADNI %>%
  mutate(
    Group = case_when(
      Zobs == 0 & Dobs == 0 ~ "Z=0, D=0 (Non-Elevated, Censored)",
      Zobs == 1 & Dobs == 0 ~ "Z=1, D=0 (Elevated, Censored)",
      Zobs == 0 & Dobs == 1 ~ "Z=0, D=1 (Non-Elevated, Event)",
      Zobs == 1 & Dobs == 1 ~ "Z=1, D=1 (Elevated, Event)",
      TRUE ~ NA_character_ 
    ),
    log_Tobs = log(Tobs) 
  ) %>%
  filter(!is.na(Group)) 

# 2. Create the "Overall" data
df_overall <- dfADNI %>%
  mutate(
    Group = "Overall",
    log_Tobs = log(Tobs)
  )

# 3. Stack them together into the final df_plot
df_plot <- bind_rows(df_overall, df_subgroups)

# 2. Generate the plot
EDPqrl_DENSITYlogT <- ggplot(df_plot, aes(x = log_Tobs, color = Group, linetype = Group)) +
  # --- 1. Fix the Legend Squares ---
  geom_density(linewidth = 1.2, key_glyph = "path") + 
  coord_cartesian(xlim = c(-6, 3), ylim = c(0, 0.5)) +
  # --- 2. Color Updates (Red/Blue Palette) & Labels ---
  scale_color_manual(
    values = c(
      "Overall" = "black", 
      "Z=0, D=0 (Non-Elevated, Censored)" = "lightcoral",  # Light Red
      "Z=1, D=0 (Elevated, Censored)" = "skyblue",       # Light Blue
      "Z=0, D=1 (Non-Elevated, Event)" = "red",          # Solid Red
      "Z=1, D=1 (Elevated, Event)" = "blue"              # Solid Blue
    ),
    labels = c(
      "Overall" = "Overall", 
      "Z=0, D=0 (Non-Elevated, Censored)" = "Non-Elevated, Censored (Z=0, D=0)", 
      "Z=1, D=0 (Elevated, Censored)" = "Elevated, Censored (Z=1, D=0)", 
      "Z=0, D=1 (Non-Elevated, Event)" = "Non-Elevated, Event (Z=0, D=1)", 
      "Z=1, D=1 (Elevated, Event)" = "Elevated, Event (Z=1, D=1)"
    )
  ) +
  # --- Distinct Linetypes & Updated Labels ---
  scale_linetype_manual(
    values = c(
      "Overall" = "solid", 
      "Z=0, D=0 (Non-Elevated, Censored)" = "dashed", 
      "Z=1, D=0 (Elevated, Censored)" = "dashed", 
      "Z=0, D=1 (Non-Elevated, Event)" = "dotted", 
      "Z=1, D=1 (Elevated, Event)" = "dotted"
    ),
    labels = c(
      "Overall" = "Overall", 
      "Z=0, D=0 (Non-Elevated, Censored)" = "Non-Elevated, Censored (Z=0, D=0)", 
      "Z=1, D=0 (Elevated, Censored)" = "Elevated, Censored (Z=1, D=0)", 
      "Z=0, D=1 (Non-Elevated, Event)" = "Non-Elevated, Event (Z=0, D=1)", 
      "Z=1, D=1 (Elevated, Event)" = "Elevated, Event (Z=1, D=1)"
    )
  ) +
  theme_classic() + 
  # --- Labels ---
  labs(
    title = "Density of Log-Transformed Observed Survival Time",
    x = "Log-Transformed Time to Dementia", 
    y = "Density",
    color = "Subgroup",
    linetype = "Subgroup" 
  ) +
  # --- 3. Legend Placement (Top-Left) ---
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.02, 0.98), 
    legend.justification = c(0, 1),    
    legend.background = element_rect(fill = alpha("white", 0.85), color = "black", linewidth = 0.5),
    legend.key.width = unit(3, "lines"), 
    legend.text = element_text(size = 9)
  ) +
  guides(
    color = guide_legend(ncol = 1),
    linetype = guide_legend(ncol = 1)
  )
# print(EDPqrl_DENSITYlogT)

# quartz(type = "pdf", file = "EDPqrl_DENSITYlogT.pdf", width = 8, height = 6)
# print(EDPqrl_DENSITYlogT)
# dev.off()

# ------------------------------------------------------------------------
# Constants --------------------------------------------------------------
# ------------------------------------------------------------------------
gibbs_thin = 5e1
gibbs_iter = 2e4
gibbs_burnin = 2e4

# # ------------------------------------------------------------------------
# # MCMC sampling ----------------------------------------------------------
# # ------------------------------------------------------------------------
# # EDPMM
# EDPMM_results_MCMC = SURVIVAL_EDPMM_MCMC(dfADNI, gibbs_iter, gibbs_burnin, gibbs_thin)
# # DPMM
# DPMM_results_MCMC = SURVIVAL_DPMM_MCMC(dfADNI, gibbs_iter, gibbs_burnin, gibbs_thin)
# 
# save.image(paste0("~/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNImcmc.RData"))

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
load(paste0("~/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNImcmc.RData"))
setwd("/Users/WooJung/Documents/Rproject/EDPqrl/source")
source("EDPqrl_r.R")
sourceCpp("EDPqrl_cpp.cpp")
setwd("/Users/woojung/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA")

# ------------------------------------------------------------------------
# Average Effects & Sensitivity Analysis ---------------------------------
# ------------------------------------------------------------------------
M = 2e3

# Sensitivity Analysis
eta_z0 = c(-0.25, 0, 0.25)
eta_z1 = c(-0.25, 0, 0.25)

# censoring = 1 - mean(dfADNI$Dobs) = 85.83397
# summary(dfADNI$Tobs)
# min:  0.0027
# max: 14.8419
# t_grid = seq(1e-3, 20, by = 0.1)  => extrapolation from 14.9 - 20.0
# nu = c(1:200)/10                  => extrapolation from 14.9 - 20.0
# rho = c(1:30)/100                 => extrapolation from 0.15 - 0.20
nu = seq(0, 20, by = 0.1)
rho = seq(0.01, 0.30, by = 0.01)
t_grid = seq(1e-3, 20, by = 0.1)

# EDPMM
EDPMM_results_POST =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC, M,
                               nu, rho, eta_z0, eta_z1)
# EDPMM_results_POST =
#   SURVIVAL_EDPMM_POST(EDPMM_results_MCMC, M,
#                       nu, rho, eta_z0, eta_z1)
# DPMM
DPMM_results_POST =
  SURVIVAL_DPMM_POST_parallel(DPMM_results_MCMC, M,
                              nu, rho, eta_z0, eta_z1)
# DPMM_results_POST =
#   SURVIVAL_DPMM_POST(DPMM_results_MCMC, M,
#                      nu, rho, eta_z0, eta_z1)
# print(EDPMM_results_POST$E_rho_nu_diff_result, max=1e5)
# print(DPMM_results_POST$E_rho_nu_diff_result, max=1e5)

SURVIVAL_EDPMM_POST_ftn_results =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                          condX = NULL, p_condX = NULL, esttype = "mean")
SURVIVAL_DPMM_POST_ftn_results =
  SURVIVAL_DPMM_POST_ftn(DPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                         condX = NULL, p_condX = NULL, esttype = "mean")

# SURVIVAL_EDPMM_POST_ftn_results$S_z0_result
# SURVIVAL_EDPMM_POST_ftn_results$S_z1_result
# SURVIVAL_EDPMM_POST_ftn_results$S_z0_store
# SURVIVAL_EDPMM_POST_ftn_results$S_z1_store

# SURVIVAL_DPMM_POST_ftn_results$S_z0_result
# SURVIVAL_DPMM_POST_ftn_results$S_z1_result
# SURVIVAL_DPMM_POST_ftn_results$S_z0_store
# SURVIVAL_DPMM_POST_ftn_results$S_z1_store

save.image(paste0("~/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNImcmc+post.RData"))

# ------------------------------------------------------------------------
# Conditional Effects & Sensitivity Analysis -----------------------------
# ------------------------------------------------------------------------
# Baseline_DX = 1 # Baseline_DX = ifelse(Baseline_DX == "MCI", 1, 0)
# EDPMM
# SURVIVAL_EDPMM_POST_condX1 =
#   SURVIVAL_EDPMM_POST(EDPMM_results_MCMC, M,
#                       nu, rho, eta_z0, eta_z1,
#                       condX = 1, p_condX = 1)
# SURVIVAL_EDPMM_POST_condX0 =
#   SURVIVAL_EDPMM_POST(EDPMM_results_MCMC, M,
#                      nu, rho, eta_z0, eta_z1,
#                      condX = 0, p_condX = 1)
SURVIVAL_EDPMM_POST_condX1 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC, M,
                               nu, rho, eta_z0, eta_z1,
                               condX = 1, p_condX = 1)
SURVIVAL_EDPMM_POST_condX0 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC, M,
                               nu, rho, eta_z0, eta_z1,
                               condX = 0, p_condX = 1)
# SURVIVAL_EDPMM_POST_condX1$E_rho_nu_diff_result
# SURVIVAL_EDPMM_POST_condX0$E_rho_nu_diff_result

# DPMM
# SURVIVAL_DPMM_POST_condX1 =
#   SURVIVAL_DPMM_POST(DPMM_results_MCMC, M,
#                      nu, rho, eta_z0, eta_z1,
#                      condX = 1, p_condX = 1)
# SURVIVAL_DPMM_POST_condX0 =
#   SURVIVAL_DPMM_POST(DPMM_results_MCMC, M,
#                      nu, rho, eta_z0, eta_z1,
#                      condX = 0, p_condX = 1)
SURVIVAL_DPMM_POST_condX1 =
  SURVIVAL_DPMM_POST_parallel(DPMM_results_MCMC, M,
                              nu, rho, eta_z0, eta_z1,
                              condX = 1, p_condX = 1)
SURVIVAL_DPMM_POST_condX0 =
  SURVIVAL_DPMM_POST_parallel(DPMM_results_MCMC, M,
                              nu, rho, eta_z0, eta_z1,
                              condX = 0, p_condX = 1)
# SURVIVAL_DPMM_POST_condX1$E_rho_nu_diff_result
# SURVIVAL_DPMM_POST_condX0$E_rho_nu_diff_result

# summary(dfADNI$Tobs)
# min:  0.0027
# max: 14.8419
SURVIVAL_EDPMM_POST_ftn_condX1_results =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                          condX = 1, p_condX = 1, esttype = "mean")
SURVIVAL_EDPMM_POST_ftn_condX0_results =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                          condX = 0, p_condX = 1, esttype = "mean")
SURVIVAL_DPMM_POST_ftn_condX1_results =
  SURVIVAL_DPMM_POST_ftn(DPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                         condX = 1, p_condX = 1, esttype = "mean")
SURVIVAL_DPMM_POST_ftn_condX0_results =
  SURVIVAL_DPMM_POST_ftn(DPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                         condX = 0, p_condX = 1, esttype = "mean")

# SURVIVAL_EDPMM_POST_ftn_condX1_results$S_z0_result
# SURVIVAL_EDPMM_POST_ftn_condX1_results$S_z1_result
# SURVIVAL_EDPMM_POST_ftn_condX1_results$S_z0_store
# SURVIVAL_EDPMM_POST_ftn_condX1_results$S_z1_store
# SURVIVAL_EDPMM_POST_ftn_condX0_results$S_z0_result
# SURVIVAL_EDPMM_POST_ftn_condX0_results$S_z1_result
# SURVIVAL_EDPMM_POST_ftn_condX0_results$S_z0_store
# SURVIVAL_EDPMM_POST_ftn_condX0_results$S_z1_store
#
# SURVIVAL_DPMM_POST_ftn_condX1_results$S_z0_result
# SURVIVAL_DPMM_POST_ftn_condX1_results$S_z1_result
# SURVIVAL_DPMM_POST_ftn_condX1_results$S_z0_store
# SURVIVAL_DPMM_POST_ftn_condX1_results$S_z1_store
# SURVIVAL_DPMM_POST_ftn_condX0_results$S_z0_result
# SURVIVAL_DPMM_POST_ftn_condX0_results$S_z1_result
# SURVIVAL_DPMM_POST_ftn_condX0_results$S_z0_store
# SURVIVAL_DPMM_POST_ftn_condX0_results$S_z1_store

save.image(paste0("~/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNImcmc+post.RData"))

# ------------------------------------------------------------------------
# Average Effects - Figure -----------------------------------------------
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
setwd("/Users/woojung/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA")
load("~/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNImcmc+post.RData")
setwd("/Users/WooJung/Documents/Rproject/EDPqrl/source")
source("EDPqrl_r.R")
sourceCpp("EDPqrl_cpp.cpp")
setwd("/Users/woojung/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA")

nu_keep <- c(0, 5, 10)
rho_keep <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

# EDPMM OSQC + sensitivity
df_ed <- as.data.frame(EDPMM_results_POST$E_rho_nu_diff_result) %>%
  mutate(
    sens = case_when(
      eta_z0 ==  0.0 & eta_z1 ==  0.0 ~ "psi0=0, psi1=0",
      eta_z0 == -0.1 & eta_z1 ==  0.1 ~ "psi0=-0.1, psi1=0.1",
      eta_z0 ==  0.1 & eta_z1 == -0.1 ~ "psi0=0.1, psi1=-0.1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sens))
df0 <- df_ed %>% filter(sens == "psi0=0, psi1=0")

# # OSQC surface (posterior mean) + contours
# # levels_osqc <- pretty(df0$estimYs, n = 8)
# EDPqrlPOST_OSQCsurface <- ggplot(df0, aes(x = nu, y = rho)) +
#   geom_raster(aes(fill = estimYs), interpolate = TRUE) +
#   # geom_contour(aes(z = estimYs), breaks = levels_osqc,
#   #              color = "white", alpha = 0.6) +
#   scale_fill_viridis_c(option = "D", name = "OSQC") +
#   labs(x = expression(nu), y = expression(rho),
#        title = "EDPMM OSQC surface",
#        subtitle = "Posterior mean") +
#   theme_bw(base_size = 12) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5)
#   )
# # quartz(type = "pdf", file = "EDPqrlPOST_OSQCsurface.pdf", width = 8, height = 6)
# # print(EDPqrlPOST_OSQCsurface)
# # dev.off()

# Sensitivity comparison (facet by sens) with viridis
df_ed <- df_ed %>%
  mutate(sens = factor(sens, levels = c("psi0=0, psi1=0",
                                        "psi0=-0.1, psi1=0.1",
                                        "psi0=0.1, psi1=-0.1")))
# levels_osqc <- pretty(df_ed$estimYs, n = 8)
EDPqrlPOST_OSQCsurface <- 
  ggplot(df_ed, aes(x = nu, y = rho)) +
  geom_raster(aes(fill = estimYs), interpolate = TRUE) +
  # geom_contour(aes(z = estimYs),
  #              breaks = levels_osqc, color = "white", alpha = 0.55) +
  facet_wrap(~ sens, ncol = 1) +
  scale_fill_viridis_c(option = "D", name = "OSQC") +
  labs(
    x = expression(nu),
    y = expression(rho),
    title = "EDPMM OSQC surface under sensitivity",
    subtitle = "Posterior mean"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCsurface.pdf", width = 8, height = 6)
print(EDPqrlPOST_OSQCsurface)
# dev.off()

# Uncertainty surface (CI length) with viridis
# levels_ci <- pretty(df0$CIlength95, n = 8)
EDPqrlPOST_CILsurface <- 
  ggplot(df0, aes(x = nu, y = rho)) +
  geom_raster(aes(fill = CIlength95), interpolate = TRUE) +
  # geom_contour(aes(z = CIlength95),
  #              breaks = levels_ci, color = "white", alpha = 0.55) +
  scale_fill_viridis_c(option = "D", name = "CI length") +
  labs(
    x = expression(nu),
    y = expression(rho),
    title = "EDPMM OSQC uncertainty",
    subtitle = "95% credible interval length"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
# quartz(type = "pdf", file = "EDPqrlPOST_CILsurface.pdf", width = 8, height = 6)
print(EDPqrlPOST_CILsurface)
# dev.off()

# Your facet-by-nu line plot, but make it readable
df_plot <- df_ed %>%
  filter(nu %in% nu_keep) %>%
  mutate(sens = factor(sens))
EDPqrlPOST_OSQCrho_nu_EDPMM <- 
  ggplot(df_plot, aes(x = rho, y = estimYs, color = sens, fill = sens)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(aes(ymin = quantile025, ymax = quantile975),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ nu, ncol = 4) +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)) +
  scale_color_viridis_d(option = "D", end = 0.9) +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(
    x = expression(rho),
    y = "OSQC",
    color = "Sensitivity",
    fill  = "Sensitivity",
    title = "OSQC vs rho at selected landmark times"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCrho_nu_EDPMM.pdf", width = 8, height = 6)
print(EDPqrlPOST_OSQCrho_nu_EDPMM)
# dev.off()

# EDPMM vs DPMM on the same plot (psi0=psi1=0), viridis colors
rho_breaks <- c(0, 0.05, 0.10, 0.15, 0.20)
df_ed0 <- as.data.frame(EDPMM_results_POST$E_rho_nu_diff_result) %>%
  filter(eta_z0 == 0, eta_z1 == 0) %>%
  mutate(model = "EDPMM")
df_dp0 <- as.data.frame(DPMM_results_POST$E_rho_nu_diff_result) %>%
  filter(eta_z0 == 0, eta_z1 == 0) %>%
  mutate(model = "DPMM")
df_both <- bind_rows(df_ed0, df_dp0) %>%
  filter(nu %in% nu_keep) %>%
  mutate(
    model = factor(model, levels = c("EDPMM", "DPMM")),
    rho   = round(as.numeric(rho), 3)   # robust
  )
EDPqrlPOST_OSQCrho_nu_EDPMMvsDPMM <-
  ggplot(df_both, aes(x = rho, y = estimYs, color = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(aes(ymin = quantile025, ymax = quantile975),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ nu, ncol = 4) +
  scale_x_continuous(
    breaks = rho_breaks,
    limits = range(rho_breaks),
    labels = c("0", "0.05", "0.10", "0.15"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  labs(
    x = expression(rho),
    y = "OSQC",
    color = "Model",
    fill = "Model",
    title = "OSQC: EDPMM vs DPMM"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCrho_nu_EDPMMvsDPMM.pdf", width = 8, height = 6)
print(EDPqrlPOST_OSQCrho_nu_EDPMMvsDPMM)
# dev.off()

# Survival curve + KM: make it cleaner and more informative
make_S_df <- function(post_ftn_res, model_name){
  df0 <- as.data.frame(post_ftn_res$S_z0_result)
  df1 <- as.data.frame(post_ftn_res$S_z1_result)
  colnames(df0) <- c("t","estimYs","sd","quantile025","quantile975","CIlength95")
  colnames(df1) <- c("t","estimYs","sd","quantile025","quantile975","CIlength95")
  
  bind_rows(
    df0 %>% mutate(Z = 0),
    df1 %>% mutate(Z = 1)
  ) %>%
    mutate(
      model = model_name,
      group = factor(ifelse(Z == 0, "Non-Elevated", "Elevated"),
                     levels = c("Non-Elevated","Elevated"))
    ) %>%
    arrange(model, group, t)  # prevents weird connecting lines
}

df_S_edp <- make_S_df(SURVIVAL_EDPMM_POST_ftn_results, "EDPMM")
# KM
km_fit <- survfit(Surv(Tobs, Dobs) ~ Zobs, data = dfADNI)
km_sum <- summary(km_fit)
df_KM <- data.frame(
  t      = km_sum$time,
  surv   = km_sum$surv,
  strata = km_sum$strata
) %>%
  mutate(
    group = factor(ifelse(grepl("=0", strata), "Non-Elevated", "Elevated"),
                   levels = c("Non-Elevated","Elevated"))
  ) %>%
  arrange(group, t)
EDPqrlPOST_OSQCsurvival <-
  ggplot() +
  # EDPMM credible band
  geom_ribbon(
    data = df_S_edp,
    aes(x = t, ymin = quantile025, ymax = quantile975, fill = group),
    alpha = 0.18, colour = NA
  ) +
  # EDPMM posterior mean
  geom_line(
    data = df_S_edp,
    aes(x = t, y = estimYs, color = group),
    linewidth = 1.15
  ) +
  # KM curve (dashed)
  geom_step(
    data = df_KM,
    aes(x = t, y = surv, group = group, color = group),
    linetype = "dashed", linewidth = 0.8, alpha = 0.8
  ) +
  labs(
    x = "Time (years)",
    y = "Probability of remaining dementia-free",
    title = "Adjusted survival curves (EDPMM) with Kaplan–Meier overlay",
    subtitle = "Solid: EDPMM posterior mean with 95% CrI; dashed: unadjusted KM"
  ) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCsurvival.pdf", width = 8, height = 6)
print(EDPqrlPOST_OSQCsurvival)
# dev.off()

rho_keep <- c(0.05, 0.10, 0.15)
eps <- 1e-8
df <- as.data.frame(EDPMM_results_POST$E_rho_nu_diff_result) %>%
  # keep only rho close to 0.05/0.10/0.15
  mutate(rho_keep_val = sapply(rho, function(r) {
    i <- which.min(abs(r - rho_keep))
    if (abs(r - rho_keep[i]) < eps) rho_keep[i] else NA_real_
  })) %>%
  filter(!is.na(rho_keep_val)) %>%
  filter(
    (eta_z0 ==  0.0 & eta_z1 ==  0.0) |
      (eta_z0 == -0.1 & eta_z1 ==  0.1) |
      (eta_z0 ==  0.1 & eta_z1 == -0.1)
  ) %>%
  mutate(
    sens = case_when(
      eta_z0 ==  0.0 & eta_z1 ==  0.0 ~ "psi0=0, psi1=0",
      eta_z0 == -0.1 & eta_z1 ==  0.1 ~ "psi0=-0.1, psi1=0.1",
      eta_z0 ==  0.1 & eta_z1 == -0.1 ~ "psi0=0.1, psi1=-0.1"
    ),
    sens = factor(
      sens,
      levels = c("psi0=0, psi1=0", "psi0=-0.1, psi1=0.1", "psi0=0.1, psi1=-0.1")
    ),
    rho = factor(rho_keep_val, levels = rho_keep, labels = paste0("rho = ", rho_keep))
  )
df_zero <- df %>%
  group_by(rho) %>%
  summarise(nu = min(nu, na.rm = TRUE), y0 = 0, .groups = "drop")
EDPqrlPOST_OSQCnu_rho <- 
  ggplot(df, aes(x = nu, y = estimYs, color = sens, fill = sens, group = sens)) +
  geom_blank(data = df_zero, aes(x = nu, y = y0), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(aes(ymin = quantile025, ymax = quantile975),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ rho, ncol = 1, scales = "free_y") +
  scale_color_viridis_d(option = "D", end = 0.9, name = "Sensitivity") +
  scale_fill_viridis_d(option = "D", end = 0.9, name = "Sensitivity") +
  labs(
    x = expression(nu),
    y = "OSQC",
    title = "EDPMM: OSQC as landmark time increases",
    subtitle = "Posterior mean and 95% credible interval"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCnu_rho.pdf", width = 8, height = 6)
print(EDPqrlPOST_OSQCnu_rho)
# dev.off()

# # DPMM: OSQC vs nu at rho = c(0.05, 0.10, 0.15) with sensitivity
# rho_keep <- c(0.05, 0.10, 0.15)
# df <- as.data.frame(DPMM_results_POST$E_rho_nu_diff_result) %>%
#   mutate(rho_num = round(as.numeric(rho), 2)) %>%   # <- key
#   filter(rho_num %in% rho_keep) %>%
#   filter(
#     (eta_z0 ==  0.0 & eta_z1 ==  0.0) |
#       (eta_z0 == -0.1 & eta_z1 ==  0.1) |
#       (eta_z0 ==  0.1 & eta_z1 == -0.1)
#   ) %>%
#   mutate(
#     sens = case_when(
#       eta_z0 ==  0.0 & eta_z1 ==  0.0 ~ "psi0=0, psi1=0",
#       eta_z0 == -0.1 & eta_z1 ==  0.1 ~ "psi0=-0.1, psi1=0.1",
#       eta_z0 ==  0.1 & eta_z1 == -0.1 ~ "psi0=0.1, psi1=-0.1"
#     ),
#     sens = factor(
#       sens,
#       levels = c("psi0=0, psi1=0", "psi0=-0.1, psi1=0.1", "psi0=0.1, psi1=-0.1")
#     ),
#     rho = factor(rho_num, levels = rho_keep, labels = paste0("rho = ", rho_keep))
#   )
# EDPqrlPOST_OSQCnu_rho_DPMM <- 
#   ggplot(df, aes(x = nu, y = estimYs, color = sens, fill = sens, group = sens)) +
#   geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
#   geom_ribbon(aes(ymin = quantile025, ymax = quantile975),
#               alpha = 0.12, colour = NA) +
#   geom_line(linewidth = 0.9) +
#   facet_wrap(~ rho, ncol = 1, scales = "free_y", drop = FALSE) +  # <- keep all panels
#   scale_color_viridis_d(option = "D", end = 0.9, name = "Sensitivity") +
#   scale_fill_viridis_d(option = "D", end = 0.9, name = "Sensitivity") +
#   labs(
#     x = expression(nu),
#     y = "OSQC",
#     title = "DPMM: OSQC as landmark time increases",
#     subtitle = "Posterior mean and 95% credible interval"
#   ) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "bottom")
# # quartz(type = "pdf", file = "EDPqrlPOST_OSQCnu_rho_DPMM.pdf", width = 8, height = 6)
# # print(EDPqrlPOST_OSQCnu_rho_DPMM)
# # dev.off()

# EDPMM vs DPMM OSQC vs nu at rho = c(0.05, 0.10, 0.15) without sensitivity
rho_keep <- c(0.05, 0.10, 0.15)
round2 <- function(x) round(x, 2)
df_ed <- as.data.frame(EDPMM_results_POST$E_rho_nu_diff_result) %>%
  mutate(rho2 = round2(rho)) %>%
  filter(eta_z0 == 0, eta_z1 == 0, rho2 %in% rho_keep) %>%
  mutate(model = "EDPMM")
df_dp <- as.data.frame(DPMM_results_POST$E_rho_nu_diff_result) %>%
  mutate(rho2 = round2(rho)) %>%
  filter(eta_z0 == 0, eta_z1 == 0, rho2 %in% rho_keep) %>%
  mutate(model = "DPMM")
df_both <- bind_rows(df_ed, df_dp) %>%
  mutate(
    model = factor(model, levels = c("EDPMM", "DPMM")),
    rho_f = factor(
      rho2,
      levels = rho_keep,
      labels = sprintf("rho = %.2f", rho_keep)
    )
  )
# ---- smooth within (model, rho_f) ----
df_smooth <- df_both %>%
  arrange(model, rho_f, nu) %>%
  group_by(model, rho_f) %>%
  group_modify(~{
    d <- .x
    spar_val <- 0.7
    s_mu <- smooth.spline(d$nu, d$estimYs,     spar = spar_val)
    s_lo <- smooth.spline(d$nu, d$quantile025, spar = spar_val)
    s_hi <- smooth.spline(d$nu, d$quantile975, spar = spar_val)
    
    d$estimYs_s     <- predict(s_mu, d$nu)$y
    d$quantile025_s <- predict(s_lo, d$nu)$y
    d$quantile975_s <- predict(s_hi, d$nu)$y
    d
  }) %>%
  ungroup()
# include y=0 in each facet but keep free_y compact
df_zero <- df_smooth %>%
  group_by(rho_f) %>%
  summarise(nu = min(nu, na.rm = TRUE), y0 = 0, .groups = "drop")
EDPqrlPOST_OSQCnu_rho_EDPMMvsDPMM <- 
  ggplot(df_smooth, aes(x = nu, y = estimYs_s, color = model, group = model)) +
  geom_blank(data = df_zero, aes(x = nu, y = y0), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(
    aes(ymin = quantile025_s, ymax = quantile975_s, fill = model),
    alpha = 0.12, colour = NA
  ) +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ rho_f, ncol = 1, scales = "free_y") +
  scale_color_viridis_d(option = "D", end = 0.85, name = "Model") +
  scale_fill_viridis_d(option = "D", end = 0.85, name = "Model") +
  labs(
    x = expression(nu),
    y = "OSQC",
    title = "OSQC vs landmark time",
    subtitle = expression(psi[0] == 0 ~ "," ~ psi[1] == 0)
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCnu_rho_EDPMMvsDPMM.pdf", width = 8, height = 6)
print(EDPqrlPOST_OSQCnu_rho_EDPMMvsDPMM)
# dev.off()

# ------------------------------------------------------------------------
# Conditional Effects - Figures ------------------------------------------
# ------------------------------------------------------------------------

# ------------------------------------------------------------
# EDPMM conditional results
# condX = 0 -> CN
# condX = 1 -> MCI
# ------------------------------------------------------------
make_cond_plot_df <- function(res_obj, subgroup_label) {
  df <- as.data.frame(res_obj$E_rho_nu_diff_result)
  # rename columns to match your plotting code
  colnames(df)[1:9] <- c(
    "eta_z0", "eta_z1", "nu", "rho",
    "estimYs", "sd", "quantile025", "quantile975", "CIlength95"
  )
  df %>%
    mutate(subgroup = subgroup_label)
}

df_edpm_cond <- bind_rows(
  make_cond_plot_df(SURVIVAL_EDPMM_POST_condX0, "CN"),
  make_cond_plot_df(SURVIVAL_EDPMM_POST_condX1, "MCI")
) %>%
  filter(
    (eta_z0 ==  0.0 & eta_z1 ==  0.0) |
      (eta_z0 == -0.1 & eta_z1 ==  0.1) |
      (eta_z0 ==  0.1 & eta_z1 == -0.1)
  ) %>%
  mutate(
    sens = case_when(
      eta_z0 ==  0.0 & eta_z1 ==  0.0 ~ "psi0 = 0, psi1 = 0",
      eta_z0 == -0.1 & eta_z1 ==  0.1 ~ "psi0 = -0.1, psi1 = 0.1",
      eta_z0 ==  0.1 & eta_z1 == -0.1 ~ "psi0 = 0.1, psi1 = -0.1"
    ),
    subgroup = factor(subgroup, levels = c("CN", "MCI"))
  )

ggplot(df_edpm_cond,
       aes(x = rho, y = estimYs,
           group = sens, color = sens, fill = sens)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(aes(ymin = quantile025, ymax = quantile975),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.3) +
  facet_grid(subgroup ~ nu, scales = "free_y") +
  labs(
    x = expression(rho),
    y = "Conditional OSQC estimate",
    color = "Sensitivity setting",
    fill = "Sensitivity setting",
    title = "EDPMM: Posterior mean and 95% credible interval for conditional OSQC",
    subtitle = "Rows: baseline diagnosis subgroup; columns: landmark time"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# ------------------------------------------------------------
# DPMM conditional results
# condX = 0 -> CN
# condX = 1 -> MCI
# ------------------------------------------------------------
df_dpmm_cond <- bind_rows(
  make_cond_plot_df(SURVIVAL_DPMM_POST_condX0, "CN"),
  make_cond_plot_df(SURVIVAL_DPMM_POST_condX1, "MCI")
) %>%
  filter(
    (eta_z0 ==  0.0 & eta_z1 ==  0.0) |
      (eta_z0 == -0.1 & eta_z1 ==  0.1) |
      (eta_z0 ==  0.1 & eta_z1 == -0.1)
  ) %>%
  mutate(
    sens = case_when(
      eta_z0 ==  0.0 & eta_z1 ==  0.0 ~ "psi0 = 0, psi1 = 0",
      eta_z0 == -0.1 & eta_z1 ==  0.1 ~ "psi0 = -0.1, psi1 = 0.1",
      eta_z0 ==  0.1 & eta_z1 == -0.1 ~ "psi0 = 0.1, psi1 = -0.1"
    ),
    subgroup = factor(subgroup, levels = c("CN", "MCI"))
  )

ggplot(df_dpmm_cond,
       aes(x = rho, y = estimYs,
           group = sens, color = sens, fill = sens)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(aes(ymin = quantile025, ymax = quantile975),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.3) +
  facet_grid(subgroup ~ nu, scales = "free_y") +
  labs(
    x = expression(rho),
    y = "Conditional OSQC estimate",
    color = "Sensitivity setting",
    fill = "Sensitivity setting",
    title = "DPMM: Posterior mean and 95% credible interval for conditional OSQC",
    subtitle = "Rows: baseline diagnosis subgroup; columns: landmark time"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# EDPMM survival curve + KM conditional CN or MCI
# Conditional EDPMM survival function results
SURVIVAL_EDPMM_POST_ftn_CN <- 
  SURVIVAL_EDPMM_POST_ftn(
    EDPMM_results_MCMC, M, t_grid,
    ftn = c("S", "f"),
    condX = 0, p_condX = 1
  )

SURVIVAL_EDPMM_POST_ftn_MCI <- 
  SURVIVAL_EDPMM_POST_ftn(
    EDPMM_results_MCMC, M, t_grid,
    ftn = c("S", "f"),
    condX = 1, p_condX = 1
  )

# -----------------------------
# EDPMM conditional survival dfs
# -----------------------------
make_surv_df <- function(res, subgroup, group) {
  df <- as.data.frame(res)
  colnames(df) <- c("t", "estimYs", "sd", "quantile025", "quantile975", "CIlength95")
  df %>%
    mutate(
      subgroup = subgroup,
      group = group
    )
}

df_edpm_S <- bind_rows(
  make_surv_df(SURVIVAL_EDPMM_POST_ftn_condX0_results$S_z0_result, "CN",  "Non-Elevated"),
  make_surv_df(SURVIVAL_EDPMM_POST_ftn_condX0_results$S_z1_result, "CN",  "Elevated"),
  make_surv_df(SURVIVAL_EDPMM_POST_ftn_condX1_results$S_z0_result, "MCI", "Non-Elevated"),
  make_surv_df(SURVIVAL_EDPMM_POST_ftn_condX1_results$S_z1_result, "MCI", "Elevated")
)

df_edpm_S <- df_edpm_S %>%
  mutate(curve = paste(subgroup, group, sep = ", "))

curve_cols <- c(
  "CN, Non-Elevated"  = "red",
  "CN, Elevated"      = "blue",
  "MCI, Non-Elevated" = "red",
  "MCI, Elevated"     = "blue"
)

curve_lts <- c(
  "CN, Non-Elevated"  = "solid",
  "CN, Elevated"      = "solid",
  "MCI, Non-Elevated" = "dashed",
  "MCI, Elevated"     = "dashed"
)

ggplot(df_edpm_S,
       aes(x = t, y = estimYs,
           color = curve, fill = curve, linetype = curve,
           group = curve)) +
  geom_ribbon(
    aes(ymin = quantile025, ymax = quantile975),
    alpha = 0.10, colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = curve_cols) +
  scale_fill_manual(values = curve_cols) +
  scale_linetype_manual(values = curve_lts) +
  labs(
    x = "Time (years)",
    y = "Probability of remaining dementia-free",
    color = "Curve",
    fill = "Curve",
    linetype = "Curve",
    title = "EDPMM conditional survival curves",
    subtitle = "One-panel display"
  ) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# -----------------------------
# DPMM conditional survival dfs
# -----------------------------
df_dpmm_S <- bind_rows(
  make_surv_df(SURVIVAL_DPMM_POST_ftn_condX0_results$S_z0_result, "CN",  "Non-Elevated"),
  make_surv_df(SURVIVAL_DPMM_POST_ftn_condX0_results$S_z1_result, "CN",  "Elevated"),
  make_surv_df(SURVIVAL_DPMM_POST_ftn_condX1_results$S_z0_result, "MCI", "Non-Elevated"),
  make_surv_df(SURVIVAL_DPMM_POST_ftn_condX1_results$S_z1_result, "MCI", "Elevated")
)

df_dpmm_S <- df_dpmm_S %>%
  mutate(curve = paste(subgroup, group, sep = ", "))

curve_cols <- c(
  "CN, Non-Elevated"  = "red",
  "CN, Elevated"      = "blue",
  "MCI, Non-Elevated" = "red",
  "MCI, Elevated"     = "blue"
)

curve_lts <- c(
  "CN, Non-Elevated"  = "solid",
  "CN, Elevated"      = "solid",
  "MCI, Non-Elevated" = "dashed",
  "MCI, Elevated"     = "dashed"
)

ggplot(df_dpmm_S,
       aes(x = t, y = estimYs,
           color = curve, fill = curve, linetype = curve,
           group = curve)) +
  geom_ribbon(
    aes(ymin = quantile025, ymax = quantile975),
    alpha = 0.10, colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = curve_cols) +
  scale_fill_manual(values = curve_cols) +
  scale_linetype_manual(values = curve_lts) +
  labs(
    x = "Time (years)",
    y = "Probability of remaining dementia-free",
    color = "Curve",
    fill = "Curve",
    linetype = "Curve",
    title = "DPMM conditional survival curves",
    subtitle = "One-panel display"
  ) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# ------------------------------------------------------------
# Helper: build conditional OSQC plotting data
# condX = 0 -> CN
# condX = 1 -> MCI
# ------------------------------------------------------------
make_cond_plot_df <- function(res_obj, subgroup_label) {
  df <- as.data.frame(res_obj$E_rho_nu_diff_result)
  colnames(df)[1:9] <- c(
    "eta_z0", "eta_z1", "nu", "rho",
    "estimYs", "sd", "quantile025", "quantile975", "CIlength95"
  )
  df %>%
    mutate(subgroup = subgroup_label)
}

# ------------------------------------------------------------
# Helper: plot conditional OSQC as nu increases
# ------------------------------------------------------------
plot_cond_osqc_nu <- function(df_plot, plot_title) {
  ggplot(
    df_plot,
    aes(x = nu, y = estimYs,
        group = sens, color = sens, fill = sens)
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
    geom_ribbon(
      aes(ymin = quantile025, ymax = quantile975),
      alpha = 0.12, colour = NA
    ) +
    geom_line(linewidth = 0.8) +
    facet_grid(subgroup ~ rho, scales = "free_y") +
    labs(
      x = expression(nu),
      y = "Conditional OSQC estimate",
      color = "Sensitivity setting",
      fill = "Sensitivity setting",
      title = plot_title,
      subtitle = "Rows: baseline diagnosis subgroup; columns: quantile level"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
}

# ------------------------------------------------------------
# EDPMM conditional results
# condX = 0 -> CN
# condX = 1 -> MCI
# ------------------------------------------------------------
df_edpm_cond_nu <- bind_rows(
  make_cond_plot_df(SURVIVAL_EDPMM_POST_condX0_nu, "CN"),
  make_cond_plot_df(SURVIVAL_EDPMM_POST_condX1_nu, "MCI")
) %>%
  filter(
    (eta_z0 ==  0.0 & eta_z1 ==  0.0) |
      (eta_z0 == -0.1 & eta_z1 ==  0.1) |
      (eta_z0 ==  0.1 & eta_z1 == -0.1)
  ) %>%
  mutate(
    sens = case_when(
      eta_z0 ==  0.0 & eta_z1 ==  0.0 ~ "psi0 = 0, psi1 = 0",
      eta_z0 == -0.1 & eta_z1 ==  0.1 ~ "psi0 = -0.1, psi1 = 0.1",
      eta_z0 ==  0.1 & eta_z1 == -0.1 ~ "psi0 = 0.1, psi1 = -0.1"
    ),
    subgroup = factor(subgroup, levels = c("CN", "MCI")),
    rho = factor(rho, levels = sort(unique(rho)),
                 labels = paste0("rho = ", sort(unique(rho))))
  )

plot_cond_osqc_nu(
  df_edpm_cond_nu,
  "EDPMM: Conditional OSQC as landmark time increases"
)

# ------------------------------------------------------------
# DPMM conditional results
# condX = 0 -> CN
# condX = 1 -> MCI
# ------------------------------------------------------------
df_dpmm_cond_nu <- bind_rows(
  make_cond_plot_df(SURVIVAL_DPMM_POST_condX0_nu, "CN"),
  make_cond_plot_df(SURVIVAL_DPMM_POST_condX1_nu, "MCI")
) %>%
  filter(
    (eta_z0 ==  0.0 & eta_z1 ==  0.0) |
      (eta_z0 == -0.1 & eta_z1 ==  0.1) |
      (eta_z0 ==  0.1 & eta_z1 == -0.1)
  ) %>%
  mutate(
    sens = case_when(
      eta_z0 ==  0.0 & eta_z1 ==  0.0 ~ "psi0 = 0, psi1 = 0",
      eta_z0 == -0.1 & eta_z1 ==  0.1 ~ "psi0 = -0.1, psi1 = 0.1",
      eta_z0 ==  0.1 & eta_z1 == -0.1 ~ "psi0 = 0.1, psi1 = -0.1"
    ),
    subgroup = factor(subgroup, levels = c("CN", "MCI")),
    rho = factor(rho, levels = sort(unique(rho)),
                 labels = paste0("rho = ", sort(unique(rho))))
  )

plot_cond_osqc_nu(
  df_dpmm_cond_nu,
  "DPMM: Conditional OSQC as landmark time increases"
)

# save.image(paste0("~/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNImcmc+post+plot.RData"))