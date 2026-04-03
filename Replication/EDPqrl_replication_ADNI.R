# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
setwd("")

# Load R code
source("EDPqrl_r.R")

# Load cpp code
sourceCpp("EDPqrl_cpp.cpp")

# seed
set.seed(1)

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
library(Rcpp)
library(RcppArmadillo)

library(survival)

library(viridis)
library(akima)
library(scales)
library(survminer)
library(grid)
library(ggplot2)
library(patchwork)

library(dplyr)
library(tidyr)
library(tibble)
library(tidyverse)
library(knitr)
library(kableExtra)
library(htmltools)

# setwd("/Users/woojung/Documents/Rproject/EDPqrl/RealDataAnalysis/DATA/ADNI")
# install.packages("ADNIMERGE2.tar.gz", repos = NULL, type = "source")
library(ADNIMERGE2)

# Due to serious data quality concerns, all data with PTIDs containing the format 
# 381_S_10### has been removed from the ADNI data repository. ADNI leadership 
# strongly advises that data with these PTIDs not be used in any data analysis, 
# submission of abstracts, or publications. Furthermore, ADNI leadership requests 
# that any abstracts or papers which have been submitted using this data be 
# appropriately revised due to concerns with data quality. 

bad_usubjid <- unique(ADSL$USUBJID[grepl("^ADNI-381-[0-9]{5}$", ADSL$USUBJID)])
length(bad_usubjid)

ADSL <- ADSL[!(ADSL$USUBJID %in% bad_usubjid), ]
ADAE <- ADAE[!(ADAE$USUBJID %in% bad_usubjid), ]
ADQS <- ADQS[!(ADQS$USUBJID %in% bad_usubjid), ]
ADRS <- ADRS[!(ADRS$USUBJID %in% bad_usubjid), ]

sum(ADSL$USUBJID %in% bad_usubjid)
sum(ADAE$USUBJID %in% bad_usubjid)
sum(ADQS$USUBJID %in% bad_usubjid)
sum(ADRS$USUBJID %in% bad_usubjid)

# ------------------------------------------------------------------------
# Load Data --------------------------------------------------------------
# ------------------------------------------------------------------------
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

##========================================================
## 0. IDs for cohort flow
##========================================================

adsl_ids <- ADSL %>%
  dplyr::distinct(USUBJID)

dx_ids <- ADRS %>%
  dplyr::filter(PARAMCD == "DX") %>%
  dplyr::distinct(USUBJID)

lost_before_merge <- adsl_ids %>%
  dplyr::anti_join(dx_ids, by = "USUBJID")

##========================================================
## 1. Define exclusion flags safely
##========================================================

merge_data_flagged <- merge_data %>%
  dplyr::mutate(
    excl_missing_amyloid = is.na(Amyloid_Positive),
    excl_baseline_dem    = dplyr::coalesce(Baseline_DX == "DEM", FALSE),
    excl_bad_time        = is.na(T_years) | dplyr::coalesce(T_years <= 0, FALSE),
    exclude_final        = excl_missing_amyloid | excl_baseline_dem | excl_bad_time
  )

##========================================================
## 2. Exclusion summary within merge_data
##========================================================

exclusion_summary <- merge_data_flagged %>%
  dplyr::summarise(
    Missing_Amyloid      = sum(is.na(Amyloid_Positive)),
    Missing_APOE         = sum(is.na(APOE)),
    Missing_Baseline_DX  = sum(is.na(Baseline_DX)),
    Baseline_DEM         = sum(dplyr::coalesce(Baseline_DX == "DEM", FALSE)),
    Time_is_Zero_or_Less = sum(dplyr::coalesce(T_years <= 0, FALSE)),
    Time_is_NA           = sum(is.na(T_years)),
    
    Missing_Amyloid_and_Base_DEM =
      sum(is.na(Amyloid_Positive) & dplyr::coalesce(Baseline_DX == "DEM", FALSE)),
    
    Missing_APOE_and_Base_DEM =
      sum(is.na(APOE) & dplyr::coalesce(Baseline_DX == "DEM", FALSE)),
    
    Missing_Both_Amyloid_APOE =
      sum(is.na(Amyloid_Positive) & is.na(APOE)),
    
    Base_DEM_with_Valid_Time =
      sum(
        dplyr::coalesce(Baseline_DX == "DEM", FALSE) &
          !is.na(T_years) &
          dplyr::coalesce(T_years > 0, FALSE)
      ),
    
    Total_Excluded_Subjects = sum(exclude_final)
  )

##========================================================
## 3. Final analytic data
##========================================================

final_data <- merge_data_flagged %>%
  dplyr::filter(!exclude_final) %>%
  dplyr::mutate(
    SEX = ifelse(SEX == "Female", 1, 0)
  ) %>%
  dplyr::select(
    T_years, Status,
    Amyloid_Positive,
    APOE4_Carrier, Baseline_DX, SEX, AGE, EDUC
  )

##========================================================
## 4. Dynamic cohort flow summary
##========================================================

N_ADSL  <- nrow(adsl_ids)
N_DX    <- nrow(dx_ids)
N_merge <- nrow(merge_data %>% dplyr::distinct(USUBJID))
N_final <- nrow(final_data)

Excluded_before_merge <- nrow(lost_before_merge)
Excluded_from_merge_to_final <- N_merge - N_final
Total_Excluded_from_ADSL <- N_ADSL - N_final

flow_summary <- tibble::tibble(
  Step = c(
    "ADSL",
    "Excluded before merge: no ADRS DX record",
    "merge_data",
    "Excluded from merge_data by analytic criteria",
    "final_data"
  ),
  N = c(
    N_ADSL,
    Excluded_before_merge,
    N_merge,
    Excluded_from_merge_to_final,
    N_final
  )
)

print(flow_summary)

## Optional checks
stopifnot(N_ADSL - Excluded_before_merge == N_merge)
stopifnot(N_merge - Excluded_from_merge_to_final == N_final)
stopifnot(N_ADSL - Total_Excluded_from_ADSL == N_final)

##========================================================
## 5. Baseline summary for final analytic cohort
##========================================================

n_event <- sum(final_data$Status == 1, na.rm = TRUE)
pct_event <- 100 * n_event / N_final

n_amyloid <- sum(final_data$Amyloid_Positive == 1, na.rm = TRUE)
pct_amyloid <- 100 * n_amyloid / N_final

n_carrier <- sum(final_data$APOE4_Carrier == 1, na.rm = TRUE)
pct_carrier <- 100 * n_carrier / N_final

age_mean <- mean(final_data$AGE, na.rm = TRUE)
age_sd   <- sd(final_data$AGE, na.rm = TRUE)

pct_female <- 100 * mean(final_data$SEX, na.rm = TRUE)

educ_mean <- mean(final_data$EDUC, na.rm = TRUE)
educ_sd   <- sd(final_data$EDUC, na.rm = TRUE)

miss_tbl <- final_data %>%
  dplyr::summarise(
    Missing_Baseline_DX = sum(is.na(Baseline_DX)),
    Missing_Amyloid     = sum(is.na(Amyloid_Positive)),
    Missing_APOE        = sum(is.na(APOE4_Carrier)),
    Missing_AGE         = sum(is.na(AGE)),
    Missing_SEX         = sum(is.na(SEX)),
    Missing_EDUC        = sum(is.na(EDUC))
  )

##========================================================
## 6. Final output
##========================================================

results_summary <- list(
  flow_summary = flow_summary,
  exclusion_summary = exclusion_summary,
  final_cohort_summary = list(
    N_final = N_final,
    events = c(n_event = n_event, pct_event = pct_event),
    amyloid_positive = c(n_amyloid = n_amyloid, pct_amyloid = pct_amyloid),
    apoe4_carrier = c(n_carrier = n_carrier, pct_carrier = pct_carrier),
    age = c(mean = age_mean, sd = age_sd),
    female = c(pct_female = pct_female),
    educ = c(mean = educ_mean, sd = educ_sd)
  ),
  missing = miss_tbl
)
print(results_summary)
dim(ADSL)

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

## pull matrix columns back out for summaries
dfADNI2 <- dfADNI %>%
  dplyr::mutate(
    Baseline_DX_MCI = Xobs[, "Baseline_DX"],
    APOE4_Carrier   = Xobs[, "APOE4_Carrier"],
    SEX             = Xobs[, "SEX"],
    AGE             = Xobs[, "AGE"],
    EDUC            = Xobs[, "EDUC"]
  )

##---------------------------------
## overall censoring
##---------------------------------
overall_censoring <- dfADNI2 %>%
  dplyr::summarise(
    N = dplyr::n(),
    n_censored = sum(Dobs == 0, na.rm = TRUE),
    censoring_rate = mean(Dobs == 0, na.rm = TRUE),
    n_event = sum(Dobs == 1, na.rm = TRUE),
    event_rate = mean(Dobs == 1, na.rm = TRUE)
  )

overall_censoring

##---------------------------------
## censoring by amyloid
##---------------------------------
censor_by_amyloid <- dfADNI2 %>%
  dplyr::mutate(
    Amyloid_Group = dplyr::case_when(
      is.na(Zobs) ~ "Missing",
      Zobs == 0 ~ "Amyloid Negative",
      Zobs == 1 ~ "Amyloid Positive"
    )
  ) %>%
  dplyr::group_by(Amyloid_Group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    n_censored = sum(Dobs == 0, na.rm = TRUE),
    censoring_rate = mean(Dobs == 0, na.rm = TRUE),
    n_event = sum(Dobs == 1, na.rm = TRUE),
    event_rate = mean(Dobs == 1, na.rm = TRUE),
    .groups = "drop"
  )

censor_by_amyloid

##---------------------------------
## censoring by baseline diagnosis (CN vs MCI only)
##---------------------------------
censor_by_dx <- dfADNI2 %>%
  dplyr::mutate(
    Baseline_DX_Group = dplyr::case_when(
      is.na(Baseline_DX_MCI) ~ "Missing",
      Baseline_DX_MCI == 0 ~ "CN",
      Baseline_DX_MCI == 1 ~ "MCI"
    )
  ) %>%
  dplyr::group_by(Baseline_DX_Group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    n_censored = sum(Dobs == 0, na.rm = TRUE),
    censoring_rate = mean(Dobs == 0, na.rm = TRUE),
    n_event = sum(Dobs == 1, na.rm = TRUE),
    event_rate = mean(Dobs == 1, na.rm = TRUE),
    .groups = "drop"
  )

censor_by_dx

##---------------------------------
## censoring by APOE4 carrier
##---------------------------------
censor_by_apoe4 <- dfADNI2 %>%
  dplyr::mutate(
    APOE4_Group = dplyr::case_when(
      is.na(APOE4_Carrier) ~ "Missing",
      APOE4_Carrier == 0 ~ "Non-carrier",
      APOE4_Carrier == 1 ~ "Carrier"
    )
  ) %>%
  dplyr::group_by(APOE4_Group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    n_censored = sum(Dobs == 0, na.rm = TRUE),
    censoring_rate = mean(Dobs == 0, na.rm = TRUE),
    n_event = sum(Dobs == 1, na.rm = TRUE),
    event_rate = mean(Dobs == 1, na.rm = TRUE),
    .groups = "drop"
  )

censor_by_apoe4

##---------------------------------
## censoring by amyloid x CN/MCI
##---------------------------------
censor_by_amyloid_dx <- dfADNI2 %>%
  dplyr::mutate(
    Amyloid_Group = dplyr::case_when(
      is.na(Zobs) ~ "Missing",
      Zobs == 0 ~ "Amyloid Negative",
      Zobs == 1 ~ "Amyloid Positive"
    ),
    Baseline_DX_Group = dplyr::case_when(
      is.na(Baseline_DX_MCI) ~ "Missing",
      Baseline_DX_MCI == 0 ~ "CN",
      Baseline_DX_MCI == 1 ~ "MCI"
    )
  ) %>%
  dplyr::group_by(Amyloid_Group, Baseline_DX_Group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    n_censored = sum(Dobs == 0, na.rm = TRUE),
    censoring_rate = mean(Dobs == 0, na.rm = TRUE),
    n_event = sum(Dobs == 1, na.rm = TRUE),
    event_rate = mean(Dobs == 1, na.rm = TRUE),
    .groups = "drop"
  )

censor_by_amyloid_dx

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
  title = "Kaplan–Meier Estimates of Remaining Dementia-Free by Baseline Amyloid Status",
  xlab = "Time (Years)",
  ylab = "Probability of Remaining Dementia-Free",
  legend.title = "Amyloid Status",
  legend.labs = c("Non-Elevated (Z=0)", "Elevated (Z=1)")
)
print(km_plot)
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

# 1. Collapse event/censored into amyloid-status groups
df_plot_main <- df_plot %>%
  mutate(
    Group_main = case_when(
      Group == "Overall" ~ "Overall",
      grepl("^Z=0", Group) ~ "Non-Elevated (Z=0)",
      grepl("^Z=1", Group) ~ "Elevated (Z=1)"
    ),
    Group_main = factor(
      Group_main,
      levels = c("Overall", "Non-Elevated (Z=0)", "Elevated (Z=1)")
    )
  )

# 2. Generate the plot on the original scale
EDPqrl_DENSITY_T <- ggplot(
  df_plot_main,
  aes(x = Tobs, color = Group_main, linetype = Group_main)
) +
  geom_density(linewidth = 1.2, key_glyph = "path") +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 0.17)) +
  scale_color_manual(
    values = c(
      "Overall" = "black",
      "Non-Elevated (Z=0)" = "red",
      "Elevated (Z=1)" = "blue"
    ),
    breaks = c("Overall", "Non-Elevated (Z=0)", "Elevated (Z=1)")
  ) +
  scale_linetype_manual(
    values = c(
      "Overall" = "solid",
      "Non-Elevated (Z=0)" = "dashed",
      "Elevated (Z=1)" = "dotted"
    ),
    breaks = c("Overall", "Non-Elevated (Z=0)", "Elevated (Z=1)")
  ) +
  theme_classic() +
  labs(
    title = "Density of Observed Dementia-Free Follow-up Time",
    x = "Time to Dementia (Years)",
    y = "Density",
    color = "Amyloid Status",
    linetype = "Amyloid Status"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.70, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      color = "black",
      linewidth = 0.5
    ),
    legend.key.width = unit(3, "lines"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  ) +
  guides(
    color = guide_legend(ncol = 1),
    linetype = guide_legend(ncol = 1)
  )
print(EDPqrl_DENSITY_T)
# quartz(type = "pdf", file = "EDPqrl_DENSITY_T.pdf", width = 8, height = 6)
# print(EDPqrl_DENSITY_T)
# dev.off()

# ------------------------------------------------------------------------
# Constants --------------------------------------------------------------
# ------------------------------------------------------------------------
gibbs_thin = 2e1
gibbs_iter = 2e4
gibbs_burnin = 2e4

# ------------------------------------------------------------------------
# MCMC sampling ----------------------------------------------------------
# ------------------------------------------------------------------------
# EDPMM
EDPMM_results_MCMC = SURVIVAL_EDPMM_MCMC(dfADNI, gibbs_iter, gibbs_burnin, gibbs_thin)
EDPMM_results_MCMC_p025 = SURVIVAL_EDPMM_MCMC(dfADNI, gibbs_iter, gibbs_burnin, gibbs_thin,
                                              location = 0.25)
EDPMM_results_MCMC_m025 = SURVIVAL_EDPMM_MCMC(dfADNI, gibbs_iter, gibbs_burnin, gibbs_thin,
                                              location = -0.25)

# ------------------------------------------------------------------------
# alpha level = 0.01
# ------------------------------------------------------------------------
level = 0.01

# ------------------------------------------------------------------------
# Average Effects & Sensitivity Analysis ---------------------------------
# ------------------------------------------------------------------------
M = 1e3

# Sensitivity Analysis
eta_z0 = 0
eta_z1 = 0

# censoring = 1 - mean(dfADNI$Dobs) = 85.83397
# summary(dfADNI$Tobs)
# min:  0.0027
# max: 14.8419
# t_grid = seq(1e-3, 20, by = 0.1)  => extrapolation from 14.9 - 20.0
# nu = c(1:200)/10                  => extrapolation from 14.9 - 20.0
# rho = c(1:30)/100                 => extrapolation from 0.15 - 0.20
nu = seq(0, 15, by = 0.5)
rho = seq(0.05, 0.30, by = 0.05)
t_grid = c(1e-3, seq(0.5, 20, by = 0.5))

# ------------------------------------------------------------------------
# Average effects: EDPMM
# ------------------------------------------------------------------------
EDPMM_results_POST_p025 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC_p025, M,
                               nu, rho, eta_z0, eta_z1,
                               level = level)

EDPMM_results_POST_m025 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC_m025, M,
                               nu, rho, eta_z0, eta_z1,
                               level = level)

SURVIVAL_EDPMM_POST_ftn_results_p025 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC_p025, M, t_grid,
                          ftn = c("S", "f"),
                          condX = NULL, p_condX = NULL,
                          esttype = "mean",
                          level = level)

SURVIVAL_EDPMM_POST_ftn_results_m025 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC_m025, M, t_grid,
                          ftn = c("S", "f"),
                          condX = NULL, p_condX = NULL,
                          esttype = "mean",
                          level = level)

# ------------------------------------------------------------------------
# Conditional Effects & Sensitivity Analysis -----------------------------
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Conditional effects: EDPMM, location = +0.25
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_condX1_p025 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC_p025, M,
                               nu, rho, eta_z0, eta_z1,
                               condX = 1, p_condX = 1,
                               level = level)

SURVIVAL_EDPMM_POST_condX0_p025 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC_p025, M,
                               nu, rho, eta_z0, eta_z1,
                               condX = 0, p_condX = 1,
                               level = level)

SURVIVAL_EDPMM_POST_ftn_condX1_results_p025 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC_p025, M, t_grid,
                          ftn = c("S", "f"),
                          condX = 1, p_condX = 1,
                          esttype = "mean",
                          level = level)

SURVIVAL_EDPMM_POST_ftn_condX0_results_p025 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC_p025, M, t_grid,
                          ftn = c("S", "f"),
                          condX = 0, p_condX = 1,
                          esttype = "mean",
                          level = level)

# ------------------------------------------------------------------------
# Conditional effects: EDPMM, location = -0.25
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_condX1_m025 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC_m025, M,
                               nu, rho, eta_z0, eta_z1,
                               condX = 1, p_condX = 1,
                               level = level)

SURVIVAL_EDPMM_POST_condX0_m025 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC_m025, M,
                               nu, rho, eta_z0, eta_z1,
                               condX = 0, p_condX = 1,
                               level = level)

SURVIVAL_EDPMM_POST_ftn_condX1_results_m025 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC_m025, M, t_grid,
                          ftn = c("S", "f"),
                          condX = 1, p_condX = 1,
                          esttype = "mean",
                          level = level)

SURVIVAL_EDPMM_POST_ftn_condX0_results_m025 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC_m025, M, t_grid,
                          ftn = c("S", "f"),
                          condX = 0, p_condX = 1,
                          esttype = "mean",
                          level = level)

# ------------------------------------------------------------------------
# Conditional effects: DPMM, location = +0.25
# ------------------------------------------------------------------------
SURVIVAL_DPMM_POST_condX1_p025 =
  SURVIVAL_DPMM_POST_parallel(DPMM_results_MCMC_p025, M,
                              nu, rho, eta_z0, eta_z1,
                              condX = 1, p_condX = 1,
                              level = level)

SURVIVAL_DPMM_POST_condX0_p025 =
  SURVIVAL_DPMM_POST_parallel(DPMM_results_MCMC_p025, M,
                              nu, rho, eta_z0, eta_z1,
                              condX = 0, p_condX = 1,
                              level = level)

SURVIVAL_DPMM_POST_ftn_condX1_results_p025 =
  SURVIVAL_DPMM_POST_ftn(DPMM_results_MCMC_p025, M, t_grid,
                         ftn = c("S", "f"),
                         condX = 1, p_condX = 1,
                         esttype = "mean",
                         level = level)

SURVIVAL_DPMM_POST_ftn_condX0_results_p025 =
  SURVIVAL_DPMM_POST_ftn(DPMM_results_MCMC_p025, M, t_grid,
                         ftn = c("S", "f"),
                         condX = 0, p_condX = 1,
                         esttype = "mean",
                         level = level)

# ------------------------------------------------------------------------
# alpha level = 0.01
# ------------------------------------------------------------------------
level = 0.01

# ------------------------------------------------------------------------
# Average Effects & Sensitivity Analysis ---------------------------------
# ------------------------------------------------------------------------
M = 1e3

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
nu = seq(0, 15, by = 0.5) # => 30
rho = seq(0.05, 0.30, by = 0.05) # => 6
t_grid = c(1e-3, seq(0.5, 20, by = 0.5)) # => 41

EDPMM_results_POST_001 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC, M,
                               nu, rho, eta_z0, eta_z1,
                               level = level)

SURVIVAL_EDPMM_POST_ftn_results_001 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                          condX = NULL, p_condX = NULL, esttype = "mean",
                          level = level)

# ------------------------------------------------------------------------
# Conditional Effects & Sensitivity Analysis -----------------------------
# ------------------------------------------------------------------------
# Baseline_DX = 1 # Baseline_DX = ifelse(Baseline_DX == "MCI", 1, 0)
SURVIVAL_EDPMM_POST_condX1_001 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC, M,
                               nu, rho, eta_z0, eta_z1,
                               condX = 1, p_condX = 1,
                               level = level)
SURVIVAL_EDPMM_POST_condX0_001 =
  SURVIVAL_EDPMM_POST_parallel(EDPMM_results_MCMC, M,
                               nu, rho, eta_z0, eta_z1,
                               condX = 0, p_condX = 1,
                               level = level)

# summary(dfADNI$Tobs)
# min:  0.0027
# max: 14.8419
SURVIVAL_EDPMM_POST_ftn_condX1_results_001 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                          condX = 1, p_condX = 1, esttype = "mean",
                          level = level)
SURVIVAL_EDPMM_POST_ftn_condX0_results_001 =
  SURVIVAL_EDPMM_POST_ftn(EDPMM_results_MCMC, M, t_grid, ftn = c("S", "f"),
                          condX = 0, p_condX = 1, esttype = "mean",
                          level = level)

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

# ============================================================
# SETTINGS
# ============================================================

rho_keep_tbl <- c(0.10, 0.20, 0.30)
nu_keep_tbl  <- c(0, 5, 10)

subgroup_order <- c("Overall", "MCI", "CN")

preferred_other_pairs <- tibble::tribble(
  ~eta_z0, ~eta_z1,
  0.00,  -0.25,
  -0.25,   0.00,
  -0.25,  -0.25,
  
  0.00,   0.25,
  0.25,   0.00,
  0.25,   0.25,
  
  0.25,  -0.25,  
  -0.25,   0.25
)

# ============================================================
# CLEAN ENVIRONMENT
# ============================================================

rm(list = intersect(
  c(
    "df_osqc_master",
    "df_overall", "df_mci", "df_cn",
    "tab1_df", "tab1_tex",
    "tab2_source", "tab2_completed", "tab2_viewer_df", "tab2_latex_df", "tab2_tex"
  ),
  ls()
))

# ============================================================
# HELPERS
# ============================================================

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x) | is.nan(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_eta <- function(x, digits = 2) {
  out <- formatC(x, format = "f", digits = digits)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  out
}

safe_view <- function(x) {
  if (interactive()) {
    View(x)
  } else {
    print(as_tibble(x), n = Inf)
  }
}

safe_print_all <- function(x) {
  x %>% as_tibble() %>% print(n = Inf)
}

label_pair_viewer <- function(eta0, eta1) {
  if (isTRUE(all.equal(eta0, 0)) && isTRUE(all.equal(eta1, 0))) {
    "OSQC (psi0=0, psi1=0)"
  } else {
    paste0("psi0=", fmt_eta(eta0), ", psi1=", fmt_eta(eta1))
  }
}

label_pair_latex <- function(eta0, eta1) {
  if (isTRUE(all.equal(eta0, 0)) && isTRUE(all.equal(eta1, 0))) {
    "OSQC ($\\psi_{0}=0,\\ \\psi_{1}=0$)"
  } else {
    paste0(
      "$\\psi_{0}=", fmt_eta(eta0),
      ",\\ \\psi_{1}=", fmt_eta(eta1), "$"
    )
  }
}

build_setting_levels <- function(df_pairs, preferred_pairs) {
  pairs_all <- df_pairs %>%
    distinct(eta_z0, eta_z1)
  
  pairs_primary <- pairs_all %>%
    filter(eta_z0 == 0, eta_z1 == 0)
  
  pairs_other <- pairs_all %>%
    filter(!(eta_z0 == 0 & eta_z1 == 0))
  
  pairs_preferred <- preferred_pairs %>%
    semi_join(pairs_other, by = c("eta_z0", "eta_z1"))
  
  pairs_remaining <- pairs_other %>%
    anti_join(pairs_preferred, by = c("eta_z0", "eta_z1")) %>%
    arrange(eta_z0, eta_z1)
  
  ordered_pairs <- bind_rows(pairs_primary, pairs_preferred, pairs_remaining)
  
  list(
    viewer = mapply(label_pair_viewer, ordered_pairs$eta_z0, ordered_pairs$eta_z1),
    latex  = mapply(label_pair_latex,  ordered_pairs$eta_z0, ordered_pairs$eta_z1)
  )
}

format_osqc_df <- function(res_obj) {
  df <- as.data.frame(res_obj$E_rho_nu_diff_result)
  colnames(df)[1:9] <- c(
    "eta_z0", "eta_z1", "nu", "rho",
    "estimYs", "sd", "quantile025", "quantile975", "CIlength95"
  )
  
  df %>%
    mutate(
      eta_z0 = as.numeric(eta_z0),
      eta_z1 = as.numeric(eta_z1),
      nu     = as.numeric(nu),
      rho    = round(as.numeric(rho), 2)
    )
}

make_panel_df <- function(res_obj, panel_label) {
  format_osqc_df(res_obj) %>%
    mutate(subgroup = panel_label)
}

# ============================================================
# BUILD MASTER FROM CURRENT OBJECTS ONLY
# ============================================================

df_overall <- make_panel_df(EDPMM_results_POST_001, "Overall")
df_mci     <- make_panel_df(SURVIVAL_EDPMM_POST_condX1_001, "MCI")
df_cn      <- make_panel_df(SURVIVAL_EDPMM_POST_condX0_001, "CN")

df_osqc_master <- bind_rows(df_overall, df_mci, df_cn) %>%
  filter(
    rho %in% rho_keep_tbl,
    nu  %in% nu_keep_tbl
  ) %>%
  distinct(subgroup, eta_z0, eta_z1, nu, rho, .keep_all = TRUE)

setting_levels <- build_setting_levels(
  df_pairs = df_osqc_master %>% distinct(eta_z0, eta_z1),
  preferred_pairs = preferred_other_pairs
)

df_osqc_master <- df_osqc_master %>%
  mutate(
    subgroup = factor(subgroup, levels = subgroup_order),
    is_primary = eta_z0 == 0 & eta_z1 == 0,
    setting_viewer = factor(
      mapply(label_pair_viewer, eta_z0, eta_z1),
      levels = setting_levels$viewer
    ),
    setting_latex = factor(
      mapply(label_pair_latex, eta_z0, eta_z1),
      levels = setting_levels$latex
    ),
    rho_lab = factor(
      paste0("rho = ", formatC(rho, format = "f", digits = 1)),
      levels = paste0("rho = ", formatC(rho_keep_tbl, format = "f", digits = 1))
    ),
    nu_lab = factor(
      paste0("nu = ", nu),
      levels = paste0("nu = ", nu_keep_tbl)
    ),
    mean_str = fmt_num(estimYs, 1),
    lwr_str  = fmt_num(quantile025, 1),
    upr_str  = fmt_num(quantile975, 1)
  ) %>%
  arrange(setting_viewer, subgroup, rho_lab, nu_lab)

# ============================================================
# CHECK: this should be empty for Overall
# ============================================================

safe_print_all(
  df_osqc_master %>%
    filter(
      subgroup == "Overall",
      nu == 0,
      (eta_z0 > 0 | eta_z1 > 0)
    ) %>%
    select(subgroup, eta_z0, eta_z1, rho, nu, estimYs, quantile025, quantile975)
)

# ============================================================
# TABLE HELPERS
# ============================================================

make_osqc_subcol_table <- function(df, row_vars) {
  nu_labs <- paste0("nu = ", nu_keep_tbl)
  
  df %>%
    mutate(
      subgroup = factor(as.character(subgroup), levels = subgroup_order),
      rho_lab  = factor(
        as.character(rho_lab),
        levels = paste0("rho = ", formatC(rho_keep_tbl, format = "f", digits = 1))
      ),
      nu_lab = factor(as.character(nu_lab), levels = nu_labs)
    ) %>%
    arrange(across(all_of(row_vars)), nu_lab) %>%
    select(all_of(row_vars), nu_lab, mean_str, lwr_str, upr_str) %>%
    distinct() %>%
    pivot_wider(
      names_from  = nu_lab,
      values_from = c(mean_str, lwr_str, upr_str),
      names_glue  = "{nu_lab}___{.value}"
    ) %>%
    select(
      all_of(row_vars),
      `nu = 0___mean_str`,  `nu = 0___lwr_str`,  `nu = 0___upr_str`,
      `nu = 5___mean_str`,  `nu = 5___lwr_str`,  `nu = 5___upr_str`,
      `nu = 10___mean_str`, `nu = 10___lwr_str`, `nu = 10___upr_str`
    )
}

show_osqc_subcol_table_viewer <- function(df, row_vars, table_title,
                                          font_size = 11,
                                          collapse_cols = integer(0)) {
  row_name_map <- c(
    subgroup       = "Group",
    rho_lab        = "rho",
    setting_viewer = "Setting",
    setting_latex  = "Setting"
  )
  
  col_names <- c(
    unname(row_name_map[row_vars]),
    rep(c("Mean", "0.5% CrI", "99.5% CrI"), 3)
  )
  
  out_tbl <- kbl(
    df,
    format = "html",
    escape = FALSE,
    align = "c",
    col.names = col_names,
    caption = table_title
  ) %>%
    kable_styling(
      bootstrap_options = c("hover", "condensed"),
      full_width = FALSE,
      font_size = font_size
    ) %>%
    add_header_above(c(" " = length(row_vars), "nu = 0" = 3, "nu = 5" = 3, "nu = 10" = 3)) %>%
    column_spec(seq_len(length(row_vars)), background = "white") %>%
    column_spec((length(row_vars) + 1):(length(row_vars) + 3), background = "grey95") %>%
    column_spec((length(row_vars) + 4):(length(row_vars) + 6), background = "grey88") %>%
    column_spec((length(row_vars) + 7):(length(row_vars) + 9), background = "grey95")
  
  if (length(collapse_cols) > 0) {
    out_tbl <- out_tbl %>%
      collapse_rows(columns = collapse_cols, valign = "top")
  }
  
  html_print(
    HTML(as.character(out_tbl)),
    viewer = getOption("viewer")
  )
  
  invisible(out_tbl)
}

blank_repeats_for_latex <- function(df, row_vars) {
  out <- df %>%
    mutate(across(all_of(row_vars), as.character))
  
  for (v in row_vars) {
    out[[paste0(v, "__orig")]] <- out[[v]]
  }
  
  for (i in seq_along(row_vars)) {
    cur_var  <- row_vars[i]
    cur_orig <- paste0(cur_var, "__orig")
    
    if (i == 1) {
      out[[cur_var]][duplicated(out[[cur_orig]])] <- ""
    } else {
      prev_orig <- paste0(row_vars[seq_len(i - 1)], "__orig")
      
      out <- out %>%
        group_by(across(all_of(prev_orig))) %>%
        mutate(
          !!cur_var := replace(
            .data[[cur_var]],
            duplicated(.data[[cur_orig]]),
            ""
          )
        ) %>%
        ungroup()
    }
  }
  
  out %>%
    select(-all_of(paste0(row_vars, "__orig")))
}

make_latex_linesep <- function(df, row_vars) {
  df_cmp <- df %>%
    mutate(across(all_of(row_vars), as.character))
  
  n_rows     <- nrow(df_cmp)
  total_cols <- ncol(df_cmp)
  linesep    <- rep("", n_rows)
  
  if (n_rows <= 1) {
    return(linesep)
  }
  
  for (i in seq_len(n_rows - 1)) {
    cur_vals  <- unlist(df_cmp[i, row_vars, drop = FALSE], use.names = FALSE)
    next_vals <- unlist(df_cmp[i + 1, row_vars, drop = FALSE], use.names = FALSE)
    
    diff_idx <- which(cur_vals != next_vals)
    
    if (length(diff_idx) > 0) {
      first_change <- diff_idx[1]
      linesep[i] <- sprintf("\\cmidrule{%d-%d}", first_change, total_cols)
    }
  }
  
  linesep
}

make_osqc_subcol_table_latex <- function(df, row_vars, table_title, label,
                                         font_size = 8) {
  row_name_map <- c(
    subgroup       = "Group",
    rho_lab        = "$\\rho$",
    setting_viewer = "Setting",
    setting_latex  = "Setting"
  )
  
  col_names <- c(
    unname(row_name_map[row_vars]),
    rep(c("Mean", "0.5\\% CrI", "99.5\\% CrI"), 3)
  )
  
  df_latex <- blank_repeats_for_latex(df, row_vars)
  linesep  <- make_latex_linesep(df, row_vars)
  
  out_tbl <- kbl(
    df_latex,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c(rep("l", length(row_vars)), rep("c", 9)),
    col.names = col_names,
    caption = table_title,
    label = label,
    linesep = linesep
  ) %>%
    add_header_above(
      c(" " = length(row_vars), "$\\nu=0$" = 3, "$\\nu=5$" = 3, "$\\nu=10$" = 3),
      escape = FALSE
    ) %>%
    kable_styling(
      latex_options = c("scale_down"),
      font_size = font_size
    )
  
  out_tbl <- as.character(out_tbl)
  out_tbl <- sub(
    "\\\\begin\\{table\\}",
    "\\\\begin\\{table\\}[hbp]",
    out_tbl
  )
  
  out_tbl
}

# ============================================================
# TABLE 1: PRIMARY OSQC
# ============================================================

tab1_df <- df_osqc_master %>%
  filter(is_primary) %>%
  make_osqc_subcol_table(
    row_vars = c("subgroup", "rho_lab")
  )

# safe_view(tab1_df)

if (interactive()) {
  show_osqc_subcol_table_viewer(
    df = tab1_df,
    row_vars = c("subgroup", "rho_lab"),
    table_title = "Posterior mean OSQC estimates under the primary analysis",
    collapse_cols = 1
  )
}

tab1_tex <- make_osqc_subcol_table_latex(
  df = tab1_df,
  row_vars = c("subgroup", "rho_lab"),
  # table_title = "Posterior mean OSQC estimates under the primary analysis (OSQC; $\\psi_{0}=0$, $\\psi_{1}=0$).",
  table_title = "",
  label = "tab:osqc_primary",
  font_size = 8
)

cat(tab1_tex, sep = "\n")

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(patchwork)

nu_keep  <- c(0, 5, 10)
rho_keep <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

# ------------------------------------------------------------------------
# Shared settings and helpers --------------------------------------------
# ------------------------------------------------------------------------

panel_levels <- c("Overall", "CN", "MCI")
model_levels <- c("EDPMM", "DPMM")

# ------------------------------------------------------------------------
# Consistent palette (print-safe / restrained)
# ------------------------------------------------------------------------
# Amyloid status
amyloid_cols <- c(
  "Non-Elevated" = "#2C7BB6",  # bluish green
  "Elevated"     = "#C44E52"   # vermillion
)

# Sensitivity settings
sens_levels <- c(
  "psi0=0, psi1=0",
  "psi0=-0.25, psi1=0.25",
  "psi0=0.25, psi1=-0.25"
)

sens_cols <- c(
  "psi0=0, psi1=0"         = "#000000", # baseline = black
  "psi0=-0.25, psi1=0.25" = "#2C7BB6",  # blue
  "psi0=0.25, psi1=-0.25" = "#C44E52"   # vermillion
)

sens_labels <- c(
  "psi0=0, psi1=0"         = "psi[0]==0*','~~psi[1]==0",
  "psi0=-0.25, psi1=0.25" = "psi[0]==-0.25*','~~psi[1]==0.25",
  "psi0=0.25, psi1=-0.25" = "psi[0]==0.25*','~~psi[1]==-0.25"
)

nu_labeller <- as_labeller(
  c(
    "0"  = "nu==0",
    "5"  = "nu==5",
    "10" = "nu==10"
  ),
  label_parsed
)

psi0_labeller <- as_labeller(
  c(
    "0.25"  = "psi[0]==0.25",
    "0"     = "psi[0]==0",
    "-0.25" = "psi[0]==-0.25"
  ),
  label_parsed
)

psi1_labeller <- as_labeller(
  c(
    "-0.25" = "psi[1]==-0.25",
    "0"     = "psi[1]==0",
    "0.25"  = "psi[1]==0.25"
  ),
  label_parsed
)

# -----------------------------
# Common OSQC result formatter
# -----------------------------
format_osqc_df <- function(res_obj) {
  df <- as.data.frame(res_obj$E_rho_nu_diff_result)
  colnames(df)[1:9] <- c(
    "eta_z0", "eta_z1", "nu", "rho",
    "estimYs", "sd", "quantile025", "quantile975", "CIlength95"
  )
  df
}

# -----------------------------
# Sensitivity labels
# -----------------------------
add_sens_label <- function(df, spaced = FALSE) {
  if (spaced) {
    df %>%
      mutate(
        sens = case_when(
          eta_z0 ==  0.00 & eta_z1 ==  0.00 ~ "psi0 = 0, psi1 = 0",
          eta_z0 == -0.25 & eta_z1 ==  0.25 ~ "psi0 = -0.25, psi1 = 0.25",
          eta_z0 ==  0.25 & eta_z1 == -0.25 ~ "psi0 = 0.25, psi1 = -0.25",
          TRUE ~ NA_character_
        )
      )
  } else {
    df %>%
      mutate(
        sens = case_when(
          eta_z0 ==  0.00 & eta_z1 ==  0.00 ~ "psi0=0, psi1=0",
          eta_z0 == -0.25 & eta_z1 ==  0.25 ~ "psi0=-0.25, psi1=0.25",
          eta_z0 ==  0.25 & eta_z1 == -0.25 ~ "psi0=0.25, psi1=-0.25",
          TRUE ~ NA_character_
        )
      )
  }
}

# -----------------------------
# Generic builders
# -----------------------------
make_panel_df <- function(res_obj, panel_label) {
  format_osqc_df(res_obj) %>%
    mutate(panel = panel_label)
}

make_subgroup_df <- function(res_obj, subgroup_label) {
  format_osqc_df(res_obj) %>%
    mutate(subgroup = subgroup_label)
}

make_model_df <- function(res_obj, panel_label, model_label) {
  format_osqc_df(res_obj) %>%
    mutate(
      panel = panel_label,
      model = model_label
    )
}

# -----------------------------
# Survival helpers
# -----------------------------
make_S_df <- function(post_ftn_res) {
  df0 <- as.data.frame(post_ftn_res$S_z0_result)
  df1 <- as.data.frame(post_ftn_res$S_z1_result)
  
  colnames(df0) <- c("t","estimYs","sd","quantile025","quantile975","CIlength95")
  colnames(df1) <- c("t","estimYs","sd","quantile025","quantile975","CIlength95")
  
  bind_rows(
    df0 %>% mutate(amyloid = "Non-Elevated"),
    df1 %>% mutate(amyloid = "Elevated")
  ) %>%
    mutate(
      amyloid = factor(amyloid, levels = c("Non-Elevated", "Elevated"))
    ) %>%
    arrange(amyloid, t)
}

make_cond_surv_df <- function(res, subgroup_label, amyloid_label) {
  df <- as.data.frame(res)
  colnames(df) <- c("t","estimYs","sd","quantile025","quantile975","CIlength95")
  
  df %>%
    mutate(
      subgroup = subgroup_label,
      amyloid  = amyloid_label
    )
}

plot_survival_estimator <- function(df_post, df_km, panel_subtitle) {
  ggplot() +
    geom_ribbon(
      data = df_post,
      aes(x = t, ymin = quantile025, ymax = quantile975, fill = amyloid),
      alpha = 0.18, colour = NA
    ) +
    geom_line(
      data = df_post,
      aes(x = t, y = estimYs, color = amyloid),
      linewidth = 1.15
    ) +
    geom_step(
      data = df_km,
      aes(x = t, y = surv, color = amyloid),
      linetype = "dashed", linewidth = 0.85, alpha = 0.8
    ) +
    scale_color_manual(values = amyloid_cols, name = "Amyloid status") +
    scale_fill_manual(values = amyloid_cols,  name = "Amyloid status") +
    labs(
      x = "Time (years)",
      y = "Probability of remaining dementia-free",
      subtitle = panel_subtitle
    ) +
    coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
    theme_bw(base_size = 12) +
    theme(
      plot.subtitle = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
}

# ------------------------------------------------------------------------
# 1. Dementia-free survival by amyloid status ----------------------------
# ------------------------------------------------------------------------

df_S_overall <- make_S_df(SURVIVAL_EDPMM_POST_ftn_results_001)

df_S_cond <- bind_rows(
  make_cond_surv_df(SURVIVAL_EDPMM_POST_ftn_condX0_results_001$S_z0_result, "CN",  "Non-Elevated"),
  make_cond_surv_df(SURVIVAL_EDPMM_POST_ftn_condX0_results_001$S_z1_result, "CN",  "Elevated"),
  make_cond_surv_df(SURVIVAL_EDPMM_POST_ftn_condX1_results_001$S_z0_result, "MCI", "Non-Elevated"),
  make_cond_surv_df(SURVIVAL_EDPMM_POST_ftn_condX1_results_001$S_z1_result, "MCI", "Elevated")
) %>%
  mutate(
    subgroup = factor(subgroup, levels = c("CN", "MCI")),
    amyloid  = factor(amyloid, levels = c("Non-Elevated", "Elevated"))
  ) %>%
  arrange(subgroup, amyloid, t)

km_fit_overall <- survfit(Surv(Tobs, Dobs) ~ Zobs, data = dfADNI)
km_sum_overall <- summary(km_fit_overall)

df_KM_overall <- data.frame(
  t      = km_sum_overall$time,
  surv   = km_sum_overall$surv,
  strata = km_sum_overall$strata
) %>%
  mutate(
    amyloid = factor(
      ifelse(grepl("=0", strata), "Non-Elevated", "Elevated"),
      levels = c("Non-Elevated", "Elevated")
    )
  ) %>%
  arrange(amyloid, t)

df_km_cond <- data.frame(
  Tobs = dfADNI$Tobs,
  Dobs = dfADNI$Dobs,
  Baseline_DX = dfADNI$Xobs[, "Baseline_DX"],
  Zobs = dfADNI$Zobs
) %>%
  mutate(
    subgroup = factor(Baseline_DX, levels = c(0, 1), labels = c("CN", "MCI")),
    amyloid  = factor(Zobs, levels = c(0, 1), labels = c("Non-Elevated", "Elevated"))
  ) %>%
  filter(!is.na(subgroup), !is.na(amyloid), !is.na(Tobs), !is.na(Dobs))

km_fit_cond <- survfit(Surv(Tobs, Dobs) ~ subgroup + amyloid, data = df_km_cond)
km_sum_cond <- summary(km_fit_cond)

df_KM_cond <- data.frame(
  t      = km_sum_cond$time,
  surv   = km_sum_cond$surv,
  strata = km_sum_cond$strata
) %>%
  mutate(
    strata   = as.character(strata),
    subgroup = ifelse(grepl("subgroup=CN", strata), "CN", "MCI"),
    amyloid  = ifelse(grepl("amyloid=Non-Elevated", strata), "Non-Elevated", "Elevated"),
    subgroup = factor(subgroup, levels = c("CN", "MCI")),
    amyloid  = factor(amyloid, levels = c("Non-Elevated", "Elevated"))
  ) %>%
  arrange(subgroup, amyloid, t)

EDPqrlPOST_survival_overall_001 <-
  plot_survival_estimator(
    df_post = df_S_overall,
    df_km   = df_KM_overall,
    panel_subtitle = "Overall"
  )

EDPqrlPOST_survival_CN_001 <-
  plot_survival_estimator(
    df_post = df_S_cond %>% filter(subgroup == "CN"),
    df_km   = df_KM_cond %>% filter(subgroup == "CN"),
    panel_subtitle = "CN"
  )

EDPqrlPOST_survival_MCI_001 <-
  plot_survival_estimator(
    df_post = df_S_cond %>% filter(subgroup == "MCI"),
    df_km   = df_KM_cond %>% filter(subgroup == "MCI"),
    panel_subtitle = "MCI"
  )

EDPqrlPOST_survival_all3_001 <-
  (EDPqrlPOST_survival_overall_001 +
     EDPqrlPOST_survival_MCI_001 +
     EDPqrlPOST_survival_CN_001) +
  plot_layout(ncol = 3, guides = "collect", axis_titles = "collect") +
  plot_annotation() &
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
print(EDPqrlPOST_survival_all3_001)
# quartz(type = "pdf", file = "EDPqrlPOST_survival_all3_001.pdf", width = 12, height = 4.5)
# print(EDPqrlPOST_survival_all3_001)
# dev.off()

# ------------------------------------------------------------------------
# 2. MAIN TEXT: OSQC vs rho for Overall + MCI only -----------------------
# ------------------------------------------------------------------------

df_plot_main <- bind_rows(
  make_panel_df(EDPMM_results_POST_001, "Overall"),
  make_panel_df(SURVIVAL_EDPMM_POST_condX1_001, "MCI")
) %>%
  add_sens_label(spaced = FALSE) %>%
  filter(!is.na(sens), nu %in% nu_keep) %>%
  mutate(
    sens = factor(sens, levels = sens_levels),
    panel = factor(panel, levels = c("Overall", "MCI")),
    nu = factor(nu, levels = c(0, 5, 10))
  )

ylim_main <- range(
  c(df_plot_main$quantile025, df_plot_main$quantile975),
  na.rm = TRUE
)

EDPqrlPOST_OSQCrho_nu_overall_MCI_001 <-
  ggplot(df_plot_main, aes(x = rho, y = estimYs, color = sens, fill = sens)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(aes(ymin = quantile025, ymax = quantile975),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(panel ~ nu, labeller = labeller(nu = nu_labeller)) +
  coord_cartesian(ylim = ylim_main) +
  scale_x_continuous(breaks = rho_keep) +
  scale_color_manual(
    values = sens_cols,
    breaks = sens_levels,
    labels = function(x) parse(text = sens_labels[x])
  ) +
  scale_fill_manual(
    values = sens_cols,
    breaks = sens_levels,
    labels = function(x) parse(text = sens_labels[x])
  ) +
  labs(
    x = expression(rho),
    y = "OSQC",
    color = NULL,
    fill  = NULL,
  ) +
  guides(fill = "none") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank()
  )
print(EDPqrlPOST_OSQCrho_nu_overall_MCI_001)
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCrho_nu_overall_MCI_001.pdf", width = 12, height = 4.5)
# print(EDPqrlPOST_OSQCrho_nu_overall_MCI_001)
# dev.off()

# ------------------------------------------------------------------------
# 3. SUPPLEMENT: OSQC vs rho for CN only ---------------------------------
# ------------------------------------------------------------------------

df_plot_CN <- make_panel_df(SURVIVAL_EDPMM_POST_condX0_001, "CN") %>%
  add_sens_label(spaced = FALSE) %>%
  filter(!is.na(sens), nu %in% nu_keep) %>%
  mutate(
    sens = factor(sens, levels = sens_levels),
    nu   = factor(nu, levels = c(0, 5, 10))
  )

ylim_CN <- range(
  c(df_plot_CN$quantile025, df_plot_CN$quantile975),
  na.rm = TRUE
)

EDPqrlPOST_OSQCrho_nu_CNonly_001 <-
  ggplot(df_plot_CN, aes(x = rho, y = estimYs, color = sens, fill = sens)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(
    aes(ymin = quantile025, ymax = quantile975),
    alpha = 0.12, colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ nu, nrow = 1, labeller = nu_labeller) +
  coord_cartesian(ylim = ylim_CN) +
  scale_x_continuous(breaks = rho_keep) +
  scale_color_manual(
    values = sens_cols,
    breaks = sens_levels,
    labels = function(x) parse(text = sens_labels[x])
  ) +
  scale_fill_manual(
    values = sens_cols,
    breaks = sens_levels,
    labels = function(x) parse(text = sens_labels[x])
  ) +
  labs(
    x = expression(rho),
    y = "OSQC",
    color = NULL,
    fill  = NULL,
  ) +
  guides(fill = "none") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank()
  )
print(EDPqrlPOST_OSQCrho_nu_CNonly_001)
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCrho_nu_CNonly_001.pdf", width = 12, height = 4.5)
# print(EDPqrlPOST_OSQCrho_nu_CNonly_001)
# dev.off()

# ------------------------------------------------------------------------
# 4. MCI only, nu = 5 only: 3 x 3 sensitivity grid -----------------------
# Each panel shows:
#   - panel-specific EDPMM curve/ribbon
#   - baseline (psi0,psi1)=(0,0) curve/ribbon overlaid for comparison
# ------------------------------------------------------------------------

df_mci_nu5 <- make_panel_df(SURVIVAL_EDPMM_POST_condX1_001, "MCI") %>%
  mutate(
    psi0 = round(eta_z0, 2),
    psi1 = round(eta_z1, 2)
  ) %>%
  filter(
    nu == 5,
    psi0 %in% c(-0.25, 0, 0.25),
    psi1 %in% c(-0.25, 0, 0.25)
  ) %>%
  add_sens_label(spaced = FALSE)

grid_all <- tidyr::crossing(
  psi0 = c(-0.25, 0, 0.25),
  psi1 = c(-0.25, 0, 0.25)
)

df_actual <- df_mci_nu5 %>%
  mutate(source = "Panel-specific")

df_baseline <- df_mci_nu5 %>%
  filter(psi0 == 0, psi1 == 0) %>%
  select(-psi0, -psi1) %>%
  tidyr::crossing(grid_all) %>%
  mutate(source = "Baseline (0,0)")

df_plot_MCI_nu5_grid <- bind_rows(df_actual, df_baseline) %>%
  mutate(
    psi0_f = factor(as.character(psi0), levels = c("0.25", "0", "-0.25")),
    psi1_f = factor(as.character(psi1), levels = c("-0.25", "0", "0.25")),
    source = factor(source, levels = c("Panel-specific", "Baseline (0,0)")),
    sens = factor(sens, levels = sens_levels)
  )

ylim_MCI_nu5_grid <- range(
  c(df_plot_MCI_nu5_grid$quantile025, df_plot_MCI_nu5_grid$quantile975),
  na.rm = TRUE
)

# fixed colors
col_baseline  <- "#1A9850"
fill_baseline <- "#1A9850"

col_panel  <- "#222222"
fill_panel <- "gray85"

EDPqrlPOST_OSQCrho_MCI_nu5_grid_compare0_001 <-
  ggplot(df_plot_MCI_nu5_grid, aes(x = rho, y = estimYs)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  
  # baseline ribbon repeated in every panel
  geom_ribbon(
    data = subset(df_plot_MCI_nu5_grid, source == "Baseline (0,0)"),
    aes(ymin = quantile025, ymax = quantile975),
    alpha = 0.30,
    fill = fill_baseline,
    colour = NA
  ) +
  
  # panel-specific sensitivity ribbon
  geom_ribbon(
    data = subset(df_plot_MCI_nu5_grid, source == "Panel-specific"),
    aes(ymin = quantile025, ymax = quantile975),
    alpha = 0.70,
    fill = fill_panel,
    colour = NA
  ) +
  
  # baseline line repeated in every panel
  geom_line(
    data = subset(df_plot_MCI_nu5_grid, source == "Baseline (0,0)"),
    color = col_baseline,
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  # panel-specific sensitivity line
  geom_line(
    data = subset(df_plot_MCI_nu5_grid, source == "Panel-specific"),
    color = col_panel,
    linetype = "solid",
    linewidth = 1.0
  ) +
  
  facet_grid(
    psi0_f ~ psi1_f,
    labeller = labeller(
      psi0_f = psi0_labeller,
      psi1_f = psi1_labeller
    ),
    switch = "y"
  ) +
  coord_cartesian(ylim = ylim_MCI_nu5_grid) +
  scale_x_continuous(breaks = rho_keep) +
  labs(
    x = expression(rho),
    y = "OSQC"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank()
  )
print(EDPqrlPOST_OSQCrho_MCI_nu5_grid_compare0_001)
# quartz(type = "pdf", file = "EDPqrlPOST_OSQCrho_MCI_nu5_grid_compare0_001.pdf", width = 12, height = 4.5)
# print(EDPqrlPOST_OSQCrho_MCI_nu5_grid_compare0_001)
# dev.off()

# ------------------------------------------------------------------------
# Figures: EDPMM, phi = 0.25 / 0 / -0.25
# Main text: Overall + MCI
# Supplement: CN only
# ------------------------------------------------------------------------

nu_keep  <- c(0, 5, 10)
rho_keep <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

# ------------------------------------------------------------------------
# Helper: extract OSQC summary table
# ------------------------------------------------------------------------
format_osqc_df <- function(res_obj) {
  df <- as.data.frame(res_obj$E_rho_nu_diff_result)
  colnames(df)[1:9] <- c(
    "eta_z0", "eta_z1", "nu", "rho",
    "estimYs", "sd", "quantile025", "quantile975", "CIlength95"
  )
  df
}

# ------------------------------------------------------------------------
# Helper: build plotting df
# zero_eta_only = TRUE is used ONLY for the old *_001 objects
# ------------------------------------------------------------------------
make_phi_df <- function(res_obj, panel_label, phi_label, zero_eta_only = FALSE) {
  df <- format_osqc_df(res_obj)
  
  if (zero_eta_only) {
    df <- df %>%
      filter(eta_z0 == 0, eta_z1 == 0)
  }
  
  df %>%
    arrange(nu, rho) %>%
    mutate(
      panel = panel_label,
      phi   = phi_label
    )
}

# ------------------------------------------------------------------------
# Facet / legend settings
# ------------------------------------------------------------------------
nu_labeller <- as_labeller(
  c(
    "0"  = "nu==0",
    "5"  = "nu==5",
    "10" = "nu==10"
  ),
  label_parsed
)

phi_levels <- c("0.25", "0", "-0.25")

phi_cols <- c(
  "0.25"  = "#C44E52",  # red
  "0"     = "#000000",  # black
  "-0.25" = "#2C7BB6"   # blue
)

phi_labels <- c(
  "0.25"  = "phi==0.25",
  "0"     = "phi==0",
  "-0.25" = "phi==-0.25"
)

# ------------------------------------------------------------------------
# Build data: MAIN (Overall + MCI)
#   phi = 0.25   -> use *_p025 objects directly
#   phi = 0      -> use *_001 objects, but filter eta_z0 = 0, eta_z1 = 0
#   phi = -0.25  -> use *_m025 objects directly
# ------------------------------------------------------------------------
df_main_phi <- bind_rows(
  make_phi_df(EDPMM_results_POST_p025,         "Overall", "0.25",  zero_eta_only = FALSE),
  make_phi_df(EDPMM_results_POST_001,          "Overall", "0",     zero_eta_only = TRUE),
  make_phi_df(EDPMM_results_POST_m025,         "Overall", "-0.25", zero_eta_only = FALSE),
  
  make_phi_df(SURVIVAL_EDPMM_POST_condX1_p025, "MCI",     "0.25",  zero_eta_only = FALSE),
  make_phi_df(SURVIVAL_EDPMM_POST_condX1_001,  "MCI",     "0",     zero_eta_only = TRUE),
  make_phi_df(SURVIVAL_EDPMM_POST_condX1_m025, "MCI",     "-0.25", zero_eta_only = FALSE)
) %>%
  filter(nu %in% nu_keep) %>%
  mutate(
    panel = factor(panel, levels = c("Overall", "MCI")),
    nu    = factor(nu, levels = c(0, 5, 10)),
    phi   = factor(phi, levels = phi_levels)
  )

ylim_main <- range(
  c(df_main_phi$quantile025, df_main_phi$quantile975),
  na.rm = TRUE
)

# ------------------------------------------------------------------------
# Main figure
# ------------------------------------------------------------------------
EDPMM_OSQCrho_nu_overall_MCI_phi <-
  ggplot(
    df_main_phi,
    aes(x = rho, y = estimYs, color = phi, fill = phi, group = phi)
  ) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(
    aes(ymin = quantile025, ymax = quantile975),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  facet_grid(panel ~ nu, labeller = labeller(nu = nu_labeller)) +
  coord_cartesian(ylim = ylim_main) +
  scale_x_continuous(breaks = rho_keep) +
  scale_color_manual(
    values = phi_cols,
    breaks = phi_levels,
    labels = function(x) parse(text = phi_labels[x])
  ) +
  scale_fill_manual(
    values = phi_cols,
    breaks = phi_levels,
    labels = function(x) parse(text = phi_labels[x])
  ) +
  labs(
    x = expression(rho),
    y = "OSQC",
    color = NULL,
    fill  = NULL,
  ) +
  guides(fill = "none") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank()
  )
print(EDPMM_OSQCrho_nu_overall_MCI_phi)
# quartz(type = "pdf",
#        file = "EDPMM_OSQCrho_nu_overall_MCI_phi.pdf",
#        width = 12, height = 4.5)
# print(EDPMM_OSQCrho_nu_overall_MCI_phi)
# dev.off()

# ------------------------------------------------------------------------
# Build data: SUPPLEMENT (CN only)
# ------------------------------------------------------------------------
df_CN_phi <- bind_rows(
  make_phi_df(SURVIVAL_EDPMM_POST_condX0_p025, "CN", "0.25",  zero_eta_only = FALSE),
  make_phi_df(SURVIVAL_EDPMM_POST_condX0_001,  "CN", "0",     zero_eta_only = TRUE),
  make_phi_df(SURVIVAL_EDPMM_POST_condX0_m025, "CN", "-0.25", zero_eta_only = FALSE)
) %>%
  filter(nu %in% nu_keep) %>%
  mutate(
    nu  = factor(nu, levels = c(0, 5, 10)),
    phi = factor(phi, levels = phi_levels)
  )

ylim_CN <- range(
  c(df_CN_phi$quantile025, df_CN_phi$quantile975),
  na.rm = TRUE
)

# ------------------------------------------------------------------------
# Supplement figure
# ------------------------------------------------------------------------
EDPMM_OSQCrho_nu_CN_phi <-
  ggplot(
    df_CN_phi,
    aes(x = rho, y = estimYs, color = phi, fill = phi, group = phi)
  ) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray40") +
  geom_ribbon(
    aes(ymin = quantile025, ymax = quantile975),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ nu, nrow = 1, labeller = nu_labeller) +
  coord_cartesian(ylim = ylim_CN) +
  scale_x_continuous(breaks = rho_keep) +
  scale_color_manual(
    values = phi_cols,
    breaks = phi_levels,
    labels = function(x) parse(text = phi_labels[x])
  ) +
  scale_fill_manual(
    values = phi_cols,
    breaks = phi_levels,
    labels = function(x) parse(text = phi_labels[x])
  ) +
  labs(
    x = expression(rho),
    y = "OSQC",
    color = NULL,
    fill  = NULL,
  ) +
  guides(fill = "none") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank()
  )
print(EDPMM_OSQCrho_nu_CN_phi)
# quartz(type = "pdf",
#        file = "EDPMM_OSQCrho_nu_CN_phi.pdf",
#        width = 12, height = 4.5)
# print(EDPMM_OSQCrho_nu_CN_phi)
# dev.off()
