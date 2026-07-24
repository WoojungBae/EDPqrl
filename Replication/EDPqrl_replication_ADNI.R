# ------------------------------------------------------------------------
# Directories ------------------------------------------------------------
# ------------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
data_dir = ""
source_dir = ""

# ------------------------------------------------------------------------
# Packages: load once only -----------------------------------------------
# ------------------------------------------------------------------------
library(Rcpp)
library(RcppArmadillo)

library(tidyverse)
library(survival)
library(viridis)
library(akima)
library(scales)
library(survminer)
library(grid)
library(patchwork)

library(knitr)
library(kableExtra)
library(htmltools)

library(ADNIMERGE2)

# ------------------------------------------------------------------------
# Source files: source once only -----------------------------------------
# ------------------------------------------------------------------------
if (nzchar(source_dir)) {
  setwd(source_dir)
}

source("EDPqrl_r.R")
sourceCpp("EDPqrl_cpp.cpp")

if (nzchar(data_dir)) {
  setwd(data_dir)
}

set.seed(1)

# ------------------------------------------------------------------------
# Checkpoint files -------------------------------------------------------
# ------------------------------------------------------------------------
checkpoint_dir = if (nzchar(data_dir)) data_dir else "."

checkpoint_mcmc  = file.path(checkpoint_dir, "ADNImcmc.RData")
checkpoint_post  = file.path(checkpoint_dir, "ADNImcmc+post.RData")
checkpoint_final = file.path(checkpoint_dir, "ADNI_full_analysis.RData")

# This script runs continuously from top to bottom.
#
# To resume from a checkpoint in a new R session:
#   1. Load the checkpoint.
#   2. Re-run the package/source-file block above.
#
# Re-sourcing EDPqrl_r.R and EDPqrl_cpp.cpp replaces stale Rcpp wrappers
# that may have been stored in the .RData checkpoint.

# ==============================================================================
# 1. Prepare ADNI analysis data
# ==============================================================================
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


# ==============================================================================
# 2. MCMC sampling
# ==============================================================================

# ------------------------------------------------------------------------
# Constants --------------------------------------------------------------
# ------------------------------------------------------------------------
gibbs_thin   = 2e1
gibbs_iter   = 4e4
gibbs_burnin = 2e4

# ------------------------------------------------------------------------
# EDPMM ------------------------------------------------------------------
# ------------------------------------------------------------------------
EDPMM_results_MCMC =
  SURVIVAL_EDPMM_MCMC(
    dfADNI,
    gibbs_iter,
    gibbs_burnin,
    gibbs_thin
  )

EDPMM_results_MCMC_e09_e11 =
  SURVIVAL_EDPMM_MCMC(
    dfADNI,
    gibbs_iter,
    gibbs_burnin,
    gibbs_thin,
    scale0 = 0.9,
    scale1 = 1.1
  )

EDPMM_results_MCMC_e11_e09 =
  SURVIVAL_EDPMM_MCMC(
    dfADNI,
    gibbs_iter,
    gibbs_burnin,
    gibbs_thin,
    scale0 = 1.1,
    scale1 = 0.9
  )

# ------------------------------------------------------------------------
# DPMM -------------------------------------------------------------------
# ------------------------------------------------------------------------
DPMM_results_MCMC =
  SURVIVAL_DPMM_MCMC(
    dfADNI,
    gibbs_iter,
    gibbs_burnin,
    gibbs_thin
  )

# ------------------------------------------------------------------------
# Checkpoint after MCMC --------------------------------------------------
# ------------------------------------------------------------------------
# save.image(file = checkpoint_mcmc)
load(file = checkpoint_mcmc)
message("Saved MCMC checkpoint: ", checkpoint_mcmc)

# ==============================================================================
# 3. POST-processing
# ==============================================================================
# ------------------------------------------------------------------------
# Settings ---------------------------------------------------------------
# ------------------------------------------------------------------------
level = 0.01
M = 1e3

# censoring = 1 - mean(dfADNI$Dobs) = 85.83397
# summary(dfADNI$Tobs)
# min: 0.0027
# max: 14.8419
# t_grid includes extrapolation beyond 14.9

nu = seq(0, 15, by = 0.5)
rho = seq(0.05, 0.30, by = 0.05)
t_grid = c(1e-3, seq(0.5, 20, by = 0.5))

# ============================================================================
# 0. Main benchmark analysis
# ============================================================================
# phi0 = 0, phi1 = 0
# psi0 = 0, psi1 = 0
# kappa = 0

psi_z0 = 0
psi_z1 = 0

# ------------------------------------------------------------------------
# Average PSQC
# ------------------------------------------------------------------------
EDPMM_results_POST_main =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = NULL,
    p_condX = NULL,
    level = level
  )

# ------------------------------------------------------------------------
# Average survival / density curves
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_ftn_results_main =
  SURVIVAL_EDPMM_POST_ftn(
    object = EDPMM_results_MCMC,
    M = M,
    t_grid = t_grid,
    ftn = c("S", "f"),
    condX = NULL,
    p_condX = NULL,
    esttype = "mean",
    level = level
  )

# ------------------------------------------------------------------------
# Conditional PSQC: Baseline_DX = 1
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_condX1_main =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 1,
    p_condX = 1,
    level = level
  )

# ------------------------------------------------------------------------
# Conditional PSQC: Baseline_DX = 0
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_condX0_main =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 0,
    p_condX = 1,
    level = level
  )

# ------------------------------------------------------------------------
# Conditional survival / density curves: Baseline_DX = 1
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_ftn_condX1_results_main =
  SURVIVAL_EDPMM_POST_ftn(
    object = EDPMM_results_MCMC,
    M = M,
    t_grid = t_grid,
    ftn = c("S", "f"),
    condX = 1,
    p_condX = 1,
    esttype = "mean",
    level = level
  )

# ------------------------------------------------------------------------
# Conditional survival / density curves: Baseline_DX = 0
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_ftn_condX0_results_main =
  SURVIVAL_EDPMM_POST_ftn(
    object = EDPMM_results_MCMC,
    M = M,
    t_grid = t_grid,
    ftn = c("S", "f"),
    condX = 0,
    p_condX = 1,
    esttype = "mean",
    level = level
  )

# ============================================================================
# 1. Unmeasured Confounding sensitivity
# ============================================================================
# phi0 = 0, phi1 = 0
# psi0, psi1 vary in post-processing
# kappa = 0

psi_z0 = c(-0.25, 0, 0.25)
psi_z1 = c(-0.25, 0, 0.25)

# ------------------------------------------------------------------------
# Average PSQC
# ------------------------------------------------------------------------
EDPMM_results_POST_uc =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = NULL,
    p_condX = NULL,
    level = level
  )

# ------------------------------------------------------------------------
# Conditional PSQC: Baseline_DX = 1
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_condX1_uc =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 1,
    p_condX = 1,
    level = level
  )

# ------------------------------------------------------------------------
# Conditional PSQC: Baseline_DX = 0
# ------------------------------------------------------------------------
SURVIVAL_EDPMM_POST_condX0_uc =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 0,
    p_condX = 1,
    level = level
  )

# ============================================================================
# 2. Sensitivity Analysis for Cross-World Dependence
# ============================================================================
# phi0 = 0, phi1 = 0
# psi0 = 0, psi1 = 0
# kappa varies in post-processing
#
# kappa = NULL corresponds to kappa = 0 and is already covered in main.
# Here we run two sensitivity settings:
#   kappa = "unif01"
#   kappa = "unifm11"

psi_z0 = 0
psi_z1 = 0

# ------------------------------------------------------------------------
# 2a. kappa ~ Uniform(0, 1)
# ------------------------------------------------------------------------
EDPMM_results_POST_kappa_unif01 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = "unif01",
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = NULL,
    p_condX = NULL,
    level = level
  )

SURVIVAL_EDPMM_POST_condX1_kappa_unif01 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = "unif01",
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 1,
    p_condX = 1,
    level = level
  )

SURVIVAL_EDPMM_POST_condX0_kappa_unif01 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = "unif01",
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 0,
    p_condX = 1,
    level = level
  )

# ------------------------------------------------------------------------
# 2b. kappa ~ Uniform(-1, 1)
# ------------------------------------------------------------------------
EDPMM_results_POST_kappa_unifm11 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = "unifm11",
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = NULL,
    p_condX = NULL,
    level = level
  )

SURVIVAL_EDPMM_POST_condX1_kappa_unifm11 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = "unifm11",
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 1,
    p_condX = 1,
    level = level
  )

SURVIVAL_EDPMM_POST_condX0_kappa_unifm11 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC,
    M = M,
    nu = nu,
    rho = rho,
    kappa = "unifm11",
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 0,
    p_condX = 1,
    level = level
  )

# ============================================================================
# 3. Informative censoring sensitivity via eta / scale
#    phi0 = phi1 = 0
#
#    eta0 = 0.90, eta1 = 0.90
#    eta0 = 0.90, eta1 = 1.10
#    eta0 = 1.10, eta1 = 0.90
#    eta0 = 1.10, eta1 = 1.10
# ============================================================================

psi_z0 = 0
psi_z1 = 0

# ------------------------------------------------------------------------
# 3a. eta0 = 0.90, eta1 = 1.10
# ------------------------------------------------------------------------

EDPMM_results_POST_eta_e09_e11 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC_e09_e11,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = NULL,
    p_condX = NULL,
    level = level
  )

SURVIVAL_EDPMM_POST_condX1_eta_e09_e11 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC_e09_e11,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 1,
    p_condX = 1,
    level = level
  )

SURVIVAL_EDPMM_POST_condX0_eta_e09_e11 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC_e09_e11,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 0,
    p_condX = 1,
    level = level
  )

# ------------------------------------------------------------------------
# 3b. eta0 = 1.10, eta1 = 0.90
# ------------------------------------------------------------------------

EDPMM_results_POST_eta_e11_e09 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC_e11_e09,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = NULL,
    p_condX = NULL,
    level = level
  )

SURVIVAL_EDPMM_POST_condX1_eta_e11_e09 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC_e11_e09,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 1,
    p_condX = 1,
    level = level
  )

SURVIVAL_EDPMM_POST_condX0_eta_e11_e09 =
  SURVIVAL_EDPMM_POST_parallel(
    object = EDPMM_results_MCMC_e11_e09,
    M = M,
    nu = nu,
    rho = rho,
    kappa = NULL,
    psi_z0 = psi_z0,
    psi_z1 = psi_z1,
    condX = 0,
    p_condX = 1,
    level = level
  )

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checkpoint after POST-processing ---------------------------------------
# ------------------------------------------------------------------------
# save.image(file = checkpoint_post)
load(file = checkpoint_post)
message("Saved POST checkpoint: ", checkpoint_post)

# ==============================================================================
# 4. Tables
# ==============================================================================
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Tables
#
# Manuscript order:
#   - Main-text Table 1 is the simulation table and is generated by the
#     separate simulation replication script.
#   - Main-text Table 2 below is the only ADNI analysis table in the PDF.
#   - Sensitivity analyses are retained in post-processing for the figures,
#     but no additional sensitivity tables are generated.
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

# ============================================================
# 1. Settings
# ============================================================

rho_keep_tbl <- c(0.10, 0.20, 0.30)
nu_keep_tbl  <- c(0, 2.5, 5)

subgroup_order <- c("Overall", "MCI", "CN")

preferred_uc_pairs <- tibble::tribble(
  ~psi_z0, ~psi_z1,
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
# 2. Clean only table objects
# ============================================================

rm(list = intersect(
  c(
    "df_primary_master",
    "df_uc_master",
    "df_kappa_master",
    "df_phi_master",
    
    "tab1_primary_df",
    "tab1_primary_tex",
    
    "tab2_uc_df",
    "tab2_uc_latex_df",
    "tab2_uc_tex",
    
    "tab3_kappa_df",
    "tab3_kappa_latex_df",
    "tab3_kappa_tex",
    
    "tab4_phi_df",
    "tab4_phi_latex_df",
    "tab4_phi_tex",
    "tab4_phi_stability_check"
  ),
  ls()
))

# ============================================================
# 3. Helpers
# ============================================================

fmt_num <- function(x, digits = 1) {
  ifelse(
    is.na(x) | is.nan(x) | !is.finite(x),
    "NA",
    formatC(x, format = "f", digits = digits)
  )
}

fmt_num_phi <- function(x, digits = 1, huge_cut = 1e4) {
  ifelse(
    is.na(x) | is.nan(x) | !is.finite(x),
    "NA",
    ifelse(
      abs(x) > huge_cut,
      "unstable",
      formatC(x, format = "f", digits = digits)
    )
  )
}

fmt_psi <- function(x, digits = 2) {
  out <- formatC(x, format = "f", digits = digits)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  out
}

label_pair_viewer <- function(psi0, psi1) {
  if (isTRUE(all.equal(psi0, 0)) && isTRUE(all.equal(psi1, 0))) {
    "psi0=0, psi1=0"
  } else {
    paste0("psi0=", fmt_psi(psi0), ", psi1=", fmt_psi(psi1))
  }
}

label_pair_latex <- function(psi0, psi1) {
  if (isTRUE(all.equal(psi0, 0)) && isTRUE(all.equal(psi1, 0))) {
    "$\\psi_{0}=0,\\ \\psi_{1}=0$"
  } else {
    paste0(
      "$\\psi_{0}=", fmt_psi(psi0),
      ",\\ \\psi_{1}=", fmt_psi(psi1), "$"
    )
  }
}

format_PSQC_df <- function(res_obj) {
  if (!is.list(res_obj) || !"E_rho_nu_diff_result" %in% names(res_obj)) {
    stop("Object does not contain E_rho_nu_diff_result.")
  }
  
  df <- as.data.frame(res_obj$E_rho_nu_diff_result)
  
  if (ncol(df) < 9) {
    stop("E_rho_nu_diff_result has fewer than 9 columns.")
  }
  
  colnames(df)[1:9] <- c(
    "psi_z0", "psi_z1", "nu", "rho",
    "estimYs", "sd", "quantile025", "quantile975", "CIlength95"
  )
  
  df %>%
    mutate(
      psi_z0 = as.numeric(psi_z0),
      psi_z1 = as.numeric(psi_z1),
      nu     = as.numeric(nu),
      rho    = round(as.numeric(rho), 2)
    )
}

make_panel_df <- function(res_obj, subgroup_label) {
  format_PSQC_df(res_obj) %>%
    mutate(subgroup = subgroup_label)
}

check_objects_exist <- function(obj_names) {
  missing_objs <- setdiff(obj_names, ls(envir = .GlobalEnv))
  
  if (length(missing_objs) > 0) {
    stop(
      "Missing required object(s):\n",
      paste(missing_objs, collapse = "\n"),
      "\n\nCannot build this table because these objects do not exist in the current R session."
    )
  }
  
  invisible(TRUE)
}

add_common_table_vars <- function(df, formatter = fmt_num) {
  df %>%
    filter(
      rho %in% rho_keep_tbl,
      nu  %in% nu_keep_tbl
    ) %>%
    mutate(
      subgroup = factor(subgroup, levels = subgroup_order),
      
      rho_lab = factor(
        paste0("rho = ", formatC(rho, format = "f", digits = 1)),
        levels = paste0("rho = ", formatC(rho_keep_tbl, format = "f", digits = 1))
      ),
      
      nu_lab = factor(
        paste0("nu = ", nu),
        levels = paste0("nu = ", nu_keep_tbl)
      ),
      
      mean_str = formatter(estimYs, 1),
      lwr_str  = formatter(quantile025, 1),
      upr_str  = formatter(quantile975, 1)
    )
}

make_PSQC_subcol_table <- function(df, row_vars) {
  expected_cols <- c(
    "nu = 0___mean_str",   "nu = 0___lwr_str",   "nu = 0___upr_str",
    "nu = 2.5___mean_str", "nu = 2.5___lwr_str", "nu = 2.5___upr_str",
    "nu = 5___mean_str",   "nu = 5___lwr_str",   "nu = 5___upr_str"
  )
  
  out <- df %>%
    arrange(across(all_of(row_vars)), nu_lab) %>%
    select(all_of(row_vars), nu_lab, mean_str, lwr_str, upr_str) %>%
    distinct() %>%
    pivot_wider(
      names_from  = nu_lab,
      values_from = c(mean_str, lwr_str, upr_str),
      names_glue  = "{nu_lab}___{.value}"
    )
  
  missing_cols <- setdiff(expected_cols, names(out))
  
  if (length(missing_cols) > 0) {
    for (cc in missing_cols) {
      out[[cc]] <- "NA"
    }
  }
  
  out %>%
    select(
      all_of(row_vars),
      `nu = 0___mean_str`,   `nu = 0___lwr_str`,   `nu = 0___upr_str`,
      `nu = 2.5___mean_str`, `nu = 2.5___lwr_str`, `nu = 2.5___upr_str`,
      `nu = 5___mean_str`,   `nu = 5___lwr_str`,   `nu = 5___upr_str`
    )
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

make_PSQC_subcol_table_latex <- function(df, row_vars, label, font_size = 8) {
  row_name_map <- c(
    subgroup      = "Group",
    rho_lab       = "$\\rho$",
    setting_latex = "Setting",
    kappa_latex   = "$\\kappa$",
    phi_latex     = "$\\phi$"
  )
  
  col_names <- c(
    unname(row_name_map[row_vars]),
    rep(c("Mean", "0.5\\% CrI", "99.5\\% CrI"), 3)
  )
  
  df_latex <- blank_repeats_for_latex(df, row_vars)
  
  first_data_col <- length(row_vars) + 1
  data_cols <- names(df_latex)[first_data_col:ncol(df_latex)]
  
  for (j in 1:3) {
    df_latex[[data_cols[j]]] <- kableExtra::cell_spec(
      df_latex[[data_cols[j]]],
      format = "latex",
      background = "gray!5"
    )
  }
  
  for (j in 4:6) {
    df_latex[[data_cols[j]]] <- kableExtra::cell_spec(
      df_latex[[data_cols[j]]],
      format = "latex",
      background = "gray!10"
    )
  }
  
  for (j in 7:9) {
    df_latex[[data_cols[j]]] <- kableExtra::cell_spec(
      df_latex[[data_cols[j]]],
      format = "latex",
      background = "gray!15"
    )
  }
  
  out_tbl <- kableExtra::kbl(
    df_latex,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c(rep("l", length(row_vars)), rep("c", 9)),
    col.names = col_names,
    caption = "",
    label = label,
    linesep = ""
  ) %>%
    kableExtra::add_header_above(
      c(" " = length(row_vars), "$\\nu=0$" = 3, "$\\nu=2.5$" = 3, "$\\nu=5$" = 3),
      escape = FALSE
    ) %>%
    kableExtra::kable_styling(
      latex_options = c("scale_down"),
      font_size = font_size
    )
  
  out_tbl <- as.character(out_tbl)
  
  out_tbl <- sub(
    "\\\\begin\\{table\\}",
    "\\\\begin{table}[hbp]",
    out_tbl
  )
  
  out_tbl
}

show_table_viewer <- function(df, row_vars, table_title,
                              font_size = 11,
                              collapse_cols = integer(0)) {
  row_name_map <- c(
    subgroup       = "Group",
    rho_lab        = "rho",
    setting_viewer = "Setting",
    kappa_viewer   = "kappa",
    phi_viewer     = "phi"
  )
  
  col_names <- c(
    unname(row_name_map[row_vars]),
    rep(c("Mean", "0.5% CrI", "99.5% CrI"), 3)
  )
  
  out_tbl <- kableExtra::kbl(
    df,
    format = "html",
    escape = FALSE,
    align = "c",
    col.names = col_names,
    caption = table_title
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "condensed"),
      full_width = FALSE,
      font_size = font_size
    ) %>%
    kableExtra::add_header_above(
      c(" " = length(row_vars), "nu = 0" = 3, "nu = 2.5" = 3, "nu = 5" = 3)
    )
  
  if (length(collapse_cols) > 0) {
    out_tbl <- out_tbl %>%
      kableExtra::collapse_rows(columns = collapse_cols, valign = "top")
  }
  
  htmltools::html_print(
    htmltools::HTML(as.character(out_tbl)),
    viewer = getOption("viewer")
  )
  
  invisible(out_tbl)
}

# ============================================================
# Main-text Table 2: Primary PSQC under the benchmark assumptions
# ============================================================

primary_objs <- c(
  "EDPMM_results_POST_main",
  "SURVIVAL_EDPMM_POST_condX1_main",
  "SURVIVAL_EDPMM_POST_condX0_main"
)

check_objects_exist(primary_objs)

df_primary_master <- bind_rows(
  make_panel_df(EDPMM_results_POST_main, "Overall"),
  make_panel_df(SURVIVAL_EDPMM_POST_condX1_main, "MCI"),
  make_panel_df(SURVIVAL_EDPMM_POST_condX0_main, "CN")
) %>%
  add_common_table_vars() %>%
  distinct(subgroup, psi_z0, psi_z1, nu, rho, .keep_all = TRUE) %>%
  arrange(subgroup, rho_lab, nu_lab)

tab1_primary_df <- df_primary_master %>%
  make_PSQC_subcol_table(
    row_vars = c("subgroup", "rho_lab")
  )

tab1_primary_tex <- make_PSQC_subcol_table_latex(
  df = tab1_primary_df,
  row_vars = c("subgroup", "rho_lab"),
  label = "tab:PSQC_primary",
  font_size = 8
)

cat("\n\n================ MAIN-TEXT TABLE 2: PRIMARY DATA ================\n")
print(as_tibble(tab1_primary_df), n = Inf)

cat("\n\n================ MAIN-TEXT TABLE 2: PRIMARY LATEX ================\n")
cat(tab1_primary_tex, sep = "\n")

if (interactive()) {
  show_table_viewer(
    df = tab1_primary_df,
    row_vars = c("subgroup", "rho_lab"),
    table_title = "Primary PSQC estimates",
    collapse_cols = 1
  )
}

# ------------------------------------------------------------------------
# No additional sensitivity tables are generated.
# The corresponding post-processing results are used below to create
# Main-text Figure 3 and Supplementary Figures S2 and S3.
# ------------------------------------------------------------------------

# ==============================================================================
# 5. Figures
# ==============================================================================
# Manuscript figure order generated by this script:
#   - Main-text Figure 1a: unadjusted Kaplan--Meier estimates
#   - Main-text Figure 1b: observed follow-up-time densities
#   - Main-text Figure 2 is the causal DAG and is not generated in this R script
#   - Main-text Figure 3: MCI PSQC sensitivity analysis
#   - Supplementary Figure S1: fitted survival and Kaplan--Meier curves
#   - Supplementary Figure S2: Overall PSQC sensitivity analysis
#   - Supplementary Figure S3: CN PSQC sensitivity analysis
#
# The additional 90% survival-crossing plot is retained at the end as an
# unnumbered diagnostic figure and is not part of the current manuscript.
# ==============================================================================

# ------------------------------------------------------------------------
# Main-text Figure 1a: Unadjusted Kaplan--Meier estimates of
# dementia-free survival by baseline amyloid status
# ------------------------------------------------------------------------
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

quartz(type = "pdf", file = "EDPqrl_KM.pdf", width = 8, height = 6)
print(km_plot)
dev.off()

# ------------------------------------------------------------------------
# Main-text Figure 1b: Kernel density estimates of observed
# dementia-free follow-up time, overall and by amyloid status
# ------------------------------------------------------------------------
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

# # 2. Generate the plot
# EDPqrl_DENSITYlogT <- ggplot(df_plot, aes(x = log_Tobs, color = Group, linetype = Group)) +
#   # --- 1. Fix the Legend Squares ---
#   geom_density(linewidth = 1.2, key_glyph = "path") + 
#   coord_cartesian(xlim = c(-6, 3), ylim = c(0, 0.5)) +
#   # --- 2. Color Updates (Red/Blue Palette) & Labels ---
#   scale_color_manual(
#     values = c(
#       "Overall" = "black", 
#       "Z=0, D=0 (Non-Elevated, Censored)" = "lightcoral",  # Light Red
#       "Z=1, D=0 (Elevated, Censored)" = "skyblue",       # Light Blue
#       "Z=0, D=1 (Non-Elevated, Event)" = "red",          # Solid Red
#       "Z=1, D=1 (Elevated, Event)" = "blue"              # Solid Blue
#     ),
#     labels = c(
#       "Overall" = "Overall", 
#       "Z=0, D=0 (Non-Elevated, Censored)" = "Non-Elevated, Censored (Z=0, D=0)", 
#       "Z=1, D=0 (Elevated, Censored)" = "Elevated, Censored (Z=1, D=0)", 
#       "Z=0, D=1 (Non-Elevated, Event)" = "Non-Elevated, Event (Z=0, D=1)", 
#       "Z=1, D=1 (Elevated, Event)" = "Elevated, Event (Z=1, D=1)"
#     )
#   ) +
#   # --- Distinct Linetypes & Updated Labels ---
#   scale_linetype_manual(
#     values = c(
#       "Overall" = "solid", 
#       "Z=0, D=0 (Non-Elevated, Censored)" = "dashed", 
#       "Z=1, D=0 (Elevated, Censored)" = "dashed", 
#       "Z=0, D=1 (Non-Elevated, Event)" = "dotted", 
#       "Z=1, D=1 (Elevated, Event)" = "dotted"
#     ),
#     labels = c(
#       "Overall" = "Overall", 
#       "Z=0, D=0 (Non-Elevated, Censored)" = "Non-Elevated, Censored (Z=0, D=0)", 
#       "Z=1, D=0 (Elevated, Censored)" = "Elevated, Censored (Z=1, D=0)", 
#       "Z=0, D=1 (Non-Elevated, Event)" = "Non-Elevated, Event (Z=0, D=1)", 
#       "Z=1, D=1 (Elevated, Event)" = "Elevated, Event (Z=1, D=1)"
#     )
#   ) +
#   theme_classic() + 
#   # --- Labels ---
#   labs(
#     title = "Density of Log-Transformed Observed Survival Time",
#     x = "Log-Transformed Time to Dementia", 
#     y = "Density",
#     color = "Subgroup",
#     linetype = "Subgroup" 
#   ) +
#   # --- 3. Legend Placement (Top-Left) ---
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     legend.position = c(0.02, 0.98), 
#     legend.justification = c(0, 1),    
#     legend.background = element_rect(fill = alpha("white", 0.85), color = "black", linewidth = 0.5),
#     legend.key.width = unit(3, "lines"), 
#     legend.text = element_text(size = 9)
#   ) +
#   guides(
#     color = guide_legend(ncol = 1),
#     linetype = guide_legend(ncol = 1)
#   )
# print(EDPqrl_DENSITYlogT)
# 
# quartz(type = "pdf", file = "EDPqrl_DENSITYlogT.pdf", width = 8, height = 6)
# print(EDPqrl_DENSITYlogT)
# dev.off()

# ------------------------------------------------------------------------
# Optional save settings -------------------------------------------------
# ------------------------------------------------------------------------

save_figures <- TRUE
fig_dir <- if (nzchar(data_dir)) data_dir else "."

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_figure_pdf <- function(plot_obj, file_stub, width = 12.0, height = 7.2) {
  
  if (!save_figures) {
    message("save_figures is FALSE. Not saving: ", file_stub)
    return(invisible(NULL))
  }
  
  out_file <- file.path(fig_dir, paste0(file_stub, ".pdf"))
  
  if (file.exists(out_file)) {
    unlink(out_file)
  }
  
  message("Saving PDF to: ", normalizePath(out_file, mustWork = FALSE))
  
  grDevices::quartz(
    type = "pdf",
    file = out_file,
    width = width,
    height = height
  )
  
  print(plot_obj)
  grDevices::dev.off()
  
  if (!file.exists(out_file)) {
    stop("PDF was not created: ", out_file)
  }
  
  message("Saved PDF: ", normalizePath(out_file))
  invisible(out_file)
}

# ------------------------------------------------------------------------
# Figure settings --------------------------------------------------------
# ------------------------------------------------------------------------

nu_keep  <- c(0, 2.5, 5)
rho_keep <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

nu_labeller <- ggplot2::as_labeller(
  c(
    "0"   = "nu==0",
    "2.5" = "nu==2.5",
    "5"   = "nu==5"
  ),
  label_parsed
)

sens_type_levels <- c("UC", "kappa", "eta")

row_labeller_short <- ggplot2::as_labeller(
  c(
    "UC"    = "psi",
    "kappa" = "kappa",
    "eta"   = "eta"
  ),
  label_parsed
)

# ------------------------------------------------------------------------
# Colors and labels ------------------------------------------------------
# ------------------------------------------------------------------------

amyloid_cols <- c(
  "Non-Elevated" = "#2C7BB6",
  "Elevated"     = "#C44E52"
)

# psi sensitivity
uc_levels <- c(
  "psi0=0, psi1=0",
  "psi0=-0.25, psi1=0.25",
  "psi0=0.25, psi1=-0.25"
)

uc_cols <- c(
  "psi0=0, psi1=0"         = "#000000",
  "psi0=-0.25, psi1=0.25" = "#2C7BB6",
  "psi0=0.25, psi1=-0.25" = "#C44E52"
)

uc_labels_plain <- c(
  "psi0=0, psi1=0"         = "ψ0 = 0, ψ1 = 0",
  "psi0=-0.25, psi1=0.25" = "ψ0 = -0.25, ψ1 = 0.25",
  "psi0=0.25, psi1=-0.25" = "ψ0 = 0.25, ψ1 = -0.25"
)

# kappa sensitivity
kappa_levels <- c(
  "kappa=0",
  "kappa~U(0,1)",
  "kappa~U(-1,1)"
)

kappa_cols <- c(
  "kappa=0"       = "#000000",
  "kappa~U(0,1)"  = "#2C7BB6",
  "kappa~U(-1,1)" = "#C44E52"
)

kappa_labels_plain <- c(
  "kappa=0"       = "κ = 0",
  "kappa~U(0,1)"  = "κ ~ U(0,1)",
  "kappa~U(-1,1)" = "κ ~ U(-1,1)"
)

# eta sensitivity
eta_levels <- c(
  "eta0=1.00, eta1=1.00",
  "eta0=1.10, eta1=0.90",
  "eta0=0.90, eta1=1.10"
)

eta_cols <- c(
  "eta0=1.00, eta1=1.00" = "#000000",
  "eta0=1.10, eta1=0.90" = "#C44E52",
  "eta0=0.90, eta1=1.10" = "#2C7BB6"
)

eta_labels_plain <- c(
  "eta0=1.00, eta1=1.00" = "η0 = 1, η1 = 1",
  "eta0=1.10, eta1=0.90" = "η0 = 1.1, η1 = 0.9",
  "eta0=0.90, eta1=1.10" = "η0 = 0.9, η1 = 1.1"
)

all_setting_levels <- c(
  uc_levels,
  kappa_levels,
  eta_levels
)

all_setting_cols <- c(
  uc_cols,
  kappa_cols,
  eta_cols
)

all_setting_labels_plain <- c(
  uc_labels_plain,
  kappa_labels_plain,
  eta_labels_plain
)

# ------------------------------------------------------------------------
# Object check -----------------------------------------------------------
# ------------------------------------------------------------------------

figure_objs <- c(
  "SURVIVAL_EDPMM_POST_ftn_results_main",
  "SURVIVAL_EDPMM_POST_ftn_condX1_results_main",
  "SURVIVAL_EDPMM_POST_ftn_condX0_results_main",
  
  "EDPMM_results_POST_main",
  "SURVIVAL_EDPMM_POST_condX1_main",
  "SURVIVAL_EDPMM_POST_condX0_main",
  
  "EDPMM_results_POST_uc",
  "SURVIVAL_EDPMM_POST_condX1_uc",
  "SURVIVAL_EDPMM_POST_condX0_uc",
  
  "EDPMM_results_POST_kappa_unif01",
  "SURVIVAL_EDPMM_POST_condX1_kappa_unif01",
  "SURVIVAL_EDPMM_POST_condX0_kappa_unif01",
  "EDPMM_results_POST_kappa_unifm11",
  "SURVIVAL_EDPMM_POST_condX1_kappa_unifm11",
  "SURVIVAL_EDPMM_POST_condX0_kappa_unifm11",
  
  "EDPMM_results_POST_eta_e09_e11",
  "SURVIVAL_EDPMM_POST_condX1_eta_e09_e11",
  "SURVIVAL_EDPMM_POST_condX0_eta_e09_e11",
  "EDPMM_results_POST_eta_e11_e09",
  "SURVIVAL_EDPMM_POST_condX1_eta_e11_e09",
  "SURVIVAL_EDPMM_POST_condX0_eta_e11_e09",
  
  "dfADNI"
)

missing_figure_objs <- setdiff(figure_objs, ls(envir = .GlobalEnv))

if (length(missing_figure_objs) > 0) {
  stop(
    "Missing required figure object(s):\n",
    paste(missing_figure_objs, collapse = "\n")
  )
}

# ------------------------------------------------------------------------
# Shared helpers ---------------------------------------------------------
# ------------------------------------------------------------------------

format_PSQC_df <- function(res_obj) {
  if (!is.list(res_obj) || !"E_rho_nu_diff_result" %in% names(res_obj)) {
    stop("Object does not contain E_rho_nu_diff_result.")
  }
  
  df <- as.data.frame(res_obj$E_rho_nu_diff_result)
  
  if (ncol(df) < 9) {
    stop("E_rho_nu_diff_result has fewer than 9 columns.")
  }
  
  colnames(df)[1:9] <- c(
    "psi_z0", "psi_z1", "nu", "rho",
    "estimYs", "sd", "quantile025", "quantile975", "CIlength95"
  )
  
  df %>%
    dplyr::mutate(
      psi_z0 = as.numeric(psi_z0),
      psi_z1 = as.numeric(psi_z1),
      nu     = as.numeric(nu),
      rho    = round(as.numeric(rho), 2),
      estimYs = as.numeric(estimYs),
      sd = as.numeric(sd),
      quantile025 = as.numeric(quantile025),
      quantile975 = as.numeric(quantile975),
      CIlength95 = as.numeric(CIlength95)
    )
}

make_panel_df <- function(res_obj, panel_label) {
  format_PSQC_df(res_obj) %>%
    dplyr::mutate(panel = panel_label)
}

make_S_df <- function(post_ftn_res) {
  df0 <- as.data.frame(post_ftn_res$S_z0_result)
  df1 <- as.data.frame(post_ftn_res$S_z1_result)
  
  colnames(df0) <- c("t", "estimYs", "sd", "quantile025", "quantile975", "CIlength95")
  colnames(df1) <- c("t", "estimYs", "sd", "quantile025", "quantile975", "CIlength95")
  
  dplyr::bind_rows(
    df0 %>% dplyr::mutate(amyloid = "Non-Elevated"),
    df1 %>% dplyr::mutate(amyloid = "Elevated")
  ) %>%
    dplyr::mutate(
      t = as.numeric(t),
      estimYs = as.numeric(estimYs),
      quantile025 = as.numeric(quantile025),
      quantile975 = as.numeric(quantile975),
      amyloid = factor(amyloid, levels = c("Non-Elevated", "Elevated"))
    ) %>%
    dplyr::arrange(amyloid, t)
}

make_cond_surv_df <- function(res, subgroup_label, amyloid_label) {
  df <- as.data.frame(res)
  colnames(df) <- c("t", "estimYs", "sd", "quantile025", "quantile975", "CIlength95")
  
  df %>%
    dplyr::mutate(
      t = as.numeric(t),
      estimYs = as.numeric(estimYs),
      quantile025 = as.numeric(quantile025),
      quantile975 = as.numeric(quantile975),
      subgroup = subgroup_label,
      amyloid  = amyloid_label
    )
}

# ------------------------------------------------------------------------
# Survival plot helper ---------------------------------------------------
# ------------------------------------------------------------------------

plot_survival_estimator <- function(df_post, df_km, panel_subtitle) {
  ggplot() +
    geom_ribbon(
      data = df_post,
      aes(x = t, ymin = quantile025, ymax = quantile975, fill = amyloid),
      alpha = 0.18,
      colour = NA
    ) +
    geom_line(
      data = df_post,
      aes(x = t, y = estimYs, color = amyloid),
      linewidth = 1.15
    ) +
    geom_step(
      data = df_km,
      aes(x = t, y = surv, color = amyloid),
      linetype = "dashed",
      linewidth = 0.85,
      alpha = 0.8
    ) +
    scale_color_manual(values = amyloid_cols, name = "Amyloid status") +
    scale_fill_manual(values = amyloid_cols, name = "Amyloid status") +
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
      plot.margin = margin(5.5, 8, 5.5, 5.5),
      panel.grid.minor = element_blank()
    )
}

# ------------------------------------------------------------------------
# Data builders for sensitivity figures ---------------------------------
# ------------------------------------------------------------------------

make_uc_df <- function(main_obj, uc_obj, panel_label) {
  dplyr::bind_rows(
    make_panel_df(main_obj, panel_label) %>%
      dplyr::filter(psi_z0 == 0, psi_z1 == 0) %>%
      dplyr::mutate(sens = "psi0=0, psi1=0"),
    
    make_panel_df(uc_obj, panel_label) %>%
      dplyr::mutate(
        sens = dplyr::case_when(
          psi_z0 == -0.25 & psi_z1 ==  0.25 ~ "psi0=-0.25, psi1=0.25",
          psi_z0 ==  0.25 & psi_z1 == -0.25 ~ "psi0=0.25, psi1=-0.25",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(sens))
  ) %>%
    dplyr::mutate(
      sens = factor(sens, levels = uc_levels),
      panel = factor(panel, levels = c("Overall", "MCI", "CN"))
    )
}

make_kappa_df <- function(main_obj, kappa01_obj, kappam11_obj, panel_label) {
  dplyr::bind_rows(
    make_panel_df(main_obj, panel_label) %>%
      dplyr::mutate(kappa = "kappa=0"),
    
    make_panel_df(kappa01_obj, panel_label) %>%
      dplyr::mutate(kappa = "kappa~U(0,1)"),
    
    make_panel_df(kappam11_obj, panel_label) %>%
      dplyr::mutate(kappa = "kappa~U(-1,1)")
  ) %>%
    dplyr::filter(psi_z0 == 0, psi_z1 == 0) %>%
    dplyr::mutate(
      kappa = factor(kappa, levels = kappa_levels),
      panel = factor(panel, levels = c("Overall", "MCI", "CN"))
    )
}

make_eta_df <- function(main_obj, eta_e11_e09_obj, eta_e09_e11_obj, panel_label) {
  dplyr::bind_rows(
    make_panel_df(main_obj, panel_label) %>%
      dplyr::filter(psi_z0 == 0, psi_z1 == 0) %>%
      dplyr::mutate(eta = "eta0=1.00, eta1=1.00"),
    
    make_panel_df(eta_e11_e09_obj, panel_label) %>%
      dplyr::filter(psi_z0 == 0, psi_z1 == 0) %>%
      dplyr::mutate(eta = "eta0=1.10, eta1=0.90"),
    
    make_panel_df(eta_e09_e11_obj, panel_label) %>%
      dplyr::filter(psi_z0 == 0, psi_z1 == 0) %>%
      dplyr::mutate(eta = "eta0=0.90, eta1=1.10")
  ) %>%
    dplyr::mutate(
      eta = factor(eta, levels = eta_levels),
      panel = factor(panel, levels = c("Overall", "MCI", "CN"))
    )
}

build_compact_sens_df_for_group <- function(group_label,
                                            main_obj,
                                            uc_obj,
                                            kappa01_obj,
                                            kappam11_obj,
                                            eta_e11_e09_obj,
                                            eta_e09_e11_obj) {
  
  df_uc <- make_uc_df(
    main_obj = main_obj,
    uc_obj = uc_obj,
    panel_label = group_label
  ) %>%
    dplyr::mutate(
      sens_type = "UC",
      setting = as.character(sens)
    ) %>%
    dplyr::select(-dplyr::any_of(c("sens", "kappa", "eta")))
  
  df_kappa <- make_kappa_df(
    main_obj = main_obj,
    kappa01_obj = kappa01_obj,
    kappam11_obj = kappam11_obj,
    panel_label = group_label
  ) %>%
    dplyr::mutate(
      sens_type = "kappa",
      setting = as.character(kappa)
    ) %>%
    dplyr::select(-dplyr::any_of(c("sens", "kappa", "eta")))
  
  df_eta <- make_eta_df(
    main_obj = main_obj,
    eta_e11_e09_obj = eta_e11_e09_obj,
    eta_e09_e11_obj = eta_e09_e11_obj,
    panel_label = group_label
  ) %>%
    dplyr::mutate(
      sens_type = "eta",
      setting = as.character(eta)
    ) %>%
    dplyr::select(-dplyr::any_of(c("sens", "kappa", "eta")))
  
  dplyr::bind_rows(df_uc, df_kappa, df_eta) %>%
    dplyr::filter(
      nu %in% nu_keep,
      rho %in% rho_keep
    ) %>%
    dplyr::mutate(
      panel = group_label,
      nu = factor(nu, levels = nu_keep),
      sens_type = factor(sens_type, levels = sens_type_levels),
      setting = factor(setting, levels = all_setting_levels)
    )
}

# ------------------------------------------------------------------------
# Main sensitivity figure helper -----------------------------------------
# ------------------------------------------------------------------------

plot_compact_sens_figure <- function(df_plot, group_label = NULL) {
  
  df_plot <- df_plot %>%
    dplyr::filter(
      nu %in% nu_keep,
      rho %in% rho_keep
    ) %>%
    dplyr::mutate(
      nu        = factor(nu, levels = nu_keep),
      sens_type = factor(sens_type, levels = sens_type_levels),
      setting   = factor(setting, levels = all_setting_levels)
    )
  
  # ----------------------------------------------------------------------
  # Left side: PSQC curves
  # ----------------------------------------------------------------------
  
  p_left <- ggplot(
    df_plot,
    aes(
      x = rho,
      y = estimYs,
      color = setting,
      fill = setting,
      group = setting
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2,
      color = "gray45"
    ) +
    geom_ribbon(
      aes(ymin = quantile025, ymax = quantile975),
      alpha = 0.12,
      colour = NA
    ) +
    geom_line(linewidth = 0.95) +
    facet_grid(
      sens_type ~ nu,
      scales = "free_y",
      switch = "y",
      labeller = labeller(
        sens_type = row_labeller_short,
        nu = nu_labeller
      )
    ) +
    scale_x_continuous(
      breaks = rho_keep
    ) +
    scale_color_manual(
      values = all_setting_cols,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = all_setting_cols,
      drop = FALSE
    ) +
    labs(
      x = expression(rho),
      y = "PSQC"
    ) +
    guides(
      color = "none",
      fill = "none"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold", size = 10),
      strip.placement = "outside",
      
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.35, "lines"),
      panel.spacing.y = unit(0.35, "lines"),
      
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text = element_text(size = 10),
      
      plot.margin = margin(5.5, 0, 5.5, 5.5)
    )
  
  # ----------------------------------------------------------------------
  # Right side: labels only
  # No box, no title, no axes.
  # angle = 90, 1st / 2nd / 3rd columns.
  # ----------------------------------------------------------------------
  
  df_legend <- dplyr::bind_rows(
    tibble::tibble(
      sens_type = "UC",
      setting = uc_levels,
      ord = c(1, 2, 3),
      label = unname(uc_labels_plain[uc_levels])
    ),
    tibble::tibble(
      sens_type = "kappa",
      setting = kappa_levels,
      ord = c(1, 2, 3),
      label = unname(kappa_labels_plain[kappa_levels])
    ),
    tibble::tibble(
      sens_type = "eta",
      setting = eta_levels,
      ord = c(1, 2, 3),
      label = unname(eta_labels_plain[eta_levels])
    )
  ) %>%
    dplyr::mutate(
      sens_type = factor(sens_type, levels = sens_type_levels),
      setting   = factor(setting, levels = all_setting_levels),
      
      x = dplyr::case_when(
        ord == 1 ~ 0.22,
        ord == 2 ~ 0.50,
        ord == 3 ~ 0.78,
        TRUE ~ 0.50
      ),
      y = 0.50
    )
  
  p_right <- ggplot(df_legend) +
    geom_text(
      aes(
        x = x,
        y = y,
        label = label,
        color = setting
      ),
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      size = 2.15
    ) +
    facet_grid(
      sens_type ~ .,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = all_setting_cols,
      drop = FALSE
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      clip = "off"
    ) +
    guides(color = "none") +
    theme_void(base_size = 12) +
    theme(
      strip.text = element_blank(),
      panel.spacing.y = unit(0.35, "lines"),
      plot.margin = margin(5.5, 0, 5.5, -2)
    )
  
  # ----------------------------------------------------------------------
  # Combine
  # ----------------------------------------------------------------------
  
  p_left + p_right +
    patchwork::plot_layout(
      ncol = 2,
      widths = c(5.20, 0.32)
    )
}

# ------------------------------------------------------------------------
# Main-text Figure 3: MCI PSQC sensitivity analysis ----------------
# ------------------------------------------------------------------------

df_compact_MCI <- build_compact_sens_df_for_group(
  group_label = "MCI",
  main_obj = SURVIVAL_EDPMM_POST_condX1_main,
  uc_obj = SURVIVAL_EDPMM_POST_condX1_uc,
  kappa01_obj = SURVIVAL_EDPMM_POST_condX1_kappa_unif01,
  kappam11_obj = SURVIVAL_EDPMM_POST_condX1_kappa_unifm11,
  eta_e11_e09_obj = SURVIVAL_EDPMM_POST_condX1_eta_e11_e09,
  eta_e09_e11_obj = SURVIVAL_EDPMM_POST_condX1_eta_e09_e11
)

EDPqrlPOST_PSQCrho_MCI_all_sens <-
  plot_compact_sens_figure(
    df_plot = df_compact_MCI,
    group_label = "MCI"
  )

print(EDPqrlPOST_PSQCrho_MCI_all_sens)

save_figure_pdf(
  EDPqrlPOST_PSQCrho_MCI_all_sens,
  "EDPqrlPOST_PSQC_MCI",
  width = 12.0,
  height = 7.2
)

# ------------------------------------------------------------------------
# Supplementary Figure S1: Dementia-free survival by amyloid status ----------------------------
# ------------------------------------------------------------------------

df_S_overall <- make_S_df(SURVIVAL_EDPMM_POST_ftn_results_main)

df_S_cond <- dplyr::bind_rows(
  make_cond_surv_df(
    SURVIVAL_EDPMM_POST_ftn_condX0_results_main$S_z0_result,
    "CN",
    "Non-Elevated"
  ),
  make_cond_surv_df(
    SURVIVAL_EDPMM_POST_ftn_condX0_results_main$S_z1_result,
    "CN",
    "Elevated"
  ),
  make_cond_surv_df(
    SURVIVAL_EDPMM_POST_ftn_condX1_results_main$S_z0_result,
    "MCI",
    "Non-Elevated"
  ),
  make_cond_surv_df(
    SURVIVAL_EDPMM_POST_ftn_condX1_results_main$S_z1_result,
    "MCI",
    "Elevated"
  )
) %>%
  dplyr::mutate(
    subgroup = factor(subgroup, levels = c("CN", "MCI")),
    amyloid  = factor(amyloid, levels = c("Non-Elevated", "Elevated"))
  ) %>%
  dplyr::arrange(subgroup, amyloid, t)

km_fit_overall <- survival::survfit(
  survival::Surv(Tobs, Dobs) ~ Zobs,
  data = dfADNI
)

km_sum_overall <- summary(km_fit_overall)

df_KM_overall <- data.frame(
  t      = km_sum_overall$time,
  surv   = km_sum_overall$surv,
  strata = km_sum_overall$strata
) %>%
  dplyr::mutate(
    amyloid = factor(
      ifelse(grepl("=0", strata), "Non-Elevated", "Elevated"),
      levels = c("Non-Elevated", "Elevated")
    )
  ) %>%
  dplyr::arrange(amyloid, t)

df_km_cond <- data.frame(
  Tobs = dfADNI$Tobs,
  Dobs = dfADNI$Dobs,
  Baseline_DX = dfADNI$Xobs[, "Baseline_DX"],
  Zobs = dfADNI$Zobs
) %>%
  dplyr::mutate(
    subgroup = factor(Baseline_DX, levels = c(0, 1), labels = c("CN", "MCI")),
    amyloid  = factor(Zobs, levels = c(0, 1), labels = c("Non-Elevated", "Elevated"))
  ) %>%
  dplyr::filter(
    !is.na(subgroup),
    !is.na(amyloid),
    !is.na(Tobs),
    !is.na(Dobs)
  )

km_fit_cond <- survival::survfit(
  survival::Surv(Tobs, Dobs) ~ subgroup + amyloid,
  data = df_km_cond
)

km_sum_cond <- summary(km_fit_cond)

df_KM_cond <- data.frame(
  t      = km_sum_cond$time,
  surv   = km_sum_cond$surv,
  strata = km_sum_cond$strata
) %>%
  dplyr::mutate(
    strata = as.character(strata),
    subgroup = ifelse(grepl("subgroup=CN", strata), "CN", "MCI"),
    amyloid  = ifelse(grepl("amyloid=Non-Elevated", strata), "Non-Elevated", "Elevated"),
    subgroup = factor(subgroup, levels = c("CN", "MCI")),
    amyloid  = factor(amyloid, levels = c("Non-Elevated", "Elevated"))
  ) %>%
  dplyr::arrange(subgroup, amyloid, t)

EDPqrlPOST_survival_overall <-
  plot_survival_estimator(
    df_post = df_S_overall,
    df_km   = df_KM_overall,
    panel_subtitle = "Overall"
  )

EDPqrlPOST_survival_MCI <-
  plot_survival_estimator(
    df_post = df_S_cond %>% dplyr::filter(subgroup == "MCI"),
    df_km   = df_KM_cond %>% dplyr::filter(subgroup == "MCI"),
    panel_subtitle = "MCI"
  )

EDPqrlPOST_survival_CN <-
  plot_survival_estimator(
    df_post = df_S_cond %>% dplyr::filter(subgroup == "CN"),
    df_km   = df_KM_cond %>% dplyr::filter(subgroup == "CN"),
    panel_subtitle = "CN"
  )

EDPqrlPOST_survival_all3 <-
  (
    EDPqrlPOST_survival_overall +
      EDPqrlPOST_survival_MCI +
      EDPqrlPOST_survival_CN
  ) +
  patchwork::plot_layout(
    ncol = 3,
    guides = "collect",
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(EDPqrlPOST_survival_all3)

save_figure_pdf(
  EDPqrlPOST_survival_all3,
  "EDPqrlPOST_survival_all3_001",
  width = 12,
  height = 4.5
)

# ------------------------------------------------------------------------
# Supplementary Figure S2: Overall PSQC sensitivity analysis ------------
# ------------------------------------------------------------------------

df_compact_Overall <- build_compact_sens_df_for_group(
  group_label = "Overall",
  main_obj = EDPMM_results_POST_main,
  uc_obj = EDPMM_results_POST_uc,
  kappa01_obj = EDPMM_results_POST_kappa_unif01,
  kappam11_obj = EDPMM_results_POST_kappa_unifm11,
  eta_e11_e09_obj = EDPMM_results_POST_eta_e11_e09,
  eta_e09_e11_obj = EDPMM_results_POST_eta_e09_e11
)

EDPqrlPOST_PSQCrho_Overall_all_sens <-
  plot_compact_sens_figure(
    df_plot = df_compact_Overall,
    group_label = "Overall"
  )

print(EDPqrlPOST_PSQCrho_Overall_all_sens)

save_figure_pdf(
  EDPqrlPOST_PSQCrho_Overall_all_sens,
  "EDPqrlPOST_PSQC_Overall",
  width = 12.0,
  height = 7.2
)

# ------------------------------------------------------------------------
# Supplementary Figure S3: CN PSQC sensitivity analysis -----------------
# ------------------------------------------------------------------------

df_compact_CN <- build_compact_sens_df_for_group(
  group_label = "CN",
  main_obj = SURVIVAL_EDPMM_POST_condX0_main,
  uc_obj = SURVIVAL_EDPMM_POST_condX0_uc,
  kappa01_obj = SURVIVAL_EDPMM_POST_condX0_kappa_unif01,
  kappam11_obj = SURVIVAL_EDPMM_POST_condX0_kappa_unifm11,
  eta_e11_e09_obj = SURVIVAL_EDPMM_POST_condX0_eta_e11_e09,
  eta_e09_e11_obj = SURVIVAL_EDPMM_POST_condX0_eta_e09_e11
)

EDPqrlPOST_PSQCrho_CN_all_sens <-
  plot_compact_sens_figure(
    df_plot = df_compact_CN,
    group_label = "CN"
  )

print(EDPqrlPOST_PSQCrho_CN_all_sens)

save_figure_pdf(
  EDPqrlPOST_PSQCrho_CN_all_sens,
  "EDPqrlPOST_PSQC_CN",
  width = 12.0,
  height = 7.2
)

# ------------------------------------------------------------------------
# Additional diagnostic figure (not reported in the manuscript):
# Kaplan--Meier plot annotated at survival probability 0.90
# ------------------------------------------------------------------------
# km plot for 
# Target survival probability
target_surv <- 0.90

# =============================================================
# Function to approximate the crossing time at S(t) = target_surv
# =============================================================
get_cross_time <- function(fit, target = 0.90) {
  ss <- summary(fit)
  
  dd <- data.frame(
    time = ss$time,
    surv = ss$surv,
    strata = ss$strata,
    stringsAsFactors = FALSE
  )
  
  res <- lapply(split(dd, dd$strata), function(x) {
    x <- x[order(x$time), ]
    
    x <- rbind(
      data.frame(
        time = 0,
        surv = 1,
        strata = x$strata[1],
        stringsAsFactors = FALSE
      ),
      x
    )
    
    idx <- which(x$surv <= target)[1]
    
    if (is.na(idx)) {
      return(data.frame(strata = x$strata[1], cross_time = NA_real_))
    }
    
    if (idx == 1) {
      return(data.frame(strata = x$strata[1], cross_time = x$time[1]))
    }
    
    t1 <- x$time[idx - 1]
    s1 <- x$surv[idx - 1]
    t2 <- x$time[idx]
    s2 <- x$surv[idx]
    
    cross_time <- if (s1 == s2) {
      t2
    } else {
      t1 + (target - s1) * (t2 - t1) / (s2 - s1)
    }
    
    data.frame(
      strata = x$strata[1],
      cross_time = cross_time
    )
  })
  
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

# =============================================================
# Calculate approximate crossing times
# =============================================================
cross_df <- get_cross_time(km_fit, target = target_surv)

cross_df$group <- ifelse(
  grepl("Zobs=0", cross_df$strata),
  "Non-Elevated (Z=0)",
  "Elevated (Z=1)"
)

cross_df$label <- ifelse(
  is.na(cross_df$cross_time),
  NA,
  paste0(sprintf("%.2f", cross_df$cross_time), " years")
)

cross_plot_df <- subset(cross_df, !is.na(cross_time))
cross_red  <- subset(cross_plot_df, group == "Non-Elevated (Z=0)")
cross_blue <- subset(cross_plot_df, group == "Elevated (Z=1)")

print(cross_df)

# =============================================================
# Tuning parameters
# =============================================================
time_label_size  <- 2.5
time_label_vjust <- 3.5
blue_x_shift     <- 0
red_x_shift      <- 0

surv_label_x     <- 0
surv_label_y     <- target_surv - 0.1
surv_label_size  <- 2.5

bottom_margin    <- 28
x_title_margin   <- 18

# =============================================================
# Add horizontal line, vertical lines, points, and labels
# =============================================================
km_plot$plot <- km_plot$plot +
  geom_hline(
    yintercept = target_surv,
    color = "darkgreen",
    linetype = "longdash",
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = surv_label_x,
    y = surv_label_y,
    label = "Survival = 90%",
    color = "darkgreen",
    size = surv_label_size,
    hjust = 0
  ) +
  geom_segment(
    data = cross_red,
    aes(
      x = cross_time,
      xend = cross_time,
      y = 0,
      yend = target_surv
    ),
    inherit.aes = FALSE,
    color = "red",
    linetype = "dotted",
    linewidth = 1.0
  ) +
  geom_segment(
    data = cross_blue,
    aes(
      x = cross_time,
      xend = cross_time,
      y = 0,
      yend = target_surv
    ),
    inherit.aes = FALSE,
    color = "blue",
    linetype = "dotted",
    linewidth = 1.0
  ) +
  geom_point(
    data = cross_red,
    aes(x = cross_time, y = target_surv),
    inherit.aes = FALSE,
    color = "red",
    size = 2.5
  ) +
  geom_point(
    data = cross_blue,
    aes(x = cross_time, y = target_surv),
    inherit.aes = FALSE,
    color = "blue",
    size = 2.5
  ) +
  geom_text(
    data = cross_blue,
    aes(
      x = cross_time + blue_x_shift,
      y = -Inf,
      label = label
    ),
    inherit.aes = FALSE,
    color = "blue",
    size = time_label_size,
    hjust = 0.5,
    vjust = time_label_vjust
  ) +
  geom_text(
    data = cross_red,
    aes(
      x = cross_time + red_x_shift,
      y = -Inf,
      label = label
    ),
    inherit.aes = FALSE,
    color = "red",
    size = time_label_size,
    hjust = 0.5,
    vjust = time_label_vjust
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = sprintf("%.2f", seq(0, 1, by = 0.25))
  ) +
  coord_cartesian(
    ylim = c(0, 1.05),
    clip = "off"
  ) +
  theme(
    plot.margin = margin(t = 10, r = 20, b = bottom_margin, l = 10),
    axis.title.x = element_text(margin = margin(t = x_title_margin))
  )
# print(km_plot)
quartz(type = "pdf", file = "EDPqrl_KM_plus.pdf", width = 8, height = 6)
print(km_plot)
dev.off()

# ------------------------------------------------------------------------
# Created manuscript figure outputs --------------------------------------
# ------------------------------------------------------------------------

cat("\n\nDone. Created manuscript figure outputs:\n")
cat("Main-text Figure 1a: EDPqrl_KM.pdf\n")
cat("Main-text Figure 1b object: EDPqrl_DENSITY_T (PDF save call remains commented)\n")
cat("Main-text Figure 3:  EDPqrlPOST_PSQC_MCI.pdf\n")
cat("Supplementary Figure S1: EDPqrlPOST_survival_all3_001.pdf\n")
cat("Supplementary Figure S2: EDPqrlPOST_PSQC_Overall.pdf\n")
cat("Supplementary Figure S3: EDPqrlPOST_PSQC_CN.pdf\n")
cat("Additional diagnostic: EDPqrl_KM_plus.pdf\n")

if (save_figures) {
  cat("\nPDF files saved in:\n")
  cat(normalizePath(fig_dir), "\n")
  
  cat("\nPDF files:\n")
  print(
    list.files(
      fig_dir,
      pattern = "\\.pdf$",
      full.names = TRUE
    )
  )
} else {
  cat("\nThe save_figures switch applies to Figures 3 and S1--S3.\n")
  cat("Figure 1a and the diagnostic plot use the existing direct quartz calls.\n")
  cat("The Figure 1b PDF save call remains commented, as in the current code.\n")
}

# ------------------------------------------------------------------------
# Final checkpoint -------------------------------------------------------
# ------------------------------------------------------------------------
# save.image(file = checkpoint_final)
load(file = checkpoint_final)
message("Saved final analysis checkpoint: ", checkpoint_final)