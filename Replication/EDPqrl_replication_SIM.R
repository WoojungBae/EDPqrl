# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Load R packages
library(Rcpp)
library(RcppArmadillo)

# library(MASS)
# library(mice)
# library(locfit)
# library(latex2exp)
# library(rlang)
# library(devtools)
# library(ggplot2)
# library(grid)
# library(tidyverse)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Load R code
source("EDPqrl_r.R")
# Load cpp code
sourceCpp("EDPqrl_cpp.cpp")

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(run_ID)

# Scenario
SCNs = c(1, 2, 3, 4)
# Sample size
Ns = c(500, 1500)

# # Supp Resutls for N = 5000 for scenario 4 (80% censoring).
# # Scenario
# SCNs = c(4)
# # Sample size
# Ns = c(5000)

for (Scn in SCNs) {
  for (N in Ns) {
    # ------------------------------------------------------------------------------
    # Data Generation --------------------------------------------------------------
    # ------------------------------------------------------------------------------
    temp_data = generate_data(Scn, N)
    
    # esttype and saveall
    esttype = "mean" # "median"
    saveall = FALSE
    
    # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
    # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
    if (N > 2000) {
      # for big N (e.g. N = 5000), can use less samples (check posterior)
      gibbs_thin = 5
      gibbs_iter = 5e3
      gibbs_burnin = 5e3
    } else {
      # for small N (e.g. N = 500), need more samples
      gibbs_thin = 2e1
      gibbs_iter = 2e4
      gibbs_burnin = 2e4
    }
    
    M = 1e3
    # M = 1e4
    # num_MC = 1e4
    # num_MC_prior = 1e5
    
    nu_vec  = c(0, 1, 2, 3)
    rho_vec = c(0.30, 0.60)
    esttype = "mean" 
    
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    DPMM_results = SURVIVAL_DPMM(temp_data, M = M, nu = nu_vec, rho = rho_vec, 
                                 gibbs_iter = gibbs_iter, gibbs_burnin = gibbs_burnin, gibbs_thin = gibbs_thin, 
                                 esttype = esttype, saveall = saveall)
    EDPMM_results = SURVIVAL_EDPMM(temp_data, M = M, nu = nu_vec, rho = rho_vec, 
                                   gibbs_iter = gibbs_iter, gibbs_burnin = gibbs_burnin, gibbs_thin = gibbs_thin, 
                                   esttype = esttype, saveall = saveall)
    
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    results_list = list()
    counter = 1
    for(ii_nu in nu_vec) {
      for(ii_rho in rho_vec) {
        Y1_sub = temp_data$Y1[temp_data$Y1 > ii_nu]
        Y0_sub = temp_data$Y0[temp_data$Y0 > ii_nu]
        
        Q1 = if(length(Y1_sub) > 0) unname(quantile(Y1_sub, probs = ii_rho, na.rm = TRUE)) else NA
        Q0 = if(length(Y0_sub) > 0) unname(quantile(Y0_sub, probs = ii_rho, na.rm = TRUE)) else NA
        
        # Store full row initially
        results_list[[counter]] = data.frame(
          nu = ii_nu,
          rho = ii_rho,
          diff   = Q1 - Q0,
          res_Y1 = Q1 - ii_nu,
          res_Y0 = Q0 - ii_nu
        )
        counter = counter + 1
      }
    }
    temp_data_results = do.call(rbind, results_list)
    
    # ------------------------------------------------------------------------------
    edpmm_matrix = EDPMM_results$E_rho_nu_diff_result
    
    # ------------------------------------------------------------------------------
    dpmm_matrix  = DPMM_results$E_rho_nu_diff_result
    
    # ------------------------------------------------------------------------------
    truth_target_cols = c("nu", "rho", "diff")
    flat_truth = c(run_ID, c(t(temp_data_results[, truth_target_cols])))
    flat_edpmm = c(run_ID, c(t(edpmm_matrix)))
    flat_dpmm  = c(run_ID, c(t(dpmm_matrix)))
    
    file_true  = paste0("Results/Truth_Scn", Scn, "_N", N, ".txt")
    file_edpmm = paste0("Results/EDPMM_Scn", Scn, "_N", N, ".txt")
    file_dpmm  = paste0("Results/DPMM_Scn", Scn, "_N", N, ".txt")
    if (run_ID == 1) {
      # header_true = c("run_ID", as.vector(outer(truth_target_cols, 1:nrow(temp_data_results), paste, sep="_row")))
      post_cols = colnames(edpmm_matrix)
      header_post = c("run_ID", as.vector(outer(post_cols, 1:nrow(edpmm_matrix), paste, sep="_row")))
      
      # write.table(t(flat_truth), file = file_true, sep = "\t", row.names = FALSE, col.names = header_true, append = FALSE)
      write.table(t(flat_edpmm), file = file_edpmm, sep = "\t", row.names = FALSE, col.names = header_post, append = FALSE)
      write.table(t(flat_dpmm),  file = file_dpmm,  sep = "\t", row.names = FALSE, col.names = header_post, append = FALSE)
      
    } else {
      # write.table(t(flat_truth), file = file_true, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      write.table(t(flat_edpmm), file = file_edpmm, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      write.table(t(flat_dpmm),  file = file_dpmm,  sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}
