# -----------------------------------------------------------------------------
# Load R packages
# -----------------------------------------------------------------------------
library(Rcpp)
library(RcppArmadillo)

# -----------------------------------------------------------------------------
# Load R code
# -----------------------------------------------------------------------------
source("EDPqrl_r.R")
sourceCpp("EDPqrl_cpp.cpp")

# -----------------------------------------------------------------------------
# Run ID and seed
# -----------------------------------------------------------------------------
run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(run_ID)) run_ID = 1L
set.seed(run_ID)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# Main results: 95% credible intervals for all censoring scenarios
level = 0.05
SCNs = c(1, 2, 3, 4)
Ns = c(500, 1500)

# # Supplementary results: 99% credible intervals under 60% and 80% censoring
# level = 0.01
# SCNs = c(3, 4)
# Ns = c(500, 1500)

# # Supplementary results: N = 5000 under 60% and 80% censoring
# level = 0.05
# SCNs = c(3, 4)
# Ns = c(5000)

M = 1e3

nu_vec = c(0, 1, 2, 3)

esttype = "mean"
saveall = FALSE

dir.create("Results", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
for (Scn in SCNs) {
  for (N in Ns) {
    
    # ---------------------------------------------------------------------------
    # MCMC settings
    # ---------------------------------------------------------------------------
    gibbs_thin = 2e1
    gibbs_iter = 2e4
    gibbs_burnin = 2e4
    
    n_MCMC = floor(gibbs_iter / gibbs_thin)
    
    # ---------------------------------------------------------------------------
    # Quantile levels
    # ---------------------------------------------------------------------------
    if (Scn %in% c(1, 2)) {
      rho_vec = c(0.30, 0.60)
    } else {
      rho_vec = c(0.10, 0.20)
    }
    
    # ---------------------------------------------------------------------------
    # Data generation
    # ---------------------------------------------------------------------------
    temp_data = generate_data(Scn, N)
    
    # ---------------------------------------------------------------------------
    # Shared kappa draws for DPMM and EDPMM
    # ---------------------------------------------------------------------------
    kappa_vec = NULL
    # kappa_vec = runif(n_MCMC, 0, 1)
    
    # ---------------------------------------------------------------------------
    # Model fitting and post-processing
    # ---------------------------------------------------------------------------
    DPMM_results = SURVIVAL_DPMM(
      temp_data,
      M = M,
      nu = nu_vec,
      rho = rho_vec,
      kappa = kappa_vec,
      gibbs_iter = gibbs_iter,
      gibbs_burnin = gibbs_burnin,
      gibbs_thin = gibbs_thin,
      esttype = esttype,
      saveall = saveall,
      level = level
    )
    
    EDPMM_results = SURVIVAL_EDPMM(
      temp_data,
      M = M,
      nu = nu_vec,
      rho = rho_vec,
      kappa = kappa_vec,
      gibbs_iter = gibbs_iter,
      gibbs_burnin = gibbs_burnin,
      gibbs_thin = gibbs_thin,
      esttype = esttype,
      saveall = saveall,
      level = level
    )
    
    # ---------------------------------------------------------------------------
    # Finite-sample true PSQC in generated data
    # ---------------------------------------------------------------------------
    results_list = list()
    counter = 1L
    
    for (ii_nu in nu_vec) {
      for (ii_rho in rho_vec) {
        
        idx_ps = (temp_data$Y0 > ii_nu) & (temp_data$Y1 > ii_nu)
        
        Y0_sub = temp_data$Y0[idx_ps] - ii_nu
        Y1_sub = temp_data$Y1[idx_ps] - ii_nu
        
        Q0 = if (length(Y0_sub) > 0) {
          unname(quantile(Y0_sub, probs = ii_rho, na.rm = TRUE))
        } else {
          NA_real_
        }
        
        Q1 = if (length(Y1_sub) > 0) {
          unname(quantile(Y1_sub, probs = ii_rho, na.rm = TRUE))
        } else {
          NA_real_
        }
        
        results_list[[counter]] = data.frame(
          nu = ii_nu,
          rho = ii_rho,
          diff = Q1 - Q0,
          res_Y1 = Q1,
          res_Y0 = Q0,
          ps_prob = mean(idx_ps)
        )
        
        counter = counter + 1L
      }
    }
    
    temp_data_results = do.call(rbind, results_list)
    
    # ---------------------------------------------------------------------------
    # Extract posterior summary matrices
    # ---------------------------------------------------------------------------
    edpmm_matrix = EDPMM_results$E_rho_nu_diff_result
    dpmm_matrix = DPMM_results$E_rho_nu_diff_result
    
    # ---------------------------------------------------------------------------
    # Flatten results
    # ---------------------------------------------------------------------------
    truth_target_cols = c("nu", "rho", "diff", "res_Y1", "res_Y0", "ps_prob")
    
    flat_truth = c(run_ID, c(t(temp_data_results[, truth_target_cols])))
    flat_edpmm = c(run_ID, c(t(edpmm_matrix)))
    flat_dpmm = c(run_ID, c(t(dpmm_matrix)))
    
    # ---------------------------------------------------------------------------
    # Output file names
    #
    # Main 95% CrI results and N = 5000 results have no level suffix.
    # Only the supplementary 99% CrI results use "_level1".
    # ---------------------------------------------------------------------------
    file_suffix = if (isTRUE(all.equal(level, 0.01))) {
      "_level1"
    } else {
      ""
    }
    
    file_true = paste0(
      "Results/Truth_Scn", Scn,
      "_N", N,
      file_suffix,
      ".txt"
    )
    
    file_edpmm = paste0(
      "Results/EDPMM_Scn", Scn,
      "_N", N,
      file_suffix,
      ".txt"
    )
    
    file_dpmm = paste0(
      "Results/DPMM_Scn", Scn,
      "_N", N,
      file_suffix,
      ".txt"
    )
    
    # ---------------------------------------------------------------------------
    # Write results
    # ---------------------------------------------------------------------------
    if (run_ID == 1) {
      
      header_true = c(
        "run_ID",
        as.vector(outer(truth_target_cols, seq_len(nrow(temp_data_results)), paste, sep = "_row"))
      )
      
      post_cols = colnames(edpmm_matrix)
      header_post = c(
        "run_ID",
        as.vector(outer(post_cols, seq_len(nrow(edpmm_matrix)), paste, sep = "_row"))
      )
      
      write.table(
        t(flat_truth),
        file = file_true,
        sep = "\t",
        row.names = FALSE,
        col.names = header_true,
        append = FALSE
      )
      
      write.table(
        t(flat_edpmm),
        file = file_edpmm,
        sep = "\t",
        row.names = FALSE,
        col.names = header_post,
        append = FALSE
      )
      
      write.table(
        t(flat_dpmm),
        file = file_dpmm,
        sep = "\t",
        row.names = FALSE,
        col.names = header_post,
        append = FALSE
      )
      
    } else {
      
      write.table(
        t(flat_truth),
        file = file_true,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE
      )
      
      write.table(
        t(flat_edpmm),
        file = file_edpmm,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE
      )
      
      write.table(
        t(flat_dpmm),
        file = file_dpmm,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE
      )
    }
  }
}