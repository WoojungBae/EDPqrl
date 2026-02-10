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
# This code must be in your current directory or you can change the path.
setwd("/Users/WooJung/Documents/Rproject/EDPqrl/source")

# Load R code
source("EDPqrl_r.R")

# Load cpp code
sourceCpp("EDPqrl_cpp.cpp")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# setwd("/Users/WooJung/Documents/Rproject/EDPqrl")
# setwd("C:/Users/Woojung.Bae/OneDrive - FDA/Documents/Rporject/ICE")

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
run_ID = 1
set.seed(run_ID)

# 
Scn = 1
N = 500
# N = 1000
{
  # esttype and saveall
  esttype = "mean" #
  # esttype = "median" #
  # saveall = TRUE
  saveall = FALSE
  
  # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
  # gibbs_thin = 1e1
  # if (N > 1000) {
  #   # for big N (e.g. N = 1500), can use less samples (check posterior)
  #   gibbs_iter = 1e3
  #   gibbs_burnin = 2e3
  # } else if (N > 500) {
  #   # for mid N (e.g. N = 1000), need more samples
  #   gibbs_iter = 2e3
  #   gibbs_burnin = 2e3
  # } else {
  #   # for small N (e.g. N = 500), need more samples
  #   gibbs_iter = 2e3
  #   gibbs_burnin = 3e3
  # }
  # 
  # gibbs_thin = 1e0
  # gibbs_iter = 2e2
  # gibbs_burnin = 2e2
  # M = 5e2
  # 
  # gibbs_thin = 1e1
  # gibbs_iter = 2e3
  # gibbs_burnin = 2e3
  # M = 5e2
  # 
  # gibbs_thin = 1e1
  # gibbs_iter = 5e3
  # gibbs_burnin = 5e3
  # M = 1e3
  # 
  # gibbs_thin = 1e2
  # gibbs_iter = 2e4
  # gibbs_burnin = 2e4
  # M = 2e4
  
  # num_MC = 1e4
  # num_MC_prior = 1e5
  
  # Define number of observations for each dataset
  
  # ------------------------------------------------------------------------------
  # Data Generation --------------------------------------------------------------
  # ------------------------------------------------------------------------------
  temp_data = generate_data(Scn, N)
  
  # # ------------------------------------------------------------------------------
  # # ------------------------------------------------------------------------------
  # # ------------------------------------------------------------------------------
  # {
  #   temp_data = generate_data(Scn, N)
  #   (temp_data$censoring)
  #   
  #   plot(density(temp_data$logY),
  #        xlim = c(-12, 15),
  #        ylim = c(0, 0.35), 
  #        main = "Truth (Y) vs Observed (T)",
  #        xlab = "Log Time",
  #        col = "black",
  #        lty = 1,
  #        lwd = 1)
  #   # lines(density(temp_data$logY0), col = "blue", lty = 1, lwd = 2)
  #   # lines(density(temp_data$logY1), col = "red", lty = 1, lwd = 3)
  #   
  #   lines(density(temp_data$logT), col = "green", lty = 3, lwd = 1)
  #   lines(density(temp_data$logT[temp_data$Zobs==0]), col = "darkgreen", lty = 3, lwd = 2)
  #   lines(density(temp_data$logT[temp_data$Zobs==1]), col = "lightgreen", lty = 3, lwd = 3)
  #   legend("topright",
  #          legend = c("logY (Truth Total)", 
  #                     "logY0 (Truth Control)", 
  #                     "logY1 (Truth Treated)", 
  #                     "logT (Observed Total)", 
  #                     "logT (Observed Z=0)", 
  #                     "logT (Observed Z=1)"),
  #          col = c("black", "blue", "red", "green", "darkgreen", "lightgreen"),
  #          lty = c(1, 1, 1, 3, 3, 3),
  #          lwd = c(1, 2, 3, 1, 2, 3),
  #          cex = 0.8)
  # }
}

# plot(density(temp_data$logY))
# plot(density(temp_data$logT))

{
  # # ------------------------------------------------------------------------------
  # # ------------------------------------------------------------------------------
  # # ------------------------------------------------------------------------------
  # txt.title = paste0("Results/estimand_Scn",Scn,"_N",N,"_result.txt")
  # if (run_ID == 1) {
  #   CE_names = c("copula.DPM.Y1", "copula.DPM.Y0", "copula.DPM.E10",
  #                "copula.single.Y1", "copula.single.Y0", "copula.single.E10",
  #                "copula.separate.Y1", "copula.separate.Y0", "copula.separate.E10",
  #                "copula.single.monotone.Y1", "copula.single.monotone.Y0", "copula.single.monotone.E10",
  #                "copula.separate.monotone.Y1", "copula.separate.monotone.Y0", "copula.separate.monotone.E10",
  #                "naive.Y1", "naive.Y0", "naive.E10",
  #                "g.single.all.Y1", "g.single.all.Y0", "g.single.all.E10",
  #                "g.single.ICE.free.Y1", "g.single.ICE.free.Y0", "g.single.ICE.free.E10",
  #                "g.separate.all.Y1", "g.separate.all.Y0", "g.separate.all.E10",
  #                "g.separate.ICE.free.Y1", "g.separate.ICE.free.Y0", "g.separate.ICE.free.E10",
  #                "MI.single.Y1", "MI.single.Y0", "MI.single.E10",
  #                "MI.separate.Y1", "MI.separate.Y0", "MI.separate.E10")
  #   if (temp_data$deterministic == F) {
  #     CE_names = c(CE_names,
  #                  "IPW.single.all.Y1", "IPW.single.all.Y0", "IPW.single.all.E10",
  #                  "IPW.single.ICE.free.Y1", "IPW.single.ICE.free.Y0", "IPW.single.ICE.free.E10",
  #                  "IPW.separate.all.Y1", "IPW.separate.all.Y0", "IPW.separate.all.E10",
  #                  "IPW.separate.ICE.free.Y1", "IPW.separate.ICE.free.Y0", "IPW.separate.ICE.free.E10")
  #   }
  #   df_col_names = c(c("run_ID", "prob.ICE", "true.Y1", "true.Y0", "true.E10"), CE_names)
  #   df = data.frame(matrix(ncol = length(df_col_names), nrow = 0))
  #   colnames(df) = df_col_names
  #   write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
  # }
  # 
  # ------------------------------------------------------------------------------
  # Estimand ---------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # SURVIVAL_DPMM_MCMC_results = SURVIVAL_DPMM_MCMC(temp_data, gibbs_iter, gibbs_burnin, gibbs_thin)
  # warnings()
  gibbs_thin = 1e1
  gibbs_iter = 2e4
  gibbs_burnin = 2e4
  
  M = 2e3
  
  # M = 1e3
  # gibbs_thin = 1e1
  # gibbs_iter = 2e3
  # gibbs_burnin = 2e3
  
  EDPMM_results_MCMC = SURVIVAL_EDPMM_MCMC(temp_data, gibbs_iter, gibbs_burnin, gibbs_thin)
  DPMM_results_MCMC = SURVIVAL_DPMM_MCMC(temp_data, gibbs_iter, gibbs_burnin, gibbs_thin)
  # save.image(paste0("~/Documents/Rproject/EDPqrl/source/MCMCrestultsN",N,".RData"))
  # allinfo = c(run_ID,
  #             temp_data$prob.ICE,
  #             c(temp_data$true.Y1,
  #               temp_data$true.Y0,
  #               temp_data$true.E10),
  #             c(BNPcopula_results$Yz1.cop_result_mc[1],
  #               BNPcopula_results$Yz0.cop_result_mc[1],
  #               BNPcopula_results$CE.cop_result_mc[1]),
  #             unlist(copula_results),
  #             unlist(reference_results))
  # 
  # allinfo = data.frame(t(allinfo))
  
  # nu  rho True_Diff True_Diff.1 True_Diff.2 True_Diff.3
  # 1  0 0.25    -0.039      -0.078      -0.117      -0.156
  # 2  0 0.50     0.781       1.562       2.344       3.125
  # 3  1 0.25     0.958       1.915       2.874       3.833
  # 4  1 0.50     4.735       9.463      14.195      18.929
  # 5  2 0.25     1.325       2.651       3.976       5.302
  # 6  2 0.50     5.975      11.942      17.910      23.882
  # 7  3 0.25     1.644       3.286       4.928       6.571
  # 8  3 0.50     6.861      13.708      20.551      27.404
}

{
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # EDPMM
  # First, check the trace plot of the number of outcome clusters ($K_y$). It
  # should stabilize around 3 (or slightly higher, like 4 or 5, if it split the
  # non-normal distributions). # Plot the trace of K (Number of Outcome Clusters)
  K_trace = unlist(EDPMM_results_MCMC$MCMCposteriors$KyLists)
  plot(K_trace, type="s", main="Trace of Ky (Should be >= 3)",
       xlab="Saved Iteration", ylab="Number of Clusters")
  abline(h=3, col="red", lwd=2, lty=2) # The Truth

  # ----------------------------------------------------------------------------
  # Since the "label switching" problem makes it hard to track "Cluster 1" vs
  # "Cluster 2" over time, the best way to verify is to pool all the intercept
  # parameters from the posterior samples and plot them.
  # You should see three distinct spikes at -3.0, -0.5, and 2.0.
  # Success Indicator: You should see the blue curve matching the red lines perfectly.
  # Loop through stored samples
  # Extract intercepts from all saved iterations
  beta_samples = EDPMM_results_MCMC$MCMCposteriors$betaLists
  all_intercepts = c()
  for(i in 1:length(beta_samples)) {
    # Get the beta matrix for this iteration
    betas = beta_samples[[i]]
    # Row 1 is usually the intercept (assuming matX column 1 is intercept)
    intercepts = betas[1, ]
    all_intercepts = c(all_intercepts, intercepts)
  }
  # Plot the density of all found intercepts
  plot(density(all_intercepts), main="Posterior Density of Intercepts",
       xlim = c(-7, 7),xlab="Intercept Value", lwd=2, col="blue")
  abline(v=c(-3.0, -0.5, 2.0), col="red", lty=2, lwd=2)
  legend("topright", legend=c("Posterior", "Truth"), col=c("blue", "red"), lty=c(1,2))

  # ----------------------------------------------------------------------------
  # Check the estimated variances ($\sigma^2$). With your new prior (0.1), these
  # should be very small (close to 0), not stuck at 1.0.
  # Success Indicator: The histogram should be clustered to the left (negative
  # values on log scale), around $-3$ or $-4$. If the histogram is centered
  # around $0$ (which is $\sigma^2=1$), the model failed.
  # Extract variances
  sig2_samples = EDPMM_results_MCMC$MCMCposteriors$sig2Lists
  all_sig2 = unlist(sig2_samples)
  # Summary
  print(summary(all_sig2))
  # Plot (Log scale to see small values)
  hist(log10(all_sig2), main="Log10 Variance Posterior",
       xlab="Log10(Sigma^2)", breaks=30, col="lightgreen")
  # Truth is approx 0.0001 -> Log10 is -4
  abline(v=-4, col="red", lwd=2)
  
  # # ----------------------------------------------------------------------------
  # # ----------------------------------------------------------------------------
  # # ----------------------------------------------------------------------------
  # # DPMM
  # # First, check the trace plot of the number of outcome clusters ($K_y$). It
  # # should stabilize around 3 (or slightly higher, like 4 or 5, if it split the
  # # non-normal distributions). # Plot the trace of K (Number of Outcome Clusters)
  # K_trace = unlist(DPMM_results_MCMC$MCMCposteriors$KLists)
  # plot(K_trace, type="s", main="Trace of Ky (Should be >= 3)",
  #      xlab="Saved Iteration", ylab="Number of Clusters")
  # abline(h=3, col="red", lwd=2, lty=2) # The Truth
  # 
  # # ----------------------------------------------------------------------------
  # # Since the "label switching" problem makes it hard to track "Cluster 1" vs
  # # "Cluster 2" over time, the best way to verify is to pool all the intercept
  # # parameters from the posterior samples and plot them.
  # # You should see three distinct spikes at -3.0, -0.5, and 2.0.
  # # Success Indicator: You should see the blue curve matching the red lines perfectly.
  # # Loop through stored samples
  # # Extract intercepts from all saved iterations
  # beta_samples = DPMM_results_MCMC$MCMCposteriors$betaLists
  # all_intercepts = c()
  # for(i in 1:length(beta_samples)) {
  #   # Get the beta matrix for this iteration
  #   betas = beta_samples[[i]]
  #   # Row 1 is usually the intercept (assuming matX column 1 is intercept)
  #   intercepts = betas[1, ]
  #   all_intercepts = c(all_intercepts, intercepts)
  # }
  # # Plot the density of all found intercepts
  # plot(density(all_intercepts), main="Posterior Density of Intercepts",
  #      xlim = c(-7, 7),xlab="Intercept Value", lwd=2, col="blue")
  # abline(v=c(-3.0, -0.5, 2.0), col="red", lty=2, lwd=2)
  # legend("topright", legend=c("Posterior", "Truth"), col=c("blue", "red"), lty=c(1,2))
  # 
  # # ----------------------------------------------------------------------------
  # # Check the estimated variances ($\sigma^2$). With your new prior (0.1), these
  # # should be very small (close to 0), not stuck at 1.0.
  # # Success Indicator: The histogram should be clustered to the left (negative
  # # values on log scale), around $-3$ or $-4$. If the histogram is centered
  # # around $0$ (which is $\sigma^2=1$), the model failed.
  # # Extract variances
  # sig2_samples = DPMM_results_MCMC$MCMCposteriors$sig2Lists
  # all_sig2 = unlist(sig2_samples)
  # # Summary
  # print(summary(all_sig2))
  # # Plot (Log scale to see small values)
  # hist(log10(all_sig2), main="Log10 Variance Posterior",
  #      xlab="Log10(Sigma^2)", breaks=30, col="lightgreen", xlim = c(-1,1))
  # # Truth is approx 0.0001 -> Log10 is -4
  # abline(v=-4, col="red", lwd=2)
}

{
  # setwd("/Users/WooJung/Documents/Rproject/EDPqrl/source")
  # load(paste0("~/Documents/Rproject/EDPqrl/source/MCMCrestultsN",N,".RData"))
  source("EDPqrl_r.R")
  sourceCpp("EDPqrl_cpp.cpp")
  nu = c(0, 1, 2, 3)
  rho = c(0.25, 0.50)
  
  esttype = "mean"
  # esttype = "median"
  EDPMM_results_POST = SURVIVAL_EDPMM_POST(EDPMM_results_MCMC, nu, rho, M, esttype)
  EDPMM_results_POST$E_rho_nu_diff_result
  
  DPMM_results_POST = SURVIVAL_DPMM_POST(DPMM_results_MCMC, nu, rho, M, esttype)
  DPMM_results_POST$E_rho_nu_diff_result
  
  nu_vec  = c(0, 1, 2, 3)
  rho_vec = c(0.25, 0.50)
  results_list = list()
  counter = 1
  for(ii_nu in nu_vec) {
    for(ii_rho in rho_vec) {
      Y1_sub = temp_data$Y1[temp_data$Y1 > ii_nu]
      Y0_sub = temp_data$Y0[temp_data$Y0 > ii_nu]
      Q1 = if(length(Y1_sub) > 0) unname(quantile(Y1_sub, probs = ii_rho, na.rm = TRUE)) else NA
      Q0 = if(length(Y0_sub) > 0) unname(quantile(Y0_sub, probs = ii_rho, na.rm = TRUE)) else NA
      val_diff = Q1 - Q0
      val_Y1 = Q1 - ii_nu
      val_Y0 = Q0 - ii_nu
      results_list[[counter]] = data.frame(
        nu = ii_nu,
        rho = ii_rho,
        diff   = val_diff,
        res_Y1 = val_Y1,
        res_Y0 = val_Y0
      )
      counter = counter + 1
    }
  }
  temp_data_results = do.call(rbind, results_list)
  
  print("temp_data_results")
  print(temp_data_results)
  print("DPMM_results_POST")
  print(DPMM_results_POST$E_rho_nu_diff_result)
  # print(DPMM_results_POST$Y_rho_nu_z1x_result)
  # print(DPMM_results_POST$Y_rho_nu_z0x_result)
  print("EDPMM_results_POST")
  print(EDPMM_results_POST$E_rho_nu_diff_result)
  # print(EDPMM_results_POST$Y_rho_nu_z1x_result)
  # print(EDPMM_results_POST$Y_rho_nu_z0x_result)
}

# {
#   esttype = "mean"
#   # esttype = "median"
#   level = 0.05
#   quantile_alpha = c(level/2,1-level/2)
# 
#   n_nu = length(nu)
#   n_rho = length(rho)
# 
#   z0_tables_list   = list()
#   z1_tables_list   = list()
#   diff_tables_list = list()
#   counter = 1
#   for(ii_nu in 1:n_nu) {
#     for(ii_rho in 1:n_rho) {
#       val_nu  = nu_vec[ii_nu]
#       val_rho = rho_vec[ii_rho]
# 
#       draws_z0   = DPMM_results_POST$Y_rho_nu_z0x_store[, ii_nu, ii_rho]
#       draws_z1   = DPMM_results_POST$Y_rho_nu_z1x_store[, ii_nu, ii_rho]
#       draws_diff = DPMM_results_POST$E_rho_nu_diff_store[, ii_nu, ii_rho]
# 
#       summ_z0   = POSTsummary(draws_z0, esttype, quantile_alpha)
#       summ_z1   = POSTsummary(draws_z1, esttype, quantile_alpha)
#       summ_diff = POSTsummary(draws_diff, esttype, quantile_alpha)
# 
#       summ_z0   = cbind(nu = val_nu, rho = val_rho, summ_z0)
#       summ_z1   = cbind(nu = val_nu, rho = val_rho, summ_z1)
#       summ_diff = cbind(nu = val_nu, rho = val_rho, summ_diff)
# 
#       z0_tables_list[[counter]]   = summ_z0
#       z1_tables_list[[counter]]   = summ_z1
#       diff_tables_list[[counter]] = summ_diff
# 
#       counter = counter + 1
#     }
#   }
# 
#   Y_rho_nu_z0x_result  = do.call(rbind, z0_tables_list)
#   Y_rho_nu_z1x_result  = do.call(rbind, z1_tables_list)
#   E_rho_nu_diff_result = do.call(rbind, diff_tables_list)
# }
