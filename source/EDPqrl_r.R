# ------------------------------------------------------------------------
# Define functions -------------------------------------------------------
# ------------------------------------------------------------------------
# Define function to update mean parameter for bivariate confounders
update_conf_bin = function(X_temp, a_bin, b_bin) {
  n_temp = length(X_temp)
  sum_X_temp = sum(X_temp)
  
  a_bin_new = a_bin + sum_X_temp
  b_bin_new = b_bin + n_temp - sum_X_temp
  
  bin_prop = c(rbeta(1, a_bin_new, b_bin_new))
  
  return(bin_prop)
  # return(list(bin_prop = bin_prop, a_bin_new = a_bin_new, b_bin_new = b_bin_new))
}

# Define function to update mean parameter for categorical confounders
update_conf_cat = function(X_temp, a_cat) {
  
  n_cat = length(a_cat)
  
  a_cat_new = count_cpp(X_temp,n_cat) + a_cat
  cat_prop  = c(rdirichlet_cpp(1,a_cat_new))
  
  return(cat_prop)
}

# Define function to update mean parameter and variance 
# parameter for continuous confounders, N-Inv-chi^2 model
update_conf_con = function(X_temp, a_var, b_var, a_mean, b_mean) {
  n_temp = length(X_temp)
  sum_X_temp = sum(X_temp)
  mean_X_temp = sum_X_temp/n_temp
  
  if (n_temp == 1) { 
    sum2_temp = 0
    weighted_mean2_temp = (b_mean * n_temp/(b_mean + n_temp)) * (X_temp - a_mean)^{2}
  } else {
    sum2_temp = (n_temp - 1) * var(X_temp)
    weighted_mean2_temp = (b_mean * n_temp/(b_mean + n_temp)) * (mean_X_temp - a_mean)^{2}
  }
  
  # Update the continuous confounders variance
  a_var_new = a_var + n_temp
  b_var_new = (a_var * b_var + sum2_temp + weighted_mean2_temp)/a_var_new
  var_prop = c(rscainvchisq_cpp(1, a_var_new, b_var_new))
  
  # Update the continuous confounders mean
  b_mean_new = b_mean + n_temp
  a_mean_new = (a_mean * b_mean + sum_X_temp)/b_mean_new
  var_new = var_prop/b_mean_new
  mean_prop = c(rnorm(1, a_mean_new, sqrt(var_new)))
  
  return(list(mean_par = mean_prop, var_par = var_prop))
}

# Define function to update outcome regression coefficients
update_reg_con = function(Y_temp, matX_temp, beta_curr,
                          a_beta, Binv_beta, aBinv_beta,
                          a_sig2, b_sig2) {
  n_temp = length(Y_temp)
  p_temp = length(beta_curr)
  
  if (n_temp == 0) {
    sig2_prior = c(rscainvchisq_cpp(1, a_sig2, b_sig2))
    sig2_prior = max(sig2_prior, 1e-5)
    beta_prior = c(rmvn_cpp(1, a_beta, sig2_prior * inv_cpp(Binv_beta)))
    return(list(beta_par = beta_prior, sig2_par = sig2_prior))
  }
  
  matX_temp = matrix(matX_temp, ncol = p_temp)
  matXtX_temp = t(matX_temp) %*% matX_temp
  matXtY_temp = t(matX_temp) %*% Y_temp
  
  # Update the outcome regression parameters (sig2, beta)
  beta_diff = beta_curr - a_beta
  SSE_prior = as.numeric(t(beta_diff) %*% Binv_beta %*% beta_diff)
  residuals = Y_temp - matX_temp %*% beta_curr
  SSE_data = sum(residuals^2)
  
  a_sig2_new = a_sig2 + n_temp + p_temp 
  b_sig2_new = (a_sig2 * b_sig2 + SSE_data + SSE_prior) / a_sig2_new
  sig2_prop = c(rscainvchisq_cpp(1, a_sig2_new, b_sig2_new))
  if(sig2_prop < 1e-5) {
    sig2_prop = 1e-5
  }
  
  Binv_beta_new = (matXtX_temp + Binv_beta) / sig2_prop
  B_beta_new = tryCatch({
    inv_cpp(Binv_beta_new)
  }, error = function(e) {
    diag(1 / diag(Binv_beta_new))
  })
  a_beta_new = B_beta_new %*% ((matXtY_temp + aBinv_beta) / sig2_prop)
  beta_prop = c(rmvn_cpp(1, a_beta_new, B_beta_new))
  
  return(list(beta_par = beta_prop, sig2_par = sig2_prop))
}

# # logit
# update_reg_bin <- function(Y_temp, matX_temp, beta_curr, 
#                            a_beta, B_beta) {
#   
#   # 1. Define a function to calculate the log of the posterior density
#   # log(posterior) = log(likelihood) + log(prior)
#   log_posterior <- function(beta) {
#     # Calculate log-likelihood in a numerically stable way
#     linear_predictor <- matX_temp %*% beta
#     log_lik <- sum(Y_temp * linear_predictor - log(1 + exp(linear_predictor)))
#     # Calculate log-prior
#     log_pri <- dmvn_cpp(beta, a_beta, B_beta, TRUE)
#     return(log_lik + log_pri)
#   }
#   
#   # 2. Propose a new candidate value for beta
#   # A common choice is a random walk proposal (Normal dist centered on the current value)
#   beta_prop <- c(rmvn_cpp(1, beta_curr, B_beta))
#   
#   # 3. Calculate the acceptance probability on the log scale
#   log_acceptance_ratio <- log_posterior(beta_prop) - log_posterior(beta_curr)
#   
#   # 4. Accept or reject the proposal
#   if (log(runif(1)) < log_acceptance_ratio) {
#     # Accept the proposal
#     return(beta_prop)
#   } else {
#     # Reject and keep the current value
#     return(beta_curr)
#   }
# }

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
generate_data = function(Scn = 1, n = 1e3) {
  # Covariates: 2 binary and 3 continous covariates
  X1 = rbinom(n, 1, 0.5)
  X2 = rbinom(n, 1, 0.4 + 0.2 * X1)
  X3 = rnorm(n, 0, 1)
  X4 = rnorm(n, - 0.2 + 0.4 * X1 - 0.3 * X3, 1)
  X5 = rnorm(n,   0.2 - 0.4 * X2 + 0.3 * X4, 0.5)
  # mean: 0.5, 0.5, 0, 0, 0
  X = cbind(X1, X2, X3, X4, X5)
  
  rowSdX = apply(X, 1, sd)
  rowMeanX = apply(X, 1, mean)
  
  # sumX2 = rowSums(X^{2})/sd(rowSums(X^{2}))
  # treatment
  z0 = 0
  z1 = 1
  Z = rbinom(n, 1, pnorm(0.2 + 0.1 * X1 + 0.2 * X3 - 0.1 * X5))
  # mean: 0.6
  
  # design matrix
  matXz0 = cbind(1, z0, X)
  matXz1 = cbind(1, z1, X)
  beta1 = c(-3.0, -0.5,  0.5,  1.0, -0.5,  0.5,  1.0)
  beta2 = c(-0.5,  0.5,  1.0, -0.5,  0.5,  1.0, -0.5)
  beta3 = c( 2.0,  1.0, -0.5,  0.5,  1.0, -0.5,  0.5)
  mean1z0 = c(matXz0 %*% beta1)
  mean2z0 = c(matXz0 %*% beta2)
  mean3z0 = c(matXz0 %*% beta3)
  mean1z1 = c(matXz1 %*% beta1)
  mean2z1 = c(matXz1 %*% beta2)
  mean3z1 = c(matXz1 %*% beta3)
  sd1 = 0.10
  sd2 = 0.05
  sd3 = 0.15
  
  betaC  = c( 0, 0.5, -0.1, 0.1, -0.2,  0.1,  -0.2)
  meanCz0 = c(matXz0 %*% betaC)
  meanCz1 = c(matXz1 %*% betaC)
  sdCz0 = 0.5
  sdCz1 = 0.5
  
  # mixtur weights
  p1 = 0.3
  p2 = 0.3
  p3 = 0.4
  cluster = sapply(1:n, function(l) sample(3, 1, F, c(p1, p2, p3)))
  
  # TRUE censoring
  # Scn1 => cstar = 2.38 # for 20% censoring
  # Scn2 => cstar = 0.57 # for 40% censoring
  # Scn3 => cstar = 2.40 # for 20% censoring
  # Scn4 => cstar = 0.56 # for 40% censoring
  if (Scn == 1) {
    cstar = 2.35
  }
  if (Scn == 2) {
    cstar = 0.76
  }
  if(Scn == 3) {
    cstar = 2.36
  }
  if(Scn == 4) {
    cstar = 0.77
  }
  logC0 = rnorm(n, cstar + meanCz0, sdCz0)
  logC1 = rnorm(n, cstar + meanCz1, sdCz1)
  
  # TRUE failure time
  # Scn1 & Scn2 => mixture of normals and uniform
  if (Scn == 1 || Scn == 2){
    logY0 = numeric(n)
    logY1 = numeric(n)
    for (i in 1:n) {
      if (cluster[i] == 1){
        logY0[i] = c(rtruncnorm_cpp(1, mean1z0[i], sd1, -5, 0))
        logY1[i] = c(rtruncnorm_cpp(1, mean1z1[i], sd1, -5, 0))
      } else if (cluster[i] == 2) {
        logY0[i] = rnorm(1, mean2z0[i], sd3)
        logY1[i] = rnorm(1, mean2z1[i], sd3)
      } else if (cluster[i] == 3) {
        logY0[i] = c(rtrunct_cpp(1, mean3z0[i], sd3, 7, 0, 5))
        logY1[i] = c(rtrunct_cpp(1, mean3z1[i], sd3, 7, 0, 5))
      }
    }
  }
  # Scn3 & Scn4 => mixture of normals
  if (Scn == 3 || Scn == 4){
    logY0 = numeric(n)
    logY1 = numeric(n)
    for (i in 1:n) {
      if (cluster[i] == 1){
        logY0[i] = rnorm(1, mean1z0[i], sd1)
        logY1[i] = rnorm(1, mean1z1[i], sd1)
      } else if (cluster[i] == 2) {
        logY0[i] = rnorm(1, mean2z0[i], sd2)
        logY1[i] = rnorm(1, mean2z1[i], sd2)
      } else if (cluster[i] == 3) {
        logY0[i] = rnorm(1, mean3z0[i], sd3)
        logY1[i] = rnorm(1, mean3z1[i], sd3)
      }
    }
  }
  
  # out
  logC = ifelse(Z == 1, logC1, logC0)
  logY = ifelse(Z == 1, logY1, logY0)
  logT = sapply(1:n, function(i) min(logY[i], logC[i]))
  D = (logY <= logC) # failure indicator
  censoring = sum(1-D)/n
  
  out = list(
    logY0 = logY0, 
    logY1 = logY1,
    logC0 = logC0,
    logC1 = logC1,
    logC  = logC,
    logY  = logY,
    Y0 = exp(logY0),
    Y1 = exp(logY1),
    C0 = exp(logC0),
    C1 = exp(logC1),
    C  = exp(logC),
    Y  = exp(logY),
    censoring = censoring,
    logT  = logT, 
    Tobs = exp(logT),
    Dobs = D,
    Zobs = Z,
    Xobs = X
  )
  
  return(out)
}

# ------------------------------------------------------------------------
POSTsummary = function(Y, esttype, quantile_alpha){
  Ytmp = Y
  Yquantile = quantile(Ytmp, c(0.25,0.75), na.rm = T)
  Yiqr = diff(Yquantile)
  Ycon = 1e2
  YiqrLOWER = which(Ytmp < (Yquantile[1] - Ycon*Yiqr))
  if (length(YiqrLOWER) > 0) {
    Ytmp = Ytmp[-YiqrLOWER]
  }
  YiqrUPPER = which(Ytmp > (Yquantile[2] + Ycon*Yiqr))
  if (length(YiqrUPPER) > 0) {
    Ytmp = Ytmp[-YiqrUPPER]
  }
  
  # Estimates: mean or median
  if (esttype == "median"){
    Y_est_mc = median(Ytmp, na.rm = T)
  } else {
    Y_est_mc = mean(Ytmp, na.rm = T)
  }
  
  # Calculate median of posterior for Y
  Y_sd_mc = sd(Ytmp, na.rm = T)
  Y_quantile_mc = quantile(Ytmp, quantile_alpha, na.rm = T)
  Y_CIlength_mc = abs(diff(Y_quantile_mc))
  Y_quantile025_mc = min(Y_quantile_mc, na.rm = T)
  Y_quantile975_mc = max(Y_quantile_mc, na.rm = T)
  
  # Save the results
  Y_result_mc = cbind(Y_est_mc, Y_sd_mc, Y_quantile025_mc, Y_quantile975_mc, Y_CIlength_mc)
  colnames(Y_result_mc) = c("estimYs","sd", "quantile025", "quantile975", "CIlength95")
  rownames(Y_result_mc) = c("")
  
  return(Y_result_mc)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Define function to fit the BNP (DPMM) model
SURVIVAL_DPMM_MCMC = function(object,
                              gibbs_iter = 1e4, gibbs_burnin = 1e4, gibbs_thin = 1e1, 
                              num_MC = NULL, num_MC_prior = NULL) {
  # =============================================================================
  # 1. SETUP & INITIALIZATION
  # =============================================================================
  # object = temp_data
  Xobs = object$Xobs
  Zobs = object$Zobs
  Dobs = object$Dobs
  Tobs = object$Tobs
  
  N = length(Tobs)
  Toriginal = Tobs
  Tobs = log(Tobs)
  
  # Identify covariate types (binary vs. continuous)
  is_binary = apply(Xobs, 2, function(col) all(col %in% c(0, 1)))
  X_bin = Xobs[, is_binary, drop = FALSE]
  X_con = Xobs[, !is_binary, drop = FALSE]
  p_X1 = ncol(X_bin)
  p_X2 = ncol(X_con)
  p_X = p_X1 + p_X2
  
  # Standardize continuous confounders (important when choosing priors)
  if (p_X1>0 && p_X2>0){
    scaX_con = apply(X_con, 2, scale)
    scaXobs = cbind(X_bin, scaX_con)
  } else if (p_X2>0){
    scaXobs = apply(Xobs, 2, scale)
  } else {
    scaXobs = Xobs
  }
  
  matX = cbind(1, Zobs, scaXobs) # Design matrix for outcome model
  p_matX = ncol(matX)
  matXtX = t(matX) %*% matX
  INVmatXtX = inv_cpp(matXtX)
  
  # =============================================================================
  # 2. HYPERPARAMETERS (PRIORS)
  # =============================================================================
  # --- Priors for outcome model (Normal-Sca-InvChi2-) ---
  library(survival)
  survreg_beta = survreg(Surv(object$Tobs, Dobs) ~ matX - 1, dist = "lognormal")
  a_beta = coef(survreg_beta)
  a_beta[is.na(a_beta)] = 0
  p_beta = length(a_beta)
  B_beta = vcov(survreg_beta)[-p_beta, -p_beta]
  c_beta = N/5
  B_beta = c_beta * B_beta
  Binv_beta = inv_cpp(B_beta)
  aBinv_beta = Binv_beta %*% a_beta
  # a_sig2 = 1
  # b_sig2 = (survreg_beta$scale^2)
  # library(aftgee)
  # aft_fit = aftgee(Surv(object$Tobs, Dobs) ~ matX - 1)
  # a_beta = coef(aft_fit)
  # a_beta[is.na(a_beta)] = 0
  # B_beta = vcov(aft_fit)
  # c_beta = N/5
  # B_beta = c_beta * B_beta
  # Binv_beta = inv_cpp(B_beta)
  # aBinv_beta = Binv_beta %*% a_beta
  # a_sig2 = 1
  # b_sig2 = var(residuals_aft[Dobs==1]) * (a_sig2 + 2)
  a_sig2 = 3
  b_sig2 = 0.01
  
  # --- Priors for treatment and covariate models (eta and omega) ---
  # Prior for binary treatment: p ~ Beta(a, b)
  a_eta = 1
  b_eta = 1
  # Prior for binary covariates: p ~ Beta(a, b)
  a_pi = 1
  b_pi = 1
  # Prior for continuous covariates: (mu, tau2) ~ Normal-Sca-InvChi2
  a_tau2 = 1
  b_tau2 = 1
  b_mu = 0.5
  a_mu = 0
  
  # --- Prior for DP concentration parameter alpha ---
  a_alpha = 1
  b_alpha = 1
  
  # =============================================================================
  # 3. INITIAL VALUES
  # =============================================================================
  # Set initial values for cluster membership
  # Use K-means to find 5 whole data clusters
  # initial cluster size K = 5
  S = kmeans(cbind(Tobs, Dobs, Zobs, scaXobs), 5)$cluster
  
  # Determine unique clusters from the cluster membership variable, S
  # Sort clusters
  unique_S = cbind(sort(unique(S)))
  
  # Calculate number of clusters
  K = length(unique_S)
  
  # Initialize alpha parameters
  alpha = 2
  
  # --- Initialize outcome parameters (theta) ---
  sig2Pars = t(rscainvchisq_cpp(K, a_sig2, b_sig2))
  betaPars = sapply(1:K, function(k) rmvn_cpp(1, a_beta, sig2Pars[k] * B_beta))
  
  # ... (Initialize latent log-times y, y*) ...
  Tstar = numeric(N)
  Tstar[Dobs == 1] = Tobs[Dobs == 1]
  for (k in 1:K) {
    ind_curr = which(S == unique_S[k])
    ind_cens = ind_curr[Dobs[ind_curr] == 0]
    
    if(length(ind_cens) > 0) {
      y_full = Tstar[ind_curr] 
      X_full = matX[ind_curr, , drop=FALSE]
      
      # Update Parameters (Collapsed Gibbs step)
      Binv_n  = Binv_beta + crossprod(X_full)
      Bn_beta = inv_cpp(Binv_n)
      an_beta = Bn_beta %*% (aBinv_beta + crossprod(X_full, y_full))
      
      residuals = y_full - X_full %*% an_beta
      beta_diff = an_beta - a_beta
      
      sse_data  = sum(residuals^2)
      sse_prior = as.numeric(t(beta_diff) %*% Binv_beta %*% beta_diff)
      
      an_sig2 = a_sig2 + length(y_full)
      bn_sig2 = b_sig2 + sse_data + sse_prior
      
      # Prepare Predictive Distribution
      X_cens = matX[ind_cens, , drop=FALSE]
      lower_vec = Tobs[ind_cens]
      
      quad_form = rowSums((X_cens %*% Bn_beta) * X_cens)
      scale_val = sqrt((bn_sig2 / an_sig2) * (1 + quad_form))
      mu_val    = as.vector(X_cens %*% an_beta)
      
      # Calculate Dynamic Upper Bound (Recommendation 2)
      # Use 99.9th percentile of the current predictive t-distribution
      t_quant   = qt(0.999, df = an_sig2)
      upper_vec = mu_val + (t_quant * scale_val)
      
      # Safety: Ensure upper > lower + buffer (4*scale) to prevent numerical trapping
      upper_safe = pmax(upper_vec, lower_vec + 4 * scale_val)
      
      # Sample
      # Note: 'upper' is now a vector, so it moves out of MoreArgs
      Tstar[ind_cens] = mapply(rtrunct_cpp, 
                               mu = mu_val, 
                               scale = scale_val, 
                               lower = lower_vec, 
                               upper = upper_safe, 
                               MoreArgs = list(n = 1, df = an_sig2))
    }
  }
  
  # --- Initialize propensity parameters (eta) ---
  etaPars = matrix(nrow = 1, ncol = K)
  # --- Initialize covariate parameters (omega) ---
  piPars = matrix(nrow = p_X1, ncol = K) # Probabilities
  muPars = matrix(nrow = p_X2, ncol = K) # Means
  tau2Pars = matrix(nrow = p_X2, ncol = K) # Variances
  
  for (k in 1:K) {
    ind_S = which(S == unique_S[k])
    n_temp = length(ind_S)
    T_temp = Tstar[ind_S]
    matX_temp = matrix(matX[ind_S,], ncol = p_matX)
    
    # update parameters for true failure time
    regression_results = update_reg_con(T_temp, matX_temp, betaPars[,k],
                                        a_beta, Binv_beta, aBinv_beta, a_sig2, b_sig2)
    sig2Pars[,k] = regression_results$sig2_par
    betaPars[,k] = regression_results$beta_par
    # update parameters for binary treatment
    etaPars[,k] = update_conf_bin(Zobs[ind_S], a_eta, b_eta)
    # update parameters for binary covariates
    if (p_X1 > 0) {
      for (q in 1:p_X1) {
        # posterior is beta
        bivariate_pars = update_conf_bin(scaXobs[ind_S, q], a_pi, b_pi)
        piPars[q, k] = bivariate_pars
      }
    }
    # update parameters for continuous covariates
    if (p_X2 > 0) {
      for (q in 1:p_X2) {
        continuous_pars = update_conf_con(scaXobs[ind_S, (p_X1 + q)], a_tau2, b_tau2, a_mu, b_mu)
        tau2Pars[q, k] = continuous_pars$var_par
        muPars[q, k] = continuous_pars$mean_par
      }
    }
  }
  
  # num_MC_prior = 1e5
  # Calculate f0_x and f0_y for use in cluster function in Gibbs Sampler
  # Average confounders distribution over prior for each x_i
  ab_tau = a_tau2 * b_tau2
  a_tau_new = (a_tau2 + 1)/2
  a_tau_half = a_tau2/2
  b_mu_ratio = b_mu/(b_mu+1)
  margin_part1 = (gamma(a_tau_new)/gamma(a_tau_half)) * sqrt(b_mu_ratio/pi) * (ab_tau)^{a_tau_half}
  
  p_Z = 1
  f0_x_all = matrix(0, nrow = N, ncol = p_Z + p_X)
  # Binary treatment
  # Beta-Binomial
  for (q in 1:p_Z) {
    f0_x_all[, q] = beta(a_eta + Zobs[q], b_eta - Zobs[q] + 1) # beta(a_pi,b_pi) = 1 in this case
  }
  # Binary confounders
  # Beta-Binomial
  if (p_X1>0){
    for (q in 1:p_X1) {
      f0_x_all[, (p_Z + q)] = beta(a_pi + scaXobs[, q], b_pi - scaXobs[, q] + 1) # beta(a_pi,b_pi) = 1 in this case
    }
  }
  # Continuous confounders
  if (p_X2>0){
    for (q in 1:p_X2) {
      x_temp = scaXobs[, (p_X1 + q)]
      margin_part2 = (ab_tau + b_mu_ratio * (x_temp-a_mu)^{2})^(-a_tau_new)
      f0_x_all[, (p_Z + p_X1 + q)] = margin_part1 * margin_part2
    }
  }
  
  # Take product (confounders are assumed to be locally independent).
  # Result is vector of size N
  f0_x = apply(f0_x_all, 1, prod)
  
  # =============================================================================
  # 4. GIBBS SAMPLING LOOP
  # =============================================================================
  gibbs_total = gibbs_iter + gibbs_burnin
  n_store = floor(gibbs_iter / gibbs_thin)
  
  # Define iteration check
  iter_check = floor(gibbs_iter/20)
  
  # Make vectors to store draws from Gibbs Sampler
  n_store = floor(gibbs_iter/gibbs_thin)
  
  # Make lists to store draws from Gibbs Sampler
  alpha_draws = numeric(n_store)
  
  betaLists = list(NA)
  sig2Lists = list(NA)
  etaLists  = list(NA)
  piLists   = list(NA)
  muLists   = list(NA)
  tau2Lists = list(NA)
  n_kLists  = list(NA)
  KLists    = list(NA)
  
  count_it = 1
  # End initial values ------------------------------------------------------------
  
  # First draw each parameter for BNP model. Then calculate causal effect.
  for (gibbs_reps in 1:gibbs_total) {
    # gibbs_reps = 1
    # Update of all parameters in BNP model--------------------------------
    # Update cluster membership -------------------------------------------
    # sourceCpp("BNPsurvival_cpp.cpp")
    cluster_res = DPMMcluster_cpp(N, p_X1, p_X2, 
                                  Tstar, matX, Zobs, scaXobs,
                                  S, unique_S, alpha, 
                                  betaPars, sig2Pars, etaPars, 
                                  piPars, muPars, tau2Pars,
                                  a_beta, B_beta, Binv_beta, 
                                  a_sig2, b_sig2, 
                                  a_eta, b_eta,
                                  a_pi, b_pi, 
                                  a_mu, b_mu, 
                                  a_tau2, b_tau2)
    # Store cluster membership output from cluster function
    S = cluster_res$S
    
    betaPars = cluster_res$betaPars
    sig2Pars = cluster_res$sig2Pars
    etaPars  = cluster_res$etaPars
    piPars   = cluster_res$piPars
    muPars   = cluster_res$muPars
    tau2Pars = cluster_res$tau2Pars
    
    # Make vector of y-clusters
    unique_S = cluster_res$unique_S
    
    # Calculate the number of clusters
    K = length(unique_S)
    
    # Calculate number of subjects in each cluster and store in vector n_k.
    # Use k to index clusters.
    n_k = rep(0, K)
    for (k in 1:K) {
      ind_k = (S == unique_S[k])
      n_k[k] = sum(ind_k)
    }
    # End update of cluster membership --------------------------------------------
    
    # Update of cluster-specific parameters for confounders---------------------
    Tstar[Dobs == 1] = Tobs[Dobs == 1]
    for (k in 1:K) {
      ind_curr = which(S == unique_S[k])
      ind_cens = ind_curr[Dobs[ind_curr] == 0]
      
      if(length(ind_cens) > 0) {
        y_full = Tstar[ind_curr] 
        X_full = matX[ind_curr, , drop=FALSE]
        
        # Update Parameters (Collapsed Gibbs step)
        Binv_n  = Binv_beta + crossprod(X_full)
        Bn_beta = inv_cpp(Binv_n)
        an_beta = Bn_beta %*% (aBinv_beta + crossprod(X_full, y_full))
        
        residuals = y_full - X_full %*% an_beta
        beta_diff = an_beta - a_beta
        
        sse_data  = sum(residuals^2)
        sse_prior = as.numeric(t(beta_diff) %*% Binv_beta %*% beta_diff)
        
        an_sig2 = a_sig2 + length(y_full)
        bn_sig2 = b_sig2 + sse_data + sse_prior
        
        # Prepare Predictive Distribution
        X_cens = matX[ind_cens, , drop=FALSE]
        lower_vec = Tobs[ind_cens]
        
        quad_form = rowSums((X_cens %*% Bn_beta) * X_cens)
        scale_val = sqrt((bn_sig2 / an_sig2) * (1 + quad_form))
        mu_val    = as.vector(X_cens %*% an_beta)
        
        # Calculate Dynamic Upper Bound (Recommendation 2)
        # Use 99.9th percentile of the current predictive t-distribution
        t_quant   = qt(0.999, df = an_sig2)
        upper_vec = mu_val + (t_quant * scale_val)
        
        # Safety: Ensure upper > lower + buffer (4*scale) to prevent numerical trapping
        upper_safe = pmax(upper_vec, lower_vec + 4 * scale_val)
        
        # Sample
        # Note: 'upper' is now a vector, so it moves out of MoreArgs
        Tstar[ind_cens] = mapply(rtrunct_cpp, 
                                 mu = mu_val, 
                                 scale = scale_val, 
                                 lower = lower_vec, 
                                 upper = upper_safe, 
                                 MoreArgs = list(n = 1, df = an_sig2))
      }
    }
    
    for (k in 1:K) {
      ind_S = which(S == unique_S[k])
      n_temp = length(ind_S)
      T_temp = Tstar[ind_S]
      matX_temp = matrix(matX[ind_S,], ncol = p_matX)
      
      # update parameters for true failure time
      regression_results = update_reg_con(T_temp, matX_temp, betaPars[,k],
                                          a_beta, Binv_beta, aBinv_beta, a_sig2, b_sig2)
      sig2Pars[,k] = regression_results$sig2_par
      betaPars[,k] = regression_results$beta_par
      # update parameters for binary treatment
      etaPars[,k] = update_conf_bin(Zobs[ind_S], a_eta, b_eta)
      # update parameters for binary covariates
      if (p_X1 > 0) {
        for (q in 1:p_X1) {
          # posterior is beta
          bivariate_pars = update_conf_bin(scaXobs[ind_S, q], a_pi, b_pi)
          piPars[q, k] = bivariate_pars
        }
      }
      # update parameters for continuous covariates
      if (p_X2 > 0) {
        for (q in 1:p_X2) {
          continuous_pars = update_conf_con(scaXobs[ind_S, (p_X1 + q)], a_tau2, b_tau2, a_mu, b_mu)
          tau2Pars[q, k] = continuous_pars$var_par
          muPars[q, k] = continuous_pars$mean_par
        }
      }
    }
    # End update of cluster-specific parameters for confounders---------------------
    
    # Update concentration parameters--------------------------------------------
    alpha = update_alpha_cpp(N, K, alpha, a_alpha, b_alpha)
    # End update of concentration parameters-----------------------------------
    # End update of all parameters in BNP model--------------------------------
    
    if (gibbs_reps < gibbs_burnin) {
    } else if (gibbs_reps == gibbs_burnin) {
      cat("Bur-In End",gibbs_reps,"Time:",date(),"\n")
    } else if (gibbs_reps > gibbs_burnin) {
      if (gibbs_reps %% gibbs_thin == 0) {
        n_kLists[[count_it]]  = n_k
        KLists[[count_it]]    = K
        betaLists[[count_it]] = betaPars
        sig2Lists[[count_it]] = sig2Pars
        etaLists[[count_it]]  = etaPars
        piLists[[count_it]]   = piPars
        muLists[[count_it]]   = muPars
        tau2Lists[[count_it]] = tau2Pars
        
        alpha_draws[count_it] = alpha
        
        count_it = count_it + 1
      }
      
      if (gibbs_reps %% iter_check == 0) {
        cat("Gibbs Iteration",(gibbs_reps-gibbs_burnin),"(",(gibbs_reps-gibbs_burnin)/gibbs_iter*100,"%)","Time:",date(),"\n")
      }
    }
  }
  # =============================================================================
  # 5. RETURN RESULTS
  # =============================================================================
  
  # constants
  # wb
  constants = list(N = N, n_MCMC = n_store,
                   p_matX = p_matX, p_X1 = p_X1, p_X2 = p_X2)
  
  # priors
  priors = list(a_sig2 = a_sig2, b_sig2 = b_sig2, 
                a_beta = a_beta, B_beta = B_beta,
                a_eta = a_eta, b_eta = b_eta,
                a_pi = a_pi, b_pi = b_pi, 
                a_mu = a_mu, b_mu = b_mu, 
                a_tau2 = a_tau2, b_tau2 = b_tau2)
  
  # # Precalculated values for Post-Processing steps
  # MCpriors = list(f0_z0_mc = f0_z0_mc, 
  #                 f0_z1_mc = f0_z1_mc,
  #                 TbetaPrior_mc = TbetaPrior_mc, 
  #                 MbetaPrior_mc = MbetaPrior_mc, 
  #                 Msig2Prior_mc = Msig2Prior_mc,
  #                 VbetaPrior_mc = VbetaPrior_mc, 
  #                 Vsig2Prior_mc = Vsig2Prior_mc)
  
  # MCMC Posteriors
  MCMCposteriors = list(betaLists   = betaLists,
                        sig2Lists   = sig2Lists,
                        etaLists    = etaLists,
                        piLists     = piLists,
                        muLists     = muLists,
                        tau2Lists   = tau2Lists,
                        n_kLists    = n_kLists,
                        KLists      = KLists,
                        alpha_draws = alpha_draws)
  
  # MCMC results
  MCMCresult = list(priors = priors,
                    # MCpriors = MCpriors,
                    constants = constants,
                    MCMCposteriors = MCMCposteriors)
  
  # {
  #   plot(unlist(KLists), type="s")
  #   
  #   # Setup the base plot with the first density
  #   plot(density(Tstar),
  #        xlim = c(-12, 15),
  #        ylim = c(0, 0.2),
  #        main = "Density Comparisons",
  #        xlab = "Log Time",
  #        col = "black",
  #        lty = 1,
  #        lwd = 2);
  #   
  #   # Add the remaining densities using lines
  #   lines(density(temp_data$logY),
  #         col = "black",
  #         lty = 1,
  #         lwd = 2);
  #   
  #   lines(density(Tstar[Zobs == 0]),
  #         col = "blue",
  #         lty = 2,
  #         lwd = 2);
  #   
  #   lines(density(temp_data$logY0),
  #         col = "blue",
  #         lty = 2,
  #         lwd = 2);
  #   
  #   lines(density(Tstar[Zobs == 1]),
  #         col = "red",
  #         lty = 3,
  #         lwd = 2);
  #   
  #   lines(density(temp_data$logY1),
  #         col = "red",
  #         lty = 3,
  #         lwd = 2);
  #   
  #   # Add a legend to distinguish the lines
  #   legend("topright",
  #          legend = c("Tstar (All)", "logY (All)", "Tstar (Z=0)", "logY0 (Z=0)", "Tstar (Z=1)", "logY1 (Z=1)"),
  #          col = c("black", "black", "blue", "blue", "red", "red"),
  #          lty = c(1, 1, 2, 2, 3, 3),
  #          lwd = 2,
  #          cex = 0.8);
  # }
  
  return(MCMCresult)
}

# -----------------------------------------------------------------------------
# Define function for the post-processing under the BNP (DPMM) model
SURVIVAL_DPMM_POST = function(object, nu, rho, M, esttype = NULL){
  # -----------------------------------------------------------------------------
  # object = DPMM_results_MCMC
  
  # -----------------------------------------------------------------------------
  normalize_log_probs = function(log_prob_matrix) {
    max_log_probs = apply(log_prob_matrix, 1, max)
    log_probs_shifted = log_prob_matrix - max_log_probs
    probs_shifted = exp(log_probs_shifted)
    row_sums = rowSums(probs_shifted)
    
    # Handle rows that are all -Inf (row_sum = 0)
    zero_rows = which(row_sums == 0)
    if(length(zero_rows) > 0) {
      probs_shifted[zero_rows,] = 1 / ncol(log_prob_matrix) # Assign uniform prob
      row_sums[zero_rows] = 1
    }
    return(probs_shifted / row_sums)
  }
  
  # -----------------------------------------------------------------------------
  update_parameter = function(S_x_mc, K, 
                              betaPars, sig2Pars, etaPars,
                              piPars, muPars, tau2Pars, 
                              priors, constants) {
    # Extract variables from lists for easier reading
    a_beta = priors$a_beta
    B_beta = priors$B_beta
    a_sig2 = priors$a_sig2
    b_sig2 = priors$b_sig2
    a_eta  = priors$a_eta
    b_eta  = priors$b_eta
    a_pi   = priors$a_pi
    b_pi   = priors$b_pi
    a_mu   = priors$a_mu
    b_mu   = priors$b_mu
    a_tau2 = priors$a_tau2
    b_tau2 = priors$b_tau2
    
    p_X1 = constants$p_X1
    p_X2 = constants$p_X2
    
    # Get properties of the new cluster assignments
    unique_S_x_mc = sort(unique(S_x_mc))
    K_x_mc = length(unique_S_x_mc)
    max_K_mc = max(unique_S_x_mc)
    
    n_k_mc = numeric(K_x_mc)
    
    # ---------------------------------------------------------------------------
    # --- Step A: Check if a NEW cluster was born (index > K) ---
    # ---------------------------------------------------------------------------
    if (max_K_mc > K) {
      # We need to generate ONE set of new parameters from the priors
      
      # Outcome Parameters
      sig2_new = c(rscainvchisq_cpp(1, a_sig2, b_sig2))
      # Assuming rmvn_cpp(mean, cov) is your custom function for a single draw
      beta_new = c(rmvn_cpp(1, a_beta, B_beta)) 
      
      # Treatment Parameter
      eta_new = rbeta(1, a_eta, b_eta)
      
      # Covariate Parameters
      pi_new   = numeric(p_X1)
      mu_new   = numeric(p_X2)
      tau2_new = numeric(p_X2)
      
      if (p_X1 > 0) {
        pi_new = rbeta(p_X1, a_pi, b_pi)
      }
      if (p_X2 > 0) {
        for (q in 1:p_X2) {
          tau2_new[q] = rinvgamma_cpp(1, a_tau2, b_tau2)
          mu_new[q]   = rnorm(1, mean = a_mu, sd = sqrt(tau2_new[q] / b_mu))
        }
      }
      
      # Append new parameters to the existing matrices
      betaPars_mc = cbind(betaPars, beta_new)
      sig2Pars_mc = cbind(sig2Pars, sig2_new)
      etaPars_mc  = cbind(etaPars, eta_new)
      piPars_mc   = cbind(piPars, pi_new)
      muPars_mc   = cbind(muPars, mu_new)
      tau2Pars_mc = cbind(tau2Pars, tau2_new)
      
    } else {
      # No new cluster created, keep existing parameters temporarily
      betaPars_mc = betaPars
      sig2Pars_mc = sig2Pars
      etaPars_mc  = etaPars
      piPars_mc   = piPars
      muPars_mc   = muPars
      tau2Pars_mc = tau2Pars
    }
    
    # ---------------------------------------------------------------------------
    # --- Step B: Handle Deleted/Empty Clusters ---
    # ---------------------------------------------------------------------------
    # We only keep the columns corresponding to clusters that *actually appear*
    # Subset matrices to keep only active clusters
    betaPars_mc = matrix(betaPars_mc[, unique_S_x_mc], ncol = K_x_mc)
    sig2Pars_mc = matrix(sig2Pars_mc[, unique_S_x_mc], ncol = K_x_mc)
    etaPars_mc  = matrix(etaPars_mc[, unique_S_x_mc],  ncol = K_x_mc)
    piPars_mc   = matrix(piPars_mc[, unique_S_x_mc],   ncol = K_x_mc)
    muPars_mc   = matrix(muPars_mc[, unique_S_x_mc],   ncol = K_x_mc)
    tau2Pars_mc = matrix(tau2Pars_mc[, unique_S_x_mc], ncol = K_x_mc)
    
    # Update cluster counts for the new set
    for (k in 1:length(unique_S_x_mc)) {
      n_k_mc[k] = sum(S_x_mc == unique_S_x_mc[k])
    }
    
    # ---------------------------------------------------------------------------
    # --- Step C: Remap Cluster Indices ---
    # ---------------------------------------------------------------------------
    # Re-label S_x_mc to be contiguous (1, 2, ..., K_x_mc)
    S_x_mc_remapped = as.numeric(factor(S_x_mc, levels = unique_S_x_mc))
    
    # Return all the updated components
    return(list(
      S_x_mc = S_x_mc_remapped,
      unique_S_x_mc = 1:K_x_mc,
      n_k_mc = n_k_mc,
      K_x_mc = K_x_mc,
      betaPars_mc = betaPars_mc,
      sig2Pars_mc = sig2Pars_mc,
      etaPars_mc = etaPars_mc,
      piPars_mc = piPars_mc,
      muPars_mc = muPars_mc,
      tau2Pars_mc = tau2Pars_mc
    ))
  }
  
  # -----------------------------------------------------------------------------
  # Extract Constants and Settings
  # Ensure nu and rho are vectors
  nu_vec = as.vector(nu)
  rho_vec = as.vector(rho)
  n_nu = length(nu_vec)
  n_rho = length(rho_vec)
  
  level = 0.05
  quantile_alpha = c(level/2,1-level/2)
  z0 = 0
  z1 = 1
  
  N = object$constants$N
  num_MC = object$constants$num_MC
  n_MCMC = object$constants$n_MCMC
  p_matX = object$constants$p_matX
  p_X1 = object$constants$p_X1
  p_X2 = object$constants$p_X2
  
  a_beta = object$priors$a_beta
  B_beta = object$priors$B_beta
  a_sig2 = object$priors$a_sig2
  b_sig2 = object$priors$b_sig2
  a_eta = object$priors$a_eta
  b_eta = object$priors$b_eta
  a_pi = object$priors$a_pi
  b_pi = object$priors$b_pi
  a_mu = object$priors$a_mu
  b_mu = object$priors$b_mu
  a_tau2 = object$priors$a_tau2
  b_tau2 = object$priors$b_tau2
  
  # -----------------------------------------------------------------------------
  # Define interation check
  iter_check = floor(n_MCMC/10)
  
  Y_rho_nu_z0x_store = array(NA, dim = c(n_MCMC, n_nu, n_rho))
  Y_rho_nu_z1x_store = array(NA, dim = c(n_MCMC, n_nu, n_rho))
  for (post_reps in 1:n_MCMC) {
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    alpha = object$MCMCposteriors$alpha_draws[post_reps]
    K = object$MCMCposteriors$KLists[[post_reps]]
    n_k = object$MCMCposteriors$n_kLists[[post_reps]]
    betaPars = object$MCMCposteriors$betaLists[[post_reps]]
    sig2Pars = object$MCMCposteriors$sig2Lists[[post_reps]]
    etaPars = object$MCMCposteriors$etaLists[[post_reps]]
    piPars = object$MCMCposteriors$piLists[[post_reps]]
    muPars = object$MCMCposteriors$muLists[[post_reps]]
    tau2Pars = object$MCMCposteriors$tau2Lists[[post_reps]]
    
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Draw new clusters for x
    lambda_x_mc = n_k
    S_x_mc = sapply(1:M, function(l) rmultinom_cpp(lambda_x_mc))
    unique_S_x_mc = sort(unique(S_x_mc))
    K_mc = max(unique_S_x_mc)
    
    # -----------------------------------------------------------------------------
    # Draw confounders for each of the M observations
    n_k_mc = numeric(K)
    X_mc = matrix(nrow = M, ncol = (p_X1 + p_X2))
    for (k in 1:K) {
      ind_S_x_mc = which(S_x_mc == unique_S_x_mc[k])
      n_k_temp = length(ind_S_x_mc)
      n_k_mc[k] = n_k_temp
      
      if (p_X1 > 0) {
        for (q in 1:p_X1) {
          X_mc[ind_S_x_mc, q] = rbinom(n_k_temp, 1, piPars[q,k]);
        }
      }
      # update parameters for continuous covariates
      if (p_X2 > 0) {
        for (q in 1:p_X2) {
          X_mc[ind_S_x_mc, (p_X1 + q)] = rnorm(n_k_temp, muPars[q,k], sqrt(tau2Pars[q,k]));
        }
      }
    }
    
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    lambda_Y_z0x_mc = matrix(nrow = M, ncol = K)
    lambda_Y_z1x_mc = matrix(nrow = M, ncol = K)
    for (k in 1:K) {
      ind_S_x_mc = which(S_x_mc == unique_S_x_mc[k])
      
      log_w_k_temp = log(n_k[k]/N)
      
      prob_z0_mc = dbinom(z0, 1, etaPars[,k], log = TRUE)
      prob_z1_mc = dbinom(z1, 1, etaPars[,k], log = TRUE)
      
      prob_x_mc = numeric(M)
      if (p_X1 > 0) {
        for (q in 1:p_X1) {
          prob_x_mc = prob_x_mc + dbinom(X_mc[, q], 1, piPars[,k], log = TRUE)
        }
      }
      # update parameters for continuous covariates
      if (p_X2 > 0) {
        for (q in 1:p_X2) {
          prob_x_mc = prob_x_mc + dnorm(X_mc[, (p_X1 + q)], muPars[q,k], sqrt(tau2Pars[q,k]), log = TRUE)
        }
      }
      lambda_Y_z0x_mc[,k] = log_w_k_temp + prob_z0_mc + prob_x_mc
      lambda_Y_z1x_mc[,k] = log_w_k_temp + prob_z1_mc + prob_x_mc
    }
    lambda_Y_z0x_mc = normalize_log_probs(lambda_Y_z0x_mc)
    lambda_Y_z1x_mc = normalize_log_probs(lambda_Y_z1x_mc)
    
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Create the design matrices for all M subjects
    matX_z0x_mc = cbind(1, z0, X_mc) # Intercept, Zobs=0, Covariates
    matX_z1x_mc = cbind(1, z1, X_mc) # Intercept, Zobs=1, Covariates
    MUy_z0x_mc = (matX_z0x_mc %*% betaPars)
    MUy_z1x_mc = (matX_z1x_mc %*% betaPars)
    Ey_z0x_mc = rowSums(MUy_z0x_mc * lambda_Y_z0x_mc)
    Ey_z1x_mc = rowSums(MUy_z1x_mc * lambda_Y_z1x_mc)
    
    # -----------------------------------------------------------------------------
    lambda_Y_z0x = apply(lambda_Y_z0x_mc, 2, mean)
    lambda_Y_z1x = apply(lambda_Y_z1x_mc, 2, mean)
    MUy_z0x = apply(MUy_z0x_mc, 2, mean)
    MUy_z1x = apply(MUy_z1x_mc, 2, mean)
    Ey_z0x = mean(Ey_z0x_mc)
    Ey_z1x = mean(Ey_z1x_mc)
    
    # LOOP OVER nu
    for(ii_nu in 1:n_nu) {
      curr_nu = nu_vec[ii_nu]
      curr_log_nu = log(curr_nu)
      
      S_nu_z0x_vec = sapply(1:M, function(m) 
        S_DPMM_cpp(curr_log_nu, lambda_Y_z0x_mc[m,], 
                   MUy_z0x_mc[m,], sqrt(sig2Pars), TRUE))
      S_nu_z1x_vec = sapply(1:M, function(m) 
        S_DPMM_cpp(curr_log_nu, lambda_Y_z1x_mc[m,], 
                   MUy_z1x_mc[m,], sqrt(sig2Pars), TRUE))
      S_nu_z0x = mean(S_nu_z0x_vec)
      S_nu_z1x = mean(S_nu_z1x_vec)
      
      # LOOP OVER rho
      for(ii_rho in 1:n_rho) {
        curr_rho = rho_vec[ii_rho]
        
        target_rho_nu_z0x = S_nu_z0x + log(1 - curr_rho)
        target_rho_nu_z1x = S_nu_z1x + log(1 - curr_rho)
        
        # -----------------------------------------------------------------------------
        Y_rho_nu_z0x = S_optim_cpp(target_rho_nu_z0x,
                                      curr_log_nu,
                                      Ey_z0x,
                                      Ey_z0x - 1,
                                      Ey_z0x + 1,
                                      lambda_Y_z0x,
                                      MUy_z0x,
                                      sqrt(sig2Pars),
                                      TRUE)$optimizer
        Y_rho_nu_z1x = S_optim_cpp(target_rho_nu_z1x,
                                      curr_log_nu,
                                      Ey_z1x,
                                      Ey_z1x - 1,
                                      Ey_z1x + 1,
                                      lambda_Y_z1x,
                                      MUy_z1x,
                                      sqrt(sig2Pars),
                                      TRUE)$optimizer
        
        # -----------------------------------------------------------------------------
        Y_rho_nu_z0x_store[post_reps, ii_nu, ii_rho] = exp(Y_rho_nu_z0x) - curr_nu
        Y_rho_nu_z1x_store[post_reps, ii_nu, ii_rho] = exp(Y_rho_nu_z1x) - curr_nu
      } # End rho loop
    } # End nu loop
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  E_rho_nu_diff_store = Y_rho_nu_z1x_store - Y_rho_nu_z0x_store
  
  z0_tables_list   = list()
  z1_tables_list   = list()
  diff_tables_list = list()
  counter = 1
  for(ii_nu in 1:n_nu) {
    for(ii_rho in 1:n_rho) {
      val_nu  = nu_vec[ii_nu]
      val_rho = rho_vec[ii_rho]
      
      draws_z0   = Y_rho_nu_z0x_store[, ii_nu, ii_rho]
      draws_z1   = Y_rho_nu_z1x_store[, ii_nu, ii_rho]
      draws_diff = E_rho_nu_diff_store[, ii_nu, ii_rho]
      
      summ_z0   = POSTsummary(draws_z0, esttype, quantile_alpha)
      summ_z1   = POSTsummary(draws_z1, esttype, quantile_alpha)
      summ_diff = POSTsummary(draws_diff, esttype, quantile_alpha)
      
      summ_z0   = cbind(nu = val_nu, rho = val_rho, summ_z0)
      summ_z1   = cbind(nu = val_nu, rho = val_rho, summ_z1)
      summ_diff = cbind(nu = val_nu, rho = val_rho, summ_diff)
      
      z0_tables_list[[counter]]   = summ_z0
      z1_tables_list[[counter]]   = summ_z1
      diff_tables_list[[counter]] = summ_diff
      
      counter = counter + 1
    }
  }
  
  Y_rho_nu_z0x_result  = do.call(rbind, z0_tables_list)
  Y_rho_nu_z1x_result  = do.call(rbind, z1_tables_list)
  E_rho_nu_diff_result = do.call(rbind, diff_tables_list)
  
  POSTresult = list(
    Y_rho_nu_z0x_result = Y_rho_nu_z0x_result, 
    Y_rho_nu_z1x_result = Y_rho_nu_z1x_result, 
    E_rho_nu_diff_result = E_rho_nu_diff_result,
    Y_rho_nu_z0x_store = Y_rho_nu_z0x_store,
    Y_rho_nu_z1x_store = Y_rho_nu_z1x_store,
    E_rho_nu_diff_store = E_rho_nu_diff_store
  )
  
  return(POSTresult)
}

# ------------------------------------------------------------------------
SURVIVAL_DPMM = function(object, nu = 0, rho = 0.5, M = 2e3, 
                         gibbs_iter = 2e4, gibbs_burnin = 2e4, gibbs_thin = 1e2, 
                         esttype = "mean", saveall = FALSE){
  # object = temp_data
  MCMCresult = SURVIVAL_DPMM_MCMC(object, gibbs_iter, gibbs_burnin, gibbs_thin)
  POSTresult = SURVIVAL_DPMM_POST(MCMCresult, nu, rho, M, esttype)
  if (isTRUE(saveall)){
    result = c(MCMCresult, POSTresult)
  } else {
    result = POSTresult
  }
  return(result)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Define function to fit the BNP (EDPMM) model
SURVIVAL_EDPMM_MCMC = function(object,
                               gibbs_iter = 1e4, gibbs_burnin = 1e4, gibbs_thin = 1e1, 
                               num_MC = NULL, num_MC_prior = NULL) {
  # =============================================================================
  # 1. SETUP & INITIALIZATION
  # =============================================================================
  # object = temp_data
  Xobs = object$Xobs
  Zobs = object$Zobs
  Dobs = object$Dobs
  Tobs = object$Tobs
  
  N = length(Tobs)
  Toriginal = Tobs
  Tobs = log(Tobs)
  
  # Identify covariate types (binary vs. continuous)
  is_binary = apply(Xobs, 2, function(col) all(col %in% c(0, 1)))
  X_bin = Xobs[, is_binary, drop = FALSE]
  X_con = Xobs[, !is_binary, drop = FALSE]
  p_X1 = ncol(X_bin)
  p_X2 = ncol(X_con)
  p_X = p_X1 + p_X2
  p_Z = length(Zobs)/N
  
  # Standardize continuous confounders (important when choosing priors)
  if (p_X1>0 && p_X2>0){
    scaX_con = apply(X_con, 2, scale)
    scaXobs = cbind(X_bin, scaX_con)
  } else if (p_X2>0){
    scaXobs = apply(Xobs, 2, scale)
  } else {
    scaXobs = Xobs
  }
  
  matX = cbind(1, Zobs, scaXobs) # Design matrix for outcome model
  p_matX = ncol(matX)
  matXtX = t(matX) %*% matX
  INVmatXtX = inv_cpp(matXtX)
  
  # =============================================================================
  # 2. HYPERPARAMETERS (PRIORS)
  # =============================================================================
  # --- Priors for outcome model (Normal-Sca-InvChi2-) ---
  library(survival)
  survreg_beta = survreg(Surv(object$Tobs, Dobs) ~ matX - 1, dist = "lognormal")
  a_beta = coef(survreg_beta)
  a_beta[is.na(a_beta)] = 0
  p_beta = length(a_beta)
  B_beta = vcov(survreg_beta)[-p_beta, -p_beta]
  c_beta = N/5
  B_beta = c_beta * B_beta
  Binv_beta = inv_cpp(B_beta)
  aBinv_beta = Binv_beta %*% a_beta
  # a_sig2 = 0.01
  # b_sig2 = (survreg_beta$scale^2)
  # library(aftgee)
  # aft_fit = aftgee(Surv(object$Tobs, Dobs) ~ matX - 1)
  # a_beta = coef(aft_fit)
  # a_beta[is.na(a_beta)] = 0
  # B_beta = vcov(aft_fit)
  # c_beta = N/5
  # B_beta = c_beta * B_beta
  # Binv_beta = inv_cpp(B_beta)
  # aBinv_beta = Binv_beta %*% a_beta
  # a_sig2 = 1
  # b_sig2 = var(residuals_aft[Dobs==1]) * (a_sig2 + 2)
  a_sig2 = 3
  b_sig2 = 0.01
  
  # --- Priors for treatment and covariate models (eta and omega) ---
  # Prior for binary treatment: p ~ Beta(a, b)
  a_eta = 1
  b_eta = 1
  # Prior for binary covariates: p ~ Beta(a, b)
  a_pi = 1
  b_pi = 1
  # Prior for continuous covariates: (mu, tau2) ~ Normal-Sca-InvChi2
  a_tau2 = 1
  b_tau2 = 1
  b_mu = 0.5
  a_mu = 0
  
  # --- Prior for DP concentration parameter alpha ---
  # for concentration parameters
  # randg: shape and scale at default in Rcpp & RcppArmadillo
  # rgamma: shape and rate at default in R
  # a_theta = (p_Y + p_M)
  a_theta = 1
  b_theta = 1
  # a_omega = (p_Z + p_C)
  a_omega = 1
  b_omega = 1
  
  # =============================================================================
  # 3. INITIAL VALUES
  # =============================================================================
  # Set initial values for cluster membership
  # Use K-means to find 5 whole data clusters
  # Set initial values for cluster membership
  Syx = matrix(nrow = N, ncol = 2, 1)
  
  # Use Ky-means to find Ky =5 whole data clusters
  Ky = 5
  Sy = kmeans(cbind(Tobs, Dobs, Zobs, scaXobs), Ky)$cluster
  Syx[,1] = Sy
  
  # Make Kx = 2 x-clusters within each y-cluster
  Kx = 2
  Sx = numeric(N)
  for (k in 1:Ky) {
    ind_Sy = (Sy == k)
    Sx[ind_Sy] = kmeans(cbind(Zobs, scaXobs)[ind_Sy,], Kx)$cluster
  }
  Syx[,2] = Sx
  
  # Determine unique clusters from the cluster membership variable, Syx
  unique_Syx = unique(Syx)
  
  # Sort clusters
  unique_Syx = unique_Syx[order(unique_Syx[,1], unique_Syx[,2]), , drop = FALSE]
  
  # Make vector of y-clusters
  unique_Sy = sort(unique(Sy))
  
  # Calculate the number of y-clusters
  Ky = length(unique_Sy)
  
  # Initialize alpha parameters
  alpha_theta = 2
  alpha_omega = 2
  
  # --- Initialize outcome parameters (theta) ---
  sig2Pars = t(rscainvchisq_cpp(Ky, a_sig2, b_sig2))
  betaPars = sapply(1:Ky, function(k) rmvn_cpp(1, a_beta, sig2Pars[k] * B_beta))
  
  # ... (Initialize latent log-times y, y*) ...
  Tstar = numeric(N)
  Tstar[Dobs == 1] = Tobs[Dobs == 1]
  for (k in 1:Ky) {
    ind_curr = which(Sy == unique_Sy[k])
    ind_cens = ind_curr[Dobs[ind_curr] == 0]
    
    if(length(ind_cens) > 0) {
      y_full = Tstar[ind_curr] 
      X_full = matX[ind_curr, , drop=FALSE]
      
      # Update Parameters (Collapsed Gibbs step)
      Binv_n  = Binv_beta + crossprod(X_full)
      Bn_beta = inv_cpp(Binv_n)
      an_beta = Bn_beta %*% (aBinv_beta + crossprod(X_full, y_full))
      
      residuals = y_full - X_full %*% an_beta
      beta_diff = an_beta - a_beta
      
      sse_data  = sum(residuals^2)
      sse_prior = as.numeric(t(beta_diff) %*% Binv_beta %*% beta_diff)
      
      an_sig2 = a_sig2 + length(y_full)
      bn_sig2 = b_sig2 + sse_data + sse_prior
      
      # Prepare Predictive Distribution
      X_cens = matX[ind_cens, , drop=FALSE]
      lower_vec = Tobs[ind_cens]
      
      quad_form = rowSums((X_cens %*% Bn_beta) * X_cens)
      scale_val = sqrt((bn_sig2 / an_sig2) * (1 + quad_form))
      mu_val    = as.vector(X_cens %*% an_beta)
      
      # Calculate Dynamic Upper Bound (Recommendation 2)
      # Use 99.9th percentile of the current predictive t-distribution
      t_quant   = qt(0.999, df = an_sig2)
      upper_vec = mu_val + (t_quant * scale_val)
      
      # Safety: Ensure upper > lower + buffer (4*scale) to prevent numerical trapping
      upper_safe = pmax(upper_vec, lower_vec + 4 * scale_val)
      
      # Sample
      # Note: 'upper' is now a vector, so it moves out of MoreArgs
      Tstar[ind_cens] = mapply(rtrunct_cpp, 
                               mu = mu_val, 
                               scale = scale_val, 
                               lower = lower_vec, 
                               upper = upper_safe, 
                               MoreArgs = list(n = 1, df = an_sig2))
    }
  }
  
  for (k in unique_Sy) {
    ind_S = which(Sy == k)
    n_temp = length(ind_S)
    T_temp = Tstar[ind_S]
    matX_temp = matrix(matX[ind_S,], ncol = p_matX)
    
    # update parameters for true failure time
    regression_results = update_reg_con(T_temp, matX_temp, betaPars[,k],
                                        a_beta, Binv_beta, aBinv_beta, a_sig2, b_sig2)
    sig2Pars[,k] = regression_results$sig2_par
    betaPars[,k] = regression_results$beta_par
  }
  
  # Update confounders parameters
  # Within each cluster parameters are stored in a long vector n_k x 1
  # Binary treatment have 1 parameter (etaPars)
  # Binary confounders have 1 parameter (piPars)
  # Continuous confounders have 2 parameters (muPars, tau2Pars)
  
  # Initialize vectors n_k is the total number of clusters (x and y clusters)
  Kyx = nrow(unique(Syx))
  
  # --- Initialize propensity parameters (eta) ---
  etaPars = matrix(nrow = 1, ncol = Kyx)
  # --- Initialize covariate parameters (omega) ---
  piPars = matrix(nrow = p_X1, ncol = Kyx) # Probabilities
  muPars = matrix(nrow = p_X2, ncol = Kyx) # Means
  tau2Pars = matrix(nrow = p_X2, ncol = Kyx) # Variances
  
  # update parameters for binary treatment
  {
    for (q in 1:p_Z) {
      count = 1
      if (p_Z==1) {
        z_temp_q = Zobs
      } else if (p_Z>1) {
        z_temp_q = Zobs[,q]
      }
      
      for (k in unique_Sy) {
        ind_Sy = (Sy == k)
        
        Sx_temp = Sx[ind_Sy]
        z_temp_qk = z_temp_q[ind_Sy]
        unique_Sx_temp = sort(unique(Sx_temp))
        for (r in unique_Sx_temp) {
          z_temp_qkr = z_temp_qk[Sx_temp == r]
          
          # update parameters for binary treatment
          bivariate_pars = update_conf_bin(z_temp_qkr, a_eta, b_eta)
          etaPars[,count] = bivariate_pars
          
          count = count + 1
        }
      }
    }
  }
  
  # update parameters for binary confounders
  if (p_X1 > 0) {
    for (q in 1:p_X1) {
      # beta(1,1) prior
      count = 1
      x_temp_q = scaXobs[,q]
      for (k in unique_Sy) {
        ind_Sy = (Sy == k)
        
        Sx_temp = Sx[ind_Sy]
        x_temp_qk = x_temp_q[ind_Sy]
        unique_Sx_temp = sort(unique(Sx_temp))
        for (r in unique_Sx_temp) {
          x_temp_qkr = x_temp_qk[Sx_temp == r]
          
          # posterior is beta
          bivariate_pars = update_conf_bin(x_temp_qkr, a_pi, b_pi)
          piPars[q, count] = bivariate_pars
          
          count = count + 1
        }
      }
    }
  }
  
  # update parameters for continuous confounders
  if (p_X2 > 0) {
    for (q in 1:p_X2) {
      # beta(1,1) prior
      count = 1
      x_temp_q = scaXobs[,(p_X1 + q)]
      for (k in unique_Sy) {
        ind_Sy = (Sy == k)
        
        Sx_temp = Sx[ind_Sy]
        x_temp_qk = x_temp_q[ind_Sy]
        unique_Sx_temp = sort(unique(Sx_temp))
        for (r in unique_Sx_temp) {
          x_temp_qkr = x_temp_qk[Sx_temp == r]
          
          # posterior for mu. prior for mu|sigma^2 mean 0 prior sample size2
          continuous_pars = update_conf_con(x_temp_qkr, a_tau2, b_tau2, a_mu, b_mu)
          tau2Pars[q, count] = continuous_pars$var_par
          muPars[q, count] = continuous_pars$mean_par
          
          count = count + 1
        }
      }
    }
  }
  
  # num_MC_prior = 1e5
  # Calculate f0_x and f0_y for use in cluster function in Gibbs Sampler
  # Average confounders distribution over prior for each x_i
  ab_tau = a_tau2 * b_tau2
  a_tau_new = (a_tau2 + 1)/2
  a_tau_half = a_tau2/2
  b_mu_ratio = b_mu/(b_mu+1)
  margin_part1 = (gamma(a_tau_new)/gamma(a_tau_half)) * sqrt(b_mu_ratio/pi) * (ab_tau)^{a_tau_half}
  
  p_Z = 1
  f0_zx_all = matrix(0, nrow = N, ncol = p_Z + p_X)
  # Binary treatment
  # Beta-Binomial
  for (q in 1:p_Z) {
    f0_zx_all[, q] = beta(a_eta + Zobs[q], b_eta - Zobs[q] + 1) # beta(a_pi,b_pi) = 1 in this case
  }
  # Binary confounders
  # Beta-Binomial
  if (p_X1>0){
    for (q in 1:p_X1) {
      f0_zx_all[, (p_Z + q)] = beta(a_pi + scaXobs[, q], b_pi - scaXobs[, q] + 1) # beta(a_pi,b_pi) = 1 in this case
    }
  }
  # Continuous confounders
  if (p_X2>0){
    for (q in 1:p_X2) {
      x_temp = scaXobs[, (p_X1 + q)]
      margin_part2 = (ab_tau + b_mu_ratio * (x_temp-a_mu)^{2})^(-a_tau_new)
      f0_zx_all[, (p_Z + p_X1 + q)] = margin_part1 * margin_part2
    }
  }
  
  # Take product (confounders are assumed to be locally independent).
  # Result is vector of size N
  f0_zx = apply(f0_zx_all, 1, prod)
  
  # =============================================================================
  # 4. GIBBS SAMPLING LOOP
  # =============================================================================
  gibbs_total = gibbs_iter + gibbs_burnin
  n_MCMC = floor(gibbs_iter / gibbs_thin)
  
  # Define iteration check
  iter_check = floor(gibbs_iter/20)
  
  # Make lists to store draws from Gibbs Sampler
  betaLists = list(NA)
  sig2Lists = list(NA)
  etaLists  = list(NA)
  piLists   = list(NA)
  muLists   = list(NA)
  tau2Lists = list(NA)
  
  n_kLists       = list(NA)
  n_rkLists      = list(NA)
  KyLists        = list(NA)
  KyxLists       = list(NA)
  max_Kx_SyLists = list(NA)
  
  alpha_theta_draws = numeric(n_MCMC)
  alpha_omega_draws = numeric(n_MCMC)
  
  count_it = 1
  # End initial values ------------------------------------------------------------
  
  # First draw each parameter for BNP model. Then calculate causal effect.
  for (gibbs_reps in 1:gibbs_total) {
    # gibbs_reps = 1
    # Update of all parameters in BNP model--------------------------------
    # Update cluster membership -------------------------------------------
    # sourceCpp("BNPsurvival_cpp.cpp")
    cluster_res = EDPMMcluster_cpp(N, p_X1, p_X2, 
                                   Tstar, matX, Zobs, scaXobs,
                                   Sy, Sx, unique_Syx, 
                                   alpha_theta, alpha_omega, 
                                   betaPars, sig2Pars, etaPars, 
                                   piPars, muPars, tau2Pars,
                                   a_beta, B_beta, Binv_beta, 
                                   a_sig2, b_sig2, 
                                   a_eta, b_eta,
                                   a_pi, b_pi, 
                                   a_mu, b_mu, 
                                   a_tau2, b_tau2,
                                   f0_zx)
    # Store cluster membership output from cluster function
    Syx = cluster_res$Syx
    Sy  = Syx[,1]
    Sx  = Syx[,2]
    
    betaPars = cluster_res$betaPars
    sig2Pars = cluster_res$sig2Pars
    etaPars  = cluster_res$etaPars
    piPars   = cluster_res$piPars
    muPars   = cluster_res$muPars
    tau2Pars = cluster_res$tau2Pars
    
    # Make matrix of clusters
    unique_Syx = cluster_res$unique_Syx
    
    # Make vector of y-clusters
    unique_Sy = cluster_res$unique_Sy
    
    # Calculate the number of clusters
    Kyx = cluster_res$Kyx
    
    # Calculate the number of y-clusters
    Ky = cluster_res$Ky
    
    # Find the largest number of x clusters
    max_Sx = cluster_res$max_Sx
    
    # Store the number of x-clusters in each y-cluster in a vector
    max_Kx_Sy = cluster_res$max_Kx_Sy
    
    # Calculate the number of subjects in each y-cluster and store in vector n_k.
    # Calculate the number of subjects in each x-cluster and store in matrix n_rk.
    # Use k to index y-clusters and r to index x-clusters.
    n_k = cluster_res$n_k
    n_rk = cluster_res$n_rk
    
    # End update of cluster membership --------------------------------------------
    
    # Update of cluster-specific parameters for confounders---------------------
    # ... (Initialize latent log-times y, y*) ...
    Tstar[Dobs == 1] = Tobs[Dobs == 1]
    for (k in 1:Ky) {
      ind_curr = which(Sy == unique_Sy[k])
      ind_cens = ind_curr[Dobs[ind_curr] == 0]
      
      if(length(ind_cens) > 0) {
        y_full = Tstar[ind_curr] 
        X_full = matX[ind_curr, , drop=FALSE]
        
        # Update Parameters (Collapsed Gibbs step)
        Binv_n  = Binv_beta + crossprod(X_full)
        Bn_beta = inv_cpp(Binv_n)
        an_beta = Bn_beta %*% (aBinv_beta + crossprod(X_full, y_full))
        
        residuals = y_full - X_full %*% an_beta
        beta_diff = an_beta - a_beta
        
        sse_data  = sum(residuals^2)
        sse_prior = as.numeric(t(beta_diff) %*% Binv_beta %*% beta_diff)
        
        an_sig2 = a_sig2 + length(y_full)
        bn_sig2 = b_sig2 + sse_data + sse_prior
        
        # Prepare Predictive Distribution
        X_cens = matX[ind_cens, , drop=FALSE]
        lower_vec = Tobs[ind_cens]
        
        quad_form = rowSums((X_cens %*% Bn_beta) * X_cens)
        scale_val = sqrt((bn_sig2 / an_sig2) * (1 + quad_form))
        mu_val    = as.vector(X_cens %*% an_beta)
        
        # Calculate Dynamic Upper Bound (Recommendation 2)
        # Use 99.9th percentile of the current predictive t-distribution
        t_quant   = qt(0.999, df = an_sig2)
        upper_vec = mu_val + (t_quant * scale_val)
        
        # Safety: Ensure upper > lower + buffer (4*scale) to prevent numerical trapping
        upper_safe = pmax(upper_vec, lower_vec + 4 * scale_val)
        
        # Sample
        # Note: 'upper' is now a vector, so it moves out of MoreArgs
        Tstar[ind_cens] = mapply(rtrunct_cpp, 
                                 mu = mu_val, 
                                 scale = scale_val, 
                                 lower = lower_vec, 
                                 upper = upper_safe, 
                                 MoreArgs = list(n = 1, df = an_sig2))
      }
    }
    
    for (k in unique_Sy) {
      ind_S = which(Sy == k)
      n_temp = length(ind_S)
      T_temp = Tstar[ind_S]
      matX_temp = matrix(matX[ind_S,], ncol = p_matX)
      
      # update parameters for true failure time
      regression_results = update_reg_con(T_temp, matX_temp, betaPars[,k],
                                          a_beta, Binv_beta, aBinv_beta, a_sig2, b_sig2)
      sig2Pars[,k] = regression_results$sig2_par
      betaPars[,k] = regression_results$beta_par
    }
    
    # Update confounders parameters
    # Within each cluster parameters are stored in a long vector n_k x 1
    # Binary treatment have 1 parameter (etaPars)
    # Binary confounders have 1 parameter (piPars)
    # Continuous confounders have 2 parameters (muPars, tau2Pars)
    {
      for (q in 1:p_Z) {
        count = 1
        if (p_Z==1) {
          z_temp_q = Zobs
        } else if (p_Z>1) {
          z_temp_q = Zobs[,q]
        }
        
        for (k in unique_Sy) {
          ind_Sy = (Sy == k)
          
          Sx_temp = Sx[ind_Sy]
          z_temp_qk = z_temp_q[ind_Sy]
          unique_Sx_temp = sort(unique(Sx_temp))
          for (r in unique_Sx_temp) {
            z_temp_qkr = z_temp_qk[Sx_temp == r]
            
            # update parameters for binary treatment
            bivariate_pars = update_conf_bin(z_temp_qkr, a_eta, b_eta)
            etaPars[,count] = bivariate_pars
            
            count = count + 1
          }
        }
      }
    }
    
    if (p_X1 > 0) {
      for (q in 1:p_X1) {
        # beta(1,1) prior
        count = 1
        x_temp_q = scaXobs[,q]
        for (k in unique_Sy) {
          ind_Sy = (Sy == k)
          
          Sx_temp = Sx[ind_Sy]
          x_temp_qk = x_temp_q[ind_Sy]
          unique_Sx_temp = sort(unique(Sx_temp))
          for (r in unique_Sx_temp) {
            x_temp_qkr = x_temp_qk[Sx_temp == r]
            
            # posterior is beta
            bivariate_pars = update_conf_bin(x_temp_qkr, a_pi, b_pi)
            piPars[q, count] = bivariate_pars
            
            count = count + 1
          }
        }
      }
    }
    
    if (p_X2 > 0) {
      for (q in 1:p_X2) {
        # beta(1,1) prior
        count = 1
        x_temp_q = scaXobs[,(p_X1 + q)]
        for (k in unique_Sy) {
          ind_Sy = (Sy == k)
          
          Sx_temp = Sx[ind_Sy]
          x_temp_qk = x_temp_q[ind_Sy]
          unique_Sx_temp = sort(unique(Sx_temp))
          for (r in unique_Sx_temp) {
            x_temp_qkr = x_temp_qk[Sx_temp == r]
            
            # posterior for mu. prior for mu|sigma^2 mean 0 prior sample size2
            continuous_pars = update_conf_con(x_temp_qkr, a_tau2, b_tau2, a_mu, b_mu)
            tau2Pars[q, count] = continuous_pars$var_par
            muPars[q, count] = continuous_pars$mean_par
            
            count = count + 1
          }
        }
      }
    }
    # End update of cluster-specific parameters for confounders---------------------
    
    # Update concentration parameters--------------------------------------------
    alpha_theta = update_alpha_theta_cpp(N, Ky, alpha_theta, a_theta, b_theta)
    alpha_omega = update_alpha_omega_cpp(Kyx, Ky, Sy, unique_Sy, alpha_omega, a_omega, b_omega)
    
    # End update of concentration parameters-----------------------------------
    # End update of all parameters in BNP model--------------------------------
    
    if (gibbs_reps < gibbs_burnin) {
    } else if (gibbs_reps == gibbs_burnin) {
      cat("Bur-In End",gibbs_reps,"Time:",date(),"\n")
    } else if (gibbs_reps > gibbs_burnin) {
      if (gibbs_reps %% gibbs_thin == 0) {
        betaLists[[count_it]] = betaPars
        sig2Lists[[count_it]] = sig2Pars
        etaLists[[count_it]]  = etaPars
        piLists[[count_it]]   = piPars
        muLists[[count_it]]   = muPars
        tau2Lists[[count_it]] = tau2Pars
        
        n_kLists[[count_it]]       = n_k
        n_rkLists[[count_it]]      = n_rk
        KyLists[[count_it]]        = Ky
        KyxLists[[count_it]]       = Kyx
        max_Kx_SyLists[[count_it]] = max_Kx_Sy
        
        alpha_theta_draws[count_it] = alpha_theta
        alpha_omega_draws[count_it] = alpha_omega
        
        count_it = count_it + 1
      }
      
      if (gibbs_reps %% iter_check == 0) {
        cat("Gibbs Iteration",(gibbs_reps-gibbs_burnin),"(",(gibbs_reps-gibbs_burnin)/gibbs_iter*100,"%)","Time:",date(),"\n")
      }
    }
  }
  
  # =============================================================================
  # 5. RETURN RESULTS
  # =============================================================================
  # constants
  constants = list(N = N, n_MCMC = n_MCMC,
                   p_matX = p_matX, p_X1 = p_X1, p_X2 = p_X2)
  
  # priors
  priors = list(a_sig2 = a_sig2, b_sig2 = b_sig2, 
                a_beta = a_beta, B_beta = B_beta,
                a_eta = a_eta, b_eta = b_eta,
                a_pi = a_pi, b_pi = b_pi, 
                a_mu = a_mu, b_mu = b_mu, 
                a_tau2 = a_tau2, b_tau2 = b_tau2)
  
  # # Precalculated values for Post-Processing steps
  # MCpriors = list(f0_z0_mc = f0_z0_mc, 
  #                 f0_z1_mc = f0_z1_mc,
  #                 TbetaPrior_mc = TbetaPrior_mc, 
  #                 MbetaPrior_mc = MbetaPrior_mc, 
  #                 Msig2Prior_mc = Msig2Prior_mc,
  #                 VbetaPrior_mc = VbetaPrior_mc, 
  #                 Vsig2Prior_mc = Vsig2Prior_mc)
  
  # MCMC Posteriors
  MCMCposteriors = list(betaLists   = betaLists,
                        sig2Lists   = sig2Lists,
                        etaLists    = etaLists,
                        piLists     = piLists,
                        muLists     = muLists,
                        tau2Lists   = tau2Lists,
                        n_kLists   = n_kLists,
                        n_rkLists  = n_rkLists,
                        KyLists    = KyLists,
                        KyxLists   = KyxLists,
                        max_Kx_SyLists    = max_Kx_SyLists,
                        alpha_theta_draws = alpha_theta_draws,
                        alpha_omega_draws = alpha_omega_draws)
  
  # MCMC results
  MCMCresult = list(priors = priors,
                    # MCpriors = MCpriors,
                    constants = constants,
                    MCMCposteriors = MCMCposteriors)
  # gibbs_reps
  # unique_Syx
  # gibbs_reps
  # n_k
  # n_rk
  # Ky
  # Kyx
  # max_Kx_Sy
  # 
  # {
  #   plot(unlist(KyLists), type="s")
  # 
  #   # Setup the base plot with the first density
  #   plot(density(Tstar[Tstar<10]),
  #        xlim = c(-12, 15),
  #        ylim = c(0, 0.2),
  #        main = "Density Comparisons",
  #        xlab = "Log Time",
  #        col = "black",
  #        lty = 1,
  #        lwd = 1);
  # 
  #   # Add the remaining densities using lines
  #   lines(density(temp_data$logY),
  #         col = "black",
  #         lty = 1,
  #         lwd = 3);
  # 
  #   lines(density(Tstar[Zobs == 0]),
  #         col = "blue",
  #         lty = 2,
  #         lwd = 1);
  # 
  #   lines(density(temp_data$logY0),
  #         col = "blue",
  #         lty = 2,
  #         lwd = 3);
  # 
  #   lines(density(Tstar[Zobs == 1]),
  #         col = "red",
  #         lty = 3,
  #         lwd = 1);
  # 
  #   lines(density(temp_data$logY1),
  #         col = "red",
  #         lty = 3,
  #         lwd = 3);
  # 
  #   lines(density(temp_data$logT[Zobs==0]),
  #         col = "lightgreen",
  #         lty = 4,
  #         lwd = 1);
  # 
  #   lines(density(temp_data$logT[Zobs==1]),
  #         col = "lightgreen",
  #         lty = 4,
  #         lwd = 3);
  # 
  #   # Add a legend to distinguish the lines
  #   legend("topright",
  #          legend = c("Tstar", "logY", "Tstar[Z==0]", "logY0", "Tstar[Z==1]", "logY1", "logT[Z==0]", "logT[Z==1]"),
  #          col = c("black", "black", "blue", "blue", "red", "red", "lightgreen", "lightgreen"),
  #          lty = c(1, 1, 2, 2, 3, 3, 4, 4),
  #          lwd = c(1, 3, 1, 3, 1, 3, 1, 3),
  #          cex = 0.8);
  # }
  
  return(MCMCresult)
}

# -----------------------------------------------------------------------------
# Define function for the post-processing under the BNP (EDPMM) model
SURVIVAL_EDPMM_POST = function(object, nu, rho, M, esttype = NULL){
  # -----------------------------------------------------------------------------
  # object = EDPMM_results_MCMC
  
  # -----------------------------------------------------------------------------
  # Define LogSumExp helper for mixture calculations
  logSumExp = function(vals) {
    max_val = max(vals)
    if(max_val == -Inf) return(-Inf)
    max_val + log(sum(exp(vals - max_val)))
  }
  
  normalize_log_probs = function(log_prob_matrix) {
    max_log_probs = apply(log_prob_matrix, 1, max)
    log_probs_shifted = log_prob_matrix - max_log_probs
    probs_shifted = exp(log_probs_shifted)
    row_sums = rowSums(probs_shifted)
    
    # Handle rows that are all -Inf (row_sum = 0)
    zero_rows = which(row_sums == 0)
    if(length(zero_rows) > 0) {
      probs_shifted[zero_rows,] = 1 / ncol(log_prob_matrix) # Assign uniform prob
      row_sums[zero_rows] = 1
    }
    return(probs_shifted / row_sums)
  }
  
  # -----------------------------------------------------------------------------
  # Extract Constants and Settings
  # Ensure nu and rho are vectors
  nu_vec = as.vector(nu)
  rho_vec = as.vector(rho)
  n_nu = length(nu_vec)
  n_rho = length(rho_vec)
  
  level = 0.05
  quantile_alpha = c(level/2,1-level/2)
  z0 = 0
  z1 = 1
  
  N = object$constants$N
  num_MC = object$constants$num_MC
  n_MCMC = object$constants$n_MCMC
  p_matX = object$constants$p_matX
  p_X1 = object$constants$p_X1
  p_X2 = object$constants$p_X2
  
  a_beta = object$priors$a_beta
  B_beta = object$priors$B_beta
  a_sig2 = object$priors$a_sig2
  b_sig2 = object$priors$b_sig2
  a_eta = object$priors$a_eta
  b_eta = object$priors$b_eta
  a_pi = object$priors$a_pi
  b_pi = object$priors$b_pi
  a_mu = object$priors$a_mu
  b_mu = object$priors$b_mu
  a_tau2 = object$priors$a_tau2
  b_tau2 = object$priors$b_tau2
  
  # -----------------------------------------------------------------------------
  # Define interation check
  iter_check = floor(n_MCMC/10)
  
  Y_rho_nu_z0x_store = array(NA, dim = c(n_MCMC, n_nu, n_rho))
  Y_rho_nu_z1x_store = array(NA, dim = c(n_MCMC, n_nu, n_rho))
  for (post_reps in 1:n_MCMC) {
    # post_reps = 1
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    alpha_omega = object$MCMCposteriors$alpha_omega_draws[post_reps]
    alpha_theta = object$MCMCposteriors$alpha_theta_draws[post_reps]
    Ky = object$MCMCposteriors$KyLists[[post_reps]]
    n_k = object$MCMCposteriors$n_kLists[[post_reps]]
    n_rk = object$MCMCposteriors$n_rkLists[[post_reps]]
    max_Kx_Sy = object$MCMCposteriors$max_Kx_SyLists[[post_reps]]
    betaPars = object$MCMCposteriors$betaLists[[post_reps]]
    sig2Pars = object$MCMCposteriors$sig2Lists[[post_reps]]
    etaPars = object$MCMCposteriors$etaLists[[post_reps]]
    piPars = object$MCMCposteriors$piLists[[post_reps]]
    muPars = object$MCMCposteriors$muLists[[post_reps]]
    tau2Pars = object$MCMCposteriors$tau2Lists[[post_reps]]
    
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Draw new clusters for x
    lambda_y_mc = n_k # c(n_k, alpha)
    S_y_mc = sapply(1:M, function(l) rmultinom_cpp(lambda_y_mc))
    unique_S_y_mc = sort(unique(S_y_mc))
    Ky_mc = max(unique_S_y_mc)
    
    # -----------------------------------------------------------------------------
    # Draw confounders for each of the M observations
    S_x_mc = numeric(M)
    n_k_mc  = numeric(Ky)
    n_rk_mc = matrix(0, nrow = nrow(n_rk), ncol = ncol(n_rk))
    max_Kx_Sy_mc = numeric(Ky)
    
    X_mc = matrix(nrow = M, ncol = (p_X1 + p_X2))
    count = 1
    for (k in 1:Ky) {
      ind_S_y_mc = which(S_y_mc == k)
      n_k_temp = length(ind_S_y_mc)
      n_k_mc[k] = n_k_temp
      
      if (n_k_temp > 0) {
        max_Kx_Sy_temp = max_Kx_Sy[k]
        weights_count = n_rk[k, 1:max_Kx_Sy_temp]
        if (max_Kx_Sy_temp == 1) {
          S_x_sub = rep(1, n_k_temp)
        } else {
          S_x_sub = sample(1:max_Kx_Sy_temp, n_k_temp, replace = TRUE, prob = n_rk[k,1:max_Kx_Sy_temp])
        }
        S_x_mc[ind_S_y_mc] = S_x_sub
        max_Kx_Sy_mc[k] = max(S_x_sub)
        
        for (r in 1:max_Kx_Sy_temp) {
          ind_Sy_x_mc = ind_S_y_mc[which(S_x_sub == r)]
          num_Sy_x_mc = length(ind_Sy_x_mc)
          
          n_rk_mc[k, r] = num_Sy_x_mc
          if (num_Sy_x_mc > 0) {
            if (p_X1 > 0) {
              for (q in 1:p_X1) {
                X_mc[ind_Sy_x_mc, q] = rbinom(num_Sy_x_mc, 1, piPars[q, count])
              }
            }
            if (p_X2 > 0) {
              for (q in 1:p_X2) {
                X_mc[ind_Sy_x_mc, (p_X1 + q)] = rnorm(num_Sy_x_mc, muPars[q, count], sqrt(tau2Pars[q, count]))
              }
            }
          } # end if observations exist
          count = count + 1
        } # end count loop
      } # end if n_k_temp > 0
    } # end k loop
    
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # sum_max_Kx_Sy_mc = sum(max_Kx_Sy_mc)
    # lambda_Y_z0x_mc = matrix(nrow = M, ncol = sum_max_Kx_Sy_mc)
    # lambda_Y_z1x_mc = matrix(nrow = M, ncol = sum_max_Kx_Sy_mc)
    lambda_Y_z0x_mc = matrix(nrow = M, ncol = Ky)
    lambda_Y_z1x_mc = matrix(nrow = M, ncol = Ky)
    count = 1
    for (k in 1:Ky) {
      n_k_temp = n_k[k]
      max_Kx_Sy_temp = max_Kx_Sy[k]
      
      log_w_k_temp = log(n_k_temp/N)
      
      w_rk_temp = n_rk[k, 1:max_Kx_Sy_temp]
      log_w_rk_temp = log(w_rk_temp / sum(w_rk_temp))
      
      log_comps_z0 = matrix(-Inf, nrow = M, ncol = max_Kx_Sy_temp)
      log_comps_z1 = matrix(-Inf, nrow = M, ncol = max_Kx_Sy_temp)
      for (r in 1:max_Kx_Sy_temp) {
        prob_x_mc = numeric(M)
        prob_z0_mc = dbinom(z0, 1, etaPars[, count], log = TRUE)
        prob_z1_mc = dbinom(z1, 1, etaPars[, count], log = TRUE)
        if (p_X1 > 0) {
          for (q in 1:p_X1) {
            prob_x_mc = prob_x_mc + dbinom(X_mc[, q], 1, piPars[q, count], log = TRUE)
          }
        }
        if (p_X2 > 0) {
          for (q in 1:p_X2) {
            prob_x_mc = prob_x_mc + dnorm(X_mc[, (p_X1 + q)], muPars[q, count], sqrt(tau2Pars[q, count]), log = TRUE)
          }
        }
        log_comps_z0[, r] = log_w_rk_temp[r] + prob_z0_mc + prob_x_mc
        log_comps_z1[, r] = log_w_rk_temp[r] + prob_z1_mc + prob_x_mc
        count = count + 1
      }
      lambda_Y_z0x_mc[,k] = apply(log_comps_z0, 1, logSumExp) + log_w_k_temp
      lambda_Y_z1x_mc[,k] = apply(log_comps_z1, 1, logSumExp) + log_w_k_temp
    }
    lambda_Y_z0x_mc = normalize_log_probs(lambda_Y_z0x_mc)
    lambda_Y_z1x_mc = normalize_log_probs(lambda_Y_z1x_mc)
    
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Create the design matrices for all M subjects
    matX_z0x_mc = cbind(1, z0, X_mc) # Intercept, Zobs=0, Covariates
    matX_z1x_mc = cbind(1, z1, X_mc) # Intercept, Zobs=1, Covariates
    MUy_z0x_mc = (matX_z0x_mc %*% betaPars)
    MUy_z1x_mc = (matX_z1x_mc %*% betaPars)
    Ey_z0x_mc = rowSums(MUy_z0x_mc * lambda_Y_z0x_mc)
    Ey_z1x_mc = rowSums(MUy_z1x_mc * lambda_Y_z1x_mc)
    SDy_z0x_mc = sqrt(rowSums(MUy_z0x_mc^{2} * lambda_Y_z0x_mc^{2}))
    SDy_z1x_mc = sqrt(rowSums(MUy_z1x_mc^{2} * lambda_Y_z1x_mc^{2}))
    
    # -----------------------------------------------------------------------------
    lambda_Y_z0x = apply(lambda_Y_z0x_mc, 2, mean)
    lambda_Y_z1x = apply(lambda_Y_z1x_mc, 2, mean)
    MUy_z0x = apply(MUy_z0x_mc, 2, mean)
    MUy_z1x = apply(MUy_z1x_mc, 2, mean)
    Ey_z0x = mean(Ey_z0x_mc)
    Ey_z1x = mean(Ey_z1x_mc)
    SDy_z0x = mean(SDy_z0x_mc)
    SDy_z1x = mean(SDy_z1x_mc)
    
    # LOOP OVER nu
    for (ii_nu in 1:n_nu) {
      curr_nu = nu_vec[ii_nu]
      curr_log_nu = log(curr_nu)
      
      S_nu_z0x_vec = sapply(1:M, function(m) 
        S_DPMM_cpp(curr_log_nu, lambda_Y_z0x_mc[m,], 
                   MUy_z0x_mc[m,], sqrt(sig2Pars), TRUE))
      S_nu_z1x_vec = sapply(1:M, function(m) 
        S_DPMM_cpp(curr_log_nu, lambda_Y_z1x_mc[m,], 
                   MUy_z1x_mc[m,], sqrt(sig2Pars), TRUE))
      S_nu_z0x = mean(S_nu_z0x_vec)
      S_nu_z1x = mean(S_nu_z1x_vec)
      
      # LOOP OVER rho
      for(ii_rho in 1:n_rho) {
        curr_rho = rho_vec[ii_rho]
        
        target_rho_nu_z0x = S_nu_z0x + log(1 - curr_rho)
        target_rho_nu_z1x = S_nu_z1x + log(1 - curr_rho)
        
        # -----------------------------------------------------------------------------
        Y_rho_nu_z0x = S_optim_cpp(target_rho_nu_z0x,
                                   curr_log_nu,
                                   Ey_z0x,
                                   Ey_z0x - 6*SDy_z0x,
                                   Ey_z0x + 6*SDy_z0x,
                                   lambda_Y_z0x,
                                   MUy_z0x,
                                   sqrt(sig2Pars),
                                   TRUE)$optimizer
        Y_rho_nu_z1x = S_optim_cpp(target_rho_nu_z1x,
                                   curr_log_nu,
                                   Ey_z1x,
                                   Ey_z1x - 6*SDy_z1x,
                                   Ey_z1x + 6*SDy_z1x,
                                   lambda_Y_z1x,
                                   MUy_z1x,
                                   sqrt(sig2Pars),
                                   TRUE)$optimizer
        
        # # Stable Log-Survival Function in R
        # calc_log_survival <- function(y, lambda, mu, sig) {
        #   z <- (y - mu) / sig
        #   log_S_k <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
        #   log_terms <- log(lambda) + log_S_k
        #   max_val <- max(log_terms)
        #   result <- max_val + log(sum(exp(log_terms - max_val)))
        #   return(result)
        # }
        # obj_fun_z0 <- function(y_guess) {
        #   val <- calc_log_survival(y_guess, lambda_Y_z0x, MUy_z0x, sqrt(sig2Pars))
        #   return((val - target_rho_nu_z0x)^2)
        # }
        # obj_fun_z1 <- function(y_guess) {
        #   val <- calc_log_survival(y_guess, lambda_Y_z1x, MUy_z1x, sqrt(sig2Pars))
        #   return((val - target_rho_nu_z1x)^2)
        # }
        # 
        # # ---------------------------------------------------------
        # # OPTIMIZATION FOR Z0 (Control)
        # # ---------------------------------------------------------
        # opt_z0 <- optim(par = Ey_z0x,
        #                 fn = obj_fun_z0,
        #                 method = "L-BFGS-B",
        #                 lower = curr_log_nu + 1e-5, 
        #                 upper = Inf)
        # Y_rho_nu_z0x <- opt_z0$par
        # opt_z1 <- optim(par = Ey_z1x,
        #                 fn = obj_fun_z1,
        #                 method = "L-BFGS-B",
        #                 lower = curr_log_nu + 1e-5,
        #                 upper = Inf)
        # Y_rho_nu_z1x <- opt_z1$par
        
        # ---------------------------------------------------------
        # Store Results (Same as before)
        # ---------------------------------------------------------
        Y_rho_nu_z0x_store[post_reps, ii_nu, ii_rho] = exp(Y_rho_nu_z0x) - curr_nu
        Y_rho_nu_z1x_store[post_reps, ii_nu, ii_rho] = exp(Y_rho_nu_z1x) - curr_nu
        # -----------------------------------------------------------------------------
        Y_rho_nu_z0x_store[post_reps, ii_nu, ii_rho] = exp(Y_rho_nu_z0x) - curr_nu
        Y_rho_nu_z1x_store[post_reps, ii_nu, ii_rho] = exp(Y_rho_nu_z1x) - curr_nu
      } # End rho loop
    } # End nu loop
    # -----------------------------------------------------------------------------
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  E_rho_nu_diff_store = Y_rho_nu_z1x_store - Y_rho_nu_z0x_store
  
  z0_tables_list   = list()
  z1_tables_list   = list()
  diff_tables_list = list()
  counter = 1
  for(ii_nu in 1:n_nu) {
    for(ii_rho in 1:n_rho) {
      val_nu  = nu_vec[ii_nu]
      val_rho = rho_vec[ii_rho]
      
      draws_z0   = Y_rho_nu_z0x_store[, ii_nu, ii_rho]
      draws_z1   = Y_rho_nu_z1x_store[, ii_nu, ii_rho]
      draws_diff = E_rho_nu_diff_store[, ii_nu, ii_rho]
      
      summ_z0   = POSTsummary(draws_z0, esttype, quantile_alpha)
      summ_z1   = POSTsummary(draws_z1, esttype, quantile_alpha)
      summ_diff = POSTsummary(draws_diff, esttype, quantile_alpha)
      
      summ_z0   = cbind(nu = val_nu, rho = val_rho, summ_z0)
      summ_z1   = cbind(nu = val_nu, rho = val_rho, summ_z1)
      summ_diff = cbind(nu = val_nu, rho = val_rho, summ_diff)
      
      z0_tables_list[[counter]]   = summ_z0
      z1_tables_list[[counter]]   = summ_z1
      diff_tables_list[[counter]] = summ_diff
      
      counter = counter + 1
    }
  }
  
  Y_rho_nu_z0x_result  = do.call(rbind, z0_tables_list)
  Y_rho_nu_z1x_result  = do.call(rbind, z1_tables_list)
  E_rho_nu_diff_result = do.call(rbind, diff_tables_list)
  Y_rho_nu_z0x_result
  Y_rho_nu_z1x_result
  E_rho_nu_diff_result
  
  POSTresult = list(
    Y_rho_nu_z0x_result = Y_rho_nu_z0x_result, 
    Y_rho_nu_z1x_result = Y_rho_nu_z1x_result, 
    E_rho_nu_diff_result = E_rho_nu_diff_result,
    Y_rho_nu_z0x_store = Y_rho_nu_z0x_store,
    Y_rho_nu_z1x_store = Y_rho_nu_z1x_store,
    E_rho_nu_diff_store = E_rho_nu_diff_store
  )
  
  return(POSTresult)
}

# ------------------------------------------------------------------------
SURVIVAL_EDPMM = function(object, nu = 0, rho = 0.5, M = 2e3, 
                          gibbs_iter = 2e4, gibbs_burnin = 2e4, gibbs_thin = 1e2, 
                          esttype = "mean", saveall = FALSE){
  # object = temp_data
  MCMCresult = SURVIVAL_EDPMM_MCMC(object, gibbs_iter, gibbs_burnin, gibbs_thin)
  POSTresult = SURVIVAL_EDPMM_POST(MCMCresult, nu, rho, M, esttype)
  if (isTRUE(saveall)){
    result = c(MCMCresult, POSTresult)
  } else {
    result = POSTresult
  }
  return(result)
}

# ------------------------------------------------------------------------
# End function definitions -----------------------------------------------
# ------------------------------------------------------------------------