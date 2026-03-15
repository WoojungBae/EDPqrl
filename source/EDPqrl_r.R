# ------------------------------------------------------------------------
# Define functions -------------------------------------------------------
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
generate_data = function(Scn, n, cstar = NULL) {
  # Covariates: 2 binary and 3 continous covariates
  X1 = rbinom(n, 1, 0.5)
  X2 = rbinom(n, 1, 0.4 + 0.2 * X1)
  X3 = rnorm(n, 0, 1)
  X4 = rnorm(n, - 0.1 + 0.2 * X1 - 0.15 * X3, 1)
  X5 = rnorm(n,   0.1 - 0.2 * X2 + 0.15 * X4, 0.5)
  # mean: 0.5, 0.5, 0, 0, 0
  X = cbind(X1, X2, X3, X4, X5)
  
  rowSdX = apply(X, 1, sd)
  rowMeanX = apply(X, 1, mean)
  
  # treatment
  z0 = 0
  z1 = 1
  Z = rbinom(n, 1, pnorm(0.2 + 0.1 * X1 + 0.2 * X3 - 0.1 * X5))
  # mean: 0.6
  
  # design matrix
  matXz0 = cbind(1, z0, X)
  matXz1 = cbind(1, z1, X)
  beta1 = c( 0.3,  0.2, -0.3, -0.5,  0.6, -0.5, -0.3)
  beta2 = c( 2.1,  0.6, -0.5, -0.3,  0.2, -0.3, -0.5)
  mean1z0 = c(matXz0 %*% beta1)
  mean1z1 = c(matXz1 %*% beta1)
  mean2z0 = c(matXz0 %*% beta2)
  mean2z1 = c(matXz1 %*% beta2)
  sd1 = 0.30
  sd2 = 0.40
  df1 = 10
  df2 = 15
  
  betaC = c( 0, 0.2, -0.1, 0.1, -0.2,  0.1,  -0.2)
  meanCz0 = c(matXz0 %*% betaC)
  meanCz1 = c(matXz1 %*% betaC)
  # sdCz0 = 0.5
  # sdCz1 = 0.5
  sdCz0 = 2
  sdCz1 = 2
  
  # TRUE censoring
  # # Scn1 => cstar = 3.20 # for 20% censoring
  # # Scn1 => cstar = 2.80 # for 25% censoring
  # # Scn1 => cstar = 1.79 # for 40% censoring
  # # Scn1 => cstar = 1.15 # for 50% censoring
  # # Scn1 => cstar = 0.53 # for 60% censoring
  # # Scn1 => cstar = -0.52 # for 75% censoring
  # # Scn1 => cstar = -0.95 # for 80% censoring
  if (is.null(cstar)) {
    if (Scn == 1) {
      cstar = 3.20
    }
    if (Scn == 2) {
      cstar = 1.79
    }
    if(Scn == 3) {
      cstar = 0.53
    }
    if(Scn == 4) {
      cstar = -0.95
    }
  }
  logC0 = rnorm(n, cstar + meanCz0, sdCz0)
  logC1 = rnorm(n, cstar + meanCz1, sdCz1)
  
  # mixtur weights
  p1 = 0.4
  p2 = 0.6
  cluster = sapply(1:n, function(l) sample(2, 1, F, c(p1, p2)))
  
  # TRUE failure time
  # Scn1 - Scn4 => mixture of truncated normal and truncated t
  logY0 = numeric(n)
  logY1 = numeric(n)
  # # Scn1 => cstar = 3.20 # for 20% censoring
  # # Scn1 => cstar = 2.80 # for 25% censoring
  # # Scn1 => cstar = 1.79 # for 40% censoring
  # # Scn1 => cstar = 1.15 # for 50% censoring
  # # Scn1 => cstar = 0.53 # for 60% censoring
  # # Scn1 => cstar = -0.52 # for 75% censoring
  # # Scn1 => cstar = -0.95 # for 80% censoring
  # logY0 <- ifelse(cluster == 1,
  #                 mean1z0 + sd1 * rt(n, df = df1),
  #                 mean2z0 + sd2 * rt(n, df = df2))
  # logY1 <- ifelse(cluster == 1,
  #                 mean1z1 + sd1 * rt(n, df = df1),
  #                 mean2z1 + sd2 * rt(n, df = df2))
  for (i in 1:n) {
    if (cluster[i] == 1){
      logY0[i] = c(mean1z0[i] + sd1 * rt(1, df = df1))
      logY1[i] = c(mean1z1[i] + sd1 * rt(1, df = df1))
    } else if (cluster[i] == 2) {
      # logY0[i] = c(rtruncnorm_cpp(1, mean2z0[i], sd2, -5, 5))
      # logY1[i] = c(rtruncnorm_cpp(1, mean2z1[i], sd2, -5, 5))
      logY0[i] = c(rnorm(1, mean2z0[i], sd2))
      logY1[i] = c(rnorm(1, mean2z1[i], sd2))
    }
  }
  
  # Scn1 => cstar = 3.20 # for 20% censoring
  # Scn1 => cstar = 2.81 # for 25% censoring
  # Scn1 => cstar = 1.79 # for 40% censoring
  # Scn1 => cstar = 1.16 # for 50% censoring
  # Scn1 => cstar = 0.54 # for 60% censoring
  # Scn1 => cstar = -0.53 # for 75% censoring
  # Scn1 => cstar = -0.94 # for 80% censoring
  # for (i in 1:n) {
  #   if (cluster[i] == 1){
  #     logY0[i] = c(rtruncnorm_cpp(1, mean1z0[i], sd1, -10, 3))
  #     logY1[i] = c(rtruncnorm_cpp(1, mean1z1[i], sd1, -10, 3))
  #   } else if (cluster[i] == 2) {
  #     logY0[i] = c(rtrunct_cpp(1, mean2z0[i], sd2, 5, -5, 7))
  #     logY1[i] = c(rtrunct_cpp(1, mean2z1[i], sd2, 5, -5, 7))
  #   }
  # }
  # # Scn1 => cstar = 2.43 # for 20% censoring
  # # Scn1 => cstar = 2.22 # for 25% censoring
  # # Scn1 => cstar = 1.66 # for 40% censoring
  # # Scn1 => cstar = 1.29 # for 50% censoring
  # # Scn1 => cstar = 0.89 # for 60% censoring
  # # Scn1 => cstar = 0.12 # for 75% censoring
  # # Scn1 => cstar = -0.19 # for 80% censoring
  # logY0 <- ifelse(cluster == 1,
  #                 mean1z0 + sd1 * rt(n, df = df1),
  #                 mean2z0 + sd2 * rt(n, df = df2))
  # logY1 <- ifelse(cluster == 1,
  #                 mean1z1 + sd1 * rt(n, df = df1),
  #                 mean2z1 + sd2 * rt(n, df = df2))
  
  # # Scn = 1 & cstar = 1.57 => 25% censoring
  # # Scn = 1 & cstar = 0.47 => 50% censoring
  # # Scn = 1 & cstar = 0.03 => 60% censoring
  # # Scn = 1 & cstar = -0.68 => 75% censoring
  # # Scn = 1 & cstar = -0.95 => 80% censoring
  # # TRUE censoring
  # # Scn1 => cstar = 1.81 # for 20% censoring
  # # Scn2 => cstar = 0.92 # for 40% censoring
  # # Scn3 => cstar = 1.79 # for 20% censoring
  # # Scn4 => cstar = 0.91 # for 40% censoring
  # if (Scn == 1) {
  #   cstar = 1.81
  # }
  # if (Scn == 2) {
  #   cstar = 0.92
  # }
  # if(Scn == 3) {
  #   cstar = 1.79
  # }
  # if(Scn == 4) {
  #   cstar = 0.91
  # }
  # logC0 = rnorm(n, cstar + meanCz0, sdCz0)
  # logC1 = rnorm(n, cstar + meanCz1, sdCz1)
  # 
  # # mixtur weights
  # p1 = 0.4
  # p2 = 0.6
  # cluster = sapply(1:n, function(l) sample(2, 1, F, c(p1, p2)))
  # 
  # # TRUE failure time
  # # Scn1 & Scn2 => mixture of truncated normal and truncated t
  # if (Scn == 1 || Scn == 2){
  #   logY0 = numeric(n)
  #   logY1 = numeric(n)
  #   for (i in 1:n) {
  #     if (cluster[i] == 1){
  #       logY0[i] = c(rtruncnorm_cpp(1, mean1z0[i], sd1, -10, 3))
  #       logY1[i] = c(rtruncnorm_cpp(1, mean1z1[i], sd1, -10, 3))
  #     } else if (cluster[i] == 2) {
  #       logY0[i] = c(rtrunct_cpp(1, mean2z0[i], sd2, 5, -5, 7))
  #       logY1[i] = c(rtrunct_cpp(1, mean2z1[i], sd2, 5, -5, 7))
  #     }
  #   }
  # }
  # # Scn3 & Scn4 => mixture of normals
  # if (Scn == 3 || Scn == 4){
  #   logY0 = numeric(n)
  #   logY1 = numeric(n)
  #   for (i in 1:n) {
  #     if (cluster[i] == 1){
  #       logY0[i] = rnorm(1, mean1z0[i], sd1)
  #       logY1[i] = rnorm(1, mean1z1[i], sd1)
  #     } else if (cluster[i] == 2) {
  #       logY0[i] = rnorm(1, mean2z0[i], sd2)
  #       logY1[i] = rnorm(1, mean2z1[i], sd2)
  #     }
  #   }
  # }
  
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
  Ycon = 1e3
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
# update_reg_con = function(Y_temp, matX_temp, beta_curr,
#                           a_beta, Binv_beta, aBinv_beta,
#                           a_sig2, b_sig2) {
#   n_temp = length(Y_temp)
#   p_temp = length(beta_curr)
#   
#   if (n_temp == 0) {
#     sig2_prior = c(rscainvchisq_cpp(1, a_sig2, b_sig2))
#     sig2_prior = max(sig2_prior, 1e-5)
#     beta_prior = c(rmvn_cpp(1, a_beta, sig2_prior * inv_cpp(Binv_beta)))
#     return(list(beta_par = beta_prior, sig2_par = sig2_prior))
#   }
#   
#   matX_temp = matrix(matX_temp, ncol = p_temp)
#   matXtX_temp = t(matX_temp) %*% matX_temp
#   matXtY_temp = t(matX_temp) %*% Y_temp
#   
#   # Update the outcome regression parameters (sig2, beta)
#   beta_diff = beta_curr - a_beta
#   SSE_prior = as.numeric(t(beta_diff) %*% Binv_beta %*% beta_diff)
#   residuals = Y_temp - matX_temp %*% beta_curr
#   SSE_data = sum(residuals^2)
#   
#   a_sig2_new = a_sig2 + n_temp + p_temp 
#   b_sig2_new = (a_sig2 * b_sig2 + SSE_data + SSE_prior) / a_sig2_new
#   sig2_prop = c(rscainvchisq_cpp(1, a_sig2_new, b_sig2_new))
#   if(sig2_prop < 1e-5) {
#     sig2_prop = 1e-5
#   }
#   
#   Binv_beta_new = (matXtX_temp + Binv_beta) / sig2_prop
#   B_beta_new = tryCatch({
#     inv_cpp(Binv_beta_new)
#   }, error = function(e) {
#     diag(1 / diag(Binv_beta_new))
#   })
#   a_beta_new = B_beta_new %*% ((matXtY_temp + aBinv_beta) / sig2_prop)
#   beta_prop = c(rmvn_cpp(1, a_beta_new, B_beta_new))
#   
#   return(list(beta_par = beta_prop, sig2_par = sig2_prop))
# }

update_reg_con = function(Y_temp, matX_temp, a_beta, Binv_beta, a_sig2, b_sig2) {
  n_temp = length(Y_temp)
  p_temp = length(a_beta)
  
  if (n_temp == 0) {
    sig2_prior = c(rscainvchisq_cpp(1, a_sig2, b_sig2))
    sig2_prior = max(sig2_prior, 1e-5)
    beta_prior = c(rmvn_cpp(1, a_beta, sig2_prior * inv_cpp(Binv_beta)))
    return(list(beta_par = beta_prior, sig2_par = sig2_prior))
  }
  
  matX_temp = matrix(matX_temp, ncol = p_temp)
  matXtX_temp = t(matX_temp) %*% matX_temp
  matXtY_temp = t(matX_temp) %*% Y_temp
  
  # 1. Update Precision and Covariance (Independent of sigma^2)
  Binv_beta_new = Binv_beta + matXtX_temp
  B_beta_new = tryCatch({ inv_cpp(Binv_beta_new) }, error = function(e) diag(1 / diag(Binv_beta_new)))
  
  # 2. Update Mean vector
  a_beta_new = B_beta_new %*% (Binv_beta %*% a_beta + matXtY_temp)
  
  # 3. Update Variance shape (Marginal update: only add n_temp, NOT p_temp)
  a_sig2_new = a_sig2 + n_temp
  
  # 4. Update Variance scale (Analytical Marginal Sum of Squares)
  prior_quad = as.numeric(t(a_beta) %*% Binv_beta %*% a_beta)
  post_quad = as.numeric(t(a_beta_new) %*% Binv_beta_new %*% a_beta_new)
  yy = sum(Y_temp^2)
  
  sum_sq = a_sig2 * b_sig2 + yy + prior_quad - post_quad
  if (sum_sq < 1e-10) sum_sq = 1e-10 # Numerical safety
  b_sig2_new = sum_sq / a_sig2_new
  
  # 5. Block Draw: Draw sigma^2 first, then draw beta given sigma^2
  sig2_prop = c(rscainvchisq_cpp(1, a_sig2_new, b_sig2_new))
  if(sig2_prop < 1e-5) sig2_prop = 1e-5
  
  beta_prop = c(rmvn_cpp(1, a_beta_new, sig2_prop * B_beta_new))
  
  return(list(beta_par = beta_prop, sig2_par = sig2_prop))
}

loglik_fun <- function(Y, X_full, theta) {
  mu <- sum(X_full * theta$beta)
  dnorm(Y, mean = mu, sd = sqrt(theta$sigma2), log = TRUE)
}

X_con_data_augmentation <- function(
    Y, X_row, q_index,
    beta, sigma2,
    mu_q, tau2_q,
    rng = function(mu, sd) rnorm(1, mean = mu, sd = sd)) {
  stopifnot(is.numeric(Y), length(Y) == 1L)
  stopifnot(is.numeric(X_row))
  stopifnot(is.numeric(beta), length(beta) == length(X_row))
  stopifnot(length(q_index) == 1L, q_index >= 1, q_index <= length(X_row))
  stopifnot(is.numeric(sigma2), length(sigma2) == 1L, sigma2 > 0)
  stopifnot(is.numeric(mu_q), length(mu_q) == 1L)
  stopifnot(is.numeric(tau2_q), length(tau2_q) == 1L, tau2_q > 0)
  
  beta_q <- beta[q_index]
  
  # partial residual eps_tilde_{-q} = Y - (sum_{j != q} beta_j X_j)
  # (assumes X_row contains the observed values for j != q)
  idx_other <- setdiff(seq_along(X_row), q_index)
  eps_tilde <- Y - sum(beta[idx_other] * X_row[idx_other])
  
  tau2_tilde <- 1 / (1 / tau2_q + (beta_q^2) / sigma2)
  mu_tilde   <- tau2_tilde * (mu_q / tau2_q + beta_q * eps_tilde / sigma2)
  
  x_imp <- rng(mu_tilde, sqrt(tau2_tilde))
  
  list(x_imp = x_imp, mu_tilde = mu_tilde, tau2_tilde = tau2_tilde)
}

X_bin_data_augmentation <- function(
    Y, X_obs, pi,
    loglik_fun, theta,
    q_name = NULL,
    rng = function(p) rbinom(1, size = 1, prob = p)) {
  stopifnot(is.numeric(Y), length(Y) == 1L)
  stopifnot(is.numeric(pi), length(pi) == 1L, pi >= 0, pi <= 1)
  
  # Build X_full under Xq=1 and Xq=0
  if (!is.null(names(X_obs))) {
    if (is.null(q_name)) stop("If X_obs is named, supply `q_name` for the missing covariate.")
    X1 <- X_obs; X1[q_name] <- 1
    X0 <- X_obs; X0[q_name] <- 0
  } else {
    # If unnamed, user must ensure X_obs is in the right order and append Xq at end
    X1 <- c(X_obs, 1)
    X0 <- c(X_obs, 0)
  }
  
  # Z is safely removed from these calls, relying entirely on X_full
  ll1 <- loglik_fun(Y = Y, X_full = X1, theta = theta)
  ll0 <- loglik_fun(Y = Y, X_full = X0, theta = theta)
  
  # Stable computation of ptilde using log-sum-exp
  # ptilde = 1 / (1 + ((1-pi)/pi) * exp(ll0-ll1))
  if (pi == 0) {
    ptilde <- 0
  } else if (pi == 1) {
    ptilde <- 1
  } else {
    log_ratio <- log1p(-pi) - log(pi) + (ll0 - ll1)
    ptilde <- 1 / (1 + exp(log_ratio))
  }
  
  x_imp <- as.integer(rng(ptilde))
  list(x_imp = x_imp, p_tilde = ptilde)
}

X_data_augmentation <- function(
    matX, Ysatr, Zobs,
    Sy, Sx, unique_Sy,
    betaPars, sig2Pars,
    piPars, muPars, tau2Pars,
    miss_bin, miss_con,
    p_Z, p_X1, p_X2,
    build_sid_count,
    X_bin_data_augmentation,
    X_con_data_augmentation,
    loglik_fun = NULL,
    use_Y = FALSE) {
  stopifnot(is.matrix(matX))
  N <- nrow(matX)
  stopifnot(length(Ysatr) == N, length(Zobs) == N, length(Sy) == N, length(Sx) == N)
  stopifnot(is.numeric(p_Z), p_Z >= 1)

  has_missing_X <- any(miss_bin) || any(miss_con)
  if (!has_missing_X) {
    scaXobs <- matX[, -(1:(1 + p_Z)), drop = FALSE]
    return(list(matX = matX, scaXobs = scaXobs, sid_vec = NULL, Kyx = NULL))
  }

  # (Sy,Sx) -> sid mapping (must match the SAME count convention used in parameter updates)
  sid_obj <- build_sid_count(Sy, Sx, unique_Sy)
  sid_vec <- sid_obj$sid
  Kyx     <- sid_obj$Kyx

  # Hard safety checks (prevents silent mis-indexing)
  if (p_X1 > 0 && Kyx > ncol(piPars)) stop("Kyx exceeds ncol(piPars): count/order mismatch.")
  if (p_X2 > 0 && Kyx > ncol(muPars)) stop("Kyx exceeds ncol(muPars): count/order mismatch.")
  if (p_X2 > 0 && Kyx > ncol(tau2Pars)) stop("Kyx exceeds ncol(tau2Pars): count/order mismatch.")

  for (i in seq_len(N)) {
    has_miss_bin <- (p_X1 > 0) && any(miss_bin[i, ], na.rm = TRUE)
    has_miss_con <- (p_X2 > 0) && any(miss_con[i, ], na.rm = TRUE)
    if (!has_miss_bin && !has_miss_con) next

    # outcome-cluster params indexed by Sy
    theta_i <- list(beta = betaPars[, Sy[i]], sigma2 = sig2Pars[Sy[i]])
    sid <- sid_vec[i]

    # Named design row (required for binary updates via q_name)
    X_row <- matX[i, ]
    names(X_row) <- colnames(matX)

    # ---- binary covariates ----
    if (p_X1 > 0 && has_miss_bin) {
      miss_j <- which(miss_bin[i, ])
      for (j in miss_j) {
        q_name <- colnames(matX)[1 + p_Z + j]  # covariates start at col (1+p_Z+1)

        if (use_Y) {
          if (is.null(loglik_fun)) stop("use_Y=TRUE requires loglik_fun.")
          out <- X_bin_data_augmentation(
            Y = Ysatr[i],
            Z = Zobs[i],
            X_obs = X_row,
            pi = piPars[j, sid],
            loglik_fun = loglik_fun,
            theta = theta_i,
            q_name = q_name
          )
        } else {
          # prior-only update (no Y likelihood)
          out <- list(x_imp = as.integer(rbinom(1, 1, piPars[j, sid])))
        }

        matX[i, 1 + p_Z + j] <- out$x_imp
        X_row[q_name] <- out$x_imp
      }
    }

    # ---- continuous covariates ----
    if (p_X2 > 0 && has_miss_con) {
      miss_j <- which(miss_con[i, ])
      for (j in miss_j) {
        sca_col <- p_X1 + j
        q_index <- 1 + p_Z + sca_col

        if (use_Y) {
          out <- X_con_data_augmentation(
            Y = Ysatr[i],
            X_row = matX[i, ],
            q_index = q_index,
            beta = theta_i$beta,
            sigma2 = theta_i$sigma2,
            mu_q = muPars[j, sid],
            tau2_q = tau2Pars[j, sid]
          )
        } else {
          # prior-only update (no Y likelihood)
          out <- list(x_imp = rnorm(1, mean = muPars[j, sid], sd = sqrt(tau2Pars[j, sid])))
        }

        matX[i, q_index] <- out$x_imp
      }
    }
  }

  scaXobs <- matX[, -(1:(1 + p_Z)), drop = FALSE]
  list(matX = matX, scaXobs = scaXobs, sid_vec = sid_vec, Kyx = Kyx)
}

# Y_data_augmentation = function(
#     Tobs, Dobs, matX, S,
#     a_beta, Binv_beta, aBinv_beta = NULL,
#     a_sig2, b_sig2,
#     inv_cpp, rtrunct_cpp,
#     upper = Inf) {
#   # ---- checks ----
#   stopifnot(is.numeric(Tobs), length(Tobs) > 0)
#   N <- length(Tobs)
#   stopifnot(length(Dobs) == N, length(S) == N)
#   stopifnot(is.matrix(matX), nrow(matX) == N)
#   p <- ncol(matX)
#   stopifnot(length(a_beta) == p)
#   stopifnot(is.matrix(Binv_beta), nrow(Binv_beta) == p, ncol(Binv_beta) == p)
#   
#   if (is.null(aBinv_beta)) aBinv_beta <- as.vector(Binv_beta %*% a_beta)
#   
#   # Coerce Dobs to integer 0/1
#   Dobs01 <- as.integer(Dobs)
#   if (!all(Dobs01 %in% c(0L, 1L))) stop("Dobs must be 0/1 (or coercible).")
#   
#   # ---- initialize ----
#   Ysatr <- numeric(N)
#   Ysatr[Dobs01 == 1L] <- Tobs[Dobs01 == 1L]
#   
#   unique_S <- unique(S)
#   K <- length(unique_S)
#   
#   # Helper to recycle upper appropriately
#   get_upper <- function(m) {
#     if (length(upper) == 1L) {
#       rep(upper, m)
#     } else if (length(upper) == m) {
#       upper
#     } else {
#       stop("`upper` must be length 1 or length equal to number of censored units in cluster.")
#     }
#   }
#   
#   for (k in seq_len(K)) {
#     ind_curr <- which(S == unique_S[k])
#     if (length(ind_curr) == 0L) next
#     
#     ind_cens <- ind_curr[Dobs01[ind_curr] == 0L]
#     if (length(ind_cens) == 0L) next
#     
#     # Use current augmented times for all in cluster (events fixed; censored updated progressively)
#     y_full <- Ysatr[ind_curr]
#     X_full <- matX[ind_curr, , drop = FALSE]
#     
#     # ---- collapsed posterior (beta, sigma^2) given augmented y_full ----
#     Binv_n  <- Binv_beta + crossprod(X_full)
#     Bn_beta <- inv_cpp(Binv_n)
#     an_beta <- Bn_beta %*% (aBinv_beta + crossprod(X_full, y_full))
#     
#     residuals <- as.vector(y_full - X_full %*% an_beta)
#     beta_diff <- as.vector(an_beta - a_beta)
#     
#     sse_data  <- sum(residuals^2)
#     sse_prior <- as.numeric(t(beta_diff) %*% Binv_beta %*% beta_diff)
#     
#     an_sig2 <- a_sig2 + length(y_full)
#     bn_sig2 <- a_sig2 * b_sig2 + sse_data + sse_prior
#     
#     # ---- predictive truncated t for censored units ----
#     X_cens    <- matX[ind_cens, , drop = FALSE]
#     lower_vec <- Tobs[ind_cens]
#     
#     quad_form <- rowSums((X_cens %*% Bn_beta) * X_cens)
#     scale_val <- sqrt((bn_sig2 / an_sig2) * (1 + quad_form))
#     mu_val    <- as.vector(X_cens %*% an_beta)
#     
#     upper_vec <- get_upper(length(ind_cens))
#     
#     # sample one draw per censored obs
#     Ysatr[ind_cens] <- mapply(
#       FUN = function(mu, scale, lower, upper_i) {
#         rtrunct_cpp(n = 1, df = an_sig2, mu = mu, scale = scale, lower = lower, upper = upper_i)
#       },
#       mu = mu_val,
#       scale = scale_val,
#       lower = lower_vec,
#       upper_i = upper_vec
#     )
#   }
#   
#   return(Ysatr)
# }

Y_data_augmentation = function(Tobs, Dobs, matX, S, betaPars, sig2Pars) {
  N <- length(Tobs)
  Dobs01 <- as.integer(Dobs)
  
  Ysatr <- numeric(N)
  Ysatr[Dobs01 == 1L] <- Tobs[Dobs01 == 1L] # Keep actual events fixed
  
  ind_cens <- which(Dobs01 == 0L)
  if (length(ind_cens) == 0L) return(Ysatr)
  
  for (i in ind_cens) {
    k <- S[i]
    # Calculate the current cluster-specific mean and standard deviation
    mu_i <- sum(matX[i, ] * betaPars[, k])
    sig_i <- sqrt(sig2Pars[1, k]) 
    
    # Draw from Truncated Normal via Inverse CDF Method
    p_min <- pnorm(Tobs[i], mean = mu_i, sd = sig_i)
    u <- runif(1, min = p_min, max = 1)
    
    # Protect against u exactly 1 yielding Inf
    if (u >= 1) u <- 1 - 1e-15 
    
    y_imp <- qnorm(u, mean = mu_i, sd = sig_i)
    
    # Fallback for extreme tail probabilities
    if (is.infinite(y_imp)) y_imp <- Tobs[i] + 0.01 
    
    Ysatr[i] <- y_imp
  }
  
  return(Ysatr)
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
  
  # --- Identify covariate types (binary vs. continuous) ---
  is_binary <- apply(Xobs, 2, function(col) {
    col2 <- col[!is.na(col)]
    if (is.logical(col2)) col2 <- as.integer(col2)
    length(col2) > 0 && all(col2 %in% c(0, 1))
  })
  
  X_bin <- Xobs[, is_binary, drop = FALSE]
  X_con <- Xobs[, !is_binary, drop = FALSE]
  p_X1 <- ncol(X_bin)
  p_X2 <- ncol(X_con)
  p_X = p_X1 + p_X2
  p_Z = length(Zobs)/N
  Zobs = matrix(Zobs, ncol = p_Z)
  
  # --- Missingness indicators (fixed throughout) ---
  # Define helper function globally for safety
  build_sid_count <- function(Sy, Sx, unique_Sy) {
    N <- length(Sy)
    sid <- integer(N)
    count <- 1L
    
    for (k in unique_Sy) {
      ind_Sy <- (Sy == k)
      unique_Sx_temp <- sort(unique(Sx[ind_Sy]))
      for (r in unique_Sx_temp) {
        ind_kr <- ind_Sy & (Sx == r)
        sid[ind_kr] <- count
        count <- count + 1L
      }
    }
    
    if (any(sid == 0L)) stop("sid has zeros: some (Sy,Sx) pairs were not mapped.")
    list(sid = sid, Kyx = count - 1L)
  }
  miss_bin <- if (p_X1 > 0) is.na(X_bin) else matrix(FALSE, nrow(Xobs), 0)
  miss_con <- if (p_X2 > 0) is.na(X_con) else matrix(FALSE, nrow(Xobs), 0)
  has_missing_X <- (any(miss_bin) || any(miss_con))
  
  # --- Initial deterministic imputation (starting values only) ---
  if (p_X1 > 0) {
    for (j in seq_len(p_X1)) {
      idx <- which(miss_bin[, j])
      if (!length(idx)) next
      mj <- mean(X_bin[, j], na.rm = TRUE)
      if (is.nan(mj)) mj <- 0 
      X_bin[idx, j] <- as.integer(mj >= 0.5)
    }
  }
  
  if (p_X2 > 0) {
    for (j in seq_len(p_X2)) {
      idx <- which(miss_con[, j])
      if (!length(idx)) next
      mj <- mean(X_con[, j], na.rm = TRUE)
      if (is.nan(mj)) mj <- 0
      X_con[idx, j] <- mj
    }
  }
  
  # --- Precompute scaling ONCE (fixed for entire run) ---
  con_center <- con_scale <- NULL
  if (p_X2 > 0) {
    con_center <- colMeans(X_con, na.rm = TRUE)
    con_scale  <- apply(X_con, 2, sd, na.rm = TRUE)
    con_scale[!is.finite(con_scale) | con_scale <= 0] <- 1
    scaX_con <- sweep(X_con, 2, con_center, "-")
    scaX_con <- sweep(scaX_con, 2, con_scale, "/")
  }
  
  # --- Build scaled covariate matrix used in the model ---
  if (p_X1 > 0 && p_X2 > 0) {
    scaXobs <- cbind(X_bin, scaX_con)
  } else if (p_X2 > 0) {
    scaXobs <- scaX_con
  } else {
    scaXobs <- X_bin
  }
  
  matX <- cbind(1, Zobs, scaXobs)
  colnames(matX) <- c("int", "Z", paste0("X", seq_len(ncol(scaXobs))))
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
  B_beta = vcov(survreg_beta)[1:p_beta, 1:p_beta, drop = FALSE]
  c_beta = N
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
  # c_beta = N
  # B_beta = c_beta * B_beta
  # Binv_beta = inv_cpp(B_beta)
  # aBinv_beta = Binv_beta %*% a_beta
  # a_sig2 = 1
  # b_sig2 = var(residuals_aft[Dobs==1]) * (a_sig2 + 2)
  a_sig2 = 3
  b_sig2 = 0.1
  
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
  # Ysatr = Y_data_augmentation(Tobs, Dobs, matX, S,
  #                             a_beta, Binv_beta, aBinv_beta,
  #                             a_sig2, b_sig2,
  #                             inv_cpp, rtrunct_cpp,
  #                             upper = Inf)
  Ysatr = Y_data_augmentation(Tobs, Dobs, matX, S, betaPars, sig2Pars)
  
  # --- Initialize propensity parameters (eta) ---
  etaPars = matrix(nrow = p_Z, ncol = K)
  # --- Initialize covariate parameters (omega) ---
  piPars = matrix(nrow = p_X1, ncol = K) # Probabilities
  muPars = matrix(nrow = p_X2, ncol = K) # Means
  tau2Pars = matrix(nrow = p_X2, ncol = K) # Variances
  
  for (k in 1:K) {
    ind_S = which(S == unique_S[k])
    n_temp = length(ind_S)
    T_temp = Ysatr[ind_S]
    matX_temp = matrix(matX[ind_S,], ncol = p_matX)
    
    # update parameters for true failure time
    # regression_results = update_reg_con(T_temp, matX_temp, betaPars[,k],
    #                                     a_beta, Binv_beta, aBinv_beta, a_sig2, b_sig2)
    regression_results = update_reg_con(T_temp, matX_temp, a_beta, Binv_beta, a_sig2, b_sig2)
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
    f0_x_all[, q] = beta(a_eta + Zobs[,q], b_eta - Zobs[,q] + 1) # beta(a_pi,b_pi) = 1 in this case
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
  n_MCMC = floor(gibbs_iter / gibbs_thin)
  
  # Define iteration check
  iter_check = floor(gibbs_iter/20)
  
  # Make vectors to store draws from Gibbs Sampler
  n_MCMC = floor(gibbs_iter/gibbs_thin)
  
  alpha_draws = numeric(n_MCMC)
  betaLists = vector("list", n_MCMC)
  sig2Lists = vector("list", n_MCMC)
  etaLists  = vector("list", n_MCMC)
  piLists   = vector("list", n_MCMC)
  muLists   = vector("list", n_MCMC)
  tau2Lists = vector("list", n_MCMC)
  n_kLists  = vector("list", n_MCMC)
  KLists    = vector("list", n_MCMC)
  
  count_it = 1
  # End initial values ------------------------------------------------------------
  
  # First draw each parameter for BNP model. Then calculate causal effect.
  for (gibbs_reps in 1:gibbs_total) {
    # gibbs_reps = 1
    # Update of all parameters in BNP model--------------------------------
    # Update cluster membership -------------------------------------------
    # sourceCpp("BNPsurvival_cpp.cpp")
    scaXobs <- matX[, -(1:(1 + p_Z)), drop = FALSE]
    cluster_res = DPMMcluster_cpp(N, p_X1, p_X2, 
                                  Ysatr, matX, Zobs, scaXobs,
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
    # Ysatr = Y_data_augmentation(Tobs, Dobs, matX, S,
    #                             a_beta, Binv_beta, aBinv_beta,
    #                             a_sig2, b_sig2,
    #                             inv_cpp, rtrunct_cpp,
    #                             upper = Inf)
    Ysatr = Y_data_augmentation(Tobs, Dobs, matX, S, betaPars, sig2Pars)
    
    # Missing covariates augmentation (DPMM specific single-cluster logic)
    if (has_missing_X) {
      for (i in seq_len(N)) {
        has_miss_bin <- (p_X1 > 0) && any(miss_bin[i, ], na.rm = TRUE)
        has_miss_con <- (p_X2 > 0) && any(miss_con[i, ], na.rm = TRUE)
        if (!has_miss_bin && !has_miss_con) next
        
        k <- S[i]
        theta_i <- list(beta = betaPars[, k], sigma2 = sig2Pars[1, k])
        X_row <- matX[i, ]
        names(X_row) <- colnames(matX)
        
        if (p_X1 > 0 && has_miss_bin) {
          for (j in which(miss_bin[i, ])) {
            target_col <- 1 + p_Z + j
            q_name <- colnames(matX)[target_col]
            out <- X_bin_data_augmentation(
              Y = Ysatr[i], X_obs = X_row, pi = piPars[j, k],
              loglik_fun = loglik_fun, theta = theta_i, q_name = q_name
            )
            matX[i, target_col] <- out$x_imp
            X_row[q_name] <- out$x_imp
          }
        }
        
        if (p_X2 > 0 && has_miss_con) {
          for (j in which(miss_con[i, ])) {
            q_index <- 1 + p_Z + p_X1 + j
            out <- X_con_data_augmentation(
              Y = Ysatr[i], X_row = matX[i, ], q_index = q_index,
              beta = theta_i$beta, sigma2 = theta_i$sigma2,
              mu_q = muPars[j, k], tau2_q = tau2Pars[j, k]
            )
            matX[i, q_index] <- out$x_imp
          }
        }
      }
    }
    scaXobs <- matX[, -(1:(1 + p_Z)), drop = FALSE]
    
    for (k in 1:K) {
      ind_S = which(S == unique_S[k])
      n_temp = length(ind_S)
      T_temp = Ysatr[ind_S]
      matX_temp = matrix(matX[ind_S,], ncol = p_matX)
      
      # update parameters for true failure time
      # regression_results = update_reg_con(T_temp, matX_temp, betaPars[,k],
      #                                     a_beta, Binv_beta, aBinv_beta, a_sig2, b_sig2)
      regression_results = update_reg_con(T_temp, matX_temp, a_beta, Binv_beta, a_sig2, b_sig2)
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
      if ((gibbs_reps - gibbs_burnin) %% gibbs_thin == 0) {
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
  constants = list(N = N, n_MCMC = n_MCMC,
                   p_matX = p_matX, p_X1 = p_X1, p_X2 = p_X2,
                   con_center = con_center, con_scale = con_scale)
  
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
  
  return(MCMCresult)
}

# -----------------------------------------------------------------------------
# Define function for the post-processing under the BNP (DPMM) model
SURVIVAL_DPMM_POST = function(object, M, nu, rho, eta_z0 = NULL, eta_z1 = NULL, 
                              condX = NULL, p_condX = NULL, esttype = "mean") {
  
  # -----------------------------------------------------------------------------
  # object = DPMM_results_MCMC
  # condX = vector of values to condition on (e.g., c(1)), NULL for marginal
  # p_condX = vector of column indices in the X matrix for condX, NULL for marginal
  # eta_z0, eta_z1 = Vectors of sensitivity shift parameters
  # -----------------------------------------------------------------------------
  
  normalize_log_probs = function(log_prob_matrix) {
    max_log_probs = apply(log_prob_matrix, 1, max)
    log_probs_shifted = log_prob_matrix - max_log_probs
    probs_shifted = exp(log_probs_shifted)
    row_sums = rowSums(probs_shifted)
    
    zero_rows = which(row_sums == 0)
    if(length(zero_rows) > 0) {
      probs_shifted[zero_rows,] = 1 / ncol(log_prob_matrix) 
      row_sums[zero_rows] = 1
    }
    return(probs_shifted / row_sums)
  }
  
  # -----------------------------------------------------------------------------
  # update_parameter (Generation is handled upstream for conditional support)
  # -----------------------------------------------------------------------------
  update_parameter = function(S_x_mc, 
                              betaPars_pot, sig2Pars_pot, etaPars_pot,
                              piPars_pot, muPars_pot, tau2Pars_pot) {
    
    unique_S_x_mc = sort(unique(S_x_mc))
    K_x_mc = length(unique_S_x_mc)
    
    n_k_mc = numeric(K_x_mc)
    for (k in 1:K_x_mc) {
      n_k_mc[k] = sum(S_x_mc == unique_S_x_mc[k])
    }
    
    # Subset matrices to keep only active clusters
    betaPars_mc = matrix(betaPars_pot[, unique_S_x_mc], ncol = K_x_mc)
    sig2Pars_mc = matrix(sig2Pars_pot[, unique_S_x_mc], ncol = K_x_mc)
    etaPars_mc  = matrix(etaPars_pot[, unique_S_x_mc],  ncol = K_x_mc)
    piPars_mc   = matrix(piPars_pot[, unique_S_x_mc],   ncol = K_x_mc)
    muPars_mc   = matrix(muPars_pot[, unique_S_x_mc],   ncol = K_x_mc)
    tau2Pars_mc = matrix(tau2Pars_pot[, unique_S_x_mc], ncol = K_x_mc)
    
    # Remap indices
    S_x_mc_remapped = as.numeric(factor(S_x_mc, levels = unique_S_x_mc))
    
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
  nu_vec = as.vector(nu)
  rho_vec = as.vector(rho)
  n_nu = length(nu_vec)
  n_rho = length(rho_vec)
  
  # Build a full grid of all combinations of eta_z0 and eta_z1
  eta_grid = expand.grid(
    eta_z0 = if(is.null(eta_z0)) 0 else as.vector(eta_z0),
    eta_z1 = if(is.null(eta_z1)) 0 else as.vector(eta_z1)
  )
  eta_z0_vec = eta_grid$eta_z0
  eta_z1_vec = eta_grid$eta_z1
  n_eta = nrow(eta_grid)
  
  level = 0.05
  quantile_alpha = c(level/2, 1-level/2)
  z0 = 0
  z1 = 1
  
  N = object$constants$N
  # num_MC = object$constants$num_MC
  n_MCMC = object$constants$n_MCMC
  p_matX = object$constants$p_matX
  p_X1 = object$constants$p_X1
  p_X2 = object$constants$p_X2
  
  # Priors for generating the K+1 cluster
  a_beta = object$priors$a_beta; B_beta = object$priors$B_beta
  a_sig2 = object$priors$a_sig2; b_sig2 = object$priors$b_sig2
  a_eta  = object$priors$a_eta;  b_eta  = object$priors$b_eta
  a_pi   = object$priors$a_pi;   b_pi   = object$priors$b_pi
  a_mu   = object$priors$a_mu;   b_mu   = object$priors$b_mu
  a_tau2 = object$priors$a_tau2; b_tau2 = object$priors$b_tau2
  
  iter_check = floor(n_MCMC/10)
  
  # 4D Arrays for Sensitivity
  Y_rho_nu_z0x_store = array(NA, dim = c(n_MCMC, n_eta, n_nu, n_rho))
  Y_rho_nu_z1x_store = array(NA, dim = c(n_MCMC, n_eta, n_nu, n_rho))
  
  for (post_reps in 1:n_MCMC) {
    
    # 1. Extract Current MCMC State
    alpha = object$MCMCposteriors$alpha_draws[post_reps]
    K = object$MCMCposteriors$KLists[[post_reps]]
    n_k = object$MCMCposteriors$n_kLists[[post_reps]]
    betaPars = object$MCMCposteriors$betaLists[[post_reps]]
    sig2Pars = object$MCMCposteriors$sig2Lists[[post_reps]]
    etaPars = object$MCMCposteriors$etaLists[[post_reps]]
    piPars = object$MCMCposteriors$piLists[[post_reps]]
    muPars = object$MCMCposteriors$muLists[[post_reps]]
    tau2Pars = object$MCMCposteriors$tau2Lists[[post_reps]]
    
    # ---------------------------------------------------------------------------
    # 2. Pre-generate potential K+1 cluster parameters from prior
    # ---------------------------------------------------------------------------
    sig2_new = c(rscainvchisq_cpp(1, a_sig2, b_sig2))
    beta_new = c(rmvn_cpp(1, a_beta, B_beta)) 
    eta_new  = rbeta(1, a_eta, b_eta)
    
    pi_new   = numeric(p_X1)
    mu_new   = numeric(p_X2)
    tau2_new = numeric(p_X2)
    
    if (p_X1 > 0) pi_new = rbeta(p_X1, a_pi, b_pi)
    if (p_X2 > 0) {
      for (q in 1:p_X2) {
        tau2_new[q] = rinvgamma_cpp(1, a_tau2, b_tau2)
        mu_new[q]   = rnorm(1, mean = a_mu, sd = sqrt(tau2_new[q] / b_mu))
      }
    }
    
    # Form "potential" matrices (Current K + 1 New)
    betaPars_pot = cbind(betaPars, beta_new)
    sig2Pars_pot = cbind(sig2Pars, sig2_new)
    etaPars_pot  = cbind(etaPars, eta_new)
    piPars_pot   = cbind(piPars, pi_new)
    muPars_pot   = cbind(muPars, mu_new)
    tau2Pars_pot = cbind(tau2Pars, tau2_new)
    
    # ---------------------------------------------------------------------------
    # 3. Calculate Conditional lambda_x_mc
    # ---------------------------------------------------------------------------
    if (!is.null(condX) && !is.null(p_condX)) {
      log_lik_cond_X = numeric(K + 1)
      
      # Calculate log-likelihood of condX for all K+1 clusters
      for (k in 1:(K + 1)) {
        ll = 0
        for (i in seq_along(p_condX)) {
          idx = p_condX[i]
          val = condX[i]
          if (idx <= p_X1) {
            ll = ll + dbinom(val, 1, piPars_pot[idx, k], log = TRUE)
          } else {
            q = idx - p_X1
            ll = ll + dnorm(val, muPars_pot[q, k], sqrt(tau2Pars_pot[q, k]), log = TRUE)
          }
        }
        log_lik_cond_X[k] = ll
      }
      
      # Bayes Rule: log(Prior Weight) + log(Likelihood)
      log_lambda_x_mc = log(c(n_k, alpha)) + log_lik_cond_X
      
      # Exponentiate safely to prevent numerical underflow
      max_ll = max(log_lambda_x_mc)
      lambda_x_mc = exp(log_lambda_x_mc - max_ll)
    } else {
      # Marginal mode: purely DP prior weights
      lambda_x_mc = c(n_k, alpha)
    }
    
    # ---------------------------------------------------------------------------
    # 4. Draw clusters and update parameters
    # ---------------------------------------------------------------------------
    S_x_mc = sapply(1:M, function(l) rmultinom_cpp(lambda_x_mc))
    orig_unique_S_x_mc = sort(unique(S_x_mc))
    
    # DP predictive weights for the chosen clusters
    w_dp = (c(n_k, alpha) / (N + alpha))[orig_unique_S_x_mc]
    
    # Subsets the _pot matrices down to only the clusters that were drawn
    updated_pars = update_parameter(S_x_mc, 
                                    betaPars_pot, sig2Pars_pot, etaPars_pot, 
                                    piPars_pot, muPars_pot, tau2Pars_pot)
    
    S_x_mc = updated_pars$S_x_mc
    unique_S_x_mc = updated_pars$unique_S_x_mc
    K = updated_pars$K_x_mc
    n_k = updated_pars$n_k_mc
    
    betaPars = updated_pars$betaPars_mc
    sig2Pars = updated_pars$sig2Pars_mc
    etaPars  = updated_pars$etaPars_mc
    piPars   = updated_pars$piPars_mc
    muPars   = updated_pars$muPars_mc
    tau2Pars = updated_pars$tau2Pars_mc
    
    # -----------------------------------------------------------------------------
    # 5. Draw confounders for each of the M observations
    # -----------------------------------------------------------------------------
    X_mc = matrix(nrow = M, ncol = (p_X1 + p_X2))
    
    # Pre-assign the conditioned values
    if (!is.null(condX) && !is.null(p_condX)) {
      for (i in seq_along(p_condX)) {
        X_mc[, p_condX[i]] = condX[i]
      }
    }
    
    for (k in 1:K) {
      ind_S_x_mc = which(S_x_mc == unique_S_x_mc[k])
      n_k_temp = length(ind_S_x_mc)
      
      if (n_k_temp > 0) {
        if (p_X1 > 0) {
          for (q in 1:p_X1) {
            # Skip sampling if variable is conditioned
            if (is.null(p_condX) || !(q %in% p_condX)) {
              X_mc[ind_S_x_mc, q] = rbinom(n_k_temp, 1, piPars[q,k])
            }
          }
        }
        if (p_X2 > 0) {
          for (q in 1:p_X2) {
            col_idx = p_X1 + q
            # Skip sampling if variable is conditioned
            if (is.null(p_condX) || !(col_idx %in% p_condX)) {
              X_mc[ind_S_x_mc, col_idx] = rnorm(n_k_temp, muPars[q,k], sqrt(tau2Pars[q,k]))
            }
          }
        }
      }
    }
    
    # -----------------------------------------------------------------------------
    # 6. Calculate Y probabilities 
    # -----------------------------------------------------------------------------
    lambda_Y_z0x_mc = matrix(nrow = M, ncol = K)
    lambda_Y_z1x_mc = matrix(nrow = M, ncol = K)
    
    for (k in 1:K) {
      log_w_k_temp = log(w_dp[k])
      
      prob_z0_mc = dbinom(z0, 1, etaPars[,k], log = TRUE)
      prob_z1_mc = dbinom(z1, 1, etaPars[,k], log = TRUE)
      
      prob_x_mc = numeric(M)
      if (p_X1 > 0) {
        for (q in 1:p_X1) {
          prob_x_mc = prob_x_mc + dbinom(X_mc[, q], 1, piPars[q,k], log = TRUE)
        }
      }
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
    # 7. Create the design matrices for all M subjects
    # -----------------------------------------------------------------------------
    # Create the design matrices for all M subjects
    matX_z0x_mc = cbind(1, z0, X_mc) 
    matX_z1x_mc = cbind(1, z1, X_mc) 
    MUy_z0x_mc = (matX_z0x_mc %*% betaPars)
    MUy_z1x_mc = (matX_z1x_mc %*% betaPars)
    Ey_z0x_mc = rowSums(MUy_z0x_mc * lambda_Y_z0x_mc)
    Ey_z1x_mc = rowSums(MUy_z1x_mc * lambda_Y_z1x_mc)
    Ey2_z0x_mc = rowSums((MUy_z0x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = K)) * lambda_Y_z0x_mc)
    Ey2_z1x_mc = rowSums((MUy_z1x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = K)) * lambda_Y_z1x_mc)
    # SDy_z0x_mc = sqrt(pmax(Ey2_z0x_mc - Ey_z0x_mc^{2}, 0))
    # SDy_z1x_mc = sqrt(pmax(Ey2_z1x_mc - Ey_z1x_mc^{2}, 0))
    
    # -----------------------------------------------------------------------------
    lambda_Y_z0x = apply(lambda_Y_z0x_mc, 2, mean)
    lambda_Y_z1x = apply(lambda_Y_z1x_mc, 2, mean)
    MUy_z0x = apply(MUy_z0x_mc, 2, mean)
    MUy_z1x = apply(MUy_z1x_mc, 2, mean)
    Ey_z0x = mean(Ey_z0x_mc)
    Ey_z1x = mean(Ey_z1x_mc)
    Ey2_z0x = mean(Ey2_z0x_mc)
    Ey2_z1x = mean(Ey2_z1x_mc)
    SDy_z0x = sqrt(pmax(Ey2_z0x - Ey_z0x^{2}, 0))
    SDy_z1x = sqrt(pmax(Ey2_z1x - Ey_z1x^{2}, 0))
    
    # -----------------------------------------------------------------------------
    # SENSITIVITY AND QRL OPTIMIZATION LOOPS
    # -----------------------------------------------------------------------------
    for (ii_eta in 1:n_eta) {
      curr_eta_z0 = eta_z0_vec[ii_eta]
      curr_eta_z1 = eta_z1_vec[ii_eta]
      for (ii_nu in 1:n_nu) {
        curr_nu = nu_vec[ii_nu]
        curr_log_nu = log(curr_nu)
        # Condition to run optimization (skip and leave as NA if < 0)
        run_z0 = (curr_nu - curr_eta_z0 >= 0)
        run_z1 = (curr_nu - curr_eta_z1 >= 0)
        # Apply the sensitivity shift (psi_z) directly to the log evaluation point
        if (run_z0) {
          S_nu_z0x = S_cond_cpp(curr_log_nu - curr_eta_z0, lambda_Y_z0x_mc, MUy_z0x_mc, sqrt(sig2Pars), TRUE)
        }
        if (run_z1) {
          S_nu_z1x = S_cond_cpp(curr_log_nu - curr_eta_z1, lambda_Y_z1x_mc, MUy_z1x_mc, sqrt(sig2Pars), TRUE)
        }
        safe_curr_log_nu = ifelse(curr_nu == 0, -100, curr_log_nu)
        for(ii_rho in 1:n_rho) {
          curr_rho = rho_vec[ii_rho]
          if (run_z0) {
            target_rho_nu_z0x = S_nu_z0x + log(1 - curr_rho)
            opt_z0 <- S_optim_cond_cpp(target_rho_nu_z0x, Ey_z0x,
                                       (safe_curr_log_nu - curr_eta_z0) + 1e-10, 100,
                                       lambda_Y_z0x_mc, MUy_z0x_mc, sqrt(sig2Pars), TRUE)
            Y_rho_nu_z0x_store[post_reps, ii_eta, ii_nu, ii_rho] = exp(opt_z0$optimizer + curr_eta_z0) - curr_nu
          }
          if (run_z1) {
            target_rho_nu_z1x = S_nu_z1x + log(1 - curr_rho)
            opt_z1 <- S_optim_cond_cpp(target_rho_nu_z1x, Ey_z1x,
                                       (safe_curr_log_nu - curr_eta_z1) + 1e-10, 100,
                                       lambda_Y_z1x_mc, MUy_z1x_mc, sqrt(sig2Pars), TRUE)
            Y_rho_nu_z1x_store[post_reps, ii_eta, ii_nu, ii_rho] = exp(opt_z1$optimizer + curr_eta_z1) - curr_nu
          }
        } 
      } 
    }
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  E_rho_nu_diff_store = Y_rho_nu_z1x_store - Y_rho_nu_z0x_store
  
  z0_tables_list   = list()
  z1_tables_list   = list()
  diff_tables_list = list()
  counter = 1
  for(ii_eta in 1:n_eta) {
    val_eta_z0 = eta_z0_vec[ii_eta]
    val_eta_z1 = eta_z1_vec[ii_eta]
    
    for(ii_nu in 1:n_nu) {
      val_nu  = nu_vec[ii_nu]
      
      for(ii_rho in 1:n_rho) {
        val_rho = rho_vec[ii_rho]
        
        # 1. Extract draws from the new 4D arrays
        draws_z0   = Y_rho_nu_z0x_store[, ii_eta, ii_nu, ii_rho]
        draws_z1   = Y_rho_nu_z1x_store[, ii_eta, ii_nu, ii_rho]
        draws_diff = E_rho_nu_diff_store[, ii_eta, ii_nu, ii_rho]
        
        # 2. Safely summarize only if the optimization was run (not all NAs)
        if (!all(is.na(draws_z0))) {
          summ_z0 = POSTsummary(draws_z0, esttype, quantile_alpha)
          # Add eta_z0 to the table output so you can track the sensitivity shift
          summ_z0 = cbind(eta_z0 = val_eta_z0, nu = val_nu, rho = val_rho, summ_z0)
          z0_tables_list[[counter]] = summ_z0
        }
        
        if (!all(is.na(draws_z1))) {
          summ_z1 = POSTsummary(draws_z1, esttype, quantile_alpha)
          # Add eta_z1 to the table output
          summ_z1 = cbind(eta_z1 = val_eta_z1, nu = val_nu, rho = val_rho, summ_z1)
          z1_tables_list[[counter]] = summ_z1
        }
        
        if (!all(is.na(draws_diff))) {
          summ_diff = POSTsummary(draws_diff, esttype, quantile_alpha)
          # For the difference table, it's helpful to see both shift parameters
          summ_diff = cbind(eta_z0 = val_eta_z0, eta_z1 = val_eta_z1, nu = val_nu, rho = val_rho, summ_diff)
          diff_tables_list[[counter]] = summ_diff
        }
        
        counter = counter + 1
      }
    }
  }
  
  # 3. do.call(rbind) will automatically ignore any NULL elements in the list 
  # where the NA condition triggered a skip.
  Y_rho_nu_z0x_result  = if(length(z0_tables_list) > 0) do.call(rbind, z0_tables_list) else NULL
  Y_rho_nu_z1x_result  = if(length(z1_tables_list) > 0) do.call(rbind, z1_tables_list) else NULL
  E_rho_nu_diff_result = if(length(diff_tables_list) > 0) do.call(rbind, diff_tables_list) else NULL
  
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

SURVIVAL_DPMM_POST_parallel = function(object, M, nu, rho,
                                       eta_z0 = NULL, eta_z1 = NULL,
                                       condX = NULL, p_condX = NULL,
                                       esttype = "mean",
                                       ncores = max(1L, parallel::detectCores() - 1L),
                                       seed = NULL) {
  # -----------------------------------------------------------------------------
  # Parallel version of SURVIVAL_DPMM_POST
  # Uses mclapply() on macOS/Linux; falls back to serial on Windows or ncores <= 1
  # -----------------------------------------------------------------------------
  
  if (.Platform$OS.type == "windows" || ncores <= 1L) {
    if (.Platform$OS.type == "windows") {
      message("Windows detected: parallel::mclapply() is unavailable; falling back to SURVIVAL_DPMM_POST().")
    } else {
      message("ncores <= 1; falling back to SURVIVAL_DPMM_POST().")
    }
    return(
      SURVIVAL_DPMM_POST(
        object   = object,
        M        = M,
        nu       = nu,
        rho      = rho,
        eta_z0   = eta_z0,
        eta_z1   = eta_z1,
        condX    = condX,
        p_condX  = p_condX,
        esttype  = esttype
      )
    )
  }
  
  normalize_log_probs = function(log_prob_matrix) {
    max_log_probs = apply(log_prob_matrix, 1, max)
    log_probs_shifted = log_prob_matrix - max_log_probs
    probs_shifted = exp(log_probs_shifted)
    row_sums = rowSums(probs_shifted)
    
    zero_rows = which(row_sums == 0)
    if (length(zero_rows) > 0) {
      probs_shifted[zero_rows, ] = 1 / ncol(log_prob_matrix)
      row_sums[zero_rows] = 1
    }
    probs_shifted / row_sums
  }
  
  update_parameter = function(S_x_mc,
                              betaPars_pot, sig2Pars_pot, etaPars_pot,
                              piPars_pot, muPars_pot, tau2Pars_pot) {
    
    unique_S_x_mc = sort(unique(S_x_mc))
    K_x_mc = length(unique_S_x_mc)
    
    n_k_mc = numeric(K_x_mc)
    for (k in seq_len(K_x_mc)) {
      n_k_mc[k] = sum(S_x_mc == unique_S_x_mc[k])
    }
    
    betaPars_mc = matrix(betaPars_pot[, unique_S_x_mc], ncol = K_x_mc)
    sig2Pars_mc = matrix(sig2Pars_pot[, unique_S_x_mc], ncol = K_x_mc)
    etaPars_mc  = matrix(etaPars_pot[, unique_S_x_mc],  ncol = K_x_mc)
    piPars_mc   = matrix(piPars_pot[, unique_S_x_mc],   ncol = K_x_mc)
    muPars_mc   = matrix(muPars_pot[, unique_S_x_mc],   ncol = K_x_mc)
    tau2Pars_mc = matrix(tau2Pars_pot[, unique_S_x_mc], ncol = K_x_mc)
    
    S_x_mc_remapped = as.numeric(factor(S_x_mc, levels = unique_S_x_mc))
    
    list(
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
    )
  }
  
  # -----------------------------------------------------------------------------
  # Extract constants and settings
  # -----------------------------------------------------------------------------
  nu_vec = as.vector(nu)
  rho_vec = as.vector(rho)
  n_nu = length(nu_vec)
  n_rho = length(rho_vec)
  
  eta_grid = expand.grid(
    eta_z0 = if (is.null(eta_z0)) 0 else as.vector(eta_z0),
    eta_z1 = if (is.null(eta_z1)) 0 else as.vector(eta_z1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  eta_z0_vec = eta_grid$eta_z0
  eta_z1_vec = eta_grid$eta_z1
  n_eta = nrow(eta_grid)
  
  level = 0.05
  quantile_alpha = c(level / 2, 1 - level / 2)
  z0 = 0
  z1 = 1
  
  N = object$constants$N
  # num_MC = object$constants$num_MC
  n_MCMC = object$constants$n_MCMC
  p_matX = object$constants$p_matX
  p_X1 = object$constants$p_X1
  p_X2 = object$constants$p_X2
  
  a_beta = object$priors$a_beta
  B_beta = object$priors$B_beta
  a_sig2 = object$priors$a_sig2
  b_sig2 = object$priors$b_sig2
  a_eta  = object$priors$a_eta
  b_eta  = object$priors$b_eta
  a_pi   = object$priors$a_pi
  b_pi   = object$priors$b_pi
  a_mu   = object$priors$a_mu
  b_mu   = object$priors$b_mu
  a_tau2 = object$priors$a_tau2
  b_tau2 = object$priors$b_tau2
  
  # -----------------------------------------------------------------------------
  # Worker
  # -----------------------------------------------------------------------------
  worker_fun = function(post_reps) {
    
    # 1. Extract current MCMC state
    alpha = object$MCMCposteriors$alpha_draws[post_reps]
    K = object$MCMCposteriors$KLists[[post_reps]]
    n_k = object$MCMCposteriors$n_kLists[[post_reps]]
    betaPars = object$MCMCposteriors$betaLists[[post_reps]]
    sig2Pars = object$MCMCposteriors$sig2Lists[[post_reps]]
    etaPars = object$MCMCposteriors$etaLists[[post_reps]]
    piPars = object$MCMCposteriors$piLists[[post_reps]]
    muPars = object$MCMCposteriors$muLists[[post_reps]]
    tau2Pars = object$MCMCposteriors$tau2Lists[[post_reps]]
    
    sqrt_sig2 = as.vector(sqrt(sig2Pars))
    
    Y_rho_nu_z0x_local = array(NA_real_, dim = c(n_eta, n_nu, n_rho))
    Y_rho_nu_z1x_local = array(NA_real_, dim = c(n_eta, n_nu, n_rho))
    
    # ---------------------------------------------------------------------------
    # 2. Pre-generate potential K+1 cluster parameters from prior
    # ---------------------------------------------------------------------------
    sig2_new = c(rscainvchisq_cpp(1, a_sig2, b_sig2))
    beta_new = c(rmvn_cpp(1, a_beta, B_beta))
    eta_new  = rbeta(1, a_eta, b_eta)
    
    pi_new   = numeric(p_X1)
    mu_new   = numeric(p_X2)
    tau2_new = numeric(p_X2)
    
    if (p_X1 > 0) pi_new = rbeta(p_X1, a_pi, b_pi)
    if (p_X2 > 0) {
      for (q in seq_len(p_X2)) {
        tau2_new[q] = rinvgamma_cpp(1, a_tau2, b_tau2)
        mu_new[q]   = rnorm(1, mean = a_mu, sd = sqrt(tau2_new[q] / b_mu))
      }
    }
    
    betaPars_pot = cbind(betaPars, beta_new)
    sig2Pars_pot = cbind(sig2Pars, sig2_new)
    etaPars_pot  = cbind(etaPars, eta_new)
    piPars_pot   = cbind(piPars, pi_new)
    muPars_pot   = cbind(muPars, mu_new)
    tau2Pars_pot = cbind(tau2Pars, tau2_new)
    
    # ---------------------------------------------------------------------------
    # 3. Calculate conditional lambda_x_mc
    # ---------------------------------------------------------------------------
    if (!is.null(condX) && !is.null(p_condX)) {
      log_lik_cond_X = numeric(K + 1)
      
      for (k in seq_len(K + 1)) {
        ll = 0
        for (i in seq_along(p_condX)) {
          idx = p_condX[i]
          val = condX[i]
          if (idx <= p_X1) {
            ll = ll + dbinom(val, 1, piPars_pot[idx, k], log = TRUE)
          } else {
            q = idx - p_X1
            ll = ll + dnorm(val, muPars_pot[q, k], sqrt(tau2Pars_pot[q, k]), log = TRUE)
          }
        }
        log_lik_cond_X[k] = ll
      }
      
      log_lambda_x_mc = log(c(n_k, alpha)) + log_lik_cond_X
      max_ll = max(log_lambda_x_mc)
      lambda_x_mc = exp(log_lambda_x_mc - max_ll)
    } else {
      lambda_x_mc = c(n_k, alpha)
    }
    
    # ---------------------------------------------------------------------------
    # 4. Draw clusters and update parameters
    # ---------------------------------------------------------------------------
    S_x_mc = vapply(seq_len(M), function(l) rmultinom_cpp(lambda_x_mc), numeric(1))
    orig_unique_S_x_mc = sort(unique(S_x_mc))
    
    w_dp = (c(n_k, alpha) / (N + alpha))[orig_unique_S_x_mc]
    
    updated_pars = update_parameter(
      S_x_mc,
      betaPars_pot, sig2Pars_pot, etaPars_pot,
      piPars_pot, muPars_pot, tau2Pars_pot
    )
    
    S_x_mc = updated_pars$S_x_mc
    unique_S_x_mc = updated_pars$unique_S_x_mc
    K = updated_pars$K_x_mc
    n_k = updated_pars$n_k_mc
    
    betaPars = updated_pars$betaPars_mc
    sig2Pars = updated_pars$sig2Pars_mc
    etaPars  = updated_pars$etaPars_mc
    piPars   = updated_pars$piPars_mc
    muPars   = updated_pars$muPars_mc
    tau2Pars = updated_pars$tau2Pars_mc
    
    sqrt_sig2 = as.vector(sqrt(sig2Pars))
    
    # ---------------------------------------------------------------------------
    # 5. Draw confounders for each of the M observations
    # ---------------------------------------------------------------------------
    X_mc = matrix(nrow = M, ncol = p_X1 + p_X2)
    
    if (!is.null(condX) && !is.null(p_condX)) {
      for (i in seq_along(p_condX)) {
        X_mc[, p_condX[i]] = condX[i]
      }
    }
    
    for (k in seq_len(K)) {
      ind_S_x_mc = which(S_x_mc == unique_S_x_mc[k])
      n_k_temp = length(ind_S_x_mc)
      
      if (n_k_temp > 0) {
        if (p_X1 > 0) {
          for (q in seq_len(p_X1)) {
            if (is.null(p_condX) || !(q %in% p_condX)) {
              X_mc[ind_S_x_mc, q] = rbinom(n_k_temp, 1, piPars[q, k])
            }
          }
        }
        if (p_X2 > 0) {
          for (q in seq_len(p_X2)) {
            col_idx = p_X1 + q
            if (is.null(p_condX) || !(col_idx %in% p_condX)) {
              X_mc[ind_S_x_mc, col_idx] = rnorm(
                n_k_temp, muPars[q, k], sqrt(tau2Pars[q, k])
              )
            }
          }
        }
      }
    }
    
    # ---------------------------------------------------------------------------
    # 6. Calculate Y probabilities
    # ---------------------------------------------------------------------------
    lambda_Y_z0x_mc = matrix(nrow = M, ncol = K)
    lambda_Y_z1x_mc = matrix(nrow = M, ncol = K)
    
    for (k in seq_len(K)) {
      log_w_k_temp = log(w_dp[k])
      
      prob_z0_mc = dbinom(z0, 1, etaPars[, k], log = TRUE)
      prob_z1_mc = dbinom(z1, 1, etaPars[, k], log = TRUE)
      
      prob_x_mc = numeric(M)
      if (p_X1 > 0) {
        for (q in seq_len(p_X1)) {
          prob_x_mc = prob_x_mc + dbinom(X_mc[, q], 1, piPars[q, k], log = TRUE)
        }
      }
      if (p_X2 > 0) {
        for (q in seq_len(p_X2)) {
          prob_x_mc = prob_x_mc + dnorm(
            X_mc[, p_X1 + q], muPars[q, k], sqrt(tau2Pars[q, k]), log = TRUE
          )
        }
      }
      
      lambda_Y_z0x_mc[, k] = log_w_k_temp + prob_z0_mc + prob_x_mc
      lambda_Y_z1x_mc[, k] = log_w_k_temp + prob_z1_mc + prob_x_mc
    }
    
    lambda_Y_z0x_mc = normalize_log_probs(lambda_Y_z0x_mc)
    lambda_Y_z1x_mc = normalize_log_probs(lambda_Y_z1x_mc)
    
    # ---------------------------------------------------------------------------
    # 7. Design matrices
    # ---------------------------------------------------------------------------
    matX_z0x_mc = cbind(1, z0, X_mc) 
    matX_z1x_mc = cbind(1, z1, X_mc) 
    MUy_z0x_mc = (matX_z0x_mc %*% betaPars)
    MUy_z1x_mc = (matX_z1x_mc %*% betaPars)
    Ey_z0x_mc = rowSums(MUy_z0x_mc * lambda_Y_z0x_mc)
    Ey_z1x_mc = rowSums(MUy_z1x_mc * lambda_Y_z1x_mc)
    Ey2_z0x_mc = rowSums((MUy_z0x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = K)) * lambda_Y_z0x_mc)
    Ey2_z1x_mc = rowSums((MUy_z1x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = K)) * lambda_Y_z1x_mc)
    # SDy_z0x_mc = sqrt(pmax(Ey2_z0x_mc - Ey_z0x_mc^{2}, 0))
    # SDy_z1x_mc = sqrt(pmax(Ey2_z1x_mc - Ey_z1x_mc^{2}, 0))
    
    # -----------------------------------------------------------------------------
    lambda_Y_z0x = apply(lambda_Y_z0x_mc, 2, mean)
    lambda_Y_z1x = apply(lambda_Y_z1x_mc, 2, mean)
    MUy_z0x = apply(MUy_z0x_mc, 2, mean)
    MUy_z1x = apply(MUy_z1x_mc, 2, mean)
    Ey_z0x = mean(Ey_z0x_mc)
    Ey_z1x = mean(Ey_z1x_mc)
    Ey2_z0x = mean(Ey2_z0x_mc)
    Ey2_z1x = mean(Ey2_z1x_mc)
    SDy_z0x = sqrt(pmax(Ey2_z0x - Ey_z0x^{2}, 0))
    SDy_z1x = sqrt(pmax(Ey2_z1x - Ey_z1x^{2}, 0))
    
    # ---------------------------------------------------------------------------
    # 8. Sensitivity and QRL optimization loops
    # ---------------------------------------------------------------------------
    for (ii_eta in seq_len(n_eta)) {
      curr_eta_z0 = eta_z0_vec[ii_eta]
      curr_eta_z1 = eta_z1_vec[ii_eta]
      
      for (ii_nu in seq_len(n_nu)) {
        curr_nu = nu_vec[ii_nu]
        curr_log_nu = log(curr_nu)
        safe_curr_log_nu = ifelse(curr_nu == 0, -100, curr_log_nu)
        
        run_z0 = (curr_nu - curr_eta_z0 >= 0)
        run_z1 = (curr_nu - curr_eta_z1 >= 0)
        
        if (run_z0) {
          S_nu_z0x = S_cond_cpp(
            curr_log_nu - curr_eta_z0,
            lambda_Y_z0x_mc, MUy_z0x_mc, sqrt_sig2, TRUE
          )
        }
        if (run_z1) {
          S_nu_z1x = S_cond_cpp(
            curr_log_nu - curr_eta_z1,
            lambda_Y_z1x_mc, MUy_z1x_mc, sqrt_sig2, TRUE
          )
        }
        
        for (ii_rho in seq_len(n_rho)) {
          curr_rho = rho_vec[ii_rho]
          
          if (run_z0) {
            target_rho_nu_z0x = S_nu_z0x + log(1 - curr_rho)
            opt_z0 <- S_optim_cond_cpp(
              target_rho_nu_z0x, Ey_z0x,
              (safe_curr_log_nu - curr_eta_z0) + 1e-10, 100,
              lambda_Y_z0x_mc, MUy_z0x_mc, sqrt_sig2, TRUE
            )
            Y_rho_nu_z0x_local[ii_eta, ii_nu, ii_rho] =
              exp(opt_z0$optimizer + curr_eta_z0) - curr_nu
          }
          
          if (run_z1) {
            target_rho_nu_z1x = S_nu_z1x + log(1 - curr_rho)
            opt_z1 <- S_optim_cond_cpp(
              target_rho_nu_z1x, Ey_z1x,
              (safe_curr_log_nu - curr_eta_z1) + 1e-10, 100,
              lambda_Y_z1x_mc, MUy_z1x_mc, sqrt_sig2, TRUE
            )
            Y_rho_nu_z1x_local[ii_eta, ii_nu, ii_rho] =
              exp(opt_z1$optimizer + curr_eta_z1) - curr_nu
          }
        }
      }
    }
    
    list(z0 = Y_rho_nu_z0x_local, z1 = Y_rho_nu_z1x_local)
  }
  
  # -----------------------------------------------------------------------------
  # Safe worker wrapper
  # -----------------------------------------------------------------------------
  safe_worker_fun = function(post_reps) {
    tryCatch(
      worker_fun(post_reps),
      error = function(e) {
        list(.error = TRUE, .post_reps = post_reps, .msg = conditionMessage(e))
      }
    )
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  start_time <- Sys.time()
  print(start_time)
  
  res_list <- parallel::mclapply(
    X = seq_len(n_MCMC),
    FUN = safe_worker_fun,
    mc.cores = ncores,
    mc.preschedule = FALSE
  )
  
  end_time <- Sys.time()
  print(end_time)
  print(end_time - start_time)
  
  bad_idx <- which(vapply(res_list, function(x) is.list(x) && isTRUE(x$.error), logical(1)))
  if (length(bad_idx) > 0) {
    print(res_list[[bad_idx[1]]])
    stop("At least one worker failed. See printed error above.")
  }
  
  # -----------------------------------------------------------------------------
  # Reassemble arrays
  # -----------------------------------------------------------------------------
  Y_rho_nu_z0x_store = array(NA_real_, dim = c(n_MCMC, n_eta, n_nu, n_rho))
  Y_rho_nu_z1x_store = array(NA_real_, dim = c(n_MCMC, n_eta, n_nu, n_rho))
  
  for (post_reps in seq_len(n_MCMC)) {
    Y_rho_nu_z0x_store[post_reps, , , ] = res_list[[post_reps]]$z0
    Y_rho_nu_z1x_store[post_reps, , , ] = res_list[[post_reps]]$z1
  }
  
  E_rho_nu_diff_store = Y_rho_nu_z1x_store - Y_rho_nu_z0x_store
  
  # -----------------------------------------------------------------------------
  # Summaries
  # -----------------------------------------------------------------------------
  z0_tables_list   = list()
  z1_tables_list   = list()
  diff_tables_list = list()
  counter = 1L
  
  for (ii_eta in seq_len(n_eta)) {
    val_eta_z0 = eta_z0_vec[ii_eta]
    val_eta_z1 = eta_z1_vec[ii_eta]
    
    for (ii_nu in seq_len(n_nu)) {
      val_nu = nu_vec[ii_nu]
      
      for (ii_rho in seq_len(n_rho)) {
        val_rho = rho_vec[ii_rho]
        
        draws_z0   = Y_rho_nu_z0x_store[, ii_eta, ii_nu, ii_rho]
        draws_z1   = Y_rho_nu_z1x_store[, ii_eta, ii_nu, ii_rho]
        draws_diff = E_rho_nu_diff_store[, ii_eta, ii_nu, ii_rho]
        
        if (!all(is.na(draws_z0))) {
          summ_z0 = POSTsummary(draws_z0, esttype, quantile_alpha)
          summ_z0 = cbind(eta_z0 = val_eta_z0, nu = val_nu, rho = val_rho, summ_z0)
          z0_tables_list[[counter]] = summ_z0
        }
        
        if (!all(is.na(draws_z1))) {
          summ_z1 = POSTsummary(draws_z1, esttype, quantile_alpha)
          summ_z1 = cbind(eta_z1 = val_eta_z1, nu = val_nu, rho = val_rho, summ_z1)
          z1_tables_list[[counter]] = summ_z1
        }
        
        if (!all(is.na(draws_diff))) {
          summ_diff = POSTsummary(draws_diff, esttype, quantile_alpha)
          summ_diff = cbind(eta_z0 = val_eta_z0, eta_z1 = val_eta_z1, nu = val_nu, rho = val_rho, summ_diff)
          diff_tables_list[[counter]] = summ_diff
        }
        
        counter = counter + 1L
      }
    }
  }
  
  Y_rho_nu_z0x_result  = if (length(z0_tables_list)   > 0) do.call(rbind, z0_tables_list)   else NULL
  Y_rho_nu_z1x_result  = if (length(z1_tables_list)   > 0) do.call(rbind, z1_tables_list)   else NULL
  E_rho_nu_diff_result = if (length(diff_tables_list) > 0) do.call(rbind, diff_tables_list) else NULL
  
  POSTresult = list(
    Y_rho_nu_z0x_result  = Y_rho_nu_z0x_result,
    Y_rho_nu_z1x_result  = Y_rho_nu_z1x_result,
    E_rho_nu_diff_result = E_rho_nu_diff_result,
    Y_rho_nu_z0x_store   = Y_rho_nu_z0x_store,
    Y_rho_nu_z1x_store   = Y_rho_nu_z1x_store,
    E_rho_nu_diff_store  = E_rho_nu_diff_store
  )
  
  return(POSTresult)
}

SURVIVAL_DPMM_POST_ftn = function(object, M, t_grid, ftn = c("S", "f"),
                                condX = NULL, p_condX = NULL, esttype = "mean") {
  
  # Ensure ftn contains valid options
  ftn = match.arg(ftn, several.ok = TRUE)
  calc_S = "S" %in% ftn
  calc_f = "f" %in% ftn
  
  # -----------------------------------------------------------------------------
  normalize_log_probs = function(log_prob_matrix) {
    max_log_probs = apply(log_prob_matrix, 1, max)
    log_probs_shifted = log_prob_matrix - max_log_probs
    probs_shifted = exp(log_probs_shifted)
    row_sums = rowSums(probs_shifted)
    
    zero_rows = which(row_sums == 0)
    if(length(zero_rows) > 0) {
      probs_shifted[zero_rows,] = 1 / ncol(log_prob_matrix) 
      row_sums[zero_rows] = 1
    }
    return(probs_shifted / row_sums)
  }
  
  # -----------------------------------------------------------------------------
  # DPMM update_parameter function
  # -----------------------------------------------------------------------------
  update_parameter = function(S_x_mc, 
                              betaPars_pot, sig2Pars_pot, etaPars_pot,
                              piPars_pot, muPars_pot, tau2Pars_pot) {
    
    unique_S_x_mc = sort(unique(S_x_mc))
    K_x_mc = length(unique_S_x_mc)
    
    n_k_mc = numeric(K_x_mc)
    for (k in 1:K_x_mc) {
      n_k_mc[k] = sum(S_x_mc == unique_S_x_mc[k])
    }
    
    # Subset matrices to keep only active clusters
    betaPars_mc = matrix(betaPars_pot[, unique_S_x_mc], ncol = K_x_mc)
    sig2Pars_mc = matrix(sig2Pars_pot[, unique_S_x_mc], ncol = K_x_mc)
    etaPars_mc  = matrix(etaPars_pot[, unique_S_x_mc],  ncol = K_x_mc)
    piPars_mc   = matrix(piPars_pot[, unique_S_x_mc],   ncol = K_x_mc)
    muPars_mc   = matrix(muPars_pot[, unique_S_x_mc],   ncol = K_x_mc)
    tau2Pars_mc = matrix(tau2Pars_pot[, unique_S_x_mc], ncol = K_x_mc)
    
    # Remap indices
    S_x_mc_remapped = as.numeric(factor(S_x_mc, levels = unique_S_x_mc))
    
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
  n_t = length(t_grid)
  
  level = 0.05
  quantile_alpha = c(level/2, 1-level/2)
  z0 = 0
  z1 = 1
  
  N = object$constants$N
  # num_MC = object$constants$num_MC
  n_MCMC = object$constants$n_MCMC
  p_matX = object$constants$p_matX
  p_X1 = object$constants$p_X1
  p_X2 = object$constants$p_X2
  
  # Priors for generating the K+1 cluster
  a_beta = object$priors$a_beta; B_beta = object$priors$B_beta
  a_sig2 = object$priors$a_sig2; b_sig2 = object$priors$b_sig2
  a_eta  = object$priors$a_eta;  b_eta  = object$priors$b_eta
  a_pi   = object$priors$a_pi;   b_pi   = object$priors$b_pi
  a_mu   = object$priors$a_mu;   b_mu   = object$priors$b_mu
  a_tau2 = object$priors$a_tau2; b_tau2 = object$priors$b_tau2
  
  iter_check = floor(n_MCMC/10)
  
  # Conditionally allocate storage arrays to save memory
  if (calc_S) {
    S_z0_store = matrix(NA, nrow = n_MCMC, ncol = n_t)
    S_z1_store = matrix(NA, nrow = n_MCMC, ncol = n_t)
  }
  if (calc_f) {
    f_z0_store = matrix(NA, nrow = n_MCMC, ncol = n_t)
    f_z1_store = matrix(NA, nrow = n_MCMC, ncol = n_t)
  }
  
  for (post_reps in 1:n_MCMC) {
    
    # 1. Extract Current MCMC State
    alpha = object$MCMCposteriors$alpha_draws[post_reps]
    K = object$MCMCposteriors$KLists[[post_reps]]
    n_k = object$MCMCposteriors$n_kLists[[post_reps]]
    betaPars = object$MCMCposteriors$betaLists[[post_reps]]
    sig2Pars = object$MCMCposteriors$sig2Lists[[post_reps]]
    etaPars = object$MCMCposteriors$etaLists[[post_reps]]
    piPars = object$MCMCposteriors$piLists[[post_reps]]
    muPars = object$MCMCposteriors$muLists[[post_reps]]
    tau2Pars = object$MCMCposteriors$tau2Lists[[post_reps]]
    
    # ---------------------------------------------------------------------------
    # 2. Pre-generate potential K+1 cluster parameters from prior
    # ---------------------------------------------------------------------------
    sig2_new = c(rscainvchisq_cpp(1, a_sig2, b_sig2))
    beta_new = c(rmvn_cpp(1, a_beta, B_beta)) 
    eta_new  = rbeta(1, a_eta, b_eta)
    
    pi_new   = numeric(p_X1)
    mu_new   = numeric(p_X2)
    tau2_new = numeric(p_X2)
    
    if (p_X1 > 0) pi_new = rbeta(p_X1, a_pi, b_pi)
    if (p_X2 > 0) {
      for (q in 1:p_X2) {
        tau2_new[q] = rinvgamma_cpp(1, a_tau2, b_tau2)
        mu_new[q]   = rnorm(1, mean = a_mu, sd = sqrt(tau2_new[q] / b_mu))
      }
    }
    
    # Form "potential" matrices (Current K + 1 New)
    betaPars_pot = cbind(betaPars, beta_new)
    sig2Pars_pot = cbind(sig2Pars, sig2_new)
    etaPars_pot  = cbind(etaPars, eta_new)
    piPars_pot   = cbind(piPars, pi_new)
    muPars_pot   = cbind(muPars, mu_new)
    tau2Pars_pot = cbind(tau2Pars, tau2_new)
    
    # ---------------------------------------------------------------------------
    # 3. Calculate Conditional lambda_x_mc
    # ---------------------------------------------------------------------------
    if (!is.null(condX) && !is.null(p_condX)) {
      log_lik_cond_X = numeric(K + 1)
      
      # Calculate log-likelihood of condX for all K+1 clusters
      for (k in 1:(K + 1)) {
        ll = 0
        for (i in seq_along(p_condX)) {
          idx = p_condX[i]
          val = condX[i]
          if (idx <= p_X1) {
            ll = ll + dbinom(val, 1, piPars_pot[idx, k], log = TRUE)
          } else {
            q = idx - p_X1
            ll = ll + dnorm(val, muPars_pot[q, k], sqrt(tau2Pars_pot[q, k]), log = TRUE)
          }
        }
        log_lik_cond_X[k] = ll
      }
      
      # Bayes Rule: log(Prior Weight) + log(Likelihood)
      log_lambda_x_mc = log(c(n_k, alpha)) + log_lik_cond_X
      
      # Exponentiate safely to prevent numerical underflow
      max_ll = max(log_lambda_x_mc)
      lambda_x_mc = exp(log_lambda_x_mc - max_ll)
    } else {
      # Marginal mode: purely DP prior weights
      lambda_x_mc = c(n_k, alpha)
    }
    
    # ---------------------------------------------------------------------------
    # 4. Draw clusters and update parameters
    # ---------------------------------------------------------------------------
    S_x_mc = sapply(1:M, function(l) rmultinom_cpp(lambda_x_mc))
    orig_unique_S_x_mc = sort(unique(S_x_mc))
    
    # DP predictive weights for the chosen clusters
    w_dp = (c(n_k, alpha) / (N + alpha))[orig_unique_S_x_mc]
    
    # Subsets the _pot matrices down to only the clusters that were drawn
    updated_pars = update_parameter(S_x_mc, 
                                    betaPars_pot, sig2Pars_pot, etaPars_pot, 
                                    piPars_pot, muPars_pot, tau2Pars_pot)
    
    S_x_mc = updated_pars$S_x_mc
    unique_S_x_mc = updated_pars$unique_S_x_mc
    K = updated_pars$K_x_mc
    n_k = updated_pars$n_k_mc
    
    betaPars = updated_pars$betaPars_mc
    sig2Pars = updated_pars$sig2Pars_mc
    etaPars  = updated_pars$etaPars_mc
    piPars   = updated_pars$piPars_mc
    muPars   = updated_pars$muPars_mc
    tau2Pars = updated_pars$tau2Pars_mc
    
    # -----------------------------------------------------------------------------
    # 5. Draw confounders for each of the M observations
    # -----------------------------------------------------------------------------
    X_mc = matrix(nrow = M, ncol = (p_X1 + p_X2))
    
    if (!is.null(condX) && !is.null(p_condX)) {
      for (i in seq_along(p_condX)) {
        X_mc[, p_condX[i]] = condX[i]
      }
    }
    
    for (k in 1:K) {
      ind_S_x_mc = which(S_x_mc == unique_S_x_mc[k])
      n_k_temp = length(ind_S_x_mc)
      
      if (n_k_temp > 0) {
        if (p_X1 > 0) {
          for (q in 1:p_X1) {
            if (is.null(p_condX) || !(q %in% p_condX)) {
              X_mc[ind_S_x_mc, q] = rbinom(n_k_temp, 1, piPars[q,k])
            }
          }
        }
        if (p_X2 > 0) {
          for (q in 1:p_X2) {
            col_idx = p_X1 + q
            if (is.null(p_condX) || !(col_idx %in% p_condX)) {
              X_mc[ind_S_x_mc, col_idx] = rnorm(n_k_temp, muPars[q,k], sqrt(tau2Pars[q,k]))
            }
          }
        }
      }
    }
    
    # -----------------------------------------------------------------------------
    # 6. Calculate Y probabilities 
    # -----------------------------------------------------------------------------
    lambda_Y_z0x_mc = matrix(nrow = M, ncol = K)
    lambda_Y_z1x_mc = matrix(nrow = M, ncol = K)
    
    for (k in 1:K) {
      log_w_k_temp = log(w_dp[k])
      
      prob_z0_mc = dbinom(z0, 1, etaPars[,k], log = TRUE)
      prob_z1_mc = dbinom(z1, 1, etaPars[,k], log = TRUE)
      
      prob_x_mc = numeric(M)
      if (p_X1 > 0) {
        for (q in 1:p_X1) {
          prob_x_mc = prob_x_mc + dbinom(X_mc[, q], 1, piPars[q,k], log = TRUE)
        }
      }
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
    # 7. Create the design matrices for all M subjects
    # -----------------------------------------------------------------------------
    # Create the design matrices for all M subjects
    matX_z0x_mc = cbind(1, z0, X_mc) 
    matX_z1x_mc = cbind(1, z1, X_mc) 
    MUy_z0x_mc = (matX_z0x_mc %*% betaPars)
    MUy_z1x_mc = (matX_z1x_mc %*% betaPars)
    Ey_z0x_mc = rowSums(MUy_z0x_mc * lambda_Y_z0x_mc)
    Ey_z1x_mc = rowSums(MUy_z1x_mc * lambda_Y_z1x_mc)
    Ey2_z0x_mc = rowSums((MUy_z0x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = Ky)) * lambda_Y_z0x_mc)
    Ey2_z1x_mc = rowSums((MUy_z1x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = Ky)) * lambda_Y_z1x_mc)
    # SDy_z0x_mc = sqrt(pmax(Ey2_z0x_mc - Ey_z0x_mc^{2}, 0))
    # SDy_z1x_mc = sqrt(pmax(Ey2_z1x_mc - Ey_z1x_mc^{2}, 0))
    
    # -----------------------------------------------------------------------------
    lambda_Y_z0x = apply(lambda_Y_z0x_mc, 2, mean)
    lambda_Y_z1x = apply(lambda_Y_z1x_mc, 2, mean)
    MUy_z0x = apply(MUy_z0x_mc, 2, mean)
    MUy_z1x = apply(MUy_z1x_mc, 2, mean)
    Ey_z0x = mean(Ey_z0x_mc)
    Ey_z1x = mean(Ey_z1x_mc)
    Ey2_z0x = mean(Ey2_z0x_mc)
    Ey2_z1x = mean(Ey2_z1x_mc)
    SDy_z0x = sqrt(pmax(Ey2_z0x - Ey_z0x^{2}, 0))
    SDy_z1x = sqrt(pmax(Ey2_z1x - Ey_z1x^{2}, 0))
    
    # -----------------------------------------------------------------------------
    # SURVIVAL & DENSITY EVALUATION LOOP
    # -----------------------------------------------------------------------------
    for (ii_t in 1:n_t) {
      curr_t = t_grid[ii_t]
      
      if (curr_t <= 0) {
        if (calc_S) {
          S_z0_store[post_reps, ii_t] = 1
          S_z1_store[post_reps, ii_t] = 1
        }
        if (calc_f) {
          f_z0_store[post_reps, ii_t] = 0
          f_z1_store[post_reps, ii_t] = 0
        }
      } else {
        curr_log_t = log(curr_t)
        
        if (calc_S) {
          log_S_z0 = S_cond_cpp(curr_log_t, lambda_Y_z0x_mc, MUy_z0x_mc, sqrt(sig2Pars), TRUE)
          log_S_z1 = S_cond_cpp(curr_log_t, lambda_Y_z1x_mc, MUy_z1x_mc, sqrt(sig2Pars), TRUE)
          
          S_z0_store[post_reps, ii_t] = exp(log_S_z0)
          S_z1_store[post_reps, ii_t] = exp(log_S_z1)
        }
        
        if (calc_f) {
          log_f_z0 = f_cond_cpp(curr_log_t, lambda_Y_z0x_mc, MUy_z0x_mc, sqrt(sig2Pars), TRUE)
          log_f_z1 = f_cond_cpp(curr_log_t, lambda_Y_z1x_mc, MUy_z1x_mc, sqrt(sig2Pars), TRUE)
          
          # Jacobian applied: exp(log_f) / t yields density of survival time T
          f_z0_store[post_reps, ii_t] = exp(log_f_z0) / curr_t
          f_z1_store[post_reps, ii_t] = exp(log_f_z1) / curr_t
        }
      }
    }
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing Curve Evaluation",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  
  # -----------------------------------------------------------------------------
  # Combine outputs conditionally
  # -----------------------------------------------------------------------------
  result_list = list()
  
  if (calc_S) {
    z0_S_list = list()
    z1_S_list = list()
    for(ii_t in 1:n_t) {
      val_t = t_grid[ii_t]
      z0_S_list[[ii_t]] = cbind(t = val_t, POSTsummary(S_z0_store[, ii_t], esttype, quantile_alpha))
      z1_S_list[[ii_t]] = cbind(t = val_t, POSTsummary(S_z1_store[, ii_t], esttype, quantile_alpha))
    }
    result_list$S_z0_result = do.call(rbind, z0_S_list)
    result_list$S_z1_result = do.call(rbind, z1_S_list)
    result_list$S_z0_store  = S_z0_store
    result_list$S_z1_store  = S_z1_store
  }
  
  if (calc_f) {
    z0_f_list = list()
    z1_f_list = list()
    for(ii_t in 1:n_t) {
      val_t = t_grid[ii_t]
      z0_f_list[[ii_t]] = cbind(t = val_t, POSTsummary(f_z0_store[, ii_t], esttype, quantile_alpha))
      z1_f_list[[ii_t]] = cbind(t = val_t, POSTsummary(f_z1_store[, ii_t], esttype, quantile_alpha))
    }
    result_list$f_z0_result = do.call(rbind, z0_f_list)
    result_list$f_z1_result = do.call(rbind, z1_f_list)
    result_list$f_z0_store  = f_z0_store
    result_list$f_z1_store  = f_z1_store
  }
  
  return(result_list)
}

# ------------------------------------------------------------------------
SURVIVAL_DPMM = function(object, M = 2e3, nu = 0, rho = 0.5, 
                         gibbs_iter = 2e4, gibbs_burnin = 2e4, gibbs_thin = 1e2, 
                         esttype = "mean", saveall = FALSE){
  # object = temp_data
  MCMCresult = SURVIVAL_DPMM_MCMC(object, gibbs_iter, gibbs_burnin, gibbs_thin)
  POSTresult = SURVIVAL_DPMM_POST(MCMCresult, M, nu, rho, 
                                  esttype = esttype)
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
  
  # --- Identify covariate types (binary vs. continuous) ---
  is_binary <- apply(Xobs, 2, function(col) {
    col2 <- col[!is.na(col)]
    if (is.logical(col2)) col2 <- as.integer(col2)
    length(col2) > 0 && all(col2 %in% c(0, 1))
  })
  
  X_bin <- Xobs[, is_binary, drop = FALSE]
  X_con <- Xobs[, !is_binary, drop = FALSE]
  p_X1 <- ncol(X_bin)
  p_X2 <- ncol(X_con)
  p_X = p_X1 + p_X2
  p_Z = length(Zobs)/N
  Zobs = matrix(Zobs, ncol = p_Z)
  
  # --- Missingness indicators (fixed throughout) ---
  # Define helper function globally for safety
  build_sid_count <- function(Sy, Sx, unique_Sy) {
    N <- length(Sy)
    sid <- integer(N)
    count <- 1L

    for (k in unique_Sy) {
      ind_Sy <- (Sy == k)
      unique_Sx_temp <- sort(unique(Sx[ind_Sy]))
      for (r in unique_Sx_temp) {
        ind_kr <- ind_Sy & (Sx == r)
        sid[ind_kr] <- count
        count <- count + 1L
      }
    }

    if (any(sid == 0L)) stop("sid has zeros: some (Sy,Sx) pairs were not mapped.")
    list(sid = sid, Kyx = count - 1L)
  }
  miss_bin <- if (p_X1 > 0) is.na(X_bin) else matrix(FALSE, nrow(Xobs), 0)
  miss_con <- if (p_X2 > 0) is.na(X_con) else matrix(FALSE, nrow(Xobs), 0)
  has_missing_X <- (any(miss_bin) || any(miss_con))

  # --- Initial deterministic imputation (starting values only) ---
  if (p_X1 > 0) {
    for (j in seq_len(p_X1)) {
      idx <- which(miss_bin[, j])
      if (!length(idx)) next
      mj <- mean(X_bin[, j], na.rm = TRUE)
      if (is.nan(mj)) mj <- 0
      X_bin[idx, j] <- as.integer(mj >= 0.5)
    }
  }

  if (p_X2 > 0) {
    for (j in seq_len(p_X2)) {
      idx <- which(miss_con[, j])
      if (!length(idx)) next
      mj <- mean(X_con[, j], na.rm = TRUE)
      if (is.nan(mj)) mj <- 0
      X_con[idx, j] <- mj
    }
  }

  # --- Precompute scaling ONCE (fixed for entire run) ---
  con_center <- con_scale <- NULL
  if (p_X2 > 0) {
    con_center <- colMeans(X_con, na.rm = TRUE)
    con_scale  <- apply(X_con, 2, sd, na.rm = TRUE)
    con_scale[!is.finite(con_scale) | con_scale <= 0] <- 1
    scaX_con <- sweep(X_con, 2, con_center, "-")
    scaX_con <- sweep(scaX_con, 2, con_scale, "/")
  }

  # --- Build scaled covariate matrix used in the model ---
  if (p_X1 > 0 && p_X2 > 0) {
    scaXobs <- cbind(X_bin, scaX_con)
  } else if (p_X2 > 0) {
    scaXobs <- scaX_con
  } else {
    scaXobs <- X_bin
  }
  
  matX <- cbind(1, Zobs, scaXobs)
  colnames(matX) <- c("int", "Z", paste0("X", seq_len(ncol(scaXobs))))
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
  B_beta = vcov(survreg_beta)[1:p_beta, 1:p_beta, drop = FALSE]
  c_beta = N
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
  # c_beta = N
  # B_beta = c_beta * B_beta
  # Binv_beta = inv_cpp(B_beta)
  # aBinv_beta = Binv_beta %*% a_beta
  # a_sig2 = 1
  # b_sig2 = var(residuals_aft[Dobs==1]) * (a_sig2 + 2)
  a_sig2 = 3
  b_sig2 = 0.1
  
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
  # Ysatr = Y_data_augmentation(Tobs, Dobs, matX, Sy,
  #                             a_beta, Binv_beta, aBinv_beta,
  #                             a_sig2, b_sig2,
  #                             inv_cpp, rtrunct_cpp,
  #                             upper = Inf)
  Ysatr = Y_data_augmentation(Tobs, Dobs, matX, Sy, betaPars, sig2Pars)
  
  for (k in unique_Sy) {
    ind_S = which(Sy == k)
    n_temp = length(ind_S)
    T_temp = Ysatr[ind_S]
    matX_temp = matrix(matX[ind_S,], ncol = p_matX)
    
    # update parameters for true failure time
    # regression_results = update_reg_con(T_temp, matX_temp, betaPars[,k],
    #                                     a_beta, Binv_beta, aBinv_beta, a_sig2, b_sig2)
    regression_results = update_reg_con(T_temp, matX_temp, a_beta, Binv_beta, a_sig2, b_sig2)
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
  etaPars = matrix(nrow = p_Z, ncol = Kyx)
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
    f0_zx_all[, q] = beta(a_eta + Zobs[,q], b_eta - Zobs[,q] + 1) # beta(a_pi,b_pi) = 1 in this case
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
  betaLists = vector("list", n_MCMC)
  sig2Lists = vector("list", n_MCMC)
  etaLists  = vector("list", n_MCMC)
  piLists   = vector("list", n_MCMC)
  muLists   = vector("list", n_MCMC)
  tau2Lists = vector("list", n_MCMC)
  
  n_kLists       = vector("list", n_MCMC)
  n_rkLists      = vector("list", n_MCMC)
  KyLists        = vector("list", n_MCMC)
  KyxLists       = vector("list", n_MCMC)
  max_Kx_SyLists = vector("list", n_MCMC)
  
  alpha_theta_draws = numeric(n_MCMC)
  alpha_omega_draws = numeric(n_MCMC)
  
  count_it = 1
  # End initial values ------------------------------------------------------------
  
  # First draw each parameter for BNP model. Then calculate causal effect.
  for (gibbs_reps in 1:gibbs_total) {
    # Update of all parameters in BNP model--------------------------------
    # Update cluster membership -------------------------------------------
    scaXobs <- matX[, -(1:(1 + p_Z)), drop = FALSE]
    cluster_res = EDPMMcluster_cpp(N, p_X1, p_X2, 
                                   Ysatr, matX, Zobs, scaXobs,
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
    # Ysatr = Y_data_augmentation(Tobs, Dobs, matX, Sy,
    #                             a_beta, Binv_beta, aBinv_beta,
    #                             a_sig2, b_sig2,
    #                             inv_cpp, rtrunct_cpp,
    #                             upper = Inf)
    Ysatr = Y_data_augmentation(Tobs, Dobs, matX, Sy, betaPars, sig2Pars)
    
    # If no missing covariates, do nothing
    if (has_missing_X) {
      outX <- X_data_augmentation(
        matX = matX,
        Ysatr = Ysatr,
        Zobs = Zobs,
        Sy = Sy, Sx = Sx, unique_Sy = unique_Sy,
        betaPars = betaPars, sig2Pars = sig2Pars,
        piPars = piPars, muPars = muPars, tau2Pars = tau2Pars,
        miss_bin = miss_bin, miss_con = miss_con,
        p_Z = p_Z, p_X1 = p_X1, p_X2 = p_X2,
        build_sid_count = build_sid_count,
        X_bin_data_augmentation = X_bin_data_augmentation,
        X_con_data_augmentation = X_con_data_augmentation,
        loglik_fun = loglik_fun,
        use_Y = FALSE   # <-- per your request: no Y-likelihood for now
      )
      matX   <- outX$matX
      scaXobs <- outX$scaXobs
    }
    
    for (k in unique_Sy) {
      ind_S = which(Sy == k)
      n_temp = length(ind_S)
      T_temp = Ysatr[ind_S]
      matX_temp = matrix(matX[ind_S,], ncol = p_matX)
      
      # update parameters for true failure time
      # regression_results = update_reg_con(T_temp, matX_temp, betaPars[,k],
      #                                     a_beta, Binv_beta, aBinv_beta, a_sig2, b_sig2)
      regression_results = update_reg_con(T_temp, matX_temp, a_beta, Binv_beta, a_sig2, b_sig2)
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
      if ((gibbs_reps - gibbs_burnin) %% gibbs_thin == 0) {
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
                   p_matX = p_matX, p_X1 = p_X1, p_X2 = p_X2,
                   con_center = con_center, con_scale = con_scale)
  
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
  
  return(MCMCresult)
}

# -----------------------------------------------------------------------------
# Define function for the post-processing under the BNP (EDPMM) model
SURVIVAL_EDPMM_POST = function(object, M, nu, rho, eta_z0 = NULL, eta_z1 = NULL, 
                               condX = NULL, p_condX = NULL, esttype = "mean"){
  # -----------------------------------------------------------------------------
  # object = EDPMM_results_MCMC
  # condX = vector of values to condition on (e.g., c(1)), NULL for marginal
  # p_condX = vector of column indices in the X matrix for condX, NULL for marginal
  # eta_z0, eta_z1 = Vectors of sensitivity shift parameters
  # -----------------------------------------------------------------------------
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
      probs_shifted[zero_rows,] = 1 / ncol(log_prob_matrix) 
      row_sums[zero_rows] = 1
    }
    return(probs_shifted / row_sums)
  }
  
  # -----------------------------------------------------------------------------
  # Extract Constants and Settings
  nu_vec = as.vector(nu)
  rho_vec = as.vector(rho)
  n_nu = length(nu_vec)
  n_rho = length(rho_vec)
  
  # Build a full grid of all combinations of eta_z0 and eta_z1
  eta_grid = expand.grid(
    eta_z0 = if(is.null(eta_z0)) 0 else as.vector(eta_z0),
    eta_z1 = if(is.null(eta_z1)) 0 else as.vector(eta_z1)
  )
  eta_z0_vec = eta_grid$eta_z0
  eta_z1_vec = eta_grid$eta_z1
  n_eta = nrow(eta_grid)
  
  level = 0.05
  quantile_alpha = c(level/2,1-level/2)
  z0 = 0
  z1 = 1
  
  N = object$constants$N
  # num_MC = object$constants$num_MC
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
  iter_check = floor(n_MCMC/10)
  
  # Set up the 4D Storage Arrays for the full Sensitivity Grid
  Y_rho_nu_z0x_store = array(NA, dim = c(n_MCMC, n_eta, n_nu, n_rho))
  Y_rho_nu_z1x_store = array(NA, dim = c(n_MCMC, n_eta, n_nu, n_rho))
  for (post_reps in 1:n_MCMC) {
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
    # Calculate Conditional lambda_y_mc and S_x Probabilities
    # -----------------------------------------------------------------------------
    if (!is.null(condX) && !is.null(p_condX)) {
      log_lambda_y_mc = numeric(Ky)
      log_prob_Sx_list = list()
      
      count_lik = 1
      for (k in 1:Ky) {
        max_Kx_Sy_temp = max_Kx_Sy[k]
        w_rk_temp = n_rk[k, 1:max_Kx_Sy_temp]
        
        # Protect against empty empirical clusters
        sum_w = sum(w_rk_temp)
        if (sum_w == 0) sum_w = 1 
        log_w_rk_temp = log(w_rk_temp / sum_w)
        
        log_lik_r = numeric(max_Kx_Sy_temp)
        for (r in 1:max_Kx_Sy_temp) {
          ll = 0
          for (i in seq_along(p_condX)) {
            idx = p_condX[i]
            val = condX[i]
            if (idx <= p_X1) {
              ll = ll + dbinom(val, 1, piPars[idx, count_lik], log = TRUE)
            } else {
              q = idx - p_X1
              ll = ll + dnorm(val, muPars[q, count_lik], sqrt(tau2Pars[q, count_lik]), log = TRUE)
            }
          }
          log_lik_r[r] = log_w_rk_temp[r] + ll
          count_lik = count_lik + 1
        }
        
        # Marginal likelihood of condX given S_y = k
        log_lambda_y_mc[k] = log(n_k[k]) + logSumExp(log_lik_r)
        
        # Store the posterior probabilities of S_x given S_y = k and condX
        log_prob_Sx_list[[k]] = log_lik_r
      }
      
      max_ll_y = max(log_lambda_y_mc)
      lambda_y_mc = exp(log_lambda_y_mc - max_ll_y)
    } else {
      lambda_y_mc = n_k 
      log_prob_Sx_list = NULL
    }
    
    # -----------------------------------------------------------------------------
    # Draw new clusters for x
    S_y_mc = sapply(1:M, function(l) rmultinom_cpp(lambda_y_mc))
    unique_S_y_mc = unique(S_y_mc)
    Ky_mc = max(unique_S_y_mc)
    
    # -----------------------------------------------------------------------------
    # Draw confounders for each of the M observations
    S_x_mc = numeric(M)
    n_k_mc  = numeric(Ky)
    n_rk_mc = matrix(0, nrow = nrow(n_rk), ncol = ncol(n_rk))
    max_Kx_Sy_mc = numeric(Ky)
    
    X_mc = matrix(nrow = M, ncol = (p_X1 + p_X2))
    
    # Pre-assign condX to X_mc
    if (!is.null(condX) && !is.null(p_condX)) {
      for (i in seq_along(p_condX)) {
        X_mc[, p_condX[i]] = condX[i]
      }
    }
    
    count = 1
    for (k in 1:Ky) {
      ind_S_y_mc = which(S_y_mc == k)
      n_k_temp = length(ind_S_y_mc)
      n_k_mc[k] = n_k_temp
      max_Kx_Sy_temp = max_Kx_Sy[k]
      
      if (n_k_temp > 0) {
        if (max_Kx_Sy_temp == 1) {
          S_x_sub = rep(1, n_k_temp)
        } else {
          if (!is.null(condX) && !is.null(p_condX)) {
            # Use inner-conditional probabilities for S_x given S_y
            log_p = log_prob_Sx_list[[k]]
            p_Sx = exp(log_p - max(log_p))
            S_x_sub = sample(1:max_Kx_Sy_temp, n_k_temp, replace = TRUE, prob = p_Sx)
          } else {
            S_x_sub = sample(1:max_Kx_Sy_temp, n_k_temp, replace = TRUE, prob = n_rk[k,1:max_Kx_Sy_temp])
          }
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
                # Skip sampling if variable is conditioned
                if (is.null(p_condX) || !(q %in% p_condX)) {
                  X_mc[ind_Sy_x_mc, q] = rbinom(num_Sy_x_mc, 1, piPars[q, count])
                }
              }
            }
            if (p_X2 > 0) {
              for (q in 1:p_X2) {
                col_idx = p_X1 + q
                # Skip sampling if variable is conditioned
                if (is.null(p_condX) || !(col_idx %in% p_condX)) {
                  X_mc[ind_Sy_x_mc, col_idx] = rnorm(num_Sy_x_mc, muPars[q, count], sqrt(tau2Pars[q, count]))
                }
              }
            }
          }
          count = count + 1
        }
      } else {
        # Advance 'count' to stay aligned with arrays even if cluster is empty
        count = count + max_Kx_Sy_temp
      }
    }
    
    # -----------------------------------------------------------------------------
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
    # Create the design matrices for all M subjects
    matX_z0x_mc = cbind(1, z0, X_mc) 
    matX_z1x_mc = cbind(1, z1, X_mc) 
    MUy_z0x_mc = (matX_z0x_mc %*% betaPars)
    MUy_z1x_mc = (matX_z1x_mc %*% betaPars)
    Ey_z0x_mc = rowSums(MUy_z0x_mc * lambda_Y_z0x_mc)
    Ey_z1x_mc = rowSums(MUy_z1x_mc * lambda_Y_z1x_mc)
    Ey2_z0x_mc = rowSums((MUy_z0x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = Ky)) * lambda_Y_z0x_mc)
    Ey2_z1x_mc = rowSums((MUy_z1x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = Ky)) * lambda_Y_z1x_mc)
    # SDy_z0x_mc = sqrt(pmax(Ey2_z0x_mc - Ey_z0x_mc^{2}, 0))
    # SDy_z1x_mc = sqrt(pmax(Ey2_z1x_mc - Ey_z1x_mc^{2}, 0))
    
    # -----------------------------------------------------------------------------
    lambda_Y_z0x = apply(lambda_Y_z0x_mc, 2, mean)
    lambda_Y_z1x = apply(lambda_Y_z1x_mc, 2, mean)
    MUy_z0x = apply(MUy_z0x_mc, 2, mean)
    MUy_z1x = apply(MUy_z1x_mc, 2, mean)
    Ey_z0x = mean(Ey_z0x_mc)
    Ey_z1x = mean(Ey_z1x_mc)
    Ey2_z0x = mean(Ey2_z0x_mc)
    Ey2_z1x = mean(Ey2_z1x_mc)
    SDy_z0x = sqrt(pmax(Ey2_z0x - Ey_z0x^{2}, 0))
    SDy_z1x = sqrt(pmax(Ey2_z1x - Ey_z1x^{2}, 0))
    
    # -----------------------------------------------------------------------------
    # SENSITIVITY AND OPTIMIZATION LOOPS
    # -----------------------------------------------------------------------------
    for (ii_eta in 1:n_eta) {
      curr_eta_z0 = eta_z0_vec[ii_eta]
      curr_eta_z1 = eta_z1_vec[ii_eta]
      
      for (ii_nu in 1:n_nu) {
        curr_nu = nu_vec[ii_nu]
        curr_log_nu = log(curr_nu)
        safe_curr_log_nu = ifelse(curr_nu == 0, -100, curr_log_nu)
        
        # Condition to run optimization
        run_z0 = (curr_nu - curr_eta_z0 >= 0)
        run_z1 = (curr_nu - curr_eta_z1 >= 0)
        
        # 1. Denominator: Use the CONDITIONAL C++ functions that accept matrices
        if (run_z0) {
          S_nu_z0x = S_cond_cpp(safe_curr_log_nu - curr_eta_z0,
                                lambda_Y_z0x_mc, MUy_z0x_mc, as.vector(sqrt(sig2Pars)), TRUE)
        }
        if (run_z1) {
          S_nu_z1x = S_cond_cpp(safe_curr_log_nu - curr_eta_z1,
                                lambda_Y_z1x_mc, MUy_z1x_mc, as.vector(sqrt(sig2Pars)), TRUE)
        }
        
        for(ii_rho in 1:n_rho) {
          curr_rho = rho_vec[ii_rho]
          
          if (run_z0) {
            target_rho_nu_z0x = S_nu_z0x + log(1 - curr_rho)
            # 2. Optimize: Use S_optim_cond_cpp
            opt_z0 <- S_optim_cond_cpp(target = target_rho_nu_z0x,
                                       y_init = safe_curr_log_nu - curr_eta_z0 + 0.1,
                                       y_min = safe_curr_log_nu - curr_eta_z0,
                                       y_max = 100,
                                       lambda_mat = lambda_Y_z0x_mc,
                                       mu_mat = MUy_z0x_mc,
                                       sig = as.vector(sqrt(sig2Pars)),
                                       logt = TRUE)
            # 3. Recover residual time y (Your math here is perfectly correct)
            Y_rho_nu_z0x_store[post_reps, ii_eta, ii_nu, ii_rho] = exp(opt_z0$optimizer + curr_eta_z0) - curr_nu
          }
          if (run_z1) {
            target_rho_nu_z1x = S_nu_z1x + log(1 - curr_rho)
            # 2. Optimize: Use S_optim_cond_cpp
            opt_z1 <- S_optim_cond_cpp(target = target_rho_nu_z1x,
                                       y_init = safe_curr_log_nu - curr_eta_z1 + 0.1,
                                       y_min = safe_curr_log_nu - curr_eta_z1,
                                       y_max = 100,
                                       lambda_mat = lambda_Y_z1x_mc,
                                       mu_mat = MUy_z1x_mc,
                                       sig = as.vector(sqrt(sig2Pars)),
                                       logt = TRUE)
            # 3. Recover residual time y (Your math here is perfectly correct)
            Y_rho_nu_z1x_store[post_reps, ii_eta, ii_nu, ii_rho] = exp(opt_z1$optimizer + curr_eta_z1) - curr_nu
          }
        } 
      } 
    }
    
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
  
  for(ii_eta in 1:n_eta) {
    val_eta_z0 = eta_z0_vec[ii_eta]
    val_eta_z1 = eta_z1_vec[ii_eta]
    
    for(ii_nu in 1:n_nu) {
      val_nu  = nu_vec[ii_nu]
      
      for(ii_rho in 1:n_rho) {
        val_rho = rho_vec[ii_rho]
        
        # 1. Extract draws from the new 4D arrays
        draws_z0   = Y_rho_nu_z0x_store[, ii_eta, ii_nu, ii_rho]
        draws_z1   = Y_rho_nu_z1x_store[, ii_eta, ii_nu, ii_rho]
        draws_diff = E_rho_nu_diff_store[, ii_eta, ii_nu, ii_rho]
        
        # 2. Safely summarize only if the optimization was run (not all NAs)
        if (!all(is.na(draws_z0))) {
          summ_z0 = POSTsummary(draws_z0, esttype, quantile_alpha)
          summ_z0 = cbind(eta_z0 = val_eta_z0, nu = val_nu, rho = val_rho, summ_z0)
          z0_tables_list[[counter]] = summ_z0
        }
        
        if (!all(is.na(draws_z1))) {
          summ_z1 = POSTsummary(draws_z1, esttype, quantile_alpha)
          summ_z1 = cbind(eta_z1 = val_eta_z1, nu = val_nu, rho = val_rho, summ_z1)
          z1_tables_list[[counter]] = summ_z1
        }
        
        if (!all(is.na(draws_diff))) {
          summ_diff = POSTsummary(draws_diff, esttype, quantile_alpha)
          summ_diff = cbind(eta_z0 = val_eta_z0, eta_z1 = val_eta_z1, nu = val_nu, rho = val_rho, summ_diff)
          diff_tables_list[[counter]] = summ_diff
        }
        
        counter = counter + 1
      }
    }
  }
  
  # 3. do.call(rbind) will automatically ignore any NULL elements in the list
  Y_rho_nu_z0x_result  = if(length(z0_tables_list) > 0) do.call(rbind, z0_tables_list) else NULL
  Y_rho_nu_z1x_result  = if(length(z1_tables_list) > 0) do.call(rbind, z1_tables_list) else NULL
  E_rho_nu_diff_result = if(length(diff_tables_list) > 0) do.call(rbind, diff_tables_list) else NULL
  
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

SURVIVAL_EDPMM_POST_parallel = function(object, M, nu, rho,
                                        eta_z0 = NULL, eta_z1 = NULL,
                                        condX = NULL, p_condX = NULL,
                                        esttype = "mean",
                                        ncores = max(1L, parallel::detectCores() - 1L),
                                        seed = NULL) {
  # -----------------------------------------------------------------------------
  # Parallel version of SURVIVAL_EDPMM_POST
  # Best on macOS/Linux via mclapply()
  # -----------------------------------------------------------------------------
  
  if (.Platform$OS.type == "windows") {
    stop("This version uses parallel::mclapply() and is intended for macOS/Linux.\n",
         "If you want, I can give you a PSOCK / parLapply version for Windows.")
  }
  
  logSumExp = function(vals) {
    max_val = max(vals, na.rm = TRUE)
    if (is.na(max_val) || max_val == -Inf) return(-Inf)
    max_val + log(sum(exp(vals - max_val)))
  }
  
  normalize_log_probs = function(log_prob_matrix) {
    max_log_probs = apply(log_prob_matrix, 1, max)
    log_probs_shifted = log_prob_matrix - max_log_probs
    probs_shifted = exp(log_probs_shifted)
    row_sums = rowSums(probs_shifted)
    
    zero_rows = which(row_sums == 0)
    if (length(zero_rows) > 0) {
      probs_shifted[zero_rows, ] = 1 / ncol(log_prob_matrix)
      row_sums[zero_rows] = 1
    }
    probs_shifted / row_sums
  }
  
  # -----------------------------------------------------------------------------
  # Constants
  # -----------------------------------------------------------------------------
  nu_vec = as.vector(nu)
  rho_vec = as.vector(rho)
  n_nu = length(nu_vec)
  n_rho = length(rho_vec)
  
  eta_grid = expand.grid(
    eta_z0 = if (is.null(eta_z0)) 0 else as.vector(eta_z0),
    eta_z1 = if (is.null(eta_z1)) 0 else as.vector(eta_z1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  eta_z0_vec = eta_grid$eta_z0
  eta_z1_vec = eta_grid$eta_z1
  n_eta = nrow(eta_grid)
  
  level = 0.05
  quantile_alpha = c(level / 2, 1 - level / 2)
  z0 = 0
  z1 = 1
  
  N = object$constants$N
  # num_MC = object$constants$num_MC
  n_MCMC = object$constants$n_MCMC
  p_matX = object$constants$p_matX
  p_X1 = object$constants$p_X1
  p_X2 = object$constants$p_X2
  
  # -----------------------------------------------------------------------------
  # Worker for a single posterior draw
  # -----------------------------------------------------------------------------
  worker_fun = function(post_reps) {
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
    
    sqrt_sig2 = sqrt(sig2Pars)
    
    Y_rho_nu_z0x_local = array(NA_real_, dim = c(n_eta, n_nu, n_rho))
    Y_rho_nu_z1x_local = array(NA_real_, dim = c(n_eta, n_nu, n_rho))
    
    # ---------------------------------------------------------------------------
    # Calculate lambda_y_mc and optional inner conditional probs for S_x
    # ---------------------------------------------------------------------------
    if (!is.null(condX) && !is.null(p_condX)) {
      log_lambda_y_mc = numeric(Ky)
      log_prob_Sx_list = vector("list", Ky)
      
      count_lik = 1L
      for (k in seq_len(Ky)) {
        max_Kx_Sy_temp = max_Kx_Sy[k]
        w_rk_temp = n_rk[k, seq_len(max_Kx_Sy_temp)]
        
        sum_w = sum(w_rk_temp)
        if (sum_w == 0) sum_w = 1
        log_w_rk_temp = log(w_rk_temp / sum_w)
        
        log_lik_r = numeric(max_Kx_Sy_temp)
        for (r in seq_len(max_Kx_Sy_temp)) {
          ll = 0
          for (i in seq_along(p_condX)) {
            idx = p_condX[i]
            val = condX[i]
            if (idx <= p_X1) {
              ll = ll + dbinom(val, 1, piPars[idx, count_lik], log = TRUE)
            } else {
              q = idx - p_X1
              ll = ll + dnorm(val, muPars[q, count_lik], sqrt(tau2Pars[q, count_lik]), log = TRUE)
            }
          }
          log_lik_r[r] = log_w_rk_temp[r] + ll
          count_lik = count_lik + 1L
        }
        
        log_lambda_y_mc[k] = log(n_k[k]) + logSumExp(log_lik_r)
        log_prob_Sx_list[[k]] = log_lik_r
      }
      
      max_ll_y = max(log_lambda_y_mc)
      lambda_y_mc = exp(log_lambda_y_mc - max_ll_y)
    } else {
      lambda_y_mc = n_k
      log_prob_Sx_list = NULL
    }
    
    # ---------------------------------------------------------------------------
    # Draw S_y and S_x, then sample X_mc
    # ---------------------------------------------------------------------------
    S_y_mc = vapply(seq_len(M), function(l) rmultinom_cpp(lambda_y_mc), numeric(1))
    Ky_mc = max(S_y_mc)
    
    S_x_mc = numeric(M)
    n_k_mc = numeric(Ky)
    n_rk_mc = matrix(0, nrow = nrow(n_rk), ncol = ncol(n_rk))
    max_Kx_Sy_mc = numeric(Ky)
    
    X_mc = matrix(nrow = M, ncol = (p_X1 + p_X2))
    
    if (!is.null(condX) && !is.null(p_condX)) {
      for (i in seq_along(p_condX)) {
        X_mc[, p_condX[i]] = condX[i]
      }
    }
    
    count = 1L
    for (k in seq_len(Ky)) {
      ind_S_y_mc = which(S_y_mc == k)
      n_k_temp = length(ind_S_y_mc)
      n_k_mc[k] = n_k_temp
      max_Kx_Sy_temp = max_Kx_Sy[k]
      
      if (n_k_temp > 0) {
        if (max_Kx_Sy_temp == 1) {
          S_x_sub = rep(1, n_k_temp)
        } else {
          if (!is.null(condX) && !is.null(p_condX)) {
            log_p = log_prob_Sx_list[[k]]
            p_Sx = exp(log_p - max(log_p))
            S_x_sub = sample(seq_len(max_Kx_Sy_temp), n_k_temp, replace = TRUE, prob = p_Sx)
          } else {
            S_x_sub = sample(seq_len(max_Kx_Sy_temp), n_k_temp, replace = TRUE,
                             prob = n_rk[k, seq_len(max_Kx_Sy_temp)])
          }
        }
        
        S_x_mc[ind_S_y_mc] = S_x_sub
        max_Kx_Sy_mc[k] = max(S_x_sub)
        
        for (r in seq_len(max_Kx_Sy_temp)) {
          ind_Sy_x_mc = ind_S_y_mc[which(S_x_sub == r)]
          num_Sy_x_mc = length(ind_Sy_x_mc)
          
          n_rk_mc[k, r] = num_Sy_x_mc
          
          if (num_Sy_x_mc > 0) {
            if (p_X1 > 0) {
              for (q in seq_len(p_X1)) {
                if (is.null(p_condX) || !(q %in% p_condX)) {
                  X_mc[ind_Sy_x_mc, q] = rbinom(num_Sy_x_mc, 1, piPars[q, count])
                }
              }
            }
            if (p_X2 > 0) {
              for (q in seq_len(p_X2)) {
                col_idx = p_X1 + q
                if (is.null(p_condX) || !(col_idx %in% p_condX)) {
                  X_mc[ind_Sy_x_mc, col_idx] = rnorm(
                    num_Sy_x_mc, muPars[q, count], sqrt(tau2Pars[q, count])
                  )
                }
              }
            }
          }
          count = count + 1L
        }
      } else {
        count = count + max_Kx_Sy_temp
      }
    }
    
    # ---------------------------------------------------------------------------
    # Conditional weights lambda_Y_z0x_mc and lambda_Y_z1x_mc
    # ---------------------------------------------------------------------------
    lambda_Y_z0x_mc = matrix(nrow = M, ncol = Ky)
    lambda_Y_z1x_mc = matrix(nrow = M, ncol = Ky)
    
    count = 1L
    for (k in seq_len(Ky)) {
      n_k_temp = n_k[k]
      max_Kx_Sy_temp = max_Kx_Sy[k]
      
      log_w_k_temp = log(n_k_temp / N)
      
      w_rk_temp = n_rk[k, seq_len(max_Kx_Sy_temp)]
      sum_w = sum(w_rk_temp)
      if (sum_w == 0) sum_w = 1
      log_w_rk_temp = log(w_rk_temp / sum_w)
      
      log_comps_z0 = matrix(-Inf, nrow = M, ncol = max_Kx_Sy_temp)
      log_comps_z1 = matrix(-Inf, nrow = M, ncol = max_Kx_Sy_temp)
      
      for (r in seq_len(max_Kx_Sy_temp)) {
        prob_x_mc = numeric(M)
        prob_z0_mc = dbinom(z0, 1, etaPars[, count], log = TRUE)
        prob_z1_mc = dbinom(z1, 1, etaPars[, count], log = TRUE)
        
        if (p_X1 > 0) {
          for (q in seq_len(p_X1)) {
            prob_x_mc = prob_x_mc + dbinom(X_mc[, q], 1, piPars[q, count], log = TRUE)
          }
        }
        if (p_X2 > 0) {
          for (q in seq_len(p_X2)) {
            prob_x_mc = prob_x_mc + dnorm(
              X_mc[, p_X1 + q], muPars[q, count], sqrt(tau2Pars[q, count]), log = TRUE
            )
          }
        }
        
        log_comps_z0[, r] = log_w_rk_temp[r] + prob_z0_mc + prob_x_mc
        log_comps_z1[, r] = log_w_rk_temp[r] + prob_z1_mc + prob_x_mc
        count = count + 1L
      }
      
      lambda_Y_z0x_mc[, k] = apply(log_comps_z0, 1, logSumExp) + log_w_k_temp
      lambda_Y_z1x_mc[, k] = apply(log_comps_z1, 1, logSumExp) + log_w_k_temp
    }
    
    lambda_Y_z0x_mc = normalize_log_probs(lambda_Y_z0x_mc)
    lambda_Y_z1x_mc = normalize_log_probs(lambda_Y_z1x_mc)
    
    # ---------------------------------------------------------------------------
    # Means
    # ---------------------------------------------------------------------------
    # Create the design matrices for all M subjects
    matX_z0x_mc = cbind(1, z0, X_mc) 
    matX_z1x_mc = cbind(1, z1, X_mc) 
    MUy_z0x_mc = (matX_z0x_mc %*% betaPars)
    MUy_z1x_mc = (matX_z1x_mc %*% betaPars)
    Ey_z0x_mc = rowSums(MUy_z0x_mc * lambda_Y_z0x_mc)
    Ey_z1x_mc = rowSums(MUy_z1x_mc * lambda_Y_z1x_mc)
    Ey2_z0x_mc = rowSums((MUy_z0x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = Ky)) * lambda_Y_z0x_mc)
    Ey2_z1x_mc = rowSums((MUy_z1x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = Ky)) * lambda_Y_z1x_mc)
    # SDy_z0x_mc = sqrt(pmax(Ey2_z0x_mc - Ey_z0x_mc^{2}, 0))
    # SDy_z1x_mc = sqrt(pmax(Ey2_z1x_mc - Ey_z1x_mc^{2}, 0))
    
    # -----------------------------------------------------------------------------
    lambda_Y_z0x = apply(lambda_Y_z0x_mc, 2, mean)
    lambda_Y_z1x = apply(lambda_Y_z1x_mc, 2, mean)
    MUy_z0x = apply(MUy_z0x_mc, 2, mean)
    MUy_z1x = apply(MUy_z1x_mc, 2, mean)
    Ey_z0x = mean(Ey_z0x_mc)
    Ey_z1x = mean(Ey_z1x_mc)
    Ey2_z0x = mean(Ey2_z0x_mc)
    Ey2_z1x = mean(Ey2_z1x_mc)
    SDy_z0x = sqrt(pmax(Ey2_z0x - Ey_z0x^{2}, 0))
    SDy_z1x = sqrt(pmax(Ey2_z1x - Ey_z1x^{2}, 0))
    
    # ---------------------------------------------------------------------------
    # Main sensitivity / nu / rho loops
    # ---------------------------------------------------------------------------
    for (ii_eta in seq_len(n_eta)) {
      curr_eta_z0 = eta_z0_vec[ii_eta]
      curr_eta_z1 = eta_z1_vec[ii_eta]
      
      for (ii_nu in seq_len(n_nu)) {
        curr_nu = nu_vec[ii_nu]
        curr_log_nu = log(curr_nu)
        safe_curr_log_nu = ifelse(curr_nu == 0, -100, curr_log_nu)
        
        run_z0 = (curr_nu - curr_eta_z0 >= 0)
        run_z1 = (curr_nu - curr_eta_z1 >= 0)
        
        if (run_z0) {
          S_nu_z0x = S_cond_cpp(
            safe_curr_log_nu - curr_eta_z0,
            lambda_Y_z0x_mc, MUy_z0x_mc, as.vector(sqrt_sig2), TRUE
          )
          y_init_z0 = safe_curr_log_nu - curr_eta_z0 + 0.1
        }
        
        if (run_z1) {
          S_nu_z1x = S_cond_cpp(
            safe_curr_log_nu - curr_eta_z1,
            lambda_Y_z1x_mc, MUy_z1x_mc, as.vector(sqrt_sig2), TRUE
          )
          y_init_z1 = safe_curr_log_nu - curr_eta_z1 + 0.1
        }
        
        for (ii_rho in seq_len(n_rho)) {
          curr_rho = rho_vec[ii_rho]
          
          if (run_z0) {
            target_rho_nu_z0x = S_nu_z0x + log(1 - curr_rho)
            opt_z0 = S_optim_cond_cpp(
              target = target_rho_nu_z0x,
              y_init = y_init_z0,
              y_min = safe_curr_log_nu - curr_eta_z0,
              y_max = 100,
              lambda_mat = lambda_Y_z0x_mc,
              mu_mat = MUy_z0x_mc,
              sig = as.vector(sqrt_sig2),
              logt = TRUE
            )
            y_init_z0 = opt_z0$optimizer
            Y_rho_nu_z0x_local[ii_eta, ii_nu, ii_rho] =
              exp(opt_z0$optimizer + curr_eta_z0) - curr_nu
          }
          
          if (run_z1) {
            target_rho_nu_z1x = S_nu_z1x + log(1 - curr_rho)
            opt_z1 = S_optim_cond_cpp(
              target = target_rho_nu_z1x,
              y_init = y_init_z1,
              y_min = safe_curr_log_nu - curr_eta_z1,
              y_max = 100,
              lambda_mat = lambda_Y_z1x_mc,
              mu_mat = MUy_z1x_mc,
              sig = as.vector(sqrt_sig2),
              logt = TRUE
            )
            y_init_z1 = opt_z1$optimizer
            Y_rho_nu_z1x_local[ii_eta, ii_nu, ii_rho] =
              exp(opt_z1$optimizer + curr_eta_z1) - curr_nu
          }
        }
      }
    }
    
    list(z0 = Y_rho_nu_z0x_local, z1 = Y_rho_nu_z1x_local)
  }
  
  # -----------------------------------------------------------------------------
  # Parallel execution
  # -----------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  start_time <- Sys.time()
  print(start_time)
  res_list <- parallel::mclapply(
    X = seq_len(n_MCMC),
    FUN = worker_fun,
    mc.cores = ncores,
    mc.preschedule = FALSE
  )
  end_time <- Sys.time()
  print(end_time)
  print(end_time - start_time)
  
  # -----------------------------------------------------------------------------
  # Reassemble 4D arrays
  # -----------------------------------------------------------------------------
  Y_rho_nu_z0x_store = array(NA_real_, dim = c(n_MCMC, n_eta, n_nu, n_rho))
  Y_rho_nu_z1x_store = array(NA_real_, dim = c(n_MCMC, n_eta, n_nu, n_rho))
  
  for (post_reps in seq_len(n_MCMC)) {
    Y_rho_nu_z0x_store[post_reps, , , ] = res_list[[post_reps]]$z0
    Y_rho_nu_z1x_store[post_reps, , , ] = res_list[[post_reps]]$z1
  }
  
  E_rho_nu_diff_store = Y_rho_nu_z1x_store - Y_rho_nu_z0x_store
  
  # -----------------------------------------------------------------------------
  # Summaries
  # -----------------------------------------------------------------------------
  z0_tables_list = list()
  z1_tables_list = list()
  diff_tables_list = list()
  counter = 1L
  
  for (ii_eta in seq_len(n_eta)) {
    val_eta_z0 = eta_z0_vec[ii_eta]
    val_eta_z1 = eta_z1_vec[ii_eta]
    
    for (ii_nu in seq_len(n_nu)) {
      val_nu = nu_vec[ii_nu]
      
      for (ii_rho in seq_len(n_rho)) {
        val_rho = rho_vec[ii_rho]
        
        draws_z0 = Y_rho_nu_z0x_store[, ii_eta, ii_nu, ii_rho]
        draws_z1 = Y_rho_nu_z1x_store[, ii_eta, ii_nu, ii_rho]
        draws_diff = E_rho_nu_diff_store[, ii_eta, ii_nu, ii_rho]
        
        if (!all(is.na(draws_z0))) {
          summ_z0 = POSTsummary(draws_z0, esttype, quantile_alpha)
          summ_z0 = cbind(eta_z0 = val_eta_z0, nu = val_nu, rho = val_rho, summ_z0)
          z0_tables_list[[counter]] = summ_z0
        }
        
        if (!all(is.na(draws_z1))) {
          summ_z1 = POSTsummary(draws_z1, esttype, quantile_alpha)
          summ_z1 = cbind(eta_z1 = val_eta_z1, nu = val_nu, rho = val_rho, summ_z1)
          z1_tables_list[[counter]] = summ_z1
        }
        
        if (!all(is.na(draws_diff))) {
          summ_diff = POSTsummary(draws_diff, esttype, quantile_alpha)
          summ_diff = cbind(
            eta_z0 = val_eta_z0, eta_z1 = val_eta_z1,
            nu = val_nu, rho = val_rho, summ_diff
          )
          diff_tables_list[[counter]] = summ_diff
        }
        
        counter = counter + 1L
      }
    }
  }
  
  Y_rho_nu_z0x_result  = if (length(z0_tables_list)   > 0) do.call(rbind, z0_tables_list)   else NULL
  Y_rho_nu_z1x_result  = if (length(z1_tables_list)   > 0) do.call(rbind, z1_tables_list)   else NULL
  E_rho_nu_diff_result = if (length(diff_tables_list) > 0) do.call(rbind, diff_tables_list) else NULL
  
  POSTresult = list(
    Y_rho_nu_z0x_result  = Y_rho_nu_z0x_result,
    Y_rho_nu_z1x_result  = Y_rho_nu_z1x_result,
    E_rho_nu_diff_result = E_rho_nu_diff_result,
    Y_rho_nu_z0x_store   = Y_rho_nu_z0x_store,
    Y_rho_nu_z1x_store   = Y_rho_nu_z1x_store,
    E_rho_nu_diff_store  = E_rho_nu_diff_store
  )
  
  return(POSTresult)
}

SURVIVAL_EDPMM_POST_ftn = function(object, M, t_grid, ftn = c("S", "f"), 
                                 condX = NULL, p_condX = NULL, esttype = "mean"){
  
  # Ensure ftn contains valid options
  ftn = match.arg(ftn, several.ok = TRUE)
  calc_S = "S" %in% ftn
  calc_f = "f" %in% ftn
  
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
    
    zero_rows = which(row_sums == 0)
    if(length(zero_rows) > 0) {
      probs_shifted[zero_rows,] = 1 / ncol(log_prob_matrix) 
      row_sums[zero_rows] = 1
    }
    return(probs_shifted / row_sums)
  }
  
  # -----------------------------------------------------------------------------
  n_t = length(t_grid)
  level = 0.05
  quantile_alpha = c(level/2, 1-level/2)
  z0 = 0
  z1 = 1
  
  N = object$constants$N
  # num_MC = object$constants$num_MC
  n_MCMC = object$constants$n_MCMC
  p_matX = object$constants$p_matX
  p_X1 = object$constants$p_X1
  p_X2 = object$constants$p_X2
  
  a_beta = object$priors$a_beta; B_beta = object$priors$B_beta
  a_sig2 = object$priors$a_sig2; b_sig2 = object$priors$b_sig2
  a_eta = object$priors$a_eta;   b_eta = object$priors$b_eta
  a_pi = object$priors$a_pi;     b_pi = object$priors$b_pi
  a_mu = object$priors$a_mu;     b_mu = object$priors$b_mu
  a_tau2 = object$priors$a_tau2; b_tau2 = object$priors$b_tau2
  
  iter_check = floor(n_MCMC/10)
  
  # Conditionally allocate storage arrays to save memory
  if (calc_S) {
    S_z0_store = matrix(NA, nrow = n_MCMC, ncol = n_t)
    S_z1_store = matrix(NA, nrow = n_MCMC, ncol = n_t)
  }
  if (calc_f) {
    f_z0_store = matrix(NA, nrow = n_MCMC, ncol = n_t)
    f_z1_store = matrix(NA, nrow = n_MCMC, ncol = n_t)
  }
  
  for (post_reps in 1:n_MCMC) {
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
    if (!is.null(condX) && !is.null(p_condX)) {
      log_lambda_y_mc = numeric(Ky)
      log_prob_Sx_list = list()
      
      count_lik = 1
      for (k in 1:Ky) {
        max_Kx_Sy_temp = max_Kx_Sy[k]
        w_rk_temp = n_rk[k, 1:max_Kx_Sy_temp]
        
        sum_w = sum(w_rk_temp)
        if (sum_w == 0) sum_w = 1 
        log_w_rk_temp = log(w_rk_temp / sum_w)
        
        log_lik_r = numeric(max_Kx_Sy_temp)
        for (r in 1:max_Kx_Sy_temp) {
          ll = 0
          for (i in seq_along(p_condX)) {
            idx = p_condX[i]
            val = condX[i]
            if (idx <= p_X1) {
              ll = ll + dbinom(val, 1, piPars[idx, count_lik], log = TRUE)
            } else {
              q = idx - p_X1
              ll = ll + dnorm(val, muPars[q, count_lik], sqrt(tau2Pars[q, count_lik]), log = TRUE)
            }
          }
          log_lik_r[r] = log_w_rk_temp[r] + ll
          count_lik = count_lik + 1
        }
        log_lambda_y_mc[k] = log(n_k[k]) + logSumExp(log_lik_r)
        log_prob_Sx_list[[k]] = log_lik_r
      }
      max_ll_y = max(log_lambda_y_mc)
      lambda_y_mc = exp(log_lambda_y_mc - max_ll_y)
    } else {
      lambda_y_mc = n_k 
      log_prob_Sx_list = NULL
    }
    
    # -----------------------------------------------------------------------------
    S_y_mc = sapply(1:M, function(l) rmultinom_cpp(lambda_y_mc))
    unique_S_y_mc = unique(S_y_mc)
    Ky_mc = max(unique_S_y_mc)
    
    # -----------------------------------------------------------------------------
    S_x_mc = numeric(M)
    n_k_mc  = numeric(Ky)
    n_rk_mc = matrix(0, nrow = nrow(n_rk), ncol = ncol(n_rk))
    max_Kx_Sy_mc = numeric(Ky)
    X_mc = matrix(nrow = M, ncol = (p_X1 + p_X2))
    
    if (!is.null(condX) && !is.null(p_condX)) {
      for (i in seq_along(p_condX)) {
        X_mc[, p_condX[i]] = condX[i]
      }
    }
    
    count = 1
    for (k in 1:Ky) {
      ind_S_y_mc = which(S_y_mc == k)
      n_k_temp = length(ind_S_y_mc)
      n_k_mc[k] = n_k_temp
      max_Kx_Sy_temp = max_Kx_Sy[k]
      
      if (n_k_temp > 0) {
        if (max_Kx_Sy_temp == 1) {
          S_x_sub = rep(1, n_k_temp)
        } else {
          if (!is.null(condX) && !is.null(p_condX)) {
            log_p = log_prob_Sx_list[[k]]
            p_Sx = exp(log_p - max(log_p))
            S_x_sub = sample(1:max_Kx_Sy_temp, n_k_temp, replace = TRUE, prob = p_Sx)
          } else {
            S_x_sub = sample(1:max_Kx_Sy_temp, n_k_temp, replace = TRUE, prob = n_rk[k,1:max_Kx_Sy_temp])
          }
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
                if (is.null(p_condX) || !(q %in% p_condX)) {
                  X_mc[ind_Sy_x_mc, q] = rbinom(num_Sy_x_mc, 1, piPars[q, count])
                }
              }
            }
            if (p_X2 > 0) {
              for (q in 1:p_X2) {
                col_idx = p_X1 + q
                if (is.null(p_condX) || !(col_idx %in% p_condX)) {
                  X_mc[ind_Sy_x_mc, col_idx] = rnorm(num_Sy_x_mc, muPars[q, count], sqrt(tau2Pars[q, count]))
                }
              }
            }
          }
          count = count + 1
        }
      } else {
        count = count + max_Kx_Sy_temp
      }
    }
    
    # -----------------------------------------------------------------------------
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
    # Create the design matrices for all M subjects
    matX_z0x_mc = cbind(1, z0, X_mc) 
    matX_z1x_mc = cbind(1, z1, X_mc) 
    MUy_z0x_mc = (matX_z0x_mc %*% betaPars)
    MUy_z1x_mc = (matX_z1x_mc %*% betaPars)
    Ey_z0x_mc = rowSums(MUy_z0x_mc * lambda_Y_z0x_mc)
    Ey_z1x_mc = rowSums(MUy_z1x_mc * lambda_Y_z1x_mc)
    Ey2_z0x_mc = rowSums((MUy_z0x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = Ky)) * lambda_Y_z0x_mc)
    Ey2_z1x_mc = rowSums((MUy_z1x_mc^{2} + matrix(sig2Pars, nrow = M, ncol = Ky)) * lambda_Y_z1x_mc)
    # SDy_z0x_mc = sqrt(pmax(Ey2_z0x_mc - Ey_z0x_mc^{2}, 0))
    # SDy_z1x_mc = sqrt(pmax(Ey2_z1x_mc - Ey_z1x_mc^{2}, 0))
    
    # -----------------------------------------------------------------------------
    lambda_Y_z0x = apply(lambda_Y_z0x_mc, 2, mean)
    lambda_Y_z1x = apply(lambda_Y_z1x_mc, 2, mean)
    MUy_z0x = apply(MUy_z0x_mc, 2, mean)
    MUy_z1x = apply(MUy_z1x_mc, 2, mean)
    Ey_z0x = mean(Ey_z0x_mc)
    Ey_z1x = mean(Ey_z1x_mc)
    Ey2_z0x = mean(Ey2_z0x_mc)
    Ey2_z1x = mean(Ey2_z1x_mc)
    SDy_z0x = sqrt(pmax(Ey2_z0x - Ey_z0x^{2}, 0))
    SDy_z1x = sqrt(pmax(Ey2_z1x - Ey_z1x^{2}, 0))
    
    # -----------------------------------------------------------------------------
    # SURVIVAL & DENSITY EVALUATION LOOP
    # -----------------------------------------------------------------------------
    for (ii_t in 1:n_t) {
      curr_t = t_grid[ii_t]
      
      if (curr_t <= 0) {
        if (calc_S) {
          S_z0_store[post_reps, ii_t] = 1
          S_z1_store[post_reps, ii_t] = 1
        }
        if (calc_f) {
          f_z0_store[post_reps, ii_t] = 0
          f_z1_store[post_reps, ii_t] = 0
        }
      } else {
        curr_log_t = log(curr_t)
        
        if (calc_S) {
          log_S_z0 = S_cond_cpp(curr_log_t, lambda_Y_z0x_mc, MUy_z0x_mc, sqrt(sig2Pars), TRUE)
          log_S_z1 = S_cond_cpp(curr_log_t, lambda_Y_z1x_mc, MUy_z1x_mc, sqrt(sig2Pars), TRUE)
          S_z0_store[post_reps, ii_t] = exp(log_S_z0)
          S_z1_store[post_reps, ii_t] = exp(log_S_z1)
        }
        
        if (calc_f) {
          log_f_z0 = f_cond_cpp(curr_log_t, lambda_Y_z0x_mc, MUy_z0x_mc, sqrt(sig2Pars), TRUE)
          log_f_z1 = f_cond_cpp(curr_log_t, lambda_Y_z1x_mc, MUy_z1x_mc, sqrt(sig2Pars), TRUE)
          
          # Jacobian applied: exp(log_f) / t yields density of survival time T
          f_z0_store[post_reps, ii_t] = exp(log_f_z0) / curr_t
          f_z1_store[post_reps, ii_t] = exp(log_f_z1) / curr_t
        }
      }
    }
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing Curve Evaluation",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  
  # -----------------------------------------------------------------------------
  # Combine outputs conditionally
  # -----------------------------------------------------------------------------
  result_list = list()
  
  if (calc_S) {
    z0_S_list = list()
    z1_S_list = list()
    for(ii_t in 1:n_t) {
      val_t = t_grid[ii_t]
      z0_S_list[[ii_t]] = cbind(t = val_t, POSTsummary(S_z0_store[, ii_t], esttype, quantile_alpha))
      z1_S_list[[ii_t]] = cbind(t = val_t, POSTsummary(S_z1_store[, ii_t], esttype, quantile_alpha))
    }
    result_list$S_z0_result = do.call(rbind, z0_S_list)
    result_list$S_z1_result = do.call(rbind, z1_S_list)
    result_list$S_z0_store  = S_z0_store
    result_list$S_z1_store  = S_z1_store
  }
  
  if (calc_f) {
    z0_f_list = list()
    z1_f_list = list()
    for(ii_t in 1:n_t) {
      val_t = t_grid[ii_t]
      z0_f_list[[ii_t]] = cbind(t = val_t, POSTsummary(f_z0_store[, ii_t], esttype, quantile_alpha))
      z1_f_list[[ii_t]] = cbind(t = val_t, POSTsummary(f_z1_store[, ii_t], esttype, quantile_alpha))
    }
    result_list$f_z0_result = do.call(rbind, z0_f_list)
    result_list$f_z1_result = do.call(rbind, z1_f_list)
    result_list$f_z0_store  = f_z0_store
    result_list$f_z1_store  = f_z1_store
  }
  
  return(result_list)
}

# ------------------------------------------------------------------------
SURVIVAL_EDPMM = function(object, M = 2e3, nu = 0, rho = 0.5, 
                          gibbs_iter = 2e4, gibbs_burnin = 2e4, gibbs_thin = 1e2, 
                          esttype = "mean", saveall = FALSE){
  # object = temp_data
  MCMCresult = SURVIVAL_EDPMM_MCMC(object, gibbs_iter, gibbs_burnin, gibbs_thin)
  POSTresult = SURVIVAL_EDPMM_POST(MCMCresult, M, nu, rho, 
                                   esttype = esttype)
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