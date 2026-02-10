#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <Rmath.h>
//#include "utilFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

template <typename T>
inline bool approx_equal_cpp(const T& lhs, const T& rhs, double tol = 1e-10) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

umat unique_rows(const umat& M) {
  int n_temp = M.n_rows;
  uvec ind_temp = arma::zeros<arma::uvec>(n_temp);
  
  for (int i = 0; i < n_temp; i++) {
    for (int ii = i + 1; ii < n_temp; ii++) {
      if (approx_equal_cpp(M.row(i), M.row(ii))) { 
        ind_temp(ii) = 1; 
        break; 
      }
    }
  }
  
  ind_temp = find(ind_temp == 0);
  umat unique_M = M.rows(ind_temp);
  
  uvec unique_M1 = unique_M.col(0);
  uvec unique_M2 = unique_M.col(1);
  
  uvec sorted_unique_M1 = sort(unique_M1);
  uvec sorted_unique_M2 = sort(unique_M2);
  
  uvec M1_temp = unique(sorted_unique_M1);
  int n_temp1 = M1_temp.size();
  int count = 0;
  for (int i = 0; i < n_temp1; i++) {
    ind_temp = find(unique_M1 == M1_temp(i));
    int n_temp2 = ind_temp.size();
    uvec unique_M2_temp = sort(unique_M2(ind_temp));
    
    for (int ii = 0; ii < n_temp2; ii++) {
      sorted_unique_M2(count) = unique_M2_temp(ii);
      count++;
    }
  }
  
  umat sorted_unique_M = join_rows(sorted_unique_M1, sorted_unique_M2);
  
  return sorted_unique_M;
}

int mod(int const& a, int const& b){
  return a - floor(a/b)*b;
}

mat reorder(mat MAT, int const& k_old, int const& k_new) {
  vec MAT_new = MAT.col(k_old);
  MAT.insert_cols(k_new, MAT_new);
  
  if (k_old > k_new) {
    MAT.shed_col(k_old);
  } else {
    MAT.shed_col((k_old+1));
  }
  
  return(MAT);
}

// [[Rcpp::export]]
int rank_cpp(mat const& MAT) {
  return rank(MAT);
}

// // [[Rcpp::export]]
// mat inv_cpp(mat const& MAT) {
//   int p_temp = MAT.n_cols;
//   int rank_temp = rank(MAT);
//   
//   mat MATinv;
//   if (p_temp > rank_temp){
//     MATinv = pinv(MAT);
//   } else if (MAT.is_sympd()){
//     MATinv = inv_sympd(symmatu(MAT));
//   } else {
//     MATinv = inv(MAT);
//   }
//   
//   return MATinv;
// }

// [[Rcpp::export]]
mat inv_cpp(mat const& MAT) {
  mat MATinv;
  bool success = false;
  
  // Attempt 1: Symmetric Positive Definite Inversion (Fastest)
  // symmatu() forces symmetry using the upper triangle, fixing numerical noise
  success = inv_sympd(MATinv, symmatu(MAT));
  
  if (!success) {
    // Attempt 2: General Inversion (Handles full-rank but non-PD/asymmetric)
    // 'inv' returns false if singular, rather than crashing
    success = inv(MATinv, MAT);
  }
  
  if (!success) {
    // Attempt 3: Pseudo-Inverse (Safest, handles Rank Deficiency/Singularity)
    // This replaces your manual 'rank' check and covers all edge cases
    MATinv = pinv(MAT);
  }
  
  return MATinv;
}

// [[Rcpp::export]]
vec count_cpp(vec const& VEC, int const& p) {
  vec count(p);
  for (int q = 0; q < p; q++) {
    uvec index = find(VEC == q);
    count(q) = index.n_elem; 
  }
  
  return(count);
}

vec IQRoutlier_cpp(vec const& x_temp) { 
  vec Q_temp = {0.25,0.5,0.75};
  vec Qx_temp = quantile(x_temp, Q_temp);
  
  double med_temp = Qx_temp[1];
  double IQR_temp = Qx_temp[2] - Qx_temp[0];
  vec IQRoutlier = {med_temp - 1.5*IQR_temp, med_temp + 1.5*IQR_temp};
  return(IQRoutlier);
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
vec rtruncnorm_cpp(int n, double mu, double sigma, 
                   double lower, double upper) {
  if (lower >= upper) {
    Rcpp::stop("'lower' bound must be less than 'upper' bound.");
  }
  if (sigma <= 0) {
    Rcpp::stop("'sigma' must be positive.");
  }
  
  // Standardize the bounds
  double alpha = (lower - mu) / sigma;
  double beta = (upper - mu) / sigma;
  
  vec z(n);
  
  // --- CASE 1: FAR RIGHT TAIL (Stable) ---
  if (alpha > 5) {
    double prob_tail_lower = R::pnorm(alpha, 0.0, 1.0, false, false); 
    double prob_tail_upper = R::pnorm(beta, 0.0, 1.0, false, false);
    
    vec u = randu<vec>(n);
    vec target_prob_tail = prob_tail_upper + u * (prob_tail_lower - prob_tail_upper);
    
    for(int i = 0; i < n; ++i) {
      z[i] = R::qnorm(target_prob_tail[i], 0.0, 1.0, false, false);
    }
    
    // --- CASE 2: FAR LEFT TAIL (New Fix) ---
    // If the interval is far to the left (e.g. [-100, -90]), pnorm underflows to 0.
    // We exploit symmetry: Sample Z' from [ -beta, -alpha ] and return -Z'.
  } else if (beta < -5) {
    double alpha_sym = -beta;
    double beta_sym  = -alpha;
    
    // Use the logic from Case 1 on the positive mirror image
    double prob_tail_lower = R::pnorm(alpha_sym, 0.0, 1.0, false, false); 
    double prob_tail_upper = R::pnorm(beta_sym, 0.0, 1.0, false, false);
    
    vec u = randu<vec>(n);
    vec target_prob_tail = prob_tail_upper + u * (prob_tail_lower - prob_tail_upper);
    
    for(int i = 0; i < n; ++i) {
      // Get positive Z, then flip sign
      double z_sym = R::qnorm(target_prob_tail[i], 0.0, 1.0, false, false);
      z[i] = -z_sym; 
    }
    
    // --- CASE 3: BODY (Standard) ---
  } else {
    double prob_lower = R::pnorm(alpha, 0.0, 1.0, true, false);
    double prob_upper = R::pnorm(beta, 0.0, 1.0, true, false);
    
    vec u = randu<vec>(n);
    vec target_prob = prob_lower + u * (prob_upper - prob_lower);
    
    for(int i = 0; i < n; ++i) {
      z[i] = R::qnorm(target_prob[i], 0.0, 1.0, true, false);
    }
  }
  
  return mu + sigma * z;
}

// [[Rcpp::export]]
vec rtrunct_cpp(int n, double mu, double scale, double df, double lower, double upper) {
  if (std::isnan(mu) || std::isnan(scale) || std::isnan(lower)) {
    Rcpp::warning("rtrunct_cpp: NaN detected in inputs. Returning NaN.");
    vec out(n); out.fill(datum::nan); return out;
  }
  if (scale <= 0) {
    Rcpp::warning("rtrunct_cpp: Scale <= 0. Returning lower bound.");
    vec out(n); out.fill(lower); return out;
  }
  if (lower >= upper) {
    vec out(n); out.fill(lower); return out;
  }
  
  double alpha = (lower - mu) / scale;
  double beta  = (upper - mu) / scale;
  
  vec z(n);
  
  // --- CASE 1: FAR RIGHT TAIL ---
  if (alpha > 6) {
    double lp_lower = R::pt(alpha, df, false, true); // log(P(T > alpha))
    double lp_upper = (upper == datum::inf) ? -datum::inf : R::pt(beta, df, false, true);
    
    vec u = randu<vec>(n);
    for(int i = 0; i < n; ++i) {
      double target_lp;
      if (upper == datum::inf) {
        target_lp = std::log(u[i]) + lp_lower;
      } else {
        // Log-space interpolation
        if (std::abs(lp_lower - lp_upper) < 1e-12) {
          target_lp = lp_lower; 
        } else {
          double p_diff = std::exp(lp_lower) - std::exp(lp_upper);
          if (p_diff <= 0) target_lp = lp_lower;
          else target_lp = std::log(std::exp(lp_upper) + u[i] * p_diff);
        }
      }
      z[i] = R::qt(target_lp, df, false, true); // log.p=TRUE
      if (z[i] == datum::inf) z[i] = alpha + 1e-3; 
    }
    
    // --- CASE 2: FAR LEFT TAIL (New Fix) ---
  } else if (beta < -6) {
    // Symmetry: T(df) is symmetric. 
    // Sample T' from [-beta, -alpha], return -T'
    double alpha_sym = -beta;
    double beta_sym  = -alpha; // Note: -alpha is larger than -beta
    
    // Same logic as Case 1 but using sym bounds
    double lp_lower = R::pt(alpha_sym, df, false, true); 
    double lp_upper = (lower == -datum::inf) ? -datum::inf : R::pt(beta_sym, df, false, true);
    
    vec u = randu<vec>(n);
    for(int i = 0; i < n; ++i) {
      double target_lp;
      // Note: "upper" for the mirror image corresponds to "lower" of original
      if (lower == -datum::inf) { 
        target_lp = std::log(u[i]) + lp_lower;
      } else {
        if (std::abs(lp_lower - lp_upper) < 1e-12) {
          target_lp = lp_lower; 
        } else {
          double p_diff = std::exp(lp_lower) - std::exp(lp_upper);
          if (p_diff <= 0) target_lp = lp_lower;
          else target_lp = std::log(std::exp(lp_upper) + u[i] * p_diff);
        }
      }
      double z_sym = R::qt(target_lp, df, false, true);
      if (z_sym == datum::inf) z_sym = alpha_sym + 1e-3;
      
      z[i] = -z_sym; // FLIP SIGN
    }
    
    // --- CASE 3: BODY ---
  } else {
    double p_cdf_lower = (lower == -datum::inf) ? 0.0 : R::pt(alpha, df, true, false);
    double p_cdf_upper = (upper == datum::inf) ? 1.0 : R::pt(beta, df, true, false);
    
    vec u = randu<vec>(n);
    vec target_cdf = p_cdf_lower + u * (p_cdf_upper - p_cdf_lower);
    
    for(int i = 0; i < n; ++i) {
      z[i] = R::qt(target_cdf[i], df, true, false);
    }
  }
  
  return mu + scale * z;
}

// [[Rcpp::export]]
vec rbeta_cpp(int const& N, double const& s, double const& r) {
  vec rgamma1 = randg(N, distr_param(s,1.0));
  vec rgamma2 = randg(N, distr_param(r,1.0));
  vec rbeta = rgamma1/(rgamma1+rgamma2);
  
  return(rbeta);
}

// [[Rcpp::export]]
vec invlogit_cpp(vec const& x_temp) { 
  return((1/(1+exp(-x_temp))));
}

// [[Rcpp::export]]
vec factorial_cpp(vec const& N, bool const& logt) {
  int p = N.size();
  
  vec factorial = zeros(p);
  for (int q = 0; q < p; q++) {
    int N_temp = N(q);
    if (N_temp > 0){
      for (int i = 1; i <= N_temp; i++) {
        factorial(q) += log(i);
      }
    }
  }
  
  if (logt == false) {factorial = exp(factorial);}
  
  return(factorial);
}

// [[Rcpp::export]]
int rmultinom_cpp(vec const& prob) {
  int p = prob.size();
  
  vec cumprob = cumsum(prob/sum(prob));
  
  int rmultinom = 1;
  double rnd = randu(1)[0];
  for (int q = 0; q < p; q++) {
    if (rnd>cumprob(q)) {
      rmultinom = rmultinom+1;
    }
  }
  
  return(rmultinom);
}

// [[Rcpp::export]]
int dmultinom_cpp(vec const& X, vec const& prob, bool const& logt) {
  vec N(1);
  N(0) = X.size();
  vec log_prob = log(prob);
  log_prob.replace(-datum::inf,0);
  
  double dmultinom = factorial_cpp(N, true)[0] - sum(factorial_cpp(X, true) + X % log_prob);
  if (logt == false) {dmultinom = exp(dmultinom);}
  
  return(dmultinom);
}

// [[Rcpp::export]]
mat rdirichlet_cpp(int const& N, vec const& alpha) {
  int p = alpha.size();
  
  mat rdirichlet(N,p);
  for (int q = 0; q < p; q++) {
    rdirichlet.col(q) = randg(N,distr_param(alpha(q),1.0));
  }
  vec sum_rdirichlet = sum(rdirichlet,1);
  rdirichlet.each_col() /= sum_rdirichlet;
  
  return rdirichlet.t();
}

// [[Rcpp::export]]
double ddirichlet_cpp(vec const& X, vec const& alpha, bool const& logt) {
  int p = alpha.size();
  vec ones_p = ones(p);
  
  double ddirichlet = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - ones_p) % log(X));
  if (logt == false) {ddirichlet = exp(ddirichlet);}
  
  return ddirichlet;
}

// [[Rcpp::export]]
vec rscainvchisq_cpp(int const& N, double const& a_tau2, double const& b_tau2) {
  double log_numerator = log(a_tau2 * b_tau2);
  
  // Upper Limit: exp(+30) = 1.068647e+13
  // Lower Limit: exp(-30) = 9.357623e-14
  double log_upper_bound = 30; 
  double log_lower_bound = -30;
  
  // Identity: ChiSq(df) is equivalent to Gamma(shape = df/2, scale = 2)
  double shape_param = a_tau2 / 2.0;
  double scale_param = 2.0;
  vec chi_vals = randg(N, distr_param(shape_param, scale_param));
  
  vec rtau(N);
  for(int i = 0; i < N; ++i) {
    double val = chi_vals(i);
    if (val < 1e-100) {
      val = 1e-100;
    }
    double log_result = log_numerator - log(val);
    if (log_result > log_upper_bound) {
      log_result = log_upper_bound;
    } else if (log_result < log_lower_bound) {
      log_result = log_lower_bound;
    }
    rtau(i) = exp(log_result);
  }
  
  return rtau;
}

// // [[Rcpp::export]]
// vec rscainvchisq_cpp(int const& N, double const& a_tau2, double const& b_tau2) {
//   vec ones_N = ones(N);
//   vec ab_tau_N = a_tau2 * b_tau2 * ones_N;
//   
//   vec rtau = ab_tau_N/chi2rnd(a_tau2,N);
//   return rtau;
// }

// [[Rcpp::export]]
double dscainvchisq_cpp(double const& tau, double const& a_tau2, double const& b_tau2, bool const& logt) {
  double ab_tau = a_tau2 * b_tau2;
  double a_tau_half = a_tau2/2;
  double ab_tau_half = ab_tau/2;
  
  double dtau = a_tau_half * log(ab_tau_half) - ab_tau_half/tau - 
    lgamma(a_tau_half) - (1+a_tau_half) * log(tau);
  if (logt == false) {dtau = exp(dtau);}
  
  return dtau;
}

// randg: shape and scale at default
// rgamma: shape and rate at default 
// => USE RATE
// [[Rcpp::export]]
vec rinvgamma_cpp(int const& N, double const& a_tau2, double const& b_tau2) {
  vec ones_N = ones(N);
  vec rtau = randg(N, distr_param(a_tau2, 1/b_tau2));
  return rtau;
}

// randg: shape and scale at default
// rgamma: shape and rate at default
// => USE RATE
// [[Rcpp::export]]
double dinvgamma_cpp(double const& tau, double const& a_tau2, double const& b_tau2, bool const& logt) {
  double dtau = - a_tau2 * log(b_tau2) - lgamma(a_tau2) - (a_tau2 + 1) * log(tau) + b_tau2/tau;
  if (logt == false) {dtau = exp(dtau);}
  return dtau;
}

// [[Rcpp::export]]
mat rmvn_cpp(int const& n, vec const& MU, mat const& SIG) {
  int k = MU.n_elem;
  
  // 1. Generate Standard Normals: Z ~ N(0, 1)
  // Direct generation is faster than creating a random object
  mat Z = randn(k, n);
  
  // 2. Cholesky Decomposition
  mat U;
  // 'chol' returns false if unstable (FAST), unlike 'mvnrnd' which throws an error (SLOW)
  bool success = chol(U, SIG); 
  
  // 3. Robust Fallback (Only executes if matrix is singular)
  if (!success) {
    mat SIG_copy = SIG;
    SIG_copy.diag() += 1e-8; // Add Jitter
    
    if (!chol(U, SIG_copy)) {
      // Ultimate Fallback: Eigen Decomposition
      // Guaranteed to work even if Cholesky fails twice
      vec eigval; 
      mat eigvec;
      eig_sym(eigval, eigvec, SIG);
      
      // Force non-negative eigenvalues (remove numerical noise)
      eigval.elem(find(eigval < 0)).zeros();
      
      // Transform: X = MU + (V * sqrt(D)) * Z
      mat X = (eigvec * diagmat(sqrt(eigval))) * Z;
      X.each_col() += MU; 
      return X;
    }
  }
  
  // 4. Transform: X = MU + U' * Z
  // Using .each_col() is faster than repmat() as it avoids allocating a new mean matrix
  mat X = trans(U) * Z;
  X.each_col() += MU;
  
  return X;
}

// [[Rcpp::export]]
double dmvn_cpp(vec X_temp, vec MU, mat SIG, bool logt) {
  int p = SIG.n_cols;
  mat Chol_Sig = chol(SIG);
  vec temp1 = solve(Chol_Sig, X_temp - MU); //check
  double temp2 = sum(temp1 % temp1);
  
  double dmvn = - sum(log(Chol_Sig.diag())) - 0.5 * p * log(2 * M_PI) - 0.5 * temp2;
  if(logt==false) {dmvn = exp(dmvn);}
  
  return dmvn;
}

// [[Rcpp::export]]
mat rinvwish_cpp(double const& a_PSI, mat const& PSI) {
  return iwishrnd(PSI, a_PSI);
}

// =============================================================================
// GIBBS STEP HELPER FUNCTIONS FOR CLUSTER PARAMETERS
// =============================================================================

// For Outcome Model: Y ~ Normal(matXy * beta, sig2)
// [[Rcpp::export]]
List update_betasig_cpp(vec const& y_k, mat const& X_k,
                        vec const& beta_current,
                        vec const& a_beta, mat const& B_beta, mat const& Binv_beta,
                        double const& a_sig2, double const& b_sig2) {
  int n_k = y_k.n_elem;
  int p = beta_current.n_elem;
  
  vec residuals = y_k - X_k * beta_current;
  vec beta_diff = beta_current - a_beta;
  
  double sse_data = dot(residuals, residuals);
  double sse_prior = as_scalar(beta_diff.t() * Binv_beta * beta_diff); // Prior contribution
  
  double a_sig2_new = a_sig2 + n_k + p;
  double b_sig2_new = (a_sig2 * b_sig2 + sse_data + sse_prior) / a_sig2_new;
  double sig2_new = rscainvchisq_cpp(1, a_sig2_new, b_sig2_new)[0];
  
  mat Bn_beta_inv = (Binv_beta + X_k.t() * X_k) / sig2_new;
  mat Bn_beta = inv_cpp(Bn_beta_inv);
  vec an_beta = Bn_beta * ((Binv_beta * a_beta + X_k.t() * y_k) / sig2_new);
  
  vec beta_new = rmvn_cpp(1, an_beta, Bn_beta).col(0);
  
  return List::create(_["beta_par"] = beta_new,
                      _["sig2_par"] = sig2_new);
}
// // [[Rcpp::export]]
// List update_betasig_cpp(vec const& y_k, mat const& X_k,
//                         vec const& beta_current,
//                         vec const& a_beta, mat const& B_beta,
//                         double const& a_sig2, double const& b_sig2) {
//   int n_k = y_k.n_elem;
//   mat B_beta_inv = inv_cpp(B_beta);
//   
//   // 1. Update sig2, conditional on the CURRENT beta
//   vec residuals = y_k - X_k * beta_current;
//   double sse = dot(residuals, residuals);
//   
//   double a_sig2_new = a_sig2 + n_k;
//   double b_sig2_new = (a_sig2 * b_sig2 + sse) / a_sig2_new;
//   
//   double sig2_new = rscainvchisq_cpp(1, a_sig2_new, b_sig2_new)[0];
//   
//   // 2. Update beta, conditional on the NEWLY SAMPLED sig2
//   mat Bn_beta_inv = B_beta_inv + (X_k.t() * X_k) / sig2_new;
//   mat Bn_beta = inv_cpp(Bn_beta_inv);
//   vec an_beta = Bn_beta * (B_beta_inv * a_beta + (X_k.t() * y_k) / sig2_new);
//   
//   vec beta_new = rmvn_cpp(1, an_beta, Bn_beta);
//   
//   return List::create(_["beta_par"] = beta_new,
//                       _["sig2_par"] = sig2_new);
// }

// // For Propensity Score Model: Z ~ Logit(matXz * eta)
// // Uses Metropolis-Hastings and your requested variable names
// // [[Rcpp::export]]
// vec update_eta_cpp(vec const& Y_temp, mat const& matX_temp,
//                    vec const& beta_curr,
//                    vec const& a_beta, mat const& B_beta) {
//   
//   // 1. Propose a new beta
//   // Note: B_beta here acts as the proposal covariance matrix
//   vec beta_prop = rmvn_cpp(1, beta_curr, B_beta);
//   
//   // 2. Calculate log-posterior for PROPOSED beta
//   vec eta_prop = matX_temp * beta_prop;
//   vec p_prop = 1.0 / (1.0 + exp(-eta_prop));
//   // Ensure probabilities are not exactly 0 or 1 for stability
//   p_prop.clamp(1e-15, 1.0 - 1e-15); 
//   double log_lik_data_prop = accu(Y_temp % log(p_prop) + (1.0 - Y_temp) % log(1.0 - p_prop));
//   double log_prior_prop = dmvn_cpp(beta_prop, a_beta, B_beta, true);
//   double log_posterior_prop = log_lik_data_prop + log_prior_prop;
//   
//   // 3. Calculate log-posterior for CURRENT beta
//   vec eta_temp = matX_temp * beta_curr;
//   vec p_temp = 1.0 / (1.0 + exp(-eta_temp));
//   p_temp.clamp(1e-15, 1.0 - 1e-15);
//   double log_lik_data_temp = accu(Y_temp % log(p_temp) + (1.0 - Y_temp) % log(1.0 - p_temp));
//   double log_prior_temp = dmvn_cpp(beta_curr, a_beta, B_beta, true);
//   double log_posterior_temp = log_lik_data_temp + log_prior_temp;
//   
//   // 4. Metropolis-Hastings acceptance step
//   double acceptance_ratio = exp(log_posterior_prop - log_posterior_temp);
//   
//   if (R::runif(0, 1) < acceptance_ratio) {
//     return beta_prop; // Accept proposal
//   } else {
//     return beta_curr; // Reject proposal
//   }
// }
// 
// // [[Rcpp::export]]
// vec update_eta_latent_cpp(vec const& Zstar, mat const& X,
//                           vec const& a_beta, mat const& Binv_beta) {
//   mat Binv_beta_new = Binv_beta + X.t() * X;
//   mat B_beta_new = inv_cpp(Binv_beta_new);
//   vec a_beta_new = B_beta_new * (Binv_beta * a_beta + X.t() * Zstar);
//   vec eta_new = rmvn_cpp(1, a_beta_new, B_beta_new);
//   return eta_new;
// }

// For a single binary covariate
// [[Rcpp::export]]
double update_omega_bin_cpp(vec const& X_tmp, double const& a_bin, double const& b_bin) {
  double sum_X_tmp = sum(X_tmp);
  double n_tmp = X_tmp.n_elem;
  double a_bin_new = a_bin + sum_X_tmp;
  double b_bin_new = b_bin + n_tmp - sum_X_tmp;
  return R::rbeta(a_bin_new, b_bin_new);
}

// For a single continuous covariate
// [[Rcpp::export]]
List update_omega_con_cpp(vec const& X_tmp,
                          double const& a_mean, double const& b_mean, // Prior for mean: a_mu, b_mu
                          double const& a_var, double const& b_var) {   // Prior for variance: a_tau, b_tau
  
  // 1. Calculate sufficient statistics
  int n_tmp = X_tmp.n_elem;
  double mean_X_tmp = mean(X_tmp);
  // (n_k - 1) * s_k^2
  double ss = (n_tmp > 1) ? accu(pow(X_tmp - mean_X_tmp, 2)) : 0.0; 
  
  // --- 2. Update Variance (tau^2) Parameters ---
  // (Using the Scaled-Inverse-Chi-Square equations from your LaTeX)
  // Posterior degrees of freedom: a_tau_n = a_tau + n_k
  double a_var_new = a_var + n_tmp;
  // Posterior scale numerator
  double weighted_mean_diff_sq = (b_mean * n_tmp / (b_mean + n_tmp)) * pow(mean_X_tmp - a_mean, 2);
  double numerator = a_var * b_var + ss + weighted_mean_diff_sq;
  // Posterior scale: b_tau_n
  double b_var_new = numerator / a_var_new;
  
  // 3. Sample new variance from its posterior
  double var_prop = rscainvchisq_cpp(1, a_var_new, b_var_new)[0];
  
  // --- 4. Update Mean (mu) Parameters ---
  // Posterior mean (a_mu_n)
  double b_mean_new = b_mean + n_tmp;
  double a_mean_new = (b_mean * a_mean + n_tmp * mean_X_tmp) / b_mean_new;
  
  // 5. Sample new mean from its posterior
  double mean_prop = R::rnorm(a_mean_new, sqrt(var_prop / b_mean_new));
  
  return List::create(_["var_par"] = var_prop,
                      _["mean_par"] = mean_prop);
}

// Define function to update alpha_theta prior Gamma(a_alpha,b_alpha)
// [[Rcpp::export]]
double update_alpha_cpp(int const& n, int const& K, double const& alpha_curr,
                        double const& a_alpha, double const& b_alpha) {
  double logxi = log(R::rbeta(alpha_curr + 1, n));
  double varrho = (a_alpha + K - 1)/(n * (b_alpha - logxi) + (a_alpha + K - 1));
  varrho = R::rbinom(1, varrho);
  
  // randg: shape and scale at default
  // rgamma: shape and rate at default
  double alpha_prop;
  if(varrho > 0){
    alpha_prop = randg(distr_param(a_alpha + K, 1/(b_alpha - logxi)));
  } else {
    alpha_prop = randg(distr_param(a_alpha + K - 1, 1/(b_alpha - logxi)));
  }
  
  return alpha_prop;
}

// Define function to update alpha_theta prior Gamma(a_theta,b_theta) with shape and rate parameters
// [[Rcpp::export]]
double update_alpha_theta_cpp(int const& N, int const& Ky, double const& alpha_theta_curr,
                              double const& a_theta, double const& b_theta) {
  double logxi = log(R::rbeta(alpha_theta_curr + 1, N));
  double a_theta_new = a_theta + Ky - 1;
  double b_theta_new = b_theta - logxi;
  double binv_theta_new = 1/b_theta_new;
  
  double varrho = a_theta_new/(N * b_theta_new + a_theta_new);
  varrho = R::rbinom(1, varrho);
  
  // randg: shape and scale at default in Rcpp & RcppArmadillo
  // rgamma: shape and rate at default in R
  double alpha_theta_prop;
  if (varrho > 0){
    alpha_theta_prop = randg(distr_param(a_theta_new + 1, binv_theta_new));
  } else {
    alpha_theta_prop = randg(distr_param(a_theta_new, binv_theta_new));
  }
  
  return alpha_theta_prop;
}

// Define function to update alpha_omega using Metropolis Hastings
// [[Rcpp::export]]
double update_alpha_omega_cpp(int const& Kyx, int const& Ky, uvec const& Sy, 
                              uvec const& unique_Sy, double const& alpha_omega_curr,
                              double const& a_omega, double const& b_omega) {
  double bivn_omega = 1/b_omega;
  double alpha_omega_prop = randg(distr_param(a_omega,bivn_omega));
  
  int n_temp = 0;
  double prod_prop = 1;
  double prod_curr = 1;
  for (int k = 0; k < Ky; k++) {
    uvec ind_Sy = find(Sy == unique_Sy(k));
    n_temp = ind_Sy.n_elem; 
    
    if (n_temp>0){
      prod_prop = prod_prop * (alpha_omega_prop + n_temp) * R::beta(alpha_omega_prop + 1, n_temp);
      prod_curr = prod_curr * (alpha_omega_curr + n_temp) * R::beta(alpha_omega_curr + 1, n_temp);
    }
  }
  
  double like_prop = R::dgamma(alpha_omega_prop, a_omega, bivn_omega, false) *
    (pow(alpha_omega_prop, (Kyx - Ky))) * prod_prop;
  double like_curr = R::dgamma(alpha_omega_curr, a_omega, bivn_omega, false) *
    (pow(alpha_omega_curr, (Kyx - Ky))) * prod_curr;
  
  double ratio = (like_prop)/(like_curr);
  if (ratio < randu(1)[0]){ alpha_omega_prop = alpha_omega_curr; }
  
  return alpha_omega_prop;
}

// =============================================================================
// LIKELIHOOD CALCULATION FOR CLUSTERING
// =============================================================================
// [[Rcpp::export]]
double calculate_log_lik_cpp(int const& p_X1, int const& p_X2,
                             double Y_tmp, rowvec const& matX_tmp, 
                             int Z_tmp, rowvec const& X_tmp,
                             vec const& beta_tmp, double const& sig2_tmp, 
                             double const& eta_tmp, vec const& pi_temp,
                             vec const& mu_temp, vec const& tau2_temp) {
  double log_lik = 0.0;
  
  log_lik += R::dnorm(Y_tmp, as_scalar(matX_tmp*beta_tmp), sqrt(sig2_tmp), true);
  log_lik += R::dbinom(Z_tmp, 1, eta_tmp, true); 
  // Likelihood for covariates X_{i} (binary)
  for (int q = 0; q < p_X1; ++q) { 
    log_lik += R::dbinom(X_tmp(q), 1, pi_temp(q), true); 
  }
  // Likelihood for covariates X_{i} (continuous)
  for (int q = 0; q < p_X2; ++q) { 
    log_lik += R::dnorm(X_tmp(p_X1 + q), mu_temp(q), sqrt(tau2_temp(q)), true); 
  }
  
  return log_lik;
}

// =============================================================================
// Survival_Y_cpp
// =============================================================================
// [[Rcpp::export]]
double S_DPMM_cpp(double const& y, vec const& lambda, 
                  vec const& mu, vec const& sig,
                  bool const& logt) {
  // Z-score: (y - mu) / sig
  vec z_scores = (y - mu) / sig;
  
  if (logt) {
    // --- log-probabiliy-space (comuptationally stable) ---
    // log(S(y)) = log(P(Y > y))
    vec log_S_k(z_scores.n_elem);
    for(size_t i = 0; i < z_scores.n_elem; ++i) {
      log_S_k(i) = R::pnorm(z_scores(i), 0.0, 1.0, false, true); // (lower.tail=F, log.p=T)
    }
    
    // Log-Sum-Exp trick
    vec log_terms = log(lambda) + log_S_k;
    double max_val = max(log_terms);
    double result = max_val + log(sum(exp(log_terms - max_val)));
    return result;
    
  } else {
    // --- probabiliy-space ---
    vec S_k(z_scores.n_elem);
    for(size_t i = 0; i < z_scores.n_elem; ++i) {
      S_k(i) = R::pnorm(z_scores(i), 0.0, 1.0, false, false); // (lower.tail=F, log.p=F)
    }
    return sum(lambda % S_k);
  }
  
  // // (1 - normcdf(y, mu, sig)) = normcdf(-y, -mu, sig)
  // if (logt) {
  //   return log(sum(lambda % normcdf(-y, -mu, sig)));
  // } else {
  //   return sum(lambda % normcdf(-y, -mu, sig));
  // }
}

// [[Rcpp::export]]
List S_optim_cpp(double const& target, double const& log_nu,// target = {S(nu)} * {1-rho}
                 double const& y_init, double const& y_min, double const& y_max,
                 vec const& lambda, vec const& Ymu, vec const& Ysig,
                 bool const& logt){
  double y_star = y_init;
  double y_m    = 0.0;
  double S_m    = 0.0;
  double step   = 0.0;
  
  double y_u = y_max;
  double y_l = std::max(y_min, log_nu + 1e-8);
  
  double S_u = S_DPMM_cpp(y_u, lambda, Ymu, Ysig, logt);
  double S_l = S_DPMM_cpp(y_l, lambda, Ymu, Ysig, logt);
  
  bool cond; int it; int maxit = 1e5; double tol = 1e-10;
  
  cond = (S_u >= target); it = 1;
  while (cond&& it < maxit){
    y_u += 1.0;
    S_u = S_DPMM_cpp(y_u, lambda, Ymu, Ysig, logt);
    cond = (S_u >= target);
    it += 1;
  }
  
  cond = (S_l <= target); it = 1;
  while (cond&& it < maxit){
    y_l -= 1.0;
    if (y_l <= log_nu) {
      y_l = log_nu + 1e-8;
      S_l = S_DPMM_cpp(y_l, lambda, Ymu, Ysig, logt);
      break;
    }
    S_l = S_DPMM_cpp(y_l, lambda, Ymu, Ysig, logt);
    cond = (S_l <= target);
    it += 1;
  }
  
  if (S_l > target && S_u < target){
    cond = true; it = 1; step = y_u - y_l;
    while (cond){
      y_m = (y_l + y_u)/2.0;
      S_m = S_DPMM_cpp(y_m, lambda, Ymu, Ysig, logt);
      if (S_m > target) {
        y_l = y_m;
      } else {
        y_u = y_m;
      }
      step /= 2.0;
      if (abs(S_m-target) < tol || step < tol || it > maxit) {
        cond = false;
      }
      it += 1;
    }
    y_star = y_m;
  } else {
    // if fail to converge, return initial value
    y_star = y_init; 
    S_m = -999.0;
    step = -0.1;
    it = -1;
  }
  
  return List::create(_["optimizer"] = y_star,
                      _["S_star"] = S_m,
                      _["step"] = step,
                      _["it"] = it,
                      _["maxit"] = maxit,
                      _["tol"] = tol);
}

// =============================================================================
// MAIN CLUSTER UPDATE FUNCTION
// =============================================================================
// [[Rcpp::export]]
List DPMMcluster_cpp(int const& n, int const& p_X1, int const& p_X2,
                     vec const& Y, mat const& matX,
                     vec const& Z, mat const& X,
                     uvec S, umat unique_S, double const& alpha,
                     mat& betaPars, mat& sig2Pars, mat& etaPars,
                     mat& piPars, mat& muPars, mat& tau2Pars,
                     vec const& a_beta, mat const& B_beta, mat const& Binv_beta,
                     double const& a_sig2, double const& b_sig2,
                     double const& a_eta, double const& b_eta,
                     double a_pi, double b_pi,
                     double a_mu, double b_mu,
                     double a_tau2, double b_tau2) {
  // initialize some dummy vars to store things
  uvec ind_S, ind_unique_S;
  int num_S;
  
  for (int i = 0; i < n; ++i) {
    // ---------------------------------------------------------------------------------------
    // ---------------------- Delete i-th person / Adjust the membership ---------------------
    // ---------------------------------------------------------------------------------------
    ind_S = find(S == S(i));
    num_S = ind_S.size();
    
    if(num_S==1) { //if lone person in X-Y cluster
      //DELETE ASSOCIATED COEFFICIENTS IN Y AND X CLUSTER
      
      ind_S = find(S == S(i)); //check if only person in Y cluster too
      num_S = ind_S.size();
      
      //Rprintf("num_S: %d\n",num_S);
      //Rprintf("S(i): %d\n", S(i));
      //Rcout << S.t() << std::endl;
      //Rcout << ind_S.t() << std::endl;
      //Rcout << ind_S.size() << std::endl;
      // delete Y coef if only one in Y cluster
      if(num_S==1) {
        betaPars.shed_col(S(i)-1);
        sig2Pars.shed_col(S(i)-1);
        etaPars.shed_col(S(i)-1);
        piPars.shed_col(S(i)-1);
        muPars.shed_col(S(i)-1);
        tau2Pars.shed_col(S(i)-1);
      }
      
      //should find row in unique_S that corresponds to person i
      ind_unique_S = find(unique_S.col(0)==S(i));
      int S_size = S.size();
      int unique_S_n_rows = unique_S.n_rows;
      // relabel Y cluster (if needed)
      if(num_S==1) {
        for(int k = 0; k < S_size; k++) {
          if(S(k) > S(i)){
            S(k) = S(k) - 1;
          }
        }
        for(int k = 0; k < unique_S_n_rows; k++) {
          if(unique_S(k,0) > S(i)) {
            unique_S(k,0) = unique_S(k,0) - 1;
          }
        }
      }
      unique_S.shed_row(ind_unique_S(0));
    }
    
    // NEED TO DELETE ROW OF S
    S.shed_row(i);
    
    // ---------------------------------------------------------------------------------------
    // --------------------- Calculate the draw probability for clusters ---------------------
    // ---------------------------------------------------------------------------------------
    // Find unique clusters and counts *without* subject i
    int K = betaPars.n_cols;
    vec log_probs = zeros(K + 1);
    int nkwoi; //counts for # in appropriate cluster
    
    // ---------------------------------------------------------------------------------------
    // ------------------ Calculate Log Probabilities for Existing Clusters ------------------
    // ---------------------------------------------------------------------------------------
    for (int k = 0; k < K; ++k) {
      ind_S = find(S == (k+1));
      nkwoi = ind_S.size();
      
      double log_lik = calculate_log_lik_cpp(p_X1, p_X2,
                                             Y(i), matX.row(i), 
                                             Z(i), X.row(i),
                                             betaPars.col(k),
                                             sig2Pars(k),
                                             etaPars(k),
                                             piPars.col(k),
                                             muPars.col(k),
                                             tau2Pars.col(k));
      log_probs(k) = log(nkwoi) + log_lik;
    }
    
    // ---------------------------------------------------------------------------------------
    // --------------------- Calculate Log Probability for a New Cluster ---------------------
    // ---------------------------------------------------------------------------------------
    // // Using MC integration
    // log_probs(K) = log(alpha) + log(f0_y(i)) + log(f0_z(i)) + log(f0_x(i));
    // // Generate new parameters based on the posterior conditional on subject i's data
    // vec beta_draw = rmvn_cpp(1, a_beta, B_beta);
    // List beta_sig_list = update_betasig_cpp(Y.subvec(i, i), 
    //                                         matX.row(i),
    //                                         beta_draw,
    //                                         a_beta, B_beta, Binv_beta,
    //                                         a_sig2, b_sig2);
    // vec beta_new = as<vec>(beta_sig_list["beta_par"]);
    // double sig2_new = beta_sig_list["sig2_par"];
    // double eta_new = update_omega_bin_cpp(Z.subvec(i, i), a_eta, b_eta);
    // 
    // vec pi_new(p_X1);
    // vec mu_new(p_X2);
    // vec tau2_new(p_X2);
    // // For binary covariates
    // for (int q = 0; q < p_X1; ++q) {
    //   pi_new(q) = update_omega_bin_cpp(X.submat(i, q, i, q), a_pi, b_pi);
    // }
    // // For continuous covariates
    // for (int q = 0; q < p_X2; ++q) {
    //   List con_list = update_omega_con_cpp(X.submat(i, p_X1 + q, i, p_X1 + q),
    //                                        a_mu, b_mu, a_tau2, b_tau2);
    //   mu_new(q) = as<double>(con_list["mean_par"]);
    //   tau2_new(q) = as<double>(con_list["var_par"]);
    // }
    // Draw Auxiliary Parameters from Prior (Base Measure G0)
    vec beta_new = rmvn_cpp(1, a_beta, B_beta);
    double sig2_new = rscainvchisq_cpp(1, a_sig2, b_sig2)[0];
    double eta_new = rscainvchisq_cpp(1, a_eta, b_eta)[0];
    vec pi_new(p_X1);
    vec mu_new(p_X2);
    vec tau2_new(p_X2);
    for (int q = 0; q < p_X1; ++q) {
      pi_new(q) = R::rbeta(a_pi, b_pi);
    }
    for (int q = 0; q < p_X2; ++q) {
      tau2_new(q) = rscainvchisq_cpp(1, a_tau2, b_tau2)[0];
      mu_new(q) = R::rnorm(a_mu, sqrt(tau2_new(q) / b_mu));
    }
    
    double log_lik_new = calculate_log_lik_cpp(p_X1, p_X2,
                                               Y(i), matX.row(i), 
                                               Z(i), X.row(i),
                                               beta_new, sig2_new, eta_new,
                                               pi_new, mu_new, tau2_new);
    log_probs(K) = log(alpha) + log_lik_new;
    
    // ---------------------------------------------------------------------------------------
    // ------------------- Normalize and Sample a New Cluster for Subject i ------------------
    // ---------------------------------------------------------------------------------------
    vec probs = exp(log_probs - log_probs.max());
    probs = probs / sum(probs);
    uword k_new = rmultinom_cpp(probs);
    
    // ---------------------------------------------------------------------------------------
    // ------------------------------------ Update State -------------------------------------
    // ---------------------------------------------------------------------------------------
    if (k_new <= K) {
      // Assign to an EXISTING cluster. No parameter updates happen here.
      S.insert_rows(i,1);
      S(i) = unique_S(k_new - 1, 0);
      
    } else {
      betaPars.insert_cols(K, beta_new);
      sig2Pars.insert_cols(K, 1); 
      sig2Pars(K) = sig2_new;
      etaPars.insert_cols(K, 1);
      etaPars(K) = eta_new;
      piPars.insert_cols(K, pi_new);
      muPars.insert_cols(K, mu_new);
      tau2Pars.insert_cols(K, tau2_new);
      
      S.insert_rows(i, 1);
      S(i) = S.max() + 1;
      unique_S.insert_rows(K, 1);
      unique_S(K,0) = S(i);
    }
  }
  
  unique_S = unique(unique_S);
  
  return List::create(
    _["S"] = S,
    _["unique_S"] = unique_S,
    _["betaPars"] = betaPars,
    _["sig2Pars"] = sig2Pars,
    _["etaPars"]  = etaPars,
    _["piPars"]   = piPars,
    _["muPars"]   = muPars,
    _["tau2Pars"] = tau2Pars
  );
}

// [[Rcpp::export]]
List EDPMMcluster_cpp(int const& n, int const& p_X1, int const& p_X2,
                      vec const& Y, mat const& matX,
                      vec const& Z, mat const& X,
                      uvec Sy, uvec Sx, umat unique_Syx, 
                      double const& alpha_theta, double const& alpha_omega,
                      mat& betaPars, mat& sig2Pars, mat& etaPars,
                      mat& piPars, mat& muPars, mat& tau2Pars,
                      vec const& a_beta, mat const& B_beta, mat const& Binv_beta,
                      double const& a_sig2, double const& b_sig2,
                      double const& a_eta, double const& b_eta,
                      double const& a_pi, double const& b_pi,
                      double const& a_mu, double const& b_mu,
                      double const& a_tau2, double const& b_tau2,
                      vec const& f0_zx) {
  
  // Initialize basic counters
  int Kyx, Ky; 
  int num_Sy, num_Syx, num_Sy0;
  uvec ind_Sy, ind_Sy0, ind_Syx, ind_unique_Syx;
  double Ylik, Zlik, Xlik;
  double Y_temp, Z_temp;
  rowvec matX_temp, X_temp;
  uword Sy_temp, Sx_temp;
  
  // Proposal containers
  vec beta_prop; 
  double sig2_prop;
  double eta_prop;
  vec pi_prop(p_X1);
  vec mu_prop(p_X2);
  vec tau2_prop(p_X2);
  
  // ---------------------------------------------------------------------------
  // GIBBS SAMPLING: Loop through every person
  // ---------------------------------------------------------------------------
  
  // Loop through every person
  for (int i = 0; i < n; i++) {
    Y_temp = Y(i);
    Z_temp = Z(i);
    X_temp = X.row(i);
    matX_temp = matX.row(i);
    
    Sy_temp = Sy(i);
    Sx_temp = Sx(i);
    
    // -----------------------------------------------------------------------
    // 2. Remove i-th person / Adjust membership
    // -----------------------------------------------------------------------
    // Check if i is lone person in the specific nested cluster (Sy, Sx)
    ind_Syx = find(Sy == Sy_temp && Sx == Sx_temp);
    num_Syx = ind_Syx.size();
    
    if (num_Syx == 1) {
      // Check if i is lone person in the top-level Y cluster
      ind_Sy = find(Sy == Sy_temp);
      num_Sy = ind_Sy.size();
      
      if (num_Sy == 1) {
        // Delete Y-level parameters (associated with Sy index)
        betaPars.shed_col(Sy_temp - 1);
        sig2Pars.shed_col(Sy_temp - 1);
      }
      
      // Find row in unique_Syx corresponding to (Sy_temp, Sx_temp)
      ind_unique_Syx = find(unique_Syx.col(0) == Sy_temp && unique_Syx.col(1) == Sx_temp);
      int unique_Syx_n_rows = unique_Syx.n_rows;
      
      // Delete X-level parameters (associated with unique_Syx row)
      etaPars.shed_col(ind_unique_Syx(0));
      piPars.shed_col(ind_unique_Syx(0));
      muPars.shed_col(ind_unique_Syx(0));
      tau2Pars.shed_col(ind_unique_Syx(0));
      
      // Relabel X clusters (Sx)
      for (int ii = 0; ii < n; ii++) {
        if (Sy(ii) == Sy_temp && Sx(ii) > Sx_temp) {
          Sx(ii) = Sx(ii) - 1;
        }
      }
      for (int k = 0; k < unique_Syx_n_rows; k++) {
        if (unique_Syx(k, 0) == Sy_temp && unique_Syx(k, 1) > Sx_temp) {
          unique_Syx(k, 1) = unique_Syx(k, 1) - 1;
        }
      }
      
      // Relabel Y clusters (Sy)
      if (num_Sy == 1) {
        for (int ii = 0; ii < n; ii++) {
          if (Sy(ii) > Sy_temp) Sy(ii) = Sy(ii) - 1;
        }
        for (int k = 0; k < unique_Syx_n_rows; k++) {
          if (unique_Syx(k, 0) > Sy_temp) unique_Syx(k, 0) = unique_Syx(k, 0) - 1;
        }
      }
      
      // Remove row from unique tracker
      unique_Syx.shed_row(ind_unique_Syx(0));
    }
    
    // Delete row from membership vectors temporarily
    Sy.shed_row(i); Sx.shed_row(i);
    
    // -----------------------------------------------------------------------
    // 3. Calculate Draw Probabilities
    // -----------------------------------------------------------------------
    Ky = betaPars.n_cols;      // Number of Y clusters
    Kyx = etaPars.n_cols;      // Number of total X sub-clusters
    int K_totalposs = Kyx + Ky + 1; // Existing X + New X in Existing Y + New Y
    
    vec probs(K_totalposs);
    int count = 0;
    int n_k_woi, n_rk_woi;
    
    // // Using MC integration
    // log_probs(K) = log(alpha) + log(f0_y(i)) + log(f0_z(i)) + log(f0_x(i));
    // // Generate new parameters based on the posterior conditional on subject i's data
    // vec beta_draw_temp = rmvn_cpp(1, a_beta, B_beta); 
    // List beta_sig_prop_list = update_betasig_cpp(Y.subvec(i, i), 
    //                                              matX.row(i), 
    //                                              beta_draw_temp,
    //                                              a_beta, B_beta, Binv_beta, 
    //                                              a_sig2, b_sig2);
    // beta_prop = as<vec>(beta_sig_prop_list["beta_par"]);
    // sig2_prop = beta_sig_prop_list["sig2_par"];
    // eta_prop = update_omega_bin_cpp(Z.subvec(i, i), a_eta, b_eta);
    // for (int q = 0; q < p_X1; ++q) {
    //   pi_prop(q) = update_omega_bin_cpp(X.submat(i, q, i, q), a_pi, b_pi);
    // }
    // for (int q = 0; q < p_X2; ++q) {
    //   List con_list_prop = update_omega_con_cpp(X.submat(i, p_X1 + q, i, p_X1 + q),
    //                                             a_mu, b_mu, a_tau2, b_tau2);
    //   mu_prop(q)   = as<double>(con_list_prop["mean_par"]);
    //   tau2_prop(q) = as<double>(con_list_prop["var_par"]);
    // }
    // Draw Auxiliary Parameters from Prior (Base Measure G0)
    beta_prop = rmvn_cpp(1, a_beta, B_beta);
    sig2_prop = rscainvchisq_cpp(1, a_sig2, b_sig2)[0];
    eta_prop = R::rbeta(a_eta, b_eta);
    for (int q = 0; q < p_X1; ++q) {
      pi_prop(q) = R::rbeta(a_pi, b_pi);
    }
    for (int q = 0; q < p_X2; ++q) {
      tau2_prop(q) = rscainvchisq_cpp(1, a_tau2, b_tau2)[0];
      mu_prop(q) = R::rnorm(a_mu, sqrt(tau2_prop(q) / b_mu));
    }
    
    // --- Calculate Ylik_prop ---
    double Ylik_prop = R::dnorm(Y_temp, as_scalar(matX_temp * beta_prop), sqrt(sig2_prop), false);
    // --- Calculate XZlik_prop ---
    double XZlik_prop = R::dbinom(Z_temp, 1, eta_prop, false);
    for (int q = 0; q < p_X1; q++) {
      XZlik_prop *= R::dbinom(X_temp(q), 1, pi_prop(q), false);
    }
    for (int q = 0; q < p_X2; q++) {
      XZlik_prop *= R::dnorm(X_temp(p_X1 + q), mu_prop(q), sqrt(tau2_prop(q)), false);
    }
    
    // Loop through existing Y clusters
    for (int k = 0; k < Ky; k++) {
      // Get count in Y cluster
      ind_Sy = find(Sy == (k + 1));
      n_k_woi = ind_Sy.size();
      
      // Likelihood for Y (Outer cluster)
      // Note: sig2Pars is 1xKy, betaPars is p x Ky
      Ylik = R::dnorm(Y_temp, as_scalar(matX_temp * betaPars.col(k)), sqrt(sig2Pars(0, k)), false);
      
      // Loop through X sub-clusters within this Y cluster
      ind_Sy0 = find(unique_Syx.col(0) == (k + 1));
      int Kx_Sy = ind_Sy0.size();

      for (int r = 0; r < Kx_Sy; r++) {
        // Get count in nested X cluster
        ind_Syx = find(Sy == (k + 1) && Sx == (r + 1));
        n_rk_woi = ind_Syx.size();
        
        // Likelihood for Z (Treatment)
        Zlik = R::dbinom(Z_temp, 1, etaPars(0, count), false);
        
        // Likelihood for X (Covariates)
        Xlik = 1.0;
        int bin_idx = 0; 
        int con_idx = 0;
        for (int q = 0; q < (p_X1 + p_X2); ++q) {
          // Note:Assuming X is ordered binary then continuous based on input
          if (q < p_X1) {
            Xlik *= R::dbinom(X_temp(q), 1, piPars(bin_idx, count), false);
            bin_idx++;
          } else {
            Xlik *= R::dnorm(X_temp(q), muPars(con_idx, count), sqrt(tau2Pars(con_idx, count)), false);
            con_idx++;
          }
        }
        
        // Probability for Existing Y, Existing X
        probs(count) = (n_k_woi * (n_rk_woi / (n_k_woi + alpha_omega))) * Ylik * Zlik * Xlik;
        count++;
      }
      probs(Kyx + k) = (n_k_woi * (alpha_omega / (n_k_woi + alpha_omega))) * Ylik * XZlik_prop;
      // probs(Kyx + k) = (n_k_woi * (alpha_omega / (n_k_woi + alpha_omega))) * Ylik * f0_zx(i);
    }
    probs(Kyx + Ky) = alpha_theta * Ylik_prop * XZlik_prop;
    // probs(Kyx + Ky) = alpha_theta * f0_y(i) * f0_zx(i);
    
    // -----------------------------------------------------------------------
    // 4. Sample New Cluster
    // -----------------------------------------------------------------------
    probs = probs / sum(probs);
    int k_new = rmultinom_cpp(probs); // 1-based index
    
    // -----------------------------------------------------------------------
    // 5. Update State
    // -----------------------------------------------------------------------
    
    // CASE 1: Existing Y, Existing X
    // Assign
    if (k_new <= Kyx) {
      Sy.insert_rows(i, 1); 
      Sy(i) = unique_Syx(k_new - 1, 0);
      Sx.insert_rows(i, 1); 
      Sx(i) = unique_Syx(k_new - 1, 1);
    } else {
      if (k_new == K_totalposs){
        // New Y, New X
        Sx.insert_rows(i,1);
        Sx(i) = 1;
        
        betaPars.insert_cols(Ky, beta_prop);
        sig2Pars.insert_cols(Ky, 1);
        sig2Pars(0, Ky) = sig2_prop;
        etaPars.insert_cols(Kyx, 1);
        etaPars(0, Kyx) = eta_prop;
        piPars.insert_cols(Kyx, pi_prop);
        muPars.insert_cols(Kyx, mu_prop);
        tau2Pars.insert_cols(Kyx, tau2_prop);
        
        Sy.insert_rows(i,1);
        Sy(i) = Sy.max() + 1;
        unique_Syx.insert_rows(Kyx, 1); 
        unique_Syx(Kyx, 0) = Sy(i); 
        unique_Syx(Kyx, 1) = Sx(i);
      } else {
        // Existing Y, New X
        Sy.insert_rows(i,1);
        Sy(i) = k_new - Kyx;
        ind_Sy0 = find(unique_Syx.col(0) == Sy(i));
        Sx.insert_rows(i,1);
        Sx(i) = ind_Sy0.size() + 1;
        
        ind_Sy0 = find(unique_Syx.col(0) <= Sy(i));
        num_Sy0 = ind_Sy0.size();
        
        etaPars.insert_cols(num_Sy0, 1);
        etaPars(0, num_Sy0) = eta_prop;
        piPars.insert_cols(num_Sy0, pi_prop);
        muPars.insert_cols(num_Sy0, mu_prop);
        tau2Pars.insert_cols(num_Sy0, tau2_prop);
        
        unique_Syx.insert_rows(num_Sy0, 1); 
        unique_Syx(num_Sy0, 0) = Sy(i); 
        unique_Syx(num_Sy0, 1) = Sx(i);
      }
    }
  }
  
  uvec unique_Sy = unique(Sy);
  umat Syx = join_rows(Sy, Sx);
  unique_Syx = unique_rows(Syx);
  
  Ky = unique_Sy.size();
  Kyx = unique_Syx.n_rows;
  int max_Sx = Sx.max();
  
  vec n_k = zeros(Ky);
  mat n_rk = zeros(Ky,max_Sx);
  vec max_Kx_Sy = zeros(Ky);
  for (int k = 0; k < Ky; k++) {
    ind_Sy = find(Sy == (k+1));
    n_k(k) = ind_Sy.size();
    
    ind_Sy0 = find(unique_Syx.col(0) == (k+1));
    num_Sy0 = ind_Sy0.size();
    max_Kx_Sy(k) = num_Sy0;
    for (int r = 0; r < num_Sy0; r++) {
      ind_Syx = find(Sx(ind_Sy) == (r+1));
      n_rk(k,r) = ind_Syx.size();
    }
  }
  
  // return List::create(
  //   _["Syx"] = Syx,
  //   _["unique_Sy"] = unique_Sy,
  //   _["unique_Syx"] = unique_Syx,
  //   _["betaPars"] = betaPars,
  //   _["sig2Pars"] = sig2Pars,
  //   _["etaPars"] = etaPars,
  //   _["piPars"] = piPars,
  //   _["muPars"] = muPars,
  //   _["tau2Pars"] = tau2Pars,
  //   _["Ky"] = Ky,
  //   _["Kyx"] = Kyx,
  //   _["n_k"] = n_k,
  //   _["n_rk"] = n_rk,
  //   _["max_Sx"] = max_Sx,
  //   _["max_Kx_Sy"] = max_Kx_Sy);
  
  // Rcout <<
  //   "Ky: " << Ky << "\n" <<
  //     "Kyx: " << Kyx << "\n" <<
  //       "n_k: " << n_k << "\n" <<
  //         "n_rk: " << n_rk << "\n" <<
  //           "max_Kx_Sy: " << max_Kx_Sy <<
  //             endl;
  //
  // Rcout << "Improve Mixing (Metropolis-Hastings)" << endl;

  // ------------------------------------------------------------------------------
  // -------------------- Improve Mixing (Metropolis-Hastings) --------------------
  // ------------------------------------------------------------------------------
  // initialize some dummy vars to store things
  double ratio;

  int num_max_Kx_Sy, num_max_Kx_Sy_temp, num_unique_Syx;
  vec max_Kx_Sy_temp, cumsum_max_Kx_Sy;

  uvec ind_max_Kx_Sy, ind_max_Kx_Sy_temp, uvec_cumsum_max_Kx_Sy, unique_Sy_temp;
  umat unique_Syx_temp;

  int n_k_curr, n_k_prop, n_rk_curr;
  int Kx1_curr, Kx1_prop, Kx2_curr, Kx2_prop;
  int k_curr, r_curr, k_prop, r_prop;
  int k_curr_1, r_curr_1, k_prop_1;

  // ---------------------------------------------------------------------------
  // 0. Pre-calculate Marginal f0(Y) via Monte Carlo
  // ---------------------------------------------------------------------------
  vec f0_y = zeros(n);
  int num_MC_prior = 1e4; // You can pass this as an argument if preferred
  for (int m = 0; m < num_MC_prior; ++m) {
    double sig2_draw = rscainvchisq_cpp(1, a_sig2, b_sig2)[0];
    vec beta_draw = rmvn_cpp(1, a_beta, B_beta);
    vec mu_vec = matX * beta_draw;
    double sd_draw = sqrt(sig2_draw);
    for(int i = 0; i < n; ++i) {
      f0_y(i) += R::dnorm(Y(i), mu_vec(i), sd_draw, false);
    }
  }
  f0_y = f0_y / num_MC_prior;

  // Rcout << "1st" << endl;
  // ---------------------------------- 1st move ----------------------------------
  if (Kyx > Ky && Ky > 1){
    ind_max_Kx_Sy = find(max_Kx_Sy == 1);
    num_max_Kx_Sy = ind_max_Kx_Sy.size();

    // Kx2_curr = Kyx - num_max_Kx_Sy;
    //
    // // Choose cluster (k,r) and h
    // int choose_kr_temp = rmultinom_cpp(ones(Kx2_curr));
    // int choose_new_k_temp = rmultinom_cpp(ones(Ky-1));
    //
    // unique_Syx_temp = unique_Syx;
    // if (num_max_Kx_Sy > 0) {
    //   ind_max_Kx_Sy = find(max_Kx_Sy != 1);
    //   cumsum_max_Kx_Sy = cumsum(max_Kx_Sy) - ones(Ky);
    //   cumsum_max_Kx_Sy = cumsum_max_Kx_Sy(ind_max_Kx_Sy);
    //   uvec_cumsum_max_Kx_Sy = conv_to<uvec>::from(cumsum_max_Kx_Sy);
    //   unique_Syx_temp = unique_Syx_temp.rows(uvec_cumsum_max_Kx_Sy);
    // }

    uvec rows_to_keep;
    for (unsigned int k = 0; k < unique_Syx.n_rows; k++) {
      int k_temp = unique_Syx(k, 0);
      if (max_Kx_Sy(k_temp - 1) > 1) {
        rows_to_keep = join_cols(rows_to_keep, uvec{k});
      }
    }
    unique_Syx_temp = unique_Syx.rows(rows_to_keep);
    Kx2_curr = unique_Syx_temp.n_rows;

    int choose_kr_temp = rmultinom_cpp(ones(Kx2_curr));
    int choose_new_k_temp = rmultinom_cpp(ones(Ky - 1));

    k_curr = unique_Syx_temp(choose_kr_temp-1,0);
    r_curr = unique_Syx_temp(choose_kr_temp-1,1);
    k_curr_1 = k_curr - 1;
    r_curr_1 = r_curr - 1;

    unique_Sy_temp = unique_Sy;
    unique_Sy_temp.shed_row(k_curr_1);
    k_prop = unique_Sy_temp(choose_new_k_temp-1);

    ind_unique_Syx = find(unique_Syx.col(0) == k_prop);
    num_unique_Syx = ind_unique_Syx.size();
    r_prop = num_unique_Syx + 1;
    k_prop_1 = k_prop - 1;

    // Calculate K_{x,2+}^{*}
    max_Kx_Sy_temp = zeros(Ky+1);
    max_Kx_Sy_temp.head(Ky) = max_Kx_Sy;
    max_Kx_Sy_temp(k_curr_1) = max_Kx_Sy_temp(k_curr_1) - 1;
    max_Kx_Sy_temp(k_prop_1) = max_Kx_Sy_temp(k_prop_1) + 1;
    ind_max_Kx_Sy_temp = find(max_Kx_Sy_temp == 1);
    num_max_Kx_Sy_temp = ind_max_Kx_Sy_temp.size();

    Kx2_prop = Kyx - num_max_Kx_Sy_temp;

    // Calculate the acceptance probability
    n_k_prop = n_k(k_prop_1);

    n_k_curr = n_k(k_curr_1);
    n_rk_curr = n_rk(k_curr_1, r_curr_1);

    ind_Syx = find(Sy == k_curr && Sx == r_curr);
    num_Syx = ind_Syx.size();

    vec Y_temp = Y(ind_Syx);
    mat matX_temp = matX.rows(ind_Syx);

    vec beta_curr = betaPars.col(k_curr_1);
    vec sig_curr = sqrt(sig2Pars(k_curr_1)) * ones(num_Syx);

    vec beta_prop = betaPars.col(k_prop_1);
    vec sig_prop = sqrt(sig2Pars(k_prop_1)) * ones(num_Syx);

    ratio =
      ((lgamma(n_k_curr - n_rk_curr) + lgamma(n_k_prop + n_rk_curr)) -
      (lgamma(n_k_curr) + lgamma(n_k_prop))) +
      ((lgamma(alpha_omega + n_k_curr) + lgamma(alpha_omega + n_k_prop)) -
      (lgamma(alpha_omega + n_k_curr - n_rk_curr) + lgamma(alpha_omega + n_k_prop + n_rk_curr))) +
      as_scalar(sum(log_normpdf(Y_temp, matX_temp * beta_prop, sig_prop)) -
      sum(log_normpdf(Y_temp, matX_temp * beta_curr, sig_curr))) +
      log(Kx2_curr/Kx2_prop);
    ratio = exp(ratio);

    // If ratio > runif (1), retrun "proposed" and o.w.,do nothing, i.e. return "current".
    if (ratio > randu(1)[0]){
      // reorder omega parameters
      uvec ind_unique_Syx_curr = find(unique_Syx.col(0)==k_curr && unique_Syx.col(1)==r_curr);
      uvec ind_Sy0_prop = find(unique_Syx.col(0) <= k_prop);
      int num_Sy0_prop = ind_Sy0_prop.size();

      etaPars = reorder(etaPars,ind_unique_Syx_curr(0),num_Sy0_prop);
      piPars  = reorder(piPars,ind_unique_Syx_curr(0),num_Sy0_prop);
      muPars  = reorder(muPars,ind_unique_Syx_curr(0),num_Sy0_prop);
      tau2Pars = reorder(tau2Pars,ind_unique_Syx_curr(0),num_Sy0_prop);

      // relabel clusters
      ind_Syx = find(Sy == k_curr && Sx == r_curr);
      num_Syx = ind_Syx.size();
      for (int ii = 0; ii < num_Syx; ii++) {
        Sy(ind_Syx(ii)) = k_prop;
        Sx(ind_Syx(ii)) = r_prop;
      }

      // relabel X cluster
      ind_Syx = find(Sy == k_curr && Sx >= r_curr);
      num_Syx = ind_Syx.size();
      for (int ii = 0; ii < num_Syx; ii++) {
        Sx(ind_Syx(ii)) -= 1;
      }

      // wrap-up
      Syx = join_rows(Sy, Sx);

      // Determine unique clusters from the cluster membership variable, Syx
      unique_Syx = unique_rows(Syx);

      // Make vector of y-clusters
      unique_Sy = unique(Sy);

      // Calculate the number of clusters
      Kyx = unique_Syx.n_rows;

      // Calculate the number of y-clusters
      Ky = unique_Sy.size();

      // Find the largest number of x clusters
      max_Sx = Sx.max();

      // Calculate the number of subjects in each y-cluster and store in vector n_k.
      // Calculate the number of subjects in each x-cluster and store in matrix n_rk.
      // Use k to index y-clusters and r to index x-clusters.
      n_k.zeros(Ky);
      n_rk.zeros(Ky, max_Sx);
      max_Kx_Sy.zeros(Ky);
      for (int k = 0; k < Ky; k++) {
        ind_Sy = find(Sy == (k+1));
        n_k(k) = ind_Sy.size();

        ind_Sy0 = find(unique_Syx.col(0) == (k+1));
        num_Sy0 = ind_Sy0.size();
        max_Kx_Sy(k) = num_Sy0;
        for (int r = 0; r < num_Sy0; r++) {
          ind_Syx = find(Sx(ind_Sy) == (r+1));
          n_rk(k,r) = ind_Syx.size();
        }
      }
    }
  }

  // Rcout << "2nd or 3rd" << endl;
  // 2nd/3rd Move
  if (randu(1)[0] < 0.5) {
    // ---------------------------------- 2nd move ----------------------------------
    if (Kyx > Ky) {
      // Calculate K_{x,1} and K_{x,2+}
      ind_max_Kx_Sy = find(max_Kx_Sy == 1);
      num_max_Kx_Sy = ind_max_Kx_Sy.size();

      // Kx2_curr = Kyx - num_max_Kx_Sy;
      //
      // // Choose cluster (k,r) and h
      // int choose_kr_temp = rmultinom_cpp(ones(Kx2_curr));
      //
      // unique_Syx_temp = unique_Syx;
      // if (num_max_Kx_Sy > 0) {
      //   ind_max_Kx_Sy = find(max_Kx_Sy != 1);
      //   cumsum_max_Kx_Sy = cumsum(max_Kx_Sy) - ones(Ky);
      //   cumsum_max_Kx_Sy = cumsum_max_Kx_Sy(ind_max_Kx_Sy);
      //   uvec_cumsum_max_Kx_Sy = conv_to<uvec>::from(cumsum_max_Kx_Sy);
      //   unique_Syx_temp = unique_Syx_temp.rows(uvec_cumsum_max_Kx_Sy);
      // }

      uvec rows_to_keep;
      for (unsigned int k = 0; k < unique_Syx.n_rows; k++) {
        int k_temp = unique_Syx(k, 0);
        if (max_Kx_Sy(k_temp - 1) > 1) {
          rows_to_keep = join_cols(rows_to_keep, uvec{k});
        }
      }
      unique_Syx_temp = unique_Syx.rows(rows_to_keep);
      Kx2_curr = unique_Syx_temp.n_rows;

      int choose_kr_temp = rmultinom_cpp(ones(Kx2_curr));
      k_curr = unique_Syx_temp(choose_kr_temp-1,0);
      r_curr = unique_Syx_temp(choose_kr_temp-1,1);
      k_prop = Ky + 1;
      r_prop = 1;

      k_curr_1 = k_curr - 1;
      r_curr_1 = r_curr - 1;
      k_prop_1 = k_prop - 1;

      // Calculate K_{x,1}^{*}
      max_Kx_Sy_temp = zeros(Ky+1);
      max_Kx_Sy_temp.head(Ky) = max_Kx_Sy;
      max_Kx_Sy_temp(k_curr_1) = max_Kx_Sy_temp(k_curr_1) - 1;
      max_Kx_Sy_temp(k_prop_1) = max_Kx_Sy_temp(k_prop_1) + 1;
      ind_max_Kx_Sy_temp = find(max_Kx_Sy_temp == 1);
      num_max_Kx_Sy_temp = ind_max_Kx_Sy_temp.size();

      Kx1_prop = num_max_Kx_Sy_temp;

      // Calculate the acceptance probability
      n_k_curr = n_k(k_curr_1);
      n_rk_curr = n_rk(k_curr_1, r_curr_1);

      ind_Syx = find(Sy == k_curr && Sx == r_curr);
      num_Syx = ind_Syx.size();

      vec Y_temp = Y(ind_Syx);
      mat matX_temp = matX.rows(ind_Syx);

      vec beta_curr = betaPars.col(k_curr_1);
      vec sig_curr = sqrt(sig2Pars(k_curr_1)) * ones(num_Syx);

      ratio =
        (lgamma(n_k_curr - n_rk_curr) + lgamma(n_rk_curr) - lgamma(n_k_curr)) +
        ((lgamma(alpha_omega + n_k_curr) + lgamma(alpha_omega)) -
        (lgamma(alpha_omega + n_k_curr - n_rk_curr) + lgamma(alpha_omega + n_rk_curr))) +
        as_scalar((sum(log(f0_y(ind_Syx)))) -
        (sum(log_normpdf(Y_temp, matX_temp * beta_curr, sig_curr)))) +
        log(Kx2_curr/(Kx1_prop*Kyx)) + log(alpha_theta);
      ratio = exp(ratio);

      // If ratio > runif (1), retrun "proposed" and o.w.,do nothing, i.e. retrun "current".
      if (ratio > randu(1)[0]){
        // Add theta parameters
        beta_prop = rmvn_cpp(1, a_beta, B_beta).col(0);
        betaPars.insert_cols(Ky, beta_prop);

        sig2_prop = rscainvchisq_cpp(1, a_sig2, b_sig2)[0];
        sig2Pars.insert_cols(Ky, 1);
        sig2Pars(0, Ky) = sig2_prop;

        // reorder omega parameters
        uvec ind_unique_Syx_curr = find(unique_Syx.col(0) == k_curr && unique_Syx.col(1) == r_curr);
        etaPars  = reorder(etaPars, ind_unique_Syx_curr(0), Kyx - 1);
        piPars   = reorder(piPars, ind_unique_Syx_curr(0), Kyx - 1);
        muPars   = reorder(muPars, ind_unique_Syx_curr(0), Kyx - 1);
        tau2Pars = reorder(tau2Pars, ind_unique_Syx_curr(0), Kyx - 1);

        // relabel clusters
        ind_Syx = find(Sy == k_curr && Sx == r_curr);
        num_Syx = ind_Syx.size();
        for (int ii = 0; ii < num_Syx; ii++) {
          Sy(ind_Syx(ii)) = k_prop;
          Sx(ind_Syx(ii)) = r_prop;
        }

        // relabel X cluster
        ind_Syx = find(Sy == k_curr && Sx > r_curr);
        num_Syx = ind_Syx.size();
        for (int ii = 0; ii < num_Syx; ii++) {
          Sx(ind_Syx(ii)) -= 1;
        }

        // wrap-up
        Syx = join_rows(Sy, Sx);

        // Determine unique clusters from the cluster membership variable, Syx
        unique_Syx = unique_rows(Syx);

        // Make vector of y-clusters
        unique_Sy = unique(Sy);

        // Calculate the number of clusters
        Kyx = unique_Syx.n_rows;

        // Calculate the number of y-clusters
        Ky = unique_Sy.size();

        // Find the largest number of x clusters
        max_Sx = Sx.max();

        // Calculate the number of subjects in each y-cluster and store in vector n_k.
        // Calculate the number of subjects in each x-cluster and store in matrix n_rk.
        // Use k to index y-clusters and r to index x-clusters.
        n_k.zeros(Ky);
        n_rk.zeros(Ky, max_Sx);
        max_Kx_Sy.zeros(Ky);
        for (int k = 0; k < Ky; k++) {
          ind_Sy = find(Sy == (k+1));
          n_k(k) = ind_Sy.size();

          ind_Sy0 = find(unique_Syx.col(0) == (k+1));
          num_Sy0 = ind_Sy0.size();
          max_Kx_Sy(k) = num_Sy0;
          for (int r = 0; r < num_Sy0; r++) {
            ind_Syx = find(Sx(ind_Sy) == (r+1));
            n_rk(k,r) = ind_Syx.size();
          }
        }
      }
    }
  } else {
    // ---------------------------------- 3rd move ----------------------------------
    // Calculate K_{x,1} and K_{x,2+}
    ind_max_Kx_Sy = find(max_Kx_Sy == 1);
    num_max_Kx_Sy = ind_max_Kx_Sy.size();

    // Initialize Kx1_curr safely
    Kx1_curr = 0;

    if (num_max_Kx_Sy>0 && Ky>1){
      Kx1_curr = num_max_Kx_Sy;

      // Choose cluster (k,r) and h
      // int choose_kr_temp = rmultinom_cpp(ones(Kx1_curr));
      // int choose_new_k_temp = rmultinom_cpp(ones(Ky-1));
      //
      // cumsum_max_Kx_Sy = cumsum(max_Kx_Sy) - ones(Ky);
      // cumsum_max_Kx_Sy = cumsum_max_Kx_Sy(ind_max_Kx_Sy);
      // uvec_cumsum_max_Kx_Sy = conv_to<uvec>::from(cumsum_max_Kx_Sy);
      // unique_Syx_temp = unique_Syx.rows(uvec_cumsum_max_Kx_Sy);
      //
      // k_curr = unique_Syx_temp(choose_kr_temp-1,0);
      uvec rows_to_keep;
      for (unsigned int k = 0; k < unique_Syx.n_rows; k++) {
        int k_temp = unique_Syx(k, 0);
        if (max_Kx_Sy(k_temp - 1) == 1) {
          rows_to_keep = join_cols(rows_to_keep, uvec{k});
        }
      }
      if (rows_to_keep.n_elem > 0) {
        unique_Syx_temp = unique_Syx.rows(rows_to_keep);
        Kx1_curr = unique_Syx_temp.n_rows;

        int choose_kr_temp = rmultinom_cpp(ones(Kx1_curr));
        k_curr = unique_Syx_temp(choose_kr_temp - 1, 0);
        r_curr = 1;
        k_curr_1 = k_curr - 1;
        r_curr_1 = r_curr - 1;

        int choose_new_k_temp = rmultinom_cpp(ones(Ky - 1));
        unique_Sy_temp = unique_Sy;
        unique_Sy_temp.shed_row(k_curr_1);
        k_prop = unique_Sy_temp(choose_new_k_temp-1);

        ind_unique_Syx = find(unique_Syx.col(0) == k_prop);
        num_unique_Syx = ind_unique_Syx.size();
        r_prop = num_unique_Syx + 1;
        k_prop_1 = k_prop - 1;

        // Calculate K_{x,2+}^{*}
        max_Kx_Sy_temp = zeros(Ky+1);
        max_Kx_Sy_temp.head(Ky) = max_Kx_Sy;
        max_Kx_Sy_temp(k_curr_1) = max_Kx_Sy_temp(k_curr_1) - 1;
        max_Kx_Sy_temp(k_prop_1) = max_Kx_Sy_temp(k_prop_1) + 1;
        ind_max_Kx_Sy_temp = find(max_Kx_Sy_temp == 1);
        num_max_Kx_Sy_temp = ind_max_Kx_Sy_temp.size();

        Kx2_prop = Kyx - num_max_Kx_Sy_temp;

        // Calculate the acceptance probability
        n_k_prop = n_k(k_prop_1);

        n_k_curr = n_k(k_curr_1);
        n_rk_curr = n_rk(k_curr_1,r_curr_1);

        ind_Syx = find(Sy == k_curr && Sx == r_curr);
        num_Syx = ind_Syx.size();

        vec Y_temp = Y(ind_Syx);
        mat matX_temp = matX.rows(ind_Syx);

        vec beta_prop = betaPars.col(k_prop_1);
        vec sig_prop = sqrt(sig2Pars(k_prop_1)) * ones(num_Syx);

        ratio =
          (lgamma(n_k_prop + n_rk_curr) - (lgamma(n_rk_curr) + lgamma(n_k_prop))) +
          ((lgamma(alpha_omega + n_rk_curr) + lgamma(alpha_omega + n_k_prop)) -
          (lgamma(alpha_omega) + lgamma(alpha_omega + n_k_prop + n_rk_curr))) +
          as_scalar((sum(log_normpdf(Y_temp, matX_temp * beta_prop, sig_prop))) -
          (sum(log(f0_y(ind_Syx))))) -
          log((Kx1_curr * Kyx - 1)/Kx2_prop) - log(alpha_theta);
        ratio = exp(ratio);

        // If ratio > runif (1), retrun "proposed" and o.w.,do nothing, i.e. retrun "current".
        if (ratio > randu(1)[0]){
          // Delete theta parameters
          betaPars.shed_col(k_curr_1);
          sig2Pars.shed_col(k_curr_1);

          // reorder omega parameters
          uvec ind_unique_Syx_curr = find(unique_Syx.col(0)==k_curr && unique_Syx.col(1)==r_curr);
          uvec ind_Sy0_prop = find(unique_Syx.col(0) <= k_prop);
          int num_Sy0_prop = ind_Sy0_prop.size();

          etaPars  = reorder(etaPars, ind_unique_Syx_curr(0), num_Sy0_prop);
          piPars   = reorder(piPars,  ind_unique_Syx_curr(0), num_Sy0_prop);
          muPars   = reorder(muPars,  ind_unique_Syx_curr(0), num_Sy0_prop);
          tau2Pars = reorder(tau2Pars, ind_unique_Syx_curr(0), num_Sy0_prop);

          // relabel clusters
          ind_Syx = find(Sy == k_curr && Sx == r_curr);
          num_Syx = ind_Syx.size();
          for (int ii = 0; ii < num_Syx; ii++) {
            Sy(ind_Syx(ii)) = k_prop;
            Sx(ind_Syx(ii)) = r_prop;
          }

          // relabel Y cluster
          ind_Sy = find(Sy >= k_curr);
          num_Sy = ind_Sy.size();
          for (int ii = 0; ii < num_Sy; ii++) {
            Sy(ind_Sy(ii)) -= 1;
          }

          // wrap-up
          Syx = join_rows(Sy, Sx);

          // Determine unique clusters from the cluster membership variable, Syx
          unique_Syx = unique_rows(Syx);

          // Make vector of y-clusters
          unique_Sy = unique(Sy);

          // Calculate the number of clusters
          Kyx = unique_Syx.n_rows;

          // Calculate the number of y-clusters
          Ky = unique_Sy.size();

          // Find the largest number of x clusters
          max_Sx = Sx.max();

          // Calculate the number of subjects in each y-cluster and store in vector n_k.
          // Calculate the number of subjects in each x-cluster and store in matrix n_rk.
          // Use k to index y-clusters and r to index x-clusters.
          n_k.zeros(Ky);
          n_rk.zeros(Ky, max_Sx);
          max_Kx_Sy.zeros(Ky);
          for (int k = 0; k < Ky; k++) {
            ind_Sy = find(Sy == (k+1));
            n_k(k) = ind_Sy.size();

            ind_Sy0 = find(unique_Syx.col(0) == (k+1));
            num_Sy0 = ind_Sy0.size();
            max_Kx_Sy(k) = num_Sy0;
            for (int r = 0; r < num_Sy0; r++) {
              ind_Syx = find(Sx(ind_Sy) == (r+1));
              n_rk(k,r) = ind_Syx.size();
            }
          }
        }
      }
    }
  }

  return List::create(
    _["Syx"] = Syx,
    _["unique_Sy"] = unique_Sy,
    _["unique_Syx"] = unique_Syx,
    _["betaPars"] = betaPars,
    _["sig2Pars"] = sig2Pars,
    _["etaPars"] = etaPars,
    _["piPars"] = piPars,
    _["muPars"] = muPars,
    _["tau2Pars"] = tau2Pars,
    _["Ky"] = Ky,
    _["Kyx"] = Kyx,
    _["n_k"] = n_k,
    _["n_rk"] = n_rk,
    _["max_Sx"] = max_Sx,
    _["max_Kx_Sy"] = max_Kx_Sy);
}
