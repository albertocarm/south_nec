// =============================================================================
// Bayesian Promotion Time Cure Model (Bounded Cumulative Hazard) v2
// =============================================================================
// Soporta priors POR VARIABLE: prior_sd_beta y prior_sd_gamma son vectores
// =============================================================================

data {
  int<lower=1> N;
  int<lower=0> K_inc;
  int<lower=0> K_lat;
  matrix[N, K_inc] X;
  matrix[N, K_lat] Z;
  vector<lower=0>[N] t;
  array[N] int<lower=0, upper=1> d;
  
  // Per-variable prior SDs (vectors, one per covariate)
  vector<lower=0>[K_inc] prior_sd_beta;
  vector<lower=0>[K_lat] prior_sd_gamma;
  real<lower=0> prior_sd_intercepts;
}

parameters {
  real beta0;
  vector[K_inc] beta;
  real gamma0;
  vector[K_lat] gamma;
  real<lower=0> alpha;
}

transformed parameters {
  vector[N] log_theta;
  vector[N] log_lambda;
  
  log_theta  = beta0  + X * beta;
  log_lambda = gamma0 + Z * gamma;
}

model {
  // Intercepts
  beta0  ~ normal(0, prior_sd_intercepts);
  gamma0 ~ normal(0, prior_sd_intercepts);
  
  // Per-variable priors
  beta  ~ normal(0, prior_sd_beta);
  gamma ~ normal(0, prior_sd_gamma);
  
  // Weibull shape
  alpha ~ gamma(2, 1);
  
  // Log-likelihood
  for (i in 1:N) {
    real theta_i  = exp(log_theta[i]);
    real lambda_i = exp(log_lambda[i]);
    real neg_cum_haz_0 = -pow(t[i] / lambda_i, alpha);
    real S0_t = exp(neg_cum_haz_0);
    real log_Spop = -theta_i * (1 - S0_t);
    
    if (d[i] == 1) {
      real log_f0 = log(alpha) - log(lambda_i) 
                    + (alpha - 1) * log(t[i] / lambda_i) 
                    + neg_cum_haz_0;
      target += log_theta[i] + log_f0 + log_Spop;
    } else {
      target += log_Spop;
    }
  }
}

generated quantities {
  vector[N] cure_frac;
  vector[N] log_lik;
  
  for (i in 1:N) {
    real theta_i  = exp(log_theta[i]);
    real lambda_i = exp(log_lambda[i]);
    real neg_cum_haz_0 = -pow(t[i] / lambda_i, alpha);
    real S0_t = exp(neg_cum_haz_0);
    real log_Spop = -theta_i * (1 - S0_t);
    
    cure_frac[i] = exp(-theta_i);
    
    if (d[i] == 1) {
      real log_f0 = log(alpha) - log(lambda_i) 
                    + (alpha - 1) * log(t[i] / lambda_i) 
                    + neg_cum_haz_0;
      log_lik[i] = log_theta[i] + log_f0 + log_Spop;
    } else {
      log_lik[i] = log_Spop;
    }
  }
}
