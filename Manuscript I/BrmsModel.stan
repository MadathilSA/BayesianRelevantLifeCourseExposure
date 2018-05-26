// generated with brms 2.2.0
functions { 

  /* hurdle lognormal log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the lognormal distribution 
   *   sigma: sd parameter of the lognormal distribution
   *   hu: hurdle probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_lognormal_lpdf(real y, real mu, real sigma, real hu) { 
     if (y == 0) { 
       return bernoulli_lpmf(1 | hu); 
     } else { 
       return bernoulli_lpmf(0 | hu) +  
              lognormal_lpdf(y | mu, sigma); 
     } 
   }
  /* hurdle lognormal log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the lognormal distribution 
   *   sigma: sd parameter of the lognormal distribution
   *   hu: linear predictor for the hurdle part 
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real hurdle_lognormal_logit_lpdf(real y, real mu, real sigma, real hu) { 
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | hu); 
     } else { 
       return bernoulli_logit_lpmf(0 | hu) +  
              lognormal_lpdf(y | mu, sigma); 
     } 
   }
} 
data { 
  int<lower=1> N;  // total number of observations 
  vector[N] Y;  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // population-level design matrix 
  // data of smooth s(AgeY)
  int nb_1;  // number of bases 
  int knots_1[nb_1]; 
  matrix[N, knots_1[1]] Zs_1_1; 
  int<lower=1> K_hu;  // number of population-level effects 
  matrix[N, K_hu] X_hu;  // population-level design matrix 
  // data of smooth s(AgeY)
  int nb_hu_1;  // number of bases 
  int knots_hu_1[nb_hu_1]; 
  matrix[N, knots_hu_1[1]] Zs_hu_1_1; 
  // data for group-level effects of ID 1
  int<lower=1> J_1[N];
  int<lower=1> N_1;
  int<lower=1> M_1;
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> J_2[N];
  int<lower=1> N_2;
  int<lower=1> M_2;
  vector[N] Z_2_hu_1;
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc = K - 1; 
  matrix[N, K - 1] Xc;  // centered version of X 
  vector[K - 1] means_X;  // column means of X before centering 
  int Kc_hu = K_hu - 1; 
  matrix[N, K_hu - 1] Xc_hu;  // centered version of X_hu 
  vector[K_hu - 1] means_X_hu;  // column means of X_hu before centering 
  for (i in 2:K) { 
    means_X[i - 1] = mean(X[, i]); 
    Xc[, i - 1] = X[, i] - means_X[i - 1]; 
  } 
  for (i in 2:K_hu) { 
    means_X_hu[i - 1] = mean(X_hu[, i]); 
    Xc_hu[, i - 1] = X_hu[, i] - means_X_hu[i - 1]; 
  } 
} 
parameters { 
  vector[Kc] b;  // population-level effects 
  real temp_Intercept;  // temporary intercept 
  // parameters of smooth s(AgeY)
  vector[knots_1[1]] zs_1_1; 
  real<lower=0> sds_1_1; 
  real<lower=0> sigma;  // residual SD 
  vector[Kc_hu] b_hu;  // population-level effects 
  real temp_hu_Intercept;  // temporary intercept 
  // parameters of smooth s(AgeY)
  vector[knots_hu_1[1]] zs_hu_1_1; 
  real<lower=0> sds_hu_1_1; 
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // unscaled group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // unscaled group-level effects
} 
transformed parameters { 
  vector[knots_1[1]] s_1_1 = sds_1_1 * zs_1_1; 
  vector[knots_hu_1[1]] s_hu_1_1 = sds_hu_1_1 * zs_hu_1_1; 
  // group-level effects 
  vector[N_1] r_1_1 = sd_1[1] * (z_1[1]);
  // group-level effects 
  vector[N_2] r_2_hu_1 = sd_2[1] * (z_2[1]);
} 
model { 
  vector[N] mu = Xc * b + Zs_1_1 * s_1_1 + temp_Intercept; 
  vector[N] hu = Xc_hu * b_hu + Zs_hu_1_1 * s_hu_1_1 + temp_hu_Intercept; 
  for (n in 1:N) { 
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    hu[n] += r_2_hu_1[J_2[n]] * Z_2_hu_1[n];
  } 
  // priors including all constants 
  target += student_t_lpdf(temp_Intercept | 3, 0, 10); 
  target += normal_lpdf(zs_1_1 | 0, 1); 
  target += student_t_lpdf(sds_1_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += student_t_lpdf(sigma | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += logistic_lpdf(temp_hu_Intercept | 0, 1); 
  target += normal_lpdf(zs_hu_1_1 | 0, 1); 
  target += student_t_lpdf(sds_hu_1_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(z_1[1] | 0, 1);
  target += student_t_lpdf(sd_2 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(z_2[1] | 0, 1);
  // likelihood including all constants 
  if (!prior_only) { 
    for (n in 1:N) { 
      target += hurdle_lognormal_logit_lpdf(Y[n] | mu[n], sigma, hu[n]); 
    } 
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_Intercept = temp_Intercept - dot_product(means_X, b); 
  // actual population-level intercept 
  real b_hu_Intercept = temp_hu_Intercept - dot_product(means_X_hu, b_hu); 
} 