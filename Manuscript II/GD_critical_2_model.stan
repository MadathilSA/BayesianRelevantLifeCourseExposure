data{
  int N;
  int t;
  matrix[N,t] X;
  int<lower=0, upper=1> y[N];
  
  real True_delta;            // True value for lifetime effect
  vector[t] StC[t];             //Reference Vector of weights for Critical periods
  vector[t] StA;                //Reference vector of weights  for Accumulation period
  vector[t] StS;                //Reference vector of weights for Sensitive period
}
transformed data{
  vector[t] StC_bw[t];             //Reference Vector of parameters for Critical periods
  vector[t] StA_bw;                //Reference vector of parameters for Accumulation period
  vector[t] StS_bw;                //Reference vector of parameters for Sensitive period
  
  for(i in 1:t){
    StC_bw[i] = StC[i] * True_delta;
  } 
  StA_bw = StA * True_delta;
  StS_bw = StS * True_delta;
}
parameters{
  real delta;
}
transformed parameters{
  vector[N] xb;
  vector[t] bw;
  
  bw[2] = delta;
  bw[1] = 0;
  bw[3] = 0;

  xb = X*bw;
}
model{
  delta ~ cauchy(0,2.5);
  y ~ bernoulli_logit(xb);
}
generated quantities{
  vector[N] log_lik;              //For WAIC and LOO
  real bias_delta;                //Bias 
  vector[t+2] EucDist;            //Euclidean distance from reference vectors
  
  bias_delta = True_delta-delta;
  
  for(i in 1:t){
    EucDist[i] = distance(bw, StC_bw[i]);
  }

  EucDist[t+1] = distance(bw, StA_bw);
  EucDist[t+2] = distance(bw, StS_bw);

  for(n in 1:N){
    log_lik[n] = bernoulli_logit_lpmf(y[n]| xb[n]);
  }
}
