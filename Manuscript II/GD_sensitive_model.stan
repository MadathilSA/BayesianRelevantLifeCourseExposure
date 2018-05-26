data{
  int N;                        //Sample size
  int t;                        //Time points
  matrix[N,t] X;                //Exposure matrix
  int<lower=0, upper=1> y[N];   //Binary outcome
  
  real True_delta;         // True value of lifetime effect
  vector[t] StC[t];             //Reference Vector of weights for Critical periods
  vector[t] StA;                //Reference vector of weights for Accumulation period
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
  vector[t] bw;
}
transformed parameters{
  vector[N] xb;
  xb = X*bw;
}
model{
  for(i in 1:t){
    bw[i] ~ cauchy(0,2.5);
  }
  y ~ bernoulli_logit(xb);
}
generated quantities{
  vector[N] log_lik;              //For WAIC and LOO
  vector[t] bias_bw;              // Bias in each period effect
  vector[t+2] EucDist;            //Euclidean distance from reference vectors

  for(i in 1:t){
    bias_bw[i] = StS_bw[i]-bw[i];
  }
  
  for(i in 1:t){
    EucDist[i] = distance(bw, StC_bw[i]);
  }

  EucDist[t+1] = distance(bw, StA_bw);
  EucDist[t+2] = distance(bw, StS_bw);

  
  for(n in 1:N){
    log_lik[n] = bernoulli_logit_lpmf(y[n]| xb[n]);
  }
}
