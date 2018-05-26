data{
  int ind;
  int t;
  int tcol;
  
  matrix[t,tcol] ts;
  
  vector[ind] age;
  vector[ind] male;
  vector[ind] alc;
  vector[ind] edu;
  int<lower=0, upper=1> hpv[ind];
  int<lower=0> pos[ind];
  
  matrix[ind,t] SmkI;
  int<lower=0, upper=1> Status[ind];
}
parameters{
  vector<lower=0>[tcol] gammaN;
  vector<lower=0>[tcol] gammaP;
  real alpha;
  real bage;
  real bmale;
  real balc;
  real bedu;
  real bhpv;
}
transformed parameters{
  matrix[t,2] deltas;
  vector[ind] xb1;
  vector[ind] xb;
  
  deltas[,1] = ts*gammaN;
  deltas[,2] = ts*gammaP;
  
  for(n in 1:ind){
    xb1[n] = SmkI[n,1:pos[n]]*to_vector(deltas[1:pos[n],(hpv[n]+1)]);
    
    xb[n]=alpha+(bage*age[n])+(bmale*male[n])+(balc*alc[n])+(bedu*edu[n])+(bhpv*hpv[n])+xb1[n];
  }
}
model{
  
  alpha ~ student_t(3,0,5);

  bage  ~ student_t(3,0,2.5);
  bmale ~ student_t(3,0,2.5);
  balc  ~ student_t(3,0,2.5);
  bedu  ~ student_t(3,0,2.5);
  bhpv  ~ student_t(3,0,2.5);

  for(i in 1:tcol){
    gammaN[i] ~ student_t(3,0,2.5);
    gammaP[i] ~ student_t(3,0,2.5);
  }
  
  Status ~ bernoulli_logit(xb);
}
generated quantities{
  vector[t] hpvP_or;
  vector[t] hpvN_or;
  
  hpvN_or = exp(deltas[,1]);
  hpvP_or = exp(deltas[,2]);
}
