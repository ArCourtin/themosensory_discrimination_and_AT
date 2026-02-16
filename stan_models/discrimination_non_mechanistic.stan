//This Stan program implements a hierarchical version of the "complete habituation - relative coding" model of thermosensory discrimination
//Licence: MIT
//Author: Arthur S. Courtin

data{
  int N;
  int P;

  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  array[N] int adapting_temperature_idx;
  array[N] int<lower=0,upper=1> choice_accuracy;
  array[N] int<lower=1,upper=P> participant;
}
transformed data{
  vector[N] deviation_from_adapting_temperature = absolute_target_temperature - absolute_adapting_temperature;
  int M=3;
  int C=5;
}
parameters{
  vector[1+C*2] mu;
  vector<lower=0>[1+C*2] tau;
  matrix[1+C*2,P] z;
  cholesky_factor_corr[1+C*2] L;
}
transformed parameters{
  matrix[C,P] alpha;
  matrix[C,P] beta;
  row_vector[P] lambda;
  vector[N] theta;
  
  {
    matrix[1+C*2,P] delta_participant = diag_pre_multiply(tau, L) * z;

    alpha[1,] = exp(mu[1] + delta_participant[1,]+mu[2] + delta_participant[2,]);
    alpha[2,] = exp(mu[1] + delta_participant[1,]+mu[3] + delta_participant[3,]);
    alpha[3,] = exp(mu[1] + delta_participant[1,]);
    alpha[4,] = exp(mu[1] + delta_participant[1,]+mu[4] + delta_participant[4,]);
    alpha[5,] = exp(mu[1] + delta_participant[1,]+mu[5] + delta_participant[5,]);
    
    beta[1,] = exp(mu[6] + delta_participant[6,]+mu[7] + delta_participant[7,]);
    beta[2,] = exp(mu[6] + delta_participant[6,]+mu[8] + delta_participant[8,]);
    beta[3,] = exp(mu[6] + delta_participant[6,]);
    beta[4,] = exp(mu[6] + delta_participant[6,]+mu[9] + delta_participant[9,]);
    beta[5,] = exp(mu[6] + delta_participant[6,]+mu[10] + delta_participant[10,]);
    
    lambda = .5 * inv_logit(mu[11] + delta_participant[11,]);
    for(n in 1:N){
      real centered_stimulus = deviation_from_adapting_temperature[n] - alpha[adapting_temperature_idx[n],participant[n]];
      real weight = inv_logit(centered_stimulus*100);
      
      theta[n] = (1-weight) * 0.5 + (weight-lambda[participant[n]]) * Phi(beta[adapting_temperature_idx[n],participant[n]] * centered_stimulus);  
    }  
  }
}
model{
  //Priors
  mu ~ normal(0,1);
  mu[1] ~ normal(-2,1);
  mu[11] ~ normal(-4,1);
  
  tau ~ normal(0,1);
  
  L ~ lkj_corr_cholesky(1);

  to_vector(z) ~ std_normal();
  
  //Likelihood
  choice_accuracy ~ bernoulli(theta);
}
generated quantities{
  vector[N] log_lik;
  
  for(n in 1:N){
    log_lik[n] = bernoulli_lpmf(choice_accuracy[n]|theta[n]);
  }
}
