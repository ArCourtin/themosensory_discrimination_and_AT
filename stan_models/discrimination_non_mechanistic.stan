//This Stan program implements a hierarchical version of the non-mechanistic model of thermosensory discrimination
//Licence: MIT
//Author: Arthur S. Courtin

data{
  int N;
  int P;
  int is_cold;

  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  array[N] int adapting_temperature_idx;
  array[N] int<lower=0,upper=1> choice_accuracy;
  array[N] int<lower=1,upper=P> participant;
}
transformed data{
  vector[N] deviation_from_adapting_temperature;
  if(is_cold==1){
    deviation_from_adapting_temperature = absolute_adapting_temperature - absolute_target_temperature;
  }else{
    deviation_from_adapting_temperature = absolute_target_temperature - absolute_adapting_temperature;
  }
  
  int M=11;
  int C=5;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  matrix[C,P] alpha;
  matrix[C,P] beta;
  row_vector[P] lambda;
  vector[N] theta;
  
  {
    matrix[M,P] delta_participant = diag_pre_multiply(tau, L) * z;
   
    for(idx in 1:5){
      alpha[idx] = exp(mu[idx] + delta_participant[idx]);
      beta[idx] = exp(mu[idx+5] + delta_participant[idx+5]);
    }
    lambda = .5 * inv_logit(mu[11] + delta_participant[11]);
    
    for(n in 1:N){
      real centered_stimulus = deviation_from_adapting_temperature[n] - alpha[adapting_temperature_idx[n],participant[n]];
      real stimulus_representation = centered_stimulus * inv_logit(100*centered_stimulus);
      
      theta[n] = lambda[participant[n]] + (1-2*lambda[participant[n]]) * Phi(beta[adapting_temperature_idx[n],participant[n]] * stimulus_representation);  
    }  
  }
}
model{
  //Priors
  mu[1:5] ~ normal(-2,1);
  mu[6:10] ~ normal(0,1);
  mu[11] ~ normal(-4,1);
  
  tau ~ normal(0,1);
  
  L ~ lkj_corr_cholesky(2);

  to_vector(z) ~ std_normal();
  
  //Likelihood
  choice_accuracy ~ bernoulli(theta);
}
generated quantities{
  corr_matrix[M] cor = multiply_lower_tri_self_transpose(L);
  vector[N] log_lik;
  
  for(n in 1:N){
    log_lik[n] = bernoulli_lpmf(choice_accuracy[n]|theta[n]);
  }
}
