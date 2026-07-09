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
   
    vector[5] alpha_condition_effect;
    matrix[5,P] alpha_participant_condition_effect;
    vector[5] beta_condition_effect;
    matrix[5,P] beta_participant_condition_effect;
    
    alpha_condition_effect[1:4] = mu[2:5];
    alpha_condition_effect[5] = -sum(mu[2:5]);
    beta_condition_effect[1:4] = mu[7:10];
    beta_condition_effect[5] = -sum(mu[7:10]);
    
    alpha_participant_condition_effect[1:4] = delta_participant[2:5];
    beta_participant_condition_effect[1:4] = delta_participant[7:10];
    for(p in 1:P){
      alpha_participant_condition_effect[5,p]
        = -sum(delta_participant[2:5,p]);
    
      beta_participant_condition_effect[5,p]
        = -sum(delta_participant[7:10,p]);
    }    

    for(c in 1:5){
      alpha[c] = exp(mu[1]+delta_participant[1]+alpha_condition_effect[c]+alpha_participant_condition_effect[c]);
      beta[c] = exp(mu[6]+delta_participant[6]+beta_condition_effect[c]+beta_participant_condition_effect[c]);
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
  mu[1] ~ normal(-2,1);
  mu[2:10] ~ normal(0,1);
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
