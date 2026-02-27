//This Stan program implements a hierarchical version of the "partial habituation - mixed coding" model of thermosensory magnitude estimation
//Licence: MIT
//Author: Arthur S. Courtin

data{
  int N;
  int P;
  int is_cold;

  vector[N] recorded_baseline_temperature;
  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  vector<lower=0,upper=1>[N] rating;

  array[N] int<lower=0,upper=1> choice_accuracy;
  array[N] int<lower=1,upper=P> participant;
}
transformed data{
  vector[N] deviation_from_adapting_temperature;
  vector[N] deviation_from_adapting_temperature_centered;
  vector[N] relative_adapting_temperature;
  
  if(is_cold==1){
    deviation_from_adapting_temperature = (absolute_adapting_temperature) - absolute_target_temperature;
    deviation_from_adapting_temperature_centered = (absolute_adapting_temperature-4) - absolute_target_temperature;
    relative_adapting_temperature = recorded_baseline_temperature - absolute_adapting_temperature;
  }else{
    deviation_from_adapting_temperature = absolute_target_temperature - (absolute_adapting_temperature);
    deviation_from_adapting_temperature_centered = absolute_target_temperature - (absolute_adapting_temperature+8);
    relative_adapting_temperature = absolute_adapting_temperature - recorded_baseline_temperature;
  }
  
  
  int M=3+5;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  vector[P] alpha;
  vector[P] sigma;
  vector[P] lambda;
  vector[N] beta;
  vector[N] theta;
  
  vector[P] intercept;
  vector[P] slope;
  vector[P] lower_bound;
  vector[P] upper_bound;
  vector[P] eta;
  
  vector[N] latent_representation;
  
  {
    matrix[P,M] delta_participant = (diag_pre_multiply(tau, L) * z)';

    alpha = exp(mu[1] + delta_participant[,1]);
    sigma = exp(mu[2] + delta_participant[,2]);
    lambda = .5 * inv_logit(mu[3] + delta_participant[,3]);
    
    intercept = mu[4] + delta_participant[,4];
    slope = exp(mu[5] + delta_participant[,5]);
    lower_bound = exp(mu[6] + delta_participant[,6]);
    upper_bound = exp(mu[7] + delta_participant[,7]);
    eta = exp(mu[8] + delta_participant[,8]);
    
    vector[N] centered_stimulus = deviation_from_adapting_temperature - alpha[participant];
    beta = 3/(2+sigma[participant]+alpha[participant]-relative_adapting_temperature);
    vector[N] weight = inv_logit(centered_stimulus*100);
    
    theta = (1-weight) .* 0.5 + weight .* (1-lambda[participant]) .* Phi(beta .* centered_stimulus);  
    
    latent_representation = intercept[participant] + slope[participant] .* beta .* deviation_from_adapting_temperature_centered;
  }
  vector[N] mu_rating = inv_logit(latent_representation);
}
model{
  //Priors
  mu[1] ~ normal(-2,1);
  mu[2] ~ normal(1,1);
  mu[3] ~ normal(-4,1);
  mu[4:5] ~ normal(-2,1);
  mu[6:7] ~ normal(2,.5);
  mu[8] ~ normal(3,1);

  tau ~ normal(0,1);
  tau[6:7] ~ normal(0,.5);
  
  L ~ lkj_corr_cholesky(1);

  to_vector(z) ~ std_normal();
  
  //Likelihood
  choice_accuracy ~ bernoulli(theta);
  
  for(idx in 1:N){
    if (rating[idx]==0){
      target += log1m_inv_logit(latent_representation[idx] + lower_bound[participant[idx]]);
    }else if(rating[idx]==1){
      target += log_inv_logit(latent_representation[idx] - upper_bound[participant[idx]]);
    }else{
      target += log(inv_logit(latent_representation[idx] + lower_bound[participant[idx]]) - inv_logit(latent_representation[idx] - upper_bound[participant[idx]]));

      rating[idx] ~ beta_proportion(mu_rating[idx], eta[participant[idx]]);
    }
  }
}
generated quantities{
  corr_matrix[M] cor = multiply_lower_tri_self_transpose(L);
  vector[N] log_lik;
  
  for (idx in 1:N) {
    if (rating[idx]==0){
      log_lik[idx] = log1m_inv_logit(latent_representation[idx] + lower_bound[participant[idx]]);
    }else if(rating[idx]==1){
      log_lik[idx] = log_inv_logit(latent_representation[idx] - upper_bound[participant[idx]]);
    }else{
      log_lik[idx] = log(inv_logit(latent_representation[idx] + lower_bound[participant[idx]]) - inv_logit(latent_representation[idx] - upper_bound[participant[idx]])) 
        + beta_proportion_lpdf(rating[idx] | mu_rating[idx], eta[participant[idx]]);
    }
  }
}
