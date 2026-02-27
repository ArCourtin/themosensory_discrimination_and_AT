//This Stan program implements a hierarchical version of the "no habituation - absolute coding" model of thermosensory magnitude estimation
//Licence: MIT
//Author: Arthur S. Courtin

data{
  int N;
  int P;
  int is_cold;

  vector[N] recorded_baseline_temperature;
  vector[N] absolute_target_temperature;
  vector<lower=0,upper=1>[N] rating;
  array[N] int<lower=1,upper=P> participant;
}
transformed data{
  vector[N] target_temperature_centered_at_max_fixed;
  if(is_cold==1){
    target_temperature_centered_at_max_fixed = (recorded_baseline_temperature-2-4) - absolute_target_temperature;
  }else{
    target_temperature_centered_at_max_fixed = absolute_target_temperature - (recorded_baseline_temperature+2+8);
  }
  int M=5;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  vector[P] intercept;
  vector[P] slope;
  vector[P] lower_bound;
  vector[P] upper_bound;
  vector[P] eta;
  
  vector[N] latent_representation;
  
  {
    matrix[P,M] delta_participant = (diag_pre_multiply(tau, L) * z)';

    intercept = mu[1] + delta_participant[,1];
    slope = exp(mu[2] + delta_participant[,2]);
    lower_bound = exp(mu[3] + delta_participant[,3]);
    upper_bound = exp(mu[4] + delta_participant[,4]);
    eta = exp(mu[5] + delta_participant[,5]);
    
    latent_representation = intercept[participant] + slope[participant] .* target_temperature_centered_at_max_fixed;
  }
  vector[N] mu_rating = inv_logit(latent_representation);
}
model{
  //Priors
  mu[1:2] ~ normal(-2,1);
  mu[3:4] ~ normal(2,.5);
  mu[5] ~ normal(3,1);

  tau[1:2] ~ normal(0,1);
  tau[3:4] ~ normal(0,.5);
  tau[5] ~ normal(0,1);

  L ~ lkj_corr_cholesky(1);

  to_vector(z) ~ std_normal();
  
  //Likelihood
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
