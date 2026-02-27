//This Stan program implements a hierarchical version of the non-mechanistic model of thermosensory magnitude estimation
//Licence: MIT
//Author: Arthur S. Courtin

data{
  int N;
  int P;
  int is_cold;

  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  vector<lower=0,upper=1>[N] rating;
  
  array[N] int adapting_temperature_idx;
  array[N] int<lower=1,upper=P> participant;
}
transformed data{
  vector[N] deviation_from_adapting_temperature_centered;
  if(is_cold==1){
    deviation_from_adapting_temperature_centered = (absolute_adapting_temperature-4) - absolute_target_temperature;
  }else{
    deviation_from_adapting_temperature_centered = absolute_target_temperature - (absolute_adapting_temperature+8);
  }
  int C=5;
  int M=3+C*2;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  matrix[C,P] intercept;
  matrix[C,P] slope;
  vector[P] lower_bound;
  vector[P] upper_bound;
  vector[P] eta;
  
  vector[N] latent_representation;
  
  {
    matrix[P,M] delta_participant = (diag_pre_multiply(tau, L) * z)';

    intercept[1,] = mu[1] + delta_participant[1,]+mu[2] + delta_participant[2,];
    intercept[2,] = mu[1] + delta_participant[1,]+mu[3] + delta_participant[3,];
    intercept[3,] = mu[1] + delta_participant[1,];
    intercept[4,] = mu[1] + delta_participant[1,]+mu[4] + delta_participant[4,];
    intercept[5,] = mu[1] + delta_participant[1,]+mu[5] + delta_participant[5,];
    
    slope[1,] = exp(mu[6] + delta_participant[6,]+mu[7] + delta_participant[7,]);
    slope[2,] = exp(mu[6] + delta_participant[6,]+mu[8] + delta_participant[8,]);
    slope[3,] = exp(mu[6] + delta_participant[6,]);
    slope[4,] = exp(mu[6] + delta_participant[6,]+mu[9] + delta_participant[9,]);
    slope[5,] = exp(mu[6] + delta_participant[6,]+mu[10] + delta_participant[10,]);
    
    lower_bound = exp(mu[11] + delta_participant[,11]);
    upper_bound = exp(mu[12] + delta_participant[,12]);
    eta = exp(mu[13] + delta_participant[,13]);
    
    for(n in 1:N){
      latent_representation[n] = intercept[adapting_temperature_idx[n],participant[n]] + slope[adapting_temperature_idx[n],participant[n]] .* deviation_from_adapting_temperature_centered[n];
    }
  }
  vector[N] mu_rating = inv_logit(latent_representation);
}
model{
  //Priors
  mu[1] ~ normal(-2,1);
  mu[2:5] ~ normal(0,1);
  mu[6] ~ normal(-2,1);
  mu[7:10] ~ normal(0,1);
  mu[11:12] ~ normal(2,.5);
  mu[13] ~ normal(3,1);

  tau ~ normal(0,1);
  tau[11:12] ~ normal(0,.5);

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
