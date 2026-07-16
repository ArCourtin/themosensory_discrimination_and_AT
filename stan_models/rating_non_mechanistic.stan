//This Stan program implements a hierarchical version of the non-mechanistic model of thermosensory magnitude estimation
//The condition-specific intercept/slope profiles across the five adapting-temperature conditions are built as a
//chain of independent steps outward from the AT=baseline condition (condition 3): each step's group mean and its
//participant-level deviation are ordinary entries in the same mu/tau hierarchy used for the zero/one-inflation
//bounds and eta, so adjacent AT conditions share more of the same cumulative terms (and so are more likely to be
//similar) than distant ones, without any shared/estimated smoothness parameter.
//Licence: MIT
//Author: Arthur S. Courtin
//Edited with the assistance of Claude Code (Anthropic).

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
  vector[N] deviation_from_adapting_temperature;
  if(is_cold==1){
    deviation_from_adapting_temperature = absolute_adapting_temperature - absolute_target_temperature;
  }else{
    deviation_from_adapting_temperature = absolute_target_temperature - absolute_adapting_temperature;
  }

  int M=13;                                       // participant-hierarchy dimensions: 5 intercept, 5 slope, lower, upper, eta
  int C=5;                                         // adapting-temperature conditions
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
  row_vector[P] lower_bound;
  row_vector[P] upper_bound;
  row_vector[P] eta;

  vector[N] latent_representation;

  {
    matrix[M,P] delta_participant = diag_pre_multiply(tau, L) * z;

    // Intercept profile: mu[1] is the level at condition 3 (adapting temperature = baseline); mu[2]/mu[3]
    // are the independent outward steps toward conditions 4 and 5, mu[4]/mu[5] the outward steps toward
    // conditions 2 and 1. Each step's participant-level deviation chains the same way.
    row_vector[P] i3 = mu[1] + delta_participant[1];
    row_vector[P] i4 = i3 + mu[2] + delta_participant[2];
    row_vector[P] i5 = i4 + mu[3] + delta_participant[3];
    row_vector[P] i2 = i3 + mu[4] + delta_participant[4];
    row_vector[P] i1 = i2 + mu[5] + delta_participant[5];
    intercept[1] = i1;
    intercept[2] = i2;
    intercept[3] = i3;
    intercept[4] = i4;
    intercept[5] = i5;

    // Slope profile: same chained construction, mu[6] the condition-3 (log-scale) level, mu[7]/mu[8] the
    // outward steps toward conditions 4 and 5, mu[9]/mu[10] the outward steps toward conditions 2 and 1.
    row_vector[P] s3 = mu[6] + delta_participant[6];
    row_vector[P] s4 = s3 + mu[7] + delta_participant[7];
    row_vector[P] s5 = s4 + mu[8] + delta_participant[8];
    row_vector[P] s2 = s3 + mu[9] + delta_participant[9];
    row_vector[P] s1 = s2 + mu[10] + delta_participant[10];
    slope[1] = exp(s1);
    slope[2] = exp(s2);
    slope[3] = exp(s3);
    slope[4] = exp(s4);
    slope[5] = exp(s5);

    lower_bound = exp(mu[11] + delta_participant[11]);
    upper_bound = exp(mu[12] + delta_participant[12]);
    eta = exp(mu[13] + delta_participant[13]);

    for(n in 1:N){
      latent_representation[n] = intercept[adapting_temperature_idx[n],participant[n]] + slope[adapting_temperature_idx[n],participant[n]] .* deviation_from_adapting_temperature[n];
    }
  }
  vector[N] mu_rating = inv_logit(latent_representation);
}
model{
  //Intercept profile: mu[1] is the condition-3 (baseline) level, mu[2:5] are independent chained steps
  mu[1] ~ normal(-2,1);
  mu[2:5] ~ normal(0,0.5);

  //Slope profile: mu[6] is the condition-3 (baseline) level, mu[7:10] are independent chained steps
  mu[6] ~ normal(-2,1);
  mu[7:10] ~ normal(0,0.5);

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
