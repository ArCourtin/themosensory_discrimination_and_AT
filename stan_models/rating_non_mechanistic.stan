//This Stan program implements a hierarchical version of the non-mechanistic model of thermosensory magnitude estimation
//The condition-specific intercept/slope of the latent linear function are free across the ordered adapting-temperature
//axis, but their group-level profiles are given a second-order random-walk (RW2) prior: the cheap default is linear in
//adapting temperature, with departures from linearity penalized by an estimated smoothness scale. This keeps the model
//strictly more flexible than the mechanistic accounts while trimming the wiggle-room no plausible DGP would use, so its
//role as a competitor in model comparison is a fair one.
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
  //Group-level RW2 profile of the latent intercept across the ordered adapting-temperature axis (identity scale)
  real i0;                                        // intercept: level at condition 1
  real i1;                                        // intercept: initial slope (per AT step)
  vector[C-2] i_curv;                             // intercept: standardized curvature increments
  real<lower=0> sigma_i;                          // intercept: curvature scale (smoothness)

  //Group-level RW2 profile of the latent slope across the ordered adapting-temperature axis (log scale)
  real s0;                                        // slope: level at condition 1
  real s1;                                        // slope: initial slope (per AT step)
  vector[C-2] s_curv;                             // slope: standardized curvature increments
  real<lower=0> sigma_s;                          // slope: curvature scale (smoothness)

  real mu_lower;                                  // group-level lower zero-inflation bound (pre-transform)
  real mu_upper;                                  // group-level upper one-inflation bound (pre-transform)
  real mu_eta;                                    // group-level precision (pre-transform)

  //Participant-level non-centered hierarchy (order: 5 intercept, 5 slope, lower, upper, eta)
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  vector[C] f_int;                                 // group-level intercept profile across AT
  vector[C] f_slope;                               // group-level log-slope profile across AT
  matrix[C,P] intercept;
  matrix[C,P] slope;
  row_vector[P] lower_bound;
  row_vector[P] upper_bound;
  row_vector[P] eta;

  vector[N] latent_representation;

  //RW2 group profiles: two free anchors (level + first slope), curvature (second differences) penalized by sigma
  f_int[1] = i0;
  f_int[2] = i0 + i1;
  for(c in 3:C){
    f_int[c] = 2*f_int[c-1] - f_int[c-2] + sigma_i * i_curv[c-2];
  }
  f_slope[1] = s0;
  f_slope[2] = s0 + s1;
  for(c in 3:C){
    f_slope[c] = 2*f_slope[c-1] - f_slope[c-2] + sigma_s * s_curv[c-2];
  }

  {
    matrix[M,P] delta_participant = diag_pre_multiply(tau, L) * z;

    for(idx in 1:5){
      intercept[idx] = f_int[idx] + delta_participant[idx];
      slope[idx] = exp(f_slope[idx] + delta_participant[idx+5]);
    }

    lower_bound = exp(mu_lower + delta_participant[11]);
    upper_bound = exp(mu_upper + delta_participant[12]);
    eta = exp(mu_eta + delta_participant[13]);

    for(n in 1:N){
      latent_representation[n] = intercept[adapting_temperature_idx[n],participant[n]] + slope[adapting_temperature_idx[n],participant[n]] .* deviation_from_adapting_temperature[n];
    }
  }
  vector[N] mu_rating = inv_logit(latent_representation);
}
model{
  //Priors on the RW2 intercept profile (i0 matches the original intercept-reference prior)
  i0 ~ normal(-2,1);
  i1 ~ normal(0,0.5);
  i_curv ~ std_normal();
  sigma_i ~ normal(0,0.3);

  //Priors on the RW2 slope profile (s0 matches the original slope-reference prior)
  s0 ~ normal(-2,1);
  s1 ~ normal(0,0.5);
  s_curv ~ std_normal();
  sigma_s ~ normal(0,0.3);

  mu_lower ~ normal(2,.5);
  mu_upper ~ normal(2,.5);
  mu_eta ~ normal(3,1);

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
