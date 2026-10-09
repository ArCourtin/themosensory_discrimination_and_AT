//This Stan program implements a hierarchical version of the non-mechanistic model of thermosensory discrimination
//The condition-specific xi (hard threshold defining the silent range) / beta (precision) profiles across the five adapting-temperature conditions are built as a chain
//of independent steps outward from the AT=baseline condition (condition 3): each step's group mean and its
//participant-level deviation are ordinary entries in the same mu/tau hierarchy used for lambda (lapse rate) and kappa (interval bias), so
//adjacent AT conditions share more of the same cumulative terms (and so are more likely to be similar) than
//distant ones, without any shared/estimated smoothness parameter.
//Response-coded: the outcome is the second-interval choice of a 2IFC task, not accuracy.
//Licence: MIT
//Author: Arthur S. Courtin
//Edited with the assistance of Claude Code (Anthropic).

data{
  int N;
  int P;
  int is_cold;

  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  vector[N] interval_sign;                       // +1 if the deviating stimulus was in the second interval, -1 if in the first
  array[N] int adapting_temperature_idx;
  array[N] int<lower=0,upper=1> chose_second;    // 1 if the participant chose the second interval
  array[N] int<lower=1,upper=P> participant;
}
transformed data{
  vector[N] deviation_from_adapting_temperature;
  if(is_cold==1){
    deviation_from_adapting_temperature = absolute_adapting_temperature - absolute_target_temperature;
  }else{
    deviation_from_adapting_temperature = absolute_target_temperature - absolute_adapting_temperature;
  }

  int M=12;                                       // participant-hierarchy dimensions: 5 xi, 5 beta, lambda, kappa
  int C=5;                                         // adapting-temperature conditions

  // beta and lambda priors come from earlier yes/no thermal detection fits (beta converted to the 2IFC scale by /sqrt(2), with that conversion's uncertainty added in quadrature).
  real beta_mean = is_cold==1 ? 0.552-30*0.0157-log(sqrt(2)) : -0.786-30*0.00239-log(sqrt(2));
  real beta_sd = is_cold==1 ? sqrt(1.1^2+log(sqrt(2))^2) : sqrt(0.865^2+log(sqrt(2))^2);
  real lambda_mean = -3.85+30*0.00171;
  real lambda_sd = 0.72;
  real kappa_sd = 0.5;
  real xi_sd = 0.5;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  matrix[C,P] xi;
  matrix[C,P] beta;
  row_vector[P] lambda;
  row_vector[P] kappa;
  vector[N] theta;

  {
    matrix[M,P] delta_participant = diag_pre_multiply(tau, L) * z;

    // Alpha profile: mu[1] is the level at condition 3 (adapting temperature = baseline); mu[2]/mu[3]
    // are the independent outward steps toward conditions 4 and 5, mu[4]/mu[5] the outward steps toward
    // conditions 2 and 1. Each step's participant-level deviation chains the same way.
    row_vector[P] f_x3 = mu[1] + delta_participant[1];
    row_vector[P] f_x4 = f_x3 + mu[2] + delta_participant[2];
    row_vector[P] f_x5 = f_x4 + mu[3] + delta_participant[3];
    row_vector[P] f_x2 = f_x3 + mu[4] + delta_participant[4];
    row_vector[P] f_x1 = f_x2 + mu[5] + delta_participant[5];
    xi[1] = exp(f_x1);
    xi[2] = exp(f_x2);
    xi[3] = exp(f_x3);
    xi[4] = exp(f_x4);
    xi[5] = exp(f_x5);

    // Beta profile: same chained construction, mu[6] the condition-3 level, mu[7]/mu[8] the outward
    // steps toward conditions 4 and 5, mu[9]/mu[10] the outward steps toward conditions 2 and 1.
    row_vector[P] f_b3 = mu[6] + delta_participant[6];
    row_vector[P] f_b4 = f_b3 + mu[7] + delta_participant[7];
    row_vector[P] f_b5 = f_b4 + mu[8] + delta_participant[8];
    row_vector[P] f_b2 = f_b3 + mu[9] + delta_participant[9];
    row_vector[P] f_b1 = f_b2 + mu[10] + delta_participant[10];
    beta[1] = exp(f_b1);
    beta[2] = exp(f_b2);
    beta[3] = exp(f_b3);
    beta[4] = exp(f_b4);
    beta[5] = exp(f_b5);

    lambda = .5 * inv_logit(mu[11] + delta_participant[11]);
    kappa = mu[12] + delta_participant[12];

    for(n in 1:N){
      real centered_stimulus = deviation_from_adapting_temperature[n] - xi[adapting_temperature_idx[n],participant[n]];
      real stimulus_representation = centered_stimulus * inv_logit(100*centered_stimulus);

      theta[n] = lambda[participant[n]] + (1-2*lambda[participant[n]]) * Phi(interval_sign[n] * beta[adapting_temperature_idx[n],participant[n]] * stimulus_representation + kappa[participant[n]]);
    }
  }
}
model{
  //Alpha profile: mu[1] is the condition-3 (baseline) level, mu[2:5] are independent chained steps
  mu[1] ~ normal(-2,xi_sd);
  mu[2:5] ~ normal(0,0.5);

  //Beta profile: mu[6] is the condition-3 (baseline) level, mu[7:10] are independent chained steps
  mu[6] ~ normal(beta_mean,beta_sd);
  mu[7:10] ~ normal(0,0.5);

  mu[11] ~ normal(lambda_mean,lambda_sd);
  mu[12] ~ normal(0,kappa_sd);

  tau[1:5] ~ normal(0,xi_sd);
  tau[6:10] ~ normal(0,beta_sd);
  tau[11] ~ normal(0,lambda_sd);
  tau[12] ~ normal(0,kappa_sd);

  L ~ lkj_corr_cholesky(1);

  to_vector(z) ~ std_normal();

  //Likelihood
  chose_second ~ bernoulli(theta);
}
generated quantities{
  corr_matrix[M] cor = multiply_lower_tri_self_transpose(L);
  vector[N] log_lik;

  for(n in 1:N){
    log_lik[n] = bernoulli_lpmf(chose_second[n]|theta[n]);
  }
}
