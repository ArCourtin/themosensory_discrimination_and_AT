//This Stan program implements a hierarchical version of the non-mechanistic model of thermosensory discrimination
//The condition-specific xi (hard threshold defining the silent range) / beta (precision) / lambda (lapse rate) / kappa (interval bias) profiles across the five adapting-temperature conditions are built as a chain
//of independent steps outward from the AT=baseline condition (condition 3): each step's group mean and its
//participant-level deviation are ordinary entries in the same mu/tau hierarchy, so
//adjacent AT conditions share more of the same cumulative terms (and so are more likely to be similar) than
//distant ones, without any shared/estimated smoothness parameter.
//Response-coded: the outcome is the second-interval choice of a 2IFC task, not accuracy.
//Licence: MIT
//Author: Arthur S. Courtin
//Edited with the assistance of Claude Code (Anthropic).

functions{
  // Chained profile over the five adapting-temperature conditions: m[1]/d[1] set the level at condition 3
  // (adapting temperature = baseline), m[2:3]/d[2:3] the outward steps toward conditions 4 and 5,
  // m[4:5]/d[4:5] the outward steps toward conditions 2 and 1. Row c of the result is condition c.
  matrix chain_profile(vector m, matrix d){
    matrix[5,cols(d)] f;
    f[3] = m[1] + d[1];
    f[4] = f[3] + m[2] + d[2];
    f[5] = f[4] + m[3] + d[3];
    f[2] = f[3] + m[4] + d[4];
    f[1] = f[2] + m[5] + d[5];
    return f;
  }
}
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

  int M=20;                                       // participant-hierarchy dimensions: 5 each of xi, beta, lambda, kappa
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
  matrix[C,P] lambda;
  matrix[C,P] kappa;
  vector[N] theta;

  {
    matrix[M,P] delta_participant = diag_pre_multiply(tau, L) * z;

    // xi, beta, lambda and kappa each get a chained profile (see chain_profile): mu[1:5] xi, mu[6:10] beta,
    // mu[11:15] lambda, mu[16:20] kappa.
    xi = exp(chain_profile(mu[1:5], delta_participant[1:5]));
    beta = exp(chain_profile(mu[6:10], delta_participant[6:10]));
    lambda = .5 * inv_logit(chain_profile(mu[11:15], delta_participant[11:15]));
    kappa = chain_profile(mu[16:20], delta_participant[16:20]);

    for(n in 1:N){
      real centered_stimulus = deviation_from_adapting_temperature[n] - xi[adapting_temperature_idx[n],participant[n]];
      real stimulus_representation = centered_stimulus * inv_logit(100*centered_stimulus);

      theta[n] = lambda[adapting_temperature_idx[n],participant[n]] + (1-2*lambda[adapting_temperature_idx[n],participant[n]]) * Phi(interval_sign[n] * beta[adapting_temperature_idx[n],participant[n]] * stimulus_representation + kappa[adapting_temperature_idx[n],participant[n]]);
    }
  }
}
model{
  //Alpha profile: mu[1] is the condition-3 (baseline) level, mu[2:5] are independent chained steps
  mu[1] ~ normal(-2,xi_sd);
  mu[2:5] ~ normal(0,xi_sd/2);

  //Beta profile: mu[6] is the condition-3 (baseline) level, mu[7:10] are independent chained steps
  mu[6] ~ normal(beta_mean,beta_sd);
  mu[7:10] ~ normal(0,beta_sd/2);

  //Lambda profile: mu[11] is the condition-3 (baseline) level, mu[12:15] are independent chained steps
  mu[11] ~ normal(lambda_mean,lambda_sd);
  mu[12:15] ~ normal(0,lambda_sd/2);

  //Kappa profile: mu[16] is the condition-3 (baseline) level, mu[17:20] are independent chained steps
  mu[16] ~ normal(0,kappa_sd);
  mu[17:20] ~ normal(0,kappa_sd/2);

  tau[1:5] ~ normal(0,xi_sd);
  tau[6:10] ~ normal(0,beta_sd);
  tau[11:15] ~ normal(0,lambda_sd);
  tau[16:20] ~ normal(0,kappa_sd);

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
