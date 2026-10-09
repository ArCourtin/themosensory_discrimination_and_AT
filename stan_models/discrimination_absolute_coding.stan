//This Stan program implements a hierarchical version of the "no habituation - absolute coding" model of thermosensory discrimination
//Response-coded: the outcome is the second-interval choice of a 2IFC task, not accuracy.
//Evidence is zero at and below the adapting temperature (adaptation masks the absolute reading up to that
//point) and equal to the raw absolute (baseline+xi anchored) reading everywhere above it; the positive hard threshold xi defines the silent range
//and beta is the precision. There is no explicit comparison of
//two representations - the adapting temperature only sets where the mask releases.
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

  vector[N] recorded_baseline_temperature;
  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  vector[N] interval_sign;                       // +1 if the deviating stimulus was in the second interval, -1 if in the first
  array[N] int<lower=1,upper=5> adapting_temperature_idx;
  array[N] int<lower=0,upper=1> chose_second;    // 1 if the participant chose the second interval
  array[N] int<lower=1,upper=P> participant;
}
transformed data{
  vector[N] centered_absolute_target_temperature;
  vector[N] centered_absolute_adapting_temperature;
  if(is_cold==1){
    centered_absolute_target_temperature = recorded_baseline_temperature - absolute_target_temperature;
    centered_absolute_adapting_temperature = recorded_baseline_temperature - absolute_adapting_temperature;
  }else{
    centered_absolute_target_temperature = absolute_target_temperature - recorded_baseline_temperature;
    centered_absolute_adapting_temperature = absolute_adapting_temperature - recorded_baseline_temperature;
  }
  int M=12;                                       // participant-hierarchy dimensions: xi, beta, 5 lambda, 5 kappa
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
  vector[P] xi;
  vector[P] beta;
  matrix[C,P] lambda;
  matrix[C,P] kappa;
  vector[N] theta;

  {
    matrix[M,P] delta_participant = diag_pre_multiply(tau, L) * z;

    xi = exp(mu[1] + delta_participant[1]');
    beta = exp(mu[2] + delta_participant[2]');
    lambda = .5 * inv_logit(chain_profile(mu[3:7], delta_participant[3:7]));
    kappa = chain_profile(mu[8:12], delta_participant[8:12]);

    // Absolute coding: the target's warm-channel activation is the soft-rectified distance from the
    // absolute reference (baseline + xi) - no recentering on the adapting temperature. Adaptation does
    // not shift this reading; it masks it up to the adapting temperature itself (no free window size),
    // releasing sharply above it. So evidence is ~0 for target readings at or below the adapting
    // temperature, and ~(target - baseline - xi) above it. When the adapting temperature itself falls
    // below xi, the mask never engages and the model reduces to pure unmasked absolute coding.
    vector[N] centered_target = centered_absolute_target_temperature - xi[participant];
    vector[N] absolute_reading = centered_target .* inv_logit(100*centered_target);
    vector[N] mask_gate = inv_logit(100*(absolute_reading - (centered_absolute_adapting_temperature - xi[participant])));
    vector[N] stimulus_representation = absolute_reading .* mask_gate;

    for(n in 1:N){
      int c = adapting_temperature_idx[n];
      int p = participant[n];
      theta[n] = lambda[c,p] + (1-2*lambda[c,p]) * Phi(interval_sign[n] * beta[p] * stimulus_representation[n] + kappa[c,p]);
    }
  }
}
model{
  //Priors
  // xi is the hard threshold defining the silent range: the absolute floor below which nothing reads as warm, positive so the floor
  // sits above baseline, where the adapted sensation is neutral. The prior keeps it within the tested range.
  mu[1] ~ normal(-2,xi_sd);
  mu[2] ~ normal(beta_mean,beta_sd);
  // lambda and kappa: mu[3]/mu[8] are the condition-3 (baseline) levels, the rest independent chained steps
  mu[3] ~ normal(lambda_mean,lambda_sd);
  mu[4:7] ~ normal(0,lambda_sd/2);
  mu[8] ~ normal(0,kappa_sd);
  mu[9:12] ~ normal(0,kappa_sd/2);

  tau[1] ~ normal(0,xi_sd);
  tau[2] ~ normal(0,beta_sd);
  tau[3:7] ~ normal(0,lambda_sd);
  tau[8:12] ~ normal(0,kappa_sd);

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
