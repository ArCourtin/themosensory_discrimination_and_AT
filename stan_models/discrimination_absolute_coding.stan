//This Stan program implements a hierarchical version of the "no habituation - absolute coding" model of thermosensory discrimination
//Response-coded: the outcome is the second-interval choice of a 2IFC task, not accuracy.
//Evidence is zero at and below the adapting temperature (adaptation masks the absolute reading up to that
//point) and equal to the raw absolute (baseline+rho anchored) reading everywhere above it; there is no
//separate detection-threshold parameter (the threshold is carried by beta) and no explicit comparison of
//two representations - the adapting temperature only sets where the mask releases.
//Licence: MIT
//Author: Arthur S. Courtin
//Edited with the assistance of Claude Code (Anthropic).

data{
  int N;
  int P;
  int is_cold;

  vector[N] recorded_baseline_temperature;
  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  vector[N] interval_sign;                       // +1 if the deviating stimulus was in the second interval, -1 if in the first
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
  int M=4;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  vector[P] rho;
  vector[P] beta;
  vector[P] lambda;
  vector[P] kappa;
  vector[N] theta;

  {
    matrix[P,M] delta_participant = (diag_pre_multiply(tau, L) * z)';

    rho = mu[1] + delta_participant[,1];
    beta = exp(mu[2] + delta_participant[,2]);
    lambda = .5 * inv_logit(mu[3] + delta_participant[,3]);
    kappa = mu[4] + delta_participant[,4];

    // Absolute coding: the target's warm-channel activation is the soft-rectified distance from the
    // absolute reference (baseline + rho) - no recentering on the adapting temperature. Adaptation does
    // not shift this reading; it masks it up to the adapting temperature itself (no free window size),
    // releasing sharply above it. So evidence is ~0 for target readings at or below the adapting
    // temperature, and ~(target - baseline - rho) above it. When the adapting temperature itself falls
    // below rho, the mask never engages and the model reduces to pure unmasked absolute coding.
    vector[N] centered_target = centered_absolute_target_temperature - rho[participant];
    vector[N] absolute_reading = centered_target .* inv_logit(100*centered_target);
    vector[N] mask_gate = inv_logit(100*(absolute_reading - (centered_absolute_adapting_temperature - rho[participant])));
    vector[N] stimulus_representation = absolute_reading .* mask_gate;

    theta = lambda[participant] + (1-2*lambda[participant]) .* Phi(interval_sign .* (beta[participant] .* stimulus_representation) - kappa[participant]);
  }
}
model{
  //Priors
  // rho is the absolute floor below which nothing reads as warm, centered at baseline (0): pushing it
  // well above the tested range would place the floor above every target temperature, flattening the
  // model to chance performance on every trial.
  mu[1] ~ normal(0,2);
  mu[2] ~ normal(0,1);
  mu[3] ~ normal(-4,1);
  mu[4] ~ normal(0,1);

  tau[1] ~ normal(0,2);
  tau[2] ~ normal(0,1);
  tau[3] ~ normal(0,1);
  tau[4] ~ normal(0,1);

  L ~ lkj_corr_cholesky(2);

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
