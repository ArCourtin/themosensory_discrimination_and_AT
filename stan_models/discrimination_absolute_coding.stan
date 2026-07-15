//This Stan program implements a hierarchical version of the "no habituation - absolute coding" model of thermosensory discrimination
//Response-coded: the outcome is the second-interval choice of a 2IFC task, not accuracy.
//Evidence is the difference between the target's and the adapting temperature's absolute (baseline-anchored)
//representations; there is no separate detection-threshold parameter (the threshold is carried by beta).
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

    // Absolute coding: each interval's warm-channel activation is the soft-rectified distance from
    // the absolute reference (baseline + rho), i.e. no recentering on the adapting temperature. The
    // discrimination evidence is the difference between the target's and the adapting temperature's
    // absolute representations, so the two accounts coincide for warm-adapting conditions and diverge
    // only when the adapting temperature falls below the absolute reference (cold-adapting).
    vector[N] centered_target = centered_absolute_target_temperature - rho[participant];
    vector[N] centered_adapting = centered_absolute_adapting_temperature - rho[participant];
    vector[N] target_representation = centered_target .* inv_logit(100*centered_target);
    vector[N] adapting_representation = centered_adapting .* inv_logit(100*centered_adapting);
    vector[N] stimulus_representation = target_representation - adapting_representation;

    theta = lambda[participant] + (1-2*lambda[participant]) .* Phi(interval_sign .* (beta[participant] .* stimulus_representation) - kappa[participant]);
  }
}
model{
  //Priors
  // rho's prior is centered above the warmest tested adapting temperature (baseline+2) rather than
  // at 0: for rho <= (adapting - baseline), the absolute model collapses to being behaviorally
  // identical to relative coding (see stimulus_representation above), so a prior centered inside the
  // tested AT range would place most of the absolute account's prior mass on parameter values where
  // it isn't actually distinguishable from the model it's meant to compete against.
  mu[1] ~ normal(2,2);
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
