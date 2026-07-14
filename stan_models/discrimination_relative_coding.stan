//This Stan program implements a hierarchical version of the "complete habituation - relative coding" model of thermosensory discrimination
//Response-coded: the outcome is the second-interval choice of a 2IFC task, not accuracy.
//Evidence is the soft-rectified deviation of the target from the adapting temperature; there is no
//separate detection-threshold parameter (the threshold is carried by beta).
//Licence: MIT
//Author: Arthur S. Courtin

data{
  int N;
  int P;
  int is_cold;

  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  vector[N] interval_sign;                       // +1 if the deviating stimulus was in the second interval, -1 if in the first
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

  int M=3;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  vector[P] beta;
  vector[P] lambda;
  vector[P] kappa;
  vector[N] theta;

  {
    matrix[P,M] delta_participant = (diag_pre_multiply(tau, L) * z)';

    beta = exp(mu[1] + delta_participant[,1]);
    lambda = .5 * inv_logit(mu[2] + delta_participant[,2]);
    kappa = mu[3] + delta_participant[,3];

    // Relative coding: the representation is referenced to the adapting temperature (complete
    // habituation), so the discrimination evidence is the soft-rectified deviation of the target
    // from the adapting temperature.
    vector[N] stimulus_representation = deviation_from_adapting_temperature .* inv_logit(100*deviation_from_adapting_temperature);

    theta = lambda[participant] + (1-2*lambda[participant]) .* Phi(interval_sign .* (beta[participant] .* stimulus_representation) - kappa[participant]);
  }
}
model{
  //Priors
  mu[1] ~ normal(0,1);
  mu[2] ~ normal(-4,1);
  mu[3] ~ normal(0,1);

  tau[1] ~ normal(0,1);
  tau[2] ~ normal(0,1);
  tau[3] ~ normal(0,1);

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
