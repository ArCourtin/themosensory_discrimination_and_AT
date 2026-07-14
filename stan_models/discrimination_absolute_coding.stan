//This Stan program implements a hierarchical version of the "no habituation - absolute coding" model of thermosensory discrimination
//Response-coded: the outcome is the second-interval choice of a 2IFC task, not accuracy.
//Licence: MIT
//Author: Arthur S. Courtin

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
  int M=5;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  vector[P] alpha;
  vector[P] rho;
  vector[P] beta;
  vector[P] lambda;
  vector[P] kappa;
  vector[N] theta;

  {
    matrix[P,M] delta_participant = (diag_pre_multiply(tau, L) * z)';

    rho = mu[1] + delta_participant[,1];
    alpha = exp(mu[2] + delta_participant[,2]);
    beta = exp(mu[3] + delta_participant[,3]);
    lambda = .5 * inv_logit(mu[4] + delta_participant[,4]);
    kappa = mu[5] + delta_participant[,5];

    vector[N] centered_stimulus = centered_absolute_target_temperature - rho[participant];
    vector[N] stimulus_representation = centered_stimulus .* inv_logit(100*centered_stimulus);
    vector[N] adapted_stimulus_representation = stimulus_representation .* inv_logit(100.*(stimulus_representation-(centered_absolute_adapting_temperature+alpha[participant]-rho[participant])));

    theta = lambda[participant] + (1-2*lambda[participant]) .* Phi(interval_sign .* (beta[participant] .* adapted_stimulus_representation) - kappa[participant]);
  }
}
model{
  //Priors
  mu[1] ~ normal(0,2);
  mu[2] ~ normal(-2,1);
  mu[3] ~ normal(0,1);
  mu[4] ~ normal(-4,1);
  mu[5] ~ normal(0,1);

  tau[1] ~ normal(0,2);
  tau[2] ~ normal(0,1);
  tau[3] ~ normal(0,1);
  tau[4] ~ normal(0,1);
  tau[5] ~ normal(0,1);

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
