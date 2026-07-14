//This Stan program implements a hierarchical version of the non-mechanistic model of thermosensory discrimination
//The condition-specific alpha/beta are free across the ordered adapting-temperature axis, but their group-level
//profiles are given a second-order random-walk (RW2) prior: the cheap default is linear in adapting temperature,
//with departures from linearity penalized by an estimated smoothness scale. This keeps the model strictly more
//flexible than the mechanistic accounts while trimming the wiggle-room no plausible DGP would use, so its role as
//a competitor in model comparison is a fair one.
//Response-coded: the outcome is the second-interval choice of a 2IFC task, not accuracy.
//Licence: MIT
//Author: Arthur S. Courtin

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

  int M=12;                                       // participant-hierarchy dimensions: 5 alpha, 5 beta, lambda, kappa
  int C=5;                                         // adapting-temperature conditions
}
parameters{
  //Group-level RW2 profiles across the ordered adapting-temperature axis (log scale)
  real a0;                                        // alpha: level at condition 1
  real a1;                                        // alpha: initial slope (per AT step)
  vector[C-2] a_curv;                             // alpha: standardized curvature increments
  real<lower=0> sigma_a;                          // alpha: curvature scale (smoothness)

  real b0;                                        // beta: level at condition 1
  real b1;                                        // beta: initial slope (per AT step)
  vector[C-2] b_curv;                             // beta: standardized curvature increments
  real<lower=0> sigma_b;                          // beta: curvature scale (smoothness)

  real mu_lambda;                                 // group-level lapse (pre-transform)
  real mu_kappa;                              // group-level interval bias

  //Participant-level non-centered hierarchy (order: 5 alpha, 5 beta, lambda, kappa)
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  vector[C] f_a;                                   // group-level log-alpha profile across AT
  vector[C] f_b;                                   // group-level log-beta profile across AT
  matrix[C,P] alpha;
  matrix[C,P] beta;
  row_vector[P] lambda;
  row_vector[P] kappa;
  vector[N] theta;

  //RW2 group profiles: two free anchors (level + first slope), curvature (second differences) penalized by sigma
  f_a[1] = a0;
  f_a[2] = a0 + a1;
  for(c in 3:C){
    f_a[c] = 2*f_a[c-1] - f_a[c-2] + sigma_a * a_curv[c-2];
  }
  f_b[1] = b0;
  f_b[2] = b0 + b1;
  for(c in 3:C){
    f_b[c] = 2*f_b[c-1] - f_b[c-2] + sigma_b * b_curv[c-2];
  }

  {
    matrix[M,P] delta_participant = diag_pre_multiply(tau, L) * z;

    for(idx in 1:5){
      alpha[idx] = exp(f_a[idx] + delta_participant[idx]);
      beta[idx] = exp(f_b[idx] + delta_participant[idx+5]);
    }
    lambda = .5 * inv_logit(mu_lambda + delta_participant[11]);
    kappa = mu_kappa + delta_participant[12];

    for(n in 1:N){
      real centered_stimulus = deviation_from_adapting_temperature[n] - alpha[adapting_temperature_idx[n],participant[n]];
      real stimulus_representation = centered_stimulus * inv_logit(100*centered_stimulus);

      theta[n] = lambda[participant[n]] + (1-2*lambda[participant[n]]) * Phi(interval_sign[n] * beta[adapting_temperature_idx[n],participant[n]] * stimulus_representation - kappa[participant[n]]);
    }
  }
}
model{
  //Priors on the RW2 alpha profile (a0 matches the original alpha-mean prior)
  a0 ~ normal(-2,1);
  a1 ~ normal(0,0.5);
  a_curv ~ std_normal();
  sigma_a ~ normal(0,0.5);

  //Priors on the RW2 beta profile (b0 matches the original beta-mean prior)
  b0 ~ normal(0,1);
  b1 ~ normal(0,0.5);
  b_curv ~ std_normal();
  sigma_b ~ normal(0,0.5);

  mu_lambda ~ normal(-4,1);
  mu_kappa ~ normal(0,1);

  tau ~ normal(0,1);

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
