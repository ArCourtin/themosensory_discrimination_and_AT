//This Stan program implements a hierarchical version of the "partial habituation - mixed coding" model of thermosensory discrimination
//Licence: MIT
//Author: Arthur S. Courtin

data{
  int N;
  int P;

  vector[N] recorded_baseline_temperature;
  vector[N] absolute_target_temperature;
  vector[N] absolute_adapting_temperature;
  array[N] int<lower=0,upper=1> choice_accuracy;
  array[N] int<lower=1,upper=P> participant;
}
transformed data{
  vector[N] deviation_from_adapting_temperature = absolute_target_temperature - absolute_adapting_temperature;
  vector[N] relative_adapting_temperature = absolute_adapting_temperature - recorded_baseline_temperature;
  int M=3;
}
parameters{
  vector[M] mu;
  vector<lower=0>[M] tau;
  matrix[M,P] z;
  cholesky_factor_corr[M] L;
}
transformed parameters{
  vector[P] alpha;
  vector[P] sigma;
  vector[P] lambda;
  vector[N] beta;
  vector[N] theta;
  
  {
    matrix[P,M] delta_participant = (diag_pre_multiply(tau, L) * z)';

    alpha = exp(mu[1] + delta_participant[,1]);
    sigma = exp(mu[2] + delta_participant[,2]);
    lambda = .5 * inv_logit(mu[3] + delta_participant[,3]);
    
    vector[N] centered_stimulus = deviation_from_adapting_temperature - alpha[participant];
    beta = 3/(2+sigma[participant]+alpha[participant]-relative_adapting_temperature);
    vector[N] weight = inv_logit(centered_stimulus*100);
    
    theta = (1-weight) .* 0.5 + (weight-lambda[participant]) .* Phi(beta .* centered_stimulus);  
  }
}
model{
  //Priors
  mu[1] ~ normal(-2,1);
  mu[2] ~ normal(1,1);
  mu[3] ~ normal(-4,1);
  
  tau ~ normal(0,1);
  
  L ~ lkj_corr_cholesky(1);

  to_vector(z) ~ std_normal();
  
  //Likelihood
  choice_accuracy ~ bernoulli(theta);
}
generated quantities{
  vector[N] log_lik;
  
  for(n in 1:N){
    log_lik[n] = bernoulli_lpmf(choice_accuracy[n]|theta[n]);
  }
}
