// Hierarchical model examining how protein abundance-flux relationships covary, and explaining variance in this relationship using enzyme-level factors

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of observations
  int<lower=1> K; // number of unique enzymes
  int<lower=1,upper=K> id[N]; // vector of group indices
  vector[N] flux; // response variables
  vector[N] protein_amounts;
  vector[K] mean_protein_amount;
  vector[K] substrate_connectivity;
  vector[K] delta_g;
}

// Parameters and hyperparameters.
parameters {
  real gamma_0;
  real gamma_1; // first population level coefficient
  real gamma_2; // second population level coefficient
  real gamma_3;
  real<lower = 0> tau;            // Variance of the population model.
  real<lower = 0> sigma;          // Variance of the likelihood.
  vector[K] beta_0; // intercept for the observation level
  vector[K] beta_1; // slope for the observation level
}


model {
  // Priors
  // gamma_0 ~ normal(0, 10);
  // gamma_1 ~ normal(0, 10);
  // gamma_2 ~ normal(0, 10);
  // gamma_3 ~ normal(0, 10);
  // tau ~ normal(0, 1);
  // sigma ~ normal(0, 1);
  // beta_0 ~ normal(0, 10);
  
  // Likelihood 
  for(k in 1:K){
    beta_1[k] ~ normal(gamma_0 + gamma_1 * mean_protein_amount[k] + gamma_2 * substrate_connectivity[k] + gamma_3 * delta_g[k], tau);
  }
  
  for(n in 1:N){
    flux[n] ~ normal(beta_0[id[n]] + beta_1[id[n]]*protein_amounts[n], sigma);
  }

}


