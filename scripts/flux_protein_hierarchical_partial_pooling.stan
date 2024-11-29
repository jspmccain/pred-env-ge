// Stan model examining relationship between enzyme abundances and associated rates.

// Function for calculating the linear predictor, *without* the variance.
functions {
  real linear_predictor(real alpha, real beta, real x) {
    real out;
    out = alpha + beta*x;
    return out;
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of observations
  int<lower=1> K; // number of unique enzymes
  int<lower=1,upper=K> id[N]; // vector of enzyme indices
  vector[N] flux; // response variables
  vector[N] protein_amounts; // observed protein amounts
  int<lower=0> N_spaced; 
  vector[N_spaced] protein_amounts_spaced; // evenly spaced protein amounts for posterior predictions
  int<lower=1,upper=K> id_spaced[N_spaced];
}

// Parameters and hyperparameters.
parameters {
  //hyperparameters
  real<lower=-1,upper=1> gamma_1;     // first population level coefficient
  real<lower = 0> tau;                // Variance of the distribution of coefficients.
  real<lower = 0> sigma;              // Variance of the likelihood
  vector[K] beta_0;                   // intercept for the observation level, for each K enzyme
  vector<lower=-1,upper=1>[K] beta_1; // slope for the observation level
}

model {
  // setting the hyperpriors
  // gamma_1 ~ normal(0, 1);
  // tau ~ normal(0, 1);
  // sigma ~ normal(0, 1);
  
  beta_1 ~ normal(gamma_1, tau);
  
  for(n in 1:N){
    flux[n] ~ normal(beta_0[id[n]] + beta_1[id[n]]*protein_amounts[n], sigma);
  }
}

generated quantities {
  // Simulate replicated data for posterior predictive checks
  vector[N] y_rep; // Replicated data
  vector[N] y_rep_r_sq;
  vector[3000] group_dist;  // Group-level distribution

  for (n in 1:N) {
    y_rep[n] = normal_rng(beta_0[id[n]] + beta_1[id[n]] * protein_amounts[n], sigma);
    y_rep_r_sq[n] = linear_predictor(beta_0[id[n]], beta_1[id[n]], protein_amounts[n]);
  }
  
  for (j in 1:3000) {
      group_dist[j] = normal_rng(gamma_1, tau);
  }
}



