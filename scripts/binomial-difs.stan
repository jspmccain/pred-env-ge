// Two Binomial Distributions Stan Model

data {
    int<lower=0> N1;  // Number of observations for the first distribution
    int<lower=0> Y1;
    int<lower=0> N2;  // Number of observations for the first distribution
    int<lower=0> Y2;
}

parameters {
    real<lower=0, upper=1> theta1;  // Probability parameter for the first binomial distribution
    real<lower=0, upper=1> theta2;  // Probability parameter for the second binomial distribution
}

model {
    // Prior distributions for the parameters
    theta1 ~ beta(1, 1);  // Beta(1, 1) is the uniform prior on [0, 1]
    theta2 ~ beta(1, 1);

    target += binomial_lupmf(N1 | Y1, theta1);
    target += binomial_lupmf(N2 | Y2, theta2);
}

generated quantities {
    real theta_difference;  // Difference between theta parameters

    theta_difference = theta2 - theta1;  // Compute the difference
}
