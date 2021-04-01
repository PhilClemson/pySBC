transformed data {
  vector[2] mu;
  real sigma1;
  real sigma2;
  real rho;
  vector[2] z_sim;

  matrix[2,2] Sigma;

  mu[1] = 0.0;
  mu[2] = 3.0;

  rho = 0.5;
  sigma1 = 1.0;
  sigma2 = 2.0;

  Sigma[1][1] = sigma1 * sigma1;
  Sigma[1][2] = rho * sigma1 * sigma2;
  Sigma[2][1] = rho * sigma1 * sigma2;
  Sigma[2][2] = sigma2 * sigma2;
  
  z_sim = multi_normal_rng(mu, Sigma);
}

parameters {
  vector[2] z;
}

model {
  z ~ multi_normal(mu, Sigma);
}

generated quantities {
  // The means of these quantities will give the difference between 
  // estimated marginal variances and correlation and their true values.
  // If everything is going correctly then these values should be with
  // the respective MCMC-SE of zero.
  real delta_var1;
  real delta_var2;
  real delta_corr;
  int<lower=0,upper=1> lt_sim[2];

  delta_var1 = square(z[1] - mu[1]) - sigma1 * sigma1;
  delta_var2 = square(z[2] - mu[2]) - sigma2 * sigma2;
  delta_corr = (z[1] - mu[1]) * (z[2] - mu[2]) / (sigma1 * sigma2) - rho;
  
  lt_sim[1] = z[1] < z_sim[1];
  lt_sim[2] = z[2] < z_sim[2];
}
