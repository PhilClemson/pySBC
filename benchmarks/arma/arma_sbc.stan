data {
  int<lower=1> T; // number of observations
}

transformed data {
  real y_sim[T];
  real mu_sim;
  real phi_sim;
  real theta_sim;
  real<lower=0> sigma_sim;
  vector[T] nu_sim;  // prediction for time t
  
  mu_sim = normal_rng(0, 10);
  phi_sim = normal_rng(0, 2);
  theta_sim = normal_rng(0, 2);
  sigma_sim = fabs(cauchy_rng(0,2.5));
  
  nu_sim[1] = mu_sim + phi_sim * mu_sim; // assume err[0] == 0
  y_sim[1] = normal_rng(nu_sim[1], sigma_sim);
  for (t in 2:T) {
    nu_sim[t] = mu_sim + phi_sim * y_sim[t - 1] + theta_sim * (y_sim[t-1] - nu_sim[t-1]);
    y_sim[t] = normal_rng(nu_sim[t], sigma_sim);
  }
}

parameters {
  real mu;             // mean coefficient
  real phi;            // autoregression coefficient
  real theta;          // moving average coefficient
  real<lower=0> sigma; // noise scale
}

model {
  vector[T] nu;  // prediction for time t
  
  mu ~ normal(0, 10);
  phi ~ normal(0, 2);
  theta ~ normal(0, 2);
  sigma ~ cauchy(0, 2.5);
  
  nu[1] = mu + phi * mu; // assume err[0] == 0
  for (t in 2:T) {
    nu[t] = mu + phi * y_sim[t - 1] + theta * (y_sim[t-1] - nu[t-1]);
  }
  
  y_sim ~ normal(nu, sigma);
}

generated quantities {
    int<lower=0,upper=1> lt_sim[4];
    lt_sim[1] = mu < mu_sim;
    lt_sim[2] = phi < phi_sim;
    lt_sim[3] = theta < theta_sim;
    lt_sim[4] = sigma < sigma_sim;
}