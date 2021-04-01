data {
 int<lower = 0> N;
}

transformed data {
  vector[N] y_sim;
  vector[2] mu_sim;
  real<lower=0> sigma_sim[2];
  real<lower=0, upper=1> theta_sim;
  
  for (n in 1:2){
    sigma_sim[n] = fabs(normal_rng(0, 2));
    mu_sim[n] = normal_rng(0, 2);
  }
  theta_sim = beta_rng(5, 5);
  
  // use Bernoulli trial to determine which samples belong to each distribution
  for (n in 1:N){
    real b;
    b = bernoulli_rng(theta_sim);
    if (b == 1){
      y_sim[n] = normal_rng(mu_sim[1], sigma_sim[1]);
    }
    else{
      y_sim[n] = normal_rng(mu_sim[2], sigma_sim[2]);
    }
  }
}

parameters {
  vector[2] mu;
  real<lower=0> sigma[2];
  real<lower=0, upper=1> theta;
}

model {
 sigma ~ normal(0, 2);
 mu ~ normal(0, 2);
 theta ~ beta(5, 5);
 for (n in 1:N)
   target += log_mix(theta,
                     normal_lpdf(y_sim[n] | mu[1], sigma[1]),
                     normal_lpdf(y_sim[n] | mu[2], sigma[2]));
}

generated quantities {
    int<lower=0,upper=1> lt_sim[5];
    lt_sim[1] = mu[1] < mu_sim[1];
    lt_sim[2] = mu[2] < mu_sim[2];
    lt_sim[3] = sigma[1] < sigma_sim[1];
    lt_sim[4] = sigma[2] < sigma_sim[2];
    lt_sim[5] = theta < theta_sim;
}
