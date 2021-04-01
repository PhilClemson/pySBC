data {
  int<lower=0> T;
  real<lower=0> sigma1; 
}

transformed data {
  real y_sim[T];
  real sigma_sim[T];
  real<lower=-100, upper=100> mu_sim;
  real<lower=0, upper=100> alpha0_sim;          
  real<lower=0, upper=1> alpha1_sim;
  real<lower=0, upper=(1-alpha1_sim)> beta1_sim; 
  mu_sim = uniform_rng(-100, 100);
  alpha0_sim = uniform_rng(0, 100);
  alpha1_sim = uniform_rng(0, 1);
  beta1_sim = uniform_rng(0,1-alpha1_sim);
  sigma_sim[1] = sigma1;
  y_sim[1] = normal_rng(mu_sim, sigma_sim[1]);
  for (t in 2:T) {
    sigma_sim[t] = sqrt(  alpha0_sim
                    + alpha1_sim * square(y_sim[t - 1] - mu_sim)
                    + beta1_sim * square(sigma_sim[t - 1]));
    y_sim[t] = normal_rng(mu_sim, sigma_sim[t]);
  }
}

parameters {
  real<lower=-100, upper=100> mu; 
  real<lower=0, upper=100> alpha0;          
  real<lower=0, upper=1> alpha1;
  real<lower=0, upper=(1-alpha1)> beta1; 
}

model {
  real sigma[T];
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(  alpha0
                    + alpha1 * square(y_sim[t - 1] - mu)
                    + beta1 * square(sigma[t - 1]));

  y_sim ~ normal(mu, sigma);
}

generated quantities {
    int<lower=0,upper=1> lt_sim[4];
    lt_sim[1] = mu < mu_sim;
    lt_sim[2] = alpha0 < alpha0_sim;
    lt_sim[3] = alpha1 < alpha1_sim;
    lt_sim[4] = beta1 < beta1_sim;
}
