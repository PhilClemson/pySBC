data {
  int<lower=0> K;
  int<lower=0> T;
}

transformed data {
  real y_sim[T];
  real y_init[2*K]; // initial buffer used to generate first K values in y_sim
  real alpha_sim;
  real beta_sim[K];
  real<lower=0> sigma_sim;
  alpha_sim = normal_rng(0, 10);
  sigma_sim = fabs(cauchy_rng(0, 2.5));
  for (k in 1:K){
    beta_sim[k] = normal_rng(0, 10);
    y_init[k] = normal_rng(alpha_sim, sigma_sim);
  }
  for (t in (K+1):2*K) {
    real mu_sim;
    mu_sim = alpha_sim;
    
    for (k in 1:K)
      mu_sim = mu_sim + beta_sim[k] * y_init[t - k];
    
    y_init[t] = normal_rng(mu_sim, sigma_sim);
    y_sim[t-K] = y_init[t];
  }
  for (t in (K+1):T) {
    real mu_sim;
    mu_sim = alpha_sim;
    
    for (k in 1:K)
      mu_sim = mu_sim + beta_sim[k] * y_sim[t - k];
    
    y_sim[t] = normal_rng(mu_sim, sigma_sim);
  }
}

parameters {
  real alpha;
  real beta[K];
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ cauchy(0, 2.5);
  
  for (t in (K+1):T) {
    real mu;
    mu = alpha;
    
    for (k in 1:K)
      mu = mu + beta[k] * y_sim[t - k];
    
    y_sim[t] ~ normal(mu, sigma);
  }
}

generated quantities {
    int<lower=0,upper=1> lt_sim[K+2];
    lt_sim[1] = alpha < alpha_sim;
    for (k in 2:K+1){
        lt_sim[k] = beta[k-1] < beta_sim[k-1];
    }
    lt_sim[K+2] = sigma < sigma_sim;
}
