functions {
  real half_cauchy_rng(real lb) {
    real p = cauchy_cdf(lb, 0, 5);
    real u = uniform_rng(p, 1);
    real y = 5*tan(pi()*(u - 0.5)); // inverse cauchy_cdf
    return y;
  }
}

data {
  int<lower=0> J;
  real<lower=0> sigma[J];
}

transformed data{
  real y_sim[J];
  real mu_sim;
  real<lower=0> tau_sim;
  real theta_tilde_sim[J];
  real theta_sim[J];
  mu_sim = normal_rng(0, 5);
  tau_sim = half_cauchy_rng(0);
  for (j in 1:J){
    theta_tilde_sim[j] = normal_rng(0, 1);
    theta_sim[j] = mu_sim + tau_sim * theta_tilde_sim[j];
  }
  y_sim = normal_rng(theta_sim, sigma);
}

parameters {
  real mu;
  real<lower=0> tau;
  real theta_tilde[J];
}

transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] = mu + tau * theta_tilde[j];
}

model {
  mu ~ normal(0, 5);
  tau ~ cauchy(0, 5);
  theta_tilde ~ normal(0, 1);
  y_sim ~ normal(theta, sigma);
}

generated quantities {
    int<lower=0,upper=1> lt_sim[J+2];
    lt_sim[1] = mu < mu_sim;
    lt_sim[2] = tau < tau_sim;
    for (n in 3:J+2){
        lt_sim[n] = theta_tilde[n-2] < theta_tilde_sim[n-2];
    }
}
