data {
  int<lower=0> I;
  int<lower=0> J;
}

transformed data {
  int<lower=0, upper=1> y_sim[I, J];
  real<lower=0> sigma_theta_sim;
  vector[J] theta_sim;

  real<lower=0> sigma_a_sim;
  vector<lower=0>[I] a_sim;

  real mu_b_sim;
  real<lower=0> sigma_b_sim;
  vector[I] b_sim;
  sigma_theta_sim = fabs(cauchy_rng(0, 2));
  for (j in 1:J){
    theta_sim[j] = normal_rng(0, sigma_theta_sim);
  } 

  sigma_a_sim = fabs(cauchy_rng(0, 2));

  mu_b_sim = normal_rng(0, 5);
  sigma_b_sim = fabs(cauchy_rng(0, 2));
  
  for (i in 1:I){
    a_sim[i] = lognormal_rng(0, sigma_a_sim);
    b_sim[i] = normal_rng(mu_b_sim, sigma_b_sim);
  }

  for (i in 1:I)
    y_sim[i] = bernoulli_logit_rng(a_sim[i] * (theta_sim - b_sim[i]));
}

parameters {
  real<lower=0> sigma_theta;
  vector[J] theta;

  real<lower=0> sigma_a;
  vector<lower=0>[I] a;

  real mu_b;
  real<lower=0> sigma_b;
  vector[I] b;
}

model {
  sigma_theta ~ cauchy(0, 2);
  theta ~ normal(0, sigma_theta); 

  sigma_a ~ cauchy(0, 2);
  a ~ lognormal(0, sigma_a);

  mu_b ~ normal(0, 5);
  sigma_b ~ cauchy(0, 2);
  b ~ normal(mu_b, sigma_b);

  for (i in 1:I)
    y_sim[i] ~ bernoulli_logit(a[i] * (theta - b[i]));
}

generated quantities {
    int<lower=0,upper=1> lt_sim[J+2*I+4];
    lt_sim[1] = sigma_theta < sigma_theta_sim;
    for (j in 1:J){
      lt_sim[j+1] =  theta[j] < theta_sim[j];
    }
    lt_sim[J+2] = sigma_a < sigma_a_sim;
    for (i in 1:I){
      lt_sim[i+J+2] = a[i] < a_sim[i];
      lt_sim[i+J+I+4] = b[i] < b_sim[i];
    }
    lt_sim[J+I+3] = mu_b < mu_b_sim;
    lt_sim[J+I+4] = sigma_b < sigma_b_sim;
}
