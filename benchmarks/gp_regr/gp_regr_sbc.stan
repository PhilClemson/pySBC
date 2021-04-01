data {
  int<lower=1> N;
  real x[N];
}

transformed data {
  vector[N] y_sim;
  real<lower=0> rho_sim;
  real<lower=0> alpha_sim;
  real<lower=0> sigma_sim;
  
  rho_sim = gamma_rng(25, 4);
  alpha_sim = fabs(normal_rng(0, 2));
  sigma_sim = fabs(normal_rng(0, 1));
  matrix[N, N] cov_sim =   cov_exp_quad(x, alpha_sim, rho_sim)
                     + diag_matrix(rep_vector(sigma_sim, N));
  matrix[N, N] L_cov_sim = cholesky_decompose(cov_sim);

  y_sim = multi_normal_cholesky_rng(rep_vector(0, N), L_cov_sim);
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
                     + diag_matrix(rep_vector(sigma, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);

  rho ~ gamma(25, 4);
  alpha ~ normal(0, 2);
  sigma ~ normal(0, 1);

  y_sim ~ multi_normal_cholesky(rep_vector(0, N), L_cov);
}

generated quantities {
    int<lower=0,upper=1> lt_sim[3];
    lt_sim[1] = rho < rho_sim;
    lt_sim[2] = alpha < alpha_sim;
    lt_sim[3] = sigma < sigma_sim;
}
