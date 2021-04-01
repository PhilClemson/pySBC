data {
  int<lower=1> N;
  real x[N];
}

transformed data {
  int k_sim[N];
  real<lower=0> rho_sim;
  real<lower=0> alpha_sim;
  vector[N] f_tilde_sim;
  vector[N] f_sim;
  
  rho_sim = gamma_rng(25, 4);
  alpha_sim = fabs(normal_rng(0, 2));
  for (n in 1:N){
    f_tilde_sim[n] = normal_rng(0, 1);
  }
  
  {
    matrix[N, N] cov_sim =   cov_exp_quad(x, alpha_sim, rho_sim)
                       + diag_matrix(rep_vector(1e-10, N));
    matrix[N, N] L_cov_sim = cholesky_decompose(cov_sim);
    f_sim = L_cov_sim * f_tilde_sim;
  }
  
  for (n in 1:N){
    k_sim[n] = poisson_log_rng(f_sim[n]);
  }
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  vector[N] f_tilde;
}

transformed parameters {
  vector[N] f;
  {
    matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
                       + diag_matrix(rep_vector(1e-10, N));
    matrix[N, N] L_cov = cholesky_decompose(cov);
    f = L_cov * f_tilde;
  }
}

model {
  rho ~ gamma(25, 4);//alpba rho need to be positive
  alpha ~ normal(0, 2);
  f_tilde ~ normal(0, 1);

  k_sim ~ poisson_log(f);
}

generated quantities {
    int<lower=0,upper=1> lt_sim[N+2];
    lt_sim[1] = rho < rho_sim;
    lt_sim[2] = alpha < alpha_sim;
    for (n in 3:(N+2)){
      lt_sim[n] = f_tilde[n-2] < f_tilde_sim[n-2];
    }
}
