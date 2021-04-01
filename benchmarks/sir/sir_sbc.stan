// Simple SIR model inspired by the presentation in
// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3380087/pdf/nihms372789.pdf

functions {
  
  // theta[1] = beta, water contact rate
  // theta[2] = kappa, C_{50}
  // theta[3] = gamma, recovery rate
  // theta[4] = xi, bacteria production rate
  // theta[5] = delta, bacteria removal rate
  real[] simple_SIR(real t,
                    real[] y,
                    real[] theta,
                    real[] x_r,
                    int[] x_i) {

    real dydt[4];

    dydt[1] = - theta[1] * y[4] / (y[4] + theta[2]) * y[1];
    dydt[2] = theta[1] * y[4] / (y[4] + theta[2]) * y[1] - theta[3] * y[2];
    dydt[3] = theta[3] * y[2];
    dydt[4] = theta[4] * y[2] - theta[5] * y[4];

    return dydt;
  }
}

data {
  int<lower=0> N_t;
  real t[N_t];
  real y0[4];
}

transformed data {
  real t0 = 0;
  real<lower=0> kappa = 1000000;

  real x_r[0];
  int x_i[0];
  
  real t_sim[N_t];
  int stoi_hat_sim[N_t];
  real B_hat_sim[N_t];
  real<lower=0> beta_sim;
  real<lower=0> gamma_sim;
  real<lower=0> xi_sim;
  real<lower=0> delta_sim;
  
  int y_pos = 0; // flag to indicate if all values of y_sim are positive
  real<lower=0> y_sim[N_t, 4];
  while (y_pos == 0){
    // sample random param values until y_sim has all positive values
    beta_sim = fabs(cauchy_rng(0, 2.5));
    gamma_sim = fabs(cauchy_rng(0, 1));
    xi_sim = fabs(cauchy_rng(0, 25));
    delta_sim = fabs(cauchy_rng(0, 1));
  
    {
      real theta_sim[5] = {beta_sim, kappa, gamma_sim, xi_sim, delta_sim};
      y_sim = integrate_ode_rk45(simple_SIR, y0, t0, t, theta_sim, x_r, x_i);
    }
    
    y_pos = 1;
    for (n in 1:N_t){
      if (y_sim[n,1] < 0){
        y_pos = 0;
      }
    }
  }
  
  real rate_sim = y0[1] - y_sim[1, 1];
  if (rate_sim <= 0){
    stoi_hat_sim[1] = 0;
  }
  else{
    stoi_hat_sim[1] = poisson_rng(rate_sim);
  }
  for (n in 2:N_t){
    rate_sim = y_sim[n - 1, 1] - y_sim[n, 1];
    if (rate_sim <= 0){
      stoi_hat_sim[n] = 0;
    }
    else{
      stoi_hat_sim[n] = poisson_rng(rate_sim);
    }
  }
  
  B_hat_sim = lognormal_rng(log(col(to_matrix(y_sim), 4)), 0.15);
}

parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> xi;
  real<lower=0> delta;
}

transformed parameters {
  real<lower=0> y[N_t, 4];
  {
    real theta[5] = {beta, kappa, gamma, xi, delta};
    y = integrate_ode_rk45(simple_SIR, y0, t0, t, theta, x_r, x_i);
  }
}
  
model {
  beta ~ cauchy(0, 2.5);
  gamma ~ cauchy(0, 1);
  xi ~ cauchy(0, 25);
  delta ~ cauchy(0, 1);

  real rate = y0[1] - y[1,1];
  if (rate > 0){
    stoi_hat_sim[1] ~ poisson(rate);
  }
  for (n in 2:N_t){
    rate = y[n-1,1] - y[n, 1];
    if (rate > 0){
      stoi_hat_sim[n] ~ poisson(rate);
    }
  }

  B_hat_sim ~ lognormal(log(col(to_matrix(y), 4)), 0.15);
  
}

generated quantities {
    int<lower=0,upper=1> lt_sim[4];
    lt_sim[1] = beta < beta_sim;
    lt_sim[2] = gamma < gamma_sim;
    lt_sim[3] = xi < xi_sim;
    lt_sim[4] = delta < delta_sim;
}
