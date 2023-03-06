// correlated probabilities

data {
  int<lower=0> nsim;
  int<lower=0> narms;
  int<lower=0> N[narms];
  int z[nsim, narms];
  real a;
  real b;
  real mu;
  real sigma;
}

parameters {
  real<lower=0, upper=1> q[nsim];
  real beta[nsim];
}

transformed parameters {
  real alpha[nsim];
  
  for (i in 1:nsim) {
    alpha[i] = log(q[i]/(1 - q[i]));
  }
}

model {
  for (i in 1:nsim) {
    for (j in 1:narms) {
      z[i,j] ~ binomial_logit(N[j], alpha[i] + beta[i] * (j-1));
    }
    
    // priors
    q[i] ~ beta(a, b);            // control probability
    beta[i] ~ normal(mu, sigma);  // log-odds ratio
  }
}

generated quantities {
  real thetaRR[nsim];
  real thetaDiff[nsim];
  real<lower=0, upper=1> pt[nsim];
  real<lower=0, upper=1> pc[nsim];
  
  for (i in 1:nsim) {
    pc[i] = inv_logit(alpha[i]);
    pt[i] = inv_logit(alpha[i] + beta[i]);
    
    thetaRR[i] = pt[i]/pc[i];         // relative risk
    thetaDiff[i] = pt[i] - pc[i];     // absolute risk difference
  }
}
