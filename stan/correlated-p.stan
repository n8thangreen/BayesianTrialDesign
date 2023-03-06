//

data {
  int<lower=0> nsim;
  int<lower=0> N1;
  int<lower=0> N2;
  int z1[nsim];
  int z2[nsim];
  
  real a1;
  real b1;
  real a2;
  real b2;
}

parameters {
  real<lower=0, upper=1> p1[nsim];
  real<lower=0, upper=1> p2[nsim];
}

model {
  
  for (i in 1:nsim) {
    z1[i] ~ binomial(N1, p1[i]);
    z2[i] ~ binomial(N2, p2[i]);
    
    // priors
    p1[i] ~ beta(a1, b1);
    p2[i] ~ beta(a2, b2);
  }
}

generated quantities {
  real thetaRR[nsim];
  real thetaDiff[nsim];
  // real prob_gt_0[nsim];
  // real positives[nsim];
  // real study_power;
  
  for (i in 1:nsim) {
    thetaRR[i] = p2[i]/p1[i];        // Relative risk
    thetaDiff[i] = p2[i] - p1[i];    // Absolute risk difference
  }
}
