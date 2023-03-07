
# Bayesian trial design 
# correlated Bayesian using Stan


library(bayesplot)
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(rstan)
library(shinystan)

n.sim <- 100

# prior for control group is beta(a1, b1)
a1 <- 1
b1 <- 1
# prior for treatment group is beta(a2, b2)
a2 <- 1
b2 <- 1

# control event/total is z1/N1
p1 <- 0.5
N1 <- 20

# treatment event/total is z2/N2.
p2 <- 0.65
N2 <- N1

n_samples <- 1000#000

# number of successes
z1 <- rbinom(n=n.sim, size=N1, prob=p1)
z2 <- rbinom(n=n.sim, size=N2, prob=p2)


##########
# Stan

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat_input <-
  list(
    nsim = n.sim,
    narms = 2,
    N = c(N1,N2),
    z = cbind(z1,z2),
    a = a1,
    b = b1,
    mu = 0,
    sigma = 1)

params <- c("alpha", "beta", "pc", "pt", "thetaRR", "thetaDiff")

n_iter <- 1e3
n_burnin <- 3e1
n_thin <- 2e1 #floor((n_iter - n_burnin)/500)


##############
# run MCMC

out <- stan(data = dat_input,
            pars = params,
            file = here::here("stan", "correlated-p.stan"),
            chains = 1,
            iter = n_iter,
            warmup = n_burnin,
            thin = n_thin,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 20))

stan_output <- extract(out)

save(stan_output, file = here::here("data", "stan_output_correlated.RData"))

prob_gt_0 <- colMeans(stan_output$thetaDiff > 0)
positives <- prob_gt_0 >= 0.95             # record if trial is "positive"
study_power <- mean(as.numeric(positives))
study_power

