
# Bayesian trial design 
# independent frequentist
# by Hakim Dehbi


library(bayesplot)
library(ggplot2)
library(purrr)

n.sim <- 100

#prior for control group is beta(a1, b1)
a1 <- 1
b1 <- 1
#prior for treatment group is beta(a2, b2)
a2 <- 1
b2 <- 1

#control event/total is z1/N1
p1 <- 0.5
N1 <- 20

#treatment event/total is z2/N2.
p2 <- 0.65
N2 <- N1

n_samples <- 1000#000

z1 <- rbinom(n=n.sim, size=N1, prob=p1)
z2 <- rbinom(n=n.sim, size=N2, prob=p2)

# simulate from posterior for control group
# posterior samples by obs data samples
pc <-
  map_dfc(z1, ~rbeta(n = n_samples,
                     shape1 = a1+.x,
                     shape2 = b1+N1-.x))
pt <-
  map_dfc(z2, ~rbeta(n = n_samples, 
                     shape1 = a2+.x, 
                     shape2 = b2+N2-.x))

thetaRR_ <- pt/pc       # Relative risk
thetaDiffv <- pt - pc   # Absolute risk difference

prob_gt_0 <- colMeans(thetaDiff > 0)
positives <- prob_gt_0 >= 0.95    # record if trial is "positive"
power <- mean(as.numeric(positives))
power


### verify with bayesCT
library(bayesCT)

##TODO: where is beta_prior()?
value <-
  binomial_outcome(p_treatment = 0.15, 
                   p_control   = 0.10) |>
  study_details(total_sample_size     = 1200, 
                study_period          = 50,
                interim_look          = NULL,
                prop_loss_to_followup = 0) |>
  hypothesis(delta                 = 0, 
             futility_prob         = 0.0, 
             prob_accept_ha        = 0.95,
             expected_success_prob = 1, 
             alternative           = "greater") |>
  beta_prior(a0 = 1, 
             b0 = 1) |>
  simulate(no_of_sim = 1000)

str(value)
value$power


df <- data.frame(characteristic = rep(c("power", "type-1"), each=4),
                 sample_size = rep(c(400, 600, 800, 1200),2),
                 value = c(0.44, 0.57, 0.70, 0.83, 0.04, 0.06, 0.05, 0.06))
df

p <-
  ggplot(df, aes(x=sample_size, y=value, group=characteristic)) +
  geom_line(aes(color=characteristic)) +
  geom_point(aes(color=characteristic)) +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))
p

