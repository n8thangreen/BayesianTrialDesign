# example

library(bayesCT)

list(p_control = 0.5,
     p_treatment = 0.65,
     prob_accept_ha = 0.95,
     total_sample_size = 20,
     no_of_sim = 100,
     alternative = "greater",
     prior_p_control = beta_distn(1,1),
     prior_log_odds_ratio = normal_distn(0,1)) |>
  simulate_joint(no_of_sim = 1000)




