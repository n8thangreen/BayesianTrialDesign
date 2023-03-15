
#'
#' @export
normal_distn <- function(mu = 1, sigma = 1, .data = NULL) 
  .data$prior <- tibble::lst(mu, sigma)

#'
#' @export
beta_distn <- function(a = 1, b = 1, .data = NULL) 
  .data$prior <- tibble::lst(a, b)


#'
#' @export
simulate_joint <- function(input, no_of_sim = 1000) {
  do.call(simulate_joint_, c(input, no_of_sim))
}


#' simulate_joint_
#' 
#' @param p_control 
#' @param p_treatment 
#' @param prob_accept_ha 
#' @param total_sample_size 
#' @param no_of_sim 
#' @param alternative 
#' @param prior_p_control 
#' @param prior_log_odds_ratio 
#' @param n_iter 
#' @param n_burnin 
#' @param n_thin
#' 
#' @examples
#' 
#' list(p_control = 0.5,
#'      p_treatment = 0.65,
#'      prob_accept_ha = 0.95,
#'      total_sample_size = 20,
#'      no_of_sim = 100,
#'      alternative = "greater",
#'      prior_p_control = beta(1,1),
#'      prior_log_odds_ratio = normal(0,1)) |> 
#'   simulate_joint(no_of_sim = 1000)
#' @export
#' 
simulate_joint_ <- function(p_control = 0.5,
                            p_treatment = NA,
                            odds_ratio = NA,
                            prob_accept_ha = 0.95,
                            total_sample_size = 20,
                            randomisation_ratio = 1,
                            no_of_sim = 100,
                            alternative = "greater",
                            prior_p_control = beta_prior(1,1),
                            prior_log_odds_ratio = normal_prior(0,1),
                            n_iter = 1e3,
                            n_burnin = 3e1,
                            n_thin = 2e1) {
  
  alt <- if (alternative == "greater") {1} else {-1}
  
  if (is.na(p_treatment))
    p_treatment <- boot::inv.logit(log(odds_ratio) + boot::logit(p_control))
  
  # number of successes
  z1 <- rbinom(n = no_of_sim, size = total_sample_size, prob = p_control)
  z2 <- rbinom(n = no_of_sim, size = total_sample_size, prob = p_treatment)
  
  ##########
  # Stan
  
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  rand_prop <- randomisation_ratio/(randomisation_ratio + 1)
  
  dat_input <-
    list(
      nsim = no_of_sim,
      narms = 2,
      N = c(rand_prop*total_sample_size, (1 - rand_prop)*total_sample_size),
      z = cbind(z1, z2),
      a = prior_p_control$a,
      b = prior_p_control$b,
      mu = prior_log_odds_ratio$mu,
      sigma = prior_log_odds_ratio$sigma)
  
  params <- c("alpha", "beta", "pc", "pt", "thetaRR", "thetaDiff")
  
  ##############
  # run MCMC
  
  out <- rstan::stan(
    data = dat_input,
    pars = params,
    file = here::here("stan", "correlated-p.stan"),
    chains = 1,
    iter = n_iter,
    warmup = n_burnin,
    thin = n_thin,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 20))
  
  stan_output <- rstan::extract(out)
  
  save(stan_output, file = here::here("data", "stan_output_correlated.RData"))
  
  tibble::lst(
    prob_gt_0 = colMeans(alt*stan_output$thetaDiff > 0),
    positives = prob_gt_0 >= prob_accept_ha,        # record if trial is "positive"
    study_power = mean(as.numeric(positives)))
}



