
#' @importFrom tibble lst
#' @export
normal_distn <- function(mu = 1, sigma = 1, .data = NULL) 
  .data$prior <- tibble::lst(mu, sigma)

#' @importFrom tibble lst
#' @export
beta_distn <- function(a = 1, b = 1, .data = NULL) 
  .data$prior <- tibble::lst(a, b)


#' @importFrom tibble lst
#' @export
simulate_joint <- function(input, no_of_sim = 1000) {
  input <- modifyList(tibble::lst(no_of_sim), input)
  do.call(simulate_joint_, input)
}


#' simulate_joint_
#' 
#' @param p_control Control probability
#' @param p_treatment  Treatment probability
#' @param prob_accept_ha Probability accept alternative hypothesis
#' @param total_sample_size Total sample size
#' @param no_of_sim Number of simulations
#' @param alternative Alternative test
#' @param prior_p_control Prior for control probability
#' @param prior_log_odds_ratio Prior for log-odds ratio
#' @param n_iter Number of iterations
#' @param n_burnin Number of burn-in
#' @param n_thin Number of thinning
#' 
#' @examples
#' 
#' list(p_control = 0.5,
#'      p_treatment = 0.65,
#'      prob_accept_ha = 0.95,
#'      total_sample_size = 20,
#'      no_of_sim = 100,
#'      alternative = "greater",
#'      prior_p_control = beta_distn(1,1),
#'      prior_log_odds_ratio = normal_distn(0,1)) |> 
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
                            prior_p_control = beta_distn(1,1),
                            prior_log_odds_ratio = normal_distn(0,1),
                            n_iter = 1e3,
                            n_burnin = 3e1,
                            n_thin = 2e1) {
  
  alt <- if (alternative == "greater") {1} else {-1}
  
  if (is.na(p_treatment))
    p_treatment <- boot::inv.logit(log(odds_ratio) + boot::logit(p_control))
  
  rand_prop <- randomisation_ratio/(randomisation_ratio + 1)
  n_control <- round(rand_prop*total_sample_size, 0)
  n_treatment <- total_sample_size - n_control
  
  # number of successes
  z_control <- rbinom(n = no_of_sim, size = n_control, prob = p_control)
  z_treatment <- rbinom(n = no_of_sim, size = n_treatment, prob = p_treatment)

  ##########
  # Stan
  
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  dat_input <-
    list(
      nsim = no_of_sim,
      narms = 2,
      N = c(n_control, n_treatment),
      z = cbind(z_control, z_treatment),
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
  
  browser()
  tibble::lst(
    prob_gt_0 = colMeans(alt*stan_output$thetaDiff > 0),
    positives = prob_gt_0 >= prob_accept_ha,        # record if trial is "positive"
    study_power = mean(as.numeric(positives)))
}



