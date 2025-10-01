fit_censored_exponential_model <- function(prior_alpha, prior_beta, data_set, cens_indicators) {
  
  # Define the log-likelihood function for the exponential distribution
  log_likelihood <- function(lambda, data, cens_indicators) {
    #lambda <- exp(log_lambda)
    n <- length(data)
    n_cens <- sum(cens_indicators)
    n_not_cens <- n - n_cens
    cens_id <- which(cens_indicators == 1)
    nocens_id <- which(cens_indicators == 0)
    
    # Likelihood for the censored observations
    cens_likelihood <- ifelse(length(cens_id)>0,sum(pexp(data[cens_id], lambda, log = TRUE)),0)
    
    # Likelihood for the not censored observations
    nocens_likelihood <- ifelse(length(nocens_id)>0,sum(dexp(data[nocens_id], lambda, log = TRUE)),0)
    
    # Combine the likelihoods
    total_likelihood <- cens_likelihood + nocens_likelihood
    return(total_likelihood)
  }
  
  # Define the log-prior function for the lognormal distribution
  log_prior <- function(log_lambda, alpha, beta) {
    dgamma(log_lambda, alpha,beta, log = TRUE)
  }
  
  # Define the log-posterior function
  log_posterior <- function(lambda, args) {
    alpha <- args[[1]]
    beta <- args[[2]]
    data <- args[[3]]
    cens_indicators <- args[[4]]
    
    value <- log_likelihood(lambda, data, cens_indicators) + log_prior(lambda, alpha, beta)
    if(lambda < 0)value = -Inf
    return(value)
  }
  args <- list(prior_alpha, prior_beta, data_set, cens_indicators)
  
  # Metropolis-Hastings MCMC sampling
  initial_value <- 12716898
  proposal_sd <- 10000
  base_sd <- proposal_sd
  iterations <- 20000
  nburn <- 10000
  cur_val <- initial_value
  samples <- numeric(iterations)
  for(i in 1:iterations){
    proposed_val <- cur_val + rnorm(1,0,proposal_sd)
    metrop_alpha <- log_posterior(proposed_val,args)-
      log_posterior(cur_val,args)
    if(log(runif(1)) < metrop_alpha){
      cur_val <- proposed_val
    }
    samples[i] <- cur_val
    if(i > 100){
      proposal_sd <- ifelse(runif(1) < .05,base_sd,
                            sqrt(var(samples[1:i]))/3)
    }
  }

  # Return the MCMC samples as the posterior distribution
  return(samples[(nburn+1):iterations])
}

# Define the log-likelihood function for the exponential distribution
exp_log_likelihood <- function(lambda, data, cens_indicators) {
  #lambda <- exp(log_lambda)
  n <- length(data)
  n_cens <- sum(cens_indicators)
  n_not_cens <- n - n_cens
  cens_id <- which(cens_indicators == 1)
  nocens_id <- which(cens_indicators == 0)
  
  # Likelihood for the censored observations
  cens_likelihood <- ifelse(length(cens_id)>0,sum(pexp(data[cens_id], lambda, log = TRUE)),0)
  
  # Likelihood for the not censored observations
  nocens_likelihood <- ifelse(length(nocens_id)>0,sum(dexp(data[nocens_id], lambda, log = TRUE)),0)
  
  # Combine the likelihoods
  total_likelihood <- cens_likelihood + nocens_likelihood
  return(total_likelihood)
}

gpd_log_likelihood <- function(params, data, censored_indicator) {
  sigma <- params[1]
  xi <- params[2]
  
  # Log likelihood for fully observed data
  observed_mask <- censored_indicator == 0
  observed_data <- data[observed_mask]
  log_likelihood_observed <- dgpd(observed_data,scale = sigma,shape = xi,log = TRUE)
  
  # Log likelihood for left-censored data
  censored_mask <- censored_indicator == 1
  censored_data <- data[censored_mask]
  log_likelihood_censored <- log(pgpd(censored_data,scale = sigma,shape = xi))
  sum(log_likelihood_observed) + sum(log_likelihood_censored)
}
