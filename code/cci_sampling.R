# cci_sampling.R — paper-consistent implementation
# ------------------------------------------------------------

# =========================
# Utilities
# =========================
rhalfnorm <- function(n, sd) abs(rnorm(n, 0, sd))
dhalfnorm_log <- function(x, sd) ifelse(x > 0, log(2) + dnorm(x, 0, sd, log=TRUE), -Inf)

# =========================
# Exponential model
# =========================

# Prior calibration: choose beta s.t. E[e^{-λT}] = (β/(β+T))^α = ε  =>  β = T / (ε^{-1/α} - 1)
exp_prior_beta_from_alpha <- function(alpha, T, epsilon) {
  T / (epsilon^(-1/alpha) - 1)
}

# Log-likelihood with left-censoring
# If cens[i]==1: contribution is log F(C_i) = log(pexp(C_i, rate=lambda))
# If cens[i]==0: contribution is log f(x_i) = log(dexp(x_i, rate=lambda))
exp_loglik <- function(lambda, x, cens) {
  if (lambda <= 0) return(-Inf)
  cen <- cens == 1
  obs <- cens == 0
  ll <- 0
  if (any(cen)) ll <- ll + sum(pexp(x[cen], rate=lambda, log=TRUE))
  if (any(obs)) ll <- ll + sum(dexp(x[obs], rate=lambda, log=TRUE))
  ll
}

# Gamma(alpha, beta) prior (rate = beta)
exp_logprior <- function(lambda, alpha, beta) {
  if (lambda <= 0) return(-Inf)
  dgamma(lambda, shape=alpha, rate=beta, log=TRUE)
}

exp_logpost <- function(lambda, x, cens, alpha, beta) {
  exp_loglik(lambda, x, cens) + exp_logprior(lambda, alpha, beta)
}

# Random-walk MH for λ with robust defaults: init at prior mean (α/β)
exp_mh <- function(x, cens, alpha, beta, init=NULL, iter=15000, burn=5000, prop_sd=1000) {
  if (is.null(init)) init <- alpha / beta
  cur <- max(init, .Machine$double.eps)
  out <- numeric(iter)
  for (i in seq_len(iter)) {
    prop <- cur + rnorm(1, 0, prop_sd)
    lp_cur <- exp_logpost(cur,  x, cens, alpha, beta)
    lp_prp <- exp_logpost(prop, x, cens, alpha, beta)
    if (log(runif(1)) < (lp_prp - lp_cur)) cur <- prop
    out[i] <- cur
    # light adaptation during burn-in
    if (i %% 200 == 0 && i <= burn) {
      s <- sd(out[1:i])
      if (is.finite(s) && s > 0) prop_sd <- 0.5 * s
    }
  }
  out[(burn+1):iter]
}

# =========================
# GPD (μ = 0) model
# =========================

# CDF F(x) with exponential limit as ξ -> 0
gpd_F <- function(x, xi, sigma) {
  if (sigma <= 0) return(NaN)
  if (abs(xi) < 1e-12) {
    1 - exp(-x / sigma)
  } else {
    z <- 1 + xi * x / sigma
    ifelse(z > 0, 1 - z^(-1/xi), 0)
  }
}

# log-pdf with exponential limit as ξ -> 0
gpd_logf <- function(x, xi, sigma) {
  if (sigma <= 0) return(-Inf)
  if (any(1 + xi * x / sigma <= 0)) return(-Inf)
  if (abs(xi) < 1e-12) {
    -log(sigma) - x / sigma
  } else {
    -log(sigma) - (1/xi + 1) * log1p(xi * x / sigma)
  }
}

# Log-likelihood for mix of left-censored and observed
gpd_loglik <- function(params, x, cens) {
  sigma <- params[1]; xi <- params[2]
  if (sigma <= 0 || xi < 0) return(-Inf)
  obs <- cens == 0
  cen <- cens == 1
  ll <- 0
  if (any(obs)) ll <- ll + sum(gpd_logf(x[obs], xi, sigma))
  if (any(cen)) {
    Fc <- gpd_F(x[cen], xi, sigma)     # x[cen] are censor points C_i
    if (any(Fc <= 0)) return(-Inf)
    ll <- ll + sum(log(Fc))
  }
  ll
}

# Priors: xi ~ HalfNormal(nu_xi), sigma ~ HalfNormal(nu_sigma)
gpd_logprior <- function(params, nu_xi, nu_sigma) {
  sigma <- params[1]; xi <- params[2]
  dhalfnorm_log(xi, nu_xi) + dhalfnorm_log(sigma, nu_sigma)
}

gpd_logpost <- function(params, x, cens, nu_xi, nu_sigma) {
  gpd_loglik(params, x, cens) + gpd_logprior(params, nu_xi, nu_sigma)
}

# Transformed MH: propose on eta = log(sigma) and z = log1p(xi) for stability
gpd_mh <- function(x, cens, nu_xi, nu_sigma,
                   init_sigma=NULL, init_xi=NULL,
                   iter=20000, burn=10000, sd_eta=0.5, sd_z=0.2, adapt_every=250) {
  
  if (is.null(init_xi))    init_xi    <- max(0.03, nu_xi * sqrt(2/pi))  # small positive
  if (is.null(init_sigma)) init_sigma <- 1
  
  eta <- log(max(init_sigma, .Machine$double.eps))
  z   <- log1p(max(init_xi, 0))
  samp <- matrix(NA_real_, nrow=iter, ncol=2)
  
  logtarget <- function(eta, z) {
    sigma <- exp(eta); xi <- exp(z) - 1
    if (sigma <= 0 || xi < 0) return(-Inf)
    # posterior in (sigma, xi) + Jacobian log|J| = log(sigma) + log(1+xi) = eta + z
    gpd_logpost(c(sigma, xi), x, cens, nu_xi, nu_sigma) + eta + z
  }
  
  for (i in seq_len(iter)) {
    eta_p <- eta + rnorm(1, 0, sd_eta)
    z_p   <- z   + rnorm(1, 0, sd_z)
    
    lp_cur <- logtarget(eta, z)
    lp_prp <- logtarget(eta_p, z_p)
    if (log(runif(1)) < (lp_prp - lp_cur)) { eta <- eta_p; z <- z_p }
    
    samp[i, ] <- c(exp(eta), exp(z) - 1)
    
    # light adaptation
    if (i %% adapt_every == 0 && i <= burn) {
      win <- pmax(2L, i - adapt_every + 1L):i
      # heuristic scaling from recent window
      s_eta <- sd(log(samp[win, 1])); s_z <- sd(log1p(samp[win, 2]))
      if (is.finite(s_eta) && s_eta > 0) sd_eta <- 0.8*sd_eta + 0.2*s_eta
      if (is.finite(s_z)   && s_z   > 0) sd_z   <- 0.8*sd_z   + 0.2*s_z
      sd_eta <- ifelse(is.finite(sd_eta) && sd_eta > 1e-3, sd_eta, 0.5)
      sd_z   <- ifelse(is.finite(sd_z)   && sd_z   > 1e-3, sd_z,   0.2)
    }
  }
  samp[(burn+1):iter, , drop=FALSE]
}

# =========================
# Prior calibration for GPD: solve nu_sigma so E[1 - F(T | xi, sigma)] = epsilon
# =========================
gpd_prior_mean_exceed <- function(nu_xi, nu_sigma, T, nsim=20000L) {
  xi    <- rhalfnorm(nsim, nu_xi)
  sigma <- rhalfnorm(nsim, nu_sigma)
  val <- ifelse(xi < 1e-8, exp(-T / sigma), (1 + xi * T / sigma)^(-1/xi))
  mean(val[is.finite(val)])
}

solve_nu_sigma <- function(nu_xi, T, epsilon, nsim=30000L, bracket=c(T/100, 100*T)) {
  f <- function(nu_sigma) gpd_prior_mean_exceed(nu_xi, nu_sigma, T, nsim) - epsilon
  lo <- bracket[1]; hi <- bracket[2]
  flo <- f(lo); fhi <- f(hi); tries <- 0
  while (sign(flo) == sign(fhi) && tries < 10) {
    lo <- lo/2; hi <- hi*2; flo <- f(lo); fhi <- f(hi); tries <- tries + 1
  }
  if (sign(flo) == sign(fhi)) stop("Could not bracket root for nu_sigma.")
  uniroot(f, c(lo, hi), tol=1e-3)$root
}

# =========================
# 1) Sample size calculators (PPD rule, planning with all-L censored)
# =========================

# Helper: linear N search (N_start, then N += N_step until success or N_max)
.linear_N_search <- function(check_fun, N_start = 1L, N_step = 5L, N_max = 20000L) {
  N  <- max(1L, as.integer(N_start))
  st <- max(1L, as.integer(N_step))
  last <- NULL
  while (N <= N_max) {
    res <- check_fun(N)
    last <- res
    if (isTRUE(res$ok)) return(list(N = N, payload = res))
    N <- N + st
  }
  stop(sprintf("Reached N_max (%d) without meeting PPD criterion. Last PPD=%.4g", N_max, last$ppd %||% NA_real_))
}


####### NEW / REPLACED: Prior calibration to make Pr(success) = 0.5 #######

# --- Exponential: set median(lambda) = -log(epsilon)/T ---
# If lambda ~ Gamma(alpha, rate=beta), then median(lambda) = qgamma(0.5, shape=alpha, rate=beta).
# Using scaling, median(Gamma(alpha, rate=beta)) = qgamma(0.5, shape=alpha, rate=1) / beta.
# Solve for beta to hit target median.
exp_prior_beta_from_alpha_median <- function(alpha, T, epsilon) {
  lambda_med_target <- -log(epsilon) / T
  if (!is.finite(lambda_med_target) || lambda_med_target <= 0)
    stop("Invalid (T, epsilon) for exponential median calibration.")
  qgamma(0.5, shape = alpha, rate = 1) / lambda_med_target
}

# --- GPD: sigma_star for a given xi on the decision boundary Pr(X>T|xi,sigma)=epsilon ---
gpd_sigma_star_from_xi <- function(xi, T, epsilon) {
  # continuous at xi -> 0 (exponential limit)
  out <- ifelse(
    xi < 1e-10,
    T / (-log(epsilon)),
    xi * T / (epsilon^(-xi) - 1)
  )
  out
}

# --- GPD: Choose nu_sigma so that Pr( sigma < sigma_star(xi) ) = 0.5 under xi~HN(nu_xi), sigma~HN(nu_sigma) ---
# Since sigma ~ |N(0, nu_sigma^2)|, Pr(sigma < s*) = 2*Phi(s*/nu_sigma) - 1.  We average this over xi.
solve_nu_sigma_median <- function(nu_xi, T, epsilon, nsim = 40000L) {
  if (nu_xi <= 0) stop("nu_xi must be > 0")
  if (!(epsilon > 0 && epsilon < 1)) stop("epsilon must be in (0,1)")
  xi_draws <- abs(rnorm(nsim, 0, nu_xi))
  s_star   <- gpd_sigma_star_from_xi(xi_draws, T, epsilon)
  s_star   <- pmax(s_star, .Machine$double.eps)  # guard
  
  # Root function: mean( 2*Phi(s_star/nu_sigma) - 1 ) - 0.5 = 0
  f <- function(nu_sigma) {
    if (nu_sigma <= 0) return(1)  # force positive
    mean(2 * pnorm(s_star / nu_sigma) - 1) - 0.5
  }
  
  # Bracket: base on s_star quantiles (robust & scale-aware)
  qs <- quantile(s_star, c(0.25, 0.75), na.rm = TRUE)
  lo <- max(qs[[1]] / 4, .Machine$double.eps)
  hi <- qs[[2]] * 4
  
  # Widen if needed
  flo <- f(lo); fhi <- f(hi); tries <- 0
  while (sign(flo) == sign(fhi) && tries < 10) {
    lo <- lo / 2; hi <- hi * 2
    flo <- f(lo);  fhi <- f(hi); tries <- tries + 1
  }
  if (sign(flo) == sign(fhi)) stop("Could not bracket root for nu_sigma (median criterion).")
  
  uniroot(f, c(lo, hi), tol = 1e-3)$root
}

####### NEW / REPLACED: Linear N search helper (fixed increment) #######
.linear_N_search <- function(check_fun, N_start = 1L, N_step = 5L, N_max = 20000L) {
  N  <- max(1L, as.integer(N_start))
  st <- max(1L, as.integer(N_step))
  last <- NULL
  while (N <= N_max) {
    res <- check_fun(N); last <- res
    if (isTRUE(res$ok)) return(list(N = N, payload = res))
    N <- N + st
  }
  stop(sprintf("Reached N_max (%d) without meeting PPD criterion. Last PPD=%.6g",
               N_max, if (!is.null(last)) last$ppd else NA_real_))
}

####### REPLACED: planners use the new 0.5-prior calibrations + linear stepping #######

# Exponential planner (PPD mean gate) with median-calibrated prior
cci_exp_sample_size <- function(C, L, T, epsilon, alpha,
                                iter = 15000, burn = 5000, prop_sd = 1000,
                                N_start = 5L, N_step = 5L, N_max = 10000L,
                                init_lambda = NULL) {
  beta <- exp_prior_beta_from_alpha_median(alpha, T, epsilon)/5
  if (is.null(init_lambda)) init_lambda <- 12716898#alpha / beta  # robust default
  
  check_N <- function(N) {
    x <- rep(L, N); cens <- rep(1, N)
    post <- exp_mh(x, cens, alpha, beta, init = init_lambda, iter = iter, burn = burn, prop_sd = prop_sd)
    ppd  <- mean(1 - pexp(T, rate = post))  # PPD mean exceedance
    list(ok = (ppd < epsilon), ppd = ppd)
  }
  
  res <- .linear_N_search(check_N, N_start = N_start, N_step = N_step, N_max = N_max)
  list(N = res$N, ppd = res$payload$ppd, alpha = alpha, beta = beta)
}

# GPD planner (PPD mean gate) with median-calibrated prior (Pr(success)=0.5)
cci_gpd_sample_size <- function(C, L, T, epsilon, nu_xi,
                                iter = 20000, burn = 10000,
                                sd_eta = 0.5, sd_z = 0.2,    # proposal sds for gpd_mh (log sigma, log1p xi)
                                N_start = 5L, N_step = 5L, N_max = 20000L,
                                prior_nsim = 40000L,
                                init_sigma = NULL, init_xi = NULL) {
  nu_sigma <- 1/12*solve_nu_sigma_median(nu_xi, T, epsilon, nsim = prior_nsim)
  
  # boundary-centered, “typical” starts
  if (is.null(init_xi))    init_xi    <- qnorm(0.75) * nu_xi         # median of half-normal = 0.67449 * sd
  if (is.null(init_sigma)) {
    sigma_star <- gpd_sigma_star_from_xi(init_xi, T, epsilon)
    init_sigma <- max(sigma_star, L, C, T/20, .Machine$double.eps)
  }
  
  check_N <- function(N) {
    x <- rep(L, N); cens <- rep(1, N)
    post <- gpd_mh(x, cens, nu_xi, nu_sigma,
                   init_sigma = init_sigma, init_xi = init_xi,
                   iter = iter, burn = burn, sd_eta = sd_eta, sd_z = sd_z)
    sigma_s <- post[,1]; xi_s <- post[,2]
    ppd <- mean(ifelse(xi_s < 1e-8, exp(-T / sigma_s),
                       (1 + xi_s * T / sigma_s)^(-1/xi_s)))
    list(ok = (ppd < epsilon), ppd = ppd)
  }
  
  res <- .linear_N_search(check_N, N_start = N_start, N_step = N_step, N_max = N_max)
  list(N = res$N, ppd = res$payload$ppd, nu_xi = nu_xi, nu_sigma = nu_sigma)
}

####### REPLACED: analyses use the new 0.5-prior calibrations #######


# =========================
# 2) PASS/FAIL analysis on current data (PPD rule)
# =========================

# PASS if E[Pr(X>T | θ) | data] < ε
cci_exp_analyze <- function(data, cens, T, epsilon, alpha,
                            iter=15000, burn=5000, prop_sd=1000) {
  beta <- exp_prior_beta_from_alpha(alpha, T, epsilon)
  post <- exp_mh(data, cens, alpha, beta, init=alpha/beta, iter=iter, burn=burn, prop_sd=prop_sd)
  ppd_exceed <- mean(1 - pexp(T, rate=post))
  decision   <- ifelse(ppd_exceed < epsilon, "PASS", "FAIL")
  pr_fraction <- mean((1 - pexp(T, rate=post)) < epsilon)  # optional diagnostic
  list(decision=decision, ppd_exceed=ppd_exceed, pr_fraction=pr_fraction,
       alpha=alpha, beta=beta, post=post)
}

cci_gpd_analyze <- function(data, cens, T, epsilon, nu_xi,
                            iter=20000, burn=10000, sd_eta=0.5, sd_z=0.2,
                            prior_nsim=30000L) {
  nu_sigma <- solve_nu_sigma(nu_xi, T, epsilon, nsim=prior_nsim)
  xi0    <- max(0.03, nu_xi * sqrt(2/pi))
  sigma0 <- xi0 * T / (epsilon^(-xi0) - 1)
  
  post <- gpd_mh(data, cens, nu_xi, nu_sigma,
                 init_sigma=sigma0, init_xi=xi0,
                 iter=iter, burn=burn, sd_eta=sd_eta, sd_z=sd_z)
  sigma_s <- post[,1]; xi_s <- post[,2]
  ppd_exceed <- mean(ifelse(xi_s < 1e-8, exp(-T / sigma_s),
                            (1 + xi_s * T / sigma_s)^(-1/xi_s)))
  decision <- ifelse(ppd_exceed < epsilon, "PASS", "FAIL")
  pr_fraction <- mean(ifelse(xi_s < 1e-8, exp(-T / sigma_s),
                             (1 + xi_s * T / sigma_s)^(-1/xi_s)) < epsilon)
  list(decision=decision, ppd_exceed=ppd_exceed, pr_fraction=pr_fraction,
       nu_xi=nu_xi, nu_sigma=nu_sigma, post=post)
}
