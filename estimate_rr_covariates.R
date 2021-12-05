#' Statistical model to estimate the relative risk of hospitalization or death
#' among non-fully vaccinated individuals compared to fully vaccinated from 
#' individual-level surveillance data.
#'
#' @param vacc Indicator of vaccination status for each individual.
#' @param coverage_vacc Vaccine coverage in the population sharing the same
#' characteristics (age, sex, time, location...) with the individual.
#' @param cov Matrix of covariates.
#' @return Maximum likelihood estimate of the relative risk together with standard error 
#' and 95% confidence interval for the reference group, and change in relative
#' risk with 95% confidence interval by covariate compared to the reference group (corresponding
#' to a matrix of covariates with all zeros). 
#' The fisher information matrix is also provided in order to compute direct relative 
#' risks for each covariate.

estimate_rr_covariates <- function(vacc, coverage_vacc, cov) {
  # compute size
  n <- length(vacc)
  # compute matrix of covariates
  X <- as.matrix(cov)
  ncov <- ncol(X)
  # correct coverage_vacc for extreme values
  coverage_vacc[coverage_vacc == 0] <- .001
  coverage_vacc[coverage_vacc == 1] <- .999
  # likelihood function
  loglik_fn <- function(x) {
    beta <- x[1]
    gamma <- x[2:(2 + ncov - 1)]
    lp <- exp(beta + X %*% gamma)
    p <- coverage_vacc / (coverage_vacc + lp * (1 - coverage_vacc))
    loglik <- dbinom(vacc, prob = p, size = 1, log = TRUE) # bernoulli is equivalent to binom with n=1
    target <- -sum(loglik)
    return(target)
  }
  # optimize likelihood function
  fit <- optim(
    par = c(1, rep(1, ncov)), # starting value
    loglik_fn,
    method = "BFGS",
    hessian = TRUE
  )
  # get standard error from hessian
  fisher_info <- solve(fit$hessian)
  se <- sqrt(diag(fisher_info))
  # extract parameter estimates
  est_rr_ref <- fit$par[1]
  se_rr_ref <- se[1]
  est_rr_group <- fit$par[2:(2 + ncov - 1)]
  se_rr_group <- se[2:(2 + ncov - 1)]
  
  return(list(
    est_rr_ref = exp(est_rr_ref),
    se_rr_ref = se_rr_ref,
    lower_rr_ref = exp(est_rr_ref - qnorm(.975) * se_rr_ref),
    upper_rr_ref = exp(est_rr_ref + qnorm(.975) * se_rr_ref),
    est_rr_group = exp(est_rr_group),
    se_rr_group = se_rr_group,
    lower_rr_group = exp(est_rr_group - qnorm(.975) * se_rr_group),
    upper_rr_group = exp(est_rr_group + qnorm(.975) * se_rr_group),
    fisher_info = fisher_info
  ))
}

