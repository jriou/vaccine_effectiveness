#' Statistical model to estimate the relative risk of hospitalization or death
#' among non-fully vaccinated individuals compared to fully vaccinated from 
#' individual-level surveillance data.
#'
#' @param vacc Indicator of vaccination status for each individual.
#' @param coverage_vacc Vaccine coverage in the population sharing the same
#' characteristics (age, sex, time, location...) with the individual.
#' @return Maximum likelihood estimate of the relative risk together with 
#' 95% confidence interval.

estimate_rr <- function(vacc, coverage_vacc) {
  # compute size
  n <- length(vacc)
  # correct coverage_vacc for extreme values
  coverage_vacc[coverage_vacc == 0] <- .001
  coverage_vacc[coverage_vacc == 1] <- .999
  # likelihood function
  loglik_fn <- function(beta) {
    lp <- exp(beta)
    p <- coverage_vacc / (coverage_vacc + lp * (1 - coverage_vacc))
    loglik <- dbinom(vacc, prob = p, size = 1, log = TRUE) # bernoulli is equivalent to binom with n=1
    target <- -sum(loglik)
    return(target)
  }
  # optimize likelihood function
  fit <- optim(
    par = log(1), # starting value
    loglik_fn,
    method = "Brent",
    hessian = TRUE,
    lower = -100, upper = 100
  )
  # get standard error from hessian
  fisher_info <- solve(fit$hessian)
  se <- sqrt(diag(fisher_info))
  upper <- exp(fit$par + qnorm(.975) * se)
  lower <- exp(fit$par - qnorm(.975) * se)
  return(list(
    rr = exp(fit$par),
    rr_se = se,
    rr_lower = lower,
    rr_upper = upper
  ))
}
