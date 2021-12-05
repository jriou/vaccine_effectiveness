#' Statistical model to estimate the relative risk of hospitalization or death
#' among individuals with vaccination status 2 and 3 compared to individuals with vaccination
#' status 1 from individual-level surveillance data.
#'
#' @param vacc3 Indicator of vaccination status for each individual (0/1 dummy variables in N*3 matrix form).
#' @param coverage_vacc3 Vaccine coverage in the population sharing the same
#' characteristics (age, sex, time, location...) with the individual in N*3 matrix form.
#' @return Maximum likelihood estimate of the relative risk and 95% confidence interval for individuals with vaccination 
#' status 2 and 3 compared to status 1.

estimate_rr_3groups <- function(vacc3, coverage_vacc3) {
  # compute size
  n <- nrow(vacc3)
  # # correct coverage for extreme values
  coverage_vacc3 = coverage_vacc3 + 1e-7
  coverage_vacc3 = coverage_vacc3 / apply(coverage_vacc3,1,sum)
  # likelihood function
  loglik_fn <- function(x) {
    # eta = exp(x[1])
    gamma21 = exp(x[1])
    gamma31 = exp(x[2])
    
    p1 = coverage_vacc3[,1] / (coverage_vacc3[,1] + gamma21 * coverage_vacc3[,2] + gamma31 * coverage_vacc3[,3])
    p2 = coverage_vacc3[,2] / (1/gamma21 * coverage_vacc3[,1] + coverage_vacc3[,2] + gamma31/gamma21 * coverage_vacc3[,3])
    p3 = 1 - p1 - p2 
    p = cbind(p1,p2,p3)
    loglik = numeric(n)
    for(i in 1:n) {
      loglik[i] <- dmultinom(vacc3[i,], prob = p[i,], size = 1, log = TRUE) # categorical distribution
    }
    target <- -sum(loglik)
    return(target)
  }
  # optimize likelihood function
  fit <- optim(
    par = c(1,1), # starting value
    loglik_fn,
    method = "L-BFGS-B",
    lower=-10,
    upper=10,
    hessian = TRUE
  )
  # get standard error from hessian
  fisher_info <- solve(fit$hessian)
  se <- sqrt(diag(fisher_info))
  upper <- exp(fit$par+ 1.96 * se)
  lower <- exp(fit$par - 1.96 * se)
  return(list(
    est_rr = exp(fit$par),
    se_rr = se,
    lower_rr = lower,
    upper_rr = upper
  ))
}
