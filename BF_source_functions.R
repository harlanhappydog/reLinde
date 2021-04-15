# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Load necessary script with functions from Gronau, Ly, & Wagenmakers,
# which is available at https://osf.io/bsp6z/.
#source("informedTtest_functions.R")



ciPlusMedian_t <- function(t,
                           n1,
                           n2 = NULL,
                           independentSamples = FALSE,
                           prior.location,
                           prior.scale,
                           prior.df,
                           ci = .95,
                           type = "two-sided",
                           tol = 0.0001,
                           max.iter = 100,
                           rel.tol = .Machine$double.eps^0.25) {
  
  lower <- (1 - ci)/2
  upper <- ci + (1 - ci)/2
  med <- .5
  
  postAreaSmaller0 <- cdf_t(x = 0, t = t, n1 = n1, n2 = n2,
                            independentSamples = independentSamples,
                            prior.location = prior.location,
                            prior.scale = prior.scale,
                            prior.df = prior.df,
                            rel.tol = rel.tol)
  
  if (type == "plus-sided") {
    
    lower <- postAreaSmaller0 + (1 - postAreaSmaller0)*lower
    upper <- postAreaSmaller0 + (1 - postAreaSmaller0)*upper
    med <- postAreaSmaller0 + (1 - postAreaSmaller0)*med
    
  } else if (type == "min-sided") {
    
    lower <- postAreaSmaller0*lower
    upper <- postAreaSmaller0*upper
    med <- postAreaSmaller0*med
    
  }
  
  ciLower <- quantile_t(lower, t = t, n1 = n1, n2 = n2,
                        independentSamples = independentSamples,
                        prior.location = prior.location,
                        prior.scale = prior.scale,
                        prior.df = prior.df,
                        rel.tol = rel.tol)
  ciUpper <- quantile_t(upper, t = t, n1 = n1, n2 = n2,
                        independentSamples = independentSamples,
                        prior.location = prior.location,
                        prior.scale = prior.scale,
                        prior.df = prior.df,
                        rel.tol = rel.tol)
  median <- quantile_t(med, t = t, n1 = n1, n2 = n2,
                       independentSamples = independentSamples,
                       prior.location = prior.location,
                       prior.scale = prior.scale,
                       prior.df = prior.df,
                       rel.tol = rel.tol)
  
  return(list(ciLower = ciLower, median = median, ciUpper = ciUpper))
  
}


integrand_t <- function(delta, t, n, nu, mu.delta, gamma, kappa) {
  
  suppressWarnings(
    dt(x = t, df = nu, ncp = sqrt(n) * delta) *
    1 / gamma * dt( (delta - mu.delta) / gamma, df = kappa)
  )
  
}

posterior_t <- function(delta,
                        t,
                        n1,
                        n2 = NULL,
                        independentSamples = FALSE,
                        prior.location,
                        prior.scale,
                        prior.df,
                        rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, n1 * n2 / (n1 + n2), n1)
  nu <- ifelse(independentSamples, n1 + n2 - 2, n1 - 1)
  
  mu.delta <- prior.location
  gamma <- prior.scale
  kappa <- prior.df
  
  numerator <- suppressWarnings(
    dt(x = t, df = nu, ncp = sqrt(neff) * delta) *
    1 / gamma * dt( (delta - mu.delta) / gamma, df = kappa)
  )
  
  denominator <- integrate(integrand_t,
                           lower = -Inf, upper = Inf,
                           t = t, n = neff, nu = nu,
                           mu.delta = mu.delta,
                           gamma = gamma,
                           kappa = kappa,
                           rel.tol = rel.tol)$value
  
  out <- numerator / denominator
  out[is.na(out)] <- 0
  
  return(out)
  
}

cdf_t <- function(x,
                  t,
                  n1,
                  n2 = NULL,
                  independentSamples = FALSE,
                  prior.location,
                  prior.scale,
                  prior.df,
                  rel.tol = .Machine$double.eps^0.25) {
  
  out <- integrate(posterior_t,
                   lower = -Inf, upper = x,
                   t = t, n1 = n1, n2 = n2,
                   independentSamples = independentSamples,
                   prior.location = prior.location,
                   prior.scale = prior.scale,
                   prior.df = prior.df,
                   rel.tol = rel.tol)$value
  
  # catch numerical errors
  if (out > 1 & out < 1.001) {
    out <- 1
    warning(
      "Numerical integration yields a CDF value slightly larger than 1. The CDF value has been replaced by 1.",
      call. = FALSE
    )
  }
  
  return(out)
  
}

quantile_t <- function(q,
                       t,
                       n1,
                       n2 = NULL,
                       independentSamples = FALSE,
                       prior.location,
                       prior.scale,
                       prior.df,
                       tol = 0.0001,
                       max.iter = 100,
                       rel.tol = .Machine$double.eps^0.25) {
  
  # compute quantiles via Newton-Raphson method
  
  x.cur <- Inf
  # get reasonable starting value
  delta <- seq(-2, 2, length.out = 400)
  dens <- posterior_t(delta, t = t, n1 = n1, n2 = n2,
                      independentSamples = independentSamples,
                      prior.location = prior.location,
                      prior.scale = prior.scale,
                      prior.df = prior.df)
  x.new <- delta[which.max(dens)]
  i <- 1
  
  while (abs(x.cur - x.new) > tol && i < max.iter) {
    
    x.cur <- x.new
    x.new <- x.cur - (cdf_t(x.cur, t = t, n1 = n1, n2 = n2,
                            independentSamples = independentSamples,
                            prior.location = prior.location,
                            prior.scale = prior.scale,
                            prior.df = prior.df,
                            rel.tol = rel.tol) - q)/
      posterior_t(x.cur, t = t, n1 = n1, n2 = n2,
                  independentSamples = independentSamples,
                  prior.location = prior.location,
                  prior.scale = prior.scale,
                  prior.df = prior.df)
    i <- i + 1
    
  }
  
  return(x.new)
  
}

