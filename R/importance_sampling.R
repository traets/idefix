

#' Importance sampling MNL
#' 
#' This function samples from the posterior distribution using importance 
#' sampling, assuming a multivariate (truncated) normal prior distribution and a
#' MNL likelihood.
#' @inheritParams SeqDB
#' @param n numeric value. The number of draws. 
#' @param prior.mean Numeric vector indicating the mean of the multivariate
#'   normal distribution (prior).
#' @param y A binary response vector. \code{\link{RespondMNL}} can be used to
#'   simulate respons data.
#' @param lower Numeric vector. Vector of lower truncation points, the default
#'   is \code{rep(-Inf, length(prior.mean))}.
#' @param upper Numeric vector. Vector of upper truncation points, the default
#'   is \code{rep(Inf, length(prior.mean))}.
#' @return \item{sample}{Numeric vector with the (unweigthted) draws from the
#' posterior distribution.} \item{weights}{Numeric vector with the associated
#' weights of the draws.} \item{max}{Numeric vector with the estimated
#' mode of the posterior distribution.} \item{covar}{Matrix representing the
#' estimated variance covariance matrix.}
#' @examples 
#' ## Example 1: sample from posterior 
#' # choice design  
#' design <- example_design 
#' # Respons.
#' truePar <- c(0.7, 0.6, 0.5, -0.5, -0.7, 1.7) # some values
#' resp <- RespondMNL(par = truePar, des = design, n.alts = 2)
#' #prior
#' pm <- rep(0, ncol(design)) # mean vector 
#' pc <- diag(2, ncol(design)) # covariance matrix 
#' # draws from posterior.
#' ImpsampMNL(n = 64, prior.mean =  pm, prior.covar = pc, 
#'            des = design, n.alts = 2, y = resp)

#' ## example 2:  sample from posterior with truncated prior
#' # choice design. 
#' design <- example_design
#' # Respons.
#' truePar <- c(0.7, 0.6, 0.5, -0.5, -0.7, 1.7) # some values
#' resp <- RespondMNL(par = truePar, des = design, n.alts = 2)
#' # prior
#' pm <- c(0.5, 0.5, 0.5, -0.5, -0.5, 0.5) # mean vector
#' pc <- diag(2, ncol(design)) # covariance matrix
#' low <- rep(-Inf, length(pm)) # lower bound 
#' up <- c(Inf, Inf, Inf, 0, 0, Inf) # upper bound
#' # draws from posterior.
#' ImpsampMNL(n = 64, prior.mean =  pm, prior.covar = pc, des = design,
#'            n.alts = 2, y = resp, lower = low, upper = up)
#' @export
ImpsampMNL <- function(n, prior.mean, prior.covar, des, n.alts, y,
                        lower = rep(-Inf, length(prior.mean)), upper = rep(Inf, length(prior.mean))){
  # Error handling.
  if (length(prior.mean) != ncol(prior.covar)) {
    stop("different number of parameters in prior mean and prior covarance matrix")
  }
  if (nrow(des) != length(y)) {
    stop("response vector length differs from the expected based on design")
  }
  if(det(prior.covar) == 0) {
    stop("prior covariance matrix is not invertible")
  }
  if(!(all(lower < prior.mean) && all(prior.mean < upper))){
    stop("prior.mean is not between lower and upper bounds")
  }
  if(length(lower) != length(prior.mean)){
    stop('length lower and prior.mean does not match')
  }
  if(length(upper) != length(prior.mean)){
    stop('length upper and prior.mean does not match')
  }
  # mode imp dens 
  impdens <- maxLik::maxNR(LogPost, start = prior.mean, prior.mean = prior.mean, prior.covar = prior.covar,
                           lower = lower, upper = upper, des = des, y = y, n.alts = n.alts)
  maxest <- impdens$estimate
  # covar imp dens 
  g.covar <- -solve(impdens$hessian)
  # draws from imp dens 
  g.draws <- tmvtnorm::rtmvt(n = n, mean = maxest, sigma = g.covar, lower = lower, upper = upper, df = length(maxest))
  # prior dens
  prior <- tmvtnorm::dtmvnorm(g.draws, mean = prior.mean, sigma = prior.covar, lower = lower, upper = upper)
  # likelihood
  likh <-  apply(g.draws, 1, Lik, des = des, y = y, n.alts = n.alts)
  # imp dens
  g.dens <- tmvtnorm::dtmvt(g.draws, mean = maxest, sigma = g.covar, df = length(maxest))
  # weights of draws.
  w <- likh * prior / g.dens 
  w <- w / sum(w)
  # Return.
  return(list(sample = g.draws, weights = w, max = maxest, covar = g.covar))
}

# log Posterior
# 
# Calculates the logposterior with a normal prior density par Numeric vector
# with parametervalues.
# @param des A design matrix in which each row is a profile.
# @param y A binary response vector.
# @param n.alts The number of alternatives in each choice set.
# @param prior.mean vector containing the prior mean.
# @param prior.covar matrix containing the prior covariance.
# @param lower Numeric vector. Vector of lower truncation points, the default
#   is \code{rep(-Inf, length(prior.mean))}.
# @param upper Numeric vector. Vector of upper truncation points, the default
#   is \code{rep(Inf, length(prior.mean))}. 
# @return the logposterior probability
LogPost <- function(par, prior.mean, prior.covar, lower, upper,  des,  n.alts, y) {
  #calcultate utility alternatives
  u <- t(t(des) * par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  #calculate probability alternatives
  expu <- exp(u)
  p <- expu / rep(rowsum(expu, rep(seq(1, nrow(des)/n.alts, 1), each = n.alts)), each = n.alts)
  #loglikelihood
  ll <- sum(y * log(p))
  #logprior
  lprior <- tmvtnorm::dtmvnorm(par, prior.mean, prior.covar, lower, upper, log = TRUE)
  #logposterior 
  return(lprior + ll)
}


# Hessian
#
# @param par Numeric vector with parametervalues.
# @param des A design matrix in which each row is a profile.
# @param covar The covariance matrix.
# @param n.alts The number of alternatives in each choice set.
# @return the hessian matrix
Hessian <- function(par, des, covar, n.alts) {
  # utility
  des <- as.matrix(des)
  u <- des %*% diag(par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  # probability 
  expu <- exp(u)
  p <- expu / rep(rowsum(expu, rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)), each = n.alts)
  # information matrix
  info <- crossprod(des * p, des) - crossprod(rowsum(des * p, rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)))
  # hessian 
  hess <- (-info - solve(covar))
  return(hess)
}

# Likelihood function
#
# @param par Numeric vector with parametervalues.
# @param des A design matrix in which each row is a profile.
# @param n.alts The number of alternatives in each choice set
# @param y A binary response vector.
# @return the likelihood
Lik <- function(par, des, n.alts, y) {
  # utility
  des <- as.matrix(des)
  u <- t(t(des) * par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  # probability
  expu <- exp(u)
  p <- expu / rep(rowsum(expu, rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)), each = n.alts)
  # likelihood
  L <- prod(p^y)
  return(L)
}


# Density multivariate t-distribution
#
# @param par Numeric vector with parametervalues.
# @param g.mean vector containing the mean of the multivariate t-distribution.
# @param g.covar covariance matrix of the multivariate t-distribution.
# @return density
Gdens <- function(par, g.mean, g.covar) {
  df <- length(g.mean)
  n <- length(par)
  dif <- g.mean - par
  invcov <- solve(g.covar)
  differ <- as.numeric(t(dif) %*% invcov %*% dif)
  iMVSTd <- 1 / (det(g.covar)^(0.5)) * (1 + ((1 / df) * differ))^(-(df + length(par)) / 2)
}






