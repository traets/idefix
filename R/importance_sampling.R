

#' log Posterior
#'
#' Calculates the logposterior with a normal prior density
#' @param par Numeric vector with parametervalues.
#' @param des A design matrix in which each row is a profile.
#' @param y A binary response vector.
#' @param n.alts The number of alternatives in each choice set.
#' @param prior.mean vector containing the prior mean.
#' @param prior.covar matrix containing the prior covariance.
#' @return the logposterior probability
LogPost <- function(par, prior.mean, prior.covar, des,  n.alts, y) {
  #calcultate utility alternatives
  u <- t(t(des) * par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  #calculate probability alternatives
  expu <- exp(u)
  p <- expu / rep(rowsum(expu, rep(seq(1, nrow(des)/n.alts, 1), each = n.alts)), each = n.alts)
  #loglikelihood
  ll <- sum(y * log(p))
  #logprior
  logprior2 <- -0.5 * (par - t(prior.mean)) %*% solve(prior.covar) %*% (as.matrix(par) - prior.mean)
  #logposterior 
  logpost <- (length(prior.mean) / 2 * log(2 * pi) - 0.5 * log(det(prior.covar))) + logprior2 + ll
  return(logpost)
}


#' Hessian
#'
#' @param par Numeric vector with parametervalues.
#' @param des A design matrix in which each row is a profile.
#' @param covar The covariance matrix.
#' @param n.alts The number of alternatives in each choice set.
#' @return the hessian matrix
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

#' Likelihood function
#'
#' @param par Numeric vector with parametervalues.
#' @param des A design matrix in which each row is a profile.
#' @param n.alts The number of alternatives in each choice set
#' @param y A binary response vector.
#' @return the likelihood
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


#' Density multivariate t-distribution
#'
#' @param par Numeric vector with parametervalues.
#' @param g.mean vector containing the mean of the multivariate t-distribution.
#' @param g.covar covariance matrix of the multivariate t-distribution.
#' @return density
Gdens <- function(par, g.mean, g.covar) {
  df <- length(g.mean)
  n <- length(par)
  dif <- g.mean - par
  invcov <- solve(g.covar)
  differ <- as.numeric(t(dif) %*% invcov %*% dif)
  iMVSTd <- 1 / (det(g.covar)^(0.5)) * (1 + ((1 / df) * differ))^(-(df + length(par)) / 2)
}


#' Importance sampling
#'
#' This functions samples from an imortance density (multivariate t-distribution),
#' and gives weightes to the samples according to the posterior distribution. The prior is
#' a normal distribution.
#' @param prior.mean Numeric vector which is the mean of the multivariate normal distribution (prior).
#' @param prior.covar A matrix, the covariance matrix of the multivariate normal distribution (prior).
#' @param des A design matrix in which each row is a profile.
#' @param n.alts Numeric value indicating the number of alternatives per choice set.
#' @param y A binary response vector.
#' @param m Numeric value. Number of samples = base^m.
#' @param b Numeric value indicating the base (default = 2).
#' @return A list containing samples, their associated weights, the maximum
#'   likelihood estimates and the estimated covariance matrix.
#' @export
ImpSamp <- function(prior.mean, prior.covar, des,  n.alts, y, m, b=2, ...) {
  #error handling
  if (length(prior.mean != ncol(prior.covar))) {
    stop("Different number of parameters in prior mean and prior covarance matrix.")
  }
  if (nrow(des) != length(y)) {
    stop("Response vector length differs from the expected based on design")
  }
  if(det(prior.covar) == 0) {
    stop("prior covariance matrix is not invertible")
  }
  #prior cte
  kPrior <- (2 * pi)^(-length(prior.mean) / 2) * (det(prior.covar))^(-0.5)
  #estimate importance distributions (g) mean by finding maximum logposterior
  maxest <- maxLik::maxNR(logPost, start = prior.mean, prior.mean = prior.mean, prior.covar = prior.covar,
                          des = des, y = y, n.alts = n.alts)$estimate
  #draws from importance density
  H <- hessian(par = maxest, des = des, covar = prior.covar, n.alts = n.alts)
  g.covar <- -solve(H)
  g.draws <- lattice_mvt(mean = maxest, cvar = g.covar, df = length(maxest), m = m, ...)
  #vectors
  prior <- LK <- dens.g <- weights <- numeric(nrow(g.draws))
  # for every sample calculate prior density, likelihood and importance density 
  for (r in 1:nrow(g.draws)) {
    #prior
    prior[r] <- kPrior * exp(-0.5 * (g.draws[r, ] - prior.mean) %*% solve(prior.covar) %*% as.matrix(g.draws[r, ] - prior.mean))
    #likelihood
    LK[r] <- Lik(par = g.draws[r, ], des = des, y = y, n.alts = n.alts)
    #density of g
    dens.g[r]<-Gdens(par = g.draws[r, ], g.mean = maxest, g.covar = g.covar)
  }
  #compute the weights of samples
  w <- LK * prior / dens.g   # posterior / importance density  
  w <- w / sum(w)
  #return
  return(list(g.draws, w, maxest, g.covar))
}





