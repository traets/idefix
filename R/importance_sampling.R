

#' Importance sampling MNL
#' 
#' This function samples from the posterior distribution using importance 
#' sampling, assuming a multivariate normal prior distribution and a MNL likelihood.
#' @inheritParams SeqDB
#' @param prior.mean Numeric vector indicating the mean of the multivariate
#'   normal distribution (prior).
#' @param y A binary response vector. \code{\link{RespondMNL}} can be used to
#'   simulate respons data.
#' @param m Numeric value. Number of samples = \code{base^m}.
#' @param b Numeric value indicating the base. The default = 2.
#' @return \item{samples}{Numeric vector with the (unweigthted) samples from the
#' posterior distribution.} \item{weights}{Numeric vector with the associated
#' weights of the samples.} \item{max}{Numeric vector with the estimated
#' mode of the posterior distribution.} \item{covar}{Matrix representing the
#' estimated variance covariance matrix.}
#' @examples 
#' # Importance sampling MNL 
#' pm <- c(0.8, 0.3, 0.2, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(pm)) # Prior variance
#' cs <- Profiles(lvls = c(3, 3), coding = c("E", "E"))
#' ps <- MASS::mvrnorm(n = 10, mu = pm, Sigma = pc) # 10 Samples.
#' # Efficient design. 
#' design <- Modfed(cand.set = cs, n.sets = 8, n.alts = 2, alt.cte = c(1,0), par.draws = ps)$design
#' # Respons.
#' resp <- RespondMNL(par = c(0.7, 0.6, 0.5, -0.5, -0.7), des = design, n.alts = 2)
#' # Parameters samples from posterior.
#' ImpsampMNL(prior.mean =  pm, prior.covar = pc, des = design, n.alts = 2, y = resp, m = 6)
#'
#' # Importance sampling MNL
#' pm <- c(0.3, 0.2, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(pm)) # Prior variance
#' cs <- Profiles(lvls = c(3, 3, 2), coding = c("D", "C", "D"), c.lvls = list(c(2,4,6)))
#' ac <- c(0, 0) # No alternative specific constants. 
#' ps <- MASS::mvrnorm(n = 10, mu = pm, Sigma = pc) # 10 Samples.
#' # Efficient design. 
#' design <- Modfed(cand.set = cs, n.sets = 8, n.alts = 2, alt.cte = c(0,0), par.draws = ps)$design
#' # Respons
#' resp <- RespondMNL(par = c(0.6, 0.5, -0.5, -0.7), des = design, n.alts = 2)
#' # Parameters samples from posterior.
#' ImpsampMNL(prior.mean =  pm, prior.covar = pc, des = design, n.alts = 2, y = resp, m = 6)
#' @export
ImpsampMNL <- function(prior.mean, prior.covar, des, n.alts, y, m, b = 2) {
  # Error handling.
  if (length(prior.mean) != ncol(prior.covar)) {
    stop("different number of parameters in prior mean and prior covarance matrix.")
  }
  if (nrow(des) != length(y)) {
    stop("response vector length differs from the expected based on design")
  }
  if(det(prior.covar) == 0) {
    stop("prior covariance matrix is not invertible")
  }
  # Prior cte.
  kPrior <- (2 * pi)^(-length(prior.mean) / 2) * (det(prior.covar))^(-0.5)
  # Estimate mean importance distribution by finding maximum logposterior.
  maxest <- maxLik::maxNR(LogPost, start = prior.mean, prior.mean = prior.mean, prior.covar = prior.covar,
                          des = des, y = y, n.alts = n.alts)$estimate
  # Draws from importance density.
  hes <- Hessian(par = maxest, des = des, covar = prior.covar, n.alts = n.alts)
  g.covar <- -solve(hes)
  g.draws <- Lattice_mvt(mean = maxest, cvar = g.covar, df = length(maxest), m = m)
  # Vectors.
  prior <- likh <- dens.g <- weights <- numeric(nrow(g.draws))
  # For every sample calculate prior density, likelihood and importance density.
  spvc <- solve(prior.covar)
  for (r in 1:nrow(g.draws)) {
    # Prior.
    prior[r] <- kPrior * exp(-0.5 * (g.draws[r, ] - prior.mean) %*% spvc %*% as.matrix(g.draws[r, ] - prior.mean))
    # Likelihood.
    likh[r] <- Lik(par = g.draws[r, ], des = des, y = y, n.alts = n.alts)
    # Density of g.
    dens.g[r] <- Gdens(par = g.draws[r, ], g.mean = maxest, g.covar = g.covar)
  }
  # Compute the weights of samples.
  w <- likh * prior / dens.g   # posterior / importance density  
  w <- w / sum(w)
  # Return.
  return(list(samples = g.draws, weights = w, max = maxest, covar = g.covar))
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
# @return the logposterior probability
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






