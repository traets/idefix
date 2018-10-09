

#' Importance sampling MNL
#' 
#' This function samples from the posterior distribution using importance 
#' sampling, assuming a multivariate (truncated) normal prior distribution and a
#' MNL likelihood.
#' 
#' For the proposal distribution a t-distribution with degrees of freedom equal 
#' to the number of parameters is used. The mode is estimated using 
#' \code{\link[stats]{optim}}, and the covariance matrix is calculated as the negative 
#' inverse of the generalized Fisher information matrix. See reference for more
#' information.
#' 
#' From this distribution a lattice grid of draws is generated.
#' 
#' If truncation is present, incorrect draws are rejected and new ones are
#' generated untill \code{n} is reached. The covariance matrix is in this case
#' still calculated as if no truncation was present.
#' 
#' @inheritParams SeqDB
#' @param n.draws numeric value indicating the number of draws. 
#' @param prior.mean Numeric vector indicating the mean of the multivariate
#'   normal distribution (prior).
#' @param y A binary response vector. \code{\link{RespondMNL}} can be used to
#'   simulate respons data.
#' @param lower Numeric vector of lower truncation points, the default
#'   is \code{NULL}.
#' @param upper Numeric vector of upper truncation points, the default
#'   is \code{NULL}.
#' @return \item{sample}{Numeric vector with the (unweigthted) draws from the
#' posterior distribution.} \item{weights}{Numeric vector with the associated
#' weights of the draws.} \item{max}{Numeric vector with the estimated
#' mode of the posterior distribution.} \item{covar}{Matrix representing the
#' estimated variance covariance matrix.}
#' @references \insertRef{ju}{idefix}
#' @examples 
#' ## Example 1: sample from posterior, no constraints, no alternative specific constants 
#' # choice design  
#' design <- example_design
#' # Respons.
#' truePar <- c(0.7, 0.6, 0.5, -0.5, -0.7, 1.7) # some values
#' set.seed(123)
#' resp <- RespondMNL(par = truePar, des = design, n.alts = 2)
#' #prior
#' pm <- c(1, 1, 1, -1, -1, 1) # mean vector 
#' pc <- diag(1, ncol(design)) # covariance matrix 
#' # draws from posterior.
#' ImpsampMNL(n.draws = 100, prior.mean =  pm, prior.covar = pc,
#'            des = design, n.alts = 2, y = resp)
#' 
#' ## example 2:  sample from posterior with constraints 
#' # and alternative specific constants
#' # choice design. 
#' design <- example_design2
#' # Respons.
#' truePar <- c(0.2, 0.8, 0.7, 0.6, 0.5, -0.5, -0.7, 1.7) # some values
#' set.seed(123)
#' resp <- RespondMNL(par = truePar, des = design, n.alts = 3)
#' # prior
#' pm <- c(1, 1, 1, 1, 1, -1, -1, 1) # mean vector 
#' pc <- diag(1, ncol(design)) # covariance matrix
#' low = c(-Inf, -Inf, 0, 0, 0, -Inf, -Inf, 0)
#' up = c(Inf, Inf, Inf, Inf, Inf, 0, 0, Inf)
#' # draws from posterior.
#' ImpsampMNL(n.draws = 100, prior.mean =  pm, prior.covar = pc, des = design,
#'            n.alts = 3, y = resp, lower = low, upper = up, alt.cte = c(1, 1, 0))
#' @export
ImpsampMNL <- function(n.draws, prior.mean, prior.covar, des, n.alts, y, 
                       alt.cte = NULL, lower = NULL, upper = NULL){
  if(!is.matrix(des)){
    stop("'des' should be a matrix")
  }
  if(!isTRUE(nrow(des) %% n.alts == 0)){
    stop("'n.alts' does not seem to match with the number of rows in 'des'")
  } else {
    n.sets <- nrow(des) / n.alts
  }
  n.sets <- nrow(des) / n.alts
  if (!is.null(lower)){
    if(!is.numeric(lower)){
      stop("lower' should be a numeric vector")
    }  
  } else {
    lower <- rep(-Inf, length(prior.mean))
  }
  if (!is.null(upper)){
    if(!is.numeric(upper)){
      stop("'upper' should be a numeric vector")
    }
  } else {
    upper <- rep(Inf, length(prior.mean))
  }
  # Error handling.
  if (length(prior.mean) != ncol(prior.covar)) {
    stop("different number of parameters in 'prior.mean' and 'prior.covar' matrix")
  }
  if (nrow(des) != length(y)) {
    stop("'y' length differs from the expected based on 'des'")
  }
  if(isTRUE(all.equal(det(prior.covar), 0))) {
    stop("prior covariance matrix is not invertible")
  }
  if(!(all(lower < prior.mean) && all(prior.mean < upper))){
    stop("'prior.mean' elements are not between 'lower' and 'upper' bounds")
  }
  if(length(lower) != length(prior.mean)){
    stop("length 'lower'' and 'prior.mean' does not match")
  }
  if(length(upper) != length(prior.mean)){
    stop("'length 'upper' and 'prior.mean' does not match")
  }
  if(!is.numeric(prior.mean)){
    stop("'prior.mean' should be a numeric vector")
  }
  if(!is.null(alt.cte)){
    if (length(alt.cte) != n.alts) {
      stop("'n.alts' does not match the 'alt.cte' vector")
    }
    if (!all(alt.cte %in% c(0, 1))){
      stop("'alt.cte' should only contain zero or ones.")
    }
    # alternative specific constants
    n.cte <- length(which(alt.cte == 1L))
    if (isTRUE(all.equal(n.cte, 0L))){
      alt.cte <- NULL
      cte.des <- NULL
    }
  } 
  if(!is.null(alt.cte)){
    cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
    if(!isTRUE(all.equal(cte.des, matrix(des[ , 1:n.cte], ncol = n.cte)))){
      stop("the first column(s) of 'des' are different from what is expected based on 'alt.cte'")
    }
  } else {
    cte.des <- NULL
    n.cte <- 0
  }
  # mode imp dens 
  maxest <- stats::optim(par = prior.mean, LogPost, lower2 = lower, upper2 = upper, prior.mean = prior.mean, prior.covar = prior.covar, 
                  des = des, y = y, n.alts = n.alts, lower = lower, upper = upper, 
                  method = "L-BFGS-B", hessian = FALSE)$par
  # covar imp dens 
  hess <- Hessian(par = maxest, des = des, prior.covar = prior.covar, n.alts = n.alts)
  g.covar <- -solve(hess)
  # draws from imp dens 
  g.draws <- Lattice_trunc(n = n.draws, mean = maxest, cvar = g.covar, lower = lower, upper = upper, df = length(maxest))
  # prior dens
  prior <- tmvtnorm::dtmvnorm(g.draws, mean = prior.mean, sigma = prior.covar, lower = lower, upper = upper)
  # likelihood
  likh <-  apply(g.draws, 1, Lik, des = des, y = y, n.alts = n.alts)
  # imp dens
  g.dens <- tmvtnorm::dtmvt(g.draws, mean = maxest, sigma = g.covar, df = length(maxest))
  # weights of draws.
  w <- likh * prior / g.dens 
  w <- w / sum(w)
  #colnames draws 
  if(!is.null(cte.des)){
    des.names <- Rcnames(n.sets = n.sets, n.alts = n.alts, alt.cte = alt.cte)
    g.draws <- list(as.matrix(g.draws[ ,1:n.cte], ncol = n.cte), 
                    as.matrix(g.draws[ ,(n.cte + 1) : ncol(des)], ncol = (ncol(des) - n.cte)))
    colnames(g.draws[[1]]) <- des.names[[2]]
    colnames(g.draws[[2]]) <- colnames(des)[(n.cte + 1) : ncol(des)]
  } else {
    colnames(g.draws) <-  colnames(des)[(n.cte + 1) : ncol(des)]
  }
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
LogPost <- function(par, prior.mean, prior.covar, lower2, upper2,  des,  n.alts, y) {
  #calcultate utility alternatives
  u <- t(t(des) * par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  #calculate probability alternatives
  expu <- exp(u)
  p <- expu / rep(rowsum(expu, rep(seq(1, nrow(des)/n.alts, 1), each = n.alts)), each = n.alts)
  #loglikelihood
  ll <- sum(y * log(p))
  #logprior
  lprior <- tmvtnorm::dtmvnorm(par, prior.mean, prior.covar, lower2, upper2, log = TRUE)
  #logposterior 
  post <- lprior + ll
  return(-post)
}


# Hessian
#
# @param par Numeric vector with parametervalues.
# @param des A design matrix in which each row is a profile.
# @param covar The covariance matrix.
# @param n.alts The number of alternatives in each choice set.
# @return the hessian matrix
Hessian <- function(par, des, prior.covar, n.alts) {
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
  hess <- (-info - solve(prior.covar))
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


# # Density multivariate t-distribution
# #
# # @param par Numeric vector with parametervalues.
# # @param g.mean vector containing the mean of the multivariate t-distribution.
# # @param g.covar covariance matrix of the multivariate t-distribution.
# # @return density
# Gdens <- function(par, g.mean, g.covar) {
#   df <- length(g.mean)
#   n <- length(par)
#   dif <- g.mean - par
#   invcov <- solve(g.covar)
#   differ <- as.numeric(t(dif) %*% invcov %*% dif)
#   iMVSTd <- 1 / (det(g.covar)^(0.5)) * (1 + ((1 / df) * differ))^(-(df + length(par)) / 2)
# }






