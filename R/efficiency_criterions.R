

#' DB error
#'
#' Function to calculate the DB-error given a design, and parameter values.
#'
#' @param par.draws Numeric matrix in which each row is a draw from a multivariate parameter distribution.
#' @param des A design matrix in which each row is an alternative.
#' @param n.alts Numeric value indicating the number of alternatives per choice
#'   set.
#' @param mean A logical value indicating whether the mean (DB) error should be returned or not. Default = TRUE. 
#' @param weights A numeric vector containing weights of \code{par.draws}. The
#'   default is \code{NULL}.
#' @return Numeric value indicating the DB-error of the design given the
#'   parameter draws.
#' @examples 
#' des <- example_design
#' mu = c(-1, -1.5, -1, -1.5, 0.5, 1)
#' Sigma = diag(length(mu))
#' par.draws <- MASS::mvrnorm(100, mu = mu, Sigma = Sigma)
#' n.alts = 2
#' DBerr(par.draws = par.draws, des = des, n.alts = n.alts)
#' 
#' mu = c(-0.5, -1, -0.5, -1, 0.5, 1)
#' Sigma = diag(length(mu))
#' par.draws <- MASS::mvrnorm(100, mu = mu, Sigma = Sigma)
#' DBerr(par.draws = par.draws, des = des, n.alts = n.alts)
#' @importFrom Rdpack reprompt
#' @export
DBerr <- function(par.draws, des, n.alts, weights = NULL, mean = TRUE) {
  if(is.list(par.draws)){
    if (!isTRUE(all.equal(length(par.draws), 2))){
      stop("If 'par.draws' is a list, it should contain two components")
    }
    if(!(all(unlist(lapply(par.draws, is.matrix))))){
      stop("'par.draws' should contain two matrices")
    }
    dims <-  as.data.frame(lapply(par.draws, dim))
    if(!isTRUE(all.equal(dims[1, 1], dims[1, 2]))){ 
      stop("the number of rows in the components of 'par.draws' should be equal")
    }
    if(!identical((dims[2, 1] + dims[2, 2]), ncol(des))){ 
      stop("the sum of the number of columns in the components of 'par.draws' 
           should equal the number of columns in 'des'")
    }
    par.draws  <- do.call("cbind", par.draws)
  }
  if(!is.matrix(par.draws)){
    stop("'par.draws' should be a matrix or a list")
  }
  if(!isTRUE(all.equal(ncol(par.draws), ncol(des)))){
    stop("ncol(des) should equal ncol(par.draws)")
  }
  if(!is.matrix(des)){
    stop("'des' should be a matrix")
  }
  if (!isTRUE(nrow(des) %% n.alts == 0)) {
    stop("'n.alts' does not seem to match with the number of rows in 'des'")
  }
  
  if(!is.wholenumber(n.alts)){
    stop("'n.alts' should be an integer")
  }
  if(!is.null(weights)){
    if(!is.vector(weights)){
      stop("'weights' should be a vector")
    }
    if(!isTRUE(all.equal(length(weights), nrow(par.draws)))){
      stop("length(weights) should equal nrow(par.draws)")
    }
  } else {
    weights <- rep(1, nrow(par.draws))
  }
  if (!is.logical(mean)){
    stop("'mean' should be TRUE or FALSE") 
  }
  d.errors <- apply(par.draws, 1, Derr_ucpp, des, n.alts)
  # DB-error.
  if(isTRUE(mean)){
    error <- mean(d.errors, na.rm = TRUE)
  } else {
    error <- d.errors
  }
  return("DBerror" = error)
}

# is integer? 
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#with covariance information
#DerrC <- function(par, des, n.alts, i.cov) {
#  info.des <- InfoDes(par, des, n.alts)
#  detinfo <- det(info.des + i.cov)
#  ifelse((detinfo <= 0), return(NA), return(detinfo^(-1 / length(par))))
#}

# Sequential D-error
# 
# Function to calculate D-errors if set would be part of design.
# @inheritParams Modfed
# @param set A choice set in which each row is a profile.
# @param des A design matrix in which each row is a profile.
# @param i.cov Inverse of covariance matrix.
# @param n.par Number of parameters.
#DerrS <- function(par.draws, set, des, n.alts, i.cov, n.par) {
#  des.f <- rbind(des, set) 
#  info.d <- InfoDes(par = par.draws, des = des.f, n.alts = n.alts) 
#  d.error <- det(info.d + i.cov)^(-1 / n.par)
#  return(d.error)
#}

#for parallel 
DerrS.P <- function(par, des, n.alts, i.cov) {
  group <- rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)
  # probability
  u <- des %*% diag(par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  p <- exp(u) / rep(rowsum(exp(u), group), each = n.alts)
  # information matrix
  info <- crossprod(des * p, des) - crossprod(rowsum(des * p, group))
  d.error <- det(info + i.cov)^(-1 / length(par))
  return(d.error)
}

# Sequential DB-error
# 
# Function to calculate DB-errors for potential choice sets in combination with 
# an initial design.
# @inheritParams Modfed
# @inheritParams DerrS
# @param full.comb A matrix with on each row a possible combination of profiles.
# @param cte.des A matrix which represent the alternative specific constants. If
#   there are none it value is \code{NULL}.
# @return The DB errors of the designs in which each design is a combination 
#  of the initial design with a potential choice set.
#  DBerrS <- function(full.comb, cand.set, par.draws, des, n.alts, cte.des, i.cov, n.par, weights) {
#  # Take set.
#  set <- as.matrix(cand.set[as.numeric(full.comb), ])
#  # Add alternative specific constants if necessary
#  if (!is.null(cte.des)) {
#    set <- as.matrix(cbind(cte.des, set))
#  }
#  # For each draw calculate D-error.
#  d.errors <- apply(par.draws, 1, DerrS, set, des, n.alts, i.cov, n.par)
#  w.d.errors <- d.errors * weights
#  # DB-error. 
#  db.error <- mean(w.d.errors, na.rm = TRUE)
#  return(db.error)
#}

# function to get DB error of start designs
StartDB <- function(des, par.draws, n.alts){
  apply(par.draws, 1, Derr_ucpp, des = des,  n.alts = n.alts)
} 


#for parallel
DBerrS.P <- function(des, par.draws, n.alts, i.cov, weights) {
  # Add alternative specific constants if necessary
  # For each draw calculate D-error.
  d.errors <- apply(par.draws, 1, DerrS.P, des, n.alts, i.cov)
  w.d.errors <- d.errors * weights
  # DB-error. 
  db.error <- mean(w.d.errors, na.rm = TRUE)
  return(db.error)
}


# Fisher Information of design
# 
# Returns the Fisher Information of a design, given parameter values.
# @inheritParams Modfed
# @param par A vector containing the parameter values
# @return Fisher Information matrix.
InfoDes <- function(par, des, n.alts) {
  group <- rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)
  # probability
  u <- des %*% diag(par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  p <- exp(u) / rep(rowsum(exp(u), group), each = n.alts)
  # information matrix
  info.des <- crossprod(des * p, des) - crossprod(rowsum( des * p, group))
  return(info.des)
}

# Infodes joined with utility balance function
# InfoDes2 <- function(par, des, n.alts, utbal = NULL) {
#   group <- rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)
#   # probability
#   u <- des %*% diag(par)
#   u <- .rowSums(u, m = nrow(des), n = length(par))
#   p <- exp(u) / rep(rowsum(exp(u), group), each = n.alts)
#   if (isTRUE(utbal)) {
#     return(p)
#   } else {
#     # information matrix
#     info.des <- crossprod(des * p, des) - crossprod(rowsum( des * p, group))
#     return(info.des)
#   }
# } 

# Utility balance 
Utbal <- function(par, des, n.alts) { 
  group <- rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)
  u <- des %*% diag(par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  p <- exp(u) / rep(rowsum(exp(u), group), each = n.alts)
  # return
  return(p)
}

# KL information
#
# Calculates the Kullback-Leibler divergence for all posible choice sets, given
# @inheritParams SeqMOD
# @param full.comb A matrix in which each row is a possible combination of
#   profiles.
# @return Numeric value indicating the Kullback-Leibler divergence.
KLs <- function(full.comb, par.draws, cte.des, cand.set, weights) {
  #take set
  set <- as.matrix(cand.set[as.numeric(full.comb), ])
  #Alternative specific constants
  set <- cbind(cte.des, set)
  #calculate for all sets the KLinfo.
  kl <- KL(set, par.draws, weights)
  return(kl)
}


# Kullback-Leibler divergence for a set
#
# Calculates the KL-divergence for a choice set given parameter values.
# @inheritParams DerrS
# @param weights A vector containing the weights of the draws. Default is
#   \code{NULL}
KL <- function (set, par.draws, weights){
  # Probabilities.
  num2 <- tcrossprod(set, par.draws)
  mmat2 <- as.matrix(t(apply(num2, 2, max)))
  numm2 <- exp(sweep(num2, 2, mmat2, FUN = "-"))
  nummax <- exp(sweep(num2, 2, mmat2, FUN = "-"))
  denom <- colSums(numm2)
  ps <- sweep(nummax, 2, denom, FUN = "/")
  lgp <- log(ps)
  wprob <- sweep(ps, 2, weights, FUN="*")
  twp <- rowSums(wprob)
  lgwp <- sweep(lgp, 2, weights, FUN="*")
  tlwp <- rowSums(lgwp)
  #kullback Leibler information
  klinfo <- sum(twp * (log(twp) - tlwp))
  return (as.numeric(klinfo))
}








