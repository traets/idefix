

# D-error
# 
# Function to calculate d error given a design, and parameter values.
# @param par Vector containing parameter values.
# @param des A design matrix in which each row is a profile.
# @param n.alts Numeric value indicating the number of alternatives per choice
#   set.
# @return D-error.
Derr <- function(par, des, n.alts) {
  info.des <- InfoDes(par, des, n.alts)
  detinfo <- det(info.des)
  ifelse((detinfo <= 0), return(NA), return(detinfo^(-1 / length(par))))
}

#with covariance information
DerrC <- function(par, des, n.alts, i.cov) {
  info.des <- InfoDes(par, des, n.alts)
  detinfo <- det(info.des + i.cov)
  ifelse((detinfo <= 0), return(NA), return(detinfo^(-1 / length(par))))
}

# Sequential D-error
# 
# Function to calculate D-errors if set would be part of design.
# @inheritParams Modfed
# @param set A choice set in which each row is a profile.
# @param des A design matrix in which each row is a profile.
# @param i.cov Inverse of covariance matrix.
# @param n.par Number of parameters.
DerrS <- function(par.draws, set, des, n.alts, i.cov, n.par) {
  des.f <- rbind(des, set) 
  info.d <- InfoDes(par = par.draws, des = des.f, n.alts = n.alts) 
  d.error <- det(info.d + i.cov)^(-1 / n.par)
  return(d.error)
}

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
#   with of the initial design with a potential choice set.
DBerrS <- function(full.comb, cand.set, par.draws, des, n.alts, cte.des, i.cov, n.par, weights) {
  # Take set.
  set <- as.matrix(cand.set[as.numeric(full.comb), ])
  # Add alternative specific constants if necessary
  if (!is.null(cte.des)) {
    set <- as.matrix(cbind(cte.des, set))
  }
  # For each draw calculate D-error.
  d.errors <- apply(par.draws, 1, DerrS, set, des, n.alts, i.cov, n.par)
  w.d.errors <- d.errors * weights
  # DB-error. 
  db.error <- mean(w.d.errors, na.rm = TRUE)
  return(db.error)
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


# Utility balance 
Utbal <- function(par, des, n.alts) { 
  group <- rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)
  u <- des %*% diag(par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  p <- exp(u) / rep(rowsum(exp(u), group), each = n.alts)
  ub <- by(p, group, function(x) {max(x) - min(x)}, simplify = TRUE)
  # return
  return(ub)
}

# # KL information
# # 
# # Calculates the Kullback-Leibler divergence for all posible choice sets, given
# # @inheritParams SeqDB
# # @param full.comb A matrix in which each row is a possible combination of
# #   profiles.
# # @return Numeric value indicating the Kullback-Leibler divergence.
# KLs <- function(full.comb, par.draws, cte.des, cand.set, weights) {
#   #take set
#   set <- as.matrix(cand.set[as.numeric(full.comb), ])
#   #Alternative specific constants 
#   set <- cbind(cte.des, set)
#   #calculate for all sets the KLinfo.
#   kl <- KL(set, par.draws, weights)
#   return(kl)
# }


# # Kullback-Leibler divergence for a set
# # 
# # Calculates the KL-divergence for a choice set given parameter values.
# # @inheritParams DerrS
# # @param weights A vector containing the weights of the draws. Default is
# #   \code{NULL}
# KL <- function (set, par.draws, weights){
#   # Probabilities.
#   num2 <- tcrossprod(set, par.draws)
#   mmat2 <- as.matrix(t(apply(num2, 2, max)))
#   numm2 <- exp(sweep(num2, 2, mmat2, FUN = "-"))
#   nummax <- exp(sweep(num2, 2, mmat2, FUN = "-"))
#   denom <- colSums(numm2)
#   ps <- sweep(nummax, 2, denom, FUN = "/")
#   lgp <- log(ps)
#   wprob <- sweep(ps, 2, weights, FUN="*")
#   twp <- rowSums(wprob)
#   lgwp <- sweep(lgp, 2, weights, FUN="*")
#   tlwp <- rowSums(lgwp)
#   #kullback Leibler information
#   klinfo <- sum(twp * (log(twp) - tlwp))
#   return (as.numeric(klinfo))
# }








