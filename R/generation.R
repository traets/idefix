
#' Profiles generation
#'
#' Function to generate all possible combinations of attribute levels (i.e. all possible profiles).
#' @param lvls  A vector which contains for each attribute, the number of levels.
#' @param coding Type op coding that need to be used. See ?contrasts for more information.
#' @param intercept Logical argument indicating whether an intercept should be included. The Default is False.
#' @return A list containing 2 matrices, one contains all possible profiles with discrete levels, the other contains the coded version.
#' @export
profiles<- function (lvls, coding, intercept = FALSE) {

  lvls <- lapply(X = as.list(lvls), function(x) (1:x))
  D <- as.data.frame(expand.grid(lvls))
  cn <- names(D)
  D[, cn] <- lapply(D[, cn],  factor)
  con <- list()

  for (i in 1:length(lvls)) {
    name <- paste("Var", i, sep = "")
    con[name] <- coding
  }

  CD <- as.data.frame(model.matrix(~., D, contrasts = con))
  if (intercept == F) { CD <- CD[, -1]}

  return(list(D, as.matrix(CD)))
}

#' Random design generation
#'
#' Function to generate a random design matrix.
#' @param lvls  A vector which contains for each attribute, the number of levels.
#' @param n_sets Numeric value indicating the number of choide sets.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param coding Type op coding that need to be used. See ?contrasts for more information.
#' @param intercept Logical argument indicating whether an intercept should be included. The default is False.
#' @return A design matrix
#' @export
rdes<-function (lvls, n_sets, n_alts, coding, intercept=FALSE){

  profs<-profiles(lvls = lvls, coding = coding, intercept = intercept)[[2]]

  R<-round(runif((n_alts*n_sets), 1, nrow(profs)))
  des<-profs[R,]
  des<-as.matrix(des)

  #return design
  return(des)
}

#' Lattice multivariate standard normal distribution.
#'
#' Generates a grid of points coming from a multivariate standard normal distribution.
#' @param K Numeric value indicating the dimensionality of the grid.
#' @param b Numeric value indicating the base.
#' @param m Numeric value. Number of samples=b^m.
#' @return Matrix of lattice points drawn from a multivariate standard normal distribution. Each row is a sample.
lattice <- function (K, b, m){


  base <- function(num){

    a1 <- c1 <- rep(0, m)
    a1[m] <- c1[m] <- num %% b

    for(j in seq(m - 1, 1, by = -1)){
      tem <- (num - c1[j + 1]) / (b^(m - j))
      a1[j] <- tem %% b
      c1[j] <- c1[j + 1] + a1[j] * (b^(m - j))
    }

    return(a1)
  }

  a <- 1571
  N <- b^m # number of lattice points
  u <- runif(K)
  av <- rep(1, K)

  for(i in 2:K) {av[i] <- (a * av[i - 1]) %% N}

  e <- matrix(0, N, K)
  seq <- rep(0, m)
  kk <- -m

  for(k in 1:m) {seq[k] <- b^kk; kk <- kk + 1;}

  for(i in 1:N){
    ei <- crossprod(seq, base(i - 1)) * av + u
    e[i, ] <- ei - floor(ei)
  }

  latt <- matrix(0, N, K)
  for(i in 1:K){
    for(j in 1:N){latt[j, i] <- 1 - abs(2 * e[j, i] - 1)}
  }

  latt <- qnorm(latt)

  return(latt)
}

#' Lattice multivariate t-distribution.
#'
#' Generates a grid of points coming from a multivariate t-distribution.
#' @param mean Numeric vector indicating the multivariate mean.
#' @param cvar A matrix which specifies the covariance matrix.
#' @param df Numeric value indicating the degrees of freedom for the multivariate t-distribution.
#' @param m Numeric value. Number of samples = b^m.
#' @param b Numeric value indicating the base (default = 2).
#' @return Matrix of lattice points drawn from a multivariate t-distribution. Each row is a sample.
#' @export
lattice_mvt<- function (mean, cvar, df, m, b=2){

  dim<-length(mean)

  #gen lattice from standard normal
  lattice<- lattice(K = dim, b, m)

  mean <- t(mean)
  X <- matrix(NA, nrow(lattice), dim)
  A <- chol(cvar)

  for(i in 1:nrow(lattice))
  {

    inv<-rgamma(1, df/2)
    invw<-inv/(df/2)
    W<-1/invw

    Z<-lattice[i,]

    r<-mean + sqrt(W) * (Z %*% t(A))
    X[i, ]<-t(r)
  }

  return (X)
}

#' Lattice multivariate normal distribution.
#'
#' Generates a grid of points coming from a multivariate normal distribution.
#' @param mean Numeric vector indicating the multivariate mean.
#' @param cvar A matrix which specifies the covariance matrix.
#' @param m Numeric value. Number of samples = b^m.
#' @param b Numeric value indicating the base (default = 2).
#' @return Matrix of lattice points drawn from a multivariate normal distribution. Each row is a sample.
#' @export
lattice_mvn<-function (mean, cvar, m, b=2) {

  dim <- length(mean)
  lattice <- lattice(K = dim, b, m)
  mean <- t(mean)
  X <- matrix(NA, nrow(lattice), dim)
  A <- chol(cvar)

  for (i in 1:nrow(lattice)) {
    Z <- lattice[i, ]
    r <- mean + (Z %*% t(A))
    X[i, ] <- t(r)
  }
  return(X)
}

#' All choice sets
#'
#' Generates all possible combinations of choice sets.
#' @param candset A numeric matrix in which each row is a possible profile.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @return Matrix with all possible combinations of profiles (choice sets).
#' @export
full_sets<- function(candset, n_alts){

  fun<-function(x){return(1:x)}

  s1<-as.list(rep(nrow(candset), n_alts))
  s2<-lapply(X = s1, fun)
  fc<-as.data.frame(expand.grid(s2))

  return(fc)
}


#' Minimum difference choice sets
#'
#' Filters all possible combinations of choice sets,
#' taking into account the minimum ammount of attributelevels that need to be different per choice set.
#' @param candset A numeric matrix in which each row is a possible profile.
#' @param fcomb A matrix with all possible combinations of profiles (choice sets).
#' @return Matrix with all possible combinations of profiles, taking into account the minimum difference between profiles per choice set.
#' @export
mindiff<-function(candset,fcomb,lvls, mindiff){

  parplace<-lvls-1
  par2<-cumsum(parplace)
  par1<-c(1,par2+1)
  par1<-par1[-length(par1)]
  diff<-numeric()
  newfcomb<-numeric()
  cand<-candset

  for (s in 1:nrow(fcomb)){

    set <- cand[as.numeric(fcomb[s,]), ]

    for (pp in 1:length(par1)){
      ifelse((sum(abs(diff(set[ , par1[pp] :par2[pp]],1))) !=0), dif<-1, dif<-0)
      diff<-c(diff,dif)
    }

    if (sum(diff) > mindiff){newfcomb<-rbind(newfcomb, fcomb[s, ])}
    diff<-numeric(0)
  }

  return(newfcomb)

}

#roxygen2::roxygenise()

