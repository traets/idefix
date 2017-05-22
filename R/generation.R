

#' Profiles generation.
#'
#' Function to generate all possible combinations of attribute levels (i.e. all possible profiles).
#' @param lvls  A vector which contains for each attribute, the number of levels.
#' @param coding Type op coding that need to be used. See \code{\link[stats]{model.matrix}}
#' @param intercept Logical argument indicating whether an intercept should be included. The default is False.
#' @return A matrix which contains all possible profiles.
#' @examples 
#' # numeric vector specifing the number of levels for each attribute.
#' attribute.levels <- c(3,4,2) #3 attributes with respectively 3, 4 and 2 levels. 
#' coding.type <- "contr.sum" #effects coding will be applied.
#' #generate all coded profiles.
#' Profiles(lvls = attribute.levels, coding = coding.type)
#' # without coding
#' Profiles(lvls = c(4,4,3), coding = "none") 
#' @export
Profiles <- function(lvls, coding, intercept = FALSE) {
  #error handling
  codings.types<-c("contr.treatment", "contr.helmert", "contr.poly", "contr.sum", "none")
  if (!(coding %in% codings.types)) {
    stop("Coding argument is incorrect.")
  } 
  #error handling
  if (length(lvls) < 2 || (!(is.numeric(lvls)))){
    stop("lvls argument is incorrect")
  }
  #create all combinations of attribute levels
  levels.list <- lapply(X = as.list(lvls), function(x) (1:x))
  dgrid <- as.data.frame(expand.grid(levels.list))
  cn <- names(dgrid)
  #create factors
  dgrid[, cn] <- lapply(dgrid[, cn],  factor)
  # Change to coded version if necesarry
  if (coding != "none") {
    #create contrast vector
    con <- list()
    for (i in 1:length(lvls)) {
      name <- paste("Var", i, sep = "")
      con[name] <- coding
    }
    #create coded version
    cgrid <- as.data.frame(model.matrix(~., dgrid, contrasts = con))
    #delete intercept if necesarry
    if (intercept == F) { 
      cgrid <- cgrid[, -1]
    }
    #return coded version
    return(as.matrix(cgrid))
  #else return uncoded version
  } else {
  return(data.matrix(dgrid, rownames.force = TRUE))
  }
 
}


#' Random design generation
#'
#' Function to generate a random design matrix.
#' @inherit Profiles
#' @param n.sets Numeric value indicating the number of choide sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice set.
#' @return A design matrix
#' @export
Rdes <- function(lvls, n.sets, n.alts, coding, intercept=FALSE) {
  #generate all possible profiles
  profs <- Profiles(lvls = lvls, coding = coding, intercept = intercept)
  #draw random profiles 
  r <- round(runif((n.alts*n.sets), 1, nrow(profs)))
  des <- as.matrix(profs[r, ])
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
Lat <- function(K, b, m) {
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
  N <- b^m   #number of lattice points
  u <- runif(K)
  av <- rep(1, K)
  for (i in 2:K) {
    av[i] <- (a * av[i - 1]) %% N
  }
  e <- matrix(0, N, K)
  seq <- rep(0, m)
  kk <- -m
  for (k in 1:m) {
    seq[k] <- b^kk
    kk <- kk + 1
  }
  for (i in 1:N) {
    ei <- c(crossprod(seq, base(i - 1))) * av + u
    e[i, ] <- ei - floor(ei)
  }
  latt <- matrix(0, N, K)
  for (i in 1:K) {
    for (j in 1:N) {
      latt[j, i] <- 1 - abs(2 * e[j, i] - 1)
    }
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
lattice_mvt <- function (mean, cvar, df, m, b=2) {
  # Dimension
  dim <- length(mean)
  # Generate lattice from standard normal
  lattice <- lattice(K = dim, b, m)
  mean <- t(mean)
  X <- matrix(NA, nrow(lattice), dim)
  A <- chol(cvar)
  # Transform to multivariate t distribution
  for (i in 1:nrow(lattice)) {
    inv <- rgamma(1, df / 2)
    invw <- inv / (df / 2)
    W <- 1 / invw
    Z <- lattice[i, ]
    r <- mean + sqrt(W) * (Z %*% t(A))
    X[i, ] <- t(r)
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
lattice_mvn <- function (mean, cvar, m, b=2) {
  dim <- length(mean)
  l <- Lat(K = dim, b, m)
  mean <- t(mean)
  X <- matrix(NA, nrow(l), dim)
  A <- chol(cvar)
  for (i in 1:nrow(l)) {
    Z <- l[i, ]
    r <- mean + (Z %*% t(A))
    X[i, ] <- t(r)
  }
  return(X)
}


#' All choice sets
#'
#' Generates all possible combinations of choice sets.
#' @param cand.set A numeric matrix in which each row is a possible profile.
#' @param n.alts Numeric value indicating the number of alternatives per choice set.
#' @return Matrix with all possible combinations of profiles (choice sets).
#' @export
FullFact <- function(cand.set, n.alts) {
  fun <- function(x) {
    return(1:x)
  }
  s1 <- as.list(rep(nrow(cand.set), n.alts))
  s2 <- lapply(X = s1, fun)
  fc <- as.data.frame(expand.grid(s2))
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
mindiff <- function(candset, fcomb, lvls, mindiff) {
  parplace <- lvls-1
  par2 <- cumsum(parplace)
  par1 <- c(1, par2+1)
  par1 <- par1[-length(par1)]
  diff <- numeric()
  newfcomb <- numeric()
  cand <- candset
  for (s in 1:nrow(fcomb)) {
    set <- cand[as.numeric(fcomb[s, ]), ]
    for (pp in 1:length(par1)) {
      ifelse((sum(abs(diff(set[, par1[pp]:par2[pp]], 1))) != 0), dif <- 1, dif <- 0)
      diff <- c(diff, dif)
    }
    if (sum(diff) > mindiff) {
      newfcomb <- rbind(newfcomb, fcomb[s, ])
    }
    diff <- numeric(0)
  }
  return(newfcomb)
}


