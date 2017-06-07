

#' Profiles generation.
#' 
#' Function to generate all possible combinations of attribute levels (i.e. all
#' possible profiles).
#' @param lvls  A numeric vector which contains for each attribute, the number of
#'   levels.
#' @param coding Type op coding that needs to be used for each attribute.
#' @param c.lvls A list containing vectors with attributelevels for continuous
#'   attributes. The default is \code{NULL}
#' @return A matrix which contains all possible profiles.
#' @examples 
#' # Without continuous attributes
#' at.lvls <- c(3,4,2) # 3 Attributes with respectively 3, 4 and 2 levels. 
#' c.type <- rep("E", length(at.lvls)) # All Effect coded.
#' Profiles(lvls = at.lvls, coding = c.type) # Generate profiles.
#' 
#' # With continuous attributes 
#' at.lvls <- c(3,4,2) # 3 attributes with respectively 3, 4 and 2 levels. 
#' # First attribute is dummy coded, second and third are continuous. 
#' c.type <- c("D", "C", "C") 
#' # Levels for continuous attributes, in the same order. 
#' con.lvls <- list(c(4,6,8,10), c(7,9))
#' Profiles(lvls = at.lvls, coding = c.type, c.lvls = con.lvls)
#' @export
Profiles <- function(lvls, coding, c.lvls = NULL) {
  # Continuous attributes. 
  contins <-  which(coding == "C")
  n.contins <-  length(contins)
  # error continuous levels 
  if (!is.null(c.lvls) && !is.list(c.lvls)) { 
    stop('c.lvls should be a list.')
  }
  # Error correct coding types.
  codings.types <- c("E", "D", "C")
  if (!all(coding %in% codings.types) || (length(coding) != length(lvls))) {
    stop("coding argument is incorrect.")
  } 
  # Error lvls vector.
  if (length(lvls) < 2 || (!(is.numeric(lvls)))){
    stop("lvls argument is incorrect.")
  }
  # Error continuous specified and NULL.
  if (length(contins) > 0 && is.null(c.lvls)) {
    stop("there are no levels provided for the continuous attributes")
  }
  # Error continuous levels specification. 
  if (!is.null(c.lvls)) {
    if (length(c.lvls) != n.contins) {
      stop("length of c.lvls does not match number of specified continuous attributes in coding")
    }
    # Error c.lvls same number of levels. 
    if (!isTRUE(all.equal(lvls[contins], lengths(c.lvls)))) {
      stop("the number of continuous attribute levels provided in c.lvls does not match the expected based on lvls")
    }
  }
  # Change into correct coding. 
  coding <- dplyr::recode(coding, D = "contr.treatment", E = "contr.sum")
  # Create all combinations of attribute levels.
  levels.list <- lapply(X = as.list(lvls), function(x) (1:x))
  # Replace continuous.
  levels.list[contins] <- c.lvls
  # Create grid. 
  dgrid <- as.data.frame(expand.grid(levels.list))
  # Apply coding to non continuous. 
  cn <- names(dgrid)
  if (!is.null(c.lvls)) {
    cn <- cn[-contins]
  }
  # Create factors.
  dgrid[, cn] <- apply(dgrid[, cn, drop = FALSE], 2, factor)
  # coding 
  con <- as.list(setNames(coding, names(dgrid)))
  con[which(con == "C")] <- NULL
  cgrid <- as.data.frame(model.matrix(~., dgrid, contrasts = con))
  # Delete intercept.
  cgrid <- cgrid[, -1]
  # Return profiles.
  return(as.matrix(cgrid))
}

#' Random design generation
#'
#' Function to generate a random design matrix.
#' @inherit Profiles
#' @param n.sets Numeric value indicating the number of choide sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice set.
#' @return A design matrix.
Rdes <- function(lvls, n.sets, n.alts, coding, c.lvls = NULL) {
  #generate all possible profiles
  profs <- Profiles(lvls = lvls, coding = coding, c.lvls = c.lvls)
  #draw random profiles 
  r <- round(runif((n.alts*n.sets), 1, nrow(profs)))
  des <- as.matrix(profs[r, ])
  #return design
  return(des)
}


#' Create alternative specific coding.
Altspec <- function(alt.cte, n.sets) {
  # create matrix
  mat <- diag(length(alt.cte))
  n.zero <- which(alt.cte == 0)
  mat[n.zero, n.zero] <- 0
  # delete zero columns
  del.col <- c(which(apply(mat, 2,   function(x) all(x == 0))))
  mat <- mat[, -del.col]
  #rbind for full design 
  mat <- as.matrix(mat)
  cte.mat <- do.call(rbind, replicate(n.sets, mat, simplify = FALSE)) 
  #return
  return(cte.mat)
}


#' Create row and column names for designs 
Rcnames <- function(n.sets, n.alts, n.cte, alt.cte) {
  # rownames
  r.s <- rep(1:n.sets, each = n.alts)
  r.a <- rep(1:n.alts, n.sets)
  r.names <- paste(paste("set", r.s, sep = ""), paste("alt", r.a, sep = ""), sep = ".")
  # colnames alternative specific constants 
  cte.names <- paste(paste("alt", which(alt.cte == 1), sep = ""), ".cte", sep = "") 
  # return
  return(list(r.names, cte.names))
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
Lattice_mvt <- function (mean, cvar, df, m, b=2) {
  # Dimension
  dim <- length(mean)
  # Generate lattice from standard normal
  lattice <- Lat(K = dim, b, m)
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
Lattice_mvn <- function(mean, cvar, m, b=2) {
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


