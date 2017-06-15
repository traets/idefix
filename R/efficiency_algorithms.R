

#' Modified Federov algorithm for MNL models.
#' 
#' The algorithm swipes every profile of an initial design with candidate 
#' profiles. By doing this it tries to minimize the D(B)-error, assuming a
#' multinomial logit model.
#' 
#' The algorithm stops when an iteration occured without replacing a profile or 
#' when \code{max.iter} is reached. An iteration is a loop through all profiles 
#' from the initial design, evaluating the change in D(B)-error for every 
#' profile from \code{cand.set}.
#' 
#' By specifying a numeric vector in \code{par.samples}, the D-error will be 
#' calculated and the design will be optimised locally. By specifying a matrix, 
#' in which each row is a sample from a multivariate distribution, the DB-error
#' will be calculated, and the design will be optimised globally. The number of
#' columns should equal the number of parameters needed for \code{alt.cte} + the
#' number of parameters needed for \code{cand.set}. This is also the order in
#' which they should be sorted (first \code{alt.cte} parameters).
#' 
#' The DB-error is calculated by taking the mean over D-errors. It could be that
#' for some samples the design results in an infinite D-error. The percentage of
#' samples for which this was true for the final design can be found in the
#' output \code{inf.error}.
#' 
#' Alternative specific constants can be specified in \code{alt.cte}. The lenght
#' of this binary vector should equal \code{n.alts}, were 0 indicates no 
#' alternative specific constant and 1 the opposit.
#' 
#' @param cand.set A numeric matrix in which each row is a possible profile. The
#'   \code{\link{Profiles}} function can be used to generate this.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice 
#'   set.
#' @param alt.cte A binary vector indicating for each alternative if an
#'   alternative specific constant is desired.
#' @param par.samples A matrix in which each row is a sample from the 
#'   multivariate prior parameter distribution.
#' @param start.des A matrix in which each row is a profile. The number of rows 
#'   equals \code{n.sets * n.alts}, and the number of columns equals the number 
#'   of columns of \code{cand.set}
#' @param max.iter A numeric value indicating the maximum number allowed 
#'   iterations.
#' @return \item{design}{A numeric matrix wich represents an efficient design.} 
#' \item{error}{Numeric value indicating the D(B)-error of the design.} 
#' \item{inf.error}{Numeric value indicating the percentage of samples for which
#' the D-error was \code{Inf}.} \item{prob.diff}{Numeric value indicating the
#' difference between the alternative with the highest and the one with the
#' lowest probability for each choice set. If a sample matrix was provided this
#' is based on the mean over all samples.}
#' @examples
#' # D-efficient design
#' # 3 Attributes, 2 are dummy coded and 1 continuous (= 3 parameters).
#' cs <- Profiles(lvls = c(2, 3, 2), coding = c("D", "C", "D"), c.lvls = list(c(2,4,6)))
#' ps <- c(0.8, 0.2, -0.3) # Prior parameter vector
#' Modfed(cand.set = cs, n.sets = 8, n.alts = 2, alt.cte = c(0,0), par.samples = ps)
#' 
#' # DB-efficient design. 
#' # 3 Attributes with 2, 3 and 2 levels, all effect coded (= 4 parameters).
#' cs <- Profiles(lvls = c(2, 3, 2), coding = c("E", "E", "E")) 
#' m <- c(0.8, 0.2, -0.3, -0.2, 0.7) # Prior mean (total = 5 parameters).
#' v <- diag(length(m)) # Prior variance. 
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = v) # 10 Samples.
#' Modfed(cand.set = cs, n.sets = 8, n.alts = 2, alt.cte = c(1,0), par.samples = ps)
#' @references
#' \insertRef{federov}{mnldes} 
#' @export
Modfed <- function(cand.set, n.sets, n.alts,  alt.cte, par.samples, start.des = NULL, max.iter = Inf) {
  # Handling par.samples.
  if (!(is.matrix(par.samples))) {
    par.samples <- matrix(par.samples, nrow = 1)
  }
  # Error alternative specific constants. 
  if (length(alt.cte) != n.alts) {
    stop("n.alts does not match the alt.cte vector")
  }
  # Error identifying model.
  if (n.sets < ncol(par.samples)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  # Errors start design. 
  if (!is.null(start.des)) {
    if (ncol(start.des) != ncol(cand.set)) {
      stop("number of colums start design is different from number of columns candidate set.")
    }
    if (nrow(start.des) != (n.alts * n.sets)) {
      stop("number of rows start design is different from number of sets times number of alternatives.")
    }
  }
  # Create alternative specific design.
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
  # Error handling cte.des
  if (ncol(cand.set) + ncol(cte.des) != ncol(par.samples)) {
    stop("dimension of par.samples does not match the dimension of alt.cte + cand.set.")
  }
  # Random start design.
  if (!is.null(start.des)) {
    des <- start.des
  } else {
    r <- round(runif((n.sets*n.alts), 1, nrow(cand.set)))
    des <- data.matrix(cand.set[r, ])
  }
  # Combine with alt.spec design.
  des <- cbind(cte.des, des)
  # Starting values. 
  d.start <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
  db.start <- mean(d.start, na.rm = TRUE)
  converge <- FALSE
  change <- FALSE
  it <- 1
  n.samples <- nrow(par.samples)
  n.cte <- ncol(cte.des)
  n.par <- ncol(des)
  # start algorithm.
  while (!converge & it <= max.iter) {
    it <- it + 1
    # show progress iteration.
    pb <- txtProgressBar(min = 0, max = nrow(des), style = 3)     
    # save design before iteration.
    iter.des <- des
    # For every row in the design.
    for (r in 1:nrow(des)) {
      # Switch with everey row in candidate set. 
      for (c in 1:nrow(cand.set)) {
        des[r, (n.cte + 1) : n.par ] <- cand.set[c, ]
        # Calculate D-errors.
        d.errors <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
        # DB-error. 
        db <- mean(d.errors, na.rm = TRUE)
        # Change if lower db error.
        if (db < db.start) {
          best.row <- as.numeric(des[r, ])
          db.start <- db
          change <- TRUE
        }
        # Show progress iteration.
        setTxtProgressBar(pb, r)
      }
      # Replace with best profile if change.
      if (change) {
        des[r, ] <- best.row
      } else {
        des[r, ] <- iter.des[r, ]
      }
      # Initialize variables again. 
      change <- FALSE
      na.percentage <- 0
    }
    close(pb)  # Print progress. 
    converge <- isTRUE(all.equal(des, iter.des)) # Convergence if no profile is swapped this iteration.
  }
  # calculate percentage NA values.
  d.errors <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
  if (any(is.na(d.errors))) {
    na.percentage <- percent(sum(is.na(d.errors))/n.samples)
  } 
  # Utility balance.
  ub <- apply(par.samples, 1, Utbal, des = des,  n.alts = n.alts)
  ub <- .rowMeans(ub, m = n.sets, n = n.samples, na.rm = FALSE)
    # Rownames design. 
  des.names <- Rcnames(n.sets = n.sets, n.alts = n.alts, n.cte = n.cte, alt.cte = alt.cte)
  rownames(des) <- des.names[[1]]
  # Colnames alternative specific constants. 
  if (n.cte != 0) {
    colnames(des)[1:n.cte] <- des.names[[2]]
  }
  # Return design, D(B)error, percentage NA's, utility balance. 
  return(list("design" = des, "error" =  db, "inf.error" = na.percentage, "prob.diff" = ub))
}


#' Sequential modified federov algorithm for MNL model.
#' 
#' Selects the choice set that minimizes the DB-error when added to an initial
#' design, given (updated) parameter values.
#' 
#' This algorithm is ideally used in an adaptive context. The algorithm will 
#' select the next DB-efficient choice set given parameter values and an initial
#' design. In an adaptive context these parameter values are updated after each 
#' observed response.
#' 
#' The initial design \code{des} can be generated with \code{\link{Modfed}}. If 
#' alternative specific constants are included in the initial design, the 
#' algorithm will use the same for selecting the new choice set. Columns of 
#' \code{des} which contain ".cte" in their name are recognized as alternative 
#' specific columns.
#' 
#' The list of potential choice sets are created using 
#' \code{\link[gtools]{combinations}}. If \code{reduce} is \code{TRUE}, 
#' \code{repeats.allowed = FALSE} and vice versa. If no alternative constants
#' are used \code{reduce} should always be \code{TRUE}. When alternative
#' specific constants are used \code{reduce} can be \code{TRUE} so that the
#' algorithm will be faster, but the combinations of constants and profiles will
#' not be evaluated exhaustively.
#' 
#' The \code{weights} can be used when the \code{par.samples} have weights. This
#' is for example the case when parameter values are updated using 
#' \code{\link{ImpSamp}}.
#' @inheritParams Modfed
#' @param par.samples A matrix in which each row is a sample from the 
#'   multivariate parameter distribution. See also \code{\link{ImpsamplingMNL}}.
#' @param des A design matrix in which each row is a profile. Can be generated
#'   with \code{\link{Modfed}}
#' @param prior.covar Covariance matrix of the prior distribution.
#' @param reduce Logical value indicating whether the candidate set should be
#'   reduced or not.
#' @param weights A vector containing the weights of the samples. Default is
#'   \code{NULL}, See also \code{\link{ImpsamplingMNL}}.
#' @return \item{set}{A matrix representing a DB efficient choice set.} 
#' \item{db.error}{A numeric value indicating the DB-error of the whole design.}
#' @examples 
#' # DB efficient choice set, given a design and parameter samples. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3), coding = c("E", "E"))
#' m <- c(0.3, 0.2, -0.3, -0.2) # Prior mean (total = 5 parameters).
#' pc <- diag(length(m)) # Prior variance
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 Samples.
#' ac <- c(0, 0) # No alternative specific constants. 
#' # Initial design.
#' des <- Modfed(cand.set = cs, n.sets = 6, n.alts = 2, alt.cte = ac, par.samples = ps)$design
#' # Efficient choice set to add. 
#' SeqDB(des = des, cand.set = cs, n.alts = 2, par.samples = ps, prior.covar = pc)
#' 
#' # DB efficient choice set, given a design and parameter samples. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3), coding = c("C", "E"), c.lvls = list(c(5,3,1)))
#' m <- c(0.7, 0.3, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(m)) # Prior variance
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 Samples.
#' ac <- c(1, 0) # Alternative specific constant. 
#' # Initial design.
#' des <- Modfed(cand.set = cs, n.sets = 6, n.alts = 2, alt.cte = ac, par.samples = ps)$design
#' # Efficient choice set to add. 
#' SeqDB(des = des, cand.set = cs, n.alts = 2, par.samples = ps, prior.covar = pc)
#' @references
#' \insertRef{ju}{mnldes} 
#' @export
SeqDB <- function(des, cand.set, n.alts, par.samples, prior.covar, reduce = TRUE, weights = NULL) {
  # Initialize.
  n.sets <- nrow(des) / n.alts
  cte.des <- NULL
  # If no weights, equal weights.
  if (is.null(weights)) {
    weights <- rep(1, nrow(par.samples))
  }
  # Detect alternative specific constants
  des.f <- as.data.frame(des)
  alt.cte <- select(des.f, contains(".cte"))
  if (ncol(alt.cte) > 0) {
    cte.des <- alt.cte[1:n.alts, ]
  }
  # Handling par.samples.
  if (!(is.matrix(par.samples))) {
    par.samples <- matrix(par.samples, nrow = 1)
  }
  # Error par.samples
  if (ncol(des) != ncol(par.samples)) {
    stop("Numbers of parameters in par.samples does not match the number of parameters in the design.")
  }
  # Error identifying model.
  if (n.sets < ncol(par.samples)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  # Starting and initializing values.
  i.cov <- solve(prior.covar)
  d.start <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
  db.start <- mean(d.start, na.rm = TRUE)
  full.comb <- gtools::combinations(n = nrow(cand.set), r = n.alts, repeats.allowed = !reduce)
  n.par <- ncol(par.samples)
  # For each potential set, select best. 
  db.errors <- apply(full.comb, 1, DBerrS, cand.set, par.samples, des, n.alts, cte.des, i.cov, n.par, weights)
  comb.nr <- as.numeric(full.comb[which.min(db.errors), ])
  set <- cand.set[comb.nr, ]
  # Add alternative specific constants if necessary
  if (!is.null(cte.des)) {
    set <- cbind(cte.des, set)
  }
  row.names(set) <- NULL
  db <- min(db.errors)
  #return best set and db error design.
  return(list(set = set, db.error = db))
}


#' Sequential Kullback-Leibler based algorithm for the MNL model.
#' 
#' Selects the choice set that maximizes the Kullback-Leibler divergence between
#' prior parameter values and the expected posterior, assuming an MNL model.
#' 
#' The algorithm selects the choice set that maximizes the Kullback-Leibler 
#' divergence between prior and expected posterior. Otherwisely framed the 
#' algorithm selects the choice set that maximizes the expected information 
#' gain.
#' @inheritParams SeqDB
#' @param alt.cte A binary vector indicating for each alternative if an
#'   alternative specific constant is desired.
#' @param reduce Logical value indicating whether the candidate set should be 
#'   reduced or not.
#' @return Choice set that maximizes the expected KL divergence.
#' @references \insertRef{crabbe}{mnldes}
#' @examples 
#' # KL efficient choice set, given parameter samples. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3), coding = c("E", "E"))
#' m <- c(0.3, 0.2, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(m)) # Prior variance
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 Samples.
#' ac <- c(0, 0) # No alternative specific constants. 
#' # Efficient choice set to add. 
#' SeqKL(cand.set = cs, n.alts = 2, alt.cte = ac, par.samples = ps, weights = NULL)
#' 
#' # KL efficient choice set, given parameter samples. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3), coding = c("C", "E"), c.lvls = list(c(5,3,1)))
#' m <- c(0.7, 0.3, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(m)) # Prior variance
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 Samples.
#' ac <- c(1, 0) # Alternative specific constant. 
#' # Efficient choice set to add. 
#' SeqKL(cand.set = cs, n.alts = 2, alt.cte = ac, par.samples = ps, weights = NULL)
#' @export
SeqKL <- function(cand.set, n.alts, alt.cte, par.samples, weights, reduce = TRUE) {
  # Handling par.samples.
  if (!(is.matrix(par.samples))) {
    par.samples <- matrix(par.samples, nrow = 1)
  }
  # Error alternative specific constants. 
  if (length(alt.cte) != n.alts) {
    stop("n.alts does not match the alt.cte vector")
  }
  # Create alternative specific design.
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = 1)
  # Error handling cte.des
  if (ncol(cand.set) + ncol(cte.des) != ncol(par.samples)) {
    stop("dimension of par.samples does not match the dimension of alt.cte + cand.set.")
  }
  # All choice sets.
  full.comb <- gtools::combinations(n = nrow(cand.set), r = n.alts, repeats.allowed = !reduce)
  # If no weights, equal weights.
  if (is.null(weights)) {
    weights <- rep(1, nrow(par.samples))
  }
  # Calculate KL for each set. 
  kl.infos <- apply(full.comb, 1, KLs, par.samples, cte.des, cand.set, weights)
  # Select maximum.
  comb.nr <- as.numeric(full.comb[which.max(kl.infos), ])
  set <- cand.set[comb.nr, ]
  # Add alternative specific constants if necessary
  if (!is.null(cte.des)) {
    set <- cbind(cte.des, set)
  }
  row.names(set) <- NULL
  # return.
  return(list(set = set, kl = max(kl.infos)))
}










