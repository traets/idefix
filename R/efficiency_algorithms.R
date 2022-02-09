

#' Modified Federov algorithm for MNL models.
#' 
#' The algorithm swaps every profile of an initial start design with candidate 
#' profiles. By doing this, it tries to minimize the D(B)-error, based on a 
#' multinomial logit model. This routine is repeated for multiple starting 
#' designs.
#' 
#' Each iteration will loop through all profiles from the initial design, 
#' evaluating the change in D(B)-error for every profile from \code{cand.set}. 
#' The algorithm stops when an iteration occured without replacing a profile or 
#' when \code{max.iter} is reached.
#' 
#' By specifying a numeric vector in \code{par.draws}, the D-error will be 
#' calculated and the design will be optimised locally. By specifying a matrix, 
#' in which each row is a draw from a multivariate distribution, the DB-error 
#' will be calculated, and the design will be optimised globally. Whenever there
#' are alternative specific constants, \code{par.draws} should be a list 
#' containing two matrices: The first matrix containing the parameter draws for
#' the alternative specific constant parameters. The second matrix containing
#' the draws for the rest of the parameters.
#' 
#' The DB-error is calculated by taking the mean over D-errors. It could be that
#' for some draws the design results in an infinite D-error. The percentage of 
#' draws for which this was true for the final design can be found in the output
#' \code{inf.error}.
#' 
#' Alternative specific constants can be specified in \code{alt.cte}. The length
#' of this binary vector should equal \code{n.alts}, were \code{0} indicates the
#' absence of an alternative specific constant and \code{1} the opposite.
#' 
#' \code{start.des} is a list with one or several matrices  corresponding to 
#' initial start design(s). In each matrix each
#' row is a profile. The number of rows equals \code{n.sets * n.alts}, and the
#' number of columns equals the number of columns of \code{cand.set} + the
#' number of non-zero elements in \code{alt.cte}. If \code{start.des
#' = NULL}, \code{n.start} random initial designs will be
#' generated. If start designs are provided, \code{n.start} is ignored.
#' 
#' If \code{no.choice} is \code{TRUE}, in each choice set an alternative with
#' one alternative specific constant is added. The return value of the
#' D(B)-error is however based on the design without the no choice option.
#' 
#' When \code{parallel} is \code{TRUE}, \code{\link[parallel]{detectCores}} will
#' be used to decide upon the number of available cores. That number minus 1 
#' cores will be used to search for efficient designs. The computation time will
#' decrease significantly when \code{parallel = TRUE}.
#' 
#' @param cand.set A numeric matrix in which each row is a possible profile. The
#'   \code{\link{Profiles}} function can be used to generate this matrix.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice 
#'   set.
#' @param alt.cte A binary vector indicating for each alternative whether an 
#'   alternative specific constant is desired. The default is \code{NULL}.
#' @param par.draws A matrix or a list, depending on \code{alt.cte}.
#' @param no.choice A logical value indicating whether a no choice alternative 
#'   should be added to each choice set. The default is \code{FALSE}.
#' @param start.des A list containing one or more matrices corresponding to initial start design(s). The default is \code{NULL}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores. The default is \code{TRUE}.
#' @param max.iter A numeric value indicating the maximum number allowed 
#'   iterations. The default is \code{Inf}.
#' @param n.start A numeric value indicating the number of random start designs
#'   to use. The default is 12.
#' @param best A logical value indicating whether only the best design should be
#'   returned. The default is \code{TRUE}.
#' @return 
#'  If \code{best = TRUE} the design with the lowest D(B)-error is returned. 
#'   If \code{best = FALSE}, the results of all (provided) start designs are
#'   returned. \item{design}{A
#'   numeric matrix wich contains an efficient design.} \item{error}{Numeric
#'   value indicating the D(B)-error of the design.} \item{inf.error}{Numeric
#'   value indicating the percentage of draws for which the D-error was
#'   \code{Inf}.} \item{probs}{Numeric matrix containing the probabilities of
#'   each alternative in each choice set. If a sample matrix was provided in
#'   \code{par.draws}, this is the average over all draws.}
#' @examples
#' \dontrun{
#' # DB-efficient designs
#' # 3 Attributes, all dummy coded. 1 alternative specific constant = 7 parameters
#' cand.set <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D"))
#' mu <- c(0.5, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' p.d <- list(matrix(pd[,1], ncol = 1), pd[,2:7])
#' Modfed(cand.set = cand.set, n.sets = 8, n.alts = 2, 
#'        alt.cte = c(1, 0), parallel = FALSE, par.draws = p.d, best = FALSE)
#' 
#' # DB-efficient design with start design provided.  
#' # 3 Attributes with 3 levels, all dummy coded (= 6 parameters).
#' cand.set <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D")) 
#' mu <- c(0.8, 0.2, -0.3, -0.2, 0.7, 0.4) # Prior mean (total = 5 parameters).
#' v <- diag(length(mu)) # Prior variance.
#' sd <- list(example_design)
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' Modfed(cand.set = cand.set, n.sets = 8, n.alts = 2, 
#'        alt.cte = c(0, 0), parallel = FALSE, par.draws = ps, start.des = sd)
#'}
#' @importFrom Rdpack reprompt 
#' @importFrom MASS mvrnorm
#' @references \insertRef{idefix}{idefix}
#' @export
Modfed <- function(cand.set, n.sets, n.alts, par.draws, alt.cte = NULL, no.choice = FALSE, 
                   start.des = NULL, parallel = TRUE, max.iter = Inf, n.start = 12,
                   best = TRUE) {
  if (is.null(alt.cte)) {
    alt.cte <- rep(0L, n.alts)
  }
  #init
  n.cte <- length(which(alt.cte == 1))
  ### Errors
  if (!is.list(par.draws)) {
    if (is.vector(par.draws)) {
      par.draws <- matrix(par.draws, nrow = 1)
    }
  }
  #handling alt.cte
  if (length(alt.cte) != n.alts) {
    stop("'n.alts' does not match the 'alt.cte' vector")
  }
  if (!all(alt.cte %in% c(0, 1))) {
    stop("'alt.cte' should only contain zero or ones.")
  }
  #if no.choice
  if (!is.logical(no.choice)) {
    stop("'no.choice' should be TRUE or FALSE")
  }
  if (no.choice) {
    if (!isTRUE(all.equal(alt.cte[n.alts], 1))) {
      stop("if 'no.choice' is TRUE, alt.cte[n.alts] should equal 1.")
    }
    ncsek <- seq(n.alts, (n.sets * n.alts), n.alts) 
  } else {
    ncsek <- NULL
  }
  # Handling par.draws with alternative specific contstants.
  if (isTRUE(all.equal(n.cte, 1))) {
    if (!(is.list(par.draws))) {
      stop("par.draws should be a list")
    }
    if (!isTRUE(all.equal(length(par.draws), 2))) {
      stop("'par.draws' should contain two components")
    }
    if (is.vector(par.draws[[1]])) {
      par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1)
    }
    if (!(all(unlist(lapply(par.draws, is.matrix))))) {
      stop("'par.draws' should contain two matrices")
    }
    if (!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))) {
      stop("the first component of 'par.draws' should contain the same number 
           of columns as there are non zero elements in 'alt.cte'")
    }
    dims <-  as.data.frame(lapply(par.draws, dim))
    if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) {
      stop("the number of rows in the components of 'par.draws' should be equal")
    }
    if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
      stop("the sum of the number of columns in the components of 'par.draws' 
           should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
    }
    par.draws  <- do.call("cbind", par.draws)
  }
  if (n.cte > 1.2) {
    if (!(is.list(par.draws))) {
      stop("par.draws should be a list")
    } 
    if (!isTRUE(all.equal(length(par.draws), 2))) {
      stop("'par.draws' should contain two components")
    }
    if (!(all(unlist(lapply(par.draws, is.matrix))))) {
      stop("'par.draws' should contain two matrices")
    }
    if (!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))) {
      stop("the first component of 'par.draws' should contain the same number 
           of columns as there are non zero elements in 'alt.cte'")
    }
    dims <-  as.data.frame(lapply(par.draws, dim))
    if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) { 
      stop("the number of rows in the components of 'par.draws' should be equal")
    }
    if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
      stop("the sum of the number of columns in the components of 'par.draws' 
           should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
    }
    par.draws  <- do.call("cbind", par.draws)
  }
  # Create alternative specific design.
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
  # Error identifying model.
  if (n.sets < ncol(par.draws)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  # Handling cand.set
  if (!all(is.finite(cand.set))) {
    stop("'cand.set' contains non finite values.")
  }
  # Error handling cte.des
  if (ncol(cand.set) + ncol(cte.des) != ncol(par.draws)) {
    stop("The number of parameters in the components of 'par.draws' does not match the number 
         of non-zero parameters in 'alt.cte' + the number of parameters in 'cand.set'.")
  }
  # Random start design.
  if (!is.null(start.des)) {
    if (!is.list(start.des)) {
      stop("'start.des' should be a list")
    }
    if (!(all(unlist(lapply(start.des, is.matrix))))) {
      stop("'start.des' should contain matrices as components")
    }
    dimstart <- as.matrix(lapply(start.des, dim))
    nr.starts <- length(dimstart)
    if (nr.starts > 1.5) {
      if (!isTRUE(all.equal(length(unique(unlist(dimstart))), 2))) {
        stop("start designs have different dimensions")
      }
    }
    if (!isTRUE(all.equal(n.alts * n.sets, unique(unlist(dimstart))[1]))) {
      stop("number of rows of start design(s) does not match with 'n.alts' * 'n.sets'")
    }
    if (!isTRUE(all.equal(sum(ncol(cand.set), ncol(cte.des)), 
                          unique(unlist(dimstart))[2]))) {
      stop("number of columns of start design(s) does not match with the number
           of columns in 'cand.set' + the non zero parameters in 'alt.cte'")
    }
    d.start <- lapply(start.des, StartDB, par.draws, n.alts)
    if (!any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))) {
      stop("One or more of the provided start designs resulted in an unvalid db-error.
           Make sure the utility of any alternative cannot exceed 700.")
    }
  } 
  if (is.null(start.des)) {
    #create start designs
    nr.starts <- n.start
    start.des <- vector(mode = 'list', length = nr.starts)
    okstart <- FALSE
    while (okstart == FALSE) {
      for (i in 1:nr.starts) {
        r <- round(stats::runif((n.sets * n.alts), 1, nrow(cand.set)))
        start.des[[i]] <- cbind(cte.des, data.matrix(cand.set[r, ]))
        if (no.choice) {
          start.des[[i]][ncsek, (ncol(cte.des) + 1):(ncol(cte.des) + ncol(cand.set))] <- c(rep(0, ncol(cand.set)))
        }
      }
      d.start <- lapply(start.des, StartDB, par.draws, n.alts)
      if(!any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))){
        stop("One or more draws resulted in an unvalid d-error. 
             Make sure the utility of any alternative cannot exceed 700.")
      } else {
        okstart <- TRUE
      } 
    }
  }
  if (parallel) {
    ########
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_cores)
    parallel::clusterExport(cl, c("n.sets", "par.draws", "cand.set", "n.alts", "n.cte", "alt.cte", "no.choice", "max.iter","ncsek"), envir = environment())
    deslist <- parallel::parLapply(cl, start.des, Modfedje_ucpp, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte, no.choice, max.iter, ncsek)
    parallel::stopCluster(cl)
    ########
  } else {
    deslist <- lapply(start.des, Modfedje_ucpp, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte, no.choice, max.iter = max.iter, ncsek)
  }                                 
  bestdes <- deslist[[which.min(unlist(lapply(deslist, function(x) (x$error))))]]
  
  ifelse(best, return(bestdes), return(deslist))
}

# Core of the Modfed algorithm
Modfedje_ucpp <- function(desje, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte,
                          no.choice, max.iter, ncsek){
  converge <- FALSE
  change <- FALSE
  it <- 1
  n.samples <- nrow(par.draws)
  n.par <- ncol(desje)
  ###
  while (!converge & it <= max.iter) {
    db.start <- mean(apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts), na.rm = TRUE)
    it <- it + 1
    # save design before iteration.
    iter.des <- desje
    # For every row in the design.
    sek <- 1:nrow(desje)
    if (no.choice) {
      sek <- sek[-ncsek]
    }
    for (r in sek) {
      # Switch with everey row in candidate set. 
      db <- numeric(nrow(cand.set))
      for (c in 1:nrow(cand.set)) {
        desje[r, (n.cte + 1):n.par ] <- cand.set[c, ]
        # Calculate D-errors.
        d.errors <- apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts)
        # DB-error. 
        db[c] <- mean(d.errors, na.rm = TRUE)
      }
      pr <- which.min(db)
      db <- min(db) 
      # Change if lower db error.
      if (!is.na(db) && !is.na(db.start)) {
        if (db < db.start) {
          best.row <- as.numeric(cand.set[pr, ])
          db.start <- db
          change <- TRUE
        }
      }
      # Replace with best profile if change.
      if (change) {
        desje[r, (n.cte + 1):n.par] <- best.row
      } else {
        desje[r, ] <- iter.des[r, ]
      }
      # Initialize variables again. 
      change <- FALSE
      na.percentage <- 0
    }
    converge <- isTRUE(all.equal(desje, iter.des)) # Convergence if no profile is swapped this iteration.
  }
  # calculate percentage NA values.
  d.errors <- apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts)
  if (any(is.na(d.errors))) {
    na.percentage <- scales::percent(sum(is.na(d.errors)) / n.samples)
  } 
  # Utility balance.
  # ub <- apply(par.draws, 1, Utbal, des = desje,  n.alts = n.alts)
  # ubcpp <- apply(par.draws, 1, InfoDes_cpp, des = desje,  n_alts = n.alts, 
  #                utbal = TRUE)
  
  # Utility balance using c++ function
  ub <- apply(par.draws, 1, InfoDes_cpp, des = desje,  n_alts = n.alts, 
              utbal = TRUE)
  pmat <- matrix(rowMeans(ub), ncol = n.alts, byrow = TRUE)
  rownames(pmat) <- paste("set", 1:n.sets, sep = "")
  colnames(pmat) <- paste(paste("Pr(", paste("alt", 1:n.alts, sep = ""), 
                                sep = ""), ")", sep = "")
  if (no.choice) {
    colnames(pmat)[n.alts] <- "Pr(no choice)"
  }
  # Rownames design. 
  des.names <- Rcnames(n.sets = n.sets, n.alts = n.alts, alt.cte = alt.cte, no.choice = no.choice)
  rownames(desje) <- des.names[[1]]
  # Colnames alternative specific constants. 
  if (n.cte != 0 && !is.null(colnames(desje))) {
    colnames(desje)[1:n.cte] <- des.names[[2]]
  }
  # Return design, D(B)error, percentage NA's, utility balance. 
  return(list("design" = desje, "error" =  db.start, "inf.error" = na.percentage, "probs" = pmat))
}




#' Sequential modified federov algorithm for MNL model.
#' 
#' Selects the choice set that minimizes the DB-error when added to an initial 
#' design, given (updated) parameter values.
#' 
#' This algorithm is ideally used in an adaptive context. The algorithm will 
#' select the next DB-efficient choice set given parameter values and possible
#' previously generated choice sets. In an adaptive context these parameter
#' values are updated after each observed response.
#' 
#' Previously generated choice sets, which together form an initial design, can
#' be provided in \code{des}. When no design is provided, the algorithm will
#' select te most efficient choice set based on the fisher information of the
#' prior covariance matrix \code{prior.covar}.
#' 
#' If \code{alt.cte = NULL}, \code{par.draws} should be a matrix in which each 
#' row is a sample from the multivariate parameter distribution. In case that 
#' \code{alt.cte} is not \code{NULL}, a list containing two matrices should be 
#' provided to \code{par.draws}. The first matrix containing the parameter draws
#' for the alternative specific parameters. The second matrix containing the
#' draws for the rest of the parameters.
#' 
#' The list of potential choice sets are created using 
#' \code{\link[utils]{combn}}. If \code{reduce} is \code{TRUE}, 
#' \code{allow.rep = FALSE} and vice versa. Furthermore, the list of 
#' potential choice sets will be screaned in order to select only those choice 
#' sets with a unique information matrix. If no alternative specific constants are used, 
#' \code{reduce} should always be \code{TRUE}. When alternative specific 
#' constants are used \code{reduce} can be \code{TRUE} so that the algorithm 
#' will be faster, but the combinations of constants and profiles will not be 
#' evaluated exhaustively.
#' 
#' The \code{weights} argument can be used when the \code{par.draws} have 
#' weights. This is for example the case when parameter values are updated using
#' \code{\link{ImpsampMNL}}.
#' 
#' When \code{parallel} is \code{TRUE}, \code{\link[parallel]{detectCores}} will
#' be used to decide upon the number of available cores. That number minus 1 
#' cores will be used to search for the optimal choice set. For small problems 
#' (6 parameters), \code{parallel = TRUE} can be slower. For larger problems the
#' computation time will decrease significantly.
#' 
#' *Note:* this function is more stable than \code{\link[idefix]{SeqCEA}}, but 
#' it takes more time to get the output. This happens because this function 
#' makes an exhaustive search to get the choice set, whereas 
#' \code{\link[idefix]{SeqCEA}} makes a random search.
#' 
#' @inheritParams Modfed
#' @param par.draws A matrix or a list, depending on \code{alt.cte}. 
#' @param des A design matrix in which each row is a profile. If alternative 
#'   specific constants are present, those should be included as the first 
#'   column(s) of the design. Can be generated with \code{\link{Modfed}} or \code{\link{CEA}}.
#' @param prior.covar Covariance matrix of the prior distribution.
#' @param weights A vector containing the weights of the draws. Default is 
#'   \code{NULL}, See also \code{\link{ImpsampMNL}}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores.
#' @param no.choice An integer indicating the no choice alternative. The default
#'   is \code{NULL}.
#' @param reduce Logical value indicating whether the candidate set should be 
#'   reduced or not.
#' @param allow.rep Logical value indicating whether repeated choice sets are
#' allowed in the design.
#' @return \item{set}{A matrix representing a DB efficient choice set.} 
#'   \item{error}{A numeric value indicating the DB-error of the whole 
#'   design.}
#' @importFrom Rdpack reprompt
#' @references \insertRef{idefix}{idefix}
#' @references \insertRef{ju}{idefix}
#' @examples 
#' # DB efficient choice set, given a design and parameter draws. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D"))
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4) # mean (total = 6 parameters).
#' pc <- diag(length(m)) # covariance matrix
#' set.seed(123)
#' sample <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
#' # Initial design.
#' des <- example_design 
#' # Efficient choice set to add. 
#' SeqMOD(des = des, cand.set = cs, n.alts = 2, par.draws = sample, 
#'            prior.covar = pc, parallel = FALSE)
#' 
#' # DB efficient choice set, given parameter draws. 
#' # with alternative specific constants 
#' des <- example_design2 
#' cs <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D"))
#' ac <- c(1, 1, 0) # Alternative specific constants. 
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4, 1.8, 1.2) # mean 
#' pc <- diag(length(m)) # covariance matrix
#' pos <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
#' sample <- list(pos[ , 1:2], pos[ , 3:8])
#' # Efficient choice set. 
#' SeqMOD(des = des, cand.set = cs, n.alts = 3, par.draws = sample, alt.cte = ac, 
#'            prior.covar = pc, parallel = FALSE)
#' @export
SeqMOD <- function(des = NULL, cand.set, n.alts, par.draws, prior.covar, 
                  alt.cte = NULL, no.choice = NULL, weights = NULL, 
                  parallel = TRUE, reduce = TRUE, 
                  allow.rep = FALSE) {
  #init
  if (is.null(des)) {
    n.sets <- 1L
    allow.rep <- TRUE
  } else { 
    if (!is.matrix(des)) {
      stop("'des' should be a matrix or NULL")
    }
    if (!isTRUE(nrow(des) %% n.alts == 0)) {
      stop("'n.alts' does not seem to match with the number of rows in 'des'")
    }
    n.sets <- nrow(des) / n.alts
  }
  # if alternative constants 
  if (!is.null(alt.cte)) {
    if (length(alt.cte) != n.alts) {
      stop("'n.alts' does not match the 'alt.cte' vector")
    }
    if (!all(alt.cte %in% c(0, 1))) {
      stop("'alt.cte' should only contain zero or ones.")
    }
    # alternative specific constants
    n.cte <- length(which(alt.cte == 1L))
    if (isTRUE(all.equal(n.cte, 0L))) {
      alt.cte <- NULL
      cte.des <- NULL
    }
  } else {
    n.cte <- 0
  }
  #if no.choice
  if (!is.null(no.choice)) {
    if (!is.wholenumber(no.choice)) {
      stop("'no.choice' should be an integer or NULL")
    }
    if (any(isTRUE(no.choice > (n.alts + 0.2)), isTRUE(no.choice < 0.2))) {
      stop("'no.choice' does not indicate one of the alternatives")
    }
    if (is.null(alt.cte)) {
      stop("if there is a no choice alternative, 'alt.cte' should be specified")
    }
    if (!isTRUE(all.equal(alt.cte[no.choice], 1))) {
      stop("the no choice alternative should correspond with a 1 in 'alt.cte'")
    }
  }
  if (!is.null(alt.cte)) {
    #prior.covar
    if (!isTRUE(all.equal(ncol(prior.covar), (ncol(cand.set) + n.cte)))) {
      stop("number of columns of 'prior.covar' does not equal 
           the number of columns in 'cand.set' + nonzero elements in 'alt.cte'")
    }
    if (isTRUE(all.equal(n.cte, 1))) {
      if (!is.list(par.draws)) {stop("'par.draws' should be a list when 'alt.cte' is not NULL")}
      if (!isTRUE(all.equal(length(par.draws), 2))) {
        stop("'par.draws' should contain two components")
      }
      if (is.vector(par.draws[[1]])) {
        par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1)
      }
      if (!(all(unlist(lapply(par.draws, is.matrix))))) {
        stop("'par.draws' should contain two matrices")
      }
      if (!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))) {
        stop("the first component of 'par.draws' should contain the same number 
             of columns as there are non zero elements in 'alt.cte'")
      }
      dims <-  as.data.frame(lapply(par.draws, dim))
      if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) { 
        stop("the number of rows in the components of 'par.draws' should be equal")
      }
      if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
        stop("the sum of the number of columns in the components of 'par.draws' 
             should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
      }
      par.draws  <- do.call("cbind", par.draws)
    }
    if (n.cte > 1.2) {
      if (!(is.list(par.draws))) {stop("'par.draws' should be a list when 'alt.cte' is not NULL")} 
      if (!isTRUE(all.equal(length(par.draws), 2))) {
        stop("'par.draws' should contain two components")
      }
      if (!(all(unlist(lapply(par.draws, is.matrix))))) {
        stop("'par.draws' should contain two matrices")
      }
      if (!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))) {
        stop("the first component of 'par.draws' should contain the same number 
             of columns as there are non zero elements in 'alt.cte'")
      }
      dims <-  as.data.frame(lapply(par.draws, dim))
      if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) { 
        stop("the number of rows in the components of 'par.draws' should be equal")
      }
      if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
        stop("the sum of the number of columns in the components of 'par.draws' 
             should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
      }
      par.draws  <- do.call("cbind", par.draws)
    }
    # Create alternative specific design.
    cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
    cte.set <- matrix(cte.des[1:n.alts, ], ncol = n.cte, byrow = FALSE)
  } else {cte.des = NULL}
  # if no alternative constants 
  if (!is.matrix(par.draws)) {
    stop("'par.draws'should be a matrix when 'alt.cte' = NULL")
  }
  #init
  n.par <- ncol(par.draws)
  if (!is.null(weights)) {
    if (!isTRUE(all.equal(length(weights), nrow(par.draws)))) {
      stop("length of 'weights' does not match number total number of rows in 'par.draws'")
    }
  }
  # If no weights, equal weights.
  if (is.null(weights)) {
    weights <- rep(1L, nrow(par.draws))
  }
  ## whenever a design is supplied 
  if (!is.null(des)) {
    # Error par.draws
    if (!isTRUE(all.equal(ncol(des), n.par))) {
      stop("Numbers of columns in 'par.draws' does not match the number of columns in 'des'")
    }
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    d.start <- apply(par.draws, 1, DerrC_ucpp, des = des,  n.alts = n.alts, 
                     i.cov = i.cov)
    db.start <- mean(d.start, na.rm = TRUE)
    full.comb <- Fullsets_ucpp(cand.set = cand.set, n.alts = n.alts, 
                               no.choice = no.choice, reduce = reduce, 
                               allow.rep = allow.rep, des = des, n.cte = n.cte)
    #if alt.cte
    if (!is.null(cte.des)) {
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    }
    full.des <- lapply(full.comb, function(x) rbind(des, x))
    # For each potential set, select best.
    ##### parallel #####
    if (parallel) {
      no_cores <- parallel::detectCores() - 1L
      cl <- parallel::makeCluster(no_cores)
      # New line to copy DerrS.P from .GlobalEnv to the cluster
      # https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function
      #clusterExport(cl=cl,varlist=c("DerrS.P_ucpp"))
      db.errors <- parallel::parLapply(cl, full.des, DBerrS.P_ucpp, par.draws, 
                                       n.alts, i.cov, weights)
      parallel::stopCluster(cl)
      ##### parallel #####
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.des, DBerrS.P_ucpp, par.draws, n.alts, i.cov, weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    db <- min(dbs)
    #return best set and db error design.
    return(list("set" = set, "error" = db))
  }
  
  if (is.null(des)) {
    
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    full.comb <- Fullsets_ucpp(cand.set = cand.set, n.alts = n.alts, 
                               no.choice = no.choice, reduce = reduce, 
                               allow.rep = allow.rep, des = des)
    #if alt.cte
    if (!is.null(cte.des)) {
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    } else {
      if (!isTRUE(all.equal(ncol(cand.set), ncol(par.draws)))) {
        stop("number of  columns of 'par.draws' and 'cand.set' should be equal when 'alt.cte = NULL'")
      }
    }
    # For each potential set, select best.
    
    ##### parallel #####
    if (parallel) {
      no_cores <- parallel::detectCores() - 1L
      cl <- parallel::makeCluster(no_cores)
      # New line to copy DerrS.P from .GlobalEnv to the cluster
      # https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function
      # Becareful with cpp functions in cluster export
      # Check https://stackoverflow.com/questions/38518387/using-rcpp-functions-inside-of-rs-parapply-functions-from-the-parallel-package
      # https://stackoverflow.com/questions/25606733/using-rcpp-function-in-parlapply-on-windows/25606950
      #clusterExport(cl=cl,varlist=c("DerrS.P_ucpp"))
      db.errors <- parallel::parLapply(cl, full.comb, DBerrS.P_ucpp, par.draws, 
                                       n.alts, i.cov, weights)
      parallel::stopCluster(cl)
      ##### parallel #####
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.comb, DBerrS.P_ucpp, par.draws, n.alts, i.cov,
                          weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    db <- min(dbs)
    #return best set and db error design.
    return(list("set" = set, "error" = db))
  }
}


#' Sequential Kullback-Leibler based algorithm for the MNL model.
#' 
#' Selects the choice set that maximizes the Kullback-Leibler divergence between
#' the prior parameter values and the expected posterior, assuming a MNL model.
#' 
#' This algorithm is ideally used in an adaptive context. The algorithm selects 
#' the choice set that maximizes the Kullback-Leibler 
#' divergence between prior and expected posterior. Otherwisely framed the 
#' algorithm selects the choice set that maximizes the expected information 
#' gain.
#' 
#' If \code{alt.cte = NULL}, \code{par.draws} should be a matrix in which each 
#' row is a sample from the multivariate parameter distribution. In case that 
#' \code{alt.cte} is not \code{NULL}, a list containing two matrices should be 
#' provided to \code{par.draws}. The first matrix containing the parameter draws
#' for the alternative specific parameters. The second matrix containing the
#' draws for the rest of the parameters.
#'
#' The list of potential choice sets are created using 
#' \code{\link[utils]{combn}}. The \code{weights} argument can be used when the
#'  \code{par.draws} have 
#' weights. This is for example the case when parameter values are updated using
#' \code{\link{ImpsampMNL}}.
#' 
#' 
#' @inheritParams SeqMOD
#' @param alt.cte A binary vector indicating for each alternative if an
#'   alternative specific constant is desired.
#' @return \item{set}{Numeric matrix containing the choice set that maximizes the expected KL divergence.}
#' \item{kl}{Numeric value which is the Kullback leibler divergence.}
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{crabbe}{idefix}
#' @examples 
#' # KL efficient choice set, given parameter draws. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3), coding = c("E", "E"))
#' m <- c(0.3, 0.2, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(m)) # Prior variance
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 draws.
#' # Efficient choice set to add. 
#' SeqKL(cand.set = cs, n.alts = 2, alt.cte = NULL, par.draws = ps, weights = NULL)
#' 
#' # KL efficient choice set, given parameter draws. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3), coding = c("C", "E"), c.lvls = list(c(5,3,1)))
#' m <- c(0.7, 0.3, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(m)) # Prior variance
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 draws.
#' sample <- list(ps[ , 1], ps[ , 2:4])
#' ac <- c(1, 0) # Alternative specific constant. 
#' # Efficient choice set to add. 
#' SeqKL(cand.set = cs, n.alts = 2, alt.cte = ac, par.draws = sample, weights = NULL)
#' @export
SeqKL <- function(des = NULL, cand.set, n.alts, par.draws, alt.cte = NULL, 
                  no.choice = NULL, weights = NULL, allow.rep = FALSE) {
  # Handling error initial design
  if (is.null(des)) {
    n.sets <- 1L
    allow.rep <- TRUE
  } else { 
    if (!is.matrix(des)) {
      stop("'des' should be a matrix or NULL")
    }
    if (!isTRUE(nrow(des) %% n.alts == 0)) {
      stop("'n.alts' does not seem to match with the number of rows in 'des'")
    }
    n.sets <- nrow(des) / n.alts
  }
  
  ### Error handling for design specifications
  # No choice errors
  if (!is.null(no.choice)) {
    if (!is.wholenumber(no.choice)) {
      stop("'no.choice' should be an integer or NULL")
    }
    if (any(isTRUE(no.choice > (n.alts + 0.2)), isTRUE(no.choice < 0.2))) {
      stop("'no.choice' does not indicate one of the alternatives")
    }
    if (is.null(alt.cte)) {
      stop("if there is a no choice alternative, 'alt.cte' should be specified")
    }
    if (!isTRUE(all.equal(alt.cte[no.choice], 1))) {
      stop("the no choice alternative should correspond with a 1 in 'alt.cte'")
    }
  }
  
  # Error alternative specific constants.
  if (!is.null(alt.cte)) {
    if (length(alt.cte) != n.alts) {
      stop("n.alts does not match the alt.cte vector")
    }
    if (!all(alt.cte %in% c(0, 1))) {
      stop("'alt.cte' should only contain zero or ones.")
    }
    # alternative specific constants
    n.cte <- length(which(alt.cte == 1L))
    if (isTRUE(all.equal(n.cte, 0L))) {
      alt.cte <- NULL
      cte.des <- NULL
    }
    
    # Handling errors when there are alternative constants
    if (n.cte > 0.2) {
      if (!is.list(par.draws)) {
        stop("'par.draws' should be a list when 'alt.cte' is not NULL")
      }
      if (!isTRUE(all.equal(length(par.draws), 2))) {
        stop("'par.draws' should contain two components")
      }
      # If there is only one specific constant and is a vector, then it is
      # transformed to a matrix
      if (isTRUE(all.equal(n.cte, 1))) {
        if (is.vector(par.draws[[1]])) {
          par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1)
        }
      }
      if (!(all(unlist(lapply(par.draws, is.matrix))))) {
        stop("'par.draws' should contain two matrices")
      }
      if (!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))) {
        stop("the first component of 'par.draws' should contain the same number 
             of columns as there are non zero elements in 'alt.cte'")
      }
      dims <-  as.data.frame(lapply(par.draws, dim))
      if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) { 
        stop("the number of rows in the components of 'par.draws' should be equal")
      }
      par.draws  <- do.call("cbind", par.draws) # Transform par.draws to a matrix
    }
    
    # Create alternative specific design.
    cte.des <- Altspec(alt.cte = alt.cte, n.sets = 1)  
    cte.set <- matrix(cte.des[1:n.alts, ], ncol = n.cte, byrow = FALSE)
    # Error handling cte.des
    if (ncol(cand.set) + ncol(cte.des) != ncol(par.draws)) {
      stop("dimension of par.draws does not match the dimension of alt.cte + cand.set.")
    }
  } else {
    n.cte <- 0
    cte.des <- NULL
    # Error handling cte.des
    if (ncol(cand.set) != ncol(par.draws)) {
      stop("dimension of par.draws does not match the dimension of alt.cte + cand.set.")
    }
  }
  
  # Check number of columns in par.draws and weights
  n.par <- ncol(par.draws)
  if (!is.null(weights)) {
    if (!isTRUE(all.equal(length(weights), nrow(par.draws)))) {
      stop("length of 'weights' does not match number total number of rows in 'par.draws'")
    }
  }
  # All choice sets.
  # full.comb <- gtools::combinations(n = nrow(cand.set), r = n.alts, 
  #                                   repeats.allowed = !reduce)
  full.comb <- Fullsets_ucpp(cand.set = cand.set, n.alts = n.alts, 
                                 no.choice = no.choice, reduce = FALSE, 
                                 allow.rep = allow.rep, des = des, n.cte = n.cte)
  
  # If no weights, equal weights.
  if (is.null(weights)) {
    weights <- rep(1, nrow(par.draws))
  }
  # Add alternative specific constants if necessary
  if (!is.null(cte.des)) {
    full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
  }
  # Calculate KL for each set.
  #kl.infos <- apply(full.comb, 1, KLs, par.draws, cte.des, cand.set, weights)
  kl.infos <- lapply(full.comb, KL, par.draws, weights)
  
  # Select maximum.
  #comb.nr <- as.numeric(full.comb[which.max(kl.infos), ])
  #set <- cand.set[comb.nr, ]
  set <- full.comb[[which.max(kl.infos)]]
  # Add alternative specific constants if necessary
  #if (!is.null(cte.des)) {
  #  set <- cbind(cte.des, set)
  #}
  row.names(set) <- NULL
  # return.
  return(list(set = set, kl = max(unlist(kl.infos))))
}
