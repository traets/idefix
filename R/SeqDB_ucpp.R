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
#' provided. The first matrix containing the parameter draws for the 
#' alternative specific parameters. The second matrix containing the draws for
#' the rest of the parameters.
#' 
#' The list of potential choice sets are created using 
#' \code{\link[gtools]{combinations}}. If \code{reduce} is \code{TRUE}, 
#' \code{repeats.allowed = FALSE} and vice versa. Furthermore, the list of 
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
#' @inheritParams Modfed
#' @param par.draws A matrix or a list, dependend on \code{alt.cte}. 
#' @param des A design matrix in which each row is a profile. If alternative 
#'   specific constants are present, those should be included as the first 
#'   column(s) of the design. Can be generated with \code{\link{Modfed}}
#' @param prior.covar Covariance matrix of the prior distribution.
#' @param weights A vector containing the weights of the draws. Default is 
#'   \code{NULL}, See also \code{\link{ImpsampMNL}}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores.
#' @param reduce Logical value indicating whether the candidate set should be 
#'   reduced or not.
#' @return \item{set}{A matrix representing a DB efficient choice set.} 
#'   \item{db.error}{A numeric value indicating the DB-error of the whole 
#'   design.}
#' @importFrom Rdpack reprompt
#' @references \insertRef{ju}{mnldes}
#' @examples 
#' # DB efficient choice set, given a design and parameter draws. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3, 3), coding = c("E", "E", "E"))
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4) # Prior mean (total = 6 parameters).
#' pc <- diag(length(m)) # Prior variance
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 draws.
#' # Initial design.
#' des <- Modfed(cand.set = cs, n.sets = 6, n.alts = 2, alt.cte = c(0, 0), par.draws = ps)$design
#' # Efficient choice set to add. 
#' SeqDB(des = des, cand.set = cs, n.alts = 2, par.draws = ps, prior.covar = pc)
#' 
#' # DB efficient choice set, given parameter draws. 
#' # with alternative specific constants 
#' cs <- Profiles(lvls = c(3, 3), coding = c("E", "E"))
#' m <- c(0.7, 0.3, 0.8, -0.2, -1.2) # Prior mean (5 parameters).
#' pc <- diag(length(m)) # Prior variance
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 draws.
#' ps <- list(ps[ , 1], ps[ , 2:5])
#' ac <- c(1, 0) # Alternative specific constant. 
#' # Efficient choice set. 
#' SeqDB(cand.set = cs, n.alts = 2, par.draws = ps, alt.cte = ac, prior.covar = pc)
#' @export
SeqDB_ucpp <- function(des = NULL, cand.set, n.alts, par.draws, prior.covar, alt.cte = NULL, no.choice = NULL, weights = NULL, parallel = TRUE, reduce = TRUE) {
  #init
  if (is.null(des)) {
    n.sets <- 1L
  } else { 
    if (!isTRUE(nrow(des) %% n.alts == 0)) {
      stop("'n.alts' does not seem to match with the number of rows in 'des'")
    }
    n.sets <- nrow(des) / n.alts
  }
  # if alternative constants 
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
  #if no.choice
  if(!is.null(no.choice)){
    if(!no.choice %% 1 == 0){
      stop("'no.choice' should be an integer")
    }
    if(any(isTRUE(no.choice > (n.alts + 0.2)), isTRUE(no.choice < 0.2))){
      stop("'no.choice' does not indicate one of the alternatives")
    }
  }
  if(!is.null(alt.cte)){
    #prior.covar
    if(!isTRUE(all.equal(ncol(prior.covar), (ncol(cand.set) + n.cte)))){
      stop("number of columns of 'prior.covar' does not equal 
           the total number of parameters (including 'alt.cte')")
    }
    if(isTRUE(all.equal(n.cte, 1))){
      if(!is.list(par.draws)){stop("'par.draws' should be a list when 'alt.cte' is not NULL")}
      if (!isTRUE(all.equal(length(par.draws), 2))){
        stop("'par.draws' should contain two components")
      }
      if(is.vector(par.draws[[1]])){
        par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1)
      }
      if(!(all(unlist(lapply(par.draws, is.matrix))))){
        stop("'par.draws' should contain 2 matrices")
      }
      if(!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))){
        stop("the first component of 'par.draws' should contain the same number 
             of columns as there are non zero elements in 'alt.cte'")
      }
      dims <-  as.data.frame(lapply(par.draws, dim))
      if(!isTRUE(all.equal(dims[1, 1], dims[1, 2]))){ 
        stop("the number of rows in the components of 'par.draws' should be equal")
      }
      if(!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))){ 
        stop("the sum of the number of columns in the components of 'par.draws' 
             should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
      }
      par.draws  <- do.call("cbind", par.draws)
      }
    if(n.cte > 1.2){
      if(!(is.list(par.draws))){stop("'par.draws' should be a list when 'alt.cte' is not NULL")} 
      if (!isTRUE(all.equal(length(par.draws), 2))){
        stop("'par.draws' should contain two components")
      }
      if(!(all(unlist(lapply(par.draws, is.matrix))))){
        stop("'par.draws' should contain two matrices")
      }
      if(!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))){
        stop("the first component of 'par.draws' should contain the same number 
             of columns as there are non zero elements in 'alt.cte'")
      }
      dims <-  as.data.frame(lapply(par.draws, dim))
      if(!isTRUE(all.equal(dims[1, 1], dims[1, 2]))){ 
        stop("the number of rows in the components of 'par.draws' should be equal")
      }
      if(!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))){ 
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
  if(!is.matrix(par.draws)){
    stop("'par.draws'should be a matrix when 'alt.cte' = NULL")
  }
  #init
  n.par <- ncol(par.draws)
  if(!is.null(weights)){
    if(!isTRUE(all.equal(length(weights), nrow(par.draws)))){
      stop("length of 'weights' does not match number total number of rows in 'par.draws'")
    }
  }
  # If no weights, equal weights.
  if (is.null(weights)) {
    weights <- rep(1L, nrow(par.draws))
  }
  ## whenever a design is supplied 
  if(!is.null(des)){
    # Error par.draws
    if (!isTRUE(all.equal(ncol(des), n.par))) {
      stop("Numbers of columns in 'par.draws' does not match the number of columns in 'des'")
    }
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    d.start <- apply(par.draws, 1, DerrC_ucpp, des = des,  n.alts = n.alts, i.cov = i.cov)
    db.start <- mean(d.start, na.rm = TRUE)
    full.comb <- Fullsets_ucpp(cand.set = cand.set, n.alts = n.alts, no.choice = no.choice, reduce = reduce)
    #if alt.cte
    if(!is.null(cte.des)){
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    }
    full.des <- lapply(full.comb, function(x) rbind(des, x))
    # For each potential set, select best.
    ##### parallel #####
    if (parallel) {
      library(parallel) # This should not be here pag 42,59 r_pack wickham
      no_cores <- detectCores() - 1L
      cl <- makeCluster(no_cores)
      # New line to copy DerrS.P from .GlobalEnv to the cluster
      # https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function
      #clusterExport(cl=cl,varlist=c("DerrS.P_ucpp"))
      db.errors <- parLapply(cl, full.des, DBerrS.P_ucpp, par.draws, 
                             n.alts, i.cov, weights)
      stopCluster(cl)
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
    return(list(set = set, db.error = db))
  }
  
  if(is.null(des)){
    
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    full.comb <- Fullsets_ucpp(cand.set = cand.set, n.alts = n.alts, no.choice = no.choice, reduce = reduce)
    #if alt.cte
    if(!is.null(cte.des)){
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    } else {
      if(!isTRUE(all.equal(ncol(cand.set), ncol(par.draws)))){
        stop("number of  columns of 'par.draws' and 'cand.set' should be equal when 'alt.cte = NULL'")
      }
    }
    # For each potential set, select best.
    
    ##### parallel #####
    if (parallel) {
      library(parallel)
      no_cores <- detectCores() - 1L
      cl <- makeCluster(no_cores)
      # New line to copy DerrS.P from .GlobalEnv to the cluster
      # https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function
      # Becareful with cpp functions in cluster export
      # Check https://stackoverflow.com/questions/38518387/using-rcpp-functions-inside-of-rs-parapply-functions-from-the-parallel-package
      # https://stackoverflow.com/questions/25606733/using-rcpp-function-in-parlapply-on-windows/25606950
      #clusterExport(cl=cl,varlist=c("DerrS.P_ucpp"))
      db.errors <- parLapply(cl, full.comb, DBerrS.P_ucpp, par.draws, 
                             n.alts, i.cov, weights)
      stopCluster(cl)
      ##### parallel #####
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.comb, DBerrS.P_ucpp, par.draws, n.alts, i.cov, weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    db <- min(dbs)
    #return best set and db error design.
    return(list(set = set, db.error = db))
  }
    }