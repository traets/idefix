
#' Coordinate Exchange algorithm for MNL models.
#' 
#' The algorithm swaps every profile of an initial start design with candidate 
#' profiles. By doing this it tries to minimize the D(B)-error, based on a 
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
#' containing two matrices. The first matrix containing the parameter draws for
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
#' \code{start.des} is a list with one or several matrices. In each matrix each
#' row is a profile. The number of rows equals \code{n.sets * n.alts}, and the
#' number of columns equals the number of columns of \code{cand.set} + the
#' number of non-zero elements in \code{alt.cte}. If \code{start.des
#' = NULL}, \code{n.start} random start designs will be
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
#'   \code{\link{Profiles}} function can be used to generate this.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice 
#'   set.
#' @param alt.cte A binary vector indicating for each alternative whether an 
#'   alternative specific constant is desired. The default is \code{NULL}.
#' @param par.draws A matrix or a list, dependend on \code{alt.cte}.
#' @param no.choice A logical value indicating whether a no choice alternative 
#'   should be added to each choice set. The default is \code{NULL}.
#' @param start.des A list containing one or more matrices. The default is \code{NULL}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores. The default is \code{TRUE}.
#' @param max.iter A numeric value indicating the maximum number allowed 
#'   iterations. The default is \code{Inf}.
#' @param n.start A numeric value indicating the number of random start designs
#'   to use. The default is 12.
#' @param best A logical value indicating whether only the best design should be
#'   returned. The default is \code{TRUE}.
#' @return 
#'   If \code{best = TRUE} the design with the lowest D(B)-error. If \code{best 
#'   = FALSE}, the result of all (provided) start designs. \item{design}{A
#'   numeric matrix wich contains an efficient design.} \item{error}{Numeric
#'   value indicating the D(B)-error of the design.} \item{inf.error}{Numeric
#'   value indicating the percentage of draws for which the D-error was
#'   \code{Inf}.} \item{probs}{Numeric matrix containing the probabilities of
#'   each alternative in each choice set. If a sample matrix was provided in
#'   \code{par.draws}, this is the average over all draws.}
#' @examples
#' \donttest{
#' # DB-efficient designs
#' # 3 Attributes, all dummy coded. 1 alternative specific constant. = 7 parameters
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
#' @references \insertRef{federov}{idefix}
#' @export
CEA <- function(lvls, coding, c.lvls = NULL, n.sets, n.alts, par.draws, 
                alt.cte = NULL, no.choice = FALSE, start.des = NULL, 
                parallel = TRUE, max.iter = Inf, n.start = 12, best = TRUE) {
  
  ### Error handling for creating initial random design
  # Error lvls vector. At least 2 attributes
  if (length(lvls) < 2 || (!(is.numeric(lvls)))) {
    stop("lvls argument is incorrect.")
  }
  # Error correct coding types.
  codings.types <- c("E", "D", "C")
  if (!all(coding %in% codings.types) || (length(coding) != length(lvls))) {
    stop("coding argument is incorrect.")
  } 
  # Continuous attributes. 
  contins <-  which(coding == "C")
  n.contins <-  length(contins)
  # error continuous levels (Not a list)
  if (!is.null(c.lvls) && !is.list(c.lvls)) { 
    stop('c.lvls should be a list.')
  }
  # Error continuous specified and NULL.
  if (n.contins > 0 && is.null(c.lvls)) {
    stop("when 'coding' contains C, 'c.lvls' should be specified")
  }
  # Error continuous levels specification. 
  if (!is.null(c.lvls)) {
    if (length(c.lvls) != n.contins) {
      stop("length of 'c.lvls' does not match number of specified continuous attributes in 'coding'")
    }
    # Error c.lvls same number of levels. 
    if (!isTRUE(all.equal(lvls[contins], lengths(c.lvls)))) {
      stop("the number of levels provided in 'c.lvls' does not match the expected based on 'lvls'")
    }
  }
  
  ### Error handling for design specifications
  # If no alternative constant is given, create the variable as a vector of 0s
  if (is.null(alt.cte)) {
    alt.cte <- rep(0L, n.alts)
  }
  #init
  n.cte <- length(which(alt.cte == 1))
  
  # If only one draw is given, transform it to a matrix
  if (!is.list(par.draws)) {
    if (is.vector(par.draws)) {
      par.draws <- matrix(par.draws, nrow = 1)
    }
  }
  
  # Alternative constant errors
  if (length(alt.cte) != n.alts) {
    stop("'n.alts' does not match the 'alt.cte' vector")
  }
  if (!all(alt.cte %in% c(0, 1))) {
    stop("'alt.cte' should only contain zero or ones.")
  }
  
  # No choice errors
  if (!is.logical(no.choice)) {
    stop("'no.choice' should be TRUE or FALSE")
  }
  if (no.choice) {
    if (!isTRUE(all.equal(alt.cte[n.alts], 1))) {
      stop("if 'no.choice' is TRUE, the last alternative constant should equal 1.")
    }
    ncsek <- seq(n.alts, (n.sets * n.alts), n.alts)  # Rows in the design matrix
    # to be assigned as zeros in all attributes
  } else {
    ncsek <- NULL
  }
  
  # Handling par.draws with alternative specific contstants.
  # This conditional is when there is only one alternative constant
  if (isTRUE(all.equal(n.cte, 1))) {
    if (!(is.list(par.draws))) {
      stop("par.draws should be a list")
    }
    if (!isTRUE(all.equal(length(par.draws), 2))) {
      stop("'par.draws' should contain two components")
    }
    # If only draws for the constant are given in a vector, transform it to a
    # matrix
    if (is.vector(par.draws[[1]])) {
      par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1) 
    }
    if (!(all(unlist(lapply(par.draws, is.matrix))))) {
      stop("'par.draws' should contain two matrices")
    }
    dims <-  as.data.frame(lapply(par.draws, dim))
    if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) { 
      stop("the number of rows in the components of 'par.draws' should be equal")
    }
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # To Do: I have to compute the number of colums of the design to make this
    # conditional
    # if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
    #   stop("the sum of the number of columns in the components of 'par.draws'
    #        should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
    # }
    par.draws  <- do.call("cbind", par.draws) # Transform draws to a matrix
  } else {
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
    if(!isTRUE(all.equal(dims[1, 1], dims[1, 2]))){
      stop("the number of rows in the components of 'par.draws' should be equal")
    }
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # To Do: I have to compute the number of colums of the design to make this
    # conditional
    # if(!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))){ 
    #   stop("the sum of the number of columns in the components of 'par.draws' 
    #        should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
    # }
    par.draws  <- do.call("cbind", par.draws) # Transform draws to a matrix
  }
  
  # Error identifying model.
  if (n.sets < ncol(par.draws)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  
  ### Create all levels for each attribute
  # Change into correct coding. 
  coding <- dplyr::recode(coding, D = "contr.treatment", E = "contr.sum")
  # Create all combinations of attribute levels.
  levels.list <- lapply(X = as.list(lvls), function(x) (1:x))
  # Replace continuous.
  levels.list[contins] <- c.lvls
  # Transform categorical attributes
  categ <-  which(coding %in% c("contr.treatment", "contr.sum"))
  n.categ <- length(categ)
  if (n.categ > 0) {
    levels.list[categ] <- lapply(X = levels.list[categ], factor)
    # Apply coding
    if (n.contins > 0) {
      for (i in 1:length(lvls)) {
        if (!(i %in% contins)) {
          contrasts(levels.list[[i]]) <- coding[i]
        }
      }
    }else {
      for (i in 1:length(lvls)) {
        contrasts(levels.list[[i]]) <- coding[i]
      }
    }
    # Compute all possible values for each categorical attribute
    levels.list[categ] <- lapply(X = levels.list[categ], contrasts)
  }
  
  # Create alternative specific design.
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
  
  
  return(levels.list)
}


#   if(is.null(alt.cte)){
#     alt.cte <- rep(0L, n.alts)
#   }
#   #init
#   n.cte <- length(which(alt.cte == 1))
#   ### Errors
#   if(!is.list(par.draws)){
#     if(is.vector(par.draws)){
#       par.draws <- matrix(par.draws, nrow = 1)
#     }
#   }
#   #handling alt.cte
#   if (length(alt.cte) != n.alts) {
#     stop("'n.alts' does not match the 'alt.cte' vector")
#   }
#   if (!all(alt.cte %in% c(0, 1))){
#     stop("'alt.cte' should only contain zero or ones.")
#   }
#   #if no.choice
#   if(!is.logical(no.choice)){
#     stop("'no.choice' should be TRUE or FALSE")
#   }
#   if(no.choice){
#     if(!isTRUE(all.equal(alt.cte[n.alts], 1))){
#       stop("if 'no.choice' is TRUE, alt.cte[n.alts] should equal 1.")
#     }
#     ncsek <- seq(n.alts, (n.sets * n.alts), n.alts) 
#   } else {
#     ncsek <- NULL
#   }
#   # Handling par.draws with alternative specific contstants.
#   if(isTRUE(all.equal(n.cte, 1))){
#     if(!(is.list(par.draws))){stop("par.draws should be a list")}
#     if (!isTRUE(all.equal(length(par.draws), 2))){
#       stop("'par.draws' should contain two components")
#     }
#     if(is.vector(par.draws[[1]])){
#       par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1)
#     }
#     if(!(all(unlist(lapply(par.draws, is.matrix))))){
#       stop("'par.draws' should contain two matrices")
#     }
#     if(!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))){
#       stop("the first component of 'par.draws' should contain the same number 
#            of columns as there are non zero elements in 'alt.cte'")
#     }
#     dims <-  as.data.frame(lapply(par.draws, dim))
#     if(!isTRUE(all.equal(dims[1, 1], dims[1, 2]))){ 
#       stop("the number of rows in the components of 'par.draws' should be equal")
#     }
#     if(!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))){ 
#       stop("the sum of the number of columns in the components of 'par.draws' 
#            should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
#     }
#     par.draws  <- do.call("cbind", par.draws)
#     }
#   if(n.cte > 1.2){
#     if(!(is.list(par.draws))){stop("par.draws should be a list")} 
#     if (!isTRUE(all.equal(length(par.draws), 2))){
#       stop("'par.draws' should contain two components")
#     }
#     if(!(all(unlist(lapply(par.draws, is.matrix))))){
#       stop("'par.draws' should contain two matrices")
#     }
#     if(!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))){
#       stop("the first component of 'par.draws' should contain the same number 
#            of columns as there are non zero elements in 'alt.cte'")
#     }
#     dims <-  as.data.frame(lapply(par.draws, dim))
#     if(!isTRUE(all.equal(dims[1, 1], dims[1, 2]))){ 
#       stop("the number of rows in the components of 'par.draws' should be equal")
#     }
#     if(!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))){ 
#       stop("the sum of the number of columns in the components of 'par.draws' 
#            should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
#     }
#     par.draws  <- do.call("cbind", par.draws)
#     }
#   # Create alternative specific design.
#   cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
#   # Error identifying model.
#   if (n.sets < ncol(par.draws)) {
#     stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
#   }
#   # Handling cand.set
#   if(!all(is.finite(cand.set))){
#     stop("'cand.set' contains non finite values.")
#   }
#   # Error handling cte.des
#   if (ncol(cand.set) + ncol(cte.des) != ncol(par.draws)) {
#     stop("The number of parameters in the components of 'par.draws' does not match the number 
#          of non-zero parameters in 'alt.cte' + the number of parameters in 'cand.set'.")
#   }
#   # Random start design.
#   if (!is.null(start.des)) {
#     if(!is.list(start.des)){
#       stop("'start.des' should be a list")
#     }
#     if(!(all(unlist(lapply(start.des, is.matrix))))){
#       stop("'start.des' should contain matrices as components")
#     }
#     dimstart <- as.matrix(lapply(start.des, dim))
#     nr.starts <- length(dimstart)
#     if(nr.starts > 1.5){
#       if(!isTRUE(all.equal(length(unique(unlist(dimstart))), 2))){
#         stop("start designs have different dimensions")
#       }
#     }
#     if(!isTRUE(all.equal(n.alts * n.sets, unique(unlist(dimstart))[1]))){
#       stop("number of rows of start design(s) does not match with 'n.alts' * 'n.sets'")
#     }
#     if(!isTRUE(all.equal(sum(ncol(cand.set), ncol(cte.des)), unique(unlist(dimstart))[2]))){
#       stop("number of columns of start design(s) does not match with the number
#            of columns in 'cand.set' + the non zero parameters in 'alt.cte'")
#     }
#     d.start <- lapply(start.des, StartDB, par.draws, n.alts)
#     if(!any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))) {
#       stop("One or more of the provided start designs resulted in an unvalid db-error.")
#     }
#     } 
#   if (is.null(start.des)) {
#     #create start designs
#     nr.starts <- n.start
#     start.des <- vector(mode = 'list', length = nr.starts)
#     okstart <- FALSE
#     while(okstart == FALSE){
#       for (i in 1:nr.starts){
#         r <- round(stats::runif((n.sets * n.alts), 1, nrow(cand.set)))
#         start.des[[i]] <- cbind(cte.des, data.matrix(cand.set[r, ]))
#         if(no.choice){
#           start.des[[i]][ncsek, (ncol(cte.des) + 1) : (ncol(cte.des) + ncol(cand.set))] <- c(rep(0, ncol(cand.set)))
#         }
#       }
#       d.start <- lapply(start.des, StartDB, par.draws, n.alts)
#       if(any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))){
#         okstart <- TRUE
#       } 
#     }
#   }
#   if (parallel){
#     ########
#     no_cores <- parallel::detectCores() - 1
#     cl <- parallel::makeCluster(no_cores)
#     parallel::clusterExport(cl, c("n.sets", "par.draws", "cand.set", "n.alts", "n.cte", "alt.cte", "no.choice", "max.iter","ncsek"), envir = environment())
#     deslist <- parallel::parLapply(cl, start.des, Modfedje_ucpp, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte, no.choice, max.iter, ncsek)
#     parallel::stopCluster(cl)
#     ########
#   } else {
#     deslist <- lapply(start.des, Modfedje_ucpp, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte, no.choice, max.iter = max.iter, ncsek)
#   }                                 
#   bestdes <- deslist[[which.min(unlist(lapply(deslist, function(x) (x$error))))]]
#   
#   ifelse(best, return(bestdes), return(deslist))
#   }
# 
# # Core of the Modfed algorithm
# Modfedje_ucpp <- function(desje, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte,
#                           no.choice, max.iter, ncsek){
#   converge <- FALSE
#   change <- FALSE
#   it <- 1
#   n.samples <- nrow(par.draws)
#   n.par <- ncol(desje)
#   ###
#   while (!converge & it <= max.iter) {
#     db.start <- mean(apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts), na.rm = TRUE)
#     it <- it + 1
#     # save design before iteration.
#     iter.des <- desje
#     # For every row in the design.
#     sek <- 1 : nrow(desje)
#     if (no.choice){
#       sek <- sek[-ncsek]
#     }
#     for (r in sek) {
#       # Switch with everey row in candidate set. 
#       db <- numeric(nrow(cand.set))
#       for (c in 1:nrow(cand.set)) {
#         desje[r, (n.cte + 1) : n.par ] <- cand.set[c, ]
#         # Calculate D-errors.
#         d.errors <- apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts)
#         # DB-error. 
#         db[c] <- mean(d.errors, na.rm = TRUE)
#       }
#       pr <- which.min(db)
#       db <- min(db) 
#       # Change if lower db error.
#       if (!is.na(db) && !is.na(db.start)) {
#         if (db < db.start) {
#           best.row <- as.numeric(cand.set[pr, ])
#           db.start <- db
#           change <- TRUE
#         }
#       }
#       # Replace with best profile if change.
#       if (change) {
#         desje[r, (n.cte + 1) : n.par] <- best.row
#       } else {
#         desje[r, ] <- iter.des[r, ]
#       }
#       # Initialize variables again. 
#       change <- FALSE
#       na.percentage <- 0
#     }
#     converge <- isTRUE(all.equal(desje, iter.des)) # Convergence if no profile is swapped this iteration.
#   }
#   # calculate percentage NA values.
#   d.errors <- apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts)
#   if (any(is.na(d.errors))) {
#     na.percentage <- scales::percent(sum(is.na(d.errors)) / n.samples)
#   } 
#   # Utility balance.
#   ub <- apply(par.draws, 1, Utbal, des = desje,  n.alts = n.alts)
#   pmat <- matrix(rowMeans(ub), ncol = n.alts, byrow = TRUE)
#   rownames(pmat) <- paste("set", 1:n.sets, sep = "")
#   colnames(pmat) <- paste(paste("Pr(", paste("alt", 1:n.alts, sep = ""), sep = ""), ")", sep= "")
#   if(no.choice){
#     colnames(pmat)[n.alts] <- "Pr(no choice)"
#   }
#   # Rownames design. 
#   des.names <- Rcnames(n.sets = n.sets, n.alts = n.alts, alt.cte = alt.cte, no.choice = no.choice)
#   rownames(desje) <- des.names[[1]]
#   # Colnames alternative specific constants. 
#   if (n.cte != 0 && !is.null(colnames(desje))) {
#     colnames(desje)[1:n.cte] <- des.names[[2]]
#   }
#   # Return design, D(B)error, percentage NA's, utility balance. 
#   return(list("design" = desje, "error" =  db.start, "inf.error" = na.percentage, "probs" = pmat))
# }