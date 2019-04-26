
#' Coordinate Exchange algorithm for MNL models.
#' 
#' The algorithm improves an initial start design by considering changes on an
#' attribute-by-attribute basis. By doing this it tries to minimize the 
#' D(B)-error, based on a multinomial logit model. This routine is repeated for
#' multiple starting designs.
#' 
#' Each iteration will loop through all profiles from the initial design, 
#' evaluating the change in D(B)-error for every level in each attribute.
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
#' number of columns equals the number of columns of the design matrix + the
#' number of non-zero elements in \code{alt.cte}. Consider that for a 
#' categorical attribute with *p* levels, there are *p - 1* columns in the design
#' matrix, whereas for a continuous attribute there is only one column. If
#' \code{start.des = NULL}, \code{n.start} random start designs will be
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
#' @param lvls  A numeric vector which contains for each attribute, the number
#'   of levels.
#' @param coding Type op coding that needs to be used for each attribute.
#' @param c.lvls A list containing numeric vectors with the attribute levels for
#'   each continuous attribute. The default is \code{NULL}.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice 
#'   set.
#' @param par.draws A matrix or a list, dependend on \code{alt.cte}.
#' @param alt.cte A binary vector indicating for each alternative whether an 
#'   alternative specific constant is desired. The default is \code{NULL}.
#' @param no.choice A logical value indicating whether a no choice alternative 
#'   should be added to each choice set. The default is \code{FALSE}.
#' @param start.des A list containing one or more matrices. The default is 
#' \code{NULL}.
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
#' mu <- c(1.2, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' p.d <- list(matrix(pd[,1], ncol = 1), pd[,2:7])
#' CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = p.d,
#' n.alts = 2, n.sets = 8, parallel = FALSE, alt.cte = c(0, 1))
#' 
#' # DB-efficient design with categorical and continuous factors
#' # 2 categorical attributes with 4 and 2 levels (effect coded) and 1 
#' # continuous attribute (= 5 parameters)
#' mu <- c(0.5, 0.8, 0.2, 0.4, 0.3) 
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 3, mu = mu, Sigma = v) # 10 draws.
#' CEA(lvls = c(4, 2, 3), coding = c("E", "E", "C"), par.draws = pd,
#' c.lvls = list(c(2, 4, 6)), n.alts = 2, n.sets = 6, parallel = F)
#' 
#' # DB-efficient design with start design provided.  
#' # 3 Attributes with 3 levels, all dummy coded (= 6 parameters).
#' mu <- c(0.8, 0.2, -0.3, -0.2, 0.7, 0.4) 
#' v <- diag(length(mu)) # Prior variance.
#' sd <- list(example_design)
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = ps,
#' n.alts = 2, n.sets = 8, parallel = F, start.des = sd)
#'}
#' @importFrom Rdpack reprompt
#' @references \insertRef{cea}{idefix}
#' @references \insertRef{cea_discrete}{idefix}
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
    par.draws  <- do.call("cbind", par.draws) # Transform draws to a matrix
  } else {
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
        stop("the first component of 'par.draws' should contain the same number of columns as there are non zero elements in 'alt.cte'")
      }
      dims <-  as.data.frame(lapply(par.draws, dim))
      if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) {
        stop("the number of rows in the components of 'par.draws' should be equal")
      }
      par.draws  <- do.call("cbind", par.draws) # Transform draws to a matrix
    }
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
  # Transform to matrix
  levels.list[contins] <- lapply(X = levels.list[contins], as.matrix)
  # Transform categorical attributes
  categ <-  which(coding %in% c("contr.treatment", "contr.sum"))
  n.categ <- length(categ)
  if (n.categ > 0) {
    levels.list[categ] <- lapply(X = levels.list[categ], factor)
    # Apply coding
    if (n.contins > 0) {
      for (i in 1:length(lvls)) {
        if (!(i %in% contins)) {
          stats::contrasts(levels.list[[i]]) <- coding[i]
        }
      }
    }else {
      for (i in 1:length(lvls)) {
        stats::contrasts(levels.list[[i]]) <- coding[i]
      }
    }
    # Compute all possible values for each categorical attribute
    levels.list[categ] <- lapply(X = levels.list[categ], stats::contrasts)
  }
  # Set colnames for the design matrix
  c.nam = list()
  for (i in 1:length(lvls)) {
    if (coding[i] == "contr.treatment") {
      c.nam[[i]] <- paste("Var", i, 2:lvls[i], sep = "")
    } else {
      if (coding[i] == "contr.sum") {
        c.nam[[i]] <- paste("Var", i, 1:(lvls[i] - 1), sep = "")
      } else {
        c.nam[[i]] <- paste("Var", i, sep = "")
      }
    }
  }
  c.nam = unlist(c.nam)
  
  # Create alternative specific design.
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
  
  # Count the number of colums of the design matrix without constants
  ncol.des.noconst <- sum(unlist(lapply(levels.list,ncol)))
  #ncol.des.noconst <- sum(unlist(lapply(levels.list,function(x){
  #  if (is.matrix(x)) { ncol(x) }
  #  else{if (is.numeric(x)) { 1 }}
  #})))
  
  if (!identical(as.integer(ncol(par.draws)), 
                 as.integer(n.cte + ncol.des.noconst))) {
    stop("The sum of the number of columns in the components of 'par.draws' should equal the number of columns of design matrix (including alternative specific constants)")
  }
  
  ### Random initial design.
  if (is.null(start.des)) {
    #create start designs
    nr.starts <- n.start
    start.des <- vector(mode = 'list', length = nr.starts)
    okstart <- FALSE
    while (okstart == FALSE) {
      for (i in 1:nr.starts) {
        # r is to know which levels to take in each attribute
        r <- NULL
        start <- NULL
        for (j in 1:length(lvls)) {
          r <- round(stats::runif((n.sets * n.alts), 1, lvls[j]))
          start <- cbind(start, levels.list[[j]][r,])
        }
        colnames(start) <- c.nam
        start.des[[i]] <- cbind(cte.des, start)
        if (no.choice) {
          start.des[[i]][ncsek, (ncol(cte.des) + 1):(ncol(cte.des) + ncol.des.noconst)] <- c(rep(0, ncol.des.noconst))
        }
      }
      # Compute D-optimality for each design and each draw
      d.start <- lapply(start.des, StartDB, par.draws, n.alts)
      # If the DB-optimality of any starting desing is finite, continue
      if (any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))) {
        okstart <- TRUE
      } 
    }
  } else {
    if (!is.list(start.des)) {
      stop("'start.des' should be a list")
    }
    if (!(all(unlist(lapply(start.des, is.matrix))))) {
      stop("'start.des' should contain matrices as components")
    }
    # Save the dimension of each starting design
    dimstart <- as.matrix(lapply(start.des, dim)) 
    # Save the number of random starts given
    nr.starts <- length(dimstart)
    if (nr.starts > 1.5) {
      if (!isTRUE(all.equal(length(unique(unlist(dimstart))), 2))) {
        stop("start designs have different dimensions")
      }
    }
    if (!isTRUE(all.equal(n.alts * n.sets, unique(unlist(dimstart))[1]))) {
      stop("number of rows of start design(s) does not match with 'n.alts' * 'n.sets'")
    }
    if (!isTRUE(all.equal(as.integer(n.cte + ncol.des.noconst), 
                         unique(unlist(dimstart))[2]))) {
      stop("number of columns of start design(s) does not match with the number of columns in the design matrix")
    }
    # Compute D-optimality for each design and each draw
    d.start <- lapply(start.des, StartDB, par.draws, n.alts)
    if (!any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))) {
      stop("One or more of the provided start designs resulted in an unvalid db-error.")
    }
  }
  
  ### Improving the initial design
  if (parallel) {
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_cores)
    parallel::clusterExport(cl, c("n.sets", "par.draws", "n.alts", "n.cte", 
                                  "alt.cte", "no.choice", "max.iter","ncsek"),
                            envir = environment())
    deslist <- parallel::parLapply(cl, start.des, CEAcore_ucpp, par.draws, 
                                   levels.list, n.alts, n.sets, n.cte, alt.cte, 
                                   no.choice, max.iter, ncsek)
    parallel::stopCluster(cl)
  } else {
      deslist <- lapply(start.des, CEAcore_ucpp, par.draws, levels.list, 
                        n.alts, n.sets, n.cte, alt.cte, no.choice, 
                        max.iter = max.iter, ncsek)
  }                                 
  bestdes <- deslist[[which.min(unlist(lapply(deslist, function(x) (x$error))))]]
  
  ifelse(best, return(bestdes), return(deslist))
}

# Core of Coordinate Exchange
CEAcore_ucpp <- function(des, par.draws, levels.list, n.alts, n.sets, n.cte, 
                         alt.cte, no.choice, max.iter, ncsek) {
  converge <- FALSE # Boolean for convergence
  change <- FALSE # Boolean for a change in the attribute
  it <- 1 # Indicator of number of iterations
  n.samples <- nrow(par.draws)  # Number of samples from distribution of betas
  n.par <- ncol(des) # Number of parameters
  
  # Beginning of algorithm
  while (!converge & it <= max.iter) {
    # Compute DB-optimality of initial design
    db.start <- mean(apply(par.draws, 1, Derr_ucpp, des = des, 
                           n.alts = n.alts), na.rm = TRUE)
    it <- it + 1 # Increase iteration
    # Save design before iteration.
    iter.des <- des
    # For every row in the design.
    sek <- 1:nrow(des)
    if (no.choice) {
      sek <- sek[-ncsek] # If no.choice is given, then it is not improved
    }
    # Loop for each row of the initial design
    for (i in sek) {
      # Loop for each attribute in the row
      for (j in 1:length(levels.list)) {
        #%%%%%%%
        # This can be improved by removing mods object and just replace the 
        # initial design with the new attributes to compute the db.error and 
        # then if one of those is better than the initial db error. Replace 
        # finally the initial design (as it is in modfed)
        # Initialize mods object to track changes of the design in each attribute
        mods <- vector("list", nrow(levels.list[[j]])) 
        mods <- lapply(mods,function(x){des})
        # Initialize db object for db error of each level
        db <- numeric(nrow(levels.list[[j]]))
        # Indicator of columns modified
        ncol.ok <- sum(unlist(lapply(levels.list[(1:j - 1)],ncol)))
        # Indicator of columns to change
        ncol.mod <- ifelse( j == 1, (n.cte + j), 
                            (n.cte + 1 + ncol.ok))
        # ncol.left tracks which columns are left to be improved in the row
        ncol.left <- sum(unlist(lapply(levels.list[-(1:j)],ncol)))
        # Loop for each level in the attribute
        for (k in 1:nrow(levels.list[[j]])) {
          # Update design with new attribute
          #mods[[k]][i, (n.cte + j):(n.par - ncol.left) ] <- levels.list[[j]][k,]
          mods[[k]][i, ncol.mod:(n.par - ncol.left) ] <- levels.list[[j]][k,]
          # Calculate D-optimality for each draw
          d.errors <- apply(par.draws, 1, Derr_ucpp, des = mods[[k]], 
                            n.alts = n.alts)
          # Compute DB-optimality 
          db[k] <- mean(d.errors, na.rm = TRUE)
        }
        pr <- which.min(db)
        db <- min(db) 
        # Change if lower db error.
        if (!is.na(db) && !is.na(db.start)) {
          if (db < db.start) {
            des <- mods[[pr]]
            db.start <- db
          }
        }
      } # End loop for each attribute in the row
    } # End loop for each row of the initial design
    converge <- isTRUE(all.equal(des, iter.des)) # Convergence if no profile is swapped this iteration.
  } # End while (after convergence)
  #return(it, des)
  # calculate percentage NA values.
  na.percentage <- 0
  d.errors <- apply(par.draws, 1, Derr_ucpp, des = des,  n.alts = n.alts)
  if (any(is.na(d.errors))) {
    na.percentage <- scales::percent(sum(is.na(d.errors)) / n.samples)
  } 
  # Utility balance.
  # ub <- apply(par.draws, 1, Utbal, des = des,  n.alts = n.alts)
  # ub2 <- apply(par.draws, 1, InfoDes2, des = des,  n.alts = n.alts, utbal = T)
  # ub_ucpp <- apply(par.draws, 1, InfoDes_cpp, des = des,  n_alts = n.alts, 
  #                  utbal = T)
  
  # Utility balance using c++ function
  ub <- apply(par.draws, 1, InfoDes_cpp, des = des,  n_alts = n.alts, 
                   utbal = T)
  pmat <- matrix(rowMeans(ub), ncol = n.alts, byrow = TRUE)
  rownames(pmat) <- paste("set", 1:n.sets, sep = "")
  colnames(pmat) <- paste(paste("Pr(", paste("alt", 1:n.alts, sep = ""), 
                                sep = ""), ")", sep = "")
  if (no.choice) {
    colnames(pmat)[n.alts] <- "Pr(no choice)"
  }
  # Rownames design. 
  des.names <- Rcnames(n.sets = n.sets, n.alts = n.alts, alt.cte = alt.cte, 
                       no.choice = no.choice)
  rownames(des) <- des.names[[1]]
  # Colnames alternative specific constants. 
  if (n.cte != 0 && !is.null(colnames(des))) {
    colnames(des)[1:n.cte] <- des.names[[2]]
  }
  # Return design, D(B)error, percentage NA's, utility balance. 
  return(list("design" = des, "error" =  db.start, "inf.error" = na.percentage,
              "probs" = pmat))
}


#' Sequential coordinate exchange algorithm for MNL model.
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
#' The list of potential choice sets is created by selecting randomly a level for
#' each attribute in an alternative/profile. \code{n.cs} controls the number of
#' potential choice sets to consider. The default is \code{
#' NULL}, which means that the number of possible choice sets is the product of
#' attribute levels considered in the experiment. For instance, an experiment 
#' with 3 attribute and 3 levels each will consider 3^3 = 27 possible choice sets. 
#' 
#' If \code{reduce} is \code{TRUE}, \code{repeats.allowed = FALSE} and vice versa.
#' Furthermore, the list of potential choice sets will be screaned in order to
#' select only those choice sets with a unique information matrix. If no 
#' alternative specific constants 
#' are used, \code{reduce} should always be \code{TRUE}. When alternative specific 
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
#' @inheritParams CEA
#' @param par.draws A matrix or a list, depending on \code{alt.cte}. 
#' @param des A design matrix in which each row is a profile. If alternative 
#'   specific constants are present, those should be included as the first 
#'   column(s) of the design. Can be generated with \code{\link{Modfed}} or
#'   \code{\link{CEA}}
#' @param n.cs An integer indicating the number of possible random choice sets to 
#' consider in the search for the next best choice set possible. The default is
#'  \code{NULL}. 
#' @param prior.covar Covariance matrix of the prior distribution.
#' @param weights A vector containing the weights of the draws. Default is 
#'   \code{NULL}, See also \code{\link{ImpsampMNL}}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores.
#' @param no.choice An integer indicating the no choice alternative. The default
#'   is \code{NULL}.
#' @param reduce Logical value indicating whether the candidate set should be 
#'   reduced or not.
#' @return \item{set}{A matrix representing a DB efficient choice set.} 
#'   \item{error}{A numeric value indicating the DB-error of the whole 
#'   design.}
#' @importFrom Rdpack reprompt
#' @references \insertRef{ju}{idefix}
#' @references \insertRef{cea}{idefix}
#' @references \insertRef{cea_discrete}{idefix}
#' @examples 
#' # DB efficient choice set, given a design and parameter draws. 
#' # 3 attributes with 3 levels each
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4) # mean (total = 6 parameters).
#' pc <- diag(length(m)) # covariance matrix
#' set.seed(123)
#' sample <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
#' # Initial design.
#' des <- example_design
#' # Efficient choice set to add.
#' SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), n.alts = 2,
#'        par.draws = sample, prior.covar = pc, parallel = FALSE)
#' 
#' # DB efficient choice set, given parameter draws. 
#' # with alternative specific constants 
#' des <- example_design2
#' ac <- c(1, 1, 0) # Alternative specific constants.
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4, 1.8, 1.2) # mean
#' pc <- diag(length(m)) # covariance matrix
#' pos <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
#' sample <- list(pos[ , 1:2], pos[ , 3:8])
#' # Efficient choice set.
#' SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), n.alts = 3, 
#'       par.draws = sample, alt.cte = ac, prior.covar = pc, parallel = FALSE)
#' @export
SeqCEA <- function(des = NULL, lvls, coding, c.lvls = NULL, n.alts, par.draws, 
                   prior.covar, n.cs = NULL, alt.cte = NULL, no.choice = NULL,
                   weights = NULL, parallel = TRUE, reduce = TRUE) {
  # Error handling initial design
  if (is.null(des)) {
    n.sets <- 1L
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
  
  # Alternative constant errors
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
    cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
    cte.set <- matrix(cte.des[1:n.alts, ], ncol = n.cte, byrow = FALSE)
  } else {
      cte.des <- NULL
      n.cte <- 0
      # if no alternative constants 
      if (!is.matrix(par.draws)) {
        stop("'par.draws'should be a matrix when 'alt.cte' = NULL")
      }
    }
  
  # Weights errors
  n.par <- ncol(par.draws)
  if (!is.null(weights)) {
    if (!isTRUE(all.equal(length(weights), nrow(par.draws)))) {
      stop("length of 'weights' does not match number total number of rows in 'par.draws'")
    }
  } else {
    weights <- rep(1L, nrow(par.draws))
  }
  
  ### Create all levels for each attribute
  # Continuous attributes. 
  contins <-  which(coding == "C")
  n.contins <-  length(contins)
  # Change into correct coding. 
  coding <- dplyr::recode(coding, D = "contr.treatment", E = "contr.sum")
  # Create all combinations of attribute levels.
  levels.list <- lapply(X = as.list(lvls), function(x) (1:x))
  # Replace continuous.
  levels.list[contins] <- c.lvls
  # Transform to matrix
  levels.list[contins] <- lapply(X = levels.list[contins], as.matrix)
  # Transform categorical attributes
  categ <-  which(coding %in% c("contr.treatment", "contr.sum"))
  n.categ <- length(categ)
  if (n.categ > 0) {
    levels.list[categ] <- lapply(X = levels.list[categ], factor)
    # Apply coding
    if (n.contins > 0) {
      for (i in 1:length(lvls)) {
        if (!(i %in% contins)) {
          stats::contrasts(levels.list[[i]]) <- coding[i]
        }
      }
    }else {
      for (i in 1:length(lvls)) {
        stats::contrasts(levels.list[[i]]) <- coding[i]
      }
    }
    # Compute all possible values for each categorical attribute
    levels.list[categ] <- lapply(X = levels.list[categ], stats::contrasts)
  }
  # Set colnames for the design matrix
  c.nam = list()
  for (i in 1:length(lvls)) {
    if (coding[i] == "contr.treatment") {
      c.nam[[i]] <- paste("Var", i, 2:lvls[i], sep = "")
    } else {
      if (coding[i] == "contr.sum") {
        c.nam[[i]] <- paste("Var", i, 1:(lvls[i] - 1), sep = "")
      } else {
        c.nam[[i]] <- paste("Var", i, sep = "")
      }
    }
  }
  c.nam = unlist(c.nam)
  
  # Count the number of colums of the design matrix without constants
  ncol.des.noconst <- sum(unlist(lapply(levels.list,ncol)))
  
  # if (!identical(as.integer(n.par), as.integer(n.cte + ncol.des.noconst))) {
  #   stop("The sum of the number of columns in the components of 'par.draws' should equal the number of columns of design matrix (including alternative specific constants)")
  # }

  if (!identical(as.integer(ncol(prior.covar)), 
                            as.integer(n.cte + ncol.des.noconst))) {
    stop("number of columns of 'prior.covar' does not equal the number of columns of design matrix (including alternative specific constants)")
  }
  
  ## When a design is supplied
  if (!is.null(des)) {
    # Error par.draws
    if (!isTRUE(all.equal(ncol(des), n.par))) {
      stop("number of columns in 'par.draws' does not match the number of columns in 'des'")
    }
    # Error dimension of initial design and design matrix
    if (!identical(as.integer(ncol(des)), as.integer(n.cte + ncol.des.noconst))) {
      stop("number of columns in 'des' does not match the number of columns of design matrix (including alternative specific constants)")
      }
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    d.start <- apply(par.draws, 1, DerrC_ucpp, des = des, n.alts = n.alts, 
                     i.cov = i.cov)
    db.start <- mean(d.start, na.rm = TRUE)
    full.comb <- Newsets_ucpp(levels.list = levels.list, n.alts = n.alts, 
                              no.choice = no.choice, reduce = reduce, 
                              c.names = c.nam, n.cs = n.cs)
    # Adding alternative constants
    if (!is.null(cte.des)) {
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    }
    
    # Adding these new choice sets to the initial design
    full.des <- lapply(full.comb, function(x) rbind(des, x))
    
    # Select the next best choice set
    if (parallel) {
      no_cores <- parallel::detectCores() - 1L
      cl <- parallel::makeCluster(no_cores)
      db.errors <- parallel::parLapply(cl, full.des, DBerrS.P_ucpp, par.draws, 
                                       n.alts, i.cov, weights)
      parallel::stopCluster(cl)
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.des, DBerrS.P_ucpp, par.draws, n.alts, i.cov,
                          weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    colnames(set) <- colnames(des)
    db <- min(dbs)
    #return best set and db error design.
    return(list("set" = set, "error" = db))
  } # End if when a design is supplied
  else {
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    full.comb <- Newsets_ucpp(levels.list = levels.list, n.alts = n.alts, 
                              no.choice = no.choice, reduce = reduce, 
                              c.names = c.nam, n.cs = n.cs)
    # Adding alternative constants
    if (!is.null(cte.des)) {
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    }
    # Select the next best choice set
    if (parallel) {
      no_cores <- parallel::detectCores() - 1L
      cl <- parallel::makeCluster(no_cores)
      db.errors <- parallel::parLapply(cl, full.comb, DBerrS.P_ucpp, par.draws, 
                                       n.alts, i.cov, weights)
      parallel::stopCluster(cl)
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.comb, DBerrS.P_ucpp, par.draws, n.alts, i.cov,
                          weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    if (n.cte > 0) {
      colnames(set) <- c(paste("alt",1:n.cte,".cte", sep = ""), c.nam)
    }
    db <- min(dbs)
    #return best set and db error design.
    return(list("set" = set, "error" = db))
  }
}