

#' Modified Federov algorithm for MNL models.
#' 
#' The algorithm swipes every profile of an initial random design with candidate
#' profiles. By doing this it tries to minimize the D(B) error, assuming
#' a multinomial logit model. See reference for more information.  
#' 
#' The algorithm stops when an iteration occured without replacing a profile or 
#' when \code{max.iter} is reached. An iteration is a loop through all profiles 
#' from the initial design, testing every profile from the candidate design.
#' 
#' By specifying a numeric vector in \code{par.samples}, the D-error will be 
#' calculated and the design will be optimised locally. By specifying a matrix 
#' in which each row is a sample from a multivariate distribution in 
#' \code{par.samples}, the DB-error will be calculated, and the design will be 
#' optimised globally. The number of columns should equal the number of
#' parameters for alternative specific constants (\code{alt.cte}) + the number
#' of attribute parameters (\code{ncol(cand.des)}). This is also the order in
#' which they should appear (first \code{alt.cte} parameters).
#' 
#' Alternative specific constants can be specified in \code{alt.cte}. The lenght
#' of this binary vector should equal \code{n.alts}, were a zero indicates no
#' alternative specific constant and a one the opposit.
#' 
#' @param cand.set A numeric matrix in which each row is a possible profile. See
#'   \code{\link{Profiles}}
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice 
#'   set.
#' @param par.samples A matrix in which each row is a sample from the 
#'   multivariate prior parameter distribution.
#' @param alt.cte A binary vector indicating for each alternative if an alternative 
#'   specific constant is desired.
#' @param start.des A matrix in which each row is a profile. The number of rows
#'   equals \code{n.sets * n.alts}, and the number of columns equals the number
#'   of columns of \code{par.samples}
#' @param max.iter A numeric value indicating the maximum number allowed 
#'   iterations.
#' @return 
#' \item{design}{An efficient design.}
#' \item{error}{The D(B)-error of the design.}
#' \item{Na.percentage}{The percentage of
#'  samples for which the D-error was \code{Inf}.}
#' \item{prob.diff}{The difference between the alternative with
#' the highest and the one with the lowest probability for each choice set. If a
#' sample matrix was provided this is based on the mean over all samples.}
#' @export
Modfed <- function(cand.set, n.sets, n.alts, par.samples, alt.cte, start.des = NULL, max.iter = Inf) {
  # handling par.samples
  if (!(is.matrix(par.samples))) {
    par.samples <- matrix(par.samples, nrow = 1)
  }
  # error alternative specific constants 
  if (length(alt.cte) != n.alts) {
    stop("n.alts does not match the alt.cte vector")
  }
  # error identifying model
  if (n.sets < ncol(par.samples)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  # errors start design 
  if (!is.null(start.des)) {
    if (ncol(start.des) != ncol(cand.set))  {
      stop("Number of colums start design is different from number of columns candidate set.")
    }
    if (nrow(start.des) != (n.alts * n.sets)) {
      stop("Number of rows start design is different from number of sets times number of alternatives.")
    }
  }
  # create alternative specific design
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
  # error handling
  if (ncol(cand.set)+ ncol(cte.des) != ncol(par.samples)) {
    stop("dimension of par.samples does not match the dimension of cand.set + alt.cte")
  }
  # random start design
  if (!is.null(start.des)) {
    des <- start.des
  } else {
    r <- round(runif((n.sets*n.alts), 1, nrow(cand.set)))
    des <- data.matrix(cand.set[r, ])
  }
  # combine with alt.spec design
  des <- cbind(cte.des, des)
  # starting values 
  d.start <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
  db.start <- mean(d.start, na.rm = TRUE)
  converge <- FALSE
  change <- FALSE
  it <- 1
  n.samples <- nrow(par.samples)
  n.cte <- ncol(cte.des)
  n.par <- ncol(des)
  # start algorithm
  while (!converge & it <= max.iter) {
    it <- it + 1
    # show progress iteration
    pb <- txtProgressBar(min = 0, max = nrow(des), style = 3)     
    # save design before iteration
    iter.des <- des
    # for every row in the design
    for (r in 1:nrow(des)) {
      # switch with everey row in candidate set 
      for (c in 1:nrow(cand.set)) {
        des[r, (n.cte + 1) : n.par ] <- cand.set[c, ]
        # calculate d errors
        d.errors <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
        #db error 
        db <- mean(d.errors, na.rm = TRUE)
        # change if lower db error
        if (db < db.start) {
          best.row <- as.numeric(des[r, ])
          db.start <- db
          change <- TRUE
        }
        # show progress iteration
        setTxtProgressBar(pb, r)
      }
      # replace with best profile if change
      if (change) {
        des[r, ] <- best.row
      } else {
        des[r, ] <- iter.des[r, ]
      }
      # initialize variables again 
      change <- FALSE
      na.percentage <- 0
    }
    close(pb)  # print progress 
    converge <- isTRUE(all.equal(des, iter.des)) # convergence if no profile is swapped this iteration 
  }
  # calculate percentage NA values
  d.errors <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
  if (any(is.na(d.errors))) {
    na.percentage <- percent(sum(is.na(d.errors))/n.samples)
  } 
  # utility balance
  ub <- apply(par.samples, 1, Utbal, des = des,  n.alts = n.alts)
  ub <- .rowMeans(ub, m = n.sets, n = n.samples, na.rm = FALSE)
    # rownames design 
  des.names <- Rcnames(n.sets = n.sets, n.alts = n.alts, n.cte = n.cte, alt.cte = alt.cte)
  rownames(des) <- des.names[[1]]
  # colnames alternative specific constants 
  if (n.cte != 0) {
    colnames(des)[1:n.cte] <- des.names[[2]]
  }
  # return design, D(B)error, percentage NA's, utility balance. 
  return(list("design" = des, "error" =  db, "Na.percentage" = na.percentage, "prob.diff" = ub))
}


#' DB sequential set
#' 
#' Adds the most DB efficient choice set to a design, given (updated) parameter
#' values.
#' @param des A design matrix in which each row is a profile.
#' @param cand.set A numeric matrix in which each row is a possible profile.
#' @param n.alts Numeric value indicating the number of alternatives per choice
#'   set.
#' @param par.samples A matrix in which each row is a sample from the
#'   multivariate parameter distribution.
#' @param weights A vector containing the weights of the samples.
#' @param prior.covar Covariance matrix of the prior distribution.
#' @return The most db efficient choice set.
#' @export
seqfed.db <- function(des, cand.set, n.alts, par.samples, weights, prior.covar) {

  #start values
  db_best<-1000
  full_comb<- full_sets(cand.set, n.alts)

  #for each potential set:
  for (p in 1:nrow(cand.set)){

    set<-as.matrix(cand.set[as.numeric(full_comb[p, ]), ])

    #for each par draw calculate d-error:
    db_error<-0

    for (s in 1: nrow(par.samples)){

      info_s<-info_design(par = par.samples[s, ], des = set, n.alts = n.alts)
      info_d<-info_design(par= par.samples[s, ], des = des, n.alts = n.alts)
      inv_cov_prior<-solve(prior.covar)

      db_error<- db_error+(det(info_d + info_s + inv_cov_prior)*(-1/ncol(par.samples))*weights[s])

    }

    #select new set if smaller db error
    if (db_error< db_best){final_set<-set; db_best<- db_error}
  }

  return(final_set)

}

#' KL set selecting
#' 
#' Provides the set that maximizes the Kullback-Leibler divergence, given
#' parameter values.
#' @param lvls A vector which contains for each attribute, the number of levels.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice
#'   set.
#' @param par.samples A matrix in which each row is a sample.
#' @param weights A vector containing the weights of the samples.
#' @return Choice set that maximizes the expected KL divergence.
#' @export
KL_select <- function(lvls, n.sets, n.alts, par.samples, weights){

  #All choice sets, without same profile twice
  fp<-profiles(lvls)
  fcomb<-full_sets(cand = fp, n.alts = n.alts)
  rows<-apply(fcomb[, 1:ncol(fcomb)], 1, function(i) length(unique(i)) > ncol(fcomb)-1)
  fcomb<-fcomb[rows, ]

  #start value
  kl_start <- 0
  best_set<-0

  for ( s in 1:nrow(fcomb)){

    #take set
    set<-as.matrix(fp[as.numeric(fcomb[s, ]), ])

    #calculate for all sets the KLinfo.
    klinfo<-KL(set, par.samples, weights)

    #if better --> keep
    if (klinfo > kl_start){
      best_set<- set
      kl_start<-klinfo
    }
    print(klinfo)
  }
  return(list(best_set, kl_start))
}











