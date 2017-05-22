

#' Modified federov algorithm.
#'
#' The algorithm swipes every profile of an initial random design with candidate profiles.
#' By doing this it tries to minimize the D error or DB error assuming a Multinomial logit model.  
#' 
#' The modified federov algorithm starts with creating an initial design by randomly selecting \code{n.alts* n.sets} profiles from the candidate set.
#' Afterwards the D error or DB error of this random initial design will be calculated. Then every profile from the initial design
#' will be replaced with every profile from the candidate set, keeping the one that minimizes the error.
#' 
#' By specifying a numeric vector in \code{par.samples}, 
#' the D error will be calculated and the design will be optimised locally. 
#' By specifying a matrix in which each row is a sample from a multivariate distribution in \code{par.samples}, 
#' the mean over D-errors will be taken as the efficiency criterion (DB-error), and the design will be optimised globally. 
#' 
#'
#' @param cand.set A numeric matrix in which each row is a possible profile.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice set.
#' @param par.samples A matrix in which each row is a sample from the multivariate prior parameter distribution.
#' @param max.iter A numeric value indicating the maximum number allowed iterations.
#' @return An efficient design, and the associated db-error.
#' @export
Modfed <- function(cand.set, n.sets, n.alts, par.samples, max.iter = Inf) {
  # error handling
  if (ncol(cand.set) != ncol(par.samples)) {
    stop("Number of paramters does not match the candidate set")
  }
  if (n.sets < ncol(par.samples)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  # random start design
  r <- round(runif((n.sets*n.alts), 1, nrow(cand.set)))
  des <- data.matrix(cand.set[r, ])
  # starting values 
  d.start <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
  db.start <- mean(d.start, na.rm = TRUE)
  converge <- FALSE
  change <- FALSE
  it <- 1
  n.samples <- nrow(par.samples)
  # start algorithm
  while (!converge & it <= max.iter) {
    it <- it + 1
    #show progress iteration
    pb <- txtProgressBar(min = 0, max = nrow(des), style = 3)
    #save design before iteration
    iter.des <- des
    # for every row in the design
    for (r in 1:nrow(des)) {
      # switch with everey row in candidate set 
      for (c in 1:nrow(cand.set)) {
        des[r, ] <- cand.set[c, ]
        # calculate d errors
        d.errors <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
        #db error 
        db <- mean(d.errors, na.rm = TRUE)
        #change if lower db error
        if (db < db.start) {
          best.row <- as.numeric(des[r, ])
          db.start <- db
          change <- TRUE
        }
        #show progress iteration
        setTxtProgressBar(pb, r)
      }
      #replace with best profile if change
      if (change) {
        des[r, ] <- best.row
      } else {
        des[r, ] <- iter.des[r, ]
      }
      #initialize variables again 
      change <- FALSE
      na.percentage <- 0
    }
    close(pb)  # print progress 
    converge <- isTRUE(all.equal(des, iter.des)) # convergence if no profile is swapped this iteration 
  }
  #calculate percentage NA values
  d.errors <- apply(par.samples, 1, Derr, des = des,  n.alts = n.alts)
  if (any(is.na(d.errors))) {
    na.percentage <- sum(is.na(d.errors))/n.samples
  } 
  #return design, db error, na.values
  return(list(des, db, na.percentage))
}


#' DB sequential set
#'
#' Adds the most DB efficient choice set to a design, given (updated) parameter values.
#' @param des A design matrix in which each row is a profile.
#' @param cand.set A numeric matrix in which each row is a possible profile.
#' @param n.alts Numeric value indicating the number of alternatives per choice set.
#' @param par.samples A matrix in which each row is a sample from the multivariate parameter distribution.
#' @param weights A vector containing the weights of the samples.
#' @param prior.covar Covariance matrix of the prior distribution.
#' @return The most db efficient choice set.
#' @export
seqfed.db <-function(des, cand.set, n.alts, par.samples, weights, prior.covar){

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
#' Provides the set that maximizes the Kullback-Leibler divergence, given parameter values.
#' @param lvls A vector which contains for each attribute, the number of levels.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice set.
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











