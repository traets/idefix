

#' DB modified federov algorithm
#'
#' Modified federov algorithm. Swipes every profile of an initial design with candidate profiles.
#' By doing this it tries to minimize the DB error. The DB error is te mean D-error over samples
#' from the prior parameter distribution.
#' @param candset A numeric matrix in which each row is a possible profile.
#' @param n_sets Numeric value indicating the number of choice sets.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param par_samples A matrix in which each row is a sample from the multivariate prior parameter distribution.
#' @param max_iter A numeric value indicating the maximum number allowed iterations.
#' @return An efficient design, and the associated DB-error.
#' @export
modfed.db<-function (candset, n_sets, n_alts, par_samples, max_iter = Inf){

  #random start des
  R<-round(runif((n_sets*n_alts), 1, nrow(candset)))
  des<-data.matrix(candset[R,])

  if (ncol(des) != ncol(par_samples)) {
    stop("Number of paramters does not match the design matrix")
  }
  if ((nrow(des)/n_alts)%%1 != 0) {
    stop("number of alternatives per choice set does not match the design")
  }

  DB_start <- 99999
  converge <- FALSE
  it <- 1

  while (!converge & it <= max_iter) {
    it <- it + 1
    pb <- txtProgressBar(min = 0, max = nrow(des), style = 3)
    iter_des <- des

    for (r in 1:nrow(des)) {
      for (c in 1:nrow(candset)) {
        des[r, ] <- candset[c, ]
        D_errors <- apply(par_samples, 1, d_err, des = des,  n_alts = n_alts)

        if (any(is.na(D_errors))) {
          warning("NaN produced in DB error. Possibly too little choice sets.")
        }
        if (any(is.infinite(D_errors))) {
          warning("Infinite D-error values. Possibly too little choice sets.")
        }
        DB <- mean(D_errors, na.rm = TRUE)
        if (DB < DB_start) {
          best_row <- as.numeric(des[r, ])
          DB_start <- DB
          change <- T
        }
        setTxtProgressBar(pb, r)
      }

      if (change) {des[r, ] <- best_row
      }else {des[r, ] <- iter_des[r, ]}

      change <- F
    }
    close(pb)
    converge <- isTRUE(all.equal(des, iter_des))
  }
  return(list(des, DB))
}


#' DB sequential set
#'
#' Adds the most DB efficient choice set to a design, given (updated) parameter values.
#' @param des A design matrix in which each row is a profile.
#' @param candset A numeric matrix in which each row is a possible profile.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param par_samples A matrix in which each row is a sample from the multivariate parameter distribution.
#' @param weights A vector containing the weights of the samples.
#' @param prior_covar Covariance matrix of the prior distribution.
#' @return The most DB efficient choice set.
#' @export
seqfed.db <-function(des, candset, n_alts, par_samples, weights, prior_covar){

  #start values
  DB_best<-1000
  full_comb<- full_sets(candset, n_alts)

  #for each potential set:
  for (p in 1:nrow(candset)){

    set<-as.matrix(candset[as.numeric(full_comb[p, ]), ])

    #for each par draw calculate d-error:
    DB_error<-0

    for (s in 1: nrow(par_samples)){

      info_s<-info_design(par = par_samples[s, ], des = set, n_alts = n_alts)
      info_d<-info_design(par= par_samples[s, ], des = des, n_alts = n_alts)
      inv_cov_prior<-solve(prior_covar)

      DB_error<- DB_error+(det(info_d + info_s + inv_cov_prior)*(-1/ncol(par_samples))*weights[s])

    }

    #select new set if smaller DB error
    if (DB_error< DB_best){final_set<-set; DB_best<- DB_error}
  }

  return(final_set)

}


#' KL information
#'
#' Calculates the Kullback-Leibler divergence for a choice set, given parameter values.
#' @param set Numeric matrix in which each row is a profile.
#' @param par_samples A matrix in which each row is a sample.
#' @param weights A vector containing the weights of the samples.
#' @return The Kullback-Leibler divergence
KL <- function (set, par_samples, weights){

  #probability
  num<-set%*%t(par_samples)
  mmat<-as.matrix(t(apply(num[, 1:ncol(num)], 2, max)))
  nummax<-exp(sweep(num, MARGIN=2, mmat, FUN="-"))
  denom<-colSums(nummax)

  probs<-sweep(nummax, MARGIN=2, denom, FUN="/")
  logprobs<-log(probs)

  wprob<- sweep(probs, MARGIN=2, weights, FUN="*")
  totwprob<- rowSums(wprob)

  logwprob<- sweep(logprobs, MARGIN=2, weights, FUN="*")
  totlogwprob<- rowSums(logwprob)

  #kullback Leibler information
  klinfo<- 0
  for( a in 1:n_alts){
    klinfo<- klinfo + (totwprob[a] * (log(totwprob[a]) - totlogwprob[a]))
  }

  return (as.numeric(klinfo))
}


#' KL set selecting
#'
#' Provides the set that maximizes the Kullback-Leibler divergence, given parameter values.
#' @param lvls A vector which contains for each attribute, the number of levels.
#' @param n_sets Numeric value indicating the number of choice sets.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param par_samples A matrix in which each row is a sample.
#' @param weights A vector containing the weights of the samples.
#' @return The most efficient choice set
#' @export
KL_select <- function(lvls, n_sets, n_alts, par_samples, weights){

  #All choice sets, without same profile twice
  fp<-profiles(lvls)
  fcomb<-full_sets(cand = fp, n_alts = n_alts)
  rows<-apply(fcomb[, 1:ncol(fcomb)], 1, function(i) length(unique(i)) > ncol(fcomb)-1)
  fcomb<-fcomb[rows, ]

  #start value
  kl_start <- 0
  best_set<-0

  for ( s in 1:nrow(fcomb)){

    #take set
    set<-as.matrix(fp[as.numeric(fcomb[s, ]), ])

    #calculate for all sets the KLinfo.
    klinfo<-KL(set, par_samples, weights)

    #if better --> keep
    if (klinfo > kl_start){
      best_set<- set
      kl_start<-klinfo
    }

    print(klinfo)
  }


  return(list(best_set, kl_start))
}













