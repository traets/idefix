

#' log Posterior
#'
#' Calculates the logposterior with a normal prior density
#' @param par Numeric vector with parametervalues.
#' @param des A design matrix in which each row is a profile.
#' @param Y A binary response vector.
#' @param n_alts The number of alternatives in each choice set.
#' @param prior_mean vector containing the prior mean.
#' @param prior_covar matrix containing the prior covariance.
#' @return the logposterior probability
#' @export
logPost<-function(par, prior_mean, prior_covar, des,  n_alts, Y ){

  p<-t(t(des) * par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)

  log_L<-sum(Y*log(p))

  logprior2=-0.5*(par -t(prior_mean))%*%solve(prior_covar)%*%(as.matrix(par) - prior_mean)

  logpost<-(length(prior_mean)/2*log(2*pi)-0.5*log(det(prior_covar))) + logprior2 + log_L

  return(logpost)

}

#' Hessian
#'
#' @param par Numeric vector with parametervalues.
#' @param des A design matrix in which each row is a profile.
#' @param covar The covariance matrix.
#' @param n_alts The number of alternatives in each choice set.
#' @return the hessian matrix
hessian<-function (par, des, covar, n_alts){

  des<-as.matrix(des)
  p<- des %*% diag(par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)

  info<-crossprod(des*p, des) - crossprod(rowsum(des*p, rep(seq(1, nrow(des)/n_alts, 1), each=n_alts)))

  hess<-(-info - solve(covar))

  return(hess)

}

#' Likelihood function
#'
#' @param par Numeric vector with parametervalues.
#' @param des A design matrix in which each row is a profile.
#' @param n_alts The number of alternatives in each choice set
#' @param Y A binary response vector.
#' @return the likelihood
Lik<-function(par, des, n_alts, Y){

  des<-as.matrix(des)

  p<-t(t(des) * par)
  p<-.rowSums(p, m= nrow(des), n=length(par))
  expp<-exp(p)
  p<-expp/rep(rowsum(expp, rep(seq(1, nrow(des)/n_alts, 1), each = n_alts)), each=n_alts)

  L<-prod(p^Y)

  return(L)

}


#' Density multivariate t-distribution
#'
#' @param par Numeric vector with parametervalues.
#' @param g_mean vector containing the mean of the multivariate t-distribution.
#' @param g_covar covariance matrix of the multivariate t-distribution.
#' @return density
g_dens<-function (par, g_mean, g_covar){

  df=length(g_mean)
  n<-length(par)
  dif<-g_mean-par
  invcov<-solve(g_covar)
  differ<-as.numeric(t(dif)%*% invcov %*% dif)

  iMVSTd=1/(det(g_covar)^(0.5))*(1+((1/df)*differ))^(-(df+length(par))/2)
}


#' Importance sampling
#'
#' This functions samples from an imortance density (multivariate t-distribution),
#' and gives weightes to the samples according to the posterior distribution. The prior is
#' a normal distribution.
#' @param prior_mean Numeric vector which is the mean of the multivariate normal distribution (prior).
#' @param prior_covar A matrix, the covariance matrix of the multivariate normal distribution (prior).
#' @param des A design matrix in which each row is a profile.
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param Y A binary response vector.
#' @param m Numeric value. Number of samples = base^m.
#' @param b Numeric value indicating the base (default = 2).
#' @return A list containing samples, their associated weights, the maximum likelihood estimates and the estimated covariance matrix.
#' @export
imp_sampling <- function (prior_mean, prior_covar, des,  n_alts, Y, m, b=2, ...){

  #cte
  prior1<-(2*pi)^(-length(prior_mean)/2)*(det(prior_covar))^(-0.5)

  #estimate importance mean
  maxest<-maxLik::maxNR(logPost, start= prior_mean, prior_mean = prior_mean, prior_covar= prior_covar,
                        des = des, Y=Y, n_alts=n_alts)$estimate

  #draws from importance density
  H<-hessian(par = maxest, des = des, covar = prior_covar, n_alts = n_alts)
  g_covar<--solve(H)
  g_draws<-lattice_mvt(mean=maxest, cvar = g_covar, df=length(maxest), m=m, ...)

  #vectors
  prior=LK=dens_g=weights<- numeric(nrow(g_draws))

  for (r in 1:nrow(g_draws)){

    #prior
    prior[r]<-prior1*exp(-0.5* (g_draws[r, ]-prior_mean) %*% solve(prior_covar) %*% as.matrix(g_draws[r, ]- prior_mean))
    #likelihood
    LK[r]<-Lik(par = g_draws[r, ], des = des, Y=Y, n_alts = n_alts)
    #density of g
    dens_g[r]<-g_dens(par = g_draws[r, ], g_mean = maxest, g_covar = g_covar)

  }

  #compute the weights of samples
  w<-LK*prior/dens_g
  w<-w/sum(w)

  return(list(g_draws, w, maxest, g_covar))
}





