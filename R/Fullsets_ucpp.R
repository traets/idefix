# Fullsets using Infodes_cpp
Fullsets_ucpp <- function(cand.set, n.alts, no.choice, reduce = TRUE, allow.rep,
                           des, n.cte){
  if (!is.null(no.choice)) {
    n.alts <- n.alts - 1
  }
  full.comb <- utils::combn(1:nrow(cand.set), n.alts, 
                            FUN = function(x)  cand.set[x, ], simplify = FALSE)
  
  # If no repeated choice sets
  if (!allow.rep) { 
    n.sets <- nrow(des) / n.alts # Number of choice sets in initial design
    i <- 1 # loop counter for rows in initial design
    repe <- matrix(F, length(full.comb), n.sets) # True where the choice set is
                                            # located
    j <- 1 # counter of choice sets inside the loop
    while (i < nrow(des)) {
      cset <- des[i:(i + n.alts - 1),] # Subset a choice set from initial design
      #per <- expand.grid(1:nrow(cset),1:nrow(cset),1:nrow(cset)) # Compute permutations
      per <- do.call(expand.grid, rep(list(1:nrow(cset)), n.alts))
      per <- per[apply(per, 1, function(x) {length(unique(x)) == n.alts}),]
      cset_comb <- list()
      for (k in 1:nrow(per)) {
        cset_comb[[k]] <- cset[as.numeric(per[k,]),] # Permutations
      } # End for
      # in repe: True where the choice set is located
      if (n.cte == 0) {
        repe[, j] <- unlist(lapply(full.comb, function(x){ 
        any(unlist(lapply(cset_comb, function(y) all(x == y))))})) 
      } else {
        repe[, j] <- unlist(lapply(full.comb, function(x){ 
          any(unlist(lapply(cset_comb, function(y) {
            y <- y[, -c(1:n.cte)] # Remove alternative constants from des
            all(x == y)
            } # End function(y)
            )))
          } # End function(x)
          )) 
      }
      i <- i + n.alts
      j <- j + 1
    } # End while
    # delete: index of repeated choice sets
    delete <- apply(repe, 2, function(x) which(x == T)) 
    full.comb <- full.comb[-delete]
  }
  
  #reduce: Remove choice sets with same information matrix
  if (reduce) {
    m <- stats::rnorm(ncol(cand.set))
    inf <- list()
    for (i in 1:length(full.comb)) {
      inf[[i]] <- round(InfoDes_cpp(m, full.comb[[i]], n.alts), digits = 3)
    }
    t <- array(unlist(inf), dim = c(length(m), length(m), length(inf))) 
    full.comb <- full.comb[!duplicated(t, MARGIN = 3)]
  }
  if (!is.null(no.choice)) {
    full.comb <- lapply(full.comb, Inchoice, no.choice = no.choice)
  }
  return(full.comb)
}


# This function is the equivalent of Fullsets for CEA algorithm. Instead
# of computing all possible combinations of alternatives, this new function
# creates alternatives by randomly choosing the level of each attribute
# n.cs: is the number of random choice sets to 
# c.names: is the column names of the design matrix
Newsets_ucpp <- function(levels.list, n.alts, no.choice, reduce = TRUE, 
                         n.cs = NULL, c.names = NULL){
  if (!is.null(no.choice)) {
    n.alts <- n.alts - 1
  }
  if (is.null(n.cs)) {
    n.cs <- prod(unlist(lapply(levels.list, nrow)))
  }
  # Generation of random choice sets
  full.comb <- vector(mode = 'list', length = n.cs)
  for (i in 1:n.cs) {
    # r is to know which levels to take in each attribute
    r <- NULL
    cs <- NULL
    for (j in 1:length(levels.list)) {
      r <- round(stats::runif(n.alts, 1, nrow(levels.list[[j]])))
      cs <- cbind(cs, levels.list[[j]][r,])
    }
    colnames(cs) <- c.names
    rownames(cs) <- NULL
    full.comb[[i]] <- cs
  } # end loop i
  #reduce
  if (reduce) {
    m <- stats::rnorm(ncol(full.comb[[1]]))
    inf <- list()
    for (i in 1:length(full.comb)) {
      inf[[i]] <- round(InfoDes_cpp(m, full.comb[[i]], n.alts), digits = 3)
    }
    t <- array(unlist(inf), dim = c(length(m), length(m), length(inf))) 
    full.comb <- full.comb[!duplicated(t, MARGIN = 3)]
  }
  if (!is.null(no.choice)) {
    full.comb <- lapply(full.comb, Inchoice, no.choice = no.choice)
  }
  return(full.comb)
}





