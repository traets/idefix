# Fullsets using Infodes_cpp
Fullsets_ucpp <- function(cand.set, n.alts, no.choice, reduce = TRUE){
  
  if(!is.null(no.choice)){
    n.alts <- n.alts - 1
  }
  full.comb <- utils::combn(1:nrow(cand.set), n.alts, FUN = function(x)  cand.set[x, ], simplify = FALSE)
  #reduce
  if (reduce){
    m <- stats::rnorm(ncol(cand.set))
    inf <-list()
    for(i in 1:length(full.comb)){
      inf[[i]] <- round(InfoDes_cpp(m, full.comb[[i]], n.alts), digits = 3)
    }
    t <- array(unlist(inf), dim = c(length(m), length(m), length(inf))) 
    full.comb <- full.comb[!duplicated(t, MARGIN = 3)]
  }
  if(!is.null(no.choice)){
    full.comb <- lapply(full.comb, Inchoice, no.choice = no.choice)
  }
  return(full.comb)
}

