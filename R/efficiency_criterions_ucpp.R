# DerrS_P using InfoDes_cpp and det_cpp
DerrS.P_ucpp <- function(par, des, n.alts, i.cov) {
  info <- InfoDes_cpp(par = par, des = des, n_alts = n.alts)
  d.error <- det_cpp(info + i.cov)^(-1 / length(par))
  return(d.error)
}

# DerrC using cpp functions
DerrC_ucpp <- function(par, des, n.alts, i.cov) {
  info.des <- InfoDes_cpp(par, des, n.alts)
  detinfo <- det_cpp(info.des + i.cov)
  ifelse((detinfo <= 0), return(NA), return(detinfo^(-1 / length(par))))
}

# DBerrS.P using DerrS.P_cpp
DBerrS.P_ucpp <- function(des, par.draws, n.alts, i.cov, weights) {
  # Add alternative specific constants if necessary
  # For each draw calculate D-error.
  d.errors <- apply(par.draws, 1, DerrS.P_ucpp, des, n.alts, i.cov)
  w.d.errors <- d.errors * weights
  # DB-error. 
  db.error <- mean(w.d.errors, na.rm = TRUE)
  return(db.error)
}