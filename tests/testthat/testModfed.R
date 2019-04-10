context("Tests for Modfed")
test_that("Errors in Modfed function", {
  cand.set <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D"))
  mu <- c(0.8, 0.2, -0.3, -1.2, 1.6, 2.2) # Prior parameter vector
  v <- diag(length(mu)) # Prior variance.
  set.seed(123) 
  pd <- MASS::mvrnorm(n = 1, mu = mu, Sigma = v) # 10 draws.
  #p.d <- list(NA, pd[,2:6])
  expect_error(Modfed(cand.set = cand.set, n.sets = 8, n.alts = 2, 
         alt.cte = c(1,0), parallel = FALSE, par.draws = pd, best = T),
         "par.draws should be a list")
  
  # 1 alternative constant in a choice set with 2 alternatives
  mu <- c(0.5, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) # Prior parameter vector
  v <- diag(length(mu)) # Prior variance.
  set.seed(123) 
  pd <- MASS::mvrnorm(n = 2, mu = mu, Sigma = v) # 10 draws.
  p.d <- list(matrix(pd[,1], ncol = 1), pd[,2:7])
  expect_error(Modfed(cand.set = cand.set, n.sets = 8, n.alts = 2, 
         alt.cte = NULL, parallel = FALSE, par.draws = p.d, best = T))
})