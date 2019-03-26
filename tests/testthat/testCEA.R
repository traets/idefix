context("Tests for CEA")
test_that("Errors for creating initial design in CEA function", {
  ## Categorical attributes
  # Misspecification of coding type
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "B")), 
               "coding argument is incorrect.")
  # Less coding types than attributes
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D")), 
               "coding argument is incorrect.")
  # One attribute only
  expect_error(CEA(lvls = c(3), coding = c("D")), "lvls argument is incorrect.")
  # Non-numeric number of levels of an attribute
  expect_error(CEA(lvls = c(3, 3, "d"), coding = c("D", "D", "D")), 
               "lvls argument is incorrect.")
  # Incorrect type of coding
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "f")), 
               "coding argument is incorrect.")
  
  ## Continuous attributes
  # Continuous levels missing
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("C", "C", "C")), 
               "when 'coding' contains C, 'c.lvls' should be specified")
  # Less number of continuous levels than attributes
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("C", "C", "C"), 
      c.lvls = list(c(4, 6, 8), c(2, 4,6))), 
      "length of 'c.lvls' does not match number of specified continuous attributes in 'coding")
  # Misspecification in "c.lvls" according to the number of levels in "lvls"
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("C", "C", "C"), 
                   c.lvls = list(c(4, 6), c(2, 4, 6), c(5, 6, 7))), 
               "the number of levels provided in 'c.lvls' does not match the expected based on 'lvls'")
  
  ## Design specifications
  # More alternative constants than alternatives in a choice set
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 2, par.draws = c(0, 0, 0, 0, 0, 0), 
                   alt.cte = c(1, 0, 0)), 
                   "'n.alts' does not match the 'alt.cte' vector")
  # An alternative constant equals to 2
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 2, par.draws = c(0, 0, 0, 0, 0, 0), 
                   alt.cte = c(2, 0)), 
               "'alt.cte' should only contain zero or ones.")
  # No boolean value in nochoice parameter
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 2, alt.cte = c(1, 0), 
                   par.draws = c(0, 0, 0, 0, 0, 0), no.choice = T), 
               "if 'no.choice' is TRUE, the last alternative constant should equal 1.")
  # 1 alternative constant and not a list in the draws
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 2, alt.cte = c(1, 0), 
                   par.draws = c(0, 0, 0, 0, 0, 0), no.choice = F), 
  "par.draws should be a list")
  # 1 alternative constant and draws in a list, but a single component
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 2, alt.cte = c(1, 0), 
                   par.draws = list(c(0, 0, 0, 0, 0, 0, 0)), no.choice = F), 
  "'par.draws' should contain two components")
  # 1 alternative constant and draws for the attributes are not in a matrix
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 2, alt.cte = c(1, 0), 
                   par.draws = list(0, c(0, 0, 0, 0, 0, 0)), no.choice = F), 
  "'par.draws' should contain two matrices")
  # Different number of draws for the alternative constant and betas
  # Note: There should be the same number of draws for both components
  # All these values are random
  mu <- c(0.5, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) 
  v <- diag(length(mu)) # Prior variance.
  set.seed(123) 
  pd <- MASS::mvrnorm(n = 2, mu = mu, Sigma = v) # 10 draws.
  p.d <- list(matrix(pd[,1], ncol = 2), pd[,2:7])
  p.d[[1]] <- p.d[[1]][,-2] # Remove a draw for the constant
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 2, alt.cte = c(1, 0), 
                   par.draws = p.d, no.choice = F), 
  "the number of rows in the components of 'par.draws' should be equal")
  # 2 alternative constants and not a list in the draws
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 3, alt.cte = c(1, 1, 0), 
                   par.draws = c(0, 0, 0, 0, 0, 0, 0, 0), no.choice = F), 
               "par.draws should be a list")
  # 2 alternative constants and draws in a list, but a single component
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 3, alt.cte = c(1, 1, 0), 
                   par.draws = list(c(0, 0, 0, 0, 0, 0, 0, 0)), no.choice = F),
               "'par.draws' should contain two components")
  # 2 alternative constants and draws for the attributes are not in a matrix
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 3, alt.cte = c(1, 1, 0), 
                   par.draws = list(c(0, 0), c(0, 0, 0, 0, 0, 0)), 
                   no.choice = F), 
               "'par.draws' should contain two matrices")
  # 2 alternative constants and draws for only one of them
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 3, alt.cte = c(1, 1, 0), 
                   par.draws = list(as.matrix(c(0)), 
                                    as.matrix(c(0, 0, 0, 0, 0, 0))), 
                   no.choice = F), 
"the first component of 'par.draws' should contain the same number
           of columns as there are non zero elements in 'alt.cte'")
  # Different number of draws for the alternative constants and betas
  # Note: There should be the same number of draws for both components
  # All these values are random
  mu <- c(0.5, 0.5, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) 
  v <- diag(length(mu)) # Prior variance.
  set.seed(123) 
  pd <- MASS::mvrnorm(n = 3, mu = mu, Sigma = v) # 10 draws.
  p.d <- list(matrix(pd[,1:2], ncol = 2), pd[,3:8])
  p.d[[1]] <- p.d[[1]][-2,] # Remove a draw for the constant
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 8,
                   n.alts = 3, alt.cte = c(1, 1, 0), 
                   par.draws = p.d, no.choice = F), 
               "the number of rows in the components of 'par.draws' should be equal")
  # Number of choice sets is smaller than parameters to estimate
  mu <- c(0.5, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) 
  v <- diag(length(mu)) # Prior variance.
  set.seed(123) 
  pd <- MASS::mvrnorm(n = 3, mu = mu, Sigma = v) # 10 draws.
  p.d <- list(matrix(pd[,1], ncol = 1), pd[,2:7])
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.sets = 6,
                   n.alts = 2, alt.cte = c(1, 0), 
                   par.draws = p.d, no.choice = F),  
  "Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
})


