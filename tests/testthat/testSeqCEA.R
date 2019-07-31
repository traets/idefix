context("Tests for SeqCEA")
test_that("Errors for creating initial design in CEA function", {
  ## Initial design
  m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4) # mean (total = 6 parameters).
  pc <- diag(length(m)) # covariance matrix
  set.seed(123)
  sample <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
  # Initial design.
  des <- list(example_design)
  # Initial design should be a matrix
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = sample, prior.covar = pc, 
                      parallel = FALSE), "'des' should be a matrix or NULL")
  
  # The number of rows in the initial design has to be a multiple of the number
  # of alternatives per choice set
  des <- example_design
  # Initial design should be a matrix
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 3, par.draws = sample, prior.covar = pc, 
                      parallel = FALSE), "'n.alts' does not seem to match with the number of rows in 'des'")
  
  # The number of alternative constants should be the same as alternatives per 
  # choice set. Here one constant is removed
  m2 <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4, 1.8, 1.2) # mean 
  pc2 <- diag(length(m2)) # covariance matrix
  pos <- MASS::mvrnorm(n = 10, mu = m2, Sigma = pc2)
  sample2 <- list(pos[ , 1:2], pos[ , 3:8])
  des2 <- example_design2 
  expect_error(SeqCEA(des = des2, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 3, par.draws = sample, alt.cte = c(1, 1), 
                      prior.covar = pc, parallel = FALSE), "'n.alts' does not match the 'alt.cte' vector")
  
  # Alternative constants can only be zeros or ones. Here is set to 2.  
  expect_error(SeqCEA(des = des2, lvls = c(3, 3, 3), coding = c("E", "E", "E"),
                      n.alts = 3, par.draws = sample, alt.cte = c(1, 2, 0),
                      prior.covar = pc, parallel = FALSE), "'alt.cte' should only contain zero or ones.")
  
  # No choice should be an integer or NULL not a decimal (here)
  expect_error(SeqCEA(des = des2, lvls = c(3, 3, 3), coding = c("E", "E", "E"),
                      n.alts = 3, par.draws = sample, alt.cte = c(1, 1, 0),
                      no.choice = 1.5, prior.covar = pc, parallel = FALSE), 
               "'no.choice' should be an integer or NULL")
  
  # No choice should be the position of the alt.cte object that indicates the 
  # constant for no.choice option. Here this position is 4 and there are only 3
  expect_error(SeqCEA(des = des2, lvls = c(3, 3, 3), coding = c("E", "E", "E"),
                      n.alts = 3, par.draws = sample, alt.cte = c(1, 1, 0),
                      no.choice = 4, prior.covar = pc, parallel = FALSE), 
               "'no.choice' does not indicate one of the alternatives")
  
  # No choice should be the position of the alt.cte object. Here there is no
  # alt.cte
  expect_error(SeqCEA(des = des2, lvls = c(3, 3, 3), coding = c("E", "E", "E"),
                      n.alts = 3, par.draws = sample, alt.cte = NULL,
                      no.choice = 1, prior.covar = pc, parallel = FALSE), 
               "if there is a no choice alternative, 'alt.cte' should be specified")
  
  # No choice position of the alt.cte object should be one. Here there is zero
  expect_error(SeqCEA(des = des2, lvls = c(3, 3, 3), coding = c("E", "E", "E"),
                      n.alts = 3, par.draws = sample, alt.cte = c(1, 1, 0),
                      no.choice = 3, prior.covar = pc, parallel = FALSE), 
               "the no choice alternative should correspond with a 1 in 'alt.cte'")
  
  # When there is only one alternative constant, par draws should be a list
  m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4, 0.4) # mean (total = 6 parameters).
  pc <- diag(length(m)) # covariance matrix
  set.seed(123)
  sample <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
  # Initial design.
  des <- example_design2[,-1]
  # Initial design should be a matrix
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = sample, prior.covar = pc, 
                      alt.cte = c(1, 0),
                      parallel = FALSE), "'par.draws' should be a list when 'alt.cte' is not NULL")
  
  # When there is only one alternative constant, par draws should be a list with
  # two components, here there is only one
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = list(sample), prior.covar = pc, 
                      alt.cte = c(1, 0), parallel = FALSE), 
               "'par.draws' should contain two components")
  
  # When there is only one alternative constant, each element of par.draws should
  # be a matrix. Here the second component is a list
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = list(sample[ , 1], 
                                                   list(sample[ , 2:7])),
                      prior.covar = pc, 
                      alt.cte = c(1, 0), parallel = FALSE), 
               "'par.draws' should contain two matrices")
  
  # When there is only one alternative constant, there must only one column in 
  # the first element of par.draws. Here there are two elements
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = sample2,
                      prior.covar = pc, 
                      alt.cte = c(1, 0), parallel = FALSE), 
               "the first component of 'par.draws' should contain the same number 
               of columns as there are non zero elements in 'alt.cte'")
  
  # When there is only one alternative constant, the number of samples from the
  # distribution of beta should be equal in both components of par.draws. Here
  # the first component has 1 row less.
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = list(sample[-1 , 1], 
                                                   sample[, 2:7]),
                      prior.covar = pc, alt.cte = c(1, 0), parallel = FALSE), 
               "the number of rows in the components of 'par.draws' should be equal")
  
  # The lenght of the weigths should be the same as the number of samples from
  # the distribution of beta
  w <- runif(9)
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = list(sample[ , 1], 
                                                   sample[, 2:7]),
                      prior.covar = pc, alt.cte = c(1, 0), parallel = FALSE, 
                      weights = w), 
               "length of 'weights' does not match number total number of rows in 'par.draws'")
  
  # The number of columns in the prior covariance matrix should be equal to the 
  # number of columns in the design matrix
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = list(sample[ , 1], 
                                                   sample[, 2:6]),
                      prior.covar = pc[,-1], alt.cte = c(1, 0), parallel = FALSE), 
               "number of columns of 'prior.covar' does not equal the number of columns of design matrix \\(including alternative specific constants\\)")
  
  # The number of columns in par.draws should be the same as the number of columns
  # of the initial design
  expect_error(SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("E", "E", "E"), 
                      n.alts = 2, par.draws = list(sample[ , 1], 
                                                   sample[, 2:6]),
                      prior.covar = pc, alt.cte = c(1, 0), parallel = FALSE), 
               "number of columns in 'par.draws' does not match the number of columns in 'des'")
  
  # The number of columns in the initial design should be the same as the number
  # of columns of the design matrix
  # It seems that I dont need this test, because there is already a conditional
  # for the number of ncol(prior.covar) with the number of colums in the design
  # matrix, and also another for the number of columns in par.draws and the
  # number of columns in des
  # expect_error(SeqCEA(des = example_design, lvls = c(2, 3, 3),
  #                      coding = c("E", "E", "E"), n.alts = 2,
  #                      par.draws = sample,
  #                      prior.covar = pc, alt.cte = NULL, parallel = FALSE),
  #               "number of columns in 'des' does not match the number of columns of design matrix \\(including alternative specific constants\\)")
  
})



