test_that("Errors for creating initial design in CEA function", {
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "B")), 
               "coding argument is incorrect.")
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D")), 
               "coding argument is incorrect.")
  expect_error(CEA(lvls = c(3), coding = c("D")), "lvls argument is incorrect.")
  expect_error(CEA(lvls = c(3, 3, "d"), coding = c("D", "D", "D")), 
               "lvls argument is incorrect.")
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("D", "D", "f")), 
               "coding argument is incorrect.")
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("C", "C", "C")), 
               "when 'coding' contains C, 'c.lvls' should be specified")
  expect_error(CEA(lvls = c(3, 3, 3), coding = c("C", "C", "C"), 
      c.lvls = list(c(4, 6, 8), c(2, 4,6))), 
      "length of 'c.lvls' does not match number of specified continuous attributes in 'coding")
      expect_error(CEA(lvls = c(3, 3, 3), coding = c("C", "C", "C"), 
      c.lvls = list(c(4, 6), c(2, 4, 6), c(5, 6, 7))), 
      "the number of levels provided in 'c.lvls' does not match the expected based on 'lvls'")
})