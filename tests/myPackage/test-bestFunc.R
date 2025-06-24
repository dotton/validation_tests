# Write relevent tests for the function in here
# Consider the type of function:
#   - is it deterministic or statistic?
#   - is it worth checking for errors/warnings under particular conditions?
local_edition(3)

test_that("Test something", {
  m <- sum(dat)/length(dat)
  
  expect_equal(m, 0.5, tolerance = 1e-3) 
})
