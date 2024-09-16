# Write relevant tests for the function in here
# Consider the type of function:
#   - is it deterministic or statistic?
#   - is it worth checking for errors/warnings under particular conditions?
local_edition(3)

## using same example of an lm object as in test-lm.R
data(mtcars)
mtcars$cyl_f <- factor(mtcars$cyl)

cmod <- lm(mpg ~ cyl, data = mtcars)
fmod <- lm(mpg ~ cyl_f, data = mtcars)

## coefficients and SEs have already been checked in test-lm.R
test_that("confint returns the expected confidence bounds", {
  expect_equal(unname(confint(cmod)[2,]), 
               coef(cmod)[2] + qt(0.975, df = 30) * c(-1,1) * summary(cmod)$coefficients[2, 2])
  expect_equal(unname(confint(fmod)[2,]), 
               coef(fmod)[2] + qt(0.975, df = 29) * c(-1,1) * summary(fmod)$coefficients[2, 2])
  expect_equal(unname(confint(fmod)[3,]), 
               coef(fmod)[3] + qt(0.975, df = 29) * c(-1,1) * summary(fmod)$coefficients[3, 2])
})
