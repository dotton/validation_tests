#if(!require('myPackage')) install.packages('myPackage')

#library(myPackage)
library(testthat)

set.seed(1234)
dat <- runif(1e5)

withr::defer({
#  detach(package:myPackage)
}, teardown_env())
