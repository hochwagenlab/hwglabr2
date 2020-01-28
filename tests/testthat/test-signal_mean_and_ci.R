context("Signal mean an CI")

set.seed(42)
input <- matrix(rexp(10), ncol=2)
#input_df <- as.data.frame(input)

# Dimnames can cause the matrices to be different
strip_dimnames <- function(x) {
  dimnames(x) <- NULL
  x
}

output <- cbind(Mean=c(0.3308183, 0.8188403),
                Lower=c(0.2133957, 0.5452756),
                Upper=c(0.4345073, 1.1147480),
                SD=c(0.2421973, 0.4968484),
                SEM=c(0.1083139, 0.2221974))

output <- strip_dimnames(output)

test_that('signal_mean_and_ci produces expected output matrix', {
  set.seed(42)
  testthat::expect_equal(strip_dimnames(signal_mean_and_ci(input,
                                                           rep_bootstrap=10)),
                         output, tolerance=1e-5)
})

#test_that('signal_mean_and_ci can handle an input data frame', {
#  set.seed(42)
#  testthat::expect_equal(strip_dimnames(signal_mean_and_ci(input_df, rep_bootstrap=10)),
#               output, tolerance=1e-5)
#})

test_that('signal_mean_and_ci throws error with object non-coercible to df', {
  testthat::expect_error(signal_mean_and_ci("a"), '"signal_data" must*')
})

