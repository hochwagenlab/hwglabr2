context("Elapsed time")

test_that('elapsed_time works for sec., min., and hrs.', {
  expect_equal(elapsed_time(50, 80), "30 sec.")
  expect_equal(elapsed_time(1, 150), "2.5 min.")
  expect_equal(elapsed_time(1, 8800), "2.44 hrs.")
  expect_equal(elapsed_time(1, 10850), "3.01 hrs.")
})

test_that('elapsed_time throws error with non-numeric input', {
  expect_error(elapsed_time(1, "a"), '"t0" and "t1" must*')
  expect_error(elapsed_time("a", 10), '"t0" and "t1" must*')
  expect_error(elapsed_time("a", "b"), '"t0" and "t1" must*')
})
