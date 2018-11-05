####################softmax gradient for each ec###############
emp <- matrix(ncol = 1, nrow = 0)

##~~~~~~~~~~~~~~~~~~~~~~~~~~SingleSpeGradSM~~~~~~~~~~~~~~~~~~~~~~~~
test_that('"1, 1" | 1', {
  expect_equal(SingleSpeGradSM(list(c(1, 1), 1), list(c(1, 1), 1), list(c(1, 1), 1), c(1, 1), 0), matrix(c(0.5, 0.5), ncol = 1))
})

test_that('1, 1 | "1"', {
  expect_equal(SingleSpeGradSM(list(c(1, 1), 1), list(c(1, 1), 1), list(c(1, 1), 1), c(1, 1), 1), matrix(c(1), ncol = 1))
})

test_that('"0, 1" | 1', {
  expect_equal(SingleSpeGradSM(list(c(1, 1), 1), list(1, 1), list(1, 1), c(0.5, 1), 0), matrix(c(2/3), ncol = 1))
})

test_that('0, 1 | "1"', {
  expect_equal(SingleSpeGradSM(list(c(1, 1), 1), list(1, 1), list(1, 1), c(0.5, 1), 1), matrix(c(1), ncol = 1))
})

test_that('"1, 0" | 0', {
  expect_equal(SingleSpeGradSM(list(c(1, 1), 1), list(1, emp), list(1, emp), c(0.5, 0), 0), matrix(c(1), ncol = 1))
})

test_that('"1, 1" | 0', {
  expect_equal(SingleSpeGradSM(list(c(1, 1), 1), list(c(1, 1), emp), list(c(1, 1), emp), c(1, 0), 0), matrix(c(0.5, 0.5), ncol = 1))
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~ECGradSM func~~~~~~~~~~~~~~~~~~~~~~~~~
test_that('1, 1 | 1', {
  expect_equal(ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(c(1, 1), 1), list(c(1, 1), 1)), matrix(c(0.5, 0.5, 1), ncol = 1))
})

test_that('0, 1 | 1', {
  expect_equal(ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(1, 1), list(1, 1)), matrix(c(2/3, 1), ncol = 1))
})

test_that('1, 0 | 1', {
  expect_equal(ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(1, 1), list(1, 1)), matrix(c(2/3, 1), ncol = 1))
})

test_that('1, 0 | 0', {
  expect_equal(ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(1, emp), list(1, emp)), matrix(c(1), ncol = 1))
})

test_that('1, 1 | 0', {
  expect_equal(ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(c(1, 1), emp), list(c(1, 1), emp)), matrix(c(0.5, 0.5), ncol = 1))
})

test_that('0, 0 | 1', {
  expect_equal(ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(emp, 1), list(emp, 1)), matrix(c(1), ncol = 1))
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############################################################
