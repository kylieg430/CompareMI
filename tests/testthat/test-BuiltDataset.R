test_that("Dataset generated with desired number of samples", {
  expect_equal(nrow(generateData(100)),
               100)
})
