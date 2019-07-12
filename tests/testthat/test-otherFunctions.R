test_that("cateNamesFromThreshold works", {
  expect_identical(cateNamesFromThreshold(c(9, 7, 6, 0), "noMotif", ""),
                   c("above9", "7-9", "6-7", "0-6", "noMotif"))
  expect_identical(cateNamesFromThreshold(c(9, 7, 6, 0),
                                          "atTSS", "kb"),
                   c("above9kb", "7-9kb", "6-7kb", "0-6kb", "atTSS"))
})
