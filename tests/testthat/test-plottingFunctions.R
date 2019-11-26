test_that("barplotFromNamedVector works", {
  topTs <- c(5, 10, 15, 20) * 1e3
  colorsForGraphs <- rainbow(2 + length(topTs))
  names(colorsForGraphs) <- c(paste0("overlapTop", sort(topTs / 1e3), "k"),
                              "overlap", "specific")
  v <- c(15809, rep(5000, 4), c(9316, 4762, 4327, 3831, 5051, 36638))
  names(v) <- paste(c(rep("ChIP1", 5), rep("ChIP2", 6)),
                    c(rep(paste0("overlap", c("", "Top5k", "Top10k",
                                              "Top15k", "Top20k")), 2),
                      "specific"),
                    sep = ".")
  b <- barplotFromNamedVector(v, colorsForGraphs = colorsForGraphs)
  expect_identical(b, c(0.7, 1.9))
})
