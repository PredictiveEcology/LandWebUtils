comps <- function() {
  data.frame(
    fireSize = c(9, 35, 104, 260, 521),
    nBurned  = c(9, 35, 104, 260, 521),
    areaHa   = c(51.8, 201.6, 599.0, 1497.6, 3001.0),
    shapeObs = c(1.32, 2.86, 3.87, 4.37, 7.69),
    shapeExp = c(2.12, 2.93, 4.21, 5.94, 7.76),
    pctIslandsObs = c(0, 5.41, 4.59, 10.65, 12.73),
    pctIslandsExp = c(4, 4, 6, 9, 9)
  )
}

test_that("landmine_plot_calibShape builds", {
  p <- landmine_plot_calibShape(comps())
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("landmine_plot_calibIslands builds", {
  p <- landmine_plot_calibIslands(comps())
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("landmine_plot_calibDiscrimination builds and reports separation", {
  set.seed(1)
  scores <- data.frame(
    group = rep(c("fitted", "random"), each = 6),
    objective = c(rnorm(6, 1.9, 0.3), rnorm(6, 3.7, 0.8))
  )
  p <- landmine_plot_calibDiscrimination(scores)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_match(p$labels$subtitle, "separation: [0-9.]+ sd")
})

test_that("landmine_plot_calibDiscrimination requires the expected columns", {
  expect_snapshot(error = TRUE, landmine_plot_calibDiscrimination(data.frame(a = 1)))
})

test_that("landmine_plot_calibParams builds and normalises to the search interval", {
  lower <- c(-2, -3, -3, -3, 1, 3.5, 0.75)
  upper <- c(-0.1, -0.5, -0.5, -1, 3.5, 5, 1)
  set.seed(2)
  params <- t(replicate(8, runif(7, lower, upper)))
  colnames(params) <- paste0("par", 1:7)

  p <- landmine_plot_calibParams(params, lower, upper)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  ## normalised positions must lie in [0, 1]
  expect_true(all(p$data$value >= 0 & p$data$value <= 1))
  expect_equal(nrow(p$data), 8 * 7)
})

test_that("landmine_plot_calibConvergence builds and reports where improvement stopped", {
  ## monotone non-increasing trace that flattens at iteration 12 of 30
  bv <- c(seq(1, 0.2, length.out = 12), rep(0.2, 18))
  p <- landmine_plot_calibConvergence(bv)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_match(p$labels$subtitle, "last improvement at iteration 12 of 30")
  expect_match(p$labels$subtitle, "18 further iterations gained nothing")
})

test_that("landmine_plot_calibConvergence accepts a DEoptim-shaped list", {
  p <- landmine_plot_calibConvergence(list(member = list(bestvalit = c(1, 0.5, 0.5))))
  expect_s3_class(p, "ggplot")
})

test_that("landmine_plot_calibConvergence errors informatively on the wrong input", {
  expect_snapshot(error = TRUE, landmine_plot_calibConvergence(list(a = 1)))
})
