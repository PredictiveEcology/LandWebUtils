## Helper: a small FRI raster with `nz` zones laid out in horizontal bands.
fri_rast <- function(vals, nrow = 4, ncol = 4) {
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol, ymin = 0, ymax = nrow)
  terra::values(r) <- vals
  r
}

# ---- landmine_area_burned_by_zone -------------------------------------------------------

test_that("a FULLY burned zone is counted, not reported as zero", {
  ## The regression this locks down: counting with `table(x)[2]` assumes the levels are
  ## always {0, 1}. A zone that burns completely has only ONE level, so `[2]` is NA, and the
  ## module's `npix[is.na(npix)] <- 0` then turned a total burn into "zero hectares burned".
  fri <- fri_rast(rep(c(30, 70), each = 8))
  burn <- fri_rast(c(rep(1, 8), rep(0, 8))) ## zone 30 fully burned, zone 70 not at all

  out <- landmine_area_burned_by_zone(burn, fri, time = 5, pixelAreaHa = 5.76)

  expect_identical(out$nPixelsBurned[out$FRI == "30"], 8)
  expect_identical(out$nPixelsBurned[out$FRI == "70"], 0)
  expect_equal(out$haBurned[out$FRI == "30"], 8 * 5.76)
})

test_that("partial and empty burns are counted correctly", {
  fri <- fri_rast(rep(c(30, 70), each = 8))
  burn <- fri_rast(c(1, 1, 1, rep(0, 5), rep(0, 8))) ## 3 of zone 30, none of zone 70

  out <- landmine_area_burned_by_zone(burn, fri, time = 1, pixelAreaHa = 1)
  expect_identical(out$nPixelsBurned, c(3, 0))
  expect_identical(as.character(out$FRI), c("30", "70"))
  expect_identical(unique(out$time), 1)
})

test_that("NA FRI pixels form no zone of their own", {
  fri <- fri_rast(c(rep(30, 8), rep(NA_real_, 8)))
  burn <- fri_rast(rep(1, 16))

  out <- landmine_area_burned_by_zone(burn, fri, time = 1, pixelAreaHa = 1)
  expect_identical(nrow(out), 1L)
  expect_identical(out$nPixelsBurned, 8)
})

# ---- landmine_fri_summary ---------------------------------------------------------------

test_that("achieved FRI is the reciprocal of the mean annual burn rate, per zone", {
  fri <- fri_rast(rep(c(30, 70), each = 8))
  flam <- fri_rast(rep(1, 16))
  ## zone 30 burns at 1/30 per year, zone 70 at 1/70 -> both exactly on target
  burn <- fri_rast(c(rep(1 / 30, 8), rep(1 / 70, 8)))

  out <- landmine_fri_summary(fri, flam, burn, "TestArea")

  expect_identical(out$LTHFC, c(30, 70))
  expect_equal(out$FRI, c(30, 70))
  expect_identical(unique(out$studyArea), "TestArea")
})

test_that("non-flammable pixels are excluded from BOTH numerator and denominator", {
  ## Half of zone 30 is non-flammable and never burns. If those pixels were left in the
  ## denominator the achieved FRI would read as twice the target.
  fri <- fri_rast(rep(c(30, 70), each = 8))
  flam <- fri_rast(c(rep(1, 4), rep(0, 4), rep(1, 8)))
  burn <- fri_rast(c(rep(1 / 30, 4), rep(0, 4), rep(1 / 70, 8)))

  out <- landmine_fri_summary(fri, flam, burn, "TestArea")
  expect_equal(out$FRI, c(30, 70)) ## NOT c(60, 70)
})

test_that("a zone that never burns yields Inf (documented, not 'fixed')", {
  fri <- fri_rast(rep(c(30, 70), each = 8))
  flam <- fri_rast(rep(1, 16))
  burn <- fri_rast(c(rep(1 / 30, 8), rep(0, 8))) ## zone 70 never burns

  out <- landmine_fri_summary(fri, flam, burn, "TestArea")
  expect_equal(out$FRI[1], 30)
  expect_identical(out$FRI[2], Inf)
})

test_that("an NA burn pixel not shared by the flammable mask makes that zone NA", {
  ## Documents the implicit contract: the NA masks of the burn map and lthfc must agree.
  ## `sum()` has no `na.rm`, so one stray NA nulls the whole zone.
  fri <- fri_rast(rep(c(30, 70), each = 8))
  flam <- fri_rast(rep(1, 16))
  burn <- fri_rast(c(NA_real_, rep(1 / 30, 7), rep(1 / 70, 8)))

  out <- landmine_fri_summary(fri, flam, burn, "TestArea")
  expect_identical(out$FRI[1], NA_real_)
  expect_equal(out$FRI[2], 70)
})

# ---- landmine_draw_num_fires ------------------------------------------------------------

test_that("zone names survive the draw", {
  ## `rnbinom()` drops names, and the caller indexes fire counts by zone -- losing them
  ## silently assigns zones each other's counts.
  nf <- c("30" = 2, "70" = 0.5, "170" = 0.05)
  withr::local_seed(42)
  out <- landmine_draw_num_fires(nf)

  expect_identical(names(out), c("30", "70", "170"))
  expect_length(out, 3L)
  expect_true(all(out >= 0))
})

test_that("the draw is overdispersed enough to give a low-rate zone many zero years", {
  ## This is the property that makes a long-FRI zone episodic rather than steady.
  withr::local_seed(1)
  draws <- landmine_draw_num_fires(rep(0.05, 5000))
  expect_gt(mean(draws == 0), 0.9)
})

test_that("fireTimestep scales the mean", {
  withr::local_seed(7)
  a <- mean(landmine_draw_num_fires(rep(1, 20000), fireTimestep = 1))
  withr::local_seed(7)
  b <- mean(landmine_draw_num_fires(rep(1, 20000), fireTimestep = 5))
  expect_gt(b, 3 * a)
})

# ---- landmine_sizes_to_pixels -----------------------------------------------------------

test_that("probabilistic rounding preserves the expectation", {
  ## The entire justification for the construct: truncating would systematically lose every
  ## sub-pixel fire and shrink the mean.
  withr::local_seed(3)
  sizesHa <- rep(3.25, 1e5)
  px <- landmine_sizes_to_pixels(sizesHa, pixelAreaHa = 6.25)

  expect_equal(mean(px), 3.25 / 6.25, tolerance = 0.01)
  expect_true(all(px %in% c(0, 1)))
  expect_lt(mean(trunc(sizesHa / 6.25)), mean(px)) ## truncation would lose it entirely
})

test_that("exact pixel multiples never round up", {
  withr::local_seed(5)
  px <- landmine_sizes_to_pixels(rep(c(6.25, 12.5), 500), pixelAreaHa = 6.25)
  expect_identical(sort(unique(px)), c(1, 2))
})

test_that("zero-length input round-trips", {
  expect_length(landmine_sizes_to_pixels(numeric(0), pixelAreaHa = 6.25), 0L)
})

# ---- landmine_reburn_ceiling ------------------------------------------------------------

test_that("only zones with outstanding fires are reported", {
  dt <- data.table::data.table(
    polygonNumeric = c(30, 30, 55, 170, 100),
    N = c(2L, 2L, 3L, 1L, 0L),
    maxSize = c(100, 250, 180, 5000, NA_real_)
  )
  out <- landmine_reburn_ceiling(dt, year = 17)

  expect_identical(out$FRI, c(30, 55, 170))   ## zone 100 (N = 0) excluded
  expect_identical(out$nFires, c(2L, 3L, 1L))
  expect_identical(out$pixelsShort, c(350L, 180L, 5000L))
  expect_identical(unique(out$year), 17)
})

test_that("NULL, empty, and nothing-short inputs give zero rows, not an error", {
  expect_identical(nrow(landmine_reburn_ceiling(NULL, year = 1)), 0L)

  empty <- data.table::data.table(
    polygonNumeric = numeric(0), N = integer(0), maxSize = numeric(0)
  )
  expect_identical(nrow(landmine_reburn_ceiling(empty, year = 1)), 0L)

  none <- data.table::data.table(polygonNumeric = c(30, 70), N = c(0L, 0L), maxSize = c(0, 0))
  expect_identical(nrow(landmine_reburn_ceiling(none, year = 1)), 0L)
})

test_that("the NA-FRI row never reaches the output", {
  ## `friByPolyDT` carries an NA-FRI row; a `FRI=NA` line in the diagnostic log would be
  ## meaningless and would mislead a maxReburns decision.
  dt <- data.table::data.table(
    polygonNumeric = c(30, NA_real_),
    N = c(1L, 4L),
    maxSize = c(10, 999)
  )
  out <- landmine_reburn_ceiling(dt, year = 3)
  expect_identical(out$FRI, 30)
})
