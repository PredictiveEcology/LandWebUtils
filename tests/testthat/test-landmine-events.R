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

test_that("STUDY AREA: pixels the fire model cannot reach are excluded from the denominator", {
  ## left half burns at the target rate; the right half is unburnable because it lies outside
  ## the polygon LandMine masks `ROSmap` and `mod$spreadProb` to
  lthfc <- terra::rast(nrows = 10L, ncols = 10L, vals = 100)
  flam <- terra::setValues(terra::rast(lthfc), 1L)
  burn <- terra::setValues(terra::rast(lthfc), rep(c(rep(1 / 100, 5), rep(0, 5)), 10))
  e <- terra::ext(lthfc)
  poly <- terra::vect(terra::ext(e[1], e[1] + 0.5 * (e[2] - e[1]), e[3], e[4]),
                      crs = terra::crs(lthfc))

  expect_equal(landmine_fri_summary(lthfc, flam, burn, "t")$FRI, 200)
  expect_equal(landmine_fri_summary(lthfc, flam, burn, "t", studyArea = poly)$FRI, 100)
})

test_that("a zone lying entirely outside the study area drops out of the summary", {
  lthfc <- terra::rast(nrows = 10L, ncols = 10L, vals = rep(c(rep(50, 5), rep(170, 5)), 10))
  flam <- terra::setValues(terra::rast(lthfc), 1L)
  burn <- terra::setValues(terra::rast(lthfc), 1 / 50)
  e <- terra::ext(lthfc)
  poly <- terra::vect(terra::ext(e[1], e[1] + 0.5 * (e[2] - e[1]), e[3], e[4]),
                      crs = terra::crs(lthfc))

  expect_equal(landmine_fri_summary(lthfc, flam, burn, "t")$LTHFC, c(50, 170))
  expect_equal(landmine_fri_summary(lthfc, flam, burn, "t", studyArea = poly)$LTHFC, 50)
})

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

# ---- landmine_ignition_budget -----------------------------------------------------------

test_that("the budget is (area / meanFireSize) / FRI over FLAMMABLE pixels", {
  fri <- fri_rast(rep(c(30, 70), each = 8))
  flam <- fri_rast(rep(1, 16))
  ## 1 x 1 map units per pixel -> 1e-4 ha; use kBest/upper that give a known mean fire size
  out <- landmine_ignition_budget(fri, flam, kBest = 0.5, biggestPossibleFireSizeHa = 1e5)

  mfs <- meanTruncPareto(k = 0.5, lower = 1, upper = 1e5, alpha = 1)
  haPerPix <- 1 / 1e4
  expect_equal(unname(out$numFiresPerYear[c("30", "70")]),
               c((8 * haPerPix / mfs) / 30, (8 * haPerPix / mfs) / 70))
})

test_that("non-flammable and zero-FRI pixels are masked OUT of the raster and the budget", {
  ## Consequence: they can never be ignition locations either, since the start-cell pool is
  ## built from this same masked raster.
  fri <- fri_rast(c(rep(30, 4), rep(0, 4), rep(70, 4), rep(30, 4)))
  flam <- fri_rast(c(rep(1, 12), rep(0, 4))) ## last 4 (zone 30) non-flammable

  out <- landmine_ignition_budget(fri, flam, kBest = 0.5, biggestPossibleFireSizeHa = 1e5)
  vals <- terra::values(out$fireReturnInterval, mat = FALSE)

  expect_true(all(is.na(vals[5:8])))    ## zero-FRI masked
  expect_true(all(is.na(vals[13:16])))  ## non-flammable masked
  expect_identical(sort(unique(na.omit(vals))), c(30, 70))
  ## zone 30 keeps only its 4 flammable pixels, so its budget is half zone 70's x (70/30)
  expect_equal(unname(out$numFiresPerYear["30"] / out$numFiresPerYear["70"]), 70 / 30)
})

test_that("STUDY AREA: the budget excludes ground the fire model cannot burn", {
  ## zone 50 spans the whole grid but only its left half is inside the burnable polygon
  fri <- terra::rast(nrows = 10L, ncols = 10L, vals = 50)
  flam <- terra::setValues(terra::rast(fri), 1L)
  e <- terra::ext(fri)
  poly <- terra::vect(terra::ext(e[1], e[1] + 0.5 * (e[2] - e[1]), e[3], e[4]),
                      crs = terra::crs(fri))

  now <- landmine_ignition_budget(fri, flam, kBest = 0.731, biggestPossibleFireSizeHa = 1e6)
  fixed <- landmine_ignition_budget(fri, flam, kBest = 0.731, biggestPossibleFireSizeHa = 1e6,
                                    studyArea = poly)

  ## the budget is proportional to burnable area, so halving the area halves the fires
  expect_equal(unname(fixed$numFiresPerYear["50"]), unname(now$numFiresPerYear["50"]) / 2)

  ## and the returned raster -- which is also the start-cell pool -- is masked, so ignitions
  ## can no longer be drawn outside the polygon
  expect_equal(sum(!is.na(terra::values(fixed$fireReturnInterval, mat = FALSE))), 50L)
  expect_equal(sum(!is.na(terra::values(now$fireReturnInterval, mat = FALSE))), 100L)
})

test_that("the ordering contract survives study-area masking", {
  fri <- terra::rast(nrows = 10L, ncols = 10L, vals = rep(c(rep(50, 5), rep(170, 5)), 10))
  flam <- terra::setValues(terra::rast(fri), 1L)
  e <- terra::ext(fri)
  ## keep the left half: zone 50 survives in full, zone 170 disappears entirely
  poly <- terra::vect(terra::ext(e[1], e[1] + 0.5 * (e[2] - e[1]), e[3], e[4]),
                      crs = terra::crs(fri))
  out <- landmine_ignition_budget(fri, flam, kBest = 0.731, biggestPossibleFireSizeHa = 1e6,
                                  studyArea = poly)

  expect_equal(out$fireReturnIntervalsByPolygonNumeric, c(50, NA))
  expect_true(is.na(out$numFiresPerYear[[2]])) ## NA-VALUED, so na.omit() drops it downstream
})

test_that("THE ORDERING CONTRACT: zones ascending, NA row last, and NA-valued so na.omit drops it", {
  ## The reburn loop indexes fire counts by GROUP POSITION of a polygonNumeric-keyed table.
  ## That is valid only if all three of these hold. Break any one and zones silently receive
  ## each other's fire counts -- no error, just a wrong fire regime.
  fri <- fri_rast(c(rep(250, 4), rep(60, 4), rep(120, 4), rep(NA_real_, 4))) ## deliberately unsorted
  flam <- fri_rast(rep(1, 16))

  out <- landmine_ignition_budget(fri, flam, kBest = 0.5, biggestPossibleFireSizeHa = 1e5)
  fri_by_poly <- out$fireReturnIntervalsByPolygonNumeric

  ## 1. ascending, regardless of the order values appear in the raster
  expect_identical(fri_by_poly[!is.na(fri_by_poly)], c(60, 120, 250))
  ## 2. the NA row is LAST
  expect_identical(which(is.na(fri_by_poly)), length(fri_by_poly))
  ## 3. the NA entry's VALUE is NA, so na.omit() drops it (an NA *name* alone would not)
  nf <- out$numFiresPerYear
  expect_true(is.na(nf[length(nf)]))
  expect_identical(names(na.omit(nf)), c("60", "120", "250"))

  ## and the contract the loop actually relies on: keyed-group order == names(numFiresPerYear)
  pool <- data.table::data.table(
    pixel = seq_along(terra::values(out$fireReturnInterval, mat = FALSE)),
    polygonNumeric = terra::values(out$fireReturnInterval, mat = FALSE)
  )
  pool <- na.omit(pool)[polygonNumeric > 0]
  data.table::setkeyv(pool, "polygonNumeric")
  expect_identical(as.character(pool[, unique(polygonNumeric)]), names(na.omit(nf)))
})

test_that("a raster with no NA cells still yields a trailing NA row", {
  fri <- fri_rast(rep(c(30, 70), each = 8))
  flam <- fri_rast(rep(1, 16))
  out <- landmine_ignition_budget(fri, flam, kBest = 0.5, biggestPossibleFireSizeHa = 1e5)

  fbp <- out$fireReturnIntervalsByPolygonNumeric
  expect_identical(which(is.na(fbp)), length(fbp))
  expect_identical(names(na.omit(out$numFiresPerYear)), c("30", "70"))
})

test_that("a single-zone raster works", {
  fri <- fri_rast(rep(40, 16))
  flam <- fri_rast(rep(1, 16))
  out <- landmine_ignition_budget(fri, flam, kBest = 0.5, biggestPossibleFireSizeHa = 1e5)
  expect_identical(names(na.omit(out$numFiresPerYear)), "40")
})

test_that("an entirely non-flammable landscape yields an empty budget, not an error", {
  fri <- fri_rast(rep(c(30, 70), each = 8))
  flam <- fri_rast(rep(0, 16))
  out <- landmine_ignition_budget(fri, flam, kBest = 0.5, biggestPossibleFireSizeHa = 1e5)

  expect_length(na.omit(out$numFiresPerYear), 0L)
  expect_true(all(is.na(terra::values(out$fireReturnInterval, mat = FALSE))))
})

# ---- landmine_reburn_budget -------------------------------------------------------------

## The inline logic this replaces, kept verbatim as the equivalence oracle.
.reburn_inline <- function(tooSmallByPoly, friByPolygon, remainingSize = NULL) {
  p <- data.table::copy(tooSmallByPoly)
  p <- p[, N := .N, by = polygonNumeric]
  if (!is.null(remainingSize)) p <- p[, maxSize := remainingSize]
  p <- p[data.table::data.table(polygonNumeric = friByPolygon), on = "polygonNumeric"]
  p[is.na(N), N := 0]
  if ("pixel" %in% names(p)) data.table::set(p, NULL, "pixel", NULL)
  list(numFiresThisPeriod = p[, N[1], by = "polygonNumeric"]$V1,
       fireSizesInPixels = na.omit(p)$maxSize)
}

test_that("output order is driven by friByPolygon: zones ascending, NA row LAST", {
  ## This is what makes the caller's `numFiresThisPeriod[.GRP]` indexing correct.
  tsbp <- data.table::data.table(
    pixel = c(11L, 12L, 13L, 14L),
    polygonNumeric = c(170, 30, 170, 55),
    maxSize = c(500, 10, 700, 90)
  )
  fbp <- c(30, 55, 170, NA)

  out <- landmine_reburn_budget(tsbp, fbp)
  expect_identical(out$polysNeedMoreFires$polygonNumeric, c(30, 55, 170, 170, NA))
  expect_identical(out$numFiresThisPeriod, c(1L, 1L, 2L, 0L)) ## zone order + NA row last
})

test_that("zones with nothing short get N = 0 and are dropped from fireSizesInPixels", {
  ## na.omit() on the whole table is load-bearing: it drops both the unmatched zones and the
  ## NA-FRI row, which is what keeps fireSizesInPixels paired with the start cells.
  tsbp <- data.table::data.table(
    pixel = c(11L, 12L), polygonNumeric = c(30, 170), maxSize = c(10, 500)
  )
  fbp <- c(30, 55, 170, NA) ## zone 55 has no short fires

  out <- landmine_reburn_budget(tsbp, fbp)
  expect_identical(out$numFiresThisPeriod, c(1L, 0L, 1L, 0L))
  expect_identical(out$fireSizesInPixels, c(10, 500)) ## zone 55 and the NA row absent
})

test_that("phase 2 replaces maxSize with the REMAINING shortfall", {
  tsbp <- data.table::data.table(
    pixel = c(11L, 12L), polygonNumeric = c(30, 170), maxSize = c(100, 500)
  )
  fbp <- c(30, 170, NA)

  ph1 <- landmine_reburn_budget(tsbp, fbp)
  ph2 <- landmine_reburn_budget(tsbp, fbp, remainingSize = c(40, 120))

  expect_identical(ph1$fireSizesInPixels, c(100, 500)) ## full original targets
  expect_identical(ph2$fireSizesInPixels, c(40, 120))  ## what is LEFT to burn
})

test_that("a zero remaining shortfall survives to the caller's >0 filter", {
  tsbp <- data.table::data.table(
    pixel = c(11L, 12L), polygonNumeric = c(30, 170), maxSize = c(100, 500)
  )
  out <- landmine_reburn_budget(tsbp, c(30, 170, NA), remainingSize = c(0, 120))
  expect_identical(out$fireSizesInPixels, c(0, 120))
})

test_that("the input table is NOT modified by reference", {
  ## The inline version mutated tooSmallByPoly via `:=`.
  tsbp <- data.table::data.table(
    pixel = 11L, polygonNumeric = 30, maxSize = 100
  )
  before <- data.table::copy(tsbp)
  invisible(landmine_reburn_budget(tsbp, c(30, NA)))
  expect_identical(tsbp, before)
})

test_that("EQUIVALENCE with the inline logic over randomised multi-zone inputs", {
  ## The blast radius of an ordering mistake here is a plausible-but-wrong fire regime, so
  ## compare against the original code path rather than against hand-computed expectations.
  withr::local_seed(11)
  zones <- c(30, 45, 55, 90, 170)
  for (i in seq_len(200)) {
    nz <- sample(2:5, 1)
    zs <- sort(sample(zones, nz))
    nfire <- sample(1:12, 1)
    tsbp <- data.table::data.table(
      pixel = sample.int(1000L, nfire),
      polygonNumeric = sample(zs, nfire, replace = TRUE), ## interleaved across zones
      maxSize = sample(1:900, nfire, replace = TRUE)
    )
    data.table::setkeyv(tsbp, NULL)
    fbp <- c(zs, NA)
    rem <- if (i %% 2 == 0) sample(0:400, nfire, replace = TRUE) else NULL

    got <- landmine_reburn_budget(tsbp, fbp, remainingSize = rem)
    want <- .reburn_inline(tsbp, fbp, remainingSize = rem)

    expect_identical(got$numFiresThisPeriod, want$numFiresThisPeriod)
    expect_identical(got$fireSizesInPixels, want$fireSizesInPixels)
  }
})

## ---- the fire-size distribution fit --------------------------------------------------------------

test_that("landmine_area_share computes the share carried by the largest fires", {
  ## 5 fires; the largest 50% is anything above the median (1) -- so just the 100
  expect_equal(landmine_area_share(c(1, 1, 1, 1, 100), 0.5), 100 / 104)
  ## an all-equal set has nothing strictly above its quantile
  expect_equal(landmine_area_share(rep(5, 10), 0.9), 0)
})

test_that("landmine_area_share is NA rather than NaN for degenerate input", {
  expect_true(is.na(landmine_area_share(numeric(0))))
  expect_true(is.na(landmine_area_share(c(0, 0, 0))))
})

test_that("THE RULE OF THUMB: the fitted k puts ~95% of area in the largest fires", {
  skip_if_not_installed("VGAM")
  set.seed(42)
  k <- landmine_estimate_kBest(1e6, nDraws = 2e5)
  expect_gt(k, 0.05)
  expect_lt(k, 0.99)

  ## the property the fit exists to produce, checked on a fresh draw
  set.seed(43)
  fs <- round(VGAM::rtruncpareto(2e5, 1, upper = 1e6, shape = k))
  expect_equal(landmine_area_share(fs, 0.95), 0.95, tolerance = 0.02)
})

test_that("a different target area share moves the fit", {
  skip_if_not_installed("VGAM")
  set.seed(7)
  loose <- landmine_estimate_kBest(1e6, nDraws = 5e4, targetAreaShare = 0.5)
  set.seed(7)
  tight <- landmine_estimate_kBest(1e6, nDraws = 5e4, targetAreaShare = 0.95)
  expect_false(isTRUE(all.equal(loose, tight)))
})

test_that("the fit is reproducible under a seed, and stochastic without one", {
  skip_if_not_installed("VGAM")
  set.seed(99)
  a <- landmine_estimate_kBest(1e6, nDraws = 2e4)
  set.seed(99)
  b <- landmine_estimate_kBest(1e6, nDraws = 2e4)
  expect_identical(a, b)

  c1 <- landmine_estimate_kBest(1e6, nDraws = 2e4)
  expect_false(isTRUE(all.equal(a, c1)))
})

test_that("a degenerate search interval is refused", {
  expect_error(landmine_estimate_kBest(1e6, interval = 0.4), "interval")
  expect_error(landmine_estimate_kBest(1e6, interval = c(0.9, 0.1)), "interval")
  expect_error(landmine_estimate_kBest(1e6, topFireQuantile = 1), "topFireQuantile")
  expect_error(landmine_estimate_kBest(1e6, targetAreaShare = 0), "targetAreaShare")
})
