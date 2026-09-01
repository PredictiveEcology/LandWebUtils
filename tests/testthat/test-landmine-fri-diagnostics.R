## Fixture: `nz` FRI zones as horizontal bands, with per-zone control over how much of each
## band is flammable, how much burns, and how much carries cohorts.
.fri_fixture <- function(zones = c(50, 100), burnRate = c(1 / 50, 1 / 100),
                         flamFrac = 1, cohortFrac = 1, cohortFracInit = NULL,
                         treedFrac = NULL, nrow = 20L, ncol = 20L) {
  nz <- length(zones)
  stopifnot(nrow %% nz == 0L)
  band <- nrow / nz

  lthfc <- terra::rast(nrows = nrow, ncols = ncol, vals = rep(zones, each = band * ncol))
  flamFrac <- rep_len(flamFrac, nz)
  cohortFrac <- rep_len(cohortFrac, nz)
  burnRate <- rep_len(burnRate, nz)

  cohortFracInit <- rep_len(cohortFracInit %||% cohortFrac, nz)
  treedFrac <- rep_len(treedFrac %||% 0, nz)

  flamV <- burnV <- vegV <- veg0V <- lccV <- numeric(0)
  for (i in seq_len(nz)) {
    n <- band * ncol
    f <- rep(0L, n)
    f[seq_len(round(n * flamFrac[i]))] <- 1L
    flamV <- c(flamV, f)
    burnV <- c(burnV, rep(burnRate[i], n))
    fill <- function(frac) {
      v <- rep(NA_integer_, n)
      if (frac > 0) v[seq_len(round(n * frac))] <- 1L
      v
    }
    vegV <- c(vegV, fill(cohortFrac[i]))
    veg0V <- c(veg0V, fill(cohortFracInit[i]))
    l <- rep(1L, n) ## 1 = not treed
    if (treedFrac[i] > 0) l[seq_len(round(n * treedFrac[i]))] <- 210L
    lccV <- c(lccV, l)
  }

  rr <- function(v) terra::setValues(terra::rast(lthfc), v)
  list(
    lthfc = lthfc, flammableMap = rr(flamV), burnMap = rr(burnV),
    vegTypeMap = rr(vegV), vegTypeMapInit = rr(veg0V), lccMap = rr(lccV)
  )
}

.fri_call <- function(f, ...) {
  landmine_fri_metrics(
    lthfc = f$lthfc, flammableMap = f$flammableMap,
    meanAnnualCumulBurnMap = f$burnMap, studyAreaName = "test",
    pixelAreaHa = 1, vegTypeMap = f$vegTypeMap, ...
  )
}

## ---- the ratio and its orientation --------------------------------------------------------------

test_that("a zone burning exactly at its target has ratio 1 and status ok", {
  m <- .fri_call(.fri_fixture())
  expect_equal(m$ratio, c(1, 1))
  expect_equal(m$status, c("ok", "ok"))
})

test_that("RATIO ORIENTATION: above 1 means the zone burns too RARELY", {
  ## zone 50 burns at 1/200 -> achieved FRI 200 -> 4x under
  m <- .fri_call(.fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100)))
  expect_equal(m[LTHFC == 50]$ratio, 4)
  expect_equal(m[LTHFC == 50]$status, "under (severe)")
  expect_equal(m[LTHFC == 100]$status, "ok")
})

test_that("burning too often gives a ratio below 1", {
  m <- .fri_call(.fri_fixture(zones = c(50, 100), burnRate = c(1 / 10, 1 / 100)))
  expect_equal(m[LTHFC == 50]$ratio, 0.2)
  expect_equal(m[LTHFC == 50]$status, "over (severe)")
})

test_that("tolerance and severe boundaries are respected", {
  ## ratio 1.3 is outside the default tolerance but inside severe
  m <- .fri_call(.fri_fixture(zones = c(100, 200), burnRate = c(1 / 130, 1 / 200)))
  expect_equal(m[LTHFC == 100]$status, "under")

  ## widening tolerance clears it
  m2 <- .fri_call(.fri_fixture(zones = c(100, 200), burnRate = c(1 / 130, 1 / 200)),
                  tolerance = c(0.5, 1.5))
  expect_equal(m2[LTHFC == 100]$status, "ok")
})

test_that("tolerance must sit inside severe", {
  expect_error(.fri_call(.fri_fixture(), tolerance = c(0.4, 1.25)), "severe")
  expect_error(.fri_call(.fri_fixture(), severe = c(0.9, 1.1)), "severe")
})

test_that("a zone that never burns is Inf, and is flagged rather than dropped", {
  m <- .fri_call(.fri_fixture(zones = c(50, 100), burnRate = c(0, 1 / 100)))
  expect_equal(m[LTHFC == 50]$FRI, Inf)
  expect_equal(m[LTHFC == 50]$status, "under (severe)")
})

## ---- the structural columns ---------------------------------------------------------------------

test_that("pctFlam is measured over the zone, not the study area", {
  m <- .fri_call(.fri_fixture(flamFrac = c(0.5, 1)))
  expect_equal(m$pctFlam, c(50, 100))
  expect_equal(m$nTotal, c(200L, 200L))
  expect_equal(m$nFlam, c(100L, 200L))
})

test_that("pctNoCohorts is measured over FLAMMABLE pixels only", {
  ## zone 1: half flammable, and only a quarter of the whole band has cohorts -- so of the
  ## flammable half (the first 100 cells), half carry cohorts
  m <- .fri_call(.fri_fixture(flamFrac = c(0.5, 1), cohortFrac = c(0.25, 1)))
  expect_equal(m$pctNoCohorts, c(50, 0))
})

test_that("pctNoCohorts is NA when no vegetation map is supplied", {
  f <- .fri_fixture()
  m <- landmine_fri_metrics(
    lthfc = f$lthfc, flammableMap = f$flammableMap, meanAnnualCumulBurnMap = f$burnMap,
    studyAreaName = "test", pixelAreaHa = 1
  )
  expect_true(all(is.na(m$pctNoCohorts)))
})

test_that("areaHa scales with pixel area", {
  f <- .fri_fixture()
  m1 <- .fri_call(f)
  m2 <- landmine_fri_metrics(
    lthfc = f$lthfc, flammableMap = f$flammableMap, meanAnnualCumulBurnMap = f$burnMap,
    studyAreaName = "test", pixelAreaHa = 5.76
  )
  expect_equal(m2$areaHa, m1$areaHa * 5.76)
})

test_that("patch stats separate one porous sheet from many islands", {
  ## a fully flammable band is one patch
  solid <- .fri_call(.fri_fixture(zones = c(50, 100)))
  expect_equal(solid$nPatches, c(1L, 1L))
  expect_equal(solid$largestPatchPct, c(100, 100))

  ## a checkerboard band is many patches under queen adjacency? no -- queen keeps it connected,
  ## which is exactly the distinction the column exists to draw
  lthfc <- terra::rast(nrows = 10L, ncols = 10L, vals = 50)
  chk <- terra::setValues(terra::rast(lthfc),
                          as.integer(outer(1:10, 1:10, function(i, j) (i + j) %% 2)))
  m <- landmine_fri_metrics(
    lthfc = lthfc, flammableMap = chk,
    meanAnnualCumulBurnMap = terra::setValues(terra::rast(lthfc), 1 / 50),
    studyAreaName = "test", pixelAreaHa = 1
  )
  expect_equal(m$pctFlam, 50)
  expect_equal(m$nPatches, 1L)
  expect_equal(m$largestPatchPct, 100)
})

test_that("AN ENTIRELY NON-FLAMMABLE ZONE IS REPORTED, not silently dropped", {
  ## `landmine_fri_summary()` masks by flammability first, so such a zone never reaches its
  ## output at all. It is carried here so the configuration is visible.
  m <- .fri_call(.fri_fixture(flamFrac = c(0, 1)))
  expect_equal(NROW(m), 2L)
  expect_equal(m[LTHFC == 50]$nFlam, 0L)
  expect_equal(m[LTHFC == 50]$nPatches, 0L)
  expect_true(is.na(m[LTHFC == 50]$largestPatchPct))
  expect_true(is.na(m[LTHFC == 50]$FRI))
  expect_equal(m[LTHFC == 50]$status, "not evaluated")
  expect_equal(m$studyArea, c("test", "test"))
})

## ---- the noise guard ----------------------------------------------------------------------------

test_that("a zone with too few expected ignitions is reported but not judged", {
  f <- .fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100))
  judged <- .fri_call(f, nYears = 1000)
  expect_equal(judged[LTHFC == 50]$status, "under (severe)")

  ## the same zone over a short run has too little opportunity to draw a conclusion from
  quiet <- .fri_call(f, nYears = 1, minIgnitions = 20L)
  expect_equal(quiet[LTHFC == 50]$status, "too few ignitions to judge")
})

test_that("the noise guard is off when nYears is unknown", {
  m <- .fri_call(.fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100)))
  expect_false(any(m$status == "too few ignitions to judge"))
  expect_true(all(is.na(m$expIgnitions)))
})

## ---- the verdict --------------------------------------------------------------------------------

test_that("an all-clear verdict says so in one line", {
  v <- landmine_fri_verdict(.fri_call(.fri_fixture()))
  expect_equal(v, "FRI check: 2/2 zones within tolerance.")
})

test_that("the verdict names failing zones worst-first, with both diagnostics", {
  f <- .fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100),
                    flamFrac = c(0.5, 1), cohortFrac = c(0.25, 1))
  v <- landmine_fri_verdict(.fri_call(f))
  expect_match(v, "1/2 zones within tolerance")
  expect_match(v, "50 \\(4\\.0x under, 50% flammable, 50% no cohorts\\)")
})

test_that("the verdict reports over-burning as a factor, not a fraction", {
  f <- .fri_fixture(zones = c(50, 100), burnRate = c(1 / 10, 1 / 100))
  expect_match(landmine_fri_verdict(.fri_call(f)), "5\\.0x over")
})

test_that("the verdict counts unjudged zones separately from failures", {
  f <- .fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100))
  v <- landmine_fri_verdict(.fri_call(f, nYears = 1, minIgnitions = 1e9))
  expect_match(v, "0/0 zones within tolerance \\(2 not judged\\)")
})

test_that("the verdict caps how many zones it names", {
  f <- .fri_fixture(zones = c(10, 20, 30, 40), burnRate = rep(1 / 500, 4), nrow = 20L)
  v <- landmine_fri_verdict(.fri_call(f), maxListed = 2L)
  expect_equal(lengths(regmatches(v, gregexpr("x under", v))), 2L)
})

## ---- the figures --------------------------------------------------------------------------------

test_that("the zone map builds, with SEPARATE fill scales per panel", {
  skip_if_not_installed("tidyterra")
  skip_if_not_installed("patchwork")
  f <- .fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100))
  m <- .fri_call(f)
  p <- landmine_plot_fri_zones(f$lthfc, m, "test")
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, 1L) ## plus the top-level plot = 2 panels

  ## a shared scale would put target intervals (tens of years) and log-ratios (about one
  ## unit) on one range, collapsing panel B to a flat block of colour
  scales <- lapply(list(p[[1]], p[[2]]), function(x) {
    x$scales$get_scales("fill")
  })
  expect_false(identical(scales[[1]]$limits, scales[[2]]$limits))
  ## zone 50 achieves 200 -> ratio 4 -> log2 = 2; the scale is symmetric about ratio 1
  expect_equal(scales[[2]]$limits, c(-2, 2))
})

test_that("the map carries the verdict as its caption", {
  skip_if_not_installed("tidyterra")
  skip_if_not_installed("patchwork")
  f <- .fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100))
  m <- .fri_call(f)
  p <- landmine_plot_fri_zones(f$lthfc, m, "test")
  expect_match(p$patches$annotation$caption, "within tolerance", fixed = TRUE)
})

test_that("the drivers plot builds with either driver column alone", {
  f <- .fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100), flamFrac = c(0.5, 1))
  expect_s3_class(landmine_plot_fri_drivers(.fri_call(f), "test"), "ggplot")

  ## with no vegetation map, pctNoCohorts is all NA and must be dropped, not plotted empty
  noVeg <- landmine_fri_metrics(
    lthfc = f$lthfc, flammableMap = f$flammableMap, meanAnnualCumulBurnMap = f$burnMap,
    studyAreaName = "test", pixelAreaHa = 1
  )
  p <- landmine_plot_fri_drivers(noVeg, "test")
  expect_s3_class(p, "ggplot")
  expect_equal(nlevels(p$data$driver), 1L)
})


## ---- treed land that was never populated with cohorts -------------------------------------------

test_that("cohortGap separates 'not forest' from 'forest that was never initialised'", {
  ## zone 1: 80% treed by LCC but only 20% ever got cohorts -> a 60 point gap
  ## zone 2: 20% treed and 20% initialised -> genuinely not forest, no gap
  f <- .fri_fixture(zones = c(50, 100), treedFrac = c(0.8, 0.2),
                    cohortFrac = c(0.2, 0.2), cohortFracInit = c(0.2, 0.2))
  m <- .fri_call(f, vegTypeMapInit = f$vegTypeMapInit, lccMap = f$lccMap,
                 treeClasses = c(81L, 210L, 220L, 230L, 240L))

  expect_equal(m$pctTreedLCC, c(80, 20))
  expect_equal(m$pctCohortsInit, c(20, 20))
  expect_equal(m$cohortGap, c(60, 0))
})

test_that("cohort LOSS is distinguished from never having been initialised", {
  ## started at 80%, ended at 20% -> 60 points lost, no initialisation gap
  f <- .fri_fixture(zones = c(50, 100), treedFrac = c(0.8, 0.8),
                    cohortFracInit = c(0.8, 0.8), cohortFrac = c(0.2, 0.8))
  m <- .fri_call(f, vegTypeMapInit = f$vegTypeMapInit, lccMap = f$lccMap,
                 treeClasses = 210L)

  expect_equal(m$pctCohortsLost, c(60, 0))
  expect_equal(m$cohortGap, c(0, 0))
})

test_that("the optional columns are NA rather than absent when their input is missing", {
  m <- .fri_call(.fri_fixture())
  expect_true(all(is.na(m$pctTreedLCC)))
  expect_true(all(is.na(m$cohortGap)))
  expect_true(all(is.na(m$pctCohortsLost)))
  expect_true(all(c("pctTreedLCC", "cohortGap", "pctCohortsLost") %in% names(m)))
})

test_that("lccMap without treeClasses is refused", {
  f <- .fri_fixture()
  expect_error(.fri_call(f, lccMap = f$lccMap), "treeClasses")
})

test_that("only the named tree classes count as treed", {
  f <- .fri_fixture(zones = c(50, 100), treedFrac = c(1, 0))
  m <- .fri_call(f, lccMap = f$lccMap, treeClasses = 999L)
  expect_equal(m$pctTreedLCC, c(0, 0))
})

## ---- the geometry guard -------------------------------------------------------------------------

test_that("A RASTER ON A DIFFERENT GRID IS REFUSED, not silently recycled", {
  ## this is the failure that produced a full table of plausible, wrong numbers when the
  ## dataPrep land cover raster (a larger extent) was read alongside the simulation grid
  f <- .fri_fixture()
  bigger <- terra::rast(nrows = 30L, ncols = 20L, vals = 1L)
  expect_snapshot(error = TRUE, .fri_call(f, lccMap = bigger, treeClasses = 210L))
  expect_snapshot(error = TRUE, .fri_call(f, vegTypeMapInit = bigger))
})

test_that("the drivers plot gains a cohortGap panel when it is available", {
  f <- .fri_fixture(zones = c(50, 100), burnRate = c(1 / 200, 1 / 100),
                    treedFrac = c(0.8, 0.2), cohortFrac = c(0.2, 0.2))
  m <- .fri_call(f, vegTypeMapInit = f$vegTypeMapInit, lccMap = f$lccMap, treeClasses = 210L)
  p <- landmine_plot_fri_drivers(m, "test")
  expect_equal(nlevels(p$data$driver), 3L)
})


## ---- the burnable-area denominator ---------------------------------------------------------------

## The fixture's zones are horizontal bands; this polygon keeps only the left `frac` of the grid,
## so each zone is clipped by the same fraction.
.fri_left_poly <- function(lthfc, frac) {
  e <- terra::ext(lthfc)
  terra::vect(terra::ext(e[1], e[1] + frac * (e[2] - e[1]), e[3], e[4]), crs = terra::crs(lthfc))
}

test_that("THE DENOMINATOR: unburnable pixels outside the study area inflate the achieved FRI", {
  ## every pixel inside the polygon burns at the target rate; outside, nothing burns
  lthfc <- terra::rast(nrows = 10L, ncols = 10L, vals = 100)
  poly <- .fri_left_poly(lthfc, 0.5)
  burnV <- rep(c(rep(1 / 100, 5), rep(0, 5)), 10) ## left half burns, right half does not
  f <- list(
    lthfc = lthfc,
    flammableMap = terra::setValues(terra::rast(lthfc), 1L),
    burnMap = terra::setValues(terra::rast(lthfc), burnV)
  )
  args <- list(lthfc = f$lthfc, flammableMap = f$flammableMap,
               meanAnnualCumulBurnMap = f$burnMap, studyAreaName = "t", pixelAreaHa = 1)

  ## without the polygon the dead half doubles the apparent interval
  bad <- do.call(landmine_fri_metrics, args)
  expect_equal(bad$ratio, 2)
  expect_equal(bad$status, "under") ## ratio is exactly 2, the severe boundary, which is inclusive

  ## with it, the zone is exactly on target
  good <- do.call(landmine_fri_metrics, c(args, list(studyArea = poly)))
  expect_equal(good$ratio, 1)
  expect_equal(good$status, "ok")

  ## and the uncorrected figure is retained so the discrepancy is explainable
  expect_equal(good$ratioUnmasked, 2)
  expect_equal(good$pctInStudyArea, 50)
})

test_that("pctFlam and the patch stats also use the burnable area only", {
  lthfc <- terra::rast(nrows = 10L, ncols = 10L, vals = 100)
  ## flammable only on the left half, which is also the burnable area
  flamV <- rep(c(rep(1L, 5), rep(0L, 5)), 10)
  m <- landmine_fri_metrics(
    lthfc = lthfc, flammableMap = terra::setValues(terra::rast(lthfc), flamV),
    meanAnnualCumulBurnMap = terra::setValues(terra::rast(lthfc), 1 / 100),
    studyAreaName = "t", pixelAreaHa = 1, studyArea = .fri_left_poly(lthfc, 0.5)
  )
  expect_equal(m$pctFlam, 100) ## 50% of the grid, but 100% of the burnable part
  expect_equal(m$nTotal, 50L)
})

test_that("the verdict warns when a zone is largely outside the burnable area", {
  lthfc <- terra::rast(nrows = 10L, ncols = 10L, vals = 100)
  m <- landmine_fri_metrics(
    lthfc = lthfc, flammableMap = terra::setValues(terra::rast(lthfc), 1L),
    meanAnnualCumulBurnMap = terra::setValues(terra::rast(lthfc), 1 / 100),
    studyAreaName = "t", pixelAreaHa = 1, studyArea = .fri_left_poly(lthfc, 0.2)
  )
  expect_match(landmine_fri_verdict(m), "lie partly outside the burnable area")
  expect_match(landmine_fri_verdict(m), "20% inside")
})

test_that("no warning, and NA columns, when no study area is given", {
  m <- .fri_call(.fri_fixture())
  expect_true(all(is.na(m$pctInStudyArea)))
  expect_true(all(is.na(m$ratioUnmasked)))
  expect_false(grepl("burnable area", landmine_fri_verdict(m)))
})
