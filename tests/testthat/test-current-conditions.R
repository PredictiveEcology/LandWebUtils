## Fixtures -----------------------------------------------------------------------------------------

## 10 x 10 grid of 30 m cells in a real projected CRS, so reprojection tests are meaningful.
cc_grid <- function(vals = NA_real_, crs = "EPSG:3979") {
  terra::rast(
    nrows = 10L, ncols = 10L, xmin = -1e6, xmax = -1e6 + 300, ymin = 5e5, ymax = 5e5 + 300,
    crs = crs, vals = vals
  )
}

## Rectangle covering columns `cols` of `cc_grid()`, as a SpatVector.
cc_poly <- function(cols = 1:10, crs = "EPSG:3979") {
  terra::vect(
    terra::ext(-1e6 + 30 * (min(cols) - 1), -1e6 + 30 * max(cols), 5e5, 5e5 + 300),
    crs = crs
  )
}

## column index of every cell, in terra's cell order (row-major)
cc_col <- rep(1:10, times = 10)

## lcc2020_classes -----------------------------------------------------------------------------------

test_that("lcc2020_classes() groups are disjoint integer codes", {
  cls <- lcc2020_classes()
  expect_identical(cls$forest, c(1L, 2L, 5L, 6L))
  expect_identical(cls$nonFlammable, c(16L, 18L, 19L))
  expect_identical(cls$urban, 17L)
  all <- unlist(cls, use.names = FALSE)
  expect_type(all, "integer")
  expect_identical(anyDuplicated(all), 0L)
})

## lcc2020_remap_table --------------------------------------------------------------------------------

## VERBATIM from LandWeb_preamble 1.0.6 InitMaps(), with its local class vectors inlined. This is the
## behaviour the promoted function must reproduce exactly.
remap_inline_preamble_1.0.6 <- function(uniqueLCCClasses, LandTypeCCvals, treeClassesToReplace) {
  treeClassesCC <- c(1L, 2L, 5L, 6L)
  nonFlammClassesCC <- c(16L, 18L, 19L)
  urbanClassCC <- 17L
  remapDT <- expand.grid(
    LCC = c(NA_integer_, sort(uniqueLCCClasses)),
    CC = c(NA_integer_, sort(unique(na.omit(LandTypeCCvals))))
  ) |>
    data.table::as.data.table()
  remapDT[LCC %in% c(0, 20, 30), newLCC := NA_integer_]
  remapDT[is.na(CC) | CC == 15L, newLCC := LCC]
  remapDT[CC %in% nonFlammClassesCC, newLCC := NA_integer_]
  remapDT[CC %in% c(treeClassesCC, 8L, 10L, 11L, 12L, 13L, 14L), newLCC := LCC]
  remapDT[is.na(LCC) & CC %in% treeClassesCC, newLCC := 99]
  remapDT[CC == urbanClassCC, newLCC := 99]
  remapDT[LCC %in% treeClassesToReplace, newLCC := 99]
  remapDT
}

## plain data.frame of integer columns: compares values, not data.table internals
as_int_df <- function(x) {
  as.data.frame(lapply(as.data.frame(x), as.integer))
}

test_that("lcc2020_remap_table() reproduces the preamble's inline remap exactly", {
  scanfi <- c(20L, 30L, 40L, 50L, 80L, 81L, 100L, 210L, 220L, 230L, 240L)
  legend <- c(1L, 2L, 5L, 6L, 8L, 10:19)
  withr::local_seed(20260916)
  for (i in seq_len(200)) {
    lcc <- sort(sample(scanfi, sample(1:length(scanfi), 1)))
    cc <- sample(c(legend, NA_integer_), sample(1:30, 1), replace = TRUE)
    replace <- sample(c(scanfi, 99L), sample(0:3, 1))
    expect_identical(
      as_int_df(lcc2020_remap_table(lcc, cc, replace)),
      as_int_df(remap_inline_preamble_1.0.6(lcc, cc, replace))
    )
  }
})

test_that("lcc2020_remap_table() applies its rules in order", {
  rt <- lcc2020_remap_table(
    lccClasses = c(20L, 210L, 240L),
    ccClasses = c(1L, 14L, 15L, 17L, 18L, NA),
    treeClassesToReplace = 240L
  )
  lookup <- function(lcc, cc) {
    rt[(if (is.na(lcc)) is.na(LCC) else LCC %in% lcc) & (if (is.na(cc)) is.na(CC) else CC %in% cc), newLCC]
  }
  expect_identical(lookup(210L, 1L), 210L) ## forest CC defers to the simulation's LCC
  expect_identical(lookup(210L, 14L), 210L) ## wetland CC does not strip simulated forest
  expect_identical(lookup(210L, 18L), NA_integer_) ## water CC drops the pixel
  expect_identical(lookup(210L, 17L), 99L) ## urban CC -> reclassify
  expect_identical(lookup(NA, 1L), 99L) ## CC forest where LCC has nothing
  expect_identical(lookup(240L, 1L), 99L) ## treeClassesToReplace wins over everything
  ## Rule 1 (LCC 20 -> NA) is overwritten by rule 2 for missing/cropland CC -- kept for v2 fidelity.
  expect_identical(lookup(20L, NA), 20L)
  expect_identical(lookup(20L, 15L), 20L)
})

test_that("lcc2020_remap_table() rule 1 only bites for codes outside the LCC 2020 legend", {
  rt <- lcc2020_remap_table(lccClasses = 30L, ccClasses = c(1L, 99L), treeClassesToReplace = integer())
  expect_identical(rt[LCC %in% 30L & CC %in% 1L, newLCC], 30L)
  expect_identical(rt[LCC %in% 30L & CC %in% 99L, newLCC], NA_integer_)
})

## cc_age_composite ----------------------------------------------------------------------------------

test_that("cc_age_composite() keeps base ages and fills only empty cells, aged to the epoch", {
  base <- cc_grid(ifelse(cc_col <= 5, 100, NA))
  fill <- cc_grid(50)
  out <- terra::values(cc_age_composite(base, cc_poly(), fill, fillYear = 2022, epoch = 2025),
    mat = FALSE
  )
  expect_identical(out[cc_col <= 5], rep(100, 50)) ## never overwritten, never minimum
  expect_identical(out[cc_col > 5], rep(53, 50)) ## 50 at 2022 -> 53 at 2025
})

test_that("cc_age_composite() never reads a NoData sentinel as an age", {
  vals <- ifelse(cc_col <= 5, 100, 65535)
  fillVals <- ifelse(cc_col <= 8, 151, 255) ## 151 = ">150"; 255 = non-treed
  check <- function(base, fill) {
    out <- terra::values(cc_age_composite(base, cc_poly(), fill, fillYear = 2022), mat = FALSE)
    expect_identical(out[cc_col <= 5], rep(100, 50))
    expect_identical(out[cc_col %in% 6:8], rep(154, 30)) ## >150 is aged like any other value
    expect_true(all(is.na(out[cc_col > 8]))) ## not 65535, and not 255 + 3 = 258
  }
  ## in memory
  check(cc_grid(vals), cc_grid(fillVals))
  ## on disk, where NAflag behaves differently -- the result must not depend on storage
  fb <- withr::local_tempfile(fileext = ".tif")
  ff <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(cc_grid(vals), fb, datatype = "INT2U", NAflag = NA)
  terra::writeRaster(cc_grid(fillVals), ff, datatype = "INT1U", NAflag = NA)
  check(terra::rast(fb), terra::rast(ff))
})

test_that("cc_age_composite() leaves cells outside the study area empty even where fill has ages", {
  out <- cc_age_composite(cc_grid(NA_real_), cc_poly(1:4), cc_grid(50), fillYear = 2025)
  v <- terra::values(terra::extend(out, cc_grid()), mat = FALSE)
  expect_identical(v[cc_col <= 4], rep(50, 40))
  expect_true(all(is.na(v[cc_col > 4])))
})

test_that("cc_age_composite() reprojects a fill in another CRS", {
  base <- cc_grid(NA_real_)
  fill <- terra::project(terra::extend(cc_grid(50), 5L, fill = 50), "EPSG:3978", method = "near")
  out <- terra::values(cc_age_composite(base, cc_poly(), fill, fillYear = 2020), mat = FALSE)
  interior <- which(cc_col %in% 3:8 & rep(1:10, each = 10) %in% 3:8)
  expect_identical(out[interior], rep(55, length(interior)))
})

test_that("cc_age_composite() accepts sf study areas", {
  base <- cc_grid(ifelse(cc_col <= 5, 100, NA))
  poly <- cc_poly(2:9)
  expect_identical(
    terra::values(cc_age_composite(base, sf::st_as_sf(poly), cc_grid(7), fillYear = 2025)),
    terra::values(cc_age_composite(base, poly, cc_grid(7), fillYear = 2025))
  )
})

test_that("cc_age_composite() without a fill just crops and clears the sentinel", {
  out <- cc_age_composite(cc_grid(c(65535, rep(10, 99))), cc_poly())
  v <- terra::values(out, mat = FALSE)
  expect_true(is.na(v[1]))
  expect_identical(v[-1], rep(10, 99))
})

test_that("cc_age_composite() rejects an unusable fill year", {
  expect_snapshot(error = TRUE, cc_age_composite(cc_grid(1), cc_poly(), cc_grid(1)))
  expect_snapshot(error = TRUE, cc_age_composite(cc_grid(1), cc_poly(), cc_grid(1), fillYear = 2030))
})

## cc_age_pct_missing --------------------------------------------------------------------------------

test_that("cc_age_pct_missing() ignores cells outside the study area (the bounding-box bug)", {
  ## The polygon covers 30% of the raster, like LandWeb's groups. Every cell inside is aged and
  ## everything outside is NA from masking: the answer is 0, not 70.
  age <- terra::mask(cc_grid(80), cc_poly(1:3))
  expect_equal(as.numeric(cc_age_pct_missing(age, cc_poly(1:3), cc_grid(1L))), 0)
})

test_that("cc_age_pct_missing() ignores non-forest (the whole-study-area bug)", {
  lc <- cc_grid(ifelse(cc_col <= 6, 1L, 18L)) ## forest | water
  age <- cc_grid(ifelse(cc_col <= 6, 80, NA)) ## lakes have no age, by design
  expect_equal(as.numeric(cc_age_pct_missing(age, cc_poly(), lc)), 0)
})

test_that("cc_age_pct_missing() counts forest with no age", {
  lc <- cc_grid(ifelse(cc_col <= 2, 5L, 16L)) ## 20 forest cells, the rest barren
  age <- cc_grid(ifelse(cc_col == 1, 80, ifelse(cc_col == 2 & rep(1:10, each = 10) <= 5, 80, NA)))
  res <- cc_age_pct_missing(age, cc_poly(), lc)
  expect_equal(as.numeric(res), 25)
  expect_identical(attr(res, "nForest"), 20)
  expect_identical(attr(res, "nAged"), 15)
})

test_that("cc_age_pct_missing() counts at the age raster's own resolution (the coarse-grid bug)", {
  ## A checkerboard is 50% missing. Aggregating first would fill every coarse cell from its aged
  ## neighbour and report 0.
  checker <- ifelse((cc_col + rep(1:10, each = 10)) %% 2 == 0, 80, NA)
  expect_equal(as.numeric(cc_age_pct_missing(cc_grid(checker), cc_poly(), cc_grid(2L))), 50)
})

test_that("cc_age_pct_missing() resamples land cover from another grid and CRS", {
  lc <- terra::project(terra::extend(cc_grid(6L), 5L, fill = 6L), "EPSG:3978", method = "near")
  expect_equal(as.numeric(cc_age_pct_missing(cc_grid(80), cc_poly(3:8), lc)), 0)
})

test_that("cc_age_pct_missing() errors when there is no forest to measure", {
  expect_snapshot(error = TRUE, cc_age_pct_missing(cc_grid(80), cc_poly(), cc_grid(18L)))
})
