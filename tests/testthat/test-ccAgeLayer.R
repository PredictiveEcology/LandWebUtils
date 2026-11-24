## The current condition's age layer for the leading-type and large-patch summaries: time since fire
## (TSF, the default since 1.0.1) or stand age (SAM), by argument or by option (#3).

test_that("a current-condition map takes its time-since-fire layer by default", {
  for (f in c("CurrentConditionVTM.tif", "CurrentConditionVTM.grd")) {
    expect_identical(
      LandWebUtils:::.ccAgeFile(file.path("out", f)),
      file.path("out", "CurrentConditionTSF.tif")
    )
  }
})

test_that("ccAgeLayer = 'SAM' takes the stand-age layer", {
  expect_identical(
    LandWebUtils:::.ccAgeFile(file.path("out", "CurrentConditionVTM.tif"), "SAM"),
    file.path("out", "CurrentConditionSAM.tif")
  )
})

test_that("the LandWebUtils.ccAgeLayer option sets the default", {
  withr::local_options(LandWebUtils.ccAgeLayer = "SAM")
  expect_identical(
    LandWebUtils:::.ccAgeFile(file.path("out", "CurrentConditionVTM.tif")),
    file.path("out", "CurrentConditionSAM.tif")
  )
})

test_that("a simulated map has no current-condition age layer", {
  expect_null(LandWebUtils:::.ccAgeFile(file.path("out", "rep01", "vegTypeMap_year0700.tif")))
})

test_that("an unknown age layer is an error", {
  expect_error(LandWebUtils:::.ccAgeFile("CurrentConditionVTM.tif", "AGE"), "should be one of")
})

test_that("LargePatches reads the age layer ccAgeLayer names", {
  withr::local_options(reproducible.useCache = FALSE)
  seen <- NULL
  testthat::local_mocked_bindings(.rasterToMemory = function(x, ...) {
    seen <<- x
    stop("read")
  })
  vtm <- file.path(withr::local_tempdir(), "CurrentConditionVTM.tif")
  expect_error(LargePatches(tsf = "ignored.tif", vtm = vtm, ccAgeLayer = "SAM"), "read")
  expect_identical(basename(seen), "CurrentConditionSAM.tif")
})

test_that("LeadingVegTypeByAgeClass reads the age layer ccAgeLayer names", {
  seen <- NULL
  testthat::local_mocked_bindings(
    rast = function(x, ...) {
      seen <<- x
      stop("read")
    },
    .package = "terra"
  )
  vtm <- file.path(withr::local_tempdir(), "CurrentConditionVTM.tif")
  expect_error(
    LeadingVegTypeByAgeClass(
      tsf = "ignored.tif",
      vtm = vtm,
      ageClassCutOffs = c(0, 40),
      ccAgeLayer = "SAM"
    ),
    "read"
  )
  expect_identical(basename(seen), "CurrentConditionSAM.tif")
})
