test_that("reportingPolygonLayers returns the expected schema", {
  lyrs <- reportingPolygonLayers()
  expect_named(lyrs, c("key", "source", "id", "labelCol"))
  expect_setequal(unique(lyrs$source), c("drive", "url"))
  expect_contains(lyrs$key, "FMA Boundaries Updated")
})

test_that("landbaseLayers returns the expected schema", {
  lyrs <- landbaseLayers()
  expect_named(lyrs, c("key", "source", "id", "layer", "status_col", "applies_to"))
  expect_contains(lyrs$key, "Spray Lake C5 Landbase")
  expect_false(any(is.na(lyrs$status_col)))
})

test_that(".findLandbaseSource prefers a File Geodatabase, then a shapefile", {
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "coverage.gdb"))
  expect_match(LandWebUtils:::.findLandbaseSource(dir), "coverage\\.gdb$")

  dir2 <- withr::local_tempdir()
  file.create(file.path(dir2, "coverage.shp"))
  expect_match(LandWebUtils:::.findLandbaseSource(dir2), "coverage\\.shp$")

  expect_null(LandWebUtils:::.findLandbaseSource(withr::local_tempdir()))
})

test_that("buildLandbasePolygons skips (no download) when no landbase applies", {
  sa <- terra::vect("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
  terra::crs(sa) <- "EPSG:3857"

  out <- buildLandbasePolygons(sa, studyAreaName = "NoSuchArea", destinationPath = withr::local_tempdir())
  expect_identical(out, list())
})

test_that("buildLandbasePolygons gates to the applicable landbase by study-area name", {
  sa <- terra::vect("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
  terra::crs(sa) <- "EPSG:3857"

  ## intercept the fetch so the test needs no network: record which layers pass the gate
  fetched <- character()
  testthat::local_mocked_bindings(
    drive_download_once = function(id, path, ...) {
      fetched <<- c(fetched, path)
      stop("stop before extraction")
    },
    .package = "workflowtools"
  )
  suppressWarnings(try(
    buildLandbasePolygons(sa, studyAreaName = "SprayLake", destinationPath = withr::local_tempdir()),
    silent = TRUE
  ))
  expect_length(fetched, 1L)
  expect_match(fetched, "Spray.Lake.C5")
})
