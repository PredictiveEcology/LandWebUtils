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

test_that(".slug uses underscores, not the dots make.names() would produce", {
  expect_identical(LandWebUtils:::.slug("National Ecozones"), "National_Ecozones")
  expect_identical(LandWebUtils:::.slug("BC Biogeoclimatic zones"), "BC_Biogeoclimatic_zones")
  expect_identical(LandWebUtils:::.slug("Parks"), "Parks")
  expect_no_match(LandWebUtils:::.slug(reportingPolygonLayers()$key), "[.]", all = TRUE)
})

test_that(".slug collapses runs of non-alphanumerics and trims them", {
  expect_identical(LandWebUtils:::.slug("a  b--c"), "a_b_c")
  expect_identical(LandWebUtils:::.slug(" leading/trailing "), "leading_trailing")
})

test_that(".findVectorFile searches recursively (archives that unpack to a subdir)", {
  dir <- withr::local_tempdir()
  nested <- file.path(dir, "Ecoregions")
  dir.create(nested)
  file.create(file.path(nested, "ecoregions.shp"))

  ## a non-recursive list.files() found nothing here, so the layer was dropped silently
  expect_match(LandWebUtils:::.findVectorFile(dir), "ecoregions\\.shp$")
  expect_length(LandWebUtils:::.findVectorFile(withr::local_tempdir()), 0L)
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
  expect_match(fetched, "Spray_Lake_C5", fixed = TRUE)
})
