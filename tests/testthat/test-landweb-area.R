test_that("the LandWeb area is the outline of the fire-cycle polygons, without holes", {
  holed <- sf::st_difference(square(0, 0, 10), square(4, 4, 2))
  lthfc <- sf::st_sf(
    LTHFC = c(100, 60),
    geometry = sf::st_sfc(holed, square(10, 0, 10), crs = 3400)
  )
  a <- landweb_area(lthfc)
  expect_s3_class(a, "sfc")
  expect_length(a, 1L)
  expect_equal(as.numeric(sf::st_area(a)), 200)
})

test_that("the fitting area keeps whole the ecological units that touch the area", {
  ## shapefile field names stop at 10 characters: the national layer's is ECOPROVINC
  eco <- sf::st_sf(
    ECOPROVINC = c("A", "B", "C"),
    geometry = sf::st_sfc(square(0, 0, 10), square(10, 0, 10), square(30, 0, 10), crs = 3400)
  )
  dest <- withr::local_tempdir()
  dir.create(file.path(dest, "Ecoprovinces"))
  sf::st_write(eco, file.path(dest, "Ecoprovinces", "ecoprovinces.shp"), quiet = TRUE)
  testthat::local_mocked_bindings(
    download_once = function(...) invisible(NULL),
    archive_extract_once = function(...) invisible(NULL),
    .package = "workflowtools"
  )
  area <- sf::st_transform(sf::st_sfc(square(5, 2, 2), crs = 3400), 3978)
  out <- landweb_anpp_area(area, "ecoprovince", destinationPath = dest)
  expect_identical(out$ECOPROVINC, "A")
  expect_equal(sf::st_crs(out), sf::st_crs(area))
  expect_equal(
    as.numeric(sf::st_area(out)),
    as.numeric(sf::st_area(sf::st_transform(eco[1, ], 3978)))
  )
})

test_that("landweb_anpp_area stops on a shapefile whose .dbf is cut short", {
  testthat::local_mocked_bindings(
    download_once = function(...) invisible(NULL),
    archive_extract_once = function(...) invisible(NULL),
    .package = "workflowtools"
  )
  dest <- withr::local_tempdir()
  local_cut_shapefile(
    keep = 2L,
    cols = list(ECOPROVINC = c("A", "B", "C")),
    dir = file.path(dest, "Ecoprovinces"),
    name = "ecoprovinces"
  )
  area <- sf::st_sfc(square(0, 0, 1), crs = 3857)
  expect_error(
    landweb_anpp_area(area, "ecoprovince", destinationPath = dest),
    "2 complete records.*3 shapes"
  )
})

test_that("landweb_lthfc stops on a shapefile whose .dbf is cut short", {
  testthat::local_mocked_bindings(
    drive_download_once = function(...) invisible(NULL),
    archive_extract_once = function(...) invisible(NULL),
    .package = "workflowtools"
  )
  dest <- withr::local_tempdir()
  local_cut_shapefile(
    keep = 2L,
    cols = list(LTFC10 = c(60, 100, 170)),
    dir = file.path(dest, "lthfc")
  )
  expect_error(landweb_lthfc(dest), "2 complete records.*3 shapes")
})
