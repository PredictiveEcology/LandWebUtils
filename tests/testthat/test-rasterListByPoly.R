## one 4 x 4 raster per replicate folder, and two polygons covering parts of it
local_rep_rasters <- function(fname, env = parent.frame()) {
  tmp <- withr::local_tempdir(.local_envir = env)
  r <- terra::rast(
    nrows = 4,
    ncols = 4,
    xmin = 0,
    xmax = 4,
    ymin = 0,
    ymax = 4,
    crs = "EPSG:3857",
    vals = 1:16
  )
  files <- file.path(tmp, c("rep01", "rep10"), fname)
  fs::dir_create(dirname(files))
  for (f in files) {
    terra::writeRaster(r, f)
  }
  files
}

polyNames <- c("Unit_A.1", "Other")

test_polys <- function() {
  sf::st_sf(
    name = polyNames,
    geometry = sf::st_sfc(square(0, 0, 2), square(2, 2, 2), crs = "EPSG:3857")
  )
}

test_that("rasterListByPoly() keeps '.' and '_' in polygon names", {
  files <- local_rep_rasters("rstTimeSinceFire_year0100.tif")

  rl <- rasterListByPoly(files, test_polys(), polyNames, "name", "rstTimeSinceFire_")

  expect_length(rl, 4L)
  expect_s4_class(rl[[1]], "SpatRaster")
  expect_identical(
    names(rl),
    paste0(rep(c("rep01", "rep10"), each = 2), ".rstTimeSinceFire_year0100_", polyNames)
  )
  expect_identical(attr(rl, "reps"), c(1L, 1L, 10L, 10L))
  expect_identical(attr(rl, "times"), rep(100L, 4))
  expect_identical(attr(rl, "polyNames"), rep(polyNames, 2))
})

test_that("rasterListByPoly() stops on a file name without a year", {
  files <- local_rep_rasters("rstTimeSinceFire.tif")

  expect_error(
    rasterListByPoly(files, test_polys(), polyNames, "name", "rstTimeSinceFire_"),
    "cannot parse raster label"
  )
})
