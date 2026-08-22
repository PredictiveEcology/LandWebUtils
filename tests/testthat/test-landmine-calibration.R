test_that("landmine_optim_unpack applies the fitted-parameter convention", {
  par <- c(-0.22, -1.55, -0.65, -1.48, 3.21, 4.72, 0.856)
  out <- landmine_optim_unpack(par)

  expect_named(out, c("spawnNewActive", "sizeCutoffs", "spreadProb"))
  expect_equal(out$spawnNewActive, 10^par[1:4])
  expect_equal(out$sizeCutoffs, 10^par[5:6])
  expect_equal(out$spreadProb, par[7])
})

test_that("landmine_optim_unpack accepts a one-row data.frame of par1..par7", {
  par <- c(-0.22, -1.55, -0.65, -1.48, 3.21, 4.72, 0.856)
  df <- stats::setNames(as.data.frame(rbind(par)), paste0("par", 1:7))

  expect_equal(landmine_optim_unpack(df), landmine_optim_unpack(par))
})

test_that("landmine_optim_unpack rejects malformed input", {
  expect_snapshot(error = TRUE, landmine_optim_unpack(c(1, 2, 3)))
  expect_snapshot(error = TRUE, landmine_optim_unpack(c(-0.22, NA, -0.65, -1.48, 3.21, 4.72, 0.856)))
})

test_that("calibration parameters round-trip through the CSV", {
  f <- withr::local_tempfile(fileext = ".csv")
  par1 <- c(-0.22, -1.55, -0.65, -1.48, 3.21, 4.72, 0.856)
  par2 <- c(-0.31, -1.20, -1.10, -2.05, 2.14, 4.64, 0.844)

  landmine_optim_params_append(f, par1, pixelSize = 240, objective = "andison",
                               date = "2026-01-01")
  landmine_optim_params_append(f, par2, pixelSize = 250, objective = "andison",
                               date = "2026-01-02")

  expect_equal(unname(landmine_optim_params_read(f, rowID = 1L)), par1)
  expect_equal(unname(landmine_optim_params_read(f, rowID = 2L)), par2)
  expect_named(landmine_optim_params_read(f, rowID = 1L), paste0("par", 1:7))

  all <- landmine_optim_params_read(f)
  expect_equal(dim(all), c(2L, 7L))

  ## appending must not clobber earlier rows (one of the two former writers did)
  expect_equal(nrow(utils::read.csv(f)), 2L)
  expect_equal(utils::read.csv(f)$pixelSize, c(240, 250))
  ## the objective is recorded, because values are only comparable within one objective
  expect_equal(utils::read.csv(f)$objective, c("andison", "andison"))
})

test_that("landmine_optim_params_read rejects a bad rowID", {
  f <- withr::local_tempfile(fileext = ".csv")
  landmine_optim_params_append(f, c(-0.22, -1.55, -0.65, -1.48, 3.21, 4.72, 0.856),
                               pixelSize = 240, date = "2026-01-01")

  expect_snapshot(error = TRUE, landmine_optim_params_read(f, rowID = 99L))
  expect_snapshot(error = TRUE, landmine_optim_params_read("no-such-file.csv"))
})

test_that("landmine_optim_landscape builds a uniform square with a central ignition cell", {
  lscape <- landmine_optim_landscape(pixelSize = 240, n = 50L)

  expect_equal(terra::ncell(lscape$ros), 2500)
  expect_equal(terra::res(lscape$ros), c(240, 240))
  expect_true(file.exists(lscape$file))

  ## uniform and fully flammable: the calibration is deliberately not landscape-specific.
  ## NOTE: the logical raster round-trips through GeoTIFF as integer 1, not TRUE.
  vals <- unique(as.vector(terra::values(lscape$ros, mat = FALSE)))
  expect_length(vals, 1L)
  expect_equal(as.numeric(vals), 1)

  rc <- terra::rowColFromCell(lscape$ros, lscape$centreCell)
  expect_equal(as.vector(rc), c(25, 25))
})

test_that("landmine_optim_shapeIndex is 1 for a circle and 2/sqrt(pi) for a square", {
  expect_equal(landmine_optim_shapeIndex(2 * pi * 100, pi * 100^2), 1)
  expect_equal(landmine_optim_shapeIndex(4 * 50, 50^2), 2 / sqrt(pi))
})

test_that("landmine_optim_shapeTarget reproduces Andison's fitted curve", {
  ## the thesis text says shape rises to "around five at 1,000 hectares"
  expect_equal(landmine_optim_shapeTarget(1000), 5.091, tolerance = 1e-3)
  expect_equal(landmine_optim_shapeTarget(10), 1.811, tolerance = 1e-3)
  ## defaults must be the EMPIRICAL fit (1.770), not the model's own (1.881)
  expect_equal(landmine_optim_shapeTarget(1), 1.770)
  expect_gt(landmine_optim_shapeTarget(1000), landmine_optim_shapeTarget(100))
})

test_that("landmine_optim_islandTarget follows Andison Table 3.5 (empirical column)", {
  expect_equal(landmine_optim_islandTarget(c(499, 500, 1000, 1001)), c(4.0, 6.0, 6.0, 9.0))
})

test_that("landmine_optim_islands finds enclosed unburned patches only", {
  ## 11x11: burn everything except a 3x3 hole in the middle and a notch on the edge
  r <- terra::rast(nrows = 11, ncols = 11, xmin = 0, xmax = 1100, ymin = 0, ymax = 1100)
  terra::values(r) <- 1L                       # 100 m cells -> 1 ha each
  m <- terra::as.matrix(r, wide = TRUE)
  m[5:7, 5:7] <- NA                            # enclosed island: 9 cells = 9 ha
  m[1, 1:3] <- NA                              # edge notch: NOT an island
  terra::values(r) <- as.vector(t(m))

  isl <- landmine_optim_islands(r, minIslandSizeHa = 2)
  expect_equal(isl$nIslands, 1L)
  expect_equal(isl$islandAreaHa, 9)
  expect_equal(isl$sizesHa, 9)
  ## fire area = burned cells + enclosed island
  expect_equal(isl$fireAreaHa, (121 - 9 - 3) + 9)
  expect_equal(isl$percentIslands, 100 * 9 / ((121 - 9 - 3) + 9))
})

test_that("landmine_optim_islands honours minIslandSizeHa", {
  r <- terra::rast(nrows = 11, ncols = 11, xmin = 0, xmax = 1100, ymin = 0, ymax = 1100)
  terra::values(r) <- 1L
  m <- terra::as.matrix(r, wide = TRUE)
  m[6, 6] <- NA                                # a single 1 ha hole
  terra::values(r) <- as.vector(t(m))

  expect_equal(landmine_optim_islands(r, minIslandSizeHa = 2)$nIslands, 0L)
  expect_equal(landmine_optim_islands(r, minIslandSizeHa = 1)$nIslands, 1L)
})

test_that("landmine_optim_islands returns zeros when nothing is unburned", {
  r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 500, ymin = 0, ymax = 500)
  terra::values(r) <- 1L
  isl <- landmine_optim_islands(r)
  expect_equal(isl$nIslands, 0L)
  expect_equal(isl$percentIslands, 0)
})

test_that("landmine_optim_fireSizes scales with pixel size so areas stay fixed", {
  ha <- c(50, 200, 600, 1500, 3000)
  px100 <- landmine_optim_fireSizes(100, ha)
  px240 <- landmine_optim_fireSizes(240, ha)

  ## same physical fires at both resolutions
  expect_equal(as.integer(px100), as.integer(ha / 1))        # 1 ha per 100 m pixel
  expect_equal(as.integer(px240), as.integer(round(ha / 5.76)))
  expect_equal(attr(px100, "haPerPixel"), 1)
  expect_equal(attr(px240, "haPerPixel"), 5.76)
  ## the pixel counts differ, which is the whole point
  expect_false(identical(as.integer(px100), as.integer(px240)))
})

test_that("landmine_optim_fireSizes warns outside Andison's fitted range", {
  ## this is what the old pixel-count defaults amounted to at 240 m
  expect_snapshot(x <- landmine_optim_fireSizes(240, c(57.6, 576000)))
})

test_that("landmine_optim_fireSizes warns when a fire is only a few pixels", {
  expect_snapshot(x <- landmine_optim_fireSizes(240, c(11, 600)))
})
