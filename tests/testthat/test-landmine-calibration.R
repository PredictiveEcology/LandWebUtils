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

  landmine_optim_params_append(f, par1, pixelSize = 240, date = "2026-01-01")
  landmine_optim_params_append(f, par2, pixelSize = 250, date = "2026-01-02")

  expect_equal(unname(landmine_optim_params_read(f, rowID = 1L)), par1)
  expect_equal(unname(landmine_optim_params_read(f, rowID = 2L)), par2)
  expect_named(landmine_optim_params_read(f, rowID = 1L), paste0("par", 1:7))

  all <- landmine_optim_params_read(f)
  expect_equal(dim(all), c(2L, 7L))

  ## appending must not clobber earlier rows (one of the two former writers did)
  expect_equal(nrow(utils::read.csv(f)), 2L)
  expect_equal(utils::read.csv(f)$pixelSize, c(240, 250))
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
