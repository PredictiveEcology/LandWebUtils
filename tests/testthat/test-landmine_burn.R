test_that("landmine burns reasonably", {
  set.seed(42L)

  optimPars <- c(-0.255744847773767, -2.48874580638343, -1.62398492055945, -2.25446136444733,
                 2.2983217794681, 4.71955387550406, 0.850655694289403)

  spawnNewActive <- 10^c(optimPars[1], optimPars[2], optimPars[3], optimPars[4])
  sizeCutoffs <- 10^c(optimPars[5], optimPars[6])

  nx <- ny <- 100L

  r <- raster::raster(nrows = ny, ncols = nx,
                      xmn = -nx / 2, xmx = nx / 2,
                      ymn = -ny / 2, ymx = ny / 2)

  rcl <- cbind(1L:5L, c(21L, 6L, 12L, 0L, 30L)) ## make interior poly (4) non-flammable

  ROSmap <- SpaDES.tools::randomPolygons(r, 5) |>
    raster::reclassify(rcl)

  nonflam <- which(ROSmap[] == 0L)

  lthfc <- r
  lthfc[] <- 100L

  spreadProb <- r
  spreadProb[] <- optimPars[7]

  ## no spread on non-flammable pixels
  lthfc[nonflam] <- NA_integer_
  ROSmap[nonflam] <- NA_integer_
  spreadProb[nonflam] <- NA_integer_

  initialPixels <- c(
    (40L * ny) + 65L,
    (55L * ny) + 80L,
    (60L * ny) + 50L
  )

  maxFireSize <- 150 ## size in pixels
  fireSizes <- rep(maxFireSize, length(initialPixels))

  fires <- landmine_burn1(
    landscape = lthfc,
    startCells = initialPixels,
    fireSizes = fireSizes,
    spreadProbRel = ROSmap,
    sizeCutoffs = sizeCutoffs,
    maxRetriesPerID = 4L,
    spawnNewActive = spawnNewActive,
    spreadProb = spreadProb
  )

  if (interactive()) {
    raster::plot(ROSmap)

    burned <- ROSmap
    burned[fires$pixels] <- 50L
    # burned[initialPixels] <- 50L
    raster::plot(burned)
  }

  expect_false(any(fires$pixels %in% nonflam))

  fa <- attr(fires, "spreadState")$clusterDT
  fa_sum <- fa[, list(numPixelsBurned = sum(size),
                      expectedNumBurned = sum(maxSize),
                      proportionBurned = sum(size) / sum(maxSize))]

  expect_true(fa_sum$numPixelsBurned == length(initialPixels) * maxFireSize)
  expect_true(fa_sum$proportionBurned == 1L)
})

test_that("landmine burns reasonably for 'burny' scenarios", {
  set.seed(42L)

  optimPars <- c(-0.255744847773767, -2.48874580638343, -1.62398492055945, -2.25446136444733,
                 2.2983217794681, 4.71955387550406, 0.850655694289403)

  spawnNewActive <- 10^c(optimPars[1], optimPars[2], optimPars[3], optimPars[4])
  sizeCutoffs <- 10^c(optimPars[5], optimPars[6])

  nx <- ny <- 100L

  r <- raster::raster(nrows = ny, ncols = nx,
                      xmn = -nx / 2, xmx = nx / 2,
                      ymn = -ny / 2, ymx = ny / 2)

  rcl <- cbind(1L:5L, c(21L, 6L, 12L, 0L, 30L)) ## make interior poly (4) non-flammable

  ROSmap <- SpaDES.tools::randomPolygons(r, 5) |>
    raster::reclassify(rcl)

  nonflam <- which(ROSmap[] == 0L)

  lthfc <- r
  lthfc[] <- 100L

  spreadProb <- r
  spreadProb[] <- optimPars[7]

  ## this time, allow spread on non-flammable pixels (don't set NAs)
  # lthfc[nonflam] <- NA_integer_
  ROSmap[nonflam] <- 6L
  # spreadProb[nonflam] <- NA_integer_

  initialPixels <- c(
    (40L * ny) + 65L,
    (55L * ny) + 80L,
    (60L * ny) + 50L
  )

  maxFireSize <- 150 ## size in pixels
  fireSizes <- rep(maxFireSize, length(initialPixels))

  fires <- landmine_burn1(
    landscape = lthfc,
    startCells = initialPixels,
    fireSizes = fireSizes,
    spreadProbRel = ROSmap,
    sizeCutoffs = sizeCutoffs,
    maxRetriesPerID = 4L,
    spawnNewActive = spawnNewActive,
    spreadProb = spreadProb,
    omitPixels = nonflam
  )

  if (interactive()) {
    raster::plot(ROSmap)

    burned <- ROSmap
    burned[nonflam] <- NA_integer_
    burned[fires$pixels] <- 50L
    # burned[initialPixels] <- 50L
    raster::plot(burned)
  }

  expect_true(any(fires$pixels %in% nonflam)) ## fires can spread into NA

  fa <- attr(fires, "spreadState")$clusterDT
  fa_sum <- fa[, list(numPixelsBurned = sum(size),
                      expectedNumBurned = sum(maxSize),
                      proportionBurned = sum(size) / sum(maxSize))]

  expect_true(fa_sum$numPixelsBurned == length(initialPixels) * maxFireSize)
  expect_true(fa_sum$proportionBurned == 1L)
})
