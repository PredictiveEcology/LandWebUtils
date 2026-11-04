test_that("landmine burns reasonably", {
  set.seed(42L)

  optimPars <- c(-0.255744847773767, -2.48874580638343, -1.62398492055945, -2.25446136444733,
                 2.2983217794681, 4.71955387550406, 0.850655694289403)

  spawnNewActive <- 10^c(optimPars[1], optimPars[2], optimPars[3], optimPars[4])
  sizeCutoffs <- 10^c(optimPars[5], optimPars[6])

  nx <- ny <- 100L

  r <- terra::rast(nrows = ny, ncols = nx,
                   xmin = -nx / 2, xmax = nx / 2,
                   ymin = -ny / 2, ymax = ny / 2)

  rcl <- cbind(1L:5L, c(21L, 6L, 12L, 0L, 30L)) ## make interior poly (4) non-flammable

  ROSmap <- SpaDES.tools::randomPolygons(r, 5) |>
    terra::classify(rcl)

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
    terra::plot(ROSmap)

    burned <- ROSmap
    burned[fires$pixels] <- 50L
    # burned[initialPixels] <- 50L
    terra::plot(burned)
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

  r <- terra::rast(nrows = ny, ncols = nx,
                   xmin = -nx / 2, xmax = nx / 2,
                   ymin = -ny / 2, ymax = ny / 2)

  rcl <- cbind(1L:5L, c(21L, 6L, 12L, 0L, 30L)) ## make interior poly (4) non-flammable

  ROSmap <- SpaDES.tools::randomPolygons(r, 5) |>
    terra::classify(rcl)

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
    terra::plot(ROSmap)

    burned <- ROSmap
    burned[nonflam] <- NA_integer_
    burned[fires$pixels] <- 50L
    # burned[initialPixels] <- 50L
    terra::plot(burned)
  }

  expect_true(any(fires$pixels %in% nonflam)) ## fires can spread into NA

  fa <- attr(fires, "spreadState")$clusterDT
  fa_sum <- fa[, list(numPixelsBurned = sum(size),
                      expectedNumBurned = sum(maxSize),
                      proportionBurned = sum(size) / sum(maxSize))]

  expect_true(fa_sum$numPixelsBurned == length(initialPixels) * maxFireSize)
  expect_true(fa_sum$proportionBurned == 1L)
})

## LandMine masks the pixels a year has already burned by setting their `spreadProb` to NA, and
## leaves `spreadProbRel` alone. The first spread step used probability 1 everywhere, so it could
## burn those pixels again; the later steps could not.
test_that("no step, the first included, spreads into a pixel whose spreadProb is NA", {
  r <- terra::rast(nrows = 21, ncols = 21, xmin = 0, xmax = 21, ymin = 0, ymax = 21)
  terra::values(r) <- 1L
  start <- terra::cellFromRowCol(r, 11, 11)
  ring <- as.vector(terra::adjacent(r, start, directions = 8))
  spreadProb <- rep(0.77, terra::ncell(r))
  spreadProb[ring] <- NA

  ## `maxRetriesPerID = 4L`: no jump, which `spread2()` makes on every 10th retry
  fires <- withr::with_seed(1L, landmine_burn1(
    r, startCells = start, fireSizes = 20, maxRetriesPerID = 4L,
    spreadProbRel = rep(1, terra::ncell(r)), spreadProb = spreadProb
  ))
  expect_false(any(fires$pixels %in% ring))
  expect_identical(fires$pixels, as.integer(start)) ## ringed in: it cannot leave its start cell
})

test_that("an NA spreadProb on a pixel spreadProbRel already masks changes nothing", {
  r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40)
  terra::values(r) <- 1L
  cells <- seq_len(terra::ncell(r))
  nonflam <- cells[terra::colFromCell(r, cells) %in% 18:22 & terra::rowFromCell(r, cells) > 10]
  ros <- rep(1, terra::ncell(r))
  ros[nonflam] <- NA
  spreadProb <- rep(0.77, terra::ncell(r))
  burn <- function(sp) {
    withr::with_seed(3L, landmine_burn1(
      r, startCells = c(205L, 1230L, 1415L), fireSizes = c(150, 60, 300),
      maxRetriesPerID = 9L, spreadProbRel = ros, spreadProb = sp
    ))
  }
  spreadProbNA <- spreadProb
  spreadProbNA[nonflam] <- NA
  expect_identical(burn(spreadProbNA), burn(spreadProb))
  expect_false(any(burn(spreadProbNA)$pixels %in% nonflam))
})
