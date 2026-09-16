## LandMine fire-spread calibration: the pieces that were previously duplicated between
## `LandMine.Rmd` (which is knitted with `eval = FALSE`, so it never ran and silently
## rotted) and `LandMine.R`. Keeping them here makes them testable and gives the
## parameter convention and the CSV schema a single definition.

#' Unpack LandMine fire-spread calibration parameters
#'
#' The DEoptim calibration fits seven parameters on a transformed scale: the four
#' `spawnNewActive` probabilities and the two `sizeCutoffs` are fitted as base-10
#' logarithms, while `spreadProb` is fitted directly. This function is the single
#' definition of that convention -- it was previously repeated at eight call sites.
#'
#' @param par Numeric vector of length 7, or anything coercible to one (e.g. a
#'            one-row `data.frame` of `par1`..`par7`), ordered
#'            `log10(spawnNewActive[1:4])`, `log10(sizeCutoffs[1:2])`, `spreadProb`.
#'
#' @return A named list with elements `spawnNewActive` (length 4), `sizeCutoffs`
#'         (length 2), and `spreadProb` (length 1), all on the natural scale.
#'
#' @export
#' @rdname landmine_optim_params
#'
#' @examples
#' landmine_optim_unpack(c(-0.22, -1.55, -0.65, -1.48, 3.21, 4.72, 0.856))
landmine_optim_unpack <- function(par) {
  par <- as.numeric(unlist(unname(par)))
  if (length(par) != 7L) {
    stop("`par` must have 7 elements (4 spawnNewActive, 2 sizeCutoffs, 1 spreadProb), not ",
         length(par), ".", call. = FALSE)
  }
  if (anyNA(par)) {
    stop("`par` must not contain NA.", call. = FALSE)
  }
  list(
    spawnNewActive = 10^par[1:4],
    sizeCutoffs = 10^par[5:6],
    spreadProb = par[7]
  )
}

#' Read and append LandMine calibration parameters
#'
#' The calibration results live in a CSV (`LandMine_DEoptim_params.csv`) whose schema is
#' `date`, `pixelSize`, then `par1`..`par7`. These accessors are the single definition of
#' that schema; previously it was written in one place and read in two others, and the two
#' writers disagreed (one appended, one overwrote).
#'
#' @param file Path to the calibration CSV.
#' @param rowID Optional row number to extract. If `NULL` (default), all rows are returned.
#'
#' @return `landmine_optim_params_read()` returns a named numeric vector of length 7 when
#'         `rowID` is given, otherwise a matrix with one row per calibration run.
#'
#' @export
#' @rdname landmine_optim_params
landmine_optim_params_read <- function(file, rowID = NULL) {
  if (!file.exists(file)) {
    stop("calibration parameter file not found: ", file, call. = FALSE)
  }
  d <- utils::read.csv(file, stringsAsFactors = FALSE)
  parCols <- grepl("^par[0-9]+$", colnames(d))
  if (sum(parCols) != 7L) {
    stop("expected 7 `par*` columns in ", file, ", found ", sum(parCols), ".", call. = FALSE)
  }
  if (is.null(rowID)) {
    return(as.matrix(d[, parCols, drop = FALSE]))
  }
  if (length(rowID) != 1L || is.na(rowID) || rowID < 1L || rowID > nrow(d)) {
    stop("`rowID` must be a single row number between 1 and ", nrow(d), ".", call. = FALSE)
  }
  stats::setNames(
    as.numeric(unlist(unname(d[rowID, parCols]))),
    colnames(d)[parCols]
  )
}

#' @param par Numeric vector of length 7 to append (see `landmine_optim_unpack()`).
#' @param pixelSize Pixel resolution (m) the calibration was run at.
#' @param objective Name of the objective the fit was made against. Recorded because parameter
#'   values are only comparable within one objective: fits made against
#'   [landmine_optim_fitSN()] before its size penalty was repaired, or before fire sizes were
#'   specified in hectares, are not comparable to anything fitted since.
#' @param date Date of the run, as a string. Defaults to today.
#'
#' @return `landmine_optim_params_append()` invisibly returns the full table as written.
#'
#' @export
#' @rdname landmine_optim_params
landmine_optim_params_append <- function(file, par, pixelSize, objective = NA_character_,
                                         date = format(Sys.time(), "%Y-%m-%d")) {
  par <- as.numeric(unlist(unname(par)))
  invisible(landmine_optim_unpack(par)) ## validate via the single convention
  newRow <- data.frame(date = date, pixelSize = pixelSize, objective = objective,
                       stringsAsFactors = FALSE)
  newRow <- cbind(newRow, stats::setNames(as.data.frame(rbind(par)), paste0("par", 1:7)))
  out <- if (file.exists(file)) {
    prev <- utils::read.csv(file, stringsAsFactors = FALSE)
    if (!identical(colnames(prev), colnames(newRow))) {
      stop("column mismatch appending to ", file, ":\n  existing: ",
           paste(colnames(prev), collapse = ", "), "\n  new: ",
           paste(colnames(newRow), collapse = ", "), call. = FALSE)
    }
    rbind(prev, newRow)
  } else {
    newRow
  }
  utils::write.csv(out, file, row.names = FALSE)
  invisible(out)
}

#' Build the synthetic landscape used for LandMine calibration
#'
#' The calibration is deliberately run on a uniform, fully-flammable square rather than a
#' real study area, so that the fitted shape statistics are not confounded by a particular
#' landscape's fuel pattern. It is run at several pixel sizes to confirm the parameters
#' transfer across resolutions.
#'
#' @param pixelSize Pixel resolution, in map units (m).
#' @param n Number of pixels per side. Default 1000, i.e. a 1e6-cell landscape.
#' @param file Optional path to write the raster to. The objective function takes the
#'             landscape as a *file path*, because `SpatRaster` objects cannot be
#'             serialized to cluster workers.
#'
#' @return A list with `file` (path), `ros` (the `SpatRaster`), and `centreCell` (the cell
#'         index at the centre of the landscape, used as the ignition point).
#'
#' @export
#'
#' @examples
#' \donttest{
#' lscape <- landmine_optim_landscape(pixelSize = 240, n = 100)
#' terra::ncell(lscape$ros)
#' }
landmine_optim_landscape <- function(pixelSize, n = 1000L,
                                     file = tempfile(fileext = ".tif")) {
  stopifnot(length(pixelSize) == 1L, pixelSize > 0, length(n) == 1L, n > 0)
  ros <- terra::rast(
    terra::ext(0, pixelSize * n, 0, pixelSize * n),
    resolution = pixelSize,
    vals = 0
  )
  ros <- ros == 0 ## uniform, fully flammable
  ros <- terra::writeRaster(ros, file, overwrite = TRUE)
  list(
    file = file,
    ros = ros,
    centreCell = terra::cellFromRowCol(
      ros,
      row = terra::nrow(ros) / 2,
      col = terra::ncol(ros) / 2
    )
  )
}

#' Run the LandMine fire-spread calibration
#'
#' Wraps the DEoptim fit that was previously written inline in `LandMine.Rmd`. Defaults
#' reproduce that fit exactly.
#'
#' @note The objective is stochastic, so treat any single fit as one draw rather than a point
#'   estimate, and pass `crnSeed`. The original objective ([landmine_optim_fitSN()]) did not
#'   identify the parameters at all: in 43 repeat fits recorded in
#'   `LandMine_DEoptim_params.csv`, four of seven scattered over more than 90% of their search
#'   interval, and repeat solutions proved statistically indistinguishable from one another --
#'   and from some randomly drawn parameter vectors. `objective = "andison"` adds constraints to
#'   address that; verify identifiability again after any change to the objective rather than
#'   assuming it.
#'
#' @param objective Which objective to optimise. `"andison"` (default) uses
#'   [landmine_optim_fitAndison()], which fits the shape-vs-area relationship and the remnant
#'   island fractions Andison (1996) reported. `"sn"` uses the original
#'   [landmine_optim_fitSN()], retained so historical fits can be reproduced.
#' @param replicates Passed to [landmine_optim_fitAndison()]: evaluations to average per fire
#'   size. The island statistic is highly variable between draws, so `replicates > 1` buys
#'   discrimination at linear cost. `objective = "andison"` only.
#' @param crnSeed Passed to [landmine_optim_fitAndison()] as common random numbers. Strongly
#'   recommended: the objective is stochastic, and without it the optimiser compares parameter
#'   vectors evaluated against *different* random draws.
#' @param pixelSize,n,file Passed to `landmine_optim_landscape()`.
#' @param fireSizesHa Target fire sizes in **hectares**, converted to pixels for the current
#'   `pixelSize` by [landmine_optim_fireSizes()]. Defaults sit inside the 10--3,000 ha range
#'   Andison fitted.
#' @param fireSizes Target fire sizes in **pixels**. Normally `NULL`; supply it only to
#'   reproduce a historical fit that was specified in pixels, and note that a pixel-count
#'   vector means a different physical fire at each resolution.
#' @param desiredPerimeterArea Target perimeter-to-area ratio (m/m^2); `objective = "sn"` only.
#' @param lower,upper Search bounds for the 7 parameters.
#' @param NP DEoptim population size. Defaults to the number of cluster nodes.
#' @param itermax,VTR,strategy Passed to `DEoptim::DEoptim.control()`.
#' @param cl An existing cluster, or `NULL` to build one from `nodes`.
#' @param nodes Node specification for `landmine_optim_clusterSetup()`.
#' @param seed Optional integer seed, set before cluster creation so the per-worker RNG
#'             streams are reproducible. Note that reproducibility also requires `NP` to be
#'             held fixed, since the stream count follows the number of workers.
#' @param paramsFile Optional path to a calibration CSV; if given, the best parameters
#'                   found are appended to it.
#'
#' @return The `DEoptim` result object.
#'
#' @export
landmine_optim_calibrate <- function(pixelSize = 240, n = 1000L,
                                     file = tempfile(fileext = ".tif"),
                                     fireSizesHa = c(50, 200, 600, 1500, 3000),
                                     fireSizes = NULL,
                                     desiredPerimeterArea = 0.003,
                                     lower = c(-2, -3, -3, -3, 1, 3.5, 0.75),
                                     upper = c(-0.1, -0.5, -0.5, -1, 3.5, 5, 1),
                                     NP = NULL, itermax = 200L, VTR = 0.001, strategy = 6L,
                                     cl = NULL, nodes = NULL, seed = NULL,
                                     objective = c("andison", "sn"), crnSeed = 1L,
                                     replicates = 1L,
                                     paramsFile = NULL) {
  objective <- match.arg(objective)
  if (!requireNamespace("DEoptim", quietly = TRUE)) {
    stop("`DEoptim` is required to run the calibration; install it and retry.", call. = FALSE)
  }
  stopifnot(length(lower) == 7L, length(upper) == 7L, all(upper > lower))

  if (!is.null(seed)) {
    set.seed(seed)
  }
  ownCluster <- is.null(cl)
  if (ownCluster && !is.null(nodes)) {
    cl <- landmine_optim_clusterSetup(nodes = nodes, seed = seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  ## NP is the DEoptim population size (>= 10 per parameter is the usual heuristic; 7
  ## parameters -> 70). It must NOT track the worker count: doing so made both the population
  ## size and the per-worker RNG streams depend on the machine, so a run was not reproducible
  ## across hosts or even across differing load on one host.
  if (is.null(NP)) {
    NP <- 70L
  }

  lscape <- landmine_optim_landscape(pixelSize = pixelSize, n = n, file = file)

  ## Fire sizes are specified in HECTARES and converted here, so the same call means the same
  ## physical fires at any resolution. `fireSizes` (pixels) is honoured when given, to reproduce
  ## historical fits that were specified that way.
  if (is.null(fireSizes)) {
    fireSizes <- landmine_optim_fireSizes(pixelSize = pixelSize, areaHa = fireSizesHa)
  }

  ctrl <- DEoptim::DEoptim.control(
    VTR = VTR, NP = as.integer(NP), itermax = itermax, cluster = cl, strategy = strategy
  )
  res <- if (identical(objective, "andison")) {
    DEoptim::DEoptim(
      fn = landmine_optim_fitAndison, lower = lower, upper = upper, control = ctrl,
      ros = lscape$file, centreCell = lscape$centreCell,
      fireSizes = fireSizes, crnSeed = crnSeed, replicates = replicates
    )
  } else {
    DEoptim::DEoptim(
      fn = landmine_optim_fitSN, lower = lower, upper = upper, control = ctrl,
      ros = lscape$file, centreCell = lscape$centreCell,
      fireSizes = fireSizes, desiredPerimeterArea = desiredPerimeterArea
    )
  }

  if (!is.null(paramsFile)) {
    landmine_optim_params_append(
      file = paramsFile,
      par = res$optim$bestmem,
      pixelSize = pixelSize,
      objective = objective
    )
  }
  res
}

#' Andison's shape index
#'
#' The dimensionless shape index used by Andison (1996, §3.3.1): the patch perimeter relative
#' to that of a circle of the same area, so a circle scores 1 and more convoluted shapes score
#' higher.
#'
#' @param perimeter Patch perimeter, in metres.
#' @param area Patch area, in square metres.
#'
#' @return Numeric shape index.
#'
#' @note Perimeters measured on a raster follow cell edges, so they exceed the perimeter of the
#'   equivalent vector polygon and inflate this index relative to the GIS-derived values Andison
#'   fitted. The bias is systematic, so it shifts the intercept rather than the slope.
#'
#' @export
#' @examples
#' ## a circle scores 1
#' landmine_optim_shapeIndex(perimeter = 2 * pi * 100, area = pi * 100^2)
landmine_optim_shapeIndex <- function(perimeter, area) {
  perimeter / (2 * sqrt(pi * area))
}

#' Andison's expected shape index for a given fire area
#'
#' `SHAPE = intercept + slope * log10(AREA)^4`, fitted by Andison (1996) over fires of
#' 10--3,000 ha. The default intercept is the one fitted to the **empirical** data
#' (1.770, R^2 = 0.670, n = 928). The value 1.881 that also appears in the thesis is the fit to
#' the *model's own* output, and targeting it would be circular.
#'
#' @param areaHa Fire area, in hectares.
#' @param intercept,slope Coefficients; defaults are Andison's empirical fit.
#'
#' @return Expected shape index.
#'
#' @export
#' @examples
#' landmine_optim_shapeTarget(c(10, 100, 1000)) ## ~1.8, ~2.1, ~5.1
landmine_optim_shapeTarget <- function(areaHa, intercept = 1.770, slope = 0.041) {
  intercept + slope * log10(areaHa)^4
}

#' Remnant islands within a simulated fire
#'
#' Unburned patches fully enclosed by the fire. Their formation is an emergent by-product of the
#' spread algorithm rather than something the model targets directly, which is why Andison (1996,
#' §3.3.2) used them as a validation statistic.
#'
#' @param burnedMap A `SpatRaster` in which burned cells are non-`NA` (as returned by
#'                  `landmine_optim_burnFun()`).
#' @param minIslandSizeHa Islands smaller than this are ignored, matching Andison, who dropped
#'                        islands below 2 ha because the model could not represent them and the
#'                        inventory data did not resolve them.
#'
#' @return A list with `nIslands`, `islandAreaHa` (total), `fireAreaHa` (burned plus enclosed
#'         islands, i.e. the area within the fire perimeter), `percentIslands`, and `sizesHa`.
#'
#' @details Unburned patches are found with 4-connectivity, the topological complement of the
#'   8-connected spread used to grow the fire. A patch touching the raster edge is exterior, not
#'   an island; if the fire reaches the edge, no enclosure can be determined for that side.
#'
#' @export
landmine_optim_islands <- function(burnedMap, minIslandSizeHa = 2) {
  cellAreaHa <- prod(terra::res(burnedMap)) / 1e4
  burned <- !is.na(burnedMap)
  nBurned <- sum(terra::values(burned, mat = FALSE), na.rm = TRUE)

  unburned <- terra::ifel(burned, NA, 1L)
  ## guard BEFORE patches(): it errors on an all-NA raster, which is exactly what a fire that
  ## consumed the whole landscape produces.
  noIslands <- list(nIslands = 0L, islandAreaHa = 0, fireAreaHa = nBurned * cellAreaHa,
                    percentIslands = 0, sizesHa = numeric(0))
  if (!any(!is.na(terra::values(unburned, mat = FALSE)))) {
    return(noIslands)
  }
  pat <- terra::patches(unburned, directions = 4, zeroAsNA = FALSE, allowGaps = FALSE)
  pv <- terra::values(pat, mat = FALSE)
  if (all(is.na(pv))) {
    return(noIslands)
  }

  ## any patch reaching the raster edge is outside the fire, not enclosed by it
  m <- terra::as.matrix(pat, wide = TRUE)
  edgeIDs <- unique(stats::na.omit(c(m[1, ], m[nrow(m), ], m[, 1], m[, ncol(m)])))

  tab <- table(pv[!is.na(pv)])
  ids <- as.integer(names(tab))
  sizes <- as.integer(tab) * cellAreaHa
  keep <- !(ids %in% edgeIDs) & sizes >= minIslandSizeHa
  sizesHa <- unname(sizes[keep])

  islandAreaHa <- sum(sizesHa)
  fireAreaHa <- nBurned * cellAreaHa + islandAreaHa
  list(
    nIslands = length(sizesHa),
    islandAreaHa = islandAreaHa,
    fireAreaHa = fireAreaHa,
    percentIslands = if (fireAreaHa > 0) 100 * islandAreaHa / fireAreaHa else 0,
    sizesHa = sizesHa
  )
}

#' Andison island-remnant target
#'
#' Percent of total fire area expected to be in remnant islands, by disturbance size class
#' (Andison 1996, Table 3.5, **empirical** column from DeLong and Tanner 1995). The model column
#' of that table (3.3 / 6.2 / 10.5) is the previous model's output, not a target.
#'
#' @param fireAreaHa Fire area, in hectares.
#'
#' @return Expected percent of fire area in remnant islands.
#'
#' @export
landmine_optim_islandTarget <- function(fireAreaHa) {
  ifelse(fireAreaHa < 500, 4.0, ifelse(fireAreaHa <= 1000, 6.0, 9.0))
}

#' Andison-constrained objective for the LandMine fire-spread calibration
#'
#' A replacement for [landmine_optim_fitSN()] that fits against the validation statistics
#' Andison (1996) actually reported, rather than a single constant perimeter-to-area target.
#'
#' @details
#' [landmine_optim_fitSN()] penalises departure from one constant `desiredPerimeterArea` at
#' *every* fire size. Andison's central result, though, is that shape complexity **increases
#' with fire area** -- that is what the `log10(AREA)^4` term describes -- so the old objective
#' works against the relationship it exists to reproduce. Combined with a size penalty that
#' never fired, it reduced to a single statistic constraining seven parameters, which left the
#' fitted values unidentifiable: across 43 repeat fits, four of seven scattered over more than
#' 90% of their search range, and repeat solutions were statistically indistinguishable from one
#' another (and from some randomly drawn parameter vectors).
#'
#' This objective sums three squared relative deviations per fire size:
#'
#' \itemize{
#'   \item **shape** -- observed vs [landmine_optim_shapeTarget()];
#'   \item **size** -- relative shortfall against the target size, continuous rather than a step,
#'         so the optimiser gets a direction instead of a cliff;
#'   \item **islands** -- percent of fire area in remnant islands vs
#'         [landmine_optim_islandTarget()].
#' }
#'
#' @param par Numeric vector of length 7; see [landmine_optim_unpack()].
#' @param ros Path to the landscape raster (a path, not a `SpatRaster`: these are not
#'            serializable to cluster workers).
#' @param centreCell Ignition cell.
#' @param fireSizes Target fire sizes, in pixels.
#' @param shapeIntercept,shapeSlope Coefficients for [landmine_optim_shapeTarget()].
#' @param minIslandSizeHa Passed to [landmine_optim_islands()].
#' @param islandRangeHa Length-2 range of fire areas over which the island term is applied.
#'   Outside it the term is dropped, because the target is not supported there: Andison
#'   "did not bother simulating many fires less than 100 hectares since disturbances below this
#'   size are relatively constant and simple in shape, and generally do not produce remnant
#'   islands", and the empirical island data (DeLong and Tanner 1995) covers 159--2,239 ha. A
#'   small fire cannot satisfy the target at all -- a 9-pixel fire has no room for an interior
#'   hole, so it scores the maximum penalty regardless of the parameters, which is noise rather
#'   than signal.
#' @param replicates Number of independent evaluations to average per fire size. The island
#'   statistic in particular is highly variable between draws, so `replicates > 1` buys
#'   discrimination at linear cost.
#' @param wShape,wSize,wIsland Weights on the three components. Set one to 0 to drop it.
#'   `wIsland` defaults to [landmine_optim_islandWeight()] for the landscape's own resolution,
#'   because a coarse grid cannot form the small islands that carry most of Andison's island
#'   area. Pass an explicit value to override.
#' @param crnSeed Optional integer. If given, the RNG is seeded as `crnSeed + i` before the
#'                `i`th fire, so every parameter vector is evaluated against the **same** random
#'                draws (common random numbers). This is variance reduction, not determinism:
#'                it makes comparisons between parameter vectors far less noisy.
#'
#' @return A single numeric objective value, with the per-fire components attached as the
#'         `"components"` attribute for diagnostics.
#'
#' @seealso [landmine_optim_fitSN()] for the original objective, retained so historical fits can
#'   be reproduced.
#'
#' @export
landmine_optim_fitAndison <- function(par, ros, centreCell,
                                      fireSizes = 10^(2:5),
                                      shapeIntercept = 1.770, shapeSlope = 0.041,
                                      minIslandSizeHa = 2, islandRangeHa = c(100, Inf),
                                      replicates = 1L,
                                      wShape = 1, wSize = 1, wIsland = NULL,
                                      crnSeed = NULL) {
  stopifnot(is.character(ros))
  ros <- terra::rast(ros)
  if (is.null(wIsland)) {
    ## scale by resolution: a coarse grid cannot form the small islands that carry most of
    ## Andison's island area, so the target is partly unreachable by construction.
    wIsland <- landmine_optim_islandWeight(terra::res(ros)[1])
  }
  p <- landmine_optim_unpack(par)

  if (!is.null(crnSeed)) {
    ## leave the caller's generator as we found it
    oldRNG <- RNGkind()
    on.exit(RNGkind(kind = oldRNG[1], normal.kind = oldRNG[2], sample.kind = oldRNG[3]),
            add = TRUE)
  }

  grid <- expand.grid(i = seq_along(fireSizes), rep = seq_len(replicates))
  comp <- lapply(seq_len(nrow(grid)), function(k) {
    i <- grid$i[k]
    if (!is.null(crnSeed)) {
      ## Distinct stream per (fire size, replicate), identical across parameter vectors.
      ## The `kind` matters: `set.seed()` alone reseeds whatever generator is current, and
      ## cluster workers run L'Ecuyer-CMRG (from `clusterSetRNGStream()`) while a plain session
      ## runs Mersenne-Twister. Without pinning it, the same parameters and the same `crnSeed`
      ## gave objective values differing by ~8x depending on where they were evaluated, so a
      ## calibration result could not be reproduced outside the cluster that produced it.
      set.seed(crnSeed + i + 1000L * grid$rep[k], kind = "Mersenne-Twister")
    }
    bf <- landmine_optim_burnFun(
      ros = ros, centreCell = centreCell, fireSize = fireSizes[i],
      spawnNewActive = p$spawnNewActive, sizeCutoffs = p$sizeCutoffs,
      spreadProb = p$spreadProb
    )
    isl <- landmine_optim_islands(bf$burnedMap, minIslandSizeHa = minIslandSizeHa)

    areaHa <- bf$LM[1, "area"] / 1e4
    shapeObs <- bf$LM[1, "shapeIndex"]
    shapeExp <- landmine_optim_shapeTarget(areaHa, shapeIntercept, shapeSlope)
    islExp <- landmine_optim_islandTarget(isl$fireAreaHa)

    ## the island target is unsupported outside `islandRangeHa`; drop rather than penalise
    islandUsed <- isl$fireAreaHa >= islandRangeHa[1] && isl$fireAreaHa <= islandRangeHa[2]

    data.frame(
      rep = grid$rep[k],
      fireSize = fireSizes[i],
      nBurned = bf$LM[1, "nBurned"],
      areaHa = areaHa,
      shapeObs = shapeObs,
      shapeExp = shapeExp,
      pctIslandsObs = isl$percentIslands,
      pctIslandsExp = islExp,
      shape = ((shapeObs - shapeExp) / shapeExp)^2,
      ## continuous shortfall in [0, 1]; 0 once the fire reaches target
      size = max(0, 1 - bf$LM[1, "nBurned"] / fireSizes[i])^2,
      islandUsed = islandUsed,
      island = if (islandUsed) ((isl$percentIslands - islExp) / islExp)^2 else 0
    )
  })
  comp <- do.call(rbind, comp)

  ## average across replicates, then sum across fire sizes
  perFire <- stats::aggregate(cbind(shape, size, island) ~ fireSize, data = comp, FUN = mean)
  out <- sum(wShape * perFire$shape + wSize * perFire$size + wIsland * perFire$island)
  attr(out, "components") <- comp
  attr(out, "weights") <- c(shape = wShape, size = wSize, island = wIsland)
  out
}

#' Fire sizes for calibration, in pixels, from areas in hectares
#'
#' The calibration targets are defined in **hectares** (Andison fitted the shape relationship
#' over 10--3,000 ha), but the fire model takes sizes in **pixels**. Specifying pixel counts
#' directly therefore means different physical fires at every resolution: `c(10, ..., 1e5)`
#' pixels is 10--100,000 ha at 100 m but 57.6--576,000 ha at 240 m. That also defeats the point
#' of calibrating at several pixel sizes, which is to check the fitted parameters transfer
#' across resolutions -- only meaningful if the fires compared are the same physical size.
#'
#' @param pixelSize Pixel resolution, in metres.
#' @param areaHa Fire sizes, in hectares.
#' @param shapeRangeHa Range over which the shape target is supported; a warning is issued for
#'                     areas outside it, since the `log10(AREA)^4` term extrapolates sharply.
#'
#' @return Integer vector of fire sizes in pixels, with `areaHa` and `haPerPixel` attributes.
#'
#' @export
#' @examples
#' landmine_optim_fireSizes(240)                      ## v3 resolution
#' landmine_optim_fireSizes(100)                      ## the resolution originally calibrated at
landmine_optim_fireSizes <- function(pixelSize, areaHa = c(50, 200, 600, 1500, 3000),
                                     shapeRangeHa = c(10, 3000)) {
  stopifnot(length(pixelSize) == 1L, pixelSize > 0, all(areaHa > 0))
  haPerPx <- pixelSize^2 / 1e4
  px <- as.integer(round(areaHa / haPerPx))

  outside <- areaHa < shapeRangeHa[1] | areaHa > shapeRangeHa[2]
  if (any(outside)) {
    warning("fire size(s) ", paste(areaHa[outside], collapse = ", "), " ha fall outside the ",
            shapeRangeHa[1], "-", shapeRangeHa[2], " ha range Andison fitted; the shape target ",
            "extrapolates sharply beyond it.", call. = FALSE)
  }
  tooSmall <- px < 5L
  if (any(tooSmall)) {
    warning("fire size(s) ", paste(areaHa[tooSmall], collapse = ", "), " ha are under 5 pixels ",
            "at ", pixelSize, " m; shape and island statistics are dominated by quantization ",
            "at that size.", call. = FALSE)
  }
  structure(px, areaHa = areaHa, haPerPixel = haPerPx)
}

#' Resolution-dependent weight for the remnant-island term
#'
#' Andison (1996) rasterised at **50 m (quarter-hectare) pixels** -- his spatial database was
#' 260 x 422 = 114,920 quarter-hectare pixels. LandWeb runs much coarser, and remnant islands
#' are the statistic that suffers most, because a coarse grid cannot form small islands at all:
#'
#' | resolution | ha/pixel | smallest possible island |
#' | ---------- | -------- | ------------------------ |
#' | Andison 50 m | 0.25 | 0.25 ha |
#' | LandWeb 100 m | 1.00 | 1 ha |
#' | LandWeb 120 m | 1.44 | 1.44 ha |
#' | LandWeb 240 m | 5.76 | 5.76 ha |
#'
#' Andison's Table 3.6 puts **65%** of island area in the 2--5 ha class for fires under 1,000 ha
#' (24% for larger fires). At 240 m a single pixel is already 5.76 ha, so that entire class is
#' unreachable by construction and the island target cannot be met however the parameters are
#' set. Weighting the term down with coarsening resolution stops the optimiser chasing a target
#' the grid forbids, at the expense of statistics it can actually reproduce.
#'
#' @param pixelSize Pixel resolution, in metres.
#' @param refPixelSize Resolution the target statistics were measured at; defaults to Andison's
#'                     50 m.
#'
#' @return A weight in `(0, 1]`: `min(1, refPixelSize / pixelSize)`. 1.0 at 50 m, 0.42 at 120 m,
#'         0.21 at 240 m.
#'
#' @note This is a deliberate, documented judgement, not a derived quantity. The linear ratio is
#'   a pragmatic choice; the representability argument above justifies *down-weighting*, not this
#'   particular functional form. Record the weight used alongside any fit.
#'
#' @export
#' @examples
#' landmine_optim_islandWeight(c(50, 100, 120, 240))
landmine_optim_islandWeight <- function(pixelSize, refPixelSize = 50) {
  stopifnot(all(pixelSize > 0), refPixelSize > 0)
  pmin(1, refPixelSize / pixelSize)
}
