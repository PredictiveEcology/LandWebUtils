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
#' @param date Date of the run, as a string. Defaults to today.
#'
#' @return `landmine_optim_params_append()` invisibly returns the full table as written.
#'
#' @export
#' @rdname landmine_optim_params
landmine_optim_params_append <- function(file, par, pixelSize,
                                         date = format(Sys.time(), "%Y-%m-%d")) {
  par <- as.numeric(unlist(unname(par)))
  invisible(landmine_optim_unpack(par)) ## validate via the single convention
  newRow <- data.frame(date = date, pixelSize = pixelSize, stringsAsFactors = FALSE)
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
#' @note The objective is stochastic and evaluated once per parameter vector, with no
#'   common random numbers, so repeated runs return *similar but different* values. In 43
#'   repeat fits recorded in `LandMine_DEoptim_params.csv`, four of the seven parameters
#'   scattered across more than 90% of their search interval, which indicates the objective
#'   does not identify them individually. Treat a single fit as one draw from a broad
#'   optimum rather than as a point estimate.
#'
#' @param pixelSize,n,file Passed to `landmine_optim_landscape()`.
#' @param fireSizes Target fire sizes (pixels) the objective fits against.
#' @param desiredPerimeterArea Target perimeter-to-area ratio (m/m^2).
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
                                     fireSizes = c(10, 100, 1000, 10000, 100000),
                                     desiredPerimeterArea = 0.003,
                                     lower = c(-2, -3, -3, -3, 1, 3.5, 0.75),
                                     upper = c(-0.1, -0.5, -0.5, -1, 3.5, 5, 1),
                                     NP = NULL, itermax = 200L, VTR = 0.001, strategy = 6L,
                                     cl = NULL, nodes = NULL, seed = NULL,
                                     paramsFile = NULL) {
  if (!requireNamespace("DEoptim", quietly = TRUE)) {
    stop("`DEoptim` is required to run the calibration; install it and retry.", call. = FALSE)
  }
  stopifnot(length(lower) == 7L, length(upper) == 7L, all(upper > lower))

  if (!is.null(seed)) {
    set.seed(seed)
  }
  ownCluster <- is.null(cl)
  if (ownCluster && !is.null(nodes)) {
    cl <- landmine_optim_clusterSetup(nodes = nodes)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  if (is.null(NP)) {
    NP <- if (is.null(cl)) 70L else length(cl)
  }

  lscape <- landmine_optim_landscape(pixelSize = pixelSize, n = n, file = file)

  res <- DEoptim::DEoptim(
    fn = landmine_optim_fitSN,
    lower = lower,
    upper = upper,
    control = DEoptim::DEoptim.control(
      VTR = VTR,
      NP = as.integer(NP),
      itermax = itermax,
      cluster = cl,
      strategy = strategy
    ),
    ros = lscape$file,
    centreCell = lscape$centreCell,
    fireSizes = fireSizes,
    desiredPerimeterArea = desiredPerimeterArea
  )

  if (!is.null(paramsFile)) {
    landmine_optim_params_append(
      file = paramsFile,
      par = res$optim$bestmem,
      pixelSize = pixelSize
    )
  }
  res
}
