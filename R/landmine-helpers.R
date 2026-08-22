utils::globalVariables(c(
  "area", "centreCell", "fireSize", "perimeter", "sizeCutoffs", "spawnNewActive"
))

#' Calculate the mean of a truncated Pareto distribution
#'
#' @param k TODO: description needed
#' @param lower TODO: description needed
#' @param upper TODO: description needed
#' @param alpha TODO: description needed
#'
#' @return TODO: description needed
#'
#' @export
meanTruncPareto <- function(k, lower, upper, alpha) {
  k * lower^k * (upper^(1 - k) - alpha^(1 - k)) / ((1 - k) * (1 - (alpha / upper)^k))
}

#' LandMine burn optimization function
#'
#' @param ros             `SpatRaster` of LandMine Rate Of Spread values.
#' @param centreCell      See `startCells` in [landmine_burn1()].
#' @param fireSize        See `fireSizes` in [landmine_burn1()].
#' @inheritParams landmine_burn1
#'
#' @return named list of length 2 containing:
#'  - `burnedMap`: `SpatRaster` of burned pixels;
#'  - `LM`: `data.frame` of patch statistics from `landscapemetrics`.
#'
#' @export
landmine_optim_burnFun <- function(
  ros,
  centreCell,
  fireSize,
  spawnNewActive,
  sizeCutoffs,
  spreadProb
) {
  stopifnot(requireNamespace("landscapemetrics", quietly = TRUE))

  burned <- landmine_burn1(
    landscape = ros,
    startCells = centreCell,
    fireSizes = fireSize,
    spreadProb = spreadProb,
    ## materialise once: `spread2()` re-reads a raster `spreadProbRel` on every spread step,
    ## which is O(ncell) per step. Bit-identical, ~2x faster per objective evaluation.
    spreadProbRel = as.numeric(terra::values(ros, mat = FALSE)),
    spawnNewActive = spawnNewActive,
    sizeCutoffs = sizeCutoffs
  )
  burnedMap <- terra::rast(ros) |> terra::as.int()
  burnedMap[] <- NA
  burnedMap[burned$pixels] <- burned$initialPixels

  ## NOTE: SDMTools is orphaned and no longer maintained;
  ## now using `landscapemetrics` instead to calculate the perimeter to area ratio of patches.
  ## Note too, that SDMTools calculations were sometimes incorrect:
  ## - double-counts internal edges (though not used here);
  ## - over-counts perimeter edges in some cases (inflating perim and perim:area).
  ##
  ## LM.orig <- SDMTools::PatchStat(burnedMap, cellsize = res(burnedMap)[1])

  ## ensure crs set, units in metres;
  if (terra::crs(burnedMap) == "") {
    terra::crs(burnedMap) <- "+proj=cea" ## cylindrical equal area projection

    # terra::crs(burnedMap, describe = TRUE)
    # terra::cellSize(burnedMap)
    # prod(res(burnedMap)) ## should be same as cellSize
    # landscapemetrics::check_landscape(burnedMap) ## verify it's ok
  }

  ## NOTE: to emulate SDMTools::PatchStat, need total perim and total area;
  ## not all fires are necessarily contiguous/connected and may result in >1 cluster ("patch");
  ## but note that perim and area are same whether clusters touch on diagonal or are
  ## separated by arbitrary number of pixels.

  # terra::trim(burnedMap) |> terra::plot()

  LM <- data.frame(
    patchID = landscapemetrics::lsm_p_para(burnedMap)$class |> unique(),
    perimeter = landscapemetrics::lsm_p_perim(burnedMap)$value |> sum(), ## m,
    area = (landscapemetrics::lsm_p_area(burnedMap)$value * 1e4) |> sum() ## ha * 1e4 = m^2
  ) |>
    dplyr::mutate(perim.area.ratio = perimeter / area) ## need m/m^2

  ## NOTE: `burnedMap` holds `initialPixels` (the ignition cell index) at every burned cell, so
  ## `sum(burnedMap[])` sums pixel IDs and is NOT a cell count. Use `nBurned` to ask whether a
  ## fire reached its target size.
  LM$nBurned <- sum(!is.na(terra::values(burnedMap, mat = FALSE)))
  LM$shapeIndex <- landmine_optim_shapeIndex(LM$perimeter, LM$area)

  list(burnedMap = burnedMap, LM = LM)
}

#' Setup a cluster for LandMine optimization
#'
#' @param nodes positive integer of length 1 specifying the number of threads
#'              to use on the current machine (`localhost`), or a character vector of
#'              hostnames on which to run worker copies.
#'
#' @param seed Optional integer passed to [parallel::clusterSetRNGStream()], so the per-worker
#'             RNG streams are reproducible. Without it the seed is drawn from ambient RNG
#'             state, which made calibration runs irreproducible unless the caller happened to
#'             call `set.seed()` first. Note reproducibility also requires a fixed worker count,
#'             since one stream is assigned per worker.
#'
#' @return a cluster object
#'
#' @export
landmine_optim_clusterSetup <- function(nodes = NULL, seed = NULL) {
  stopifnot(
    requireNamespace("parallel", quietly = TRUE),
    requireNamespace("parallelly", quietly = TRUE),
    !is.null(nodes)
  )

  nnodes <- if (is.numeric(nodes) && length(nodes) == 1) {
    as.integer(nodes)
  } else if (is.character(nodes)) {
    length(nodes)
  } else {
    stop("nodes must be an integer of length 1 or a character vector of nodenames")
  }

  ## NOTE: this used to also require `nnodes >= 70` and `nnodes %% 7 == 0`. Those are
  ## constraints on DEoptim's POPULATION SIZE (>= 10 per parameter, 7 parameters), not on the
  ## number of workers -- `NP` and the worker count are independent, and DEoptim distributes a
  ## population of any size across whatever workers exist. Conflating them made the calibration
  ## unrunnable on a 64-core host, and made both `NP` and the per-worker RNG streams depend on
  ## `availableCores()`, i.e. on the machine and its current load. Set `NP` explicitly instead.
  ## R can't create more than ~125 socket connections.
  stopifnot(
    nnodes >= 1,
    nnodes <= parallelly::availableCores(constraints = c("connections"))
  )

  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(nnodes)
  } else if (nnodes > 1) {
    cl <- parallelly::makeClusterPSOCK(nodes, autoStop = TRUE)
  } else {
    cl <- parallel::makeCluster(nnodes, type = "FORK")
  }

  ## seeding from `sample(1e8, 1)` drew from ambient RNG state, so streams were not
  ## reproducible unless the caller happened to set.seed() first.
  parallel::clusterSetRNGStream(cl, if (is.null(seed)) sample(1e8, 1) else seed)

  return(cl)
}

#' Export objects used for optimization and load packages on cluster
#'
#' @param cl    a cluster object or `NULL`
#' @param objs  character vector of names of objects to export
#' @param pkgs  character vector of packages to pre-load on the cluster nodes
#'
#' @return `NULL`. Invoked for its side effects.
#' @export
landmine_optim_clusterExport <- function(cl = NULL, objs = NULL, pkgs = NULL) {
  stopifnot(requireNamespace("parallel", quietly = TRUE))

  if (!is.null(cl)) {
    objs <- c(objs)
    parallel::clusterExport(cl, objs, envir = parent.frame())
    parallel::clusterExport(cl, "pkgs", envir = parent.frame())
    parallel::clusterEvalQ(cl, {
      lapply(pkgs, library, character.only = TRUE)
    })
  }

  data.table::setDTthreads(1L)

  return(NULL)
}

#' Wrapper function to setup cluster, export objects and load packages
#'
#' @inheritParams landmine_optim_clusterSetup
#' @inheritParams landmine_optim_clusterExport
#' @param reps integer. number of replicates to run.
#'
#' @return named list of length 2 containing:
#'  `cl`, a cluster object;
#'   `out`, a list of burn maps (aka `burnMapList`).
#' @export
landmine_optim_clusterWrap <- function(cl = NULL, nodes, reps, objs, pkgs) {
  stopifnot(requireNamespace("parallel", quietly = TRUE))

  if (is.null(cl)) {
    cl <- landmine_optim_clusterSetup(nodes)
  }
  landmine_optim_clusterExport(cl, objs, pkgs)
  objs <- mget(objs, envir = parent.frame())

  burnMapList <- parallel::clusterApplyLB(cl, reps, function(r) {
    ros <- terra::rast(ros)
    out <- LandWebUtils::landmine_optim_burnFun(
      ros = ros,
      centreCell = centreCell,
      fireSize = fireSize,
      spawnNewActive = spawnNewActive,
      sizeCutoffs = sizeCutoffs,
      spreadProb = spreadProb
    )
    out$burnedMap <- terra::wrap(out$burnedMap)

    return(out)
  })

  return(list(cl = cl, out = burnMapList))
}

#' LandMine objective functions
#'
#' `landmine_fitSN()` is used for the module.
#'
#' `landmine_fitSN2()` is an alternative version tries the optimization using fewer parameters,
#' to test whether a simpler version gets better/different results.
#' Although this version was not used for the final module, we preserve it here for posterity.
#'
#' @seealso [landmine_burn1()], [landmine_optim_burnFun()]
#'
#' @param sna (i.e., `spawnNewActive`) A numeric vector of length 4.
#'  These are the probabilities of creating spreading to 2 neighbours
#'  instead of the 1 default neighbour, each time step.
#'  The 4 values are for 4 different fire size conditions.
#'  See details in [landmine_optim_burnFun()].
#'
#' @param ros Character, specifying the file path to raster of LandMine Rate Of Spread values.
#'
#' @param centreCell Integer id of the centre (start) cell of `ros` raster.
#'  See `startCells` in [landmine_burn1()].
#'
#' @param fireSizes A numeric vector indicating the final size of each of the fires.
#'  See [landmine_burn1()].
#'
#' @param desiredPerimeterArea Numeric target perimeter-area ratio.
#'
#' @return Summary `data.table` of fit results.
#'
#' @export
#' @rdname landmine_fitSN
landmine_optim_fitSN <- function(
  sna,
  ros,
  centreCell,
  fireSizes = 10^(2:5),
  desiredPerimeterArea = 0.004
) {
  stopifnot(is.character(ros))

  ros <- terra::rast(ros)

  sizeCutoffs <- 10^sna[5:6]
  spreadProb <- sna[7]
  sna <- c(10^(sna[1]), 10^(sna[2]), 10^(sna[3]), 10^(sna[4]))
  bfs1 <- lapply(fireSizes, function(fireSize) {
    landmine_optim_burnFun(
      ros = ros,
      centreCell = centreCell,
      fireSize = fireSize,
      spawnNewActive = sna,
      sizeCutoffs = sizeCutoffs,
      spreadProb = spreadProb
    )
  })
  res <- lapply(seq(bfs1), function(bfCount) {
    abs(log(bfs1[[bfCount]]$LM[1, "perim.area.ratio"]) - log(desiredPerimeterArea)) +
      100 * (bfs1[[bfCount]]$LM[1, "nBurned"] < fireSizes[bfCount])
    ## NOTE: this previously read `sum(burnedMap[], na.rm = TRUE)`, which sums the ignition-cell
    ## INDEX stored at each burned cell (~5e9), not the number of cells burned -- so the penalty
    ## never fired in any calibration run. Fixing it CHANGES this objective: rows in
    ## `LandMine_DEoptim_params.csv` fitted before this are not comparable to rows fitted after.
  })
  a <- sum(unlist(res)) # * log10(fireSizes)) # weigh larger ones more
  attr(a, "bfs1") <- bfs1
  a
}

#' @param par         parameter vector of length 5
#' @param spreadProb  spread probability
#'
#' @export
#' @rdname landmine_fitSN
landmine_optim_fitSN2 <- function(
  par,
  ros,
  centreCell,
  fireSizes = 10^(2:5),
  desiredPerimeterArea = 0.003,
  spreadProb = 0.9
) {
  sizeCutoffs <- 10^c(par[4], par[5])
  bfs1 <- lapply(fireSizes, function(fireSize) {
    sna <- min(-0.15, par[1] + par[2] * log10(fireSize))
    sna <- 10^c(sna * par[3], sna * 2 * par[3], sna * 3 * par[3], sna * 4 * par[3])
    # sna <- -1
    landmine_optim_burnFun(ros, centreCell, fireSize, sna, sizeCutoffs, spreadProb)
  })
  res <- lapply(seq(bfs1), function(bfCount) {
    abs(log(bfs1[[bfCount]]$LM[1, "perim.area.ratio"]) - log(desiredPerimeterArea)) +
      100 * (bfs1[[bfCount]]$LM[1, "nBurned"] < fireSizes[bfCount])
    ## NOTE: this previously read `sum(burnedMap[], na.rm = TRUE)`, which sums the ignition-cell
    ## INDEX stored at each burned cell (~5e9), not the number of cells burned -- so the penalty
    ## never fired in any calibration run. Fixing it CHANGES this objective: rows in
    ## `LandMine_DEoptim_params.csv` fitted before this are not comparable to rows fitted after.
  })
  a <- sum(unlist(res))
  attr(a, "bfs1") <- bfs1
  a
}
