utils::globalVariables(c("area", "perimeter"))

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
  k *
    lower^k *
    (upper^(1 - k) - alpha^(1 - k)) /
    ((1 - k) * (1 - (alpha / upper)^k))
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
landmine_optim_burnFun <- function(ros, centreCell, fireSize, spawnNewActive,
                                   sizeCutoffs, spreadProb) {
  stopifnot(requireNamespace("landscapemetrics", quietly = TRUE))

  burned <- landmine_burn1(
    landscape = ros,
    startCells = centreCell,
    fireSizes = fireSize,
    spreadProb = spreadProb,
    spreadProbRel = ros,
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

  list(burnedMap = burnedMap, LM = LM)
}

#' Setup a cluster for LandMine optimization
#'
#' @param nodes positive integer of length 1 specifying the number of threads
#'              to use on the current machine (`localhost`), or a character vector of
#'              hostnames on which to run worker copies.
#'
#' @return a cluster object
#'
#' @export
landmine_optim_clusterSetup <- function(nodes = NULL) {
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
    stop(
      "nodes must be an integer of length 1 or a character vector of nodenames"
    )
  }

  ## check the number of nodes used for the cluster:
  ## 1. at least 10 populations per parameter; 7 params;
  ## 2. ensure it's a multiple of the number of params (7);
  ## 3. R can't create more than ~125 socket connections.
  stopifnot(
    nnodes >= 70,
    nnodes %% 7 == 0,
    nnodes <= parallelly::availableCores(constraints = c("connections"))
  )

  if (Sys.info()["sysname"] == "Windows") {
    cl <- parallel::makeCluster(nnodes)
  } else if (nnodes > 1) {
    cl <- parallelly::makeClusterPSOCK(nodes, autoStop = TRUE)
  } else {
    cl <- parallel::makeCluster(nnodes, type = "FORK")
  }

  parallel::clusterSetRNGStream(cl, sample(1e8, 1))

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
    parallel::clusterExport(cl, objs)
    # env <- environment()
    parallel::clusterExport(cl, "pkgs", envir = parent.frame())
    parallel::clusterEvalQ(cl, {
      lapply(pkgs, library, character.only = TRUE)
      data.table::setDTthreads(1L)
    })
  } else {
    data.table::setDTthreads(1L)
  }

  return(NULL)
}

#' Wrapper function to setup cluster, export objects and load packages
#'
#' @inheritParams landmine_optim_clusterSetup
#' @inheritParams landmine_optim_clusterExport
#' @param reps integer. number of replicates to run.
#'
#' @return named list of length 2 containing:
#'   `cl` - a cluster object;
#'   `out` - a list of burn maps (aka `burnMapList`)
#' @export
landmine_optim_clusterWrap <- function(cl = NULL, nodes, reps, objs, pkgs) {
  stopifnot(requireNamespace("parallel", quietly = TRUE))

  if (is.null(cl)) {
    cl <- landmine_optim_clusterSetup(nodes)
  }
  landmine_optim_clusterExport(cl, objs, pkgs)
  objs <- mget(objs, envir = parent.frame())

  burnMapList <- parallel::clusterApplyLB(cl, reps, function(r) {
    do.call("landmine_optim_burnFun", objs)
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
#' @param ros `SpatRaster` of LandMine Rate Of Spread values.
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
landmine_optim_fitSN <- function(sna, ros, centreCell, fireSizes = 10^(2:5),
                                 desiredPerimeterArea = 0.004) {
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
    abs(
      log(bfs1[[bfCount]]$LM[1, "perim.area.ratio"]) - log(desiredPerimeterArea)
    ) +
      100 * (sum(bfs1[[bfCount]]$burnedMap[], na.rm = TRUE) < fireSizes[bfCount])
    ## it needs to get to above 90,000 HA for it to count
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
landmine_optim_fitSN2 <- function(par, ros, centreCell, fireSizes = 10^(2:5),
                                  desiredPerimeterArea = 0.003, spreadProb = 0.9) {
  sizeCutoffs <- 10^c(par[4], par[5])
  bfs1 <- lapply(fireSizes, function(fireSize) {
    sna <- min(-0.15, par[1] + par[2] * log10(fireSize))
    sna <- 10^c(
      sna * par[3],
      sna * 2 * par[3],
      sna * 3 * par[3],
      sna * 4 * par[3]
    )
    # sna <- -1
    landmine_optim_burnFun(
      ros,
      centreCell,
      fireSize,
      sna,
      sizeCutoffs,
      spreadProb
    )
  })
  res <- lapply(seq(bfs1), function(bfCount) {
    abs(log(bfs1[[bfCount]]$LM[1, "perim.area.ratio"]) - log(desiredPerimeterArea)) +
      100 * (sum(bfs1[[bfCount]]$burnedMap[], na.rm = TRUE) < fireSizes[bfCount])
    ## it needs to get to above 90,000 HA for it to count
  })
  a <- sum(unlist(res))
  attr(a, "bfs1") <- bfs1
  a
}
