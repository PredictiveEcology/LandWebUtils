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
  k * lower^k * (upper^(1 - k) - alpha^(1 - k)) / ((1 - k) * (1 - (alpha/upper)^k))
}

#' LandMine burn optimization function
#'
#' @param ros             `RasterLayer` of LandMine Raster Of Spread values
#' @param centreCell      TODO
#' @param fireSize        TODO
#' @param spawnNewActive  TODO
#' @param sizeCutoffs     TODO
#' @param spreadProb      TODO
#'
#' @return named list of length 2 containing:
#'         `burnedMap`: `rasterLayer` of burned pixels;
#'         `LM`: `data.frame` of patch statistics from `SDMTools::PatchStats()`.
#'
#' @export
#' @importFrom raster raster res
landmine_optim_burnFun <- function(ros, centreCell, fireSize, spawnNewActive, sizeCutoffs, spreadProb) {
  burned <- landmine_burn1(
    landscape = ros,
    startCells = centreCell,
    fireSizes = fireSize,
    spreadProb = spreadProb,
    spreadProbRel = ros,
    spawnNewActive = spawnNewActive,
    sizeCutoffs = sizeCutoffs
  )
  burnedMap <- raster(ros)
  burnedMap[] <- NA

  burnedMap[burned$pixels] <- burned$initialPixels

  ## TODO: SDMTools is orphaned and no longer maintained;
  ## use landscapemetrics package instead to calculate the perimiter to area ratio of patches
  LM <- SDMTools::PatchStat(burnedMap, cellsize = res(burnedMap)[1])
  list(burnedMap = burnedMap, LM = LM)
}

#' Setup a cluster for LandMine optimization
#'
#' @param nodes positive integer of length 1 specifying the number of threads
#'              to use on the current machine (`localhost`), or a character vector of
#'              hostnames on which to run worker copies.
#'
#' @return a cluster objct
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
    stop("nodes must be an integer of length 1 or a character vector of nodenames")
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
    cl <- parallel::makeCluster(nodes, homogeneous = FALSE, verbose = TRUE)
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
landmine_optim_clusterWrap <- function(nodes, reps, objs, pkgs) {
  stopifnot(requireNamespace("parallel", quietly = TRUE))

  cl <- landmine_optim_clusterSetup(nodes)
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
#' @param sna                   TODO
#' @param ros                   TODO
#' @param centreCell            TODO
#' @param fireSizes             TODO
#' @param desiredPerimeterArea  TODO
#'
#' @return `data.table` (TODO)
#'
#' @export
#' @rdname landmine_fitSN
landmine_optim_fitSN <- function(sna, ros, centreCell, fireSizes = 10^(2:5),
                                 desiredPerimeterArea = 0.004) {
  sizeCutoffs <- 10^sna[5:6]
  spreadProb <- sna[7]
  sna <- c(10^(sna[1]), 10^(sna[2]), 10^(sna[3]), 10^(sna[4]))
  bfs1 <- lapply(fireSizes, function(fireSize) {
    landmine_optim_burnFun(ros, centreCell, fireSize,
                           sna, sizeCutoffs, spreadProb = spreadProb)
  })
  res <- lapply(seq(bfs1), function(bfCount) {
    abs(log(bfs1[[bfCount]]$LM[1, "perim.area.ratio"]) - log(desiredPerimeterArea)) +
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
    sna <- 10^c(sna*par[3], sna*2*par[3], sna*3*par[3], sna*4*par[3])
    #sna <- -1
    landmine_optim_burnFun(ros, centreCell, fireSize, sna, sizeCutoffs, spreadProb)
  })
  res <- lapply(seq(bfs1), function(bfCount) {
    abs(log(bfs1[[bfCount]]$LM[1, "perim.area.ratio"]) - log(desiredPerimeterArea)) +
      100*(sum(bfs1[[bfCount]]$burnedMap[], na.rm = TRUE) < fireSizes[bfCount])
    ## it needs to get to above 90,000 HA for it to count
  })
  a <- sum(unlist(res))
  attr(a, "bfs1") <- bfs1
  a
}
