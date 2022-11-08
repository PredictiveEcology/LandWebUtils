#' Create a list of rasters in each rep, at each timestep, for each polygon area
#'
#' @note uses `future_lapply` internally to loop across `files`;
#'       set e.g., option `future.availableCores.fallback` appropriately for your system.
#'
#' @param files character vector giving paths to raster files
#' @param polys polygon object of class `sf`
#' @param names character vector giving the names of each of the subpolygons in `poly`
#' @param col character string giving the column name in `poly` to use
#' @param filter regex string giving partial filename in `files` to be filtered out
#'               (i.e., when extracting reps, times, etc. from filenames)
#'               (e.g., `'rstTimeSinceFire_'`, `'vegTypeMap_'`)
#'
#' @return list of `RasterLayer` objects with attributes `reps`, `times`, `polyNames`
#'
#' @export
#' @importFrom future.apply future_lapply
#' @importFrom purrr transpose
#' @importFrom raster crop mask raster
#' @importFrom tools file_path_sans_ext
rasterListByPoly <- function(files, polys, names, col, filter) {
  rlbp <- future.apply::future_lapply(files, function(f) {
    byPoly <- lapply(names, function(polyName) {
      subpoly <- polys[polys[[col]] == polyName, ]
      r <- raster::raster(f)
      rc <- raster::crop(r, subpoly)
      rcm <- raster::mask(rc, subpoly)
      rcm
    })
    names(byPoly) <- paste(tools::file_path_sans_ext(basename(f)), names , sep = "_") ## <filter>_yearXXXX_polyName

    byPoly
  }, future.packages = c("raster", "sp", "sf"))
  names(rlbp) <- basename(dirname(files)) ## repXX
  rlbp <- unlist(rlbp, recursive = FALSE, use.names = TRUE)

  labels <- purrr::transpose(strsplit(names(rlbp), "[.]"))
  labels1 <- unlist(labels[[1]])
  labels2 <- gsub(filter, "", unlist(labels[[2]]))
  labels2a <- purrr::transpose(strsplit(labels2, "_{1}"))
  labels2a1 <- unlist(labels2a[[1]])
  labels2a2 <- unlist(labels2a[[2]])

  attr(rlbp, "reps") <- as.integer(gsub("rep", "", labels1))
  attr(rlbp, "times") <- as.integer(gsub("year", "", labels2a1))
  attr(rlbp, "polyNames") <- labels2a2

  rlbp
}
