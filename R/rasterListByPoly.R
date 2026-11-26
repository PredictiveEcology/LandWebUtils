#' Create a list of rasters in each rep, at each timestep, for each polygon area
#'
#' @note uses `future_lapply` internally to loop across `files`;
#'       set e.g., option `future.availableCores.fallback` appropriately for your system.
#'
#' @param files character vector giving paths to raster files
#' @param polys polygon object of class `sf`
#' @param names character vector giving the names of each of the subpolygons in `poly`
#' @param col character string giving the column name in `poly` to use
#' @param filter unused; kept for compatibility. The replicate, time and polygon name of each
#'   raster are parsed from its folder (`repXX`) and file name (`<prefix>_year<YYYY>`).
#'
#' @return list of `SpatRaster` objects, named `<repXX>.<prefix>_year<YYYY>_<polyName>`,
#'   with attributes `reps`, `times`, `polyNames`
#'
#' @export
rasterListByPoly <- function(files, polys, names, col, filter) {
  rlbp <- future.apply::future_lapply(
    files,
    function(f) {
      byPoly <- lapply(names, function(polyName) {
        subpoly <- polys[polys[[col]] == polyName, ]
        terra::rast(f) |> terra::crop(subpoly, mask = TRUE)
      })
      ## <prefix>_yearXXXX_polyName
      names(byPoly) <- paste(tools::file_path_sans_ext(basename(f)), names, sep = "_")

      byPoly
    },
    future.packages = c("terra", "sf")
  )
  names(rlbp) <- basename(dirname(files)) ## repXX
  rlbp <- unlist(rlbp, recursive = FALSE, use.names = TRUE)

  ## one anchored match per label, so a polygon name may contain "." or "_"
  labels <- names(rlbp)
  m <- regmatches(labels, regexec("^(.*?)\\.(.*)_year([0-9]+)_(.*)$", labels))
  bad <- lengths(m) != 5L
  if (any(bad)) {
    stop(
      "cannot parse raster label(s): ",
      paste(utils::head(labels[bad], 3L), collapse = ", "),
      " (expected '<repXX>.<prefix>_year<YYYY>_<polyName>')",
      call. = FALSE
    )
  }
  part <- function(i) vapply(m, `[[`, character(1L), i)

  attr(rlbp, "reps") <- as.integer(sub("^rep", "", part(2L)))
  attr(rlbp, "times") <- as.integer(part(4L))
  attr(rlbp, "polyNames") <- part(5L)

  rlbp
}
