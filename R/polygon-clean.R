#' Clean up the LandWeb LTHFCs
#'
#' @param poly A polygon or character string identifying the path to polygon
#'
#' @param minFRI Numeric or integer, the minimum fire return interval kept.
#'               Intervals strictly below this become `NA`; an interval equal
#'               to `minFRI` is kept.
#' @export
.cleanLandWebStudyArea <- function(poly, minFRI = 40) {
  if (is.character(poly)) {
    pemisc::createPrjFile(poly)
    poly <- sf::st_read(poly)
  }

  stopifnot(any(c("LTHFC", "LTHRC") %in% names(poly)))

  ## Apparently, sometimes it is LTHFC, sometimes LTHRC. Renamed with names<-, which sf and terra both
  ## support: dplyr::rename() on a SpatVector needs tidyterra loaded to find its method.
  nms <- names(poly)
  if (all(c("LTHFC", "LTHRC") %in% nms)) {
    stop("the polygons have both an LTHFC and an LTHRC column", call. = FALSE)
  }
  nms[nms %in% c("LTHFC", "LTHRC")] <- "fireReturnInterval"
  names(poly) <- nms

  ## fires with Fire Return Interval 30 years are not correctly simulated; remove.
  ## Strictly BELOW `minFRI`, as documented: this was `<=`, which at the default of 40 also dropped
  ## every FRI-40 polygon -- 27 of them, 91,797 km2 of LTHFC v10, much of it NW Alberta.
  poly$fireReturnInterval[poly$fireReturnInterval < minFRI] <- NA

  return(poly)
}

#' Do an arbitrary set of operations on a polygon
#'
#' @param poly A polygon object, or a character string identifying the shapefile
#'             path to load, and clean.
#'
#' @param fn   A function identifying the type of cleaning to do.
#'
#' @param type If `fn` is not known, an character string can be specified to
#'             identify which `fn` to use.
#'             This must be a known type for this function.
#'
#' @param ...  Passed to `fn`
#'
#' @export
polygonClean <- function(poly, fn = NULL, type = NULL, ...) {
  if (is.null(fn)) {
    if (is.null(type)) {
      stop("Either fn or type must be specified")
    } else {
      if (type == "LandWeb") {
        fn <- .cleanLandWebStudyArea
      } else {
        stop("Unknown type")
      }
    }
  }
  poly <- fn(poly, ...)
}
