#' Clean up the LandWeb LTHFCs
#'
#' @param poly A polygon or character string identifying the path to polygon
#'
#' @param minFRI Numeric or integer, indicating the minimum fire return interval
#'               that will be part of the cleanup of polygon. Anything below
#'               this will be `NA`.
#' @export
.cleanLandWebStudyArea <- function(poly, minFRI = 40) {
  if (is.character(poly)) {
    pemisc::createPrjFile(poly)
    poly <- sf::st_read(poly)
  }

  stopifnot(any(c("LTHFC", "LTHRC") %in% names(poly)))

  ## Apparently, sometimes it is LTHFC, sometimes LTHRC; use LTHFC
  if (isTRUE("LTHRC" %in% names(poly))) {
    poly <- dplyr::rename(poly, LTHFC = "LTHRC")
  }

  poly <- dplyr::rename(poly, fireReturnInterval = "LTHFC")

  ## fires with Fire Return Interval 30 years are not correctly simulated; remove
  poly$fireReturnInterval[poly$fireReturnInterval <= minFRI] <- NA

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
