#' Join reporting polygons and intersect their features
#'
#' Join two reporting polygons, preserving their features;
#' concatenate the `Name.*` fields into a single `Name` field.
#' E.g., if `x` and `y` each contain 2 features, the resulting object will contain 4
#' features (corresponding to `x1.y1`, `x1.y2`, `x2.y1`, and `x2.y2`).
#'
#' @param x,y a `sf`
#'
#' @return an `sf` polygons object
#'
#' @export
joinReportingPolygons <- function(x, y) {
  if (is.null(x[["Name"]]) && !is.null(x[["Name.1"]]) && !is.null(x[["Name.2"]])) {
    z <- x

    z[["Name"]] <- paste(z[["Name.2"]], z[["Name.1"]])
    z[["Name.1"]] <- z[["Name.2"]] <- NULL

  } else {
    if (!is(x, "sf")) {
      x <- sf::st_as_sf(x)
    }
    if (!is(y, "sf")) {
      y <- sf::st_as_sf(y)
    }

    x <- sf::st_set_precision(x, 1e5) |> reproducible::fixErrors()
    y <- sf::st_set_precision(y, 1e5) |> reproducible::fixErrors()
    z <- sf::st_intersection(x, y)

    ## sfc_GEOMETRY may itself contain points, so filter them out
    z <- suppressWarnings(sf::st_collection_extract(z, "POLYGON"))

    ## ensure polygon name not duplicated in ACTIVE/PASSIVE x Caribou polygon names
    z[["Name"]] <- gsub("^(ACTIVE|PASSIVE).*", "\\1", z[["Name"]], ignore.case = TRUE)
    z[["Name.1"]] <- gsub("^(ACTIVE|PASSIVE).*", "\\1", z[["Name.1"]], ignore.case = TRUE)

    ## concatenate polygon names
    z[["Name"]] <- paste(z[["Name"]], z[["Name.1"]])
    z[["Name.1"]] <- NULL

    z <- as(z, "Spatial")
  }

  return(z)
}
