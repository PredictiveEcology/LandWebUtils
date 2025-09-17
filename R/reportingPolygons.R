#' Target CRS (projection) to use with LandWeb
#'
#' @note needs to be character, not `CRS` class, for downstream use with `data.table`
#'
#' @export
LandWebCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                    "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0") ## TODO: use SCANFI

#' Prepare reporting polygons
#'
#' - LandWeb FMA boundaries (FMAs);
#' - Alberta FMU boundaries (FMUs);
#' - Alberta Natural Subregions (ANSRs);
#'
#' @param destinationPath character specifying path to destination directory
#' @param targetCRS character specifying the CRS to reproject polygons to
#'
#' @export
#' @rdname prepReportingPolygons
prepFMAs <- function(destinationPath, targetCRS = LandWebCRS) {
  reproducible::prepInputs(
    url = "https://drive.google.com/file/d/1yCbq8rcRXCfUKHJGg-Fzlnrjl48LJfCO", ## 2020
    destinationPath = destinationPath,
    projectTo = targetCRS
  )
}

#' @export
#' @rdname prepReportingPolygons
prepFMUs <- function(destinationPath, targetCRS = LandWebCRS) {
  reproducible::prepInputs(
    url = "https://drive.google.com/open?id=1OH3b5pwjumm1ToytDBDI6jthVe2pp0tS", ## 2024-08 added C5
    destinationPath = destinationPath,
    projectTo = targetCRS
  )
}

#' @export
#' @rdname prepReportingPolygons
prepANSRs <- function(destinationPath, targetCRS = LandWebCRS) {
  reproducible::prepInputs(
    url = "https://drive.google.com/file/d/1hW6zy0CpUBdk-K2IAjzW4INjVl1J4aLJ",
    destinationPath = destinationPath,
    projectTo = targetCRS
  )
}

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
