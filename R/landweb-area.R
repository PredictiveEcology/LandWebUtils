#' The LandWeb long-term historic fire cycle (LTHFC) map
#'
#' Downloads the v10 LTHFC map once, reads it, projects it, and keeps its fire cycle as `LTHFC` and each
#' polygon's area in m^2 as `area`. Prepare it with [polygonClean()]; [landweb_area()] takes its outline.
#'
#' The download needs Google Drive credentials (the LandWeb service account in a pipeline session).
#'
#' @param destinationPath directory under which the archive is downloaded and extracted (in `lthfc/`).
#' @param targetCRS coordinate reference system to project to.
#'
#' @return `SpatVector`.
#'
#' @export
landweb_lthfc <- function(destinationPath, targetCRS = LandWebCRS) {
  lthfcID <- "176yAq5NCfZZ5ZQX36zHcu0w3uh-V9qvf" ## landweb_ltfc_v10 (Google Drive)
  lthfcDir <- fs::dir_create(file.path(destinationPath, "lthfc"))
  lthfcZip <- file.path(lthfcDir, "landweb_ltfc_v10.zip")
  workflowtools::drive_download_once(googledrive::as_id(lthfcID), lthfcZip)
  workflowtools::archive_extract_once(lthfcZip, dir = lthfcDir)

  lthfc <- terra::vect(list.files(lthfcDir, "\\.shp$", full.names = TRUE)[[1]]) |>
    terra::project(targetCRS)
  lthfc <- lthfc[, "LTFC10"] ## v10 renamed the fire cycle column
  names(lthfc) <- "LTHFC"
  lthfc$area <- terra::expanse(lthfc, unit = "m")
  lthfc
}

#' The outline of the LandWeb area
#'
#' @param lthfc the cleaned LTHFC polygons: [polygonClean()] of [landweb_lthfc()].
#'
#' @return `sfc` polygon: the union of `lthfc`, made valid, without holes.
#'
#' @export
landweb_area <- function(lthfc) {
  sf::st_as_sf(lthfc) |> sf::st_union() |> sf::st_make_valid() |> nngeo::st_remove_holes()
}

#' Ecological units that touch an area, for fitting growth curves
#'
#' The polygons of the national ecological framework, at ecozone or ecoprovince level, that intersect
#' `area`, kept whole. `LandWeb_preamble` builds each study area's `studyAreaANPP` this way from the
#' buffered study area; the shared growth-curve fit uses the ecoprovinces that touch the LandWeb area.
#'
#' @param area `sf`, `sfc` or `SpatVector` polygons.
#' @param ecoLevel `"ecoprovince"` or `"ecozone"`.
#' @param destinationPath directory the framework's shapefile is downloaded to.
#'
#' @return `sf` of the intersecting polygons, in the coordinate reference system of `area`.
#'
#' @export
landweb_anpp_area <- function(area, ecoLevel = c("ecoprovince", "ecozone"), destinationPath) {
  ecoLevel <- match.arg(ecoLevel)
  url <- switch(
    ecoLevel,
    ecoprovince = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/province/ecoprovince_shp.zip",
    ecozone = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"
  )
  if (inherits(area, "SpatVector")) {
    area <- sf::st_as_sf(area)
  }
  eco <- reproducible::prepInputs(
    url = url,
    destinationPath = destinationPath,
    fun = "sf::st_read",
    overwrite = TRUE
  )
  eco <- sf::st_transform(eco, sf::st_crs(area))
  eco[which(lengths(sf::st_intersects(eco, area)) > 0), ]
}
