#' Target CRS (projection) to use with LandWeb
#'
#' As of September 2025, this corresponds to the SCANFI's `NAD_1983_Canada_Lambert` projection.
#'
#' @note needs to be character, not `CRS` class, for downstream use with `data.table`
#'
#' @export
LandWebCRS <- paste(
  "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77",
  "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
)

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
  ## v10: FMAs_LandWebltfc_map_v10_LCC1 (supersedes FMA_Boundary_Updated_2024)
  .prepDriveVector("1c5vkNPG81DB5wkAVT3CP_h2jFCUDqCoL", destinationPath, targetCRS)
}

#' @export
#' @rdname prepReportingPolygons
prepFMUs <- function(destinationPath, targetCRS = LandWebCRS) {
  .prepDriveVector("1OH3b5pwjumm1ToytDBDI6jthVe2pp0tS", destinationPath, targetCRS)
}

#' @export
#' @rdname prepReportingPolygons
prepANSRs <- function(destinationPath, targetCRS = LandWebCRS) {
  .prepDriveVector("1hW6zy0CpUBdk-K2IAjzW4INjVl1J4aLJ", destinationPath, targetCRS)
}

#' @export
#' @rdname prepReportingPolygons
prepEcoregionLayer <- function(destinationPath, targetCRS = LandWebCRS) {
  .prepUrlVector(
    "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
    destinationPath, targetCRS,
    subdir = "Ecoregions"
  )
}

#' @export
#' @rdname prepReportingPolygons
prepEcoprovinceLayer <- function(destinationPath, targetCRS = LandWebCRS) {
  .prepUrlVector(
    "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/province/ecoprovince_shp.zip",
    destinationPath, targetCRS,
    subdir = "Ecoprovinces"
  )
}

## Download (googledrive direct, bypassing reproducible's broken Drive path /
## reproducible #447), extract, load, and reproject a Drive-hosted vector. Returns
## `sf` to match the legacy map/postProcess/joinReportingPolygons machinery.
.prepDriveVector <- function(id, destinationPath, targetCRS) {
  dir <- reproducible::checkPath(file.path(destinationPath, id), create = TRUE)
  zip <- file.path(dir, paste0(id, ".zip"))
  workflowtools::drive_download_once(googledrive::as_id(id), zip)
  workflowtools::archive_extract_once(zip, dir = dir)
  shp <- list.files(dir, "\\.shp$", full.names = TRUE)[[1]]
  v <- terra::project(terra::makeValid(terra::vect(shp)), targetCRS)
  sf::st_as_sf(v)
}

## URL analog of `.prepDriveVector()`: download (once), extract (once), load, and
## reproject a zipped shapefile from a plain URL (e.g. the AAFC ecostrat
## ecoregion/ecoprovince layers), composing the idempotent `workflowtools` helpers.
.prepUrlVector <- function(url, destinationPath, targetCRS, subdir) {
  dir <- reproducible::checkPath(file.path(destinationPath, subdir), create = TRUE)
  zip <- file.path(dir, basename(url))
  workflowtools::download_once(url, zip)
  workflowtools::archive_extract_once(zip, dir = dir)
  shp <- list.files(dir, "\\.shp$", full.names = TRUE)[[1L]]
  v <- terra::project(terra::makeValid(terra::vect(shp)), targetCRS)
  sf::st_as_sf(v)
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

#' Candidate reporting-polygon source layers
#'
#' The reporting-polygon layers considered for LandWeb/NRV summaries. Each is
#' fetched, clipped to the study area with [spatialutils::prep_vector()], and
#' kept only if it intersects (see [buildReportingPolygons()]). `source` is
#' either `"drive"` (a Google Drive file id) or `"url"` (a direct URL);
#' `labelCol` is the attribute used for sub-region labels.
#'
#' @return A `tibble` with columns `key`, `source`, `id`, and `labelCol`.
#' @export
reportingPolygonLayers <- function() {
  tibble::tribble(
    ~key                               , ~source , ~id                                                                   , ~labelCol    ,
    ## NB: the source FMA shapefile's `Name` column is CORRUPT -- every row holds the
    ## deparsed whole-column vector as a literal string ("c(NA, NA, ..., \"Fort Providence\",
    ## ...)"), a bug in whatever built FMAs_LandWebltfc_map_v10. Use the clean per-row
    ## `FMA_NAME` (FMA holder) instead.
    "FMA Boundaries Updated"           , "drive" , "1c5vkNPG81DB5wkAVT3CP_h2jFCUDqCoL"                                   , "FMA_NAME"   ,
    "Caribou Ranges"                   , "drive" , "1gwqq3TO-vKfTR3Za7bf7GWW7BobbAusQ"                                   , "RANGE_NAME" ,
    "Parks"                            , "drive" , "10-TpJaEOUCN6MNISWhpU-CRLpadyHDS2"                                   , "Name"       ,
    "Alberta Natural Subregions"       , "drive" , "1hW6zy0CpUBdk-K2IAjzW4INjVl1J4aLJ"                                   , "Name"       ,
    "BC Biogeoclimatic zones"          , "drive" , "1NS15Gd7dHEhvPOy-Ol_LBtf-4Ch6mPnS"                                   , "ZONE_NAME"  ,
    "Northwest Territories Ecoregions" , "drive" , "1iRAQfARkmS6-XVHFnTkB-iltzMNPAczC"                                   , "ECO4_NAM_1" ,
    "National Ecozones"                , "url"   , "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"     , "ZONE_NAME"  ,
    "National Ecoregions"              , "url"   , "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip" , "REGION_NAM"
  )
}

#' Build the reporting-polygons named list
#'
#' Fetch each candidate reporting layer ([reportingPolygonLayers()]), clip it to
#' `studyArea` with [spatialutils::prep_vector()], and keep only those that
#' intersect -- a layer that does not overlap the study area is dropped
#' (`nrow == 0`), which replaces the former per-study-area hardcoding. Downloads
#' use the idempotent `workflowtools` `*_once` helpers (bypassing
#' `reproducible::prepInputs`). The result is a named list keyed by layer name,
#' suitable as the `reportingPolygons` input to `NRV_summary`; each layer's
#' `labelCol` is copied to a common `Name` column. The `"CC SAM"`/`"CC TSF"`
#' (current condition) and `"ecoregionLayer"` entries that `NRV_summary` also
#' expects come from the simulation, not from here, and are merged in separately.
#'
#' @param studyArea A `SpatVector` (or source readable by [terra::vect()]).
#' @param destinationPath Directory for downloads/extraction.
#' @param targetCRS Target CRS (default [LandWebCRS]).
#' @param layers Candidate-layer table (default [reportingPolygonLayers()]).
#'
#' @return A named `list` of `SpatVector`s, one per intersecting layer.
#' @export
buildReportingPolygons <- function(
  studyArea,
  destinationPath,
  targetCRS = LandWebCRS,
  layers = reportingPolygonLayers()
) {
  if (!inherits(studyArea, "SpatVector")) {
    studyArea <- terra::vect(studyArea)
  }

  out <- lapply(seq_len(nrow(layers)), function(i) {
    lyr <- layers[i, ]
    dir <- reproducible::checkPath(file.path(destinationPath, make.names(lyr$key)), create = TRUE)
    if (identical(lyr$source, "drive")) {
      dest <- file.path(dir, paste0(make.names(lyr$key), ".zip"))
      workflowtools::drive_download_once(googledrive::as_id(lyr$id), dest)
    } else {
      dest <- file.path(dir, basename(lyr$id))
      workflowtools::download_once(lyr$id, dest)
    }
    workflowtools::archive_extract_once(dest, dir = dir)

    shp <- list.files(dir, "\\.(shp|gpkg)$", full.names = TRUE)
    if (length(shp) == 0L) {
      return(NULL)
    }
    v <- spatialutils::prep_vector(terra::vect(shp[[1]]), studyArea, crs = targetCRS)
    if (nrow(v) == 0L) {
      return(NULL)
    }
    if (!is.null(lyr$labelCol) && lyr$labelCol %in% names(v)) {
      ## NB: v[[col]] returns a 1-column data.frame; as.character() on that deparses
      ## the whole column into one string recycled to every row. Use `[[1]]` to pull
      ## the vector so each feature keeps its own label.
      v$Name <- as.character(v[[lyr$labelCol]][[1]])
    }
    v
  })

  names(out) <- layers$key
  out[!vapply(out, is.null, logical(1))]
}

#' Candidate active/passive landbase-status source layers
#'
#' The per-FMA landbase-status ("active/passive", a.k.a.
#' "contributing/non-contributing") coverages, each fetched, dissolved by status,
#' and clipped to the study area by [buildLandbasePolygons()]. This generalises
#' the per-study-area landbase processing formerly hard-coded across the
#' `LandWeb_preamble` study-area helpers (`SprayLake.R`, `WestFraser.R`,
#' `Tolko.R`, `SundreFP.R`, `Manning.R`).
#'
#' Each landbase applies to only a subset of study areas (not every FMA has a
#' landbase coverage), so `applies_to` is a regular expression matched against
#' the study-area name: [buildLandbasePolygons()] fetches a source **only** when
#' its `applies_to` matches, avoiding the download of large (multi-GB) coverages
#' that could not intersect the study area anyway. `status_col` is the source
#' attribute holding the landbase status; `layer` names the layer for
#' multi-layer File Geodatabases (`NA` = the source's only/first layer).
#'
#' @return A `tibble` with columns `key`, `source`, `id`, `layer`, `status_col`,
#'   and `applies_to`.
#' @export
landbaseLayers <- function() {
  tibble::tribble(
    ~key                             , ~source , ~id                                 , ~layer            , ~status_col    , ~applies_to           ,
    "Spray Lake C5 Landbase"         , "drive" , "1FpMg6dJ4eblMjMkEHdcoFAYSDyo1OvTv" , "lb_20230901_tsa" , "f_active"     , "SprayLake"           ,
    "West Fraser Blue Ridge Landbase", "drive" , "1Mk5L6287sKFGLY4ZfwWIUAczF5AGqAOV" , NA_character_     , "LBC_LBStatus" , "BlueRidge"           ,
    "West Fraser N CLS S17 Landbase" , "drive" , "1XrF9ygQruC2FsUulWhDxR-nD3Cd4eu7B" , NA_character_     , "F_CONDITIO"   , "WestFraser_N"        ,
    "West Fraser N CLS S20 Landbase" , "drive" , "17fZw80w3n2jIKRP1-X6gq8tyOjSWh0ky" , NA_character_     , "F_CONDITIO"   , "WestFraser_N"        ,
    "West Fraser N CLS S21 Landbase" , "drive" , "1akMUL-lRumTfmmWG7WF7-9KDrFxTPN3Z" , NA_character_     , "F_CONDITIO"   , "WestFraser_N"        ,
    "Tolko AB North Landbase"        , "drive" , "1qzWhMR-nP_2N2KLL84ZHCVGNNlWxd2bZ" , NA_character_     , "LBC_Landbase" , "Tolko_AB_N|tolko_AB_N",
    "Sundre FP Landbase"             , "drive" , "1oh_w9nALKufCQXb1PR3VIlChPHPysHF4" , NA_character_     , "LBC_LBStatus" , "Sundre"              ,
    "Manning Landbase"               , "drive" , "1lY0p6Ms84paja9p1lmGXz5jaCgv2_VkY" , NA_character_     , "LBC_LBStat"   , "Manning"
  )
}

## Locate a readable vector source inside an extracted download directory:
## a File Geodatabase directory (`*.gdb`) if present, else a shapefile/GeoPackage.
.findLandbaseSource <- function(dir) {
  gdb <- list.files(dir, pattern = "\\.gdb$", full.names = TRUE, include.dirs = TRUE)
  gdb <- gdb[dir.exists(gdb)]
  if (length(gdb) > 0L) {
    return(gdb[[1]])
  }
  vec <- list.files(dir, pattern = "\\.(shp|gpkg)$", full.names = TRUE, recursive = TRUE)
  if (length(vec) > 0L) {
    return(vec[[1]])
  }
  NULL
}

#' Build the active/passive landbase-status reporting polygons
#'
#' For each candidate landbase source ([landbaseLayers()]) whose `applies_to`
#' matches `studyAreaName`, download it (idempotent `workflowtools` `*_once`
#' helpers), then dissolve it by landbase status and clip it to `studyArea` with
#' [spatialutils::prep_landbase()] -- which pushes the column and spatial
#' filters down to the read, so a large coverage is reduced to the study area
#' before any geometry work. Sources that do not intersect the study area (or
#' whose `applies_to` does not match) are skipped; the result is a named list of
#' `SpatVector`s keyed by layer name, each with a `Name` column holding the
#' landbase status, suitable for merging into the `reportingPolygons` input to
#' `NRV_summary`.
#'
#' @param studyArea A `SpatVector` (or a source readable by [terra::vect()]).
#' @param studyAreaName Character. The study-area name, matched against each
#'   layer's `applies_to` to decide which landbases to fetch.
#' @param destinationPath Directory for downloads/extraction.
#' @param targetCRS Target CRS (default [LandWebCRS]).
#' @param layers Candidate-layer table (default [landbaseLayers()]).
#'
#' @return A named `list` of `SpatVector`s, one per applicable, intersecting
#'   landbase; empty when none apply.
#' @seealso [buildReportingPolygons()], [landbaseLayers()]
#' @export
buildLandbasePolygons <- function(
  studyArea,
  studyAreaName,
  destinationPath,
  targetCRS = LandWebCRS,
  layers = landbaseLayers()
) {
  if (!inherits(studyArea, "SpatVector")) {
    studyArea <- terra::vect(studyArea)
  }

  ## gate: only consider landbases applicable to this study area
  applicable <- vapply(
    layers$applies_to,
    function(p) any(grepl(p, studyAreaName)),
    logical(1)
  )
  layers <- layers[applicable, , drop = FALSE]
  if (nrow(layers) == 0L) {
    return(list())
  }

  out <- lapply(seq_len(nrow(layers)), function(i) {
    lyr <- layers[i, ]
    dir <- reproducible::checkPath(file.path(destinationPath, make.names(lyr$key)), create = TRUE)
    dest <- file.path(dir, paste0(make.names(lyr$key), ".zip"))
    if (identical(lyr$source, "drive")) {
      workflowtools::drive_download_once(googledrive::as_id(lyr$id), dest)
    } else {
      workflowtools::download_once(lyr$id, dest)
    }
    workflowtools::archive_extract_once(dest, dir = dir)

    src <- .findLandbaseSource(dir)
    if (is.null(src)) {
      return(NULL)
    }

    layerName <- if (is.na(lyr$layer)) NULL else lyr$layer
    v <- spatialutils::prep_landbase(
      src,
      studyArea,
      status_col = lyr$status_col,
      layer = layerName,
      crs = targetCRS
    )
    if (nrow(v) == 0L) {
      return(NULL)
    }
    v$Name <- as.character(v$lbstatus) ## $ gives the vector; [["lbstatus"]] would deparse (see above)
    v
  })

  names(out) <- layers$key
  out[!vapply(out, is.null, logical(1))]
}
