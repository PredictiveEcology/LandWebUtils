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

## Filesystem-safe slug for generated directory/file names. NB: deliberately NOT
## `make.names()`, which maps spaces to DOTS ("National Ecozones" ->
## "National.Ecozones"); dot-laden names are awkward on Windows (and for Windows
## users), so collapse any run of non-alphanumerics to a single underscore.
.slug <- function(x) {
  gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
}

## Locate an extracted vector file. Searches RECURSIVELY: several source archives
## unpack into a subdirectory (e.g. the AAFC national ecoregion/ecozone zips extract
## to `<dir>/Ecoregions/ecoregions.shp`), and a non-recursive `list.files()` silently
## found nothing there -- which callers could not distinguish from a layer that
## legitimately does not intersect the study area.
.findVectorFile <- function(dir, pattern = "\\.(shp|gpkg|geojson)$") {
  list.files(dir, pattern, full.names = TRUE, recursive = TRUE)
}

## `source` gates which fetch path runs, and anything unrecognised would be treated as a zipped URL --
## so a typo ("geojsn") would silently take the archive path and surface later as an opaque "no vector
## file found". Fail at the table instead. "assembly" means the layer is built by a dedicated assembler
## rather than fetched from one file (see `.ASSEMBLERS`).
.LAYER_SOURCES <- c("drive", "url", "geojson", "assembly")

.validateLayerSources <- function(source) {
  bad <- setdiff(source, .LAYER_SOURCES)
  if (length(bad)) {
    stop(
      "unknown source(s) ", paste(sQuote(bad), collapse = ", "), "; must be one of ",
      paste(sQuote(.LAYER_SOURCES), collapse = ", "), ".", call. = FALSE
    )
  }
  invisible(source)
}

## Fetch one layer to disk and return the vector file to read. Shared by `buildReportingPolygons()` and
## `buildCaribouRanges()` so the three transports behave identically in both.
##   "drive"   -- Google Drive id, zipped shapefile
##   "url"     -- direct URL to a zipped shapefile
##   "geojson" -- URL returning GeoJSON directly (an ArcGIS REST `query` or an OGC WFS
##                `outputFormat=application/json`). No archive step, and the download CANNOT be named
##                from `basename(url)`: a query URL's basename is the query string itself.
.fetchVectorFile <- function(source, id, dir, slug) {
  if (identical(source, "drive")) {
    dest <- file.path(dir, paste0(slug, ".zip"))
    workflowtools::drive_download_once(googledrive::as_id(id), dest)
    workflowtools::archive_extract_once(dest, dir = dir)
  } else if (identical(source, "geojson")) {
    dest <- file.path(dir, paste0(slug, ".geojson"))
    workflowtools::download_once(id, dest)
  } else {
    dest <- file.path(dir, basename(id))
    workflowtools::download_once(id, dest)
    workflowtools::archive_extract_once(dest, dir = dir)
  }
  .findVectorFile(dir)
}

## Download (googledrive direct, bypassing reproducible's broken Drive path /
## reproducible #447), extract, load, and reproject a Drive-hosted vector. Returns
## `sf` to match the legacy map/postProcess/joinReportingPolygons machinery.
.prepDriveVector <- function(id, destinationPath, targetCRS) {
  dir <- reproducible::checkPath(file.path(destinationPath, id), create = TRUE)
  zip <- file.path(dir, paste0(id, ".zip"))
  workflowtools::drive_download_once(googledrive::as_id(id), zip)
  workflowtools::archive_extract_once(zip, dir = dir)
  shp <- .findVectorFile(dir, "\\.shp$")
  if (length(shp) == 0L) {
    stop("no shapefile found under '", dir, "' (Drive id '", id, "').", call. = FALSE)
  }
  v <- terra::project(terra::makeValid(terra::vect(shp[[1L]])), targetCRS)
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
  shp <- .findVectorFile(dir, "\\.shp$")
  if (length(shp) == 0L) {
    stop("no shapefile found under '", dir, "' (from '", url, "').", call. = FALSE)
  }
  v <- terra::project(terra::makeValid(terra::vect(shp[[1L]])), targetCRS)
  sf::st_as_sf(v)
}

#' Join reporting polygons and intersect their features
#'
#' Join two reporting polygons, preserving their features;
#' concatenate the `Name.*` fields into a single `Name` field.
#' E.g., if `x` and `y` each contain 2 features, the resulting object will contain 4
#' features (corresponding to `x1.y1`, `x1.y2`, `x2.y1`, and `x2.y2`).
#'
#' The concatenated `Name` is what makes the crossings work downstream: both
#' `nrvtools::calculateLandWebMetrics()` and `nrvtools::subregion_forested_area()`
#' group by **name**, not by feature, so a sub-region split into several disjoint
#' parts inside a tenure is tallied as one unit for free.
#'
#' @section Geometry repair and the intersection engine:
#' Repair goes through [spatialutils::repair_geoms()] (`terra::makeValid()` on just
#' the invalid subset), **not** [sf::st_make_valid()]: the latter collapses the
#' reversed ring-winding-order polygons present in the v10 FMA layer (e.g. Prince
#' Albert, -31,600 km^2 -> a 4 km^2 sliver) that terra correctly re-orients.
#'
#' The intersection itself uses [terra::intersect()] rather than v2's
#' [sf::st_intersection()]. On the real v10 layers `st_intersection()` throws GEOS
#' `TopologyException: side location conflict` on some already-repaired tenure x
#' Caribou pairs (e.g. Dawson Creek TSA), which silently cost those units their
#' reporting polygons; `terra::intersect()` handles the same geometries. Using one
#' engine (the same one the repair already runs through) also keeps the result
#' deterministic, rather than depending on which pairs happen to trip GEOS.
#'
#' @param x,y an `sf` or `SpatVector` (anything [terra::vect()] accepts).
#'
#' @return a polygons object of the same class as `x` (`sf` unless `x` was a
#'   `SpatVector`).
#'
#' @export
joinReportingPolygons <- function(x, y) {
  asVect <- inherits(x, "SpatVector")

  if (is.null(x[["Name"]]) && !is.null(x[["Name.1"]]) && !is.null(x[["Name.2"]])) {
    ## already-crossed input that carries the two source names but no combined one
    ## (the triple-crossing path): just re-concatenate, no geometry work.
    z <- x
    z$Name <- paste(.vecCol(z, "Name.2"), .vecCol(z, "Name.1"))
    return(.dropCols(z, c("Name.1", "Name.2")))
  }

  ## two layers that simply do not overlap is a normal, expected outcome here (most tenures
  ## meet most sub-region layers nowhere), so terra's "[intersect] no intersection" is muffled
  ## -- narrowly, so a genuine geometry warning still surfaces.
  z <- withCallingHandlers(
    terra::intersect(spatialutils::repair_geoms(x), spatialutils::repair_geoms(y)),
    warning = function(w) {
      if (grepl("no intersection", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  ## terra disambiguates the clashing `Name` columns as `Name_1` (from x) / `Name_2` (from y);
  ## restore the sf-style `Name` / `Name.1` the rest of this function is written against.
  if (all(c("Name_1", "Name_2") %in% names(z))) {
    names(z)[names(z) == "Name_1"] <- "Name"
    names(z)[names(z) == "Name_2"] <- "Name.1"
  }
  if (nrow(z) == 0L) {
    return(if (asVect) z else sf::st_as_sf(z))
  }

  ## an intersection can produce points/lines along shared edges; keep polygons only.
  if (!identical(terra::geomtype(z), "polygons")) {
    z <- z[terra::geomtype(z) == "polygons", ]
  }

  ## Ensure the tenure name is not repeated when crossing an ACTIVE/PASSIVE landbase (already
  ## "<status> <tenure>") with another already-tenure-crossed layer (e.g. "<range> <tenure>"):
  ## reduce the landbase side back to the bare status so the concatenation reads
  ## "<status> <range> <tenure>" rather than repeating the tenure twice.
  nm <- gsub("^(ACTIVE|PASSIVE).*", "\\1", .vecCol(z, "Name"), ignore.case = TRUE)
  nm2 <- .vecCol(z, "Name.1")
  if (!is.null(nm2)) {
    ## concatenate polygon names
    nm <- paste(nm, gsub("^(ACTIVE|PASSIVE).*", "\\1", nm2, ignore.case = TRUE))
    z <- .dropCols(z, "Name.1")
  }
  z$Name <- nm

  if (asVect) z else sf::st_as_sf(z)
}

## Pull one attribute as a plain vector, for `sf` and `SpatVector` alike. NB: for a
## `SpatVector`, `v[["col"]]` returns a 1-COLUMN DATA.FRAME -- passing that to `gsub()`/`paste()`
## deparses the whole column into one string and recycles it to every feature, which silently
## turned every crossed name into `c("A", "B") c("T1", "T1")`.
.vecCol <- function(v, nm) {
  x <- v[[nm]]
  if (is.data.frame(x)) x[[1L]] else x
}

## Drop attribute column(s), for `sf` and `SpatVector` alike (`v$col <- NULL` does not remove a
## `SpatVector` column). Geometry is preserved.
.dropCols <- function(v, nms) {
  keep <- setdiff(names(v), nms)
  if (inherits(v, "SpatVector")) v[, keep] else v[, c(keep, attr(v, "sf_column"))]
}

#' Candidate reporting-polygon source layers
#'
#' The reporting-polygon layers considered for LandWeb/NRV summaries. Each is
#' fetched, clipped to the study area with [spatialutils::prep_vector()], and
#' kept only if it intersects (see [buildReportingPolygons()]). `source` is one of
#' `"drive"` (a Google Drive file id), `"url"` (a direct URL to a zipped
#' shapefile), or `"geojson"` (a URL returning GeoJSON directly, e.g. an ArcGIS
#' REST `query` endpoint -- fetched with no archive-extraction step);
#' `labelCols` names the attribute(s) holding the sub-region label.
#'
#' `labelCols` is a **comma-separated priority list**, coalesced left to right, because the
#' v10 layers were assembled by merging per-jurisdiction sources and several name their
#' features in different columns depending on the province. The Caribou ranges are the sharp
#' case: `RANGE_NAME` is populated for only 21 of the 74 ranges (Ontario and Manitoba) ---
#' British Columbia uses `HERD_NAME`, Alberta `LOCALRANGE`, Saskatchewan `CONUNIT`/`RNGEUNIT`
#' --- so labelling on `RANGE_NAME` alone left 53 ranges unnamed. Unnamed features are dropped
#' downstream, which silently removed every tenure x Caribou reporting unit outside ON/MB.
#'
#' `NAME_SHORT` is the curated short token that names the layer everywhere it is
#' keyed: it is the name of the layer's entry in the [buildReportingPolygons()]
#' list, the second half of a crossed unit's name
#' (`"<tenure NAME_SHORT> <layer NAME_SHORT>"`, see
#' [buildCrossedReportingPolygons()]), and --- via `.slug()` --- the `refCode` that
#' keys the parquet aggregates and figure files. `key` is retained as the human
#' title for reports. See [validateShortNames()] for the invariants enforced.
#'
#' `isTenure` marks the **tenure** layer --- the forest-management boundaries every
#' other layer is crossed against ([buildCrossedReportingPolygons()]). Exactly one
#' layer may carry it.
#'
#' `cross` marks the layers crossed with each tenure by
#' [buildCrossedReportingPolygons()], reproducing the v2 per-tenure reporting units.
#' Every non-tenure layer is crossed: reporting is `studyAreaReporting` intersected
#' by tenure, so a stand-alone layer covering ground *outside* the tenures has
#' nothing to report against.
#'
#' A reporting layer may be assembled from **several sources**: rows sharing a
#' `NAME_SHORT` are built independently and then merged into one layer (reduced to the
#' common `Name` column). `where` optionally keeps only the features whose label matches
#' a regular expression, so a supplementary source can contribute just the part the
#' primary source is missing. Removing a supplement later is a one-row deletion.
#'
#' @return A `tibble` with columns `key`, `NAME_SHORT`, `isTenure`, `cross`,
#'   `source`, `id`, `labelCols`, and `where`.
#' @export
reportingPolygonLayers <- function() {
  out <- tibble::tribble(
    ~key                               , ~NAME_SHORT     , ~isTenure , ~cross , ~source , ~id                                                                   , ~labelCols   , ~where,
    ## NB: the source FMA shapefile's `Name` column is CORRUPT for most rows -- they hold the
    ## deparsed whole-column vector as a literal string ("c(NA, NA, ..., \"Fort Providence\",
    ## ...)"), a bug in whatever built FMAs_LandWebltfc_map_v10. `labelCols` is a FALLBACK here:
    ## the tenure layer is labelled by `.fmaMemberIdentity()`, which coalesces the per-jurisdiction
    ## identity fields (AB/NT `FMA_NAME`, ON `FMU_NAME`, BC `TSA_NUMB_1`, SK `FOREST_NAM`,
    ## MB `FML_NAME`) -- `FMA_NAME` alone is NA for every non-AB member, which would have left the
    ## tenure layer empty for the BC/SK/ON/MB study-area groups.
    "FMA Boundaries Updated"           , "FMA"           , TRUE      , FALSE  , "drive" , "1c5vkNPG81DB5wkAVT3CP_h2jFCUDqCoL"                                   , "FMA_NAME"   , NA_character_,
    ## Caribou ranges are ASSEMBLED from the six jurisdictional sources rather than fetched from one
    ## file -- see `caribouRangeLayers()` / `buildCaribouRanges()` for the sources and the reasoning.
    ## This replaces (a) Julie's pre-combined v10 layer, whose manual combine had dropped the NWT
    ## ranges entirely and left range names scattered across five columns, and (b) the interim
    ## GNWT-only supplement that patched the NWT hole. Reporting goes to *jurisdictional* partners,
    ## who each want their own management-unit names, so the jurisdictional sources are authoritative
    ## here; the ECCC national layer is kept as a documented comparison only (see
    ## `scripts/make_caribou_reference.R`), NOT as a reporting layer.
    "Caribou Ranges"                   , "Caribou"       , FALSE     , TRUE   , "assembly", "caribou"                                                            , NA_character_, NA_character_,
    ## NB: no `Parks` layer. Reporting is always `studyAreaReporting` INTERSECTED BY TENURE, and
    ## protected areas are carved OUT of forest tenures, so a parks layer has essentially nothing
    ## to report against here: for WesternAlbertaUpland, 81 of 732 parks merely abut the FMAs and
    ## the whole intersection is 0.05 km^2 (vs 6,538 km^2 against the buffered `studyArea`). v2
    ## never reported parks either. Parks reporting only becomes meaningful at province-wide
    ## scale, and even that does not resolve it -- provincial parks sit inside one province, but
    ## NATIONAL parks cross provincial boundaries, so summarising "all parks" needs a domain
    ## wider than any single province. Excluded pending that review (LandWeb#118).
    "Alberta Natural Subregions"       , "ANSR"          , FALSE     , TRUE   , "drive" , "1hW6zy0CpUBdk-K2IAjzW4INjVl1J4aLJ"                                   , "Name"       , NA_character_,
    "BC Biogeoclimatic zones"          , "BEC"           , FALSE     , TRUE   , "drive" , "1NS15Gd7dHEhvPOy-Ol_LBtf-4Ch6mPnS"                                   , "ZONE_NAME"  , NA_character_,
    "Northwest Territories Ecoregions" , "NT_Ecoregion"  , FALSE     , TRUE   , "drive" , "1iRAQfARkmS6-XVHFnTkB-iltzMNPAczC"                                   , "ECO4_NAM_1" , NA_character_,
    "National Ecozones"                , "Ecozone"       , FALSE     , TRUE   , "url"   , "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"     , "ZONE_NAME"  , NA_character_,
    "National Ecoregions"              , "Ecoregion"     , FALSE     , TRUE   , "url"   , "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip" , "REGION_NAM" , NA_character_
  )

  ## flags must agree across the rows that make up one layer, and exactly one layer is the tenure
  byShort <- split(out, out$NAME_SHORT)
  stopifnot(
    all(vapply(byShort, function(d) length(unique(d$isTenure)) == 1L, logical(1))),
    all(vapply(byShort, function(d) length(unique(d$cross)) == 1L, logical(1))),
    sum(vapply(byShort, function(d) d$isTenure[[1L]], logical(1))) == 1L,
    !any(out$isTenure & out$cross)
  )

  .validateLayerSources(out$source)

  ## validate the LAYER set, not the rows: repeated NAME_SHORTs are an explicit multi-source
  ## merge (see above), not the accidental collision the validator guards against.
  validateShortNames(unique(out$NAME_SHORT), what = "reportingPolygonLayers()$NAME_SHORT")
  out
}

## ---- GNWT NWT boreal caribou ------------------------------------------------------------------
## The v10 Caribou layer has no NWT ranges. The Government of the Northwest Territories publishes
## them through the ArcGIS MapServer below, which supports GeoJSON export, so the layer is fetched
## live rather than shipped as a file. Portal page (what Julie's `DataSources.xlsx` cites for
## `CaribouRrange_NWT_LCC1`): https://www.geomatics.gov.nt.ca/en/boreal-caribou-range-planning-data
##
##   * layer 96 -- "NT1 Boreal Caribou Range (GNWT 2016 version)": the NT1 range as a SINGLE
##     442,920 km^2 polygon (compiled J. Hodson & A. Smith, ENR/GNWT 2015). Same granularity as
##     the retired stopgap, so it would buy correctness/currency but not detail.
##   * layer 97 -- "Boreal Caribou Range Planning Regions": NT1 subdivided into 6 NAMED regions.
##     IN USE (see `reportingPolygonLayers()`): it subdivides NT1 the way `LOCALRANGE` /
##     `HERD_NAME` subdivide AB/BC, so an NWT tenure gets a real sub-region breakdown. Both
##     LandWeb NWT tenures (Fort Providence, Fort Resolution) fall in "Southern NWT".
##   * layer 98 -- "Other Canada Boreal Caribou Range Boundaries": the rest-of-Canada ranges, for
##     cross-checking against the v10 layer.
##
## Verified live 2026-08-12: layer 97 returns exactly 6 features, all with valid geometries, and
## `geoJSON` is in `supportedQueryFormats`. NB the service's native SR is EPSG:102002 (Canada
## Lambert Conformal Conic) but the `f=geojson` export arrives as EPSG:4326, per the GeoJSON spec
## -- immaterial here, since `prep_vector()` reprojects to `targetCRS` either way.
## Region areas from `AREA_HA` -- Southern NWT 16,241,765 ha; Sahtu 14,901,479; Wek'eezhii
## 4,950,506; Gwich'in 3,866,210; Inuvialuit 3,439,298; Yukon 892,790. The layer also carries
## disturbance attributes we do not use (`FFPpct`, `BHDpct`, `TDISTpct` over `FIRE_YR_PD`
## 1984-2023, `HDYR` 2020), which is a useful sign it is actively maintained.
##
## Open with Julie: whether NWT gets folded into the maintained LandWeb caribou layer, which would
## make this row deletable. Until then LandWeb carries it as a second, NWT-only source.
.gnwtBiologicEcologic <- paste0(
  "https://www.apps.geomatics.gov.nt.ca/arcgis/rest/services/GNWT/",
  "BiologicEcologic_LCC/MapServer"
)

#' Build the reporting-polygons named list
#'
#' Fetch each candidate reporting layer ([reportingPolygonLayers()]), clip it to
#' `studyArea` with [spatialutils::prep_vector()], and keep only those that
#' intersect -- a layer that does not overlap the study area is dropped
#' (`nrow == 0`), which replaces the former per-study-area hardcoding. Downloads
#' use the idempotent `workflowtools` `*_once` helpers (bypassing
#' `reproducible::prepInputs`). The result is a named list keyed by each layer's
#' curated `NAME_SHORT` (which is what `refCode` is slugged from, so the on-disk
#' key matches the label in the figures), suitable as the `reportingPolygons`
#' input to `NRV_summary`; each layer's `labelCols` label is copied to a common `Name`
#' column. The `"CC SAM"`/`"CC TSF"` (current condition) and `"ecoregionLayer"`
#' entries that `NRV_summary` also expects come from the simulation, not from
#' here, and are merged in separately.
#'
#' The **tenure** layer (`isTenure`) is labelled differently: its `Name` is the
#' curated tenure short name ([tenureShortNames()]) resolved from the coalesced
#' v10 member identity, not the raw `labelCols`. That both keeps the long (up to
#' 97-character) `FMA_NAME`s out of the output keys and gives the non-Alberta
#' members --- whose `FMA_NAME` is `NA` --- a label at all.
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
    dir <- reproducible::checkPath(file.path(destinationPath, .slug(lyr$key)), create = TRUE)
    ## an "assembly" layer is not one file: hand off to its assembler, which does its own fetching,
    ## per-source cleaning and merging, and returns an already-clipped layer.
    if (identical(lyr$source, "assembly")) {
      v <- .ASSEMBLERS[[lyr$id]](
        studyArea = studyArea, destinationPath = dir, targetCRS = targetCRS
      )
      return(if (is.null(v) || nrow(v) == 0L) NULL else v)
    }

    shp <- .fetchVectorFile(lyr$source, lyr$id, dir, .slug(lyr$key))
    if (length(shp) == 0L) {
      ## NB: warn rather than drop silently -- a missing file and a layer that
      ## legitimately does not intersect the study area both returned NULL, making
      ## genuine data loss indistinguishable from expected behaviour.
      warning(
        "no vector file found under '", dir, "' for reporting layer '", lyr$key,
        "'; the layer will be MISSING from the reporting polygons.",
        call. = FALSE
      )
      return(NULL)
    }
    v <- spatialutils::prep_vector(terra::vect(shp[[1]]), studyArea, crs = targetCRS)
    if (nrow(v) == 0L) {
      return(NULL)
    }
    v <- .labelReportingLayer(v, lyr)
    ## `where`: keep only the features this source is meant to contribute (a supplementary
    ## source filling a gap in the primary one -- see reportingPolygonLayers()).
    if (!is.null(lyr$where) && !is.na(lyr$where)) {
      v <- v[!is.na(v$Name) & grepl(lyr$where, v$Name), ]
    }
    if (nrow(v) == 0L) {
      return(NULL)
    }
    v
  })

  names(out) <- layers$NAME_SHORT
  out <- out[!vapply(out, is.null, logical(1))]
  .mergeLayerSources(out)
}

## Merge the built sources that make up one reporting layer (rows sharing a `NAME_SHORT`).
## Single-source layers are returned untouched, keeping all their attributes; a merged layer is
## reduced to the common `Name` column first, since the sources have unrelated schemas and `Name`
## is all the downstream summaries use.
.mergeLayerSources <- function(out) {
  if (!length(out) || !anyDuplicated(names(out))) {
    return(out)
  }
  merged <- lapply(split(out, factor(names(out), levels = unique(names(out)))), function(parts) {
    if (length(parts) == 1L) {
      return(parts[[1L]])
    }
    ## unname(): the parts share a name, which do.call() would pass as an argument name and
    ## terra's rbind() would then bind to `deparse.level`.
    do.call(rbind, unname(lapply(parts, function(v) v[, "Name"])))
  })
  merged
}

## Set the common `Name` label column on one clipped reporting layer.
## Non-tenure layers take the first non-missing `labelCols` value. The tenure layer instead resolves each feature's
## coalesced v10 identity (`.fmaMemberIdentity()`) to its curated short name, and DROPS members
## with no curated code -- an uncurated member must not silently key its outputs off a long,
## collision-prone name, so it is reported loudly and left out (add it to `tenureShortNames()`).
.labelReportingLayer <- function(v, lyr) {
  ## guard rather than `lyr$isTenure`: this is also called from `buildCaribouRanges()`, whose source
  ## table has no `isTenure` column (no caribou source is the tenure layer), and tibble's `$` warns
  ## "Unknown or uninitialised column" for a missing name.
  if ("isTenure" %in% names(lyr) && isTRUE(lyr$isTenure)) {
    ident <- .fmaMemberIdentity(as.data.frame(v))
    v$Name <- shortNameFor(ident$member)
    unknown <- unique(ident$member[is.na(v$Name) & !is.na(ident$member)])
    if (length(unknown)) {
      warning(
        "no curated short name for tenure(s): ", paste(sQuote(unknown), collapse = ", "),
        "; they are DROPPED from the reporting polygons. Add them to `tenureShortNames()`.",
        call. = FALSE
      )
    }
    return(v[!is.na(v$Name), ])
  }
  cols <- .labelColList(lyr$labelCols)
  cols <- cols[cols %in% names(v)]
  if (length(cols)) {
    ## NB: v[[col]] returns a 1-column data.frame; as.character() on that deparses
    ## the whole column into one string recycled to every row. Use `[[1]]` to pull
    ## the vector so each feature keeps its own label.
    v$Name <- do.call(.coalesceChr, lapply(cols, function(k) .blankToNA(v[[k]][[1]])))
  }
  v
}

## Split a comma-separated `labelCols` entry into its priority list of column names.
.labelColList <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) {
    return(character(0))
  }
  trimws(strsplit(as.character(x[[1L]]), ",", fixed = TRUE)[[1L]])
}

## Source placeholders that mean "no value here" and must not become a reporting-unit name.
## The Saskatchewan caribou rows, for instance, carry a literal "Not Applicable" in the
## columns that do not apply to them, which would otherwise label a real polygon.
.blankToNA <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & (trimws(x) == "" | tolower(trimws(x)) %in%
                   c("not applicable", "n/a", "na", "none", "unknown"))] <- NA_character_
  x
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
#' A landbase coverage belongs to **one tenure** (not every tenure has one), so
#' `tenure` holds that tenure's curated short name ([tenureShortNames()]):
#' [buildLandbasePolygons()] fetches a source **only** when its tenure is present
#' in the study area, avoiding the download of large (multi-GB) coverages that
#' could not intersect it anyway. `status_col` is the source attribute holding the
#' landbase status; `layer` names the layer for multi-layer File Geodatabases
#' (`NA` = the source's only/first layer). `NAME_SHORT` discriminates a tenure's
#' several coverages (West Fraser's three CLS supply blocks); it need only be
#' unique **within** a tenure, since the reporting-unit name is
#' `"<tenure> <NAME_SHORT>"`.
#'
#' @note Gating was formerly a regex (`applies_to`) matched against the
#'   *study-area name*, which silently stopped matching anything once study areas
#'   became ecoregion groups (`"WesternAlbertaUpland"` matches no company name), so
#'   no landbase was ever fetched. Gating on the tenures actually present restores it.
#'
#' @return A `tibble` with columns `key`, `NAME_SHORT`, `tenure`, `source`, `id`,
#'   `layer`, and `status_col`.
#' @export
landbaseLayers <- function() {
  out <- tibble::tribble(
    ~key                             , ~NAME_SHORT , ~tenure            , ~source , ~id                                 , ~layer            , ~status_col    ,
    "Spray Lake C5 Landbase"         , "LB_C5"     , "SprayLake"        , "drive" , "1FpMg6dJ4eblMjMkEHdcoFAYSDyo1OvTv" , "lb_20230901_tsa" , "f_active"     ,
    "West Fraser Blue Ridge Landbase", "LB"        , "BlueRidge"        , "drive" , "1Mk5L6287sKFGLY4ZfwWIUAczF5AGqAOV" , NA_character_     , "LBC_LBStatus" ,
    ## "West Fraser N" is the v2 name for the v10 Slave Lake consortium (Tolko + Vanderwell +
    ## West Fraser Mills Ltd. (Slave Lake)); S17/S20/S21 are its three CLS supply blocks.
    "West Fraser N CLS S17 Landbase" , "LB_S17"    , "Tolko_Vand_WF_SL" , "drive" , "1XrF9ygQruC2FsUulWhDxR-nD3Cd4eu7B" , NA_character_     , "F_CONDITIO"   ,
    "West Fraser N CLS S20 Landbase" , "LB_S20"    , "Tolko_Vand_WF_SL" , "drive" , "17fZw80w3n2jIKRP1-X6gq8tyOjSWh0ky" , NA_character_     , "F_CONDITIO"   ,
    "West Fraser N CLS S21 Landbase" , "LB_S21"    , "Tolko_Vand_WF_SL" , "drive" , "1akMUL-lRumTfmmWG7WF7-9KDrFxTPN3Z" , NA_character_     , "F_CONDITIO"   ,
    "Tolko AB North Landbase"        , "LB"        , "Tolko_Norbord_LC" , "drive" , "1qzWhMR-nP_2N2KLL84ZHCVGNNlWxd2bZ" , NA_character_     , "LBC_Landbase" ,
    "Sundre FP Landbase"             , "LB"        , "Sundre"           , "drive" , "1oh_w9nALKufCQXb1PR3VIlChPHPysHF4" , NA_character_     , "LBC_LBStatus" ,
    "Manning Landbase"               , "LB"        , "Manning"          , "drive" , "1lY0p6Ms84paja9p1lmGXz5jaCgv2_VkY" , NA_character_     , "LBC_LBStat"
  )

  ## unique WITHIN a tenure (the reporting-unit name carries the tenure prefix), and every
  ## `tenure` must name a curated tenure -- a typo here would silently disable that landbase.
  validateShortNames(out$NAME_SHORT, what = "landbaseLayers()$NAME_SHORT", by = out$tenure)
  unknown <- setdiff(out$tenure, tenureShortNames()$NAME_SHORT)
  if (length(unknown)) {
    stop(
      "landbaseLayers(): unknown tenure(s) ", paste(sQuote(unknown), collapse = ", "),
      "; must be a `tenureShortNames()$NAME_SHORT`.", call. = FALSE
    )
  }
  out
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
#' For each candidate landbase source ([landbaseLayers()]) whose `tenure` is
#' present in `tenures`, download it (idempotent `workflowtools` `*_once`
#' helpers), then dissolve it by landbase status and clip it to `studyArea` with
#' [spatialutils::prep_landbase()] -- which pushes the column and spatial
#' filters down to the read, so a large coverage is reduced to the study area
#' before any geometry work. Sources that do not intersect the study area (or
#' whose tenure is absent) are skipped; the result is a named list of
#' `SpatVector`s keyed `"<tenure> <NAME_SHORT>"`, each with a `Name` column
#' holding the landbase status, suitable for merging into the
#' `reportingPolygons` input to `NRV_summary`.
#'
#' @param studyArea A `SpatVector` (or a source readable by [terra::vect()]).
#' @param tenures Character vector of the curated tenure short names present in
#'   this study area (e.g. `reportingPolygons[["FMA"]]$Name`, or the study-area
#'   group's `name_short`s from [studyAreaCrosswalk()]). Only landbases belonging
#'   to one of these are fetched.
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
  tenures,
  destinationPath,
  targetCRS = LandWebCRS,
  layers = landbaseLayers()
) {
  if (!inherits(studyArea, "SpatVector")) {
    studyArea <- terra::vect(studyArea)
  }

  ## gate: only consider landbases belonging to a tenure present in this study area
  layers <- layers[layers$tenure %in% as.character(tenures), , drop = FALSE]
  if (nrow(layers) == 0L) {
    return(list())
  }

  out <- lapply(seq_len(nrow(layers)), function(i) {
    lyr <- layers[i, ]
    dir <- reproducible::checkPath(file.path(destinationPath, .slug(lyr$key)), create = TRUE)
    dest <- file.path(dir, paste0(.slug(lyr$key), ".zip"))
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

  names(out) <- paste(layers$tenure, layers$NAME_SHORT)
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) list() else out ## a plain empty list, matching the no-gate early return
}
