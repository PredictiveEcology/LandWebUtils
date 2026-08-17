#' Jurisdictional caribou range source layers
#'
#' The six provincial/territorial sources assembled into LandWeb's caribou reporting
#' layer by [buildCaribouRanges()]. LandWeb reports to **jurisdictional** partners, who
#' each want their own management-unit names, so these are authoritative for reporting
#' rather than any national layer.
#'
#' This replaces two earlier approaches. Julie Duval's pre-combined v10 layer
#' (`CaribouRange_Regional_ltfc_v10_LCC1`) was assembled by hand, and the combine had
#' **dropped the Northwest Territories ranges** --- 0 of its 74 features intersect either
#' NWT tenure --- while leaving range names scattered across five per-jurisdiction
#' columns. A GNWT-only supplement then patched the NWT hole. Assembling from the
#' sources here removes both problems and the whole class of defect they came from: a bad
#' merge becomes a one-row fix instead of a re-delivery.
#'
#' `source` uses the [reportingPolygonLayers()] transports (`"drive"`, `"url"`,
#' `"geojson"`). `labelCols` is a comma-separated priority list coalesced left to right,
#' as there. `statusCols` names the candidate column(s) carrying herd status and
#' `extirpated` the regular expression matching values to **exclude**; both `NA` where a
#' source ships no status field.
#'
#' @section Why the ECCC national layer is not here:
#' Measured against the assembled layer, ECCC adds no coverage: with designatable-unit
#' scope matched (boreal + mountain), the two footprints agree to **0.7%** --- the
#' assembly divides the same ground into 74 units where the federal layer has 52. The
#' difference is granularity and naming, not extent. ECCC's Southern Mountain units are
#' aggregate Local Population Units that *merge* the finer provincial herds, and for
#' ranges wholly inside one province the provincial polygon is geometrically identical to
#' ECCC's (AB IoU = 1.0000 for 9 of 12 boreal ranges). It is retained as a documented
#' comparison in `scripts/make_caribou_reference.R`, not as a reporting layer.
#'
#' @section Per-source notes:
#' * **AB** --- labelled on `SUBUNIT` before `LOCALRANGE`. v2 reported AB's mountain
#'   contribution at subunit granularity: `A La Peche`, `Narraway` and
#'   `Redrock-Prairie Creek` were three separate v2 units but are all subunits of the one
#'   `LOCALRANGE` `"West Central"`, so labelling on `LOCALRANGE` alone would collapse
#'   them. NB `Bischto` is Alberta's own misspelling of `Bistcho`; normalised by
#'   [.caribouNameFixes()].
#' * **BC** --- fetched live from the BC Data Catalogue WFS rather than a cached export.
#'   The WFS release is **newer and materially different** from the v2-era shapefile: 55
#'   features against 72, and it flags **5** extirpated herds where the old export flagged
#'   1. NB the WFS uses unabbreviated field names (`HERD_STATUS`, `RISK_STATUS`) while a
#'   shapefile export truncates them to 10 characters (`HERD_STAT`), hence the candidate
#'   list in `statusCols`.
#' * **SK** --- `CONUNIT` first, falling back to `RNGEUNIT`: SK1 carries the literal
#'   `"Not Applicable"` in `CONUNIT`, which [buildReportingPolygons()]'s placeholder
#'   handling treats as missing. Same footprint as ECCC (IoU 0.9983), subdivided finer.
#' * **MB** --- by-request from the Government of Manitoba, so Drive-hosted; there is no
#'   public endpoint to refresh from. This 2015 re-delineation is **newer than any ECCC
#'   boundary set** and 11.9% of it falls outside every ECCC range, so it is not
#'   substitutable. Worth re-requesting a current copy.
#' * **NWT** --- GNWT MapServer layer 97, which subdivides NT1 into 6 named planning
#'   regions the way `LOCALRANGE`/`HERD_NAME` subdivide AB/BC. Labelled on `REGION`, not
#'   `NAME`: both exist, but `NAME` carries diacritics and apostrophes (`Sahtú`,
#'   `Wek'èezhìi`, `Gwich'in`) while `REGION` is ASCII-safe, and this label is slugged into
#'   the `refCode` that keys parquet aggregates and figure filenames.
#' * **ON** --- the LIO direct download. Required, and easy to overlook: `RANGE_NAME` is
#'   populated for exactly the ON and MB features of the old combined layer, so omitting
#'   Ontario silently loses every ON tenure x Caribou reporting unit.
#'
#' @return A `tibble` with columns `juris`, `key`, `source`, `id`, `labelCols`,
#'   `statusCols` and `extirpated`.
#' @seealso [buildCaribouRanges()], [reportingPolygonLayers()]
#' @export
caribouRangeLayers <- function() {
  out <- tibble::tribble(
    ~juris , ~key                              , ~source   , ~id                                                                                      , ~labelCols          , ~statusCols            , ~extirpated,
    "AB"   , "Alberta caribou ranges"          , "url"     , "https://extranet.gov.ab.ca/srd/geodiscover/srd_pub/LAT/FWDSensitivity/CaribouRange.zip"  , "SUBUNIT,LOCALRANGE", "STATUS"               , "^Extirp",
    "BC"   , "BC caribou herd locations"       , "geojson" , .bcCaribouWFS                                                                             , "HERD_NAME"         , "HERD_STATUS,HERD_STAT", "^Extirpated$",
    "SK"   , "Saskatchewan caribou ranges"     , "geojson" , paste0("https://geohub.saskatchewan.ca/api/download/v1/items/",
                                                                    "61c26278ee28415abc9bce5018e6ec10/geojson?layers=0")                               , "CONUNIT,RNGEUNIT"  , NA_character_          , NA_character_,
    "MB"   , "Manitoba caribou ranges"         , "drive"   , "1Y_Qi3twoU3fHaNgMzF5QEl1CosGmGyha"                                                       , "RANGE_NAME"        , NA_character_          , NA_character_,
    "NWT"  , "Northwest Territories caribou"   , "geojson" , paste0(.gnwtBiologicEcologic,
                                                                    "/97/query?where=1%3D1&outFields=REGION,NAME&returnGeometry=true&f=geojson")       , "REGION"            , NA_character_          , NA_character_,
    "ON"   , "Ontario caribou range boundary"  , "url"     , "https://ws.gisetl.lrc.gov.on.ca/fmedatadownload/Packages/CARIBRNG.zip"                    , "RANGE_NAME,RANGE"  , NA_character_          , NA_character_
  )
  .validateLayerSources(out$source)
  stopifnot(
    "caribouRangeLayers(): `juris` must be unique" = !anyDuplicated(out$juris),
    "caribouRangeLayers(): `key` must be unique" = !anyDuplicated(out$key)
  )
  out
}

## BC Data Catalogue WFS for GCPB_CARIBOU_POPULATION_SP. Preferred over the cached order-based
## BCGW export (`BCGW_<order>_...zip`), whose URL is not reproducible and whose content is stale.
.bcCaribouWFS <- paste0(
  "https://openmaps.gov.bc.ca/geo/pub/WHSE_WILDLIFE_INVENTORY.GCPB_CARIBOU_POPULATION_SP/ows",
  "?service=WFS&version=2.0.0&request=GetFeature",
  "&typeName=pub:WHSE_WILDLIFE_INVENTORY.GCPB_CARIBOU_POPULATION_SP",
  "&outputFormat=application/json"
)

## Source spellings that disagree across jurisdictions for the SAME range. Left unfixed these go
## straight into partner-facing figures, and they also block the name-grouping that stitches a
## cross-border range back into one reporting unit (the summaries group by NAME, not by feature) --
## which is how `Chinchaga` and `Narraway` already come out whole across the AB/BC border.
## Keyed by the wrong spelling; values are the accepted form.
.caribouNameFixes <- function() {
  c(
    "Bischto" = "Bistcho",          ## Alberta's own typo
    "Snake-Sahtaneh" = "Snake-Sahtahneh" ## BC drops the second 'h'
  )
}

## Boreal vs mountain caribou, per source. Only BC ships an explicit field.
##   * BC `ECOTYPE` -- "Northern" and "Mountain" are BOTH mountain caribou (the Northern Mountain and
##     Southern Mountain designatable units); the split that matters here is boreal vs mountain.
##   * AB -- `Banff`, `Jasper` and `West Central` are its mountain ranges (`West Central`'s subunits
##     being A La Peche, Narraway and Redrock-Prairie Creek); every other LOCALRANGE is boreal.
##   * SK, MB, NWT, ON -- boreal caribou only within this AOI.
## Barren-ground and Peary caribou are out of scope for LandWeb and are not represented by any of
## these sources (the jurisdictions publish barren-ground separately -- GNWT has its own layer, and
## BC's ECOTYPE has no barren-ground class).
.caribouEcotype <- function(x, juris) {
  nm <- names(x)
  if (identical(juris, "BC") && "ECOTYPE" %in% nm) {
    e <- trimws(as.character(x$ECOTYPE))
    return(ifelse(!is.na(e) & e == "Boreal", "Boreal", "Mountain"))
  }
  if (identical(juris, "AB") && "LOCALRANGE" %in% nm) {
    lr <- trimws(as.character(x$LOCALRANGE))
    return(ifelse(lr %in% c("Banff", "Jasper", "West Central"), "Mountain", "Boreal"))
  }
  rep("Boreal", nrow(x))
}

## Drop locally extirpated herds. `statusCols` is a candidate list because a source's field name
## differs between its WFS (`HERD_STATUS`) and a shapefile export truncated to 10 chars
## (`HERD_STAT`); the first column present wins.
.dropExtirpated <- function(x, statusCols, extirpated) {
  if (is.na(statusCols) || is.na(extirpated)) {
    return(x)
  }
  cols <- .labelColList(statusCols)
  cols <- cols[cols %in% names(x)]
  if (!length(cols)) {
    return(x)
  }
  v <- trimws(as.character(x[[cols[[1L]]]][[1L]]))
  x[is.na(v) | !grepl(extirpated, v), ]
}

#' Assemble the caribou range reporting layer
#'
#' Fetch each jurisdictional source in [caribouRangeLayers()], drop locally extirpated
#' herds, label it from its `labelCols`, tag its `juris` and boreal/mountain `ecotype`,
#' clip it to `studyArea` with [spatialutils::prep_vector()], and bind the intersecting
#' pieces into one layer. Sources that do not reach the study area are dropped, so this
#' works unchanged for a study area in any jurisdiction.
#'
#' Range names are normalised via an explicit spelling table (see the source of
#' [caribouRangeLayers()]): a range that two provinces spell differently would otherwise
#' appear as two reporting units, since the summaries group by **name** rather than by
#' feature. That same name-grouping is what stitches a cross-border range whose two
#' provincial shares carry the same name --- `Chinchaga` and `Narraway` --- back into a
#' single reporting unit.
#'
#' @param studyArea A `SpatVector` (or anything [terra::vect()] accepts).
#' @param destinationPath Directory for downloads/extraction.
#' @param targetCRS Target CRS (default [LandWebCRS]).
#' @param layers Candidate-source table (default [caribouRangeLayers()]).
#'
#' @return A `SpatVector` with columns `Name`, `juris` and `ecotype`, or `NULL` if no
#'   source intersects `studyArea`.
#' @seealso [caribouRangeLayers()], [buildReportingPolygons()]
#' @export
buildCaribouRanges <- function(
  studyArea,
  destinationPath,
  targetCRS = LandWebCRS,
  layers = caribouRangeLayers()
) {
  if (!inherits(studyArea, "SpatVector")) {
    studyArea <- terra::vect(studyArea)
  }
  fixes <- .caribouNameFixes()

  parts <- lapply(seq_len(nrow(layers)), function(i) {
    lyr <- layers[i, ]
    dir <- reproducible::checkPath(
      file.path(destinationPath, paste0("caribou_", lyr$juris)), create = TRUE
    )
    f <- .fetchVectorFile(lyr$source, lyr$id, dir, paste0("caribou_", lyr$juris))
    if (length(f) == 0L) {
      ## warn rather than drop silently: a missing file and a source that legitimately does not
      ## intersect the study area both return NULL, making genuine data loss invisible.
      warning(
        "no vector file found under '", dir, "' for caribou source '", lyr$key,
        "'; ", lyr$juris, " will be MISSING from the caribou ranges.",
        call. = FALSE
      )
      return(NULL)
    }
    v <- terra::vect(f[[1L]])
    v <- .dropExtirpated(v, lyr$statusCols, lyr$extirpated)
    if (nrow(v) == 0L) {
      return(NULL)
    }
    v$ecotype <- .caribouEcotype(v, lyr$juris)
    v <- .labelReportingLayer(v, lyr)
    v <- v[!is.na(v$Name), ]
    if (nrow(v) == 0L) {
      return(NULL)
    }
    hit <- v$Name %in% names(fixes)
    if (any(hit)) {
      v$Name[hit] <- unname(fixes[v$Name[hit]])
    }
    v$juris <- lyr$juris
    v <- spatialutils::prep_vector(v[, c("juris", "Name", "ecotype")], studyArea, crs = targetCRS)
    if (nrow(v) == 0L) NULL else v
  })

  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) {
    return(NULL)
  }
  do.call(rbind, parts)
}

## Layers assembled from several sources rather than fetched from one file, keyed by the `id` of
## their `source == "assembly"` row in `reportingPolygonLayers()`. Each takes
## (studyArea, destinationPath, targetCRS) and returns an already-clipped SpatVector or NULL.
.ASSEMBLERS <- list(caribou = function(studyArea, destinationPath, targetCRS) {
  buildCaribouRanges(studyArea, destinationPath = destinationPath, targetCRS = targetCRS)
})
