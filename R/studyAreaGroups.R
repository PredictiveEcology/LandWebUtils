utils::globalVariables(c("fma_name", "eco_unit", "province", "company", "group"))

## Element-wise coalesce of character vectors: first non-missing, non-empty value.
.coalesceChr <- function(...) {
  vs <- lapply(list(...), as.character)
  out <- rep(NA_character_, length(vs[[1L]]))
  for (v in vs) {
    take <- is.na(out) & !is.na(v) & nzchar(v)
    out[take] <- v[take]
  }
  out
}

## Turn an ecological-unit name ("Western Alberta Upland") into a study-area group
## token ("WesternAlbertaUpland") safe for use as a target name / output directory.
## Preserves the source casing (the eco layer's names are already title-cased); only
## strips punctuation and whitespace.
.groupToken <- function(x) {
  gsub("\\s+", "", trimws(gsub("[^A-Za-z0-9 ]", "", x)))
}

## sf-friendly adapter for spatialutils::repair_geoms() (terra::makeValid on just the
## invalid subset). Used instead of sf::st_make_valid(), which COLLAPSES polygons with
## reversed ring winding order (present in the v10 FMA layer, e.g. Prince Albert:
## -31,600 km^2 -> a 4 km^2 sliver) that terra correctly re-orients. The repair logic
## lives in spatialutils; this only adapts sf <-> SpatVector.
.repairGeom <- function(x) {
  sf::st_as_sf(spatialutils::repair_geoms(x))
}

## The v10 combined-FMA layer identifies areas in different fields per jurisdiction
## (AB/NT: FMA_NAME; ON: FMU_NAME; BC: TSA_NUMB_1; SK: FOREST_NAM; MB: FML_NAME).
## Coalesce them into a single stable member identity, and derive the province.
.fmaMemberIdentity <- function(d) {
  member <- .coalesceChr(
    d[["FMA_NAME"]], d[["FMU_NAME"]], d[["TSA_NUMB_1"]],
    d[["FOREST_NAM"]], d[["FML_NAME"]], d[["Name"]]
  )
  province <- ifelse(!is.na(d[["TSA_NUMB_1"]]), "BC",
    ifelse(!is.na(d[["FOREST_NAM"]]), "SK",
      ifelse(!is.na(d[["FML_NAME"]]), "MB",
        ifelse(!is.na(d[["FMU_NAME"]]) & is.na(d[["FMA_NAME"]]), "ON", "AB"))))
  ## NWT areas carry an FMA_NAME (so default to "AB"); re-tag by name.
  province[grepl("Fort Providence|Fort Resolution", member)] <- "NT"
  company <- .coalesceChr(d[["LICENSEE_N"]], member)
  data.frame(member = member, province = province, company = company, stringsAsFactors = FALSE)
}

#' Build the LandWeb study-area grouping crosswalk
#'
#' Groups the LandWeb FMA/TSA/FML boundary polygons into ecologically-coherent
#' study areas by assigning each to its **dominant ecological unit** (an ecoregion by
#' default) via spatial overlap, then names the group after that unit. Small,
#' scattered forest tenures (e.g. a single small FMA) are thereby run *inside* a
#' larger, ecologically homogeneous area, so the downstream `studyAreaANPP` PSP
#' catchment is both large enough for trait estimation and ecologically appropriate
#' (species traits can differ across ecoregions, so PSP pools should not span
#' drastically different ecology).
#'
#' The grouping unit is a parameter: pass a coarser layer (ecoprovinces) for fewer,
#' larger groups, or a finer one (ecodistricts) to subdivide groups that prove too
#' large computationally. Membership is by **verified spatial overlap against the
#' current FMA layer** --- never by company/tenure labels, which drift over time (a
#' polygon labelled for one licensee in an older file may belong to another now).
#'
#' @param fmas `sf` of the combined LandWeb FMA boundaries (e.g. [prepFMAs()]).
#' @param eco `sf` of ecological units to group by (e.g. an ecoregions layer). Supply
#'   an ecoprovince or ecodistrict layer to change the grouping granularity.
#' @param eco_field character naming the column in `eco` that holds each unit's name
#'   (default `"REGION_NAM"`, the ecoregion name field).
#' @param res_m numeric simulation pixel size in metres, used only for the
#'   `mpix` size estimate (default `240`).
#' @param min_area_km2 numeric; member tenures smaller than this are dropped as
#'   sliver artifacts (with a message listing them), since they cannot form a
#'   meaningful study area. Default `100`. Set to `0` to keep everything.
#'
#' @return A `data.frame` crosswalk, one row per member tenure, with columns:
#'   `group` (study-area group token), `fma_name` (the member's coalesced
#'   FMA/TSA/FML identity --- the join key back to `fmas`), `company` (licensee, for
#'   company-driven reruns), `province`, `eco_unit` (the dominant ecological unit),
#'   `area_km2`, and `mpix` (approx. pixels at `res_m`). Sort/split by `group` to get
#'   each study area's member set.
#'
#' @export
build_studyarea_crosswalk <- function(fmas, eco, eco_field = "REGION_NAM", res_m = 240,
                                      min_area_km2 = 100) {
  stopifnot(inherits(fmas, "sf"), inherits(eco, "sf"))
  if (!eco_field %in% names(eco)) {
    stop("`eco_field` '", eco_field, "' is not a column of `eco`.")
  }
  fmas <- .repairGeom(fmas) # fix reversed-winding polygons before any area/overlap
  ident <- .fmaMemberIdentity(sf::st_drop_geometry(fmas))
  fmas$fma_name <- ident$member
  fmas$province <- ident$province
  fmas$company <- ident$company
  fmas <- fmas[!is.na(fmas$fma_name), ]

  ## dissolve multi-part tenures to one geometry per member (removes sliver-fragment
  ## mis-assignments, e.g. a 4 km^2 outlier polygon landing in the wrong ecoregion).
  fmad <- stats::aggregate(
    fmas[, c("province", "company")],
    by = list(fma_name = fmas$fma_name),
    FUN = function(x) x[[1L]]
  )
  fmad <- .repairGeom(fmad)
  eco <- .repairGeom(sf::st_transform(eco, sf::st_crs(fmad)))

  dom <- vapply(seq_len(nrow(fmad)), function(i) {
    g <- sf::st_geometry(fmad[i, ])
    hit <- which(sf::st_intersects(eco, g, sparse = FALSE)[, 1L])
    if (!length(hit)) {
      return(NA_character_)
    }
    a <- vapply(hit, function(j) {
      sum(as.numeric(sf::st_area(sf::st_intersection(sf::st_geometry(eco[j, ]), g))))
    }, numeric(1))
    as.character(eco[[eco_field]][hit[which.max(a)]])
  }, character(1))

  pxkm2 <- (res_m^2) / 1e6
  out <- data.frame(
    group = .groupToken(dom),
    fma_name = fmad$fma_name,
    company = fmad$company,
    province = fmad$province,
    eco_unit = dom,
    area_km2 = round(as.numeric(sf::st_area(fmad)) / 1e6),
    stringsAsFactors = FALSE
  )
  out$mpix <- round(out$area_km2 / pxkm2 / 1e6, 3)

  tiny <- out[out$area_km2 < min_area_km2, , drop = FALSE]
  if (nrow(tiny)) {
    message(
      "Dropping ", nrow(tiny), " sub-threshold tenure(s) (< ", min_area_km2,
      " km^2, likely sliver artifacts): ",
      paste(sprintf("%s [%s, %g km2]", tiny$fma_name, tiny$province, tiny$area_km2), collapse = "; ")
    )
    out <- out[out$area_km2 >= min_area_km2, , drop = FALSE]
  }

  out <- out[order(out$group, out$fma_name), ]
  rownames(out) <- NULL
  out
}

#' Build (or load a cached) LandWeb study-area grouping crosswalk
#'
#' Convenience wrapper that assembles [build_studyarea_crosswalk()]'s inputs from
#' the packaged prep helpers --- [prepFMAs()] for the FMA boundaries and, by
#' default, [prepEcoregions()] for the grouping layer --- and caches the result
#' under `destinationPath` so repeated study-area setups within a run do not rebuild
#' it (the overlay is the expensive part).
#'
#' @inheritParams prepStudyArea
#' @param eco a function of `(destinationPath, targetCRS)` returning the ecological
#'   grouping layer as `sf` (default [prepEcoregions]); pass [prepEcoprovinces] for
#'   coarser groups or a custom prep for a different eco\* layer.
#' @param eco_field character naming the grouping-unit column in that layer
#'   (default `"REGION_NAM"`, the ecoregion name).
#' @param cache logical; reuse/save `file.path(destinationPath, "studyAreaCrosswalk.rds")`.
#'
#' @return the crosswalk `data.frame` (see [build_studyarea_crosswalk()]).
#' @export
studyAreaCrosswalk <- function(destinationPath, targetCRS = LandWebCRS,
                               eco = prepEcoregions, eco_field = "REGION_NAM",
                               cache = TRUE) {
  cachePath <- file.path(destinationPath, "studyAreaCrosswalk.rds")
  if (cache && file.exists(cachePath)) {
    return(readRDS(cachePath))
  }
  fmas <- prepFMAs(destinationPath, targetCRS)
  ecoLayer <- eco(destinationPath, targetCRS)
  cw <- build_studyarea_crosswalk(fmas, ecoLayer, eco_field = eco_field)
  if (cache) {
    saveRDS(cw, cachePath)
  }
  cw
}

## Union the member FMA polygons of one crosswalk `group` into a single study-area
## `sf` boundary. Members are matched to `fmas` by the same coalesced identity that
## built the crosswalk.
.extractStudyAreaGroup <- function(fmas, group, crosswalk) {
  members <- crosswalk[["fma_name"]][crosswalk[["group"]] == group]
  fmas <- .repairGeom(fmas)
  ident <- .fmaMemberIdentity(sf::st_drop_geometry(fmas))
  sel <- fmas[ident$member %in% members, ]
  if (!nrow(sel)) {
    stop("no FMA polygons matched study-area group '", group, "'.")
  }
  sf::st_sf(studyAreaName = group, geometry = sf::st_union(sf::st_geometry(sel)))
}
