## Growth-curve traits fitted once for all of LandWeb (Biomass_speciesParameters run on every species in
## a shared pipeline stage), and given to each study area's simulated units.

growthTraitCols <- c(
  "growthcurve",
  "mortalityshape",
  "mANPPproportion",
  "inflationFactor",
  "longevity"
)

#' Species table for the shared growth-curve fit
#'
#' [landweb_species_sppEquiv()] with the hybrid spruce counted as Engelmann spruce, so that
#' `Biomass_speciesParameters` relabels BC provincial plot records of white spruce outside the boreal
#' BEC zones as the hybrid (its `relabelHybridSprucePSP()` requires the hybrid row's `LandR` code to be
#' `Pice_eng`). The hybrid row also takes Engelmann spruce's `LANDIS_traits` code, so that Engelmann
#' spruce's traits do not change by merging it.
#'
#' @param sppEquiv `data.table` from [landweb_species_sppEquiv()].
#'
#' @return A copy of `sppEquiv`.
#'
#' @export
landweb_growth_sppEquiv <- function(sppEquiv) {
  eq <- data.table::copy(data.table::as.data.table(sppEquiv))
  hybrid <- which(eq[["LandR"]] %in% "Pice_eng_gla")
  parent <- which(eq[["LandR"]] %in% "Pice_eng")
  if (length(hybrid) && length(parent)) {
    data.table::set(eq, hybrid, c("LandR", "LandWeb"), list("Pice_eng", "Pice_eng"))
    if ("LANDIS_traits" %in% names(eq)) {
      data.table::set(eq, hybrid, "LANDIS_traits", eq[["LANDIS_traits"]][parent[1L]])
    }
  }
  eq
}

#' The growth-curve traits of each species, from the shared fit
#'
#' @param species the shared growth-curve stage's species table: one row per species, after
#'   `Biomass_speciesParameters` has fitted the species it can and given the others their hardwood or
#'   softwood means (`growthCurveSource` `"imputed"`).
#'
#' @return `data.table` with `species`, the five traits a growth curve sets (`growthcurve`,
#'   `mortalityshape`, `mANPPproportion`, `inflationFactor`, `longevity`) and `growthCurveSource`.
#'
#' @export
landweb_growth_traits <- function(species) {
  sp <- data.table::as.data.table(species)
  cols <- c("species", growthTraitCols, "growthCurveSource")
  miss <- setdiff(cols, names(sp))
  if (length(miss)) {
    stop("the species table has no ", paste(miss, collapse = ", "), call. = FALSE)
  }
  out <- sp[, cols, with = FALSE]
  if (anyNA(out[, growthTraitCols, with = FALSE])) {
    stop(
      "growth-curve traits are missing for: ",
      paste(
        out[["species"]][!stats::complete.cases(out[, growthTraitCols, with = FALSE])],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  out[]
}

#' Give each simulated unit the growth-curve traits of a member species
#'
#' Each unit of a study area takes, from the shared fit, the traits of its member with most cover among
#' those whose curve was fitted (`growthCurveSource` `"estimated"`); a unit none of whose members was
#' fitted takes its dominant member's traits, which are then hardwood or softwood means. Members follow
#' the unit's `ranked` order (species in `neverSplit` last); without it, the dominant member is used. A
#' merged fir unit dominated by unfitted balsam fir thus takes subalpine fir's fitted curve rather than
#' the conifer mean. The traits' source is recorded in `growthTraitSource`, e.g.
#' `"estimated (Abie_las)"`. Apply the result to `speciesEcoregion` with
#' `LandR::modifySpeciesAndSpeciesEcoregionTable()`.
#'
#' @param species a study area's species table (one row per unit), from data preparation.
#' @param traits from [landweb_growth_traits()].
#' @param units the study area's units, from [landweb_species_units()].
#'
#' @return A copy of `species` with the traits set.
#'
#' @export
landweb_unit_growth_traits <- function(species, traits, units) {
  sp <- data.table::copy(data.table::as.data.table(species))
  tr <- data.table::as.data.table(traits)
  dominantOf <- stats::setNames(units[["dominant"]], units[["unit"]])
  ranked <- if ("ranked" %in% names(units)) {
    stats::setNames(units[["ranked"]], units[["unit"]])
  } else {
    stats::setNames(as.list(units[["dominant"]]), units[["unit"]])
  }
  fitted <- tr[["species"]][tr[["growthCurveSource"]] %in% "estimated"]
  src <- vapply(
    as.character(sp[["species"]]),
    function(u) {
      firstFitted <- intersect(ranked[[u]], fitted)
      if (length(firstFitted)) firstFitted[[1L]] else unname(dominantOf[u])
    },
    character(1),
    USE.NAMES = FALSE
  )
  if (anyNA(src)) {
    stop(
      "no dominant species for unit(s): ",
      paste(sp[["species"]][is.na(src)], collapse = ", "),
      call. = FALSE
    )
  }
  noTraits <- setdiff(src, tr[["species"]])
  if (length(noTraits)) {
    stop("no growth-curve traits for: ", paste(noTraits, collapse = ", "), call. = FALSE)
  }
  m <- tr[match(src, tr[["species"]])]
  for (cl in growthTraitCols) {
    v <- m[[cl]]
    if (cl %in% c("mortalityshape", "longevity")) {
      v <- as.integer(round(v))
    } else {
      v <- as.numeric(v)
    }
    data.table::set(sp, j = cl, value = v)
  }
  data.table::set(
    sp,
    j = "growthTraitSource",
    value = paste0(m[["growthCurveSource"]], " (", src, ")")
  )
  sp[]
}

#' BC BEC zones over the BC plots of a growth-curve fitting area
#'
#' Fetches the biogeoclimatic (BEC) zone polygons from the BC Data Catalogue over the BC provincial plots
#' inside `area`, for `Biomass_speciesParameters`' hybrid spruce relabelling (its input `BECzonesBC`).
#' The module fetches them itself when they are not supplied, but a failed fetch there only warns and
#' skips the relabelling, which changes the spruce curves; here it stops.
#'
#' @param PSPgis `sf` of plot locations, with `OrigPlotID1` (from `PSPclean::getPSP()`).
#' @param PSPmeasure `data.table` of tree measurements, with `OrigPlotID1` and `source`.
#' @param area the fitting area (`sf` or `sfc`).
#'
#' @return `sf` of BEC polygons with a `ZONE` column, or `NULL` if no BC plot lies in `area`.
#'
#' @export
landweb_bec_zones <- function(PSPgis, PSPmeasure, area) {
  if (!requireNamespace("bcdata", quietly = TRUE)) {
    stop("landweb_bec_zones() needs the bcdata package", call. = FALSE)
  }
  bcPlots <- unique(PSPmeasure[["OrigPlotID1"]][PSPmeasure[["source"]] %in% "BC"])
  pts <- PSPgis[PSPgis[["OrigPlotID1"]] %in% bcPlots, ]
  if (nrow(pts) > 0L) {
    pts <- pts[lengths(sf::st_intersects(sf::st_transform(pts, sf::st_crs(area)), area)) > 0L, ]
  }
  if (nrow(pts) == 0L) {
    return(NULL)
  }
  ## the box is built before the query: bcdata cannot translate a geometry built inside filter()
  box <- sf::st_as_sfc(sf::st_bbox(sf::st_transform(pts, 3005)))
  bcdata::bcdc_query_geodata("WHSE_FOREST_VEGETATION.BEC_BIOGEOCLIMATIC_POLY") |>
    bcdata::filter(bcdata::INTERSECTS(box)) |>
    bcdata::select("ZONE") |>
    bcdata::collect()
}
