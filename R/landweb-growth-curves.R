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

#' Damage agent codes to exclude from growth-curve plots, by plot source
#'
#' LandWeb simulates neither bark beetles nor defoliators, so a growth curve fitted on trees they
#' damaged would read their damage as growth and senescence. BC and the NFI code damage agents alike: a
#' leading `I` for insects, then `B` for bark beetles and `D` for defoliators; BC lists each species
#' (`IBM` mountain pine beetle, `IDE` spruce budworm, ...) and the budworm genus as `CHX`. The BC codes
#' are those of BC's damage agent table as PSPclean reads it (`BCForestry_DamageAgentCodes.csv`).
#' Alberta records damage causes as numbers (1 spruce budworm, 2 defoliator, 3 mountain pine beetle);
#' Saskatchewan records only the cause of a tree's death, with all insects as 3, so either class
#' excludes every tree that insects killed there.
#'
#' @param agents damage agent classes: any of `"barkBeetles"` and `"defoliators"`.
#'
#' @return list named by plot source (`BC`, `AB`, `SK`, `NFI`), for `PSPclean::getPSP()`'s
#'   `codesToExclude`.
#'
#' @export
landweb_damage_codes <- function(agents = c("barkBeetles", "defoliators")) {
  agents <- match.arg(agents, several.ok = TRUE)
  codes <- list(
    barkBeetles = list(
      BC = c("IB", "IBB", "IBD", "IBI", "IBM", "IBP", "IBS", "IBT", "IBW"),
      AB = 3L,
      SK = 3L,
      NFI = "IB"
    ),
    defoliators = list(
      BC = c(
        "ID",
        "IDA",
        "IDB",
        "IDC",
        "IDD",
        "IDE",
        "IDF",
        "IDG",
        "IDH",
        "IDI",
        "IDL",
        "IDM",
        "IDN",
        "IDP",
        "IDR",
        "IDS",
        "IDT",
        "IDU",
        "IDV",
        "IDW",
        "IDX",
        "IDZ",
        "CHX"
      ),
      AB = c(1L, 2L),
      SK = 3L,
      NFI = "ID"
    )
  )[agents]
  sources <- c("BC", "AB", "SK", "NFI")
  stats::setNames(
    lapply(sources, function(src) unique(unlist(lapply(codes, `[[`, src), use.names = FALSE))),
    sources
  )
}

#' Resample the plots of a growth-curve fit, whole plots with replacement
#'
#' A bootstrap sample of the plots inside the fitting area, for refitting the shared growth curves to
#' see how much the fitted traits depend on which plots happen to have been measured. As many plots as
#' lie inside `area` are drawn with replacement, each with all its measurements and trees. A plot drawn
#' more than once is copied, and each copy's plot and measurement IDs get a suffix (`_r1`, `_r2`, ...),
#' so that `Biomass_speciesParameters` counts the copies as separate plots. `area` selects plots as the
#' module does (`PSPgis[studyAreaANPP, ]`), so plots the fit would not use are not drawn.
#'
#' @param psp list with `PSPmeasure`, `PSPplot` and `PSPgis`, from `PSPclean::getPSP()`.
#' @param area the fitting area (`sf` or `sfc`), given to `Biomass_speciesParameters` as
#'   `studyAreaANPP`.
#' @param seed integer seed for the draw. The global random number generator is left as it was.
#'
#' @return `psp`, with the drawn plots only.
#'
#' @export
landweb_resample_psp <- function(psp, area, seed) {
  gis <- psp[["PSPgis"]]
  area <- sf::st_transform(sf::st_as_sf(area), sf::st_crs(gis))
  inside <- unique(gis[area, ][["OrigPlotID1"]])
  if (length(inside) == 0L) {
    stop("no plot lies inside the fitting area", call. = FALSE)
  }
  draw <- withr::with_seed(seed, inside[sample.int(length(inside), replace = TRUE)])
  suffix <- paste0("_r", stats::ave(seq_along(draw), draw, FUN = seq_along))
  resampleRows <- function(x) {
    rows <- split(seq_len(nrow(x)), x[["OrigPlotID1"]])[draw]
    out <- x[unlist(rows, use.names = FALSE), ]
    copyOfRow <- rep(suffix, lengths(rows))
    for (cl in intersect(c("OrigPlotID1", "MeasureID"), names(out))) {
      data.table::set(out, j = cl, value = paste0(out[[cl]], copyOfRow))
    }
    out
  }
  psp[["PSPmeasure"]] <- resampleRows(data.table::as.data.table(psp[["PSPmeasure"]]))
  psp[["PSPplot"]] <- resampleRows(data.table::as.data.table(psp[["PSPplot"]]))
  psp[["PSPgis"]] <- resampleRows(gis)
  psp
}

#' How often refits on resampled plots chose each species' growth-curve traits
#'
#' Summarises the refits of the shared growth curves on bootstrap samples of plots
#' ([landweb_resample_psp()]): for each species, each distinct set of the five traits a growth curve
#' sets, with how many refits chose it and whether the fit to all plots chose it too.
#'
#' @param resampled `data.table` of [landweb_growth_traits()] from every refit, with the refit's number
#'   in `resample`.
#' @param fitted [landweb_growth_traits()] of the fit to all plots.
#'
#' @return `data.table` with one row per species and set of traits: `species`, the five traits,
#'   `growthCurveSource`, `n` (refits that chose the set), `share` (`n` over all refits) and `fullFit`
#'   (`TRUE` for the set the fit to all plots chose). Sorted by species, then by decreasing `n`.
#'
#' @export
landweb_growth_trait_frequency <- function(resampled, fitted) {
  rs <- data.table::as.data.table(resampled)
  keys <- c("species", growthTraitCols, "growthCurveSource")
  miss <- setdiff(c("resample", keys), names(rs))
  if (length(miss)) {
    stop("the refits' traits have no ", paste(miss, collapse = ", "), call. = FALSE)
  }
  out <- rs[, list(n = .N), by = keys]
  data.table::set(out, j = "share", value = out[["n"]] / data.table::uniqueN(rs[["resample"]]))
  full <- unique(data.table::as.data.table(fitted)[, keys, with = FALSE])
  data.table::set(full, j = "fullFit", value = TRUE)
  out <- merge(out, full, by = keys, all.x = TRUE, sort = FALSE)
  data.table::set(out, i = which(is.na(out[["fullFit"]])), j = "fullFit", value = FALSE)
  data.table::setorderv(out, c("species", "n"), order = c(1L, -1L))
  out[]
}
