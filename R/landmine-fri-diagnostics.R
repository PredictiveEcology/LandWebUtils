utils::globalVariables(c(
  "areaHa", "expIgnitions", "FRI", "largestPatchPct", "LTHFC", "nFlam", "pctFlam",
  "pctNoCohorts", "ratio", "status", "studyArea", "zone", "driver", "lyr", "pct",
  "pctTreedLCC", "pctCohortsInit", "cohortGap", "pctInStudyArea", "ratioUnmasked",
  "FRIunmasked", "pctInStudyArea"
))

## Diagnostics for FRI zones that miss their target fire return interval. The point of these is
## that a user should never have to read a simulation log to find out a zone is failing, nor run a
## bespoke investigation to find out why: the metrics below are exactly the ones that separated the
## WesternAlbertaUpland failure modes from one another.

#' Per-zone fire-return-interval diagnostics
#'
#' Extends [landmine_fri_summary()] with the structural measurements needed to tell
#' *why* a fire-return-interval zone misses its target, and classifies each zone
#' against a tolerance.
#'
#' @param lthfc `SpatRaster` of target fire return intervals (the "LTHFC").
#' @param flammableMap `SpatRaster`; `1` is flammable, `0` and `NA` are not.
#' @param meanAnnualCumulBurnMap `SpatRaster` of mean annual burn counts per pixel.
#' @param studyAreaName character, recorded in the returned table.
#' @param pixelAreaHa area of one pixel, in hectares.
#' @param vegTypeMap optional `SpatRaster` of leading vegetation type. When supplied, a
#'   `pctNoCohorts` column reports the share of each zone's flammable pixels that carry
#'   no cohorts at all -- non-forest land, not young forest. See Details.
#' @param vegTypeMapInit optional `SpatRaster` of leading vegetation type at the START
#'   of the run. Supplying it separates land that never had cohorts from land that lost
#'   them, via `pctCohortsInit` and `pctCohortsLost`.
#' @param studyArea optional `SpatVector`/`sf` polygon of the area the fire model was
#'   actually allowed to burn. **Supplying this is strongly recommended**; see Details.
#' @param lccMap optional `SpatRaster` of land cover class, on the same grid as `lthfc`.
#' @param treeClasses integer vector of `lccMap` values that count as treed. Required
#'   with `lccMap`.
#' @param tolerance length-2 numeric; achieved/target ratios outside this range are
#'   flagged. Symmetric in log space by default.
#' @param severe length-2 numeric; ratios outside this range are flagged `"severe"`.
#' @param minIgnitions integer; a zone expected to see fewer ignitions than this over
#'   the whole run is reported but not judged, because its achieved FRI is dominated
#'   by sampling noise. Set `NA` to judge every zone.
#' @param nYears length of the run in years, used with `minIgnitions`. `NA` skips the
#'   noise check.
#'
#' @details
#' `ratio` is **achieved / target**, so a value above 1 means the zone burns *less*
#' often than intended (under-burning) and below 1 means more often. That orientation
#' is what the diverging fill in [landmine_plot_fri_zones()] is centred on.
#'
#' # Why `studyArea` matters
#'
#' LandMine masks both its rate-of-spread map and its spread probability to
#' `sim$studyArea`, so pixels outside that polygon cannot be ignited or spread into. The
#' fire-return-interval raster is not similarly masked, so it commonly extends past the
#' polygon -- and [landmine_fri_summary()] divides burns by **every** flammable pixel
#' carrying an interval, unburnable ones included. Each such pixel inflates a zone's
#' achieved interval without contributing any opportunity to burn.
#'
#' This is not a small correction. On WesternAlbertaUpland, 988,437 of 3,361,240
#' flammable zone pixels (29.4%) lay outside the polygon and burned at a mean rate of
#' 0.51 against 16.59 inside. Uncorrected, two zones appeared to under-burn by 3.5x and
#' 12.2x; with the denominator restricted to the burnable area, every zone came in
#' between 0.90 and 1.07 -- there was no attainment problem at all. Supply `studyArea`
#' and both figures are reported, `ratio` from the burnable area and `ratioUnmasked`
#' the way the module's own summary computes it, so the discrepancy is visible rather
#' than mistaken for a fire-model failure.
#'
#' The structural columns exist because the failure modes are not distinguishable from
#' the ratio alone:
#'
#' * `pctFlam` -- the fraction of the zone that can carry fire at all. `spread2()`
#'   treats a non-flammable neighbour as a hard mask, so a porous zone loses perimeter
#'   growth every iteration and its fires stall short of their target size. On
#'   WesternAlbertaUpland the two failing zones were the only ones below 93% flammable
#'   (54.6% and 57.2%), which no previously-reported metric would have shown -- the
#'   landscape-wide figure was 89%.
#' * `nPatches` / `largestPatchPct` -- these separate *porous* from *fragmented*. A zone
#'   can be half non-flammable and still be one connected sheet (zone 170's largest
#'   patch held 86.8% of its flammable pixels), which is a very different problem from
#'   a zone broken into islands smaller than the fire-size draw.
#' * `pctNoCohorts` -- on WesternAlbertaUpland this tracked the shortfall more closely
#'   than anything else: every zone below 5% hit its target (ratios 0.96-1.00), while
#'   the two worst zones were at 87.0% and 99.9%. Note it means *non-forest land*, not
#'   young forest: those same zones had the OLDEST stands in the study area (median time
#'   since fire 777 and 1067 years), so a young-fuel explanation is ruled out -- high
#'   stand age is a consequence of under-burning, not its cause.
#' * `cohortGap` = `pctTreedLCC - pctCohortsInit` -- treed land that never received
#'   cohorts when the run was initialised. This separates "the zone genuinely is not
#'   forest" from "the zone is forest that was never populated", which the other columns
#'   cannot tell apart, and which imply completely different fixes. On
#'   WesternAlbertaUpland the gap was 0.1-2.1 points in every zone that met its target
#'   and 10.7-33.1 points in every zone that did not; cohort *loss* over 1000 years was
#'   at most 13.4% anywhere and ~0% in the two worst zones, so the shortfall is an
#'   initialisation property, not regeneration failure.
#'
#' @return A `data.table`, one row per zone, sorted by `LTHFC`, with columns
#'   `studyArea`, `LTHFC`, `FRI`, `ratio`, `status`, `areaHa`, `nTotal`, `nFlam`,
#'   `pctFlam`, `nPatches`, `largestPatchPct`, `pctNoCohorts`, `pctCohortsInit`,
#'   `pctCohortsLost`, `pctTreedLCC`, `cohortGap`, `pctInStudyArea`, `FRIunmasked`,
#'   `ratioUnmasked` and `expIgnitions`. The columns
#'   depending on an optional input are `NA` when it is not supplied.
#'
#' @seealso [landmine_fri_verdict()], [landmine_plot_fri_zones()]
#'
#' @export
#' @importFrom data.table data.table setorderv
#' @importFrom stats na.omit
#' @importFrom terra compareGeom deepcopy freq ncell patches rasterize values vect
landmine_fri_metrics <- function(lthfc, flammableMap, meanAnnualCumulBurnMap, studyAreaName,
                                 pixelAreaHa, studyArea = NULL,
                                 vegTypeMap = NULL, vegTypeMapInit = NULL,
                                 lccMap = NULL, treeClasses = NULL,
                                 tolerance = c(0.8, 1.25), severe = c(0.5, 2),
                                 minIgnitions = 20L, nYears = NA_real_) {
  stopifnot(length(tolerance) == 2L, length(severe) == 2L,
            tolerance[[1]] < tolerance[[2]], severe[[1]] <= tolerance[[1]],
            severe[[2]] >= tolerance[[2]])
  if (!is.null(lccMap) && is.null(treeClasses)) {
    stop("`treeClasses` must be supplied alongside `lccMap`.")
  }

  ## Every raster must be on the SAME GRID. Without this check a mismatched layer is read as a
  ## shorter vector and silently recycled, which produces a complete table of plausible, wrong
  ## numbers -- the dataPrep land-cover raster is on a larger extent than the simulation grid, so
  ## this is a live hazard rather than a hypothetical one.
  .fri_check_geom(lthfc, list(
    flammableMap = flammableMap, meanAnnualCumulBurnMap = meanAnnualCumulBurnMap,
    vegTypeMap = vegTypeMap, vegTypeMapInit = vegTypeMapInit, lccMap = lccMap
  ))

  ## `inSA` is TRUE where the fire model was permitted to burn. Everything downstream -- the
  ## achieved interval, the flammable fraction, the patch statistics -- is computed over that
  ## area only, because a pixel the model cannot reach is not evidence about the model.
  inSA <- if (is.null(studyArea)) {
    NULL
  } else {
    sv <- if (inherits(studyArea, "SpatVector")) studyArea else terra::vect(studyArea)
    !is.na(terra::values(terra::rasterize(sv, lthfc, field = 1L), mat = FALSE))
  }

  achieved <- landmine_fri_summary(
    lthfc = if (is.null(inSA)) lthfc else .fri_mask(lthfc, inSA),
    flammableMap = flammableMap,
    meanAnnualCumulBurnMap = meanAnnualCumulBurnMap, studyAreaName = studyAreaName
  )
  unmasked <- if (is.null(inSA)) {
    NULL
  } else {
    landmine_fri_summary(
      lthfc = lthfc, flammableMap = flammableMap,
      meanAnnualCumulBurnMap = meanAnnualCumulBurnMap, studyAreaName = studyAreaName
    )[, list(LTHFC, FRIunmasked = FRI)]
  }

  lthfcVals <- terra::values(lthfc, mat = FALSE)
  flamVals <- terra::values(flammableMap, mat = FALSE)
  isFlam <- !is.na(flamVals) & flamVals == 1L
  vegVals <- if (is.null(vegTypeMap)) NULL else terra::values(vegTypeMap, mat = FALSE)
  vegInitVals <- if (is.null(vegTypeMapInit)) NULL else terra::values(vegTypeMapInit, mat = FALSE)
  isTreed <- if (is.null(lccMap)) {
    NULL
  } else {
    lccVals <- terra::values(lccMap, mat = FALSE)
    !is.na(lccVals) & lccVals %in% treeClasses
  }

  ## Enumerate zones from the UNMASKED raster. `landmine_fri_summary()` reports only zones with
  ## flammable pixels, so an entirely non-flammable zone silently disappears from it -- which is
  ## exactly the configuration a diagnostic exists to surface. Such a zone is carried here with an
  ## NA achieved FRI and reported as "not evaluated".
  zones <- sort(unique(stats::na.omit(lthfcVals)))

  struct <- lapply(zones, function(z) {
    inZoneAll <- !is.na(lthfcVals) & lthfcVals == z
    inZone <- if (is.null(inSA)) inZoneAll else inZoneAll & inSA
    zoneFlam <- inZone & isFlam
    nFlamZ <- sum(zoneFlam)

    p <- .fri_patch_stats(lthfc, zoneFlam, nFlamZ)

    data.table::data.table(
      LTHFC = z,
      areaHa = sum(inZone) * pixelAreaHa,
      nTotal = sum(inZone),
      nFlam = nFlamZ,
      pctFlam = if (sum(inZone) > 0L) 100 * nFlamZ / sum(inZone) else NA_real_,
      nPatches = p$nPatches,
      largestPatchPct = p$largestPct,
      pctNoCohorts = .fri_pct(vegVals, zoneFlam, nFlamZ, function(v) is.na(v)),
      pctCohortsInit = .fri_pct(vegInitVals, zoneFlam, nFlamZ, function(v) !is.na(v)),
      pctCohortsLost = if (is.null(vegVals) || is.null(vegInitVals) || nFlamZ == 0L) {
        NA_real_
      } else {
        100 * sum(!is.na(vegInitVals[zoneFlam]) & is.na(vegVals[zoneFlam])) / nFlamZ
      },
      pctTreedLCC = .fri_pct(isTreed, zoneFlam, nFlamZ, function(v) v),
      pctInStudyArea = if (is.null(inSA)) NA_real_ else 100 * sum(inZone) / sum(inZoneAll)
    )
  })
  struct <- data.table::rbindlist(struct)

  out <- achieved[struct, on = "LTHFC"]
  if (!is.null(unmasked)) {
    out <- unmasked[out, on = "LTHFC"]
  } else {
    out[, FRIunmasked := NA_real_]
  }
  out[is.na(studyArea), studyArea := studyAreaName]
  out[, ratio := FRI / LTHFC]

  ## expected ignitions over the whole run, on the same (area / meanFireSize) / FRI footing the
  ## ignition budget uses -- but meanFireSize is not available here, so this is the per-zone
  ## *relative* opportunity: flammable area / FRI, scaled by run length.
  out[, expIgnitions := if (is.na(nYears)) NA_real_ else (nFlam * pixelAreaHa / LTHFC) * nYears]

  out[, ratioUnmasked := FRIunmasked / LTHFC]
  out[, cohortGap := pctTreedLCC - pctCohortsInit]
  out[, status := .fri_status(ratio, tolerance, severe)]
  if (!is.na(minIgnitions) && !is.na(nYears)) {
    out[expIgnitions < minIgnitions, status := "too few ignitions to judge"]
  }

  data.table::setorderv(out, "LTHFC")
  out[]
}

#' Set cells outside a logical mask to NA
#'
#' @param r `SpatRaster`.
#' @param keep logical vector, one element per cell.
#'
#' @return `SpatRaster` with `!keep` cells set to `NA`.
#'
#' @noRd
.fri_mask <- function(r, keep) {
  out <- terra::deepcopy(r)
  out[!keep] <- NA
  out
}

#' Percentage of a zone's flammable pixels satisfying a predicate
#'
#' @param vals vector of per-cell values, or `NULL`.
#' @param zoneFlam logical vector selecting the zone's flammable cells.
#' @param nFlamZ number of flammable cells in the zone.
#' @param f predicate applied to the selected values.
#'
#' @return numeric percentage, or `NA_real_` when `vals` is `NULL` or the zone is empty.
#'
#' @noRd
.fri_pct <- function(vals, zoneFlam, nFlamZ, f) {
  if (is.null(vals) || nFlamZ == 0L) {
    return(NA_real_)
  }
  100 * sum(f(vals[zoneFlam])) / nFlamZ
}

#' Require every supplied raster to share the reference grid
#'
#' @param ref reference `SpatRaster`.
#' @param others named list of `SpatRaster`s or `NULL`s.
#'
#' @return `TRUE`, invisibly; stops otherwise.
#'
#' @noRd
.fri_check_geom <- function(ref, others) {
  for (nm in names(others)) {
    r <- others[[nm]]
    if (is.null(r)) next
    if (!isTRUE(terra::compareGeom(ref, r, stopOnError = FALSE))) {
      stop(
        "`", nm, "` is not on the same grid as `lthfc` (", terra::ncell(r), " vs ",
        terra::ncell(ref), " cells). Align it first -- reading rasters of different ",
        "lengths into one table recycles the shorter one and yields a complete set of ",
        "plausible, wrong numbers."
      )
    }
  }
  invisible(TRUE)
}

#' Patch statistics for one zone's flammable pixels
#'
#' @param template `SpatRaster` giving the grid.
#' @param zoneFlam logical vector, one element per cell.
#' @param nFlamZ number of flammable cells in the zone.
#'
#' @return list with `nPatches` and `largestPct`.
#'
#' @noRd
.fri_patch_stats <- function(template, zoneFlam, nFlamZ) {
  if (nFlamZ == 0L) {
    return(list(nPatches = 0L, largestPct = NA_real_))
  }
  r <- terra::rast(template)
  terra::values(r) <- ifelse(zoneFlam, 1L, NA_integer_)
  f <- terra::freq(terra::patches(r, directions = 8, zeroAsNA = TRUE))
  if (!NROW(f)) {
    return(list(nPatches = 0L, largestPct = NA_real_))
  }
  list(nPatches = NROW(f), largestPct = 100 * max(f$count) / nFlamZ)
}

#' Classify achieved/target ratios against a tolerance
#'
#' @param ratio numeric vector of achieved/target ratios.
#' @param tolerance,severe length-2 numeric ranges.
#'
#' @return character vector of statuses.
#'
#' @noRd
.fri_status <- function(ratio, tolerance, severe) {
  out <- rep("ok", length(ratio))
  out[is.na(ratio)] <- "not evaluated"
  out[!is.na(ratio) & ratio > tolerance[[2]]] <- "under"
  out[!is.na(ratio) & ratio < tolerance[[1]]] <- "over"
  out[!is.na(ratio) & ratio > severe[[2]]] <- "under (severe)"
  out[!is.na(ratio) & ratio < severe[[1]]] <- "over (severe)"
  out
}

#' One-line verdict on fire-return-interval attainment
#'
#' Summarises [landmine_fri_metrics()] as a single sentence, for the end of a run's
#' log and for a figure caption.
#'
#' @param metrics the `data.table` from [landmine_fri_metrics()].
#' @param maxListed maximum number of failing zones to name.
#'
#' @details
#' Failing zones are named worst-first, each with its ratio and its flammable
#' fraction, because those two numbers together are usually enough to tell a
#' porosity-limited zone from one that is genuinely mis-parameterised. When
#' `pctNoCohorts` is available it is named too, since it has been the better predictor
#' of the two.
#'
#' @return A length-1 character vector.
#'
#' @export
landmine_fri_verdict <- function(metrics, maxListed = 5L) {
  judged <- metrics[status != "too few ignitions to judge" & status != "not evaluated"]
  nOK <- sum(judged$status == "ok")
  skipped <- NROW(metrics) - NROW(judged)

  head <- sprintf("FRI check: %d/%d zones within tolerance", nOK, NROW(judged))
  if (skipped > 0L) {
    head <- paste0(head, sprintf(" (%d not judged)", skipped))
  }

  ## a zone mostly outside the burnable area is the single most common cause of an apparent
  ## shortfall, and the one most likely to be mistaken for a fire-model failure
  clipped <- metrics[!is.na(pctInStudyArea) & pctInStudyArea < 90]
  note <- if (NROW(clipped)) {
    sprintf(
      " NOTE: %d zone(s) lie partly outside the burnable area (as little as %.0f%% inside); ratios use the burnable area only.",
      NROW(clipped), min(clipped$pctInStudyArea)
    )
  } else {
    ""
  }

  bad <- judged[status != "ok"]
  if (!NROW(bad)) {
    return(paste0(head, ".", note))
  }

  data.table::setorderv(bad, "ratio", order = -1L)
  bad <- utils::head(bad, maxListed)
  detail <- paste(sprintf(
    "%g (%.1fx %s, %.0f%% flammable%s)",
    bad$LTHFC, ifelse(bad$ratio > 1, bad$ratio, 1 / bad$ratio),
    ifelse(bad$ratio > 1, "under", "over"), bad$pctFlam,
    ifelse(is.na(bad$pctNoCohorts), "",
           sprintf(", %.0f%% no cohorts", bad$pctNoCohorts))
  ), collapse = ", ")

  paste0(head, ". OFF TARGET: ", detail, ".", note)
}

#' Map fire-return-interval zones and their attainment
#'
#' Two panels sharing one layout: the study area's zones filled by their **target**
#' fire return interval, and the same geometry filled by **achieved/target ratio** on
#' a diverging scale centred on 1. Zones outside tolerance are outlined so they read
#' in greyscale.
#'
#' @param lthfc `SpatRaster` of target fire return intervals.
#' @param metrics the `data.table` from [landmine_fri_metrics()].
#' @param studyAreaName character, used in the title.
#' @param caption optional character; defaults to [landmine_fri_verdict()].
#' @param maxCells passed to [tidyterra::geom_spatraster()] as the display resampling
#'   budget; the default keeps a full-resolution study area readable without rendering
#'   every cell.
#'
#' @details
#' Panel B is the one that answers "which zones, and how badly". The fill is
#' `log2(ratio)` so that burning twice as often and half as often are equally far from
#' centre; the legend is labelled in plain ratios.
#'
#' @return A `patchwork`-free `ggplot` built with [ggplot2::facet_wrap()], so no extra
#'   dependency is needed to place the two panels side by side.
#'
#' @seealso [landmine_fri_metrics()]
#'
#' @export
#' @importFrom data.table data.table melt rbindlist
#' @importFrom ggplot2 aes coord_sf element_blank element_text facet_wrap ggplot
#' @importFrom ggplot2 labs scale_fill_gradientn theme theme_minimal
#' @importFrom scales rescale
#' @importFrom terra classify
landmine_plot_fri_zones <- function(lthfc, metrics, studyAreaName, caption = NULL,
                                    maxCells = 5e5) {
  if (!requireNamespace("tidyterra", quietly = TRUE)) {
    stop("`tidyterra` is required to map FRI zones; add it to the calling context.")
  }
  if (is.null(caption)) {
    caption <- landmine_fri_verdict(metrics)
  }

  ratioRast <- terra::classify(lthfc, cbind(metrics$LTHFC, log2(metrics$ratio)),
                               others = NA_real_)
  targetRast <- terra::classify(lthfc, cbind(metrics$LTHFC, metrics$LTHFC),
                                others = NA_real_)

  panels <- c(targetRast, ratioRast)
  names(panels) <- c("A. Target fire return interval (years)",
                     "B. Achieved / target (red = burns too rarely)")

  ## one shared fill would put years and log-ratios on the same scale, so each panel is
  ## drawn with its own and they are stacked by `facet_wrap` on the layer name.
  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = panels, maxcell = maxCells) +
    ggplot2::facet_wrap(~lyr, ncol = 2L) +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166ac", "#f7f7f7", "#b2182b"),
      values = scales::rescale(c(-1, 0, 1)),
      na.value = "transparent", name = NULL
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(
      title = paste("Fire return interval attainment:", studyAreaName),
      caption = caption
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(hjust = 0, face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9),
      axis.title = ggplot2::element_blank()
    )
}

#' Plot what separates the failing zones from the rest
#'
#' Achieved/target ratio against each structural measurement that is available, one
#' point per zone. This is the figure that makes a porosity-limited zone, a genuinely
#' non-forest one, and one that is forest but was never populated with cohorts
#' distinguishable at a glance -- three failure modes the ratio alone cannot separate.
#'
#' @param metrics the `data.table` from [landmine_fri_metrics()].
#' @param studyAreaName character, used in the title.
#'
#' @return A `ggplot`.
#'
#' @export
#' @importFrom ggplot2 aes element_text facet_wrap geom_hline geom_point geom_text
#' @importFrom ggplot2 ggplot labs scale_y_continuous theme theme_bw
landmine_plot_fri_drivers <- function(metrics, studyAreaName) {
  cols <- c(pctFlam = "% of zone flammable",
            pctNoCohorts = "% of zone without cohorts",
            cohortGap = "treed by LCC but never initialised (pts)")
  cols <- cols[intersect(names(cols), names(metrics))]
  cols <- cols[vapply(names(cols), function(x) any(!is.na(metrics[[x]])), logical(1))]
  if (!length(cols)) {
    stop("`metrics` has none of `pctFlam`, `pctNoCohorts` or `cohortGap` to plot against.")
  }

  long <- data.table::melt(
    metrics[, c("LTHFC", "ratio", "status", names(cols)), with = FALSE],
    id.vars = c("LTHFC", "ratio", "status"),
    variable.name = "driver", value.name = "pct"
  )
  long[, driver := factor(cols[as.character(driver)], levels = unname(cols))]

  ggplot2::ggplot(long, ggplot2::aes(x = pct, y = ratio)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
    ggplot2::geom_point(ggplot2::aes(colour = status), size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(label = LTHFC), hjust = -0.4, size = 3) +
    ggplot2::facet_wrap(~driver, scales = "free_x") +
    ggplot2::scale_y_continuous(transform = "log2") +
    ggplot2::labs(
      x = NULL, y = "achieved / target fire return interval",
      colour = NULL,
      title = paste("What separates the off-target FRI zones:", studyAreaName),
      caption = paste(
        "Points are FRI zones, labelled by target interval. Above the dotted line the zone",
        "burns too rarely.\nThe y axis is log2, so equal distances are equal factors."
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8))
}
