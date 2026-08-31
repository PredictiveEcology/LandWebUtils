utils::globalVariables(c("maxSize", "N", "nFires", "polygonNumeric"))

## Pure helpers extracted from the LandMine module's event code, so they can be unit-tested.
## The module has no working test suite (its `tests/` are the untouched SpaDES stub), so every
## one of these was previously exercised only by running a 1000-year simulation -- where the
## failure modes below all produce plausible-looking numbers rather than errors.

#' Summarize achieved fire return interval by FRI zone
#'
#' For each distinct value of `lthfc` (the target fire return interval, held as a
#' raster of per-pixel FRI values), computes the achieved FRI as the reciprocal of
#' the mean annual burn rate over that zone's **flammable** pixels.
#'
#' @param lthfc `SpatRaster` of target fire return intervals (the "LTHFC").
#' @param flammableMap `SpatRaster`; pixels that are `NA` or `0` are non-flammable
#'   and are excluded from both the numerator and the denominator.
#' @param meanAnnualCumulBurnMap `SpatRaster` of mean annual burn counts per pixel.
#' @param studyAreaName character, recorded in the returned table.
#'
#' @details
#' This was duplicated verbatim in the LandMine module's single- and multi-mode
#' summary functions, so the two could silently diverge.
#'
#' Two edge cases are **preserved as-is** rather than "fixed", because changing them
#' would change published outputs:
#'
#' * a zone in which nothing burned yields `1/0` = `Inf`, not `NA` and not a large
#'   finite number;
#' * the burn total is a plain `sum()` with no `na.rm`, so an `NA` burn pixel inside
#'   a zone makes that whole zone's achieved FRI `NA`. Masking all three rasters by
#'   `flammableMap` is what normally prevents this -- an implicit contract that the
#'   `NA` masks of `lthfc` and `meanAnnualCumulBurnMap` agree. Both behaviours are
#'   covered by tests so a future change to either is deliberate.
#'
#' @return A `data.table` with columns `studyArea`, `LTHFC`, `FRI`, sorted by `LTHFC`.
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom terra values
landmine_fri_summary <- function(lthfc, flammableMap, meanAnnualCumulBurnMap, studyAreaName) {
  nonFlammable <- which(is.na(flammableMap[]) | flammableMap[] == 0)
  if (length(nonFlammable) > 0) {
    lthfc[nonFlammable] <- NA
    meanAnnualCumulBurnMap[nonFlammable] <- NA
  }

  lthfcVals <- terra::values(lthfc, mat = FALSE)
  expFRIs <- sort(na.omit(unique(lthfcVals)))

  simFRIs <- vapply(expFRIs, function(fri) {
    pixIds <- which(lthfcVals == fri)
    1 / (sum(meanAnnualCumulBurnMap[pixIds]) / length(pixIds))
  }, numeric(1))

  data.table::data.table(
    studyArea = studyAreaName,
    LTHFC = expFRIs,
    FRI = simFRIs,
    stringsAsFactors = FALSE
  )
}

#' Area burned per FRI zone in a single year
#'
#' Counts burned pixels of `currentBurn` within each FRI zone and converts to hectares.
#'
#' @param currentBurn `SpatRaster` of this year's burn; burned pixels are `1`, unburned `0`.
#' @param fireReturnInterval `SpatRaster` of per-pixel FRI values (the zone labels).
#' @param time numeric, the simulation time recorded in the returned table.
#' @param pixelAreaHa numeric, the area of one pixel in hectares.
#'
#' @details
#' **Fixes a bug.** This previously counted burned pixels with
#' `unname(table(currBurn[ids])[2])` -- positional indexing into a `table()`, which assumes
#' the levels are always `{0, 1}` in that order. When a zone burns **completely** there is
#' only one level, so `[2]` is `NA`, and the following `npix[is.na(npix)] <- 0` (whose comment
#' says it handles the *no-burn* case) silently converted it to its exact opposite: a
#' fully-burned zone was recorded as **zero** hectares burned. Counting with
#' `sum(vals == 1, na.rm = TRUE)` is correct in all three cases.
#'
#' @return A `data.frame` with columns `time`, `nPixelsBurned`, `haBurned`, `FRI` (a factor).
#'
#' @export
#' @importFrom terra values
landmine_area_burned_by_zone <- function(currentBurn, fireReturnInterval, time, pixelAreaHa) {
  friVals <- terra::values(fireReturnInterval, mat = FALSE)
  burnVals <- terra::values(currentBurn, mat = FALSE)
  fris <- sort(na.omit(unique(friVals)))

  npix <- vapply(fris, function(x) {
    sum(burnVals[which(friVals == x)] == 1, na.rm = TRUE)
  }, numeric(1))

  data.frame(
    time = as.numeric(time),
    nPixelsBurned = npix,
    haBurned = npix * pixelAreaHa,
    FRI = as.factor(fris)
  )
}

#' Draw the number of fires per FRI zone for one period
#'
#' @param numFiresPerYear numeric, expected fires per year for each zone.
#' @param fireTimestep numeric, years per burn event.
#' @param size numeric, the negative-binomial dispersion. The default of `1.3765` is
#'   LandMine's calibrated value (lowered from `1.8765` in 2018 because the draw was
#'   otherwise too constant); it is deliberately overdispersed, so a zone with a small
#'   `numFiresPerYear` sees long runs of zero-fire years punctuated by bursts.
#'
#' @return An integer vector, one count per element of `numFiresPerYear`, with names preserved.
#'
#' @export
#' @importFrom stats rnbinom
landmine_draw_num_fires <- function(numFiresPerYear, fireTimestep = 1, size = 1.3765) {
  out <- stats::rnbinom(
    length(numFiresPerYear),
    mu = numFiresPerYear * fireTimestep,
    size = size
  )
  ## `rnbinom()` drops names; the caller indexes fire counts by zone, so keep them.
  names(out) <- names(numFiresPerYear)
  out
}

#' Convert fire sizes in hectares to whole pixels, preserving the expectation
#'
#' @param sizesHa numeric, fire sizes in hectares.
#' @param pixelAreaHa numeric, the area of one pixel in hectares.
#'
#' @details
#' Truncating to whole pixels would systematically lose every fire smaller than one
#' pixel -- at a 6.25 ha pixel, a 3.25 ha fire is "not detectable". Rounding the
#' fractional part **probabilistically** instead keeps `E[pixels] == sizesHa/pixelAreaHa`,
#' so the fire-size distribution's mean survives the discretisation even though any
#' individual small fire may round to zero. Fires that do round to zero are dropped by
#' the caller.
#'
#' Exact multiples of `pixelAreaHa` never round up, because the fractional part is `0`
#' and `0 > runif(1)` is `FALSE`.
#'
#' @return A numeric vector of whole pixel counts, the same length as `sizesHa`.
#'
#' @export
#' @importFrom stats runif
landmine_sizes_to_pixels <- function(sizesHa, pixelAreaHa) {
  inPixels <- sizesHa / pixelAreaHa
  truncVals <- trunc(inPixels)
  decimalVals <- unname(inPixels - truncVals) > stats::runif(length(inPixels))
  truncVals + decimalVals
}

#' Summarize FRI zones still short of target at the reburn ceiling
#'
#' @param polysNeedMoreFires a `data.table` with columns `polygonNumeric` (the FRI zone),
#'   `N` (fires still needed) and `maxSize` (pixels still to burn).
#' @param year numeric, the simulation year.
#'
#' @details
#' Emitted only when the reburn loop exits on the `maxReburns` ceiling rather than because
#' every fire reached its target -- the two are otherwise indistinguishable in a log, and
#' the second means the year's area-burned target was abandoned unmet.
#'
#' `pixelsShort` means different things in the two reburn phases: in phase 1 `maxSize` is a
#' fire's *full* target, in phase 2 it is the *remaining shortfall*. The caller is
#' responsible for passing the table from the phase it means.
#'
#' @return A `data.table` with columns `year`, `FRI`, `nFires`, `pixelsShort`, containing only
#'   zones with at least one outstanding fire; zero rows when nothing is short.
#'
#' @export
#' @importFrom data.table as.data.table
landmine_reburn_ceiling <- function(polysNeedMoreFires, year) {
  if (is.null(polysNeedMoreFires) || NROW(polysNeedMoreFires) == 0L) {
    return(data.table::data.table(
      year = numeric(0), FRI = numeric(0),
      nFires = integer(0), pixelsShort = integer(0)
    ))
  }
  dt <- data.table::as.data.table(polysNeedMoreFires)
  out <- dt[
    , list(nFires = N[1], pixelsShort = sum(maxSize, na.rm = TRUE)),
    by = "polygonNumeric"
  ]
  out <- out[!is.na(polygonNumeric) & nFires > 0]
  data.table::data.table(
    year = as.numeric(year),
    FRI = as.numeric(out$polygonNumeric),
    nFires = as.integer(out$nFires),
    pixelsShort = as.integer(out$pixelsShort)
  )
}

#' Ignition budget per FRI zone
#'
#' Masks ineligible pixels out of the fire-return-interval raster, tabulates pixels per FRI
#' zone, and converts that to an expected number of fires per year per zone.
#'
#' @param fireReturnInterval `SpatRaster` whose pixel values are the target fire return
#'   interval. Zero-valued and non-flammable pixels are set to `NA`.
#' @param flammableMap `SpatRaster`; pixels that are `NA` or `0` are non-flammable.
#' @param kBest numeric, the fitted truncated-Pareto shape for the fire-size distribution.
#' @param biggestPossibleFireSizeHa numeric, the distribution's upper bound, in hectares.
#'
#' @details
#' The budget is `(zone area / mean fire size) / FRI`, so the expected area burned per year in a
#' zone equals `area / FRI` by construction -- provided each fire ignites in the zone, reaches
#' its drawn size, and burns inside the zone.
#'
#' **This function creates an ordering contract that the reburn loop then consumes
#' positionally.** Downstream, `numFiresThisPeriod[.GRP]` indexes the per-zone fire counts by
#' the *group position* of a `polygonNumeric`-keyed table, which is valid only because:
#'
#' 1. `terra::freq()` returns rows sorted ascending by value;
#' 2. the `value = NA` row is appended **last**;
#' 3. dividing by the `NA` return interval makes that entry's **value** `NA`, so a later
#'    `na.omit()` drops it -- an `NA` *name* alone would not.
#'
#' Break any one and zones are handed each other's fire counts, with no error. All three are
#' pinned by tests here, at the point where the contract is created.
#'
#' Non-flammable pixels are masked out of the raster *before* tabulation, so the budget is a
#' flammable-area budget and matches how achieved FRI is measured in
#' [landmine_fri_summary()]. It also means such pixels can never be ignition locations, since
#' the start-cell pool is built from this same masked raster.
#'
#' @return A list with:
#'   \describe{
#'     \item{`fireReturnInterval`}{the masked raster.}
#'     \item{`fireReturnIntervalsByPolygonNumeric`}{FRI per zone, ascending, with a trailing `NA`.}
#'     \item{`numFiresPerYear`}{expected fires per year per zone, named by FRI, in the same order.}
#'   }
#'
#' @export
#' @importFrom terra freq res
landmine_ignition_budget <- function(fireReturnInterval, flammableMap, kBest,
                                     biggestPossibleFireSizeHa) {
  ## fireReturnInterval should have no zeros
  zeros <- fireReturnInterval[] == 0L
  if (any(zeros, na.rm = TRUE)) {
    fireReturnInterval[zeros] <- NA_integer_
  }

  ## 2023-09: exclude non-flammable pixels for FRI calculations
  nonFlammable <- which(is.na(flammableMap[]) | flammableMap[] == 0)
  if (length(nonFlammable) > 0) {
    fireReturnInterval[nonFlammable] <- NA
  }

  ## NOTE: the NA-FRI count is kept (as a trailing row) because the burn function needs it.
  ## The module previously also built a `value = seq_len(NROW(.))` column and `order()`ed by
  ## it -- a trivially-identity permutation, since the column IS the row index. Dropped.
  freqDT <- stats::setNames(
    rbind(
      terra::freq(fireReturnInterval, bylayer = FALSE),
      terra::freq(fireReturnInterval, value = NA, bylayer = FALSE)
    ),
    c("fri", "count")
  )

  friByPolygon <- freqDT[, "fri"]
  numPixels <- stats::setNames(freqDT[, "count"], friByPolygon)

  numHa <- numPixels * (prod(terra::res(fireReturnInterval)) / 1e4)
  meanFireSizeHa <- meanTruncPareto(
    k = kBest, lower = 1, upper = biggestPossibleFireSizeHa, alpha = 1
  )

  list(
    fireReturnInterval = fireReturnInterval,
    fireReturnIntervalsByPolygonNumeric = friByPolygon,
    numFiresPerYear = (numHa / meanFireSizeHa) / friByPolygon
  )
}
