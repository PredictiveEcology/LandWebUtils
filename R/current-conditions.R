utils::globalVariables(c("CC", "LCC", "newLCC"))

## Current-condition helpers promoted from `LandWeb_preamble`. Each one encodes a decision that
## used to live inline in `InitMaps()`, where the only thing that exercised it was a full preamble
## run -- and where three separate defects (an unflagged NoData sentinel, a bounding-box
## denominator, and a coarse-grid denominator that loosened a threshold 3-20x) went unnoticed
## until someone happened to measure the output.

#' Canada LCC 2020 class groups used by LandWeb
#'
#' The single definition of which Land Cover of Canada 2020 (NALCMS level II, 30 m) codes LandWeb
#' treats as forest, non-flammable, urban, and so on. The same codes are needed by the
#' land-cover remap ([lcc2020_remap_table()]), the flammability map, and the current-condition age
#' check ([cc_age_pct_missing()]); defining them once keeps those from drifting apart.
#'
#' @return A named list of integer vectors:
#'   \describe{
#'     \item{`forest`}{`1, 2, 5, 6` -- temperate/sub-polar needleleaf, sub-polar taiga needleleaf,
#'       temperate/sub-polar broadleaf deciduous, mixed forest (LandWeb v2 CC classes 0-2).}
#'     \item{`shrubGrassMoss`}{`8, 10, 11, 12, 13` -- shrubland, grassland and the sub-polar
#'       lichen/moss classes: no vegetation dynamics, but they burn (v2 CC 3).}
#'     \item{`wetland`}{`14` -- burns (v2 CC 3).}
#'     \item{`cropland`}{`15` -- treated as grassland for fire (v2 CC 5).}
#'     \item{`nonFlammable`}{`16, 18, 19` -- barren, water, snow/ice (v2 CC 4).}
#'     \item{`urban`}{`17` -- no v2 equivalent; reclassified to the nearest type for simulation.}
#'   }
#'   The groups are disjoint.
#'
#' @export
#' @examples
#' lcc2020_classes()$forest
lcc2020_classes <- function() {
  list(
    forest = c(1L, 2L, 5L, 6L),
    shrubGrassMoss = c(8L, 10L, 11L, 12L, 13L),
    wetland = 14L,
    cropland = 15L,
    nonFlammable = c(16L, 18L, 19L),
    urban = 17L
  )
}

#' Remap table for overlaying LCC 2020 onto the simulation land cover
#'
#' Builds the `remapTable` that `LandR::overlayLCCs()` uses to decide, for every combination of a
#' simulation land-cover code (`LCC`, SCANFI-derived) and a current-condition LCC 2020 code (`CC`),
#' which code the pixel ends up with. Promoted verbatim from `LandWeb_preamble`, which mirrors
#' LandWeb v2's rule ordering with LCC 2020 codes substituted for v2's CC 0-5.
#'
#' @param lccClasses integer vector of the simulation land-cover codes present (`NA` is added).
#' @param ccClasses integer vector of the LCC 2020 codes present; `NA`s are dropped, then `NA` is
#'   added back as its own row.
#' @param treeClassesToReplace integer vector of `LCC` codes that must be reclassified (set to 99).
#' @param classes class groups, as returned by [lcc2020_classes()].
#'
#' @details
#' **Rule order matters: later rules overwrite earlier ones.** In order:
#'
#' 1. `LCC` 0, 20, 30 -> `NA`;
#' 2. `CC` missing or cropland -> keep `LCC`;
#' 3. `CC` non-flammable (barren/water/snow-ice) -> `NA`;
#' 4. `CC` forest, shrub/grass/moss or wetland -> keep `LCC`;
#' 5. `LCC` missing but `CC` forest -> 99;
#' 6. `CC` urban -> 99;
#' 7. `LCC` in `treeClassesToReplace` -> 99.
#'
#' `CC` does **not** override the simulation's forest determination: every `CC` class except
#' barren/water/snow-ice defers to `LCC`, so LCC 2020 calling a treed wetland "wetland" does not
#' strip forest from the simulation.
#'
#' Note that rule 1 is overwritten by rules 2-6 for **every** code in the LCC 2020 legend, so it only
#' takes effect for a `CC` value outside that legend. It is kept for fidelity with v2 rather than
#' removed, and the behaviour is pinned by a test.
#'
#' @return A `data.table` with integer columns `LCC`, `CC` and `newLCC`, one row per combination.
#'
#' @export
#' @examples
#' lcc2020_remap_table(
#'   lccClasses = c(20L, 30L, 210L, 220L),
#'   ccClasses = c(1L, 14L, 17L, 18L),
#'   treeClassesToReplace = 240L
#' )
lcc2020_remap_table <- function(lccClasses, ccClasses, treeClassesToReplace,
                                classes = lcc2020_classes()) {
  remapDT <- expand.grid(
    LCC = c(NA_integer_, sort(as.integer(lccClasses))),
    CC = c(NA_integer_, sort(unique(as.integer(stats::na.omit(ccClasses)))))
  ) |>
    data.table::as.data.table()
  remapDT[LCC %in% c(0L, 20L, 30L), newLCC := NA_integer_]
  remapDT[is.na(CC) | CC %in% classes$cropland, newLCC := LCC]
  remapDT[CC %in% classes$nonFlammable, newLCC := NA_integer_]
  remapDT[CC %in% c(classes$forest, classes$shrubGrassMoss, classes$wetland), newLCC := LCC]
  remapDT[is.na(LCC) & CC %in% classes$forest, newLCC := 99L]
  remapDT[CC %in% classes$urban, newLCC := 99L]
  remapDT[LCC %in% treeClassesToReplace, newLCC := 99L]
  remapDT[]
}

#' Composite current-condition stand age, filling gaps from a second source
#'
#' Crops the primary age raster to the study area, converts its NoData sentinel to `NA`, and
#' optionally fills the cells it leaves empty from a second age raster aged forward to the same
#' epoch. Promoted from `LandWeb_preamble`, where `base` is fRI Research's `age_in2025` and `fill`
#' is NTEMS forest age.
#'
#' @param base `SpatRaster` of stand age at `epoch`; authoritative wherever it has a value.
#' @param studyArea `SpatVector` or `sf` polygon(s). Both rasters are cropped to it **before** any
#'   other work (see Details).
#' @param fill optional `SpatRaster` of stand age at `fillYear`, in any CRS.
#' @param fillYear the year `fill` represents. Required with `fill`.
#' @param epoch the year `base` represents; `fill` is aged forward by `epoch - fillYear`.
#' @param baseNAflag,fillNAflag the NoData sentinel of each source (`age_in2025`: 65535; NTEMS: 255,
#'   meaning non-treed). `NA` if the source has none.
#'
#' @details
#' # Fill, never minimum
#' `base` stays authoritative: `fill` only supplies cells `base` left empty. A cell-wise minimum
#' would instead let a *modelled* age drag *inventory* ages younger, and because a minimum is
#' one-directional that is a systematic young bias into exactly the areas where the data is best.
#'
#' # Sentinels
#' A sentinel left unconverted is read as a real age -- 65,535 years, or a 255-year stand on every
#' non-treed pixel -- and silently poisons every downstream aggregate. `terra` does not read
#' `age_in2025`'s NoData tag, so the sentinel must be supplied. It is converted with
#' `terra::classify()` after cropping rather than with `NAflag<-`: for an in-memory raster,
#' `NAflag<-` converts the values but `NAflag()` still reports `NaN`, so an assertion on the flag
#' fails on correct data, and whether the flag "took" depends on how the raster is stored.
#'
#' `fill`'s sentinel is converted **before** it is reprojected and aged, so it can never surface as
#' `sentinel + (epoch - fillYear)` (e.g. 258).
#'
#' # Cost
#' Cropping first matters: `age_in2025` alone is about 6.5e9 cells, so aligning a fill at full
#' extent would move billions of cells to produce one study area's output. `fill` is reprojected
#' with `method = "near"`, so the result contains only the source's own ages.
#'
#' @return A `SpatRaster` on `base`'s grid, cropped and masked to `studyArea`.
#'
#' @seealso [cc_age_pct_missing()] to check the result is complete enough to use.
#' @export
cc_age_composite <- function(base, studyArea, fill = NULL, fillYear = NULL, epoch = 2025L,
                             baseNAflag = 65535, fillNAflag = 255) {
  stopifnot(inherits(base, "SpatRaster"))
  sa <- .asSpatVector(studyArea)

  age <- terra::crop(base, terra::project(sa, terra::crs(base)), mask = TRUE)
  age <- .sentinelToNA(age, baseNAflag)

  if (!is.null(fill)) {
    stopifnot(inherits(fill, "SpatRaster"))
    if (is.null(fillYear) || length(fillYear) != 1L || is.na(fillYear)) {
      stop("`fillYear` is required when `fill` is supplied.", call. = FALSE)
    }
    offset <- as.integer(epoch - fillYear)
    if (offset < 0L) {
      stop("`fillYear` (", fillYear, ") is after `epoch` (", epoch, ").", call. = FALSE)
    }
    src <- terra::crop(fill, terra::project(sa, terra::crs(fill)), mask = TRUE)
    src <- .sentinelToNA(src, fillNAflag)
    age <- terra::cover(age, terra::project(src, age, method = "near") + offset)
  }
  age
}

#' Percent of forest with no current-condition stand age
#'
#' The completeness check for a current-condition age raster: of the forest inside the study area,
#' what share has no age? Promoted from `LandWeb_preamble`'s `ccAgeMaxMissing` gate.
#'
#' @param age `SpatRaster` of stand age, e.g. from [cc_age_composite()].
#' @param studyArea `SpatVector` or `sf` polygon(s).
#' @param landcover `SpatRaster` of land cover, in any CRS and at any resolution; it is cropped
#'   and reprojected onto `age`'s grid with `method = "near"`.
#' @param forestClasses the `landcover` codes that count as forest. Defaults to the LCC 2020 forest
#'   classes from [lcc2020_classes()].
#'
#' @details
#' The denominator is the part of the study area that is forest. It took three attempts to get
#' there, and each wrong version is pinned by a test:
#'
#' 1. **Not the raster.** Cropping with `mask = TRUE` sets every cell outside the polygon to `NA`,
#'    indistinguishable from a cell that genuinely lacks an age, so `mean(is.na(age))` measures the
#'    shape of the bounding box. LandWeb's study-area groups fill only 30-44% of theirs:
#'    WesternAlbertaUpland read 71.4% missing against a true 3.9%.
#' 2. **Not the whole study area.** Stand age is only defined for forest, and a treed-only source
#'    such as NTEMS leaves water and wetland empty by design, so the polygon-wide figure mostly
#'    measures lakes. ChurchillRiverUpland read 39.4% missing after filling; 76% of that was water,
#'    barren or wetland.
#' 3. **Not on a coarser grid.** Aggregating with `average` fills a coarse cell from any one of its
#'    children, so a cell only counts as missing when *all* of them are. At 240 m that loosened the
#'    figure 3-20x and let an unfilled study area pass.
#'
#' Counting uses `terra::global()`, which works in chunks, so no full-resolution vector is held in
#' memory.
#'
#' @return A single number, the percent (0-100) of forest cells with no age, with attributes
#'   `nForest` and `nAged` giving the counts behind it.
#'
#' @export
cc_age_pct_missing <- function(age, studyArea, landcover,
                               forestClasses = lcc2020_classes()$forest) {
  stopifnot(inherits(age, "SpatRaster"), inherits(landcover, "SpatRaster"))
  sa <- terra::project(.asSpatVector(studyArea), terra::crs(age))

  inStudyArea <- terra::rasterize(sa, age, field = 1L)
  lc <- terra::crop(landcover, terra::project(sa, terra::crs(landcover))) |>
    terra::project(age, method = "near")
  ## `classify()`, not `ifel(lc %in% forestClasses, ...)`: `%in%` is a plain base closure, so terra's
  ## SpatRaster method is only found when terra is ATTACHED. Called from a package namespace (or any
  ## session that uses `terra::` without `library(terra)`) it silently resolves to base `%in%`,
  ## which returns a logical vector instead of a raster.
  forest <- terra::classify(lc, cbind(as.numeric(forestClasses), 1), others = NA) |>
    terra::mask(inStudyArea)

  nForest <- terra::global(forest, "notNA")[[1L]]
  if (!isTRUE(nForest > 0)) {
    stop(
      "No forest (classes ", paste(forestClasses, collapse = ", "), ") in the study area, so ",
      "the share of forest missing an age is undefined.",
      call. = FALSE
    )
  }
  nAged <- terra::global(terra::mask(age, forest), "notNA")[[1L]]
  structure(100 * (1 - nAged / nForest), nForest = nForest, nAged = nAged)
}

## `sf`/`sfc` -> `SpatVector`; a `SpatVector` passes through.
.asSpatVector <- function(x) {
  if (inherits(x, "SpatVector")) {
    return(x)
  }
  if (inherits(x, c("sf", "sfc"))) {
    return(terra::vect(sf::st_as_sf(x)))
  }
  stop("`studyArea` must be a SpatVector or an sf object.", call. = FALSE)
}

## Replace a NoData sentinel with NA by value, independent of how the raster is stored.
.sentinelToNA <- function(x, flag) {
  if (is.null(flag) || is.na(flag)) {
    return(x)
  }
  terra::classify(x, cbind(as.numeric(flag), NA))
}
