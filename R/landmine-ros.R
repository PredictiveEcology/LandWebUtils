utils::globalVariables(c(
  "age", "LandMine", "LandWeb", "leading", "pixelValue", "ros", "species", "used", "V2"
))

## The rate-of-spread calculation lifted out of the LandMine module's `fireROS()`, so it can be
## unit-tested without running a simulation. In the module it read eight values off `sim`/`mod`/
## `P(sim)`, which made every one of its failure modes -- all of which produce a plausible-looking
## ROS map rather than an error -- reachable only through a multi-hour run.

#' Leading-species to fuel-type lookup used by LandMine
#'
#' The mapping from `LandWeb`-column species codes to the four LandMine fuel types
#' (`spruce`, `pine`, `decid`, `softwood`) that [landmine_fire_ros()] joins the
#' `ROSTable` on.
#'
#' @details
#' The LandMine module previously defined this inline in `init()` and used it in two
#' places, one of them inside `fireROS()`. Both now call this, so the fuel-type
#' assignment cannot drift between them.
#'
#' Note the keys are `LandWeb`-column codes, **not** `LandR`-column codes.
#'
#' @return A named character vector; names are species codes, values are fuel types.
#'
#' @export
#' @examples
#' landmine_known_species()[["Pice_mar"]] ## "spruce"
landmine_known_species <- function() {
  c(
    Abie_spp = "softwood",
    Lari_spp = "decid",
    Pice_gla = "spruce",
    Pice_mar = "spruce",
    Pinu_spp = "pine",
    Popu_spp = "decid",
    Pseu_men = "softwood",
    Thuj_pli = "softwood",
    Tsug_het = "softwood"
  )
}

#' Build the pixel-value to rate-of-spread lookup
#'
#' Maps each raster attribute table (RAT) entry of a vegetation type map onto a
#' LandMine fuel type, then joins the `ROSTable` to give one rate of spread per
#' (pixel value, age class).
#'
#' @param vegTypes `data.table` of the vegetation type map's RAT: column 1 the
#'   raster value, column 2 the label.
#' @param sppEquiv `data.table` with (at least) columns `LandMine` and `LandWeb`.
#' @param sppEquivCol character; the `sppEquiv` column holding the labels used in
#'   `vegTypes`.
#' @param ROSTable `data.table` with columns `age`, `leading`, `ros`.
#' @param knownSpecies named character vector, as from [landmine_known_species()].
#'
#' @return A `data.table` keyed on `used` (one of `"mature"`, `"immature"`,
#'   `"young"`, or `"no"`) with columns `leading`, `age`, `ros`, `pixelValue`, `used`.
#'
#' @noRd
#' @importFrom data.table as.data.table data.table rbindlist setkeyv setnames
#' @importFrom stats na.omit
.ros_lookup <- function(vegTypes, sppEquiv, sppEquivCol, ROSTable, knownSpecies) {
  sppNames <- equivalentName(as.character(vegTypes[[2]]), sppEquiv, sppEquivCol)

  ## NOTE: these `grep()` calls return POSITIONS in `sppNames` (i.e. RAT row numbers), which are
  ## then used as raster VALUES ("pixelValue"). That is only correct when the RAT's value column
  ## is 1..nrow, which `.ros_check_rat()` asserts.
  suppressWarnings({
    onRaster <- data.table::rbindlist(list(
      list("mixed", which(is.na(sppNames))),
      list("spruce", grep(sppNames, pattern = "Pice")),
      list("pine", grep(sppNames, pattern = "Pinu")),
      list("decid", grep(sppNames, pattern = "Popu")),
      list("softwood", grep(sppNames, pattern = "Pice|Pinu|Popu", invert = TRUE))
    ))
  })

  ## remove duplicates of softwood, which is NA
  onRaster <- stats::na.omit(unique(onRaster, by = "V2")) |>
    data.table::setnames(old = 1:2, new = c("leading", "pixelValue")) |>
    data.table::setkeyv("pixelValue")

  ## `sppNames` is assigned positionally onto a table that has been deduplicated and re-sorted by
  ## key. That only aligns because every RAT position appears exactly once (the "softwood" group is
  ## the complement of the other three, so `unique(by = "V2")` keeps one row each) and the key sorts
  ## them back into 1..n order.
  if (NROW(onRaster) != length(sppNames)) {
    stop(
      "internal error building the ROS lookup: ", NROW(onRaster), " rows for ",
      length(sppNames), " vegetation types; `species` would be recycled, silently ",
      "mis-assigning rates of spread."
    )
  }
  onRaster[, species := sppNames]

  ## `pixelValue` has so far held attribute-table ROW POSITIONS, because that is what `which()` and
  ## `grep()` return. Translate them into the raster VALUES they denote. The two coincide when the
  ## table is in value order -- which is how `LandR::vegTypeMapGenerator` builds it, `ID =
  ## seq_along(levels(leading))` -- but NOT after a `writeRaster()`/`rast()` round-trip, which can
  ## return the rows in a different order while keeping each ID paired with its own label. Reading
  ## a saved `vegTypeMap` back off disk and relying on positions therefore assigns whole species
  ## the wrong fuel type, silently.
  onRaster[, pixelValue := vegTypes[[1]][pixelValue]]
  data.table::setkeyv(onRaster, "pixelValue")

  sppEquiv <- sppEquiv[, c("LandMine", "LandWeb")][, leading := knownSpecies[LandWeb]] |>
    stats::na.omit(on = "LandMine")
  sppEquiv <- sppEquiv[onRaster, on = c("LandMine" = "leading", "LandWeb" = "species")] |>
    unique()

  sppEquivHere <- unique(stats::na.omit(sppEquiv$LandWeb))
  haveAllKnown <- sppEquivHere %in% names(knownSpecies)
  if (!all(haveAllKnown)) {
    stop(
      "LandMine only has rate of spread burn rates for\n",
      paste(names(knownSpecies), collapse = ", "),
      "\nMissing rate of spread for ", paste(sppEquivHere[!haveAllKnown], collapse = ", ")
    )
  }

  ## NOTE: `knownSpecies` maps species onto spruce/pine/decid/softwood only -- never "mixed", which
  ## `onRaster` invents for RAT entries that match no species. The join above therefore maps "mixed"
  ## onto `LandMine` and leaves `leading` NA, and the `ROSTable` join below (on `leading`) drops the
  ## mixed row entirely. Mixedwood pixels consequently never receive the `ROSTable`'s mixed rates and
  ## fall through to `ROSother`. Preserved as-is -- see the "Mixedwood" section of
  ## [landmine_fire_ros()] -- and pinned by a test.
  sppEquiv <- unique(sppEquiv, by = c("LandMine", "leading", "pixelValue"))
  sppEquiv <- sppEquiv[ROSTable, on = "leading", allow.cartesian = TRUE, nomatch = NULL]
  sppEquiv <- sppEquiv[, c("leading", "age", "ros", "pixelValue")]
  sppEquiv <- unique(sppEquiv, by = c("age", "leading", "pixelValue"))

  ## a compound age label such as "immature_young" is claimed by the OLDEST class it names, so the
  ## order of these three assignments is load-bearing.
  sppEquiv[, used := "no"]
  sppEquiv[(used == "no") & grepl("(^|_)mature", age), used := "mature"]
  sppEquiv[(used == "no") & grepl("(^|_)immature", age), used := "immature"]
  sppEquiv[(used == "no") & grepl("(^|_)young", age), used := "young"]
  data.table::setkeyv(sppEquiv, "used")

  sppEquiv
}

#' Build the three age-class masks
#'
#' @param lookup the `data.table` returned by `.ros_lookup()`.
#' @param tsf integer vector of time since fire, one element per pixel.
#' @param vegType integer vector of vegetation type raster values, one per pixel.
#' @param vegTypeIDs the RAT's value column.
#' @param ageCutoffs length-2 numeric, the young/immature and immature/mature
#'   boundaries.
#' @param youngFilter see [landmine_fire_ros()].
#'
#' @return A named list of three logical vectors: `mature`, `immature`, `young`.
#'
#' @noRd
.ros_age_cuts <- function(lookup, tsf, vegType, vegTypeIDs, ageCutoffs, youngFilter) {
  cuts <- list()

  ## a `ROSTable` that already merges mature with another class needs no upper age cut: every
  ## non-NA pixel is eligible, and the species filter below does the discriminating.
  if (!any(grepl("_mature$|^mature_|_mature_", lookup$age))) {
    cuts[[1]] <- tsf > ageCutoffs[[2]]
  } else {
    cuts[[1]] <- !is.na(tsf)
  }

  if (!any(grepl("_immature$|^immature_|_immature_", lookup$age))) {
    cuts[[2]] <- tsf > ageCutoffs[[1]] & tsf <= ageCutoffs[[2]]
  } else {
    cuts[[2]] <- tsf <= ageCutoffs[[2]]
  }

  cuts[[3]] <- tsf <= ageCutoffs[[1]]

  ## Now go through from mature through immature through young.
  ##
  ## Each filter restricts a class to the vegetation types that actually have a rate in it. That
  ## only changes the answer for the YOUNG class, because young is a strict subset of immature:
  ## an unfiltered young pixel whose type has no young rate is assigned `NA`, overwriting the
  ## immature rate already assigned to it. Mature and immature are disjoint from what precedes
  ## them, so an `NA` there is just an `NA`.
  ##
  ## The two guards below are therefore equivalent to applying the filter unconditionally: every
  ## `pixelValue` in the lookup came from the attribute table, so `all(... %in% vegTypeIDs)` is
  ## TRUE whenever the class has any rows at all, and FALSE only for the `NA` row a keyed subset
  ## returns when it has none -- in which case the filter empties the class, which is right.
  if (!all(lookup["mature"]$pixelValue %in% vegTypeIDs)) {
    cuts[[1]] <- cuts[[1]] & vegType %in% lookup["mature"]$pixelValue
  }

  if (!all(lookup["immature"]$pixelValue %in% vegTypeIDs)) {
    cuts[[2]] <- cuts[[2]] & vegType %in% lookup["immature"]$pixelValue
  }

  applyYoung <- identical(youngFilter, "always") ||
    all(lookup["young"]$pixelValue %in% vegTypeIDs)
  if (applyYoung) {
    cuts[[3]] <- cuts[[3]] & vegType %in% lookup["young"]$pixelValue
  }

  stats::setNames(cuts, c("mature", "immature", "young"))
}

#' Verify that a derived ROS map agrees with the lookup that produced it
#'
#' Re-derives, for every (pixel value, age class) pair present in the output, the
#' rate of spread actually assigned, and compares it against the lookup.
#'
#' @param ROS integer vector of assigned rates of spread.
#' @param vegType integer vector of vegetation type raster values.
#' @param tsf integer vector of time since fire.
#' @param cuts the list returned by `.ros_age_cuts()`.
#' @param lookup the `data.table` returned by `.ros_lookup()`.
#' @param ageCutoffs length-2 numeric, as for [landmine_fire_ros()].
#'
#' @return `TRUE`, invisibly; stops on disagreement.
#'
#' @noRd
.ros_self_test <- function(ROS, vegType, tsf, cuts, lookup, ageCutoffs) {
  ## NOTE: `cut()` is right-closed, so `tsf == 0` and `tsf > 999` both land in `NA` and are dropped
  ## below, whereas `.ros_age_cuts()` puts `tsf == 0` in the young class and anything above the
  ## upper cutoff in mature. Those pixels have never been covered by this check; `test-landmine-
  ## ros.R` covers them instead.
  dt <- data.table::data.table(
    ROS = ROS,
    pixelValue = vegType,
    age = cut(tsf,
      breaks = c(0, ageCutoffs[[1]], ageCutoffs[[2]], 999),
      labels = c("young", "immature", "mature")
    ),
    data.table::as.data.table(cuts)
  )
  dt <- stats::na.omit(dt, cols = c("ROS", "age"))
  dtSumm <- dt[, list(derivedROS = unique(ROS)), by = c("pixelValue", "age")]
  dtSumm <- dtSumm[lookup, on = c("pixelValue", "age" = "used"), nomatch = NULL]

  if (!identical(dtSumm$derivedROS, dtSumm$ros)) {
    mismatched <- dtSumm[dtSumm$derivedROS != dtSumm$ros]
    stop(
      "fireROS failed its test: ", NROW(mismatched), " of ", NROW(dtSumm),
      " (pixel value, age class) combinations were assigned a rate of spread that ",
      "disagrees with the lookup table."
    )
  }

  invisible(TRUE)
}

#' Assign a rate of spread to every pixel
#'
#' Combines a vegetation type map, a stand-age (time since fire) map and a table of
#' rates of spread by fuel type and age class into one rate-of-spread vector, then
#' fills the remaining flammable and non-flammable pixels with their defaults.
#'
#' @param vegTypeMap `SpatRaster` of leading vegetation type, **with** a raster
#'   attribute table.
#' @param rstTimeSinceFire `SpatRaster` of time since fire, or a numeric vector with
#'   one element per pixel of `vegTypeMap`.
#' @param flammableMap `SpatRaster` (or numeric vector) in which `1` marks flammable
#'   pixels and `0` or `NA` marks non-flammable ones.
#' @param ROSTable `data.table` with columns `age`, `leading`, `ros`. The `age`
#'   values may be `"mature"`, `"immature"`, `"young"`, or compound versions such as
#'   `"immature_young"` where two classes share one rate.
#' @param sppEquiv `data.table` of species equivalencies, with (at least) columns
#'   `LandMine`, `LandWeb`, and `sppEquivCol`.
#' @param sppEquivCol character; the `sppEquiv` column whose values appear as the
#'   labels of `vegTypeMap`'s attribute table.
#' @param ROSother integer; the rate of spread for flammable pixels that are not
#'   forested (grassland, lichen, shrub). Must lie within the range of
#'   `ROSTable$ros` and within 5% of mature spruce.
#' @param knownSpecies named character vector mapping species codes to fuel types;
#'   defaults to [landmine_known_species()].
#' @param ROStype character; `"default"` leaves non-flammable pixels as `NA`,
#'   `"burny"` gives them the young-deciduous rate so that fire can cross
#'   discontinuous fuels.
#' @param ageCutoffs length-2 numeric giving the young/immature and immature/mature
#'   stand-age boundaries, in years. These are the only place the `ROSTable`'s age
#'   **labels** acquire a numeric meaning -- the table itself carries no ages.
#' @param youngFilter character; whether the young age class is always restricted to
#'   the vegetation types that have a young rate (`"always"`), or only under the
#'   module's original condition (`"legacy"`, the default). See Details.
#' @param assertions logical; run `.ros_self_test()` and the attribute-table check.
#'
#' @details
#' # The `youngFilter` argument
#'
#' The module applied the species filter to the mature and immature classes when
#' `!all(pixelValue %in% vegTypeIDs)`, but to the young class when
#' `all(pixelValue %in% vegTypeIDs)` -- the opposite condition. The git history shows
#' this is an unfinished copy-paste: the line was introduced (`ad5ddb3`) as a copy of
#' the mature guard that still read `sppEquiv["mature"]` in its condition while
#' assigning to the young class, and the follow-up (`b6d600c`) repaired the subset key
#' but not the dropped `!`.
#'
#' Restoring the `!` is nonetheless **not** the fix. Because every `pixelValue` in the
#' lookup comes from the attribute table, the condition is only ever asking whether the
#' class has any rows at all, so the correct behaviour -- for all three classes -- is to
#' apply the filter unconditionally. That is what `"always"` does.
#'
#' `"legacy"` and `"always"` agree except in one case: when no vegetation type present
#' in the landscape has a `"young"`-classified age row. Then `"legacy"` skips the filter,
#' every young pixel is assigned `NA` from an empty lookup, and that `NA` overwrites the
#' immature rate the pixel had already been given -- so young stands silently fall
#' through to `ROSother` instead of their `immature_young` rate. With LandMine's default
#' `ROSTable`, `"young"` appears only for pine, so this bites on any landscape with no
#' pine in it, and not otherwise.
#'
#' The default is `"legacy"` so that promoting this function out of the module changed
#' no results.
#'
#' # Attribute-table row order
#'
#' The lookup is built from `which()` and `grep()` over the attribute table, which give
#' row *positions*; those are translated into raster *values* via the table's own value
#' column. Row position and raster value coincide as `LandR::vegTypeMapGenerator` builds
#' the map, but not necessarily after a `terra::writeRaster()`/`terra::rast()` round
#' trip, which may return the rows in a different order. Code that reloads a saved
#' `vegTypeMap` and indexes by position -- as this function's predecessor in the LandMine
#' module did -- assigns whole species the wrong fuel type with no error. Measured on a
#' saved WesternAlbertaUpland `vegTypeMap`, that meant aspen and pine, 36% of the
#' landscape between them, taking each other's and softwood's rates of spread.
#'
#' # Mixedwood stands never use the table's mixed rates
#'
#' `knownSpecies` maps species onto `spruce`, `pine`, `decid` and `softwood` only.
#' The fifth fuel type, `mixed`, is invented for attribute-table entries that match
#' no species, so it exists on one side of the species join and not the other: the
#' join leaves its `leading` value `NA`, and the subsequent join onto `ROSTable`
#' (which is keyed on `leading`) drops it. Mixedwood pixels therefore never receive
#' the `ROSTable`'s `mixed` rates at any age, and instead fall through to `ROSother`.
#'
#' With LandMine's default table this means mixedwood burns at 30 rather than 12
#' (young/immature) or 17 (mature) -- i.e. at the mature-spruce rate -- and both
#' `mixed` rows of the table are dead entries. This is preserved here rather than
#' fixed, because changing it changes the fire regime; a test pins the behaviour so
#' that any future change is deliberate.
#'
#' @return An integer vector of rates of spread, one element per pixel of
#'   `vegTypeMap`.
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom SpaDES.tools inRange
#' @importFrom terra levels ncell values
landmine_fire_ros <- function(vegTypeMap, rstTimeSinceFire, flammableMap, ROSTable,
                              sppEquiv, sppEquivCol, ROSother,
                              knownSpecies = landmine_known_species(),
                              ROStype = c("default", "burny"),
                              ageCutoffs = c(young = 40, immature = 120),
                              youngFilter = c("legacy", "always"),
                              assertions = getOption("LandR.assertions", TRUE)) {
  ROStype <- match.arg(ROStype)
  youngFilter <- match.arg(youngFilter)

  stopifnot(length(ageCutoffs) == 2L, ageCutoffs[[1]] < ageCutoffs[[2]])

  nCells <- terra::ncell(vegTypeMap)
  ROS <- rep(NA_integer_, nCells)

  vegType <- terra::values(vegTypeMap, mat = FALSE) ## vector, not matrix (else pixelValue col malformed)
  rat <- terra::levels(vegTypeMap)
  rat <- if (length(rat) >= 1L) rat[[1L]] else NULL
  if (is.null(rat) || !is.data.frame(rat) || NROW(rat) == 0L) {
    stop(
      "`vegTypeMap` has no raster attribute table, so its values cannot be mapped to ",
      "species. `terra::writeRaster()` stores the attribute table in an `.aux.xml` ",
      "sidecar -- copying the `.tif` alone silently drops it."
    )
  }
  vegTypes <- data.table::data.table(rat)

  tsf <- .as_cell_vector(rstTimeSinceFire, nCells, "rstTimeSinceFire")
  flammable <- .as_cell_vector(flammableMap, nCells, "flammableMap")

  lookup <- .ros_lookup(vegTypes, sppEquiv, sppEquivCol, ROSTable, knownSpecies)

  if (isTRUE(assertions)) {
    .ros_check_rat(vegTypes)
  }

  cuts <- .ros_age_cuts(lookup, tsf, vegType, vegTypes[[1]], ageCutoffs, youngFilter)

  mature <- which(cuts[["mature"]])
  immature <- which(cuts[["immature"]])
  young <- which(cuts[["young"]])

  if (length(mature)) {
    ROS[mature] <- lookup["mature"]$ros[match(vegType[mature], lookup["mature"]$pixelValue)]
  }

  if (length(immature)) {
    ROS[immature] <- lookup["immature"]$ros[match(vegType[immature], lookup["immature"]$pixelValue)]
  }

  if (length(young)) {
    ROS[young] <- lookup["young"]$ros[match(vegType[young], lookup["young"]$pixelValue)]
  }

  if (isTRUE(assertions)) {
    .ros_self_test(ROS, vegType, tsf, cuts, lookup, ageCutoffs)
  }

  ## Other vegetation (flammable, non-forest; e.g., grasslands, lichen, shrub)
  ## The original default value is the same as that of mature spruce stands (30L)
  ROSotherRef <- ROSTable[leading == "spruce" & age == "mature", ros]
  if (length(ROSotherRef) != 1L) {
    stop(
      "`ROSTable` must have exactly one mature-spruce row to validate `ROSother` ",
      "against; found ", length(ROSotherRef), "."
    )
  }

  ## Non-flammable cover types
  ## 2023-02: discontinuous fuels (e.g., shield) may require spread through non-flammable pixels;
  ##          use same value as young deciduous (6L), per Dave's text messages.
  ROSnonflam <- switch(ROStype,
    burny = ROSTable[leading == "decid" & age == "immature_young", ros],
    NA_integer_
  )
  if (length(ROSnonflam) != 1L) {
    stop(
      "`ROStype = \"burny\"` needs exactly one immature_young deciduous row in ",
      "`ROSTable` to use for non-flammable pixels; found ", length(ROSnonflam), "."
    )
  }

  ## the second of these previously compared `ROSother` against +/-5% of itself, which is
  ## vacuously true; it is meant to hold `ROSother` near the mature-spruce rate its default matches.
  if (!isTRUE(SpaDES.tools::inRange(ROSother, min(ROSTable$ros), max(ROSTable$ros)))) {
    stop(
      "`ROSother` (", ROSother, ") is outside the range of `ROSTable$ros` (",
      min(ROSTable$ros), "-", max(ROSTable$ros), ")."
    )
  }
  ## TODO: tweak this to allow greater range
  if (!isTRUE(SpaDES.tools::inRange(ROSother, 0.95 * ROSotherRef, 1.05 * ROSotherRef))) {
    stop(
      "`ROSother` (", ROSother, ") is not within 5% of the mature spruce rate of ",
      "spread (", ROSotherRef, ")."
    )
  }

  ROS[flammable == 1L & is.na(ROS)] <- as.integer(ROSother) ## flammable
  ROS[flammable == 0L | is.na(flammable)] <- as.integer(ROSnonflam) ## nonflammable

  ROS
}

#' Coerce a raster or vector of per-pixel values to a plain vector
#'
#' @param x `SpatRaster` or vector.
#' @param nCells expected length.
#' @param what name of the argument, for error messages.
#'
#' @return A vector of length `nCells`.
#'
#' @noRd
.as_cell_vector <- function(x, nCells, what) {
  out <- if (inherits(x, "SpatRaster")) terra::values(x, mat = FALSE) else x
  if (length(out) != nCells) {
    stop(
      "`", what, "` has ", length(out), " values but `vegTypeMap` has ", nCells,
      " cells; they must be on the same grid."
    )
  }
  out
}

#' Check a raster attribute table's value column
#'
#' @param vegTypes `data.table` of the RAT.
#'
#' @return `TRUE`, invisibly; stops otherwise.
#'
#' @noRd
.ros_check_rat <- function(vegTypes) {
  ids <- vegTypes[[1]]
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop(
      "`vegTypeMap`'s attribute table has ",
      if (anyNA(ids)) "missing " else "duplicate ",
      "values, so a raster value cannot be mapped to one vegetation type."
    )
  }
  invisible(TRUE)
}
