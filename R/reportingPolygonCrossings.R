## ---- tenure x sub-region reporting-polygon crossings -------------------------------------------
## v2 reported each sub-region layer *per tenure*: an "<FMA> ANSR" unit was the Alberta Natural
## Subregions clipped to THAT FMA, so a partner saw their own share of each subregion. v3's
## `buildReportingPolygons()` clips to `studyArea`, which is now an ecoregion group unioning
## several tenures, so the per-tenure breakdown was not produced at all.
##
## ORDERING MATTERS: cross first, group second. Both nrvtools consumers group features by NAME, so
## clipping ANSR to the whole group and labelling features "Central Mixedwood" pools every tenure's
## share of Central Mixedwood into one number -- the very reporting error this restores. Crossing
## first gives distinctly-named units ("Central Mixedwood <tenure>") that then group correctly, and
## the same name-grouping tallies a subregion's disjoint parts within a tenure for free.

utils::globalVariables("Name")

## Dissolve the tenure layer to ONE geometry per tenure short name. A tenure can be several
## features in the v10 layer; crossing against a single part would report only that part.
##
## Tenures contributing less than `min_area_km2` are dropped: a neighbouring tenure that merely
## grazes the study-area boundary leaves a sliver in the clipped layer, and crossing it would
## mint a full set of reporting units (one per sub-region layer, each with its own refCode,
## aggregate directory, and figures) over a handful of pixels. Mirrors the sliver threshold in
## `build_studyarea_crosswalk()`.
.tenureUnits <- function(tenure, min_area_km2 = 1) {
  if (!inherits(tenure, "SpatVector")) {
    tenure <- terra::vect(tenure)
  }
  tenure <- tenure[!is.na(tenure$Name), "Name"]
  if (nrow(tenure) == 0L) {
    return(list())
  }
  d <- terra::aggregate(tenure, by = "Name")
  d$Name <- as.character(d$Name) ## aggregate() can return a factor

  km2 <- terra::expanse(d, unit = "km")
  tiny <- km2 < min_area_km2
  if (any(tiny)) {
    message(
      "Not crossing ", sum(tiny), " sliver tenure(s) contributing < ", min_area_km2,
      " km^2 to this study area: ",
      paste(sprintf("%s [%.2f km2]", d$Name[tiny], km2[tiny]), collapse = "; ")
    )
    d <- d[!tiny, ]
  }
  if (nrow(d) == 0L) {
    return(list())
  }
  stats::setNames(lapply(seq_len(nrow(d)), function(i) d[i, "Name"]), d$Name)
}

## Collapse a tenure token repeated in a crossed name down to its last occurrence
## ("ACTIVE SprayLake Montane SprayLake" -> "ACTIVE Montane SprayLake"). The ACTIVE|PASSIVE
## reduction in `joinReportingPolygons()` normally prevents this in the triple crossings, but it
## only fires when the landbase status literally starts with "ACTIVE"/"PASSIVE" -- source coverages
## label the status differently (e.g. contributing/non-contributing), and a doubled tenure name
## would then be baked into a partner-facing label.
.dropRepeatedToken <- function(x, token) {
  vapply(as.character(x), function(s) {
    parts <- strsplit(s, " ", fixed = TRUE)[[1L]]
    hits <- which(parts == token)
    if (length(hits) > 1L) {
      parts <- parts[-utils::head(hits, -1L)]
    }
    paste(parts, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}

## Cross one sub-region layer with one tenure geometry, returning NULL when they do not overlap
## (or the layer has no labelled features left). Unlabelled features are dropped BEFORE the
## concatenation, which would otherwise mint a literal "NA <tenure>" reporting unit.
.crossOne <- function(lyr, tenureGeom, unit = "<unit>") {
  if (!inherits(lyr, "SpatVector")) {
    lyr <- terra::vect(lyr)
  }
  lyr <- lyr[!is.na(lyr$Name), "Name"]
  if (nrow(lyr) == 0L) {
    return(NULL)
  }
  ## A failed crossing means a MISSING reporting unit, which is the class of silent bug this
  ## whole exercise is about -- so name the unit that was lost rather than reporting a bare
  ## geometry error. (Two of these, GEOS `TopologyException`s from the former
  ## `sf::st_intersection()` engine, cost Dawson Creek TSA its Caribou units unnoticed.)
  z <- tryCatch(joinReportingPolygons(lyr, tenureGeom), error = function(e) {
    warning(
      "REPORTING UNIT MISSING: crossing '", unit, "' failed: ", conditionMessage(e),
      call. = FALSE, immediate. = TRUE
    )
    NULL
  })
  if (is.null(z) || nrow(z) == 0L) {
    return(NULL)
  }
  z
}

#' Build the tenure x sub-region crossed reporting polygons
#'
#' Restores the v2 per-tenure reporting units on top of v3's study-area-clipped layers.
#' Each tenure in the tenure layer (dissolved to one geometry per tenure, see
#' [reportingPolygonLayers()]`$isTenure`) is crossed --- via [joinReportingPolygons()] ---
#' with:
#'
#' 1. every sub-region layer marked `cross` in `layers` (Caribou ranges, Alberta Natural
#'    Subregions, BC biogeoclimatic zones, NWT/national ecoregions, national ecozones);
#' 2. its own active/passive landbase coverage(s), if any ([buildLandbasePolygons()]);
#' 3. and, where a tenure has both, the landbase **and** Caribou crossing together --- the
#'    triple (`tenure x landbase x Caribou`) v2 produced for the West Fraser CLS blocks.
#'
#' Products are named `"<tenure> <layer NAME_SHORT>"` (e.g. `"SprayLake ANSR"`,
#' `"Tolko_Vand_WF_SL LB_S17 Caribou"`), so each is a distinct `reportingPolygons` entry
#' with its own `refCode`, parquet aggregate directory, and figure folder. Their features
#' are named `"<sub-region> <tenure>"` by the concatenation in [joinReportingPolygons()].
#' Crossings that do not intersect are dropped.
#'
#' @param polys Named list of study-area-clipped reporting layers, keyed by `NAME_SHORT`
#'   (the output of [buildReportingPolygons()]); must contain the tenure layer.
#' @param landbases Named list of per-tenure landbase coverages keyed
#'   `"<tenure> <NAME_SHORT>"` (the output of [buildLandbasePolygons()]); may be empty.
#' @param layers Candidate-layer table (default [reportingPolygonLayers()]), supplying
#'   `isTenure` and `cross`.
#' @param members optional character vector of the study area's **member** tenure short
#'   names (a study-area group's `name_short`s from [studyAreaCrosswalk()]). When supplied,
#'   only members are crossed. Tenure polygons overlap in the v10 layer, so clipping to a
#'   group boundary also catches neighbouring tenures that belong to a *different* group ---
#'   and since the crosswalk assigns each tenure to exactly one group, crossing them here
#'   would report a fragment of a tenure that is reported in full elsewhere. `NULL` (the
#'   default) crosses every tenure present.
#' @param min_area_km2 numeric; tenures contributing less than this to the study area are
#'   not crossed (default `1`). A neighbouring tenure that merely grazes the study-area
#'   boundary would otherwise mint a full set of reporting units over a handful of pixels.
#'
#' @return A named `list` of crossed layers (`SpatVector`s), ready to be appended to the
#'   `reportingPolygons` input to `NRV_summary`. The tenure-crossed landbases are included
#'   here, so the raw `landbases` list should **not** also be appended.
#' @seealso [buildReportingPolygons()], [buildLandbasePolygons()], [joinReportingPolygons()]
#' @export
buildCrossedReportingPolygons <- function(
  polys,
  landbases = list(),
  layers = reportingPolygonLayers(),
  members = NULL,
  min_area_km2 = 1
) {
  tenureKey <- layers$NAME_SHORT[layers$isTenure][[1L]]
  if (is.null(polys[[tenureKey]])) {
    warning(
      "no tenure layer ('", tenureKey, "') among the reporting polygons; ",
      "no tenure x sub-region crossings will be produced.", call. = FALSE
    )
    return(list())
  }

  units <- .tenureUnits(polys[[tenureKey]], min_area_km2 = min_area_km2)

  ## restrict to this study area's own tenures: v10 tenure polygons OVERLAP, so a group
  ## boundary also clips neighbouring tenures assigned (by the crosswalk) to a different
  ## group. Crossing those would report a fragment of a tenure that is reported in full in
  ## its own group -- the same partner's numbers split across two study-area reports.
  if (!is.null(members)) {
    nonMembers <- setdiff(names(units), as.character(members))
    if (length(nonMembers)) {
      message(
        "Not crossing ", length(nonMembers), " non-member tenure(s) overlapping this study ",
        "area (each is reported in its own study-area group): ",
        paste(nonMembers, collapse = ", ")
      )
      units <- units[intersect(names(units), as.character(members))]
    }
  }
  if (!length(units)) {
    return(list())
  }

  crossKeys <- intersect(layers$NAME_SHORT[layers$cross], names(polys))
  out <- list()

  for (tn in names(units)) {
    geom <- units[[tn]]

    ## (1) tenure x sub-region layers
    for (k in crossKeys) {
      z <- .crossOne(polys[[k]], geom, unit = paste(tn, k))
      if (!is.null(z)) {
        out[[paste(tn, k)]] <- z
      }
    }

    ## (2) tenure x landbase -- keyed "<tenure> <lb>", so select this tenure's coverages.
    ## `names()` of an empty list is NULL, which startsWith() rejects.
    lbNames <- names(landbases)
    lbKeys <- if (is.null(lbNames)) character(0) else lbNames[startsWith(lbNames, paste0(tn, " "))]
    for (k in lbKeys) {
      z <- .crossOne(landbases[[k]], geom, unit = k)
      if (is.null(z)) {
        next
      }
      out[[k]] <- z

      ## (3) triple: tenure x landbase x Caribou. Both operands are already tenure-crossed, so
      ## joinReportingPolygons()'s ACTIVE|PASSIVE reduction is what keeps the tenure name from
      ## appearing twice in the combined name.
      caribouKey <- paste(tn, "Caribou")
      if (!is.null(out[[caribouKey]])) {
        zz <- .crossOne(z, out[[caribouKey]], unit = paste(k, "Caribou"))
        if (!is.null(zz)) {
          zz$Name <- .dropRepeatedToken(zz$Name, tn)
          out[[paste(k, "Caribou")]] <- zz
        }
      }
    }
  }

  out
}
