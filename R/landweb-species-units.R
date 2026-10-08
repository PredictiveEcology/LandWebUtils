## Species-level LandWeb: each SCANFI species is simulated on its own where it is common enough in the
## study area, and the rest stay in today's merged groups. landweb_species_map() and landweb_sppEquiv()
## (the merged groups) are unchanged.

#' SCANFI species to LandWeb species, one code per species
#'
#' Like [landweb_sppEquiv()], but the simulation column `LandWeb` holds each species' own `LandR` code,
#' and the merged group goes to `LandWebGroup`. [landweb_species_units()] then decides, per study
#' area, which species are common enough to simulate on their own.
#'
#' @details
#' The rows are those [landweb_sppEquiv()] keeps, in the same order. Hybrid white x Engelmann spruce
#' keeps its own code, `Pice_eng_gla`; [landweb_species_units()] puts it with whichever parent has more
#' cover in the study area.
#'
#' Every row of a species carries the labels and colour of its first row, so LandR names it the same
#' way whichever row it looks up. Every row of *Pinus contorta* -- shore pine, var. *latifolia* and the generic row -- takes interior
#' lodgepole pine's (`PINU_CON_LAT`) labels, colour and `LANDIS_traits` code. LandR lists shore pine
#' first, so its labels would otherwise name the species, and shore pine's trait code has no rows for
#' the ecozones LandWeb reads. Western larch shares tamarack's colour in LandR and gets its own here,
#' since `LandR::sppColors()` falls back to a generated palette when two species share a colour.
#'
#' @param sppEquiv `data.table` species equivalency table with `SCANFI` and `LandR` columns, normally
#'   `LandR::sppEquivalencies_CA`. It is not modified.
#'
#' @return A copy of the rows of `sppEquiv` that map to a LandWeb species, with columns `LandWeb`
#'   (the species' `LandR` code), `LandWebGroup` (its group in [landweb_species_map()]) and
#'   `LandWebReport` (its reporting group, [landweb_report_map()]).
#'
#' @seealso [landweb_species_units()]
#'
#' @export
landweb_species_sppEquiv <- function(sppEquiv) {
  if (!data.table::is.data.table(sppEquiv)) {
    stop("`sppEquiv` must be a data.table.", call. = FALSE)
  }
  for (col in c("SCANFI", "LandR")) {
    if (!col %in% names(sppEquiv)) {
      stop("`sppEquiv` has no `", col, "` column.", call. = FALSE)
    }
  }
  out <- data.table::copy(sppEquiv)
  rows <- .landweb_mapped_rows(out)
  keep <- c(rows$mapped, rows$sameSpp)
  data.table::set(out, j = "LandWebGroup", value = rows$groups)
  data.table::set(
    out,
    j = "LandWeb",
    value = ifelse(is.na(rows$groups), NA_character_, out[["LandR"]])
  )

  ## one row speaks for a species listed under several: its first mapped row, or for Pinus contorta
  ## interior lodgepole pine's, whose trait code it also takes
  ref <- .landweb_generic_donor()
  labelCols <- intersect(c("EN_generic_short", "EN_generic_full", "Leading", "colorHex"), names(out))
  for (sp in unique(out[["LandR"]][keep])) {
    to <- keep[out[["LandR"]][keep] %in% sp]
    from <- if (sp %in% names(ref)) match(ref[[sp]], out[["SCANFI"]]) else to[1L]
    cols <- if (sp %in% names(ref)) c(labelCols, intersect("LANDIS_traits", names(out))) else labelCols
    for (col in cols) {
      data.table::set(out, i = to, j = col, value = out[[col]][from])
    }
  }
  if ("colorHex" %in% names(out)) {
    for (sp in names(.landweb_species_colours())) {
      data.table::set(
        out,
        i = keep[out[["LandR"]][keep] %in% sp],
        j = "colorHex",
        value = .landweb_species_colours()[[sp]]
      )
    }
  }

  out <- out[keep, ]
  out <- landweb_add_report_column(out, col = "LandWeb")
  landweb_require_species(out, "sppEquiv")
}

#' Which species LandWeb simulates on its own in a study area
#'
#' Splits a species out of its merged LandWeb group where it holds at least `threshold` percent of the
#' study area's tree cover; the group's other species stay together. Returns everything a study area's
#' simulation needs to use the result: its species table, colours, and which species layers make up
#' each simulated unit.
#'
#' @details
#' Shares are measured as `Biomass_borealDataPrep` will see the layers. A candidate unit's layer is the
#' sum of its species' layers, rounded; pixels whose total cover is at most `floor` are left out;
#' cells of a unit with at most `floor` cover are dropped; and each pixel is rescaled to 100. A unit's
#' share is then its mean percent of a pixel's tree cover over `studyArea`.
#'
#' The rule, for each merged group (`LandWebGroup`):
#' 1. hybrid spruce (`Pice_eng_gla`) joins whichever parent, white or Engelmann spruce, has more cover
#'    in `studyArea`, before the floor;
#' 2. a member with a share of at least `threshold` is simulated on its own;
#' 3. the members left form the group's remainder: one unit if its share is at least `threshold`,
#'    otherwise merged into the group's largest split-out member, and otherwise (none split out)
#'    dropped, its cover going to the other species of each pixel;
#' 4. species in `neverSplit` take no part in 2-3: they join the largest unit of their group at the
#'    end, or are dropped if their group has none.
#'
#' A unit is named by its species' code, or by its group's code when it holds several species. It
#' takes its labels and colour from that species or group, and the traits of its member with most
#' cover: the other members' `LANDIS_traits` codes are blanked (see [landweb_dominant_sppEquiv()]),
#' and hybrid rows take their parent's code.
#'
#' A species unit's growth curve is fitted on that species' permanent sample plot trees only; a group
#' unit's on all its members'. `Biomass_speciesParameters` assigns plot trees to units by `Latin_full`,
#' so the rows of a species unit's merged minor members (hybrid spruce included), and of species in
#' `neverSplit`, have it blanked: their trees still count in plot biomass but do not shape the curve.
#' On WesternAlbertaUpland, merged Engelmann and hybrid spruce plots from the montane ecozone made white
#' spruce's pooled curve unidentifiable.
#'
#' @param sppEquiv `data.table` from [landweb_species_sppEquiv()].
#' @param layers `SpatRaster` of percent cover, one layer per `LandWeb` code: the speciesData stage's
#'   `speciesLayers` for the study area.
#' @param studyArea `SpatVector` or `sf` polygon(s) the shares are measured over, normally the reporting
#'   study area; `NULL` uses every cell of `layers`.
#' @param threshold numeric; percent of the study area's tree cover a species (or a group's remainder)
#'   needs to be simulated on its own.
#' @param floor numeric; `Biomass_borealDataPrep`'s `minCoverThreshold`.
#' @param neverSplit character; species never simulated on their own (western redcedar and western
#'   hemlock, whose traits are unreviewed).
#'
#' @return A list:
#'   - `sppEquiv`: the study area's species table, `LandWeb` holding each row's unit; rows of dropped
#'     species are removed;
#'   - `sppColorVect`: named colours of the units and `"Mixed"`;
#'   - `sppColorVectReport`: named colours of the reporting groups present and `"Mixed"`;
#'   - `layerMap`: named character vector, each layer of `layers` to its unit (`NA`: dropped), for
#'     [landweb_sum_layers()];
#'   - `units`: `data.table` of each unit's members, dominant member and share (percent);
#'   - `hybridParent`: the species hybrid spruce joined.
#'
#' @seealso [landweb_sum_layers()]
#'
#' @export
landweb_species_units <- function(
  sppEquiv,
  layers,
  studyArea = NULL,
  threshold = 1,
  floor = 5,
  neverSplit = c("Thuj_pli", "Tsug_het")
) {
  need <- c("LandWeb", "LandWebGroup", "LandWebReport", "LANDIS_traits")
  miss <- setdiff(need, names(sppEquiv))
  if (length(miss)) {
    stop(
      "`sppEquiv` lacks column(s) ",
      paste(miss, collapse = ", "),
      "; see landweb_species_sppEquiv().",
      call. = FALSE
    )
  }
  spp <- names(layers)
  unknown <- setdiff(spp, sppEquiv[["LandWeb"]])
  if (length(unknown)) {
    stop(
      "Layer(s) with no species in `sppEquiv`: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  ## cover matrix over the study area: pixels x species layers
  if (!is.null(studyArea)) {
    sa <- if (inherits(studyArea, "SpatVector")) studyArea else terra::vect(studyArea)
    if (!terra::same.crs(sa, layers)) {
      sa <- terra::project(sa, terra::crs(layers))
    }
    layers <- terra::mask(terra::crop(layers, sa, snap = "out"), sa)
  }
  m <- terra::values(layers, mat = TRUE)
  m <- m[rowSums(!is.na(m)) > 0, , drop = FALSE]
  m[is.na(m)] <- 0
  colnames(m) <- spp

  ## 1. hybrid spruce joins its dominant parent
  hybrid <- "Pice_eng_gla"
  hybridParent <- NA_character_
  layerSpecies <- stats::setNames(spp, spp)
  if (hybrid %in% spp) {
    parents <- intersect(c("Pice_gla", "Pice_eng"), spp)
    hybridParent <- if (length(parents)) {
      parents[which.max(colSums(m[, parents, drop = FALSE]))]
    } else {
      "Pice_gla"
    }
    if (!hybridParent %in% spp) {
      m <- cbind(m, stats::setNames(matrix(0, nrow(m), 1L), hybridParent))
      colnames(m)[ncol(m)] <- hybridParent
    }
    m[, hybridParent] <- m[, hybridParent] + m[, hybrid]
    m <- m[, colnames(m) != hybrid, drop = FALSE]
    layerSpecies[[hybrid]] <- hybridParent
  }
  species <- colnames(m)
  cover <- colSums(m)

  groupOf <- unique(sppEquiv[, c("LandWeb", "LandWebGroup"), with = FALSE])
  groupOf <- stats::setNames(groupOf[["LandWebGroup"]], groupOf[["LandWeb"]])
  rule <- .landweb_unit_rule(
    m,
    groupOf[species],
    threshold = threshold,
    floor = floor,
    neverSplit = neverSplit
  )
  unitOf <- rule$unitOf ## species -> unit; NA = dropped
  shares <- .landweb_floored_shares(
    m,
    split(names(unitOf)[!is.na(unitOf)], unitOf[!is.na(unitOf)]),
    floor
  )

  ## dominant member of each unit (most cover; species in neverSplit only if nothing else)
  units <- data.table::data.table(
    species = names(unitOf),
    unit = unname(unitOf),
    cover = cover[names(unitOf)]
  )
  units <- units[!is.na(units[["unit"]])]
  units[, "eligible" := !(species %in% neverSplit)]
  dom <- units[,
    list(
      members = paste(species, collapse = "+"),
      dominant = if (any(eligible)) {
        species[eligible][which.max(cover[eligible])]
      } else {
        species[which.max(cover)]
      }
    ),
    by = "unit"
  ]
  dom[, "share" := unname(shares[unit])]
  dom[, "kind" := unname(rule$kind[unit])]
  dom[, "group" := unname(groupOf[dominant])]

  ## the study area's species table
  eq <- data.table::copy(data.table::as.data.table(sppEquiv))
  rowSpecies <- unname(ifelse(eq[["LandWeb"]] == hybrid, hybridParent, eq[["LandWeb"]]))
  rowUnit <- unname(unitOf[rowSpecies])
  ## a hybrid row takes its parent's trait code
  isHybrid <- eq[["LandWeb"]] %in% hybrid
  if (any(isHybrid)) {
    parentCode <- eq[["LANDIS_traits"]][match(hybridParent, eq[["LandWeb"]])]
    data.table::set(eq, i = which(isHybrid), j = "LANDIS_traits", value = parentCode)
  }
  ## traits of the unit's dominant member only
  domOf <- stats::setNames(dom[["dominant"]], dom[["unit"]])
  notDominant <- !is.na(rowUnit) & rowSpecies != domOf[rowUnit]
  data.table::set(eq, i = which(notDominant), j = "LANDIS_traits", value = NA_character_)
  ## whose plot trees shape the unit's growth curve: Biomass_speciesParameters assigns PSP trees to a
  ## unit by `Latin_full`, and counts trees it cannot assign only as plot biomass. A species unit is
  ## fitted on its own species' trees, a group unit on all its members'; merged minor members of a
  ## species unit (hybrid spruce included) and species in `neverSplit` are never fitted.
  kindOf <- stats::setNames(dom[["kind"]], dom[["unit"]])
  ownSpecies <- eq[["LandWeb"]] ## before the hybrid joins its parent
  notFitted <- !is.na(rowUnit) &
    (ownSpecies %in% neverSplit | (kindOf[rowUnit] %in% "species" & ownSpecies != rowUnit))
  if ("Latin_full" %in% names(eq)) {
    data.table::set(eq, i = which(notFitted), j = "Latin_full", value = NA_character_)
  }
  data.table::set(eq, j = "LandWeb", value = rowUnit)
  eq <- eq[!is.na(rowUnit)]

  ## labels and colour of each unit, on every one of its rows
  labs <- .landweb_unit_labels(sppEquiv, dom)
  for (col in intersect(
    c("EN_generic_short", "EN_generic_full", "Leading", "colorHex"),
    names(eq)
  )) {
    data.table::set(eq, j = col, value = labs[[col]][match(eq[["LandWeb"]], labs[["unit"]])])
  }
  cols <- stats::setNames(labs[["colorHex"]], labs[["unit"]])
  if (anyDuplicated(cols)) {
    stop(
      "Two LandWeb units share a colour: ",
      paste(names(cols)[duplicated(cols) | duplicated(cols, fromLast = TRUE)], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    sppEquiv = landweb_require_species(eq, "sppEquiv"),
    sppColorVect = c(cols[order(names(cols))], Mixed = "#7D5BA6"),
    sppColorVectReport = landweb_report_colours(unique(eq[["LandWebReport"]])),
    layerMap = stats::setNames(unname(unitOf[layerSpecies[spp]]), spp),
    units = dom[order(-dom[["share"]])],
    hybridParent = hybridParent
  )
}

#' Sum species cover layers into the units LandWeb simulates
#'
#' @param layers `SpatRaster` of percent cover, one layer per species (named by its `LandWeb` code).
#' @param layerMap named character vector from [landweb_species_units()]: each layer to its unit; `NA`
#'   drops the layer.
#'
#' @return A `SpatRaster` with one layer per unit, the sum of its species' layers, named by unit. Cells
#'   that are `NA` in every layer stay `NA`. Summing these layers is exactly what `Biomass_speciesData`
#'   does when it merges species (it resamples and thresholds each SCANFI layer, then sums them, with no
#'   cap or rescaling).
#'
#' @export
landweb_sum_layers <- function(layers, layerMap) {
  nm <- names(layers)
  miss <- setdiff(nm, names(layerMap))
  if (length(miss)) {
    stop(
      "No unit given in `layerMap` for layer(s): ",
      paste(miss, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  map <- layerMap[nm]
  units <- unique(map[!is.na(map)])
  if (!length(units)) {
    stop("`layerMap` drops every layer.", call. = FALSE)
  }
  noCover <- terra::allNA(layers) ## cells NA in every layer stay NA
  out <- terra::rast(lapply(units, function(u) {
    ids <- which(map %in% u)
    r <- if (length(ids) == 1L) layers[[ids]] else sum(layers[[ids]], na.rm = TRUE)
    terra::mask(r, noCover, maskvalues = 1)
  }))
  names(out) <- units
  out
}

## Shares (percent) of each unit in `units` (named list of species -> unit members) after
## Biomass_borealDataPrep's floor: unit layers summed and rounded; pixels with total cover <= floor left
## out; unit cells <= floor dropped; pixels rescaled; mean over pixels.
.landweb_floored_shares <- function(m, units, floor) {
  u <- vapply(units, function(ss) round(rowSums(m[, ss, drop = FALSE])), numeric(nrow(m)))
  if (!is.matrix(u)) {
    u <- matrix(u, ncol = length(units), dimnames = list(NULL, names(units)))
  }
  u <- u[rowSums(u) > floor, , drop = FALSE]
  u[u <= floor] <- 0
  rs <- rowSums(u)
  u <- u[rs > 0, , drop = FALSE] / rs[rs > 0]
  if (!nrow(u)) {
    return(stats::setNames(rep(0, length(units)), names(units)))
  }
  100 * colSums(u) / nrow(u)
}

## The split rule of landweb_species_units() on a cover matrix; returns `unitOf`, species -> unit.
.landweb_unit_rule <- function(m, groupOf, threshold, floor, neverSplit) {
  species <- colnames(m)
  alone <- .landweb_floored_shares(m, stats::setNames(as.list(species), species), floor)
  units <- list()
  kind <- character(0) ## "species": one species, maybe with minor members merged in; "group": several
  for (g in unique(groupOf[setdiff(species, neverSplit)])) {
    mem <- setdiff(species[groupOf[species] == g], neverSplit)
    split <- mem[alone[mem] >= threshold]
    rest <- setdiff(mem, split)
    for (s in split) {
      units[[s]] <- s
      kind[[s]] <- "species"
    }
    if (length(rest)) {
      restName <- if (length(rest) == 1L) rest else g
      others <- setdiff(species, c(unlist(units), rest))
      trial <- c(
        units,
        stats::setNames(list(rest), restName),
        stats::setNames(as.list(others), others)
      )
      if (.landweb_floored_shares(m, trial, floor)[[restName]] >= threshold) {
        units[[restName]] <- rest
        kind[[restName]] <- if (length(rest) == 1L) "species" else "group"
      } else if (length(split)) {
        big <- split[which.max(alone[split])]
        units[[big]] <- c(units[[big]], rest)
      }
    }
  }
  for (s in intersect(neverSplit, species)) {
    inGroup <- names(units)[vapply(units, function(u) any(groupOf[u] == groupOf[[s]]), logical(1))]
    if (length(inGroup)) {
      big <- inGroup[which.max(vapply(units[inGroup], function(u) sum(alone[u]), numeric(1)))]
      units[[big]] <- c(units[[big]], s)
    }
  }
  unitOf <- stats::setNames(rep(NA_character_, length(species)), species)
  for (u in names(units)) {
    unitOf[units[[u]]] <- u
  }
  list(unitOf = unitOf, kind = kind)
}

## Labels and colour of each unit: a species unit takes its species' (from its first row in `sppEquiv`,
## as landweb_species_sppEquiv() fixed them); a group unit takes its group's.
.landweb_unit_labels <- function(sppEquiv, dom) {
  first <- sppEquiv[!duplicated(sppEquiv[["LandWeb"]])]
  grpLabels <- .landweb_group_labels()
  grpCols <- .landweb_group_colours()
  data.table::rbindlist(lapply(seq_len(nrow(dom)), function(i) {
    u <- dom[["unit"]][i]
    if (identical(dom[["kind"]][i], "species")) {
      r <- first[first[["LandWeb"]] == u]
      data.table::data.table(
        unit = u,
        EN_generic_short = r[["EN_generic_short"]],
        EN_generic_full = r[["EN_generic_full"]],
        Leading = r[["Leading"]],
        colorHex = r[["colorHex"]]
      )
    } else {
      gl <- grpLabels[[u]]
      ## a group some of whose species run on their own is labelled as the rest of it
      partial <- any(dom[["group"]][-i] == dom[["group"]][i])
      suffix <- if (partial) " (other)" else ""
      data.table::data.table(
        unit = u,
        EN_generic_short = paste0(gl[["EN_generic_short"]], suffix),
        EN_generic_full = paste0(gl[["EN_generic_full"]], suffix),
        Leading = paste0(gl[["EN_generic_full"]], suffix, " leading"),
        colorHex = grpCols[[u]]
      )
    }
  }))
}

## Colours for species whose LandR colour another LandWeb species shares.
.landweb_species_colours <- function() {
  c(Lari_occ = "#C2691E")
}

## Colours for units that merge several species, from the LandWeb figures' palette.
.landweb_group_colours <- function() {
  c(
    Abie_spp = "#19A7A0",
    Lari_spp = "#D8602E",
    Pice_gla = "#1B8A5A",
    Pinu_spp = "#E0A100",
    Popu_spp = "#9BBF3B"
  )
}
