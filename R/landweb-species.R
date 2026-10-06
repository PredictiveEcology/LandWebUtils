utils::globalVariables(c("code", "cover", "nCodes", "total", "LandWeb"))

## The LandWeb species groups lifted out of LandWeb_preamble's `InitSpecies()`, so the mapping and
## its labels have one tested definition rather than an inline table reachable only through a
## preamble run. The LandWeb project's SCANFI cover summary kept a hand copy of the same map.

#' SCANFI species to LandWeb species groups
#'
#' The mapping from SCANFI species codes to the seven species groups LandWeb
#' simulates: `Abie_spp`, `Lari_spp`, `Pice_gla`, `Pice_mar`, `Pinu_spp`,
#' `Popu_spp` and `Pseu_men`.
#'
#' @details
#' Every SCANFI species found in the LandWeb study areas is mapped, except
#' `POPU_GRA`, whose SCANFI layer is unreliable.
#'
#' Western redcedar (`THUJ_PLI`) and western hemlock (`TSUG_HET`) merge into
#' `Abie_spp`, the closest shade-tolerant softwood, rather than being simulated as
#' distinct species. Neither is among the original Silvacom current-condition
#' species groups, and in Alberta they look like SCANFI over-attribution
#' (western hemlock is about 25% of the Spray Lake FMA, well outside its range).
#' Engelmann spruce and its hybrid with white spruce merge into `Pice_gla`, and
#' paper birch into `Popu_spp`.
#'
#' @return A named character vector; names are SCANFI codes, values are LandWeb
#'   species groups.
#'
#' @seealso [landweb_sppEquiv()], which adds these groups to a species
#'   equivalency table.
#'
#' @export
#' @examples
#' landweb_species_map()[["TSUG_HET"]] ## "Abie_spp"
landweb_species_map <- function() {
  c(
    ABIE_BAL = "Abie_spp",
    ABIE_LAS = "Abie_spp",
    BETU_PAP = "Popu_spp",
    LARI_LAR = "Lari_spp",
    LARI_OCC = "Lari_spp",
    PICE_ENG = "Pice_gla",
    PICE_ENG_GLA = "Pice_gla",
    PICE_GLA = "Pice_gla",
    PICE_MAR = "Pice_mar",
    PINU_BAN = "Pinu_spp",
    PINU_CON_CON = "Pinu_spp", ## shore pine (Pinus contorta var. contorta; coastal)
    PINU_CON_LAT = "Pinu_spp", ## lodgepole pine (Pinus contorta var. latifolia; interior)
    POPU_BAL = "Popu_spp",
    POPU_TRE = "Popu_spp",
    PSEU_MEN = "Pseu_men",
    PSEU_MEN_GLA = "Pseu_men",
    ## TODO: revisit -- confirm Abie_spp is the right target (vs. dropping them, or a per-study-area
    ## rule for FMAs nearer the BC coast where they may genuinely occur).
    THUJ_PLI = "Abie_spp",
    TSUG_HET = "Abie_spp"
  )
}

#' Add the LandWeb species groups to a species equivalency table
#'
#' Adds a `LandWeb` column to a `LandR`-style species equivalency table, giving
#' each SCANFI species its LandWeb species group, and relabels the groups that
#' merge several species.
#'
#' @details
#' Groups that merge species take a group label in `EN_generic_short`,
#' `EN_generic_full` and `Leading`: `Abie_spp`, `Lari_spp`, `Pice_gla`,
#' `Pinu_spp`, `Popu_spp` and `Pseu_men`. `Abie_spp` is labelled `"Fir"`,
#' although western redcedar and western hemlock are merged into it too (see
#' [landweb_species_map()]). `Pice_mar` holds one species, so its row keeps its
#' own labels.
#'
#' Consumers such as `Biomass_core`'s leading-vegetation maps look a group's
#' label up from its first row, so a group without a label would be named after
#' whichever of its species comes first.
#'
#' A row without a SCANFI code is kept when another row of the same species
#' (the same `LandR` code) maps, and joins that row's group. `LandR` lists some
#' species under a generic name as well as by variety, and gives only one of
#' them a SCANFI code: the generic *Pinus contorta* row has none, while var.
#' *latifolia* and var. *contorta* do. Generic *Pinus contorta* takes its group
#' and colour from var. *latifolia* (`PINU_CON_LAT`, interior lodgepole pine), not
#' from whichever variety comes first. `Biomass_speciesParameters` assigns PSP
#' trees to groups by their `Latin_full` in this table and discards trees whose
#' name is missing, and the NFI records lodgepole pine as *Pinus contorta* and
#' black cottonwood as *Populus trichocarpa*. Without these rows neither
#' contributed to the growth curves fitted for `Pinu_spp` and `Popu_spp`. A row
#' with a SCANFI code that LandWeb does not use (coastal Douglas-fir,
#' `PSEU_MEN_MEN`) stays out.
#'
#' These rows add no SCANFI layer, and they come after the rows mapped by SCANFI
#' code, so the first row of each group is unchanged. One with no `colorHex`
#' takes the colour of its species' mapped row: `LandR::sppColors()` uses the
#' table's colours only when every row has one.
#'
#' @param sppEquiv `data.table` species equivalency table with `SCANFI` and
#'   `LandR` columns, normally `LandR::sppEquivalencies_CA`. It is not modified.
#'
#' @return A copy of `sppEquiv` restricted to the rows that map to a LandWeb
#'   species group, by SCANFI code or through their species (see Details), with
#'   a new `LandWeb` column. Rows mapped by SCANFI code come first. It is an
#'   error for no row to map (see [landweb_require_species()]).
#'
#' @seealso [landweb_species_map()]
#'
#' @export
landweb_sppEquiv <- function(sppEquiv) {
  if (!data.table::is.data.table(sppEquiv)) {
    stop("`sppEquiv` must be a data.table.", call. = FALSE)
  }
  if (!"SCANFI" %in% names(sppEquiv)) {
    stop("`sppEquiv` has no `SCANFI` column to map species from.", call. = FALSE)
  }
  if (!"LandR" %in% names(sppEquiv)) {
    stop("`sppEquiv` has no `LandR` column to match species by.", call. = FALSE)
  }

  out <- data.table::copy(sppEquiv)
  groups <- unname(landweb_species_map()[out[["SCANFI"]]])
  mapped <- which(!is.na(groups))

  ## rows with no SCANFI code, for a species that has a mapped row (e.g. generic Pinus contorta)
  spp <- out[["LandR"]]
  noCode <- is.na(out[["SCANFI"]]) | !nzchar(out[["SCANFI"]])
  sameSpp <- which(noCode & !is.na(spp) & nzchar(spp) & spp %in% spp[mapped])
  donor <- mapped[match(spp[sameSpp], spp[mapped])]
  ## where a species has several mapped varieties, the one LandWeb means
  code <- unname(.landweb_generic_donor()[spp[sameSpp]])
  pref <- ifelse(is.na(code), NA_integer_, match(code, out[["SCANFI"]]))
  donor[!is.na(pref)] <- pref[!is.na(pref)]
  groups[sameSpp] <- groups[donor]
  data.table::set(out, j = "LandWeb", value = groups)

  if ("colorHex" %in% names(out)) {
    noColour <- is.na(out[["colorHex"]][sameSpp]) | !nzchar(out[["colorHex"]][sameSpp])
    data.table::set(
      out,
      i = sameSpp[noColour],
      j = "colorHex",
      value = out[["colorHex"]][donor[noColour]]
    )
  }

  labels <- .landweb_group_labels()
  for (grp in names(labels)) {
    rows <- which(out[["LandWeb"]] == grp)
    for (col in names(labels[[grp]])) {
      data.table::set(out, i = rows, j = col, value = labels[[grp]][[col]])
    }
  }

  out <- out[c(mapped, sameSpp), ]
  landweb_require_species(out, "sppEquiv")
}

## SCANFI code of the variety a generic (no-SCANFI-code) row stands for, by `LandR` code: the
## NFI records interior lodgepole pine as plain *Pinus contorta*.
.landweb_generic_donor <- function() {
  c(Pinu_con = "PINU_CON_LAT")
}

#' Total SCANFI cover of each member species in a study area
#'
#' Sums the per-species SCANFI cover layers that the speciesData stage writes for
#' a study area (`SCANFI_spsCC_<CODE>_<year>_...tif`), one value per SCANFI code.
#'
#' @param dir directory holding the `SCANFI_spsCC_*.tif` layers.
#'
#' @return A named numeric vector of total cover (summed cover percent over the
#'   layer's cells); names are SCANFI codes. It is an error for `dir` to hold no
#'   such layer.
#'
#' @seealso [landweb_dominant_sppEquiv()]
#'
#' @export
landweb_member_cover <- function(dir) {
  files <- list.files(dir, pattern = "^SCANFI_spsCC_.+_[0-9]{4}_.*\\.tif$", full.names = TRUE)
  if (!length(files)) {
    stop("No SCANFI_spsCC_*.tif species cover layers in `", dir, "`.", call. = FALSE)
  }
  codes <- sub("^SCANFI_spsCC_(.+?)_[0-9]{4}_.*$", "\\1", basename(files))
  cover <- vapply(files, function(f) {
    terra::global(terra::rast(f), "sum", na.rm = TRUE)[[1]]
  }, numeric(1))
  stats::setNames(unname(cover), codes)
}

#' Give each merged LandWeb species group its dominant member's traits
#'
#' Where a LandWeb species group merges several species, keeps only the trait
#' code (`LANDIS_traits`) of the member with the most cover, and blanks the
#' others, so that `LandR::prepSpeciesTable()` and `LandR::speciesTableUpdate()`
#' take the group's traits from that one member.
#'
#' @details
#' Both LandR functions give a group the minimum of each numeric trait across its
#' members' `LANDIS_traits` codes. For LandWeb's `Pice_gla` (white, Engelmann and
#' hybrid spruce) that meant Engelmann spruce's 30 m effective seed dispersal for a
#' group that is mostly white spruce (100 m). A member's cover is summed over the
#' rows that share its `LANDIS_traits` code (the hybrid spruce shares Engelmann
#' spruce's), so rows with the dominant code keep it. A group with one trait code,
#' or with no cover at all in the study area, is left unchanged. Rows without a
#' SCANFI code (see [landweb_sppEquiv()]) contribute no cover, and keep their code
#' only if it is the dominant one.
#'
#' @param sppEquiv `data.table` from [landweb_sppEquiv()], with `SCANFI`,
#'   `LANDIS_traits` and `LandWeb` columns. It is not modified.
#' @param cover named numeric vector of total cover by SCANFI code, as from
#'   [landweb_member_cover()].
#'
#' @return A copy of `sppEquiv` with the non-dominant members' `LANDIS_traits` set
#'   to `NA`. Attribute `"dominant"` is a `data.table` of each group's dominant
#'   trait code and its share of the group's cover.
#'
#' @export
landweb_dominant_sppEquiv <- function(sppEquiv, cover) {
  need <- c("SCANFI", "LANDIS_traits", "LandWeb")
  miss <- setdiff(need, names(sppEquiv))
  if (length(miss)) {
    stop("`sppEquiv` lacks column(s) ", paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  out <- data.table::copy(data.table::as.data.table(sppEquiv))
  rowCover <- unname(cover[out[["SCANFI"]]])
  rowCover[is.na(rowCover)] <- 0
  codeCover <- data.table::data.table(LandWeb = out[["LandWeb"]], code = out[["LANDIS_traits"]], cover = rowCover)
  codeCover <- codeCover[!is.na(code) & nzchar(code), list(cover = sum(cover)), by = c("LandWeb", "code")]
  dominant <- codeCover[, list(
    code = code[which.max(cover)],
    share = if (sum(cover) > 0) max(cover) / sum(cover) else NA_real_,
    nCodes = .N,
    total = sum(cover)
  ), by = "LandWeb"]
  dominant <- dominant[nCodes > 1 & total > 0]
  for (i in seq_len(nrow(dominant))) {
    rows <- which(out[["LandWeb"]] == dominant$LandWeb[i] &
                    !is.na(out[["LANDIS_traits"]]) & out[["LANDIS_traits"]] != dominant$code[i])
    data.table::set(out, i = rows, j = "LANDIS_traits", value = NA_character_)
  }
  data.table::setattr(out, "dominant", dominant[, list(LandWeb, code, share)])
  out
}

## Labels for the LandWeb groups that merge several species.
.landweb_group_labels <- function() {
  list(
    Abie_spp = list(
      EN_generic_full = "Fir",
      EN_generic_short = "Fir",
      Leading = "Fir leading"
    ),
    Lari_spp = list(
      EN_generic_full = "Western Larch & Tamarack",
      EN_generic_short = "Larch & Tamarack",
      Leading = "Larch & Tamarack leading"
    ),
    Pice_gla = list(
      EN_generic_full = "White & Engelmann's Spruce",
      EN_generic_short = "Whi & Eng Spr",
      Leading = "White & Engelmann's Spruce leading"
    ),
    Pinu_spp = list(
      EN_generic_short = "Pine",
      EN_generic_full = "Pine",
      Leading = "Pine leading"
    ),
    Popu_spp = list(
      EN_generic_full = "Deciduous",
      EN_generic_short = "Decid",
      Leading = "Deciduous leading"
    ),
    Pseu_men = list(
      EN_generic_full = "Douglas fir",
      EN_generic_short = "Doug fir",
      Leading = "Douglas fir leading"
    )
  )
}

#' Stop when a species input is empty
#'
#' Checks that a species input handed between LandWeb pipeline stages is not
#' empty, and returns it unchanged so the check can wrap the reference inline.
#'
#' @details
#' The upstream `Biomass_*` modules and `LandR` treat a study area with no tree
#' species as valid: an empty species table, `NULL` species layers or an empty
#' `cohortData` make them skip their work and finish without error. For
#' LandWeb every study area is forested, so an empty species input can only
#' mean a broken species mapping, and a run would otherwise "succeed" with no
#' vegetation dynamics. This turns that into an error at the stage boundary.
#'
#' @param x the input to check: a `data.frame` (e.g. `sppEquiv`, `cohortData`),
#'   a `SpatRaster` (e.g. `speciesLayers`), or a named list holding one of those
#'   as element `name` (e.g. the result of `SpaDES.targets::sim_objects()`).
#' @param name character; the object's name, used to find it in a list and in
#'   the error message.
#'
#' @return `x`, unchanged.
#'
#' @export
#' @examples
#' landweb_require_species(data.frame(species = "Pice_mar"), "sppEquiv")
landweb_require_species <- function(x, name) {
  obj <- if (is.list(x) && !is.data.frame(x)) x[[name]] else x
  empty <- is.null(obj) ||
    (inherits(obj, "SpatRaster") && terra::nlyr(obj) == 0L) ||
    (is.data.frame(obj) && nrow(obj) == 0L)
  if (empty) {
    stop(
      "`", name, "` is empty: no tree species reached this stage. LandWeb study areas are ",
      "forested, so this indicates a broken species mapping (e.g. `sppEquiv`/`sppEquivCol`), ",
      "which the upstream modules would otherwise treat as a valid no-species run.",
      call. = FALSE
    )
  }
  x
}
