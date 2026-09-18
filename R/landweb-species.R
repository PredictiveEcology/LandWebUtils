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
#' @param sppEquiv `data.table` species equivalency table with a `SCANFI`
#'   column, normally `LandR::sppEquivalencies_CA`. It is not modified.
#'
#' @return A copy of `sppEquiv` restricted to the rows that map to a LandWeb
#'   species group, with a new `LandWeb` column. It is an error for no row to
#'   map (see [landweb_require_species()]).
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

  out <- data.table::copy(sppEquiv)
  groups <- unname(landweb_species_map()[out[["SCANFI"]]])
  data.table::set(out, j = "LandWeb", value = groups)

  labels <- .landweb_group_labels()
  for (grp in names(labels)) {
    rows <- which(out[["LandWeb"]] == grp)
    for (col in names(labels[[grp]])) {
      data.table::set(out, i = rows, j = col, value = labels[[grp]][[col]])
    }
  }

  out <- out[!is.na(out[["LandWeb"]]), ]
  landweb_require_species(out, "sppEquiv")
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
