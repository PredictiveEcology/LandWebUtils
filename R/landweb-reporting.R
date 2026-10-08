## The groups LandWeb reports leading vegetation and large patches by, kept apart from what it
## simulates. LandR decides a pixel's leading type from the species codes in `cohortData`, summing
## biomass within a code only, so reporting merged groups of separately simulated species needs the
## codes recoded to their group first; relabelling a species-level map afterwards is not the same.

#' LandWeb reporting groups
#'
#' The species groups LandWeb reports leading vegetation and large patches by.
#'
#' @details
#' These are the species groups of the original LandWeb current-condition data (Silvacom, 2018),
#' as listed in the "Notes" of each "Percent Cover by" table of its data dictionary:
#' White Spruce (white, Engelmann and hybrid spruce), Black Spruce (black spruce and the larches),
#' Pine (all pines), Fir (the true firs) and Deciduous (all broadleaves). Douglas-fir, which the
#' dictionary puts under Fir, is reported on its own, as in the v2 Spray Lake runs. Western redcedar
#' and western hemlock, which the dictionary does not list, are reported as Fir: LandWeb simulates
#' them with the firs (see [landweb_species_map()]).
#'
#' `label` is the v2 output label of each group. `Type` is the group's conifer or broadleaf class,
#' which `LandR::vegTypeMapGenerator()` reads to find mixedwood stands (`mixedType = 2`); every
#' member of a group has the same `Type` in `LandR::sppEquivalencies_CA`. The colours are those of
#' the LandWeb figures, checked for colour-vision deficiency; see [landweb_report_colours()].
#'
#' @return A `data.table` with one row per group and columns `code`, `label`, `Type` and `colour`.
#'
#' @seealso [landweb_report_map()], [landweb_add_report_column()], [landweb_report_colours()]
#'
#' @export
#' @examples
#' landweb_report_groups()
landweb_report_groups <- function() {
  data.table::data.table(
    code = c("Wh_Spruce", "Bl_Spruce", "Pine", "Fir", "Doug_fir", "Decid"),
    label = c("Wh Spruce", "Bl Spruce", "Pine", "Fir", "Doug fir", "Decid"),
    Type = c("Conifer", "Conifer", "Conifer", "Conifer", "Conifer", "Deciduous"),
    colour = c("#1B8A5A", "#2A6FB0", "#E0A100", "#19A7A0", "#C77CB0", "#9BBF3B")
  )
}

#' Simulated species and groups to LandWeb reporting groups
#'
#' The reporting group ([landweb_report_groups()]) of every code LandWeb simulates: single species,
#' by their `LandR` code, and the merged species groups of [landweb_species_map()].
#'
#' @return A named character vector; names are simulated codes, values are reporting-group codes.
#'
#' @seealso [landweb_add_report_column()]
#'
#' @export
#' @examples
#' landweb_report_map()[c("Lari_lar", "Pinu_spp")] ## "Bl_Spruce" "Pine"
landweb_report_map <- function() {
  c(
    Abie_bal = "Fir",
    Abie_las = "Fir",
    Abie_spp = "Fir",
    Betu_pap = "Decid",
    Lari_lar = "Bl_Spruce",
    Lari_occ = "Bl_Spruce",
    Lari_spp = "Bl_Spruce",
    Pice_eng = "Wh_Spruce",
    Pice_eng_gla = "Wh_Spruce",
    Pice_gla = "Wh_Spruce",
    Pice_mar = "Bl_Spruce",
    Pinu_ban = "Pine",
    Pinu_con = "Pine",
    Pinu_spp = "Pine",
    Popu_bal = "Decid",
    Popu_spp = "Decid",
    Popu_tre = "Decid",
    Pseu_men = "Doug_fir",
    Thuj_pli = "Fir",
    Tsug_het = "Fir"
  )
}

#' Add the LandWeb reporting groups to a species equivalency table
#'
#' Adds a column giving each row's reporting group, looked up from the code it is simulated as.
#'
#' @param sppEquiv `data.table` species equivalency table with column `col`, e.g. from
#'   [landweb_sppEquiv()]. It is not modified.
#' @param col character; the column holding the simulated codes.
#' @param reportCol character; the name of the new column.
#'
#' @return A copy of `sppEquiv` with column `reportCol`. It is an error for a row's simulated code
#'   to have no reporting group, or for `reportCol` to exist already.
#'
#' @seealso [landweb_report_map()]
#'
#' @export
landweb_add_report_column <- function(sppEquiv, col = "LandWeb", reportCol = "LandWebReport") {
  if (!data.table::is.data.table(sppEquiv)) {
    stop("`sppEquiv` must be a data.table.", call. = FALSE)
  }
  if (!col %in% names(sppEquiv)) {
    stop("`sppEquiv` has no `", col, "` column.", call. = FALSE)
  }
  if (reportCol %in% names(sppEquiv)) {
    stop("`sppEquiv` already has a `", reportCol, "` column.", call. = FALSE)
  }
  codes <- sppEquiv[[col]]
  group <- unname(landweb_report_map()[codes])
  missing <- unique(codes[is.na(group) & !is.na(codes) & nzchar(codes)])
  if (length(missing)) {
    stop("No LandWeb reporting group for: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  out <- data.table::copy(sppEquiv)
  data.table::set(out, j = reportCol, value = group)
  out
}

#' Colours for the LandWeb reporting groups
#'
#' @param groups character; reporting-group codes to include ([landweb_report_groups()]).
#'
#' @return A named character vector of colours for `groups` and `"Mixed"`, in the order of
#'   [landweb_report_groups()]: the form `LandR::vegTypeMapGenerator()` takes as `colors`, whose
#'   names must be exactly the groups in its `sppEquiv` plus `"Mixed"`.
#'
#' @export
#' @examples
#' landweb_report_colours(c("Pine", "Decid"))
landweb_report_colours <- function(groups = landweb_report_groups()[["code"]]) {
  tab <- landweb_report_groups()
  unknown <- setdiff(groups, tab[["code"]])
  if (length(unknown)) {
    stop(
      "Unknown LandWeb reporting group(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  keep <- tab[["code"]] %in% groups
  c(stats::setNames(tab[["colour"]][keep], tab[["code"]][keep]), Mixed = "#7D5BA6")
}
