## Grouping a simulation's species, for reporting and for fire. LandR::vegTypeMapGenerator() decides a
## pixel's leading type from the species codes in `cohortData`, summing biomass within a code only, so
## a map of groups needs the cohorts recoded to their groups first.

#' Group the simulated species of a species table
#'
#' Pairs each simulated species code with its group, and builds the one-row-per-group table that
#' `LandR::vegTypeMapGenerator()` takes as `sppEquiv` for cohorts recoded to their groups
#' ([recode_cohorts()]).
#'
#' @details
#' It is an error for a simulated code to have no group, or more than one, and for a group to hold
#' both conifers and broadleaves (`Type`). `LandR::vegTypeMapGenerator()` reads `Type`, one row per
#' code, to find mixedwood stands (`mixedType = 2`), so a group of both would be counted twice there.
#'
#' @param sppEquiv `data.table` species equivalency table with columns `sppEquivCol`, `groupCol` and
#'   `Type`. It is not modified.
#' @param sppEquivCol character; the column of simulated species codes.
#' @param groupCol character; the column of groups.
#'
#' @return A list: `map`, a named character vector from simulated code to group; `groups`, a
#'   `data.table` with one row per group (sorted) and columns `groupCol` and `Type`; and `col`,
#'   `groupCol`.
#'
#' @seealso [recode_cohorts()], [landweb_add_report_column()]
#'
#' @export
sppEquiv_groups <- function(sppEquiv, sppEquivCol, groupCol) {
  need <- c(sppEquivCol, groupCol, "Type")
  miss <- setdiff(need, names(sppEquiv))
  if (length(miss)) {
    stop("`sppEquiv` lacks column(s) ", paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  eq <- unique(data.table::as.data.table(sppEquiv)[, need, with = FALSE])
  eq <- eq[!is.na(eq[[sppEquivCol]]) & nzchar(eq[[sppEquivCol]])]
  pairs <- unique(eq[, c(sppEquivCol, groupCol), with = FALSE])
  noGroup <- pairs[[sppEquivCol]][is.na(pairs[[groupCol]]) | !nzchar(pairs[[groupCol]])]
  if (length(noGroup)) {
    stop(
      "No group in `",
      groupCol,
      "` for: ",
      paste(unique(noGroup), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  twoGroups <- unique(pairs[[sppEquivCol]][duplicated(pairs[[sppEquivCol]])])
  if (length(twoGroups)) {
    stop(
      "Simulated code(s) in more than one group: ",
      paste(twoGroups, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  typed <- !is.na(eq[["Type"]]) & nzchar(eq[["Type"]])
  types <- unique(eq[typed, c(groupCol, "Type"), with = FALSE])
  mixed <- unique(types[[groupCol]][duplicated(types[[groupCol]])])
  if (length(mixed)) {
    stop(
      "Group(s) holding both conifers and broadleaves: ",
      paste(mixed, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  untyped <- setdiff(unique(pairs[[groupCol]]), types[[groupCol]])
  if (length(untyped)) {
    stop("No `Type` for group(s): ", paste(untyped, collapse = ", "), ".", call. = FALSE)
  }
  list(
    map = stats::setNames(pairs[[groupCol]], pairs[[sppEquivCol]]),
    groups = types[order(types[[groupCol]])],
    col = groupCol
  )
}

#' Recode cohorts' species to their groups
#'
#' @param cohortData `data.table` of cohorts with columns `pixelGroup`, `speciesCode` and `B`. It is
#'   not modified.
#' @param map named character vector from species code to group, as `sppEquiv_groups()$map`.
#'
#' @return A new `data.table` with columns `pixelGroup`, `speciesCode` (now the group) and `B`, ready
#'   for `LandR::vegTypeMapGenerator()` with `sppEquiv_groups()$groups`. It is an error for a
#'   cohort's species to have no group.
#'
#' @seealso [sppEquiv_groups()]
#'
#' @export
recode_cohorts <- function(cohortData, map) {
  cd <- data.table::as.data.table(cohortData)[, c("pixelGroup", "speciesCode", "B"), with = FALSE]
  code <- as.character(cd[["speciesCode"]])
  group <- unname(map[code])
  if (anyNA(group)) {
    stop(
      "No group for species: ",
      paste(unique(code[is.na(group)]), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  data.table::set(cd, j = "speciesCode", value = group)
  cd
}
