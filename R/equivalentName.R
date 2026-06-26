## Vendored from LandR (GPL-3) to drop the LandR dependency. These look up
## equivalent species names in a `sppEquivalencies`-style table supplied by the
## caller (no data is bundled). Refresh from LandR if the upstream logic changes.

equivalentNameAsList <- function(value, df, multi) {
  lapply(df, function(x) {
    if (isTRUE(multi)) {
      which(x %in% as.character(value))
    } else {
      match(as.character(value), x)
    }
  })
}

equivalentNameColumn <- function(value, df, column, multi = FALSE, searchColumn = NULL) {
  out <- equivalentNameAsList(value, df, multi)
  if (is.null(searchColumn)) {
    names(which.max(unlist(lapply(out, function(x) sum(!is.na(x))))))
  } else {
    searchColumn
  }
}

equivalentName <- function(value, df, column, multi = FALSE, searchColumn = NULL) {
  out <- equivalentNameAsList(value, df, multi)
  likelyMatch <- equivalentNameColumn(
    value, df, column,
    multi = multi, searchColumn = searchColumn
  )
  df[[column]][out[[likelyMatch]]]
}
