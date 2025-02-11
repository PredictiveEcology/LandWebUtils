utils::globalVariables(c("modification_time", "path"))

#' Find LandWeb simulation output files
#'
#' @param outputDir path to LandWeb output directory.
#'
#' @param rep integer giving the replicate id, or character string in the form of `"rep01"`.
#'
#' @param before date before which files are considered outdated (default "2024-01-01").
#'
#' @return path to the file
#'
#' @export
#' @rdname simfiles
findSimFile <- function(outputDir = NULL, rep = NULL) {
  stopifnot(!is.null(outputDir), !is.null(rep))

  if (is.numeric(rep)) {
    rep <- sprintf("rep%02d", as.integer(rep))
  }
  ## TODO: check other file prefixes
  fsim <- file.path(outputDir, rep, "mySimOut_1000.qs")

  ## try alt/older names
  if (!file.exists(fsim)) {
    fsim <- file.path(outputDir, rep, "mySimOut_1000.rds")
  }

  if (!file.exists(fsim)) {
    fsim <- file.path(outputDir, rep, "mySimOut_year1000.rds")
  }

  stopifnot(file.exists(fsim))

  return(fsim)
}

#' @examples
#' \dontrun{
#'   outDir <- file.path("~/GitHub/LandWeb/outputs/Tolko_AB_S_aspenDispersal_logROS")
#'   if (dir.exists(outdir)) {
#'     oldFiles <- findOldSimFiles(outDir)      ## search all reps
#'     oldFiles <- findOldSimFiles(outDir, 1:5) ## search specific reps
#'     # fs::file_delete(oldFiles) ## double check before deleting
#'   }
#' }
#'
#' @export
#' @rdname simfiles
findOldSimFiles <- function(outputDir = NULL, rep = NULL, before = "2024-01-01") {
  stopifnot(!is.null(outputDir))

  if (is.null(rep)) {
    outputDir <- fs::dir_ls(outputDir, regexp = "rep([0-9])+", recurse = FALSE, type = "directory")
  } else if (is.numeric(rep)) {
    outputDir <- file.path(outputDir, sprintf("rep%02d", as.integer(rep)))
  }

  f1 <- fs::dir_ls(outputDir, regexp = "vegTypeMap|rstTimeSinceFire|standAgeMap",
                   recurse = FALSE, type = "file") |>
    grep("year([0-9]){3}[.](grd|gri|tif)", x = _, value = TRUE) |>
    fs::file_info() |>
    dplyr::filter(as.Date(modification_time) < as.Date("2025-01-15")) |>
    dplyr::pull(path)

  f2 <- fs::dir_ls(outputDir, regexp = paste0(c(
    "Abie_sp[.]tif",
    "Pice_gla[.]tif",
    "Pice_mar[.]tif",
    "Pinu_sp[.]tif",
    "Popu_sp[.]tif",
    "CASFRIAbie_sp[.]tif",
    "CASFRIPice_gla[.]tif",
    "CASFRIPice_mar[.]tif",
    "CASFRIPinu_sp[.]tif",
    "CASFRIPopu_sp[.]tif",
    "CHECKSUMS[.]txt",
    "mySimOut[.](qs|rds)",
    "mySimOut_0[0-9]00[.](qs|rds)",
    "rstFlammable_year1000[.]grd",
    "rstFlammable_year1000[.]gri"
  ) , collapse = "|"), recurse = FALSE, type = "file") |>
    fs::file_info() |>
    dplyr::filter(as.Date(modification_time) < as.Date(!!before)) |>
    dplyr::pull(path)

  return(c(f1, f2))
}
