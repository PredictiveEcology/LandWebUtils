#' @keywords internal
.ageClassCutOffs <- c(0L, 40L, 80L, 120L) ## keep as integer

#' @keywords internal
.ageClasses <- c("Young", "Immature", "Mature", "Old")

#' Simulation timesteps for analyses
#'
#' @param period numeric vector of length 2 corresponding to the start and end times
#'               to use for analyses.
#'
#' @param interval numeric indicating the interval between timesteps for analyses
#'
#' @export
#' @return numeric vector of timesteps for which to run analyses
analysesOutputsTimes <- function(period, interval) {
  seq(period[1], period[2], by = interval)
}

#' Extract study area name from run name
#'
#' @param area Simulated area (i.e., run) name
#'
#' @export
cleanAreaName <- Vectorize(function(area) {
  strsplit(area, "_")[[1]] %>%
    grep("Dispersal|ROS", ., invert = TRUE, value = TRUE) %>%
    paste(., collapse = "_")
})

#' Find LandWeb simulation output file
#'
#' @param outputDir path to LandWeb output directory.
#' @param area character string giving study area with run name.
#' @param rep integer giving the replicate id, or character string in the form of `"rep01"`.
#'
#' @return path to the file
#' @export
findSimFile <- function(outputDir, area, rep) {
  if (is.numeric(rep)) {
    rep <- sprintf("rep%02d", as.integer(rep))
  }

  fsim <- file.path(outputDir, area, rep, "mySimOut_1000.qs")

  ## try alt/older names
  if (!file.exists(fsim)) {
    fsim <- file.path(outputDir, area, rep, "mySimOut_1000.rds")
  }

  if (!file.exists(fsim)) {
    fsim <- file.path(outputDir, area, rep, "mySimOut_year1000.rds")
  }

  stopifnot(file.exists(fsim))

  return(fsim)
}
