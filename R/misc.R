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
