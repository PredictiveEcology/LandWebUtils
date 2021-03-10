#' @keywords internal
.ageClassCutOffs <- c(0, 40, 80, 120)

#' @keywords internal
.ageClasses <- c("Young", "Immature", "Mature", "Old")

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
