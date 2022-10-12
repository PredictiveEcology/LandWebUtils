#' @keywords internal
.ageClassCutOffs <- c(0L, 40L, 80L, 120L) ## keep as integer

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
