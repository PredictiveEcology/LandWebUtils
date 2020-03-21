#' @keywords internal
.ageClassCutOffs <- c(0, 40, 80, 120)

#' @keywords internal
.ageClasses <- c("Young", "Immature", "Mature", "Old")

#' Generate simulation file name
#'
#' Assists with saving/retrieving LandWeb simulations.
#'
#' @param name Object name (e.g., \code{"mySimOut"})
#' @param path Directory location in where the file will be located (e.g., an \code{outputPath}).
#' @param time Optional simulation time to use as filename suffix. Default \code{NULL}.
#' @param ext  The file extension to use (default \code{"rds"}).
#' @export
#' @importFrom SpaDES.core paddedFloatToChar
#' @importFrom reproducible normPath
simFile <- function(name, path, time = NULL, ext = "rds") {
  if (is.null(time))
    file.path(normPath(path), paste0(name, ".", ext))
  else {
    file.path(normPath(path), paste0(name, "_", paddedFloatToChar(time, padL = 4), ".", ext))
  }
}
