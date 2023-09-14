utils::globalVariables(c(
  "FRI", "haBurned", "LTHFC", "studyArea", "time"
))

#' LandMine diagnostic plots
#'
#' - `landmine_plot_areaBurnedOverTime()` plots the area burned over time by LTHFC polygon;
#' - `landmine_plot_LTHFC()` produces a [rasterVis::levelplot()] a map of the LTHFC polygons;
#' - `landmine_plot_FRI()` plots
#'
#' @return a `ggplot2` or `rasterVis` object; invoked for side effect of creating plots.
#'
#' @name landmine_plots
#' @rdname landmine_plots
NULL

#' @param areaBurnedOverTime Summary `data.frame` of area burned over time, containing
#'                           the following columns:
#'                           `time` (numeric) gives the simulation time (year);
#'                           `haBurned` (numeric) gives the burned area in hectares;
#'                           `FRI` (numeric) identifies the fire return interval polygon.
#'
#' @export
#' @importFrom ggplot2 aes element_text geom_area ggplot theme
#' @rdname landmine_plots
landmine_plot_areaBurnedOverTime <- function(areaBurnedOverTime) {
  ggplot(areaBurnedOverTime, aes(x = time, y = haBurned, fill = FRI, ymin = 0)) +
    # geom_line(size = 1.5) +
    geom_area() +
    theme(legend.text = element_text(size = 6))
}

#' @param lthfc long-term historic fire cycle map (raster).
#'
#' @param studyAreaName study area name (character).
#'
#' @param ... additional arguments passed to [rasterVis::levelplot()].
#'
#' @export
#' @importFrom rasterVis levelplot PuOrTheme
#' @rdname landmine_plots
landmine_plot_LTHFC <- function(lthfc, studyAreaName, ...) {
  rasterVis::levelplot(
    lthfc,
    main = paste("Long-term historic fire cycle (LTHFC) map for", studyAreaName),
    margin = FALSE,
    par.settings = PuOrTheme,
    ...
  )
}

#' @param friSummary Summary `data.frame` of simulated fire return intervals (FRI) vs.
#'                   the long-term historic fire cycles (LTHFC), containing the following columns:
#'                   `simArea` (character) gives the study area name;
#'                   `LTHFC` (numeric) is the expected (i.e., historic) FRI;
#'                   `FRI` (numeric) is the simulated FRI.
#'
#' @export
#' @importFrom ggplot2 aes geom_abline geom_point ggplot ggtitle
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw xlab ylab
#' @rdname landmine_plots
landmine_plot_FRI <- function(friSummary) {
  studyAreaName <- unique(friSummary$simArea)
  ggplot(friSummary, aes(x = LTHFC, y = FRI, col = studyArea)) +
    geom_point() +
    xlab("Expected fire return interval (years)") +
    ylab("Simulated fire return interval (years)") +
    ggtitle(paste("Expected vs. simulated fire return intervals in", studyAreaName)) +
    theme_bw() +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    geom_abline(slope = 1, lty = "dotted")
}
