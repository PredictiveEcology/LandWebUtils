#' Plot boundary polygon(s) for forest management areas
#'
#' @param x        `sf` object corresponding to the FMA to be plotted
#'
#' @param provs    `sf` object corresponding to the provincial (or territorial) boundaries to plot
#'
#' @param caribou  Optional `sf` object corresponding to caribou boundaries
#'
#' @param xsr      Optional `sf` object corresponding to a buffered `studyArea`
#'
#' @param title    Character string to use for plot title
#'
#' @param png      Optional. If non-NULL, must be a valid file path to write a png
#'
#' @export
#' @rdname plotFMA
plotFMA <- function(x, provs, caribou = NULL, xsr = NULL, title = NULL, png = NULL) {
  provs <- sf::st_transform(provs, x)

  ## regular boring old plot
  if (!is.null(png)) png(filename = png, width = 1200, height = 800)
  plot(provs)
  plot(x[, "Name"], main = title, col = "lightblue", add = TRUE)
  if (!is.null(caribou)) plot(caribou, col = "magenta", add = TRUE)
  if (!is.null(xsr)) plot(xsr, add = TRUE)
  if (!is.null(png)) dev.off()

  ## sexy ggplot version
  x_gg <- plotGG(x, provs, caribou, png, title)

  if (!is.null(png)) {
    png2 <- gsub("[.]png", "_gg.png", png)
    ggsave(png2, x_gg, width = 12, height = 8) ## a bit slow...
  }
}

#' @export
#' @rdname plotFMA
plotLandWeb <- function(x, provs, caribou = NULL, xsr = NULL, title = NULL, png = NULL) {
  provs <- sf::st_transform(provs, x)

  ## regular boring old plot
  if (!is.null(png)) png(filename = png, width = 1800, height = 1200)
  plot(provs)
  plot(x, main = title, col = "lightblue", add = TRUE)
  if (!is.null(caribou)) plot(caribou, col = "magenta", add = TRUE)
  if (!is.null(xsr)) plot(xsr, add = TRUE)
  if (!is.null(png)) dev.off()

  ## sexy ggplot version
  x_gg <- plotGG(x, provs, caribou, png, title)

  if (!is.null(png)) {
    png2 <- gsub("[.]png", "_gg.png", png)
    ggsave(png2, x_gg, width = 18, height = 12) ## a bit slow...
  }
}

#' @export
#' @rdname plotFMA
plotGG <- function(x, provs, caribou = NULL, png = NULL, title = NULL) {
  x.sf <- sf::st_as_sf(x)
  provs.sf <- sf::st_as_sf(provs)

  x_gg <- ggplot(provs.sf) +
    geom_sf() +
    geom_sf(data = x.sf, color = "white", fill = scales::hue_pal()(16)[11]) +
    coord_sf() +
    theme_bw()

  if (!is.null(title)) {
    x_gg <- x_gg + ggtitle(title)
  }

  if (!is.null(caribou)) {
    caribou.sf <- sf::st_as_sf(caribou)

    x_gg <- x_gg + geom_sf(data = caribou.sf, color = "white", fill = scales::hue_pal()(16)[15])
  }

  return(x_gg)
}
