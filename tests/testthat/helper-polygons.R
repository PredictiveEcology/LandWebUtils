## An s x s square polygon with its lower-left corner at (x0, y0).
square <- function(x0, y0, s) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
}
