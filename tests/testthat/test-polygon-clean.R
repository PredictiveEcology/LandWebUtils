lthfc_sf <- function(fri, col = "LTHFC") {
  polys <- lapply(seq_along(fri), function(i) {
    sf::st_polygon(list(rbind(c(i, 0), c(i + 1, 0), c(i + 1, 1), c(i, 1), c(i, 0))))
  })
  x <- sf::st_sf(geometry = sf::st_sfc(polys, crs = 3979))
  x[[col]] <- fri
  x
}

test_that(".cleanLandWebStudyArea() drops intervals strictly below minFRI and keeps minFRI itself", {
  out <- .cleanLandWebStudyArea(lthfc_sf(c(10, 39, 40, 41, 170)), minFRI = 40)
  expect_identical(out$fireReturnInterval, c(NA, NA, 40, 41, 170))
})

test_that(".cleanLandWebStudyArea() defaults to minFRI = 40", {
  out <- .cleanLandWebStudyArea(lthfc_sf(c(35, 40, 45)))
  expect_identical(out$fireReturnInterval, c(NA, 40, 45))
})

test_that(".cleanLandWebStudyArea() accepts the LTHRC spelling", {
  out <- .cleanLandWebStudyArea(lthfc_sf(c(20, 100), col = "LTHRC"), minFRI = 25)
  expect_identical(out$fireReturnInterval, c(NA, 100))
  expect_false("LTHRC" %in% names(out))
})

test_that(".cleanLandWebStudyArea() requires a fire-return-interval column", {
  expect_snapshot(error = TRUE, .cleanLandWebStudyArea(lthfc_sf(1, col = "other")))
})

test_that(".cleanLandWebStudyArea() works on a SpatVector, as LandWeb_preamble passes it", {
  skip_if_not_installed("tidyterra")
  loadNamespace("tidyterra") ## registers the dplyr methods for SpatVector
  v <- terra::vect(lthfc_sf(c(25, 40, 55)))
  out <- .cleanLandWebStudyArea(v, minFRI = 40)
  expect_s4_class(out, "SpatVector")
  expect_identical(out$fireReturnInterval, c(NA, 40, 55))
})

test_that("polygonClean() dispatches on type and rejects what it does not know", {
  x <- lthfc_sf(c(30, 60))
  expect_identical(
    polygonClean(x, type = "LandWeb", minFRI = 50)$fireReturnInterval,
    c(NA, 60)
  )
  expect_snapshot(error = TRUE, polygonClean(x))
  expect_snapshot(error = TRUE, polygonClean(x, type = "other"))
})
