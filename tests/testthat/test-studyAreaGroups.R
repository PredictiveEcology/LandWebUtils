mk_sq <- function(x0, x1, y0, y1) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
}

test_that("build_studyarea_crosswalk groups tenures by dominant ecoregion", {
  crs <- 3400 # NAD83 / Alberta 10-TM (metres)
  eco <- sf::st_sf(
    REGION_NAM = c("Region A", "Region B"),
    geometry = sf::st_sfc(mk_sq(0, 2e5, 0, 2e5), mk_sq(2e5, 4e5, 0, 2e5), crs = crs)
  )
  fmas <- sf::st_sf(
    FMA_NAME = c("FMA One", "FMA Two", NA, NA),
    FMU_NAME = NA_character_,
    TSA_NUMB_1 = NA_character_,
    FOREST_NAM = c(NA, NA, "Forest Three", NA),
    FML_NAME = NA_character_,
    Name = c(NA, NA, NA, "Tiny Sliver"),
    LICENSEE_N = c("Co X", "Co X", NA, NA),
    geometry = sf::st_sfc(
      mk_sq(2e4, 6e4, 2e4, 6e4), # Region A, 1600 km^2
      mk_sq(1.2e5, 1.6e5, 2e4, 6e4), # Region A, 1600 km^2
      mk_sq(2.2e5, 2.6e5, 2e4, 6e4), # Region B, 1600 km^2
      mk_sq(3e5, 3.05e5, 2e4, 2.5e4), # Region B, 25 km^2 -> dropped
      crs = crs
    )
  )
  cw <- suppressMessages(build_studyarea_crosswalk(fmas, eco, min_area_km2 = 100))

  expect_named(cw, c("group", "fma_name", "company", "province", "eco_unit", "area_km2", "mpix"))
  expect_setequal(cw$group[cw$fma_name %in% c("FMA One", "FMA Two")], "RegionA")
  expect_equal(cw$group[cw$fma_name == "Forest Three"], "RegionB")
  expect_equal(cw$province[cw$fma_name == "FMA One"], "AB")
  expect_equal(cw$province[cw$fma_name == "Forest Three"], "SK")
  expect_equal(cw$company[cw$fma_name == "FMA One"], "Co X")
  expect_equal(cw$company[cw$fma_name == "Forest Three"], "Forest Three")
})

test_that("sub-threshold tenures are dropped as slivers", {
  crs <- 3400
  eco <- sf::st_sf(
    REGION_NAM = "Region A",
    geometry = sf::st_sfc(mk_sq(0, 2e5, 0, 2e5), crs = crs)
  )
  fmas <- sf::st_sf(
    FMA_NAME = c("Big FMA", "Sliver FMA"),
    FMU_NAME = NA_character_, TSA_NUMB_1 = NA_character_, FOREST_NAM = NA_character_,
    FML_NAME = NA_character_, Name = NA_character_, LICENSEE_N = NA_character_,
    geometry = sf::st_sfc(
      mk_sq(2e4, 6e4, 2e4, 6e4), # 1600 km^2
      mk_sq(1e5, 1.05e5, 1e5, 1.05e5), # 25 km^2
      crs = crs
    )
  )
  cw <- suppressMessages(build_studyarea_crosswalk(fmas, eco, min_area_km2 = 100))
  expect_equal(cw$fma_name, "Big FMA")
})

test_that("build_studyarea_crosswalk errors on a missing eco_field", {
  crs <- 3400
  eco <- sf::st_sf(FOO = "x", geometry = sf::st_sfc(mk_sq(0, 1e5, 0, 1e5), crs = crs))
  fmas <- sf::st_sf(FMA_NAME = "A", geometry = sf::st_sfc(mk_sq(0, 1e5, 0, 1e5), crs = crs))
  expect_snapshot(build_studyarea_crosswalk(fmas, eco), error = TRUE)
})

test_that(".extractStudyAreaGroup unions a group's member FMAs", {
  crs <- 3400
  fmas <- sf::st_sf(
    FMA_NAME = c("FMA One", "FMA Two", "FMA Three"),
    FMU_NAME = NA_character_, TSA_NUMB_1 = NA_character_, FOREST_NAM = NA_character_,
    FML_NAME = NA_character_, Name = NA_character_, LICENSEE_N = NA_character_,
    geometry = sf::st_sfc(
      mk_sq(0, 4e4, 0, 4e4), # GroupA
      mk_sq(4e4, 8e4, 0, 4e4), # GroupA (adjacent)
      mk_sq(2e5, 2.4e5, 0, 4e4), # GroupB
      crs = crs
    )
  )
  cw <- data.frame(
    group = c("GroupA", "GroupA", "GroupB"),
    fma_name = c("FMA One", "FMA Two", "FMA Three"),
    stringsAsFactors = FALSE
  )
  sa <- .extractStudyAreaGroup(fmas, "GroupA", cw)
  expect_s3_class(sa, "sf")
  expect_equal(nrow(sa), 1L)
  expect_equal(sa$studyAreaName, "GroupA")
  expect_equal(round(as.numeric(sf::st_area(sa)) / 1e6), 3200) # 40km x 80km union
})

test_that("prepStudyArea rejects an unknown study area", {
  cw <- data.frame(group = "GroupA", fma_name = "FMA One", stringsAsFactors = FALSE)
  expect_snapshot(
    prepStudyArea("ZZZNonExistent", destinationPath = tempdir(), crosswalk = cw),
    error = TRUE
  )
})
