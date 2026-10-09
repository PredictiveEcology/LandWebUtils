## `LandR::sppEquivalencies_CA` (LandR 1.2.0.9047), as in test-landweb-species.R.
lr_table_growth <- function() {
  data.table::fread(
    test_path("fixtures", "sppEquivalencies_CA_LandR-1.2.0.9047.csv"),
    colClasses = "character",
    na.strings = NULL
  )
}

## the shared growth-curve stage's species table, trimmed to the columns that matter here
growth_species_table <- function() {
  data.table::data.table(
    species = c("Pice_gla", "Pinu_con", "Lari_lar"),
    growthcurve = c(0.68, 0.74, 0.7),
    mortalityshape = c(21L, 23L, 22L),
    mANPPproportion = c(5.123, 5.502, 5.34),
    inflationFactor = c(1.006, 1.009, 1.027),
    longevity = c(400L, 338L, 350L),
    growthCurveSource = c("estimated", "estimated", "imputed"),
    shadetolerance = c(2, 1, 1)
  )
}

square <- function(x0, y0, s) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
}

test_that("the growth-curve species table counts hybrid spruce as Engelmann spruce", {
  eq <- landweb_species_sppEquiv(lr_table_growth())
  g <- landweb_growth_sppEquiv(eq)
  hyb <- g[["Latin_full"]] %in% "Picea engelmannii x glauca"
  eng <- g[["Latin_full"]] %in% "Picea engelmannii"
  expect_true(any(hyb))
  ## Biomass_speciesParameters relabels BC white spruce records only if the hybrid's LandR code is Pice_eng
  expect_identical(unique(g[["LandR"]][hyb]), "Pice_eng")
  expect_identical(unique(g[["LandWeb"]][hyb]), "Pice_eng")
  expect_identical(unique(g[["LANDIS_traits"]][hyb]), unique(g[["LANDIS_traits"]][eng]))
  expect_identical(g[!hyb], data.table::as.data.table(eq)[!hyb])
  expect_identical(unique(eq[["LandR"]][hyb]), "Pice_eng_gla") ## the input is left as it was
})

test_that("growth traits keep the five traits and their source, and refuse missing ones", {
  tr <- landweb_growth_traits(growth_species_table())
  expect_named(
    tr,
    c(
      "species",
      "growthcurve",
      "mortalityshape",
      "mANPPproportion",
      "inflationFactor",
      "longevity",
      "growthCurveSource"
    )
  )
  bad <- growth_species_table()
  bad$inflationFactor[2] <- NA
  expect_error(landweb_growth_traits(bad), "Pinu_con")
  expect_error(landweb_growth_traits(growth_species_table()[, !"longevity"]), "longevity")
})

test_that("each unit takes its dominant species' growth traits, with their source", {
  tr <- landweb_growth_traits(growth_species_table())
  units <- data.table::data.table(
    unit = c("Pice_gla", "Pinu_spp", "Lari_lar"),
    dominant = c("Pice_gla", "Pinu_con", "Lari_lar")
  )
  sp <- data.table::data.table(
    species = c("Pinu_spp", "Pice_gla", "Lari_lar"),
    longevity = c(150L, 250L, 350L),
    growthcurve = c(0, 0, 0),
    mortalityshape = c(15L, 15L, 15L),
    hardsoft = "soft"
  )
  out <- landweb_unit_growth_traits(sp, tr, units)
  expect_identical(out$species, sp$species)
  expect_identical(out$longevity, c(338L, 400L, 350L))
  expect_identical(out$mortalityshape, c(23L, 21L, 22L))
  expect_identical(out$growthcurve, c(0.74, 0.68, 0.7))
  expect_identical(out$mANPPproportion, c(5.502, 5.123, 5.34))
  expect_identical(out$inflationFactor, c(1.009, 1.006, 1.027))
  expect_identical(
    out$growthTraitSource,
    c("estimated (Pinu_con)", "estimated (Pice_gla)", "imputed (Lari_lar)")
  )
  expect_identical(out$hardsoft, sp$hardsoft)
  ## the input is left as it was
  expect_identical(sp$longevity, c(150L, 250L, 350L))
  expect_false("inflationFactor" %in% names(sp))
})

test_that("a unit without a dominant species, or whose dominant has no traits, stops", {
  tr <- landweb_growth_traits(growth_species_table())
  sp <- data.table::data.table(species = c("Pice_gla", "Abie_spp"))
  expect_error(
    landweb_unit_growth_traits(
      sp,
      tr,
      data.table::data.table(unit = "Pice_gla", dominant = "Pice_gla")
    ),
    "Abie_spp"
  )
  units <- data.table::data.table(
    unit = c("Pice_gla", "Abie_spp"),
    dominant = c("Pice_gla", "Abie_bal")
  )
  expect_error(landweb_unit_growth_traits(sp, tr, units), "Abie_bal")
})

test_that("a merged unit takes the curve of its member with most cover that was fitted", {
  firs <- data.table::data.table(
    species = c("Abie_bal", "Abie_las"),
    growthcurve = c(0.72, 0.68),
    mortalityshape = c(23L, 21L),
    mANPPproportion = c(5.48, 4.915),
    inflationFactor = c(1.012, 1.041),
    longevity = c(200L, 236L),
    growthCurveSource = c("imputed", "estimated")
  )
  tr <- data.table::rbindlist(list(landweb_growth_traits(growth_species_table()), firs))
  units <- data.table::data.table(
    unit = c("Abie_spp", "Lari_lar"),
    dominant = c("Abie_bal", "Lari_lar"),
    ranked = list(c("Abie_bal", "Abie_las", "Thuj_pli"), "Lari_lar")
  )
  sp <- data.table::data.table(species = c("Abie_spp", "Lari_lar"))
  out <- landweb_unit_growth_traits(sp, tr, units)
  ## balsam fir leads by cover but was not fitted; tamarack has no fitted member, so keeps its mean
  expect_identical(out$growthTraitSource, c("estimated (Abie_las)", "imputed (Lari_lar)"))
  expect_identical(out$longevity, c(236L, 350L))
  ## without the ranking, the dominant member's traits
  expect_identical(
    landweb_unit_growth_traits(sp, tr, units[, !"ranked"])$growthTraitSource,
    c("imputed (Abie_bal)", "imputed (Lari_lar)")
  )
})

test_that("damage agent codes: bark beetles and defoliators, by plot source", {
  both <- landweb_damage_codes()
  expect_named(both, c("BC", "AB", "SK", "NFI"))
  expect_identical(both$NFI, c("IB", "ID"))
  expect_identical(both$AB, c(3L, 1L, 2L)) ## mountain pine beetle; spruce budworm, defoliator
  expect_identical(both$SK, 3L) ## death by insects, once
  expect_identical(landweb_damage_codes("defoliators")$SK, 3L)
  expect_true(all(c("IBM", "IBS", "IDE", "IDW", "IDX", "CHX") %in% both$BC))
  expect_true(all(grepl("^I[BD]|^CHX$", both$BC)))
  expect_identical(anyDuplicated(both$BC), 0L)
  beetles <- landweb_damage_codes("barkBeetles")
  expect_identical(beetles$NFI, "IB")
  expect_true(all(startsWith(beetles$BC, "IB")))
  expect_error(landweb_damage_codes("fire"))
})

test_that("the LandWeb area is the outline of the fire-cycle polygons, without holes", {
  holed <- sf::st_difference(square(0, 0, 10), square(4, 4, 2))
  lthfc <- sf::st_sf(
    LTHFC = c(100, 60),
    geometry = sf::st_sfc(holed, square(10, 0, 10), crs = 3400)
  )
  a <- landweb_area(lthfc)
  expect_s3_class(a, "sfc")
  expect_length(a, 1L)
  expect_equal(as.numeric(sf::st_area(a)), 200)
})

test_that("the fitting area keeps whole the ecological units that touch the area", {
  eco <- sf::st_sf(
    ECOPROVINCE = c("A", "B", "C"),
    geometry = sf::st_sfc(square(0, 0, 10), square(10, 0, 10), square(30, 0, 10), crs = 3400)
  )
  local_mocked_bindings(prepInputs = function(...) eco, .package = "reproducible")
  area <- sf::st_transform(sf::st_sfc(square(5, 2, 2), crs = 3400), 3978)
  out <- landweb_anpp_area(area, "ecoprovince", destinationPath = withr::local_tempdir())
  expect_identical(out$ECOPROVINCE, "A")
  expect_equal(sf::st_crs(out), sf::st_crs(area))
  expect_equal(
    as.numeric(sf::st_area(out)),
    as.numeric(sf::st_area(sf::st_transform(eco[1, ], 3978)))
  )
})

test_that("no BC plot in the fitting area: no BEC zones are fetched", {
  skip_if_not_installed("bcdata")
  gis <- sf::st_sf(
    OrigPlotID1 = c("p1", "p2"),
    geometry = sf::st_sfc(sf::st_point(c(1, 1)), sf::st_point(c(50, 50)), crs = 3400)
  )
  meas <- data.table::data.table(OrigPlotID1 = c("p1", "p2"), source = c("NFI", "BC"))
  expect_null(landweb_bec_zones(gis, meas, sf::st_sfc(square(0, 0, 10), crs = 3400)))
})
