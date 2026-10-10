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

test_that("the growth-curve species table counts hybrid spruce as Engelmann spruce", {
  eq <- landweb_species_sppEquiv(lr_table_growth())
  g <- landweb_growth_sppEquiv(eq)
  hyb <- g[["Latin_full"]] %in% "Picea engelmannii x glauca"
  eng <- g[["Latin_full"]] %in% "Picea engelmannii"
  expect_identical(sum(hyb), 1L)
  ## Biomass_speciesParameters relabels BC white spruce records only if the hybrid's LandR code is Pice_eng
  expect_identical(g[["LandR"]][hyb], "Pice_eng")
  expect_identical(g[["LandWeb"]][hyb], "Pice_eng")
  expect_identical(g[["LANDIS_traits"]][hyb], unique(g[["LANDIS_traits"]][eng]))
  expect_identical(g[!hyb], data.table::as.data.table(eq)[!hyb])
  expect_identical(eq[["LandR"]][hyb], "Pice_eng_gla") ## the input is left as it was
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
  expect_snapshot(error = TRUE, landweb_growth_traits(bad))
  expect_snapshot(error = TRUE, landweb_growth_traits(growth_species_table()[, !"longevity"]))
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
  expect_named(sp, c("species", "longevity", "growthcurve", "mortalityshape", "hardsoft"))
  expect_identical(sp$longevity, c(150L, 250L, 350L))
})

test_that("a unit without a dominant species, or whose dominant has no traits, stops", {
  tr <- landweb_growth_traits(growth_species_table())
  sp <- data.table::data.table(species = c("Pice_gla", "Abie_spp"))
  expect_snapshot(error = TRUE, {
    landweb_unit_growth_traits(
      sp,
      tr,
      data.table::data.table(unit = "Pice_gla", dominant = "Pice_gla")
    )
  })
  units <- data.table::data.table(
    unit = c("Pice_gla", "Abie_spp"),
    dominant = c("Pice_gla", "Abie_bal")
  )
  expect_snapshot(error = TRUE, landweb_unit_growth_traits(sp, tr, units))
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
  expect_contains(both$BC, c("IBM", "IBS", "IDE", "IDW", "IDX", "CHX"))
  expect_match(both$BC, "^I[BD]|^CHX$")
  expect_identical(anyDuplicated(both$BC), 0L)
  beetles <- landweb_damage_codes("barkBeetles")
  expect_identical(beetles$NFI, "IB")
  expect_match(beetles$BC, "^IB")
  expect_snapshot(error = TRUE, landweb_damage_codes("fire"))
})

## plots p1-p4 inside a 10-degree square, p5 outside; two measurements per plot, three trees each
psp_fixture <- function() {
  ids <- paste0("p", 1:5)
  plot <- data.table::data.table(
    OrigPlotID1 = rep(ids, each = 2),
    MeasureID = paste0(rep(ids, each = 2), c("_m1", "_m2")),
    MeasureYear = rep(c(2000, 2010), 5),
    source = "NFI"
  )
  measure <- plot[rep(seq_len(nrow(plot)), each = 3)]
  measure[["TreeNumber"]] <- rep(1:3, nrow(plot))
  gis <- sf::st_as_sf(
    data.table::data.table(
      OrigPlotID1 = ids,
      baseSA = 50,
      x = c(1, 2, 3, 4, 20),
      y = c(1, 2, 3, 4, 20)
    ),
    coords = c("x", "y"),
    crs = 4326
  )
  list(PSPmeasure = measure, PSPplot = plot, PSPgis = gis)
}
psp_area <- function() sf::st_transform(sf::st_sfc(square(0, 0, 10), crs = 4326), 3857)

test_that("a bootstrap draws as many whole plots as lie in the fitting area, and only those", {
  out <- landweb_resample_psp(psp_fixture(), psp_area(), seed = 1L)
  drawn <- sub("_r[0-9]+$", "", out$PSPgis$OrigPlotID1)
  expect_length(drawn, 4L)
  expect_in(drawn, paste0("p", 1:4))
  expect_identical(anyDuplicated(out$PSPgis$OrigPlotID1), 0L)
  expect_identical(anyDuplicated(out$PSPplot$MeasureID), 0L)
  expect_identical(nrow(out$PSPplot), 8L)
  expect_identical(nrow(out$PSPmeasure), 24L)
  expect_setequal(out$PSPmeasure$OrigPlotID1, out$PSPgis$OrigPlotID1)
  expect_setequal(out$PSPplot$OrigPlotID1, out$PSPgis$OrigPlotID1)
  ## Biomass_speciesParameters subsets and keys the locations as a data.table
  expect_s3_class(out$PSPgis, "sf")
  expect_s3_class(out$PSPgis, "data.table")
  area <- sf::st_as_sf(sf::st_transform(psp_area(), 4326))
  expect_identical(nrow(data.table::setkeyv(out$PSPgis[area, ], "OrigPlotID1")), 4L)
})

test_that("a plot drawn twice becomes two plots, each with all its measurements", {
  psp <- psp_fixture()
  for (s in 1:50) {
    out <- landweb_resample_psp(psp, psp_area(), seed = s)
    counts <- table(sub("_r[0-9]+$", "", out$PSPgis$OrigPlotID1))
    if (any(counts == 2L)) break
  }
  twice <- names(counts)[counts == 2L][1L]
  copies <- out$PSPplot[startsWith(out$PSPplot$OrigPlotID1, paste0(twice, "_r"))]
  expect_setequal(copies$OrigPlotID1, paste0(twice, c("_r1", "_r2")))
  expect_setequal(copies$MeasureID, paste0(twice, c("_m1_r1", "_m2_r1", "_m1_r2", "_m2_r2")))
})

test_that("a bootstrap is reproducible from its seed and leaves the global seed alone", {
  withr::local_seed(42)
  before <- get(".Random.seed", envir = globalenv())
  a <- landweb_resample_psp(psp_fixture(), psp_area(), seed = 7L)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(landweb_resample_psp(psp_fixture(), psp_area(), seed = 7L), a)
})

test_that("a bootstrap stops when no plot lies in the fitting area", {
  far <- sf::st_sfc(square(50, 50, 1), crs = 4326)
  expect_snapshot(error = TRUE, landweb_resample_psp(psp_fixture(), far, seed = 1L))
})

test_that("trait frequencies count each species' trait sets over the refits", {
  tr <- landweb_growth_traits(growth_species_table())
  alt <- data.table::copy(tr)
  alt$growthcurve[alt$species == "Pice_gla"] <- 0.8
  refits <- data.table::rbindlist(list(
    cbind(resample = 1L, tr),
    cbind(resample = 2L, tr),
    cbind(resample = 3L, alt)
  ))
  out <- landweb_growth_trait_frequency(refits, tr)
  expect_identical(out$species, c("Lari_lar", "Pice_gla", "Pice_gla", "Pinu_con"))
  expect_identical(out$growthcurve[out$species == "Pice_gla"], c(0.68, 0.8))
  expect_identical(out$n, c(3L, 2L, 1L, 3L))
  expect_equal(out$share, c(3, 2, 1, 3) / 3)
  expect_identical(out$fullFit, c(TRUE, TRUE, FALSE, TRUE))
  expect_snapshot(error = TRUE, landweb_growth_trait_frequency(refits[, !"resample"], tr))
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
