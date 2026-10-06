## Fixtures -----------------------------------------------------------------------------------------

## `LandR::sppEquivalencies_CA` (LandR 1.2.0.9047, the version LandWeb pins): every row, with the
## columns `landweb_sppEquiv()` reads or relabels. Missing values are "", as in LandR's own table.
lr_sppEquiv <- function() {
  data.table::fread(
    test_path("fixtures", "sppEquivalencies_CA_LandR-1.2.0.9047.csv"),
    colClasses = "character", na.strings = NULL
  )
}

## VERBATIM from LandWeb_preamble 1.0.7 InitSpecies(), minus its `if (FALSE)` exploration block and
## the `sim$` assignments. This is the behaviour the promoted function must reproduce exactly. The
## module ran it on the lazy-loaded LandR table, so the copy stands in for that.
sppEquiv_inline_preamble_1.0.7 <- function(sppEquiv) {
  sppEquiv <- data.table::copy(sppEquiv)

  ## Make LandWeb spp equivalencies
  sppEquiv[,
    LandWeb := c(
      ABIE_BAL = "Abie_spp",
      ABIE_LAS = "Abie_spp",
      BETU_PAP = "Popu_spp",
      LARI_LAR = "Lari_spp",
      LARI_OCC = "Lari_spp",
      PICE_ENG = "Pice_gla",
      PICE_ENG_GLA = "Pice_gla", ## TODO: confirm merge with Pice_gla
      PICE_GLA = "Pice_gla",
      PICE_MAR = "Pice_mar",
      PINU_BAN = "Pinu_spp",
      PINU_CON_CON = "Pinu_spp", ## shore pine (Pinus contorta var. contorta; coastal)
      PINU_CON_LAT = "Pinu_spp", ## lodgepole pine (Pinus contorta var. latifolia; interior)
      POPU_BAL = "Popu_spp",
      POPU_TRE = "Popu_spp",
      PSEU_MEN = "Pseu_men",
      PSEU_MEN_GLA = "Pseu_men",
      THUJ_PLI = "Abie_spp",
      TSUG_HET = "Abie_spp"
    )[SCANFI]
  ]

  sppEquiv[
    LandWeb == "Lari_spp",
    `:=`(
      EN_generic_full = "Western Larch & Tamarack",
      EN_generic_short = "Larch & Tamarack",
      Leading = "Larch & Tamarack leading"
    )
  ]

  sppEquiv[
    LandWeb == "Pice_gla",
    `:=`(
      EN_generic_full = "White & Engelmann's Spruce",
      EN_generic_short = "Whi & Eng Spr",
      Leading = "White & Engelmann's Spruce leading"
    )
  ]

  sppEquiv[
    grep("Pin", LandWeb),
    `:=`(
      EN_generic_short = "Pine",
      EN_generic_full = "Pine",
      Leading = "Pine leading"
    )
  ]

  sppEquiv[
    LandWeb == "Popu_spp",
    `:=`(
      EN_generic_full = "Deciduous",
      EN_generic_short = "Decid",
      Leading = "Deciduous leading"
    )
  ]

  sppEquiv[
    LandWeb == "Pseu_men",
    `:=`(
      EN_generic_full = "Douglas fir",
      EN_generic_short = "Doug fir",
      Leading = "Douglas fir leading"
    )
  ]

  sppEquiv[!is.na(LandWeb), ]
}

## landweb_species_map ------------------------------------------------------------------------------

test_that("landweb_species_map() maps each SCANFI code once, onto the seven LandWeb groups", {
  map <- landweb_species_map()
  expect_identical(anyDuplicated(names(map)), 0L)
  expect_setequal(
    unname(map),
    c("Abie_spp", "Lari_spp", "Pice_gla", "Pice_mar", "Pinu_spp", "Popu_spp", "Pseu_men")
  )
  expect_false("POPU_GRA" %in% names(map))
})

test_that("every LandWeb group has a LandMine fuel type", {
  expect_in(unname(landweb_species_map()), names(landmine_known_species()))
})

test_that("every mapped SCANFI code exists in LandR's species table", {
  expect_in(names(landweb_species_map()), lr_sppEquiv()$SCANFI)
})

test_that("each species with a mapped SCANFI code maps to one LandWeb group", {
  ## landweb_sppEquiv() gives a species' rows without a SCANFI code the group of its first mapped row
  map <- landweb_species_map()
  coded <- lr_sppEquiv()[SCANFI %in% names(map)][, LandWeb := map[SCANFI]]
  expect_identical(unique(coded[, .(n = data.table::uniqueN(LandWeb)), by = LandR]$n), 1L)
})

## landweb_sppEquiv ---------------------------------------------------------------------------------

test_that("landweb_sppEquiv() reproduces the preamble's inline table, plus the Abie_spp label", {
  expected <- sppEquiv_inline_preamble_1.0.7(lr_sppEquiv())
  ## the one deliberate change: 1.0.7 gave Abie_spp no group label
  expected[LandWeb == "Abie_spp", `:=`(
    EN_generic_full = "Fir",
    EN_generic_short = "Fir",
    Leading = "Fir leading"
  )]
  ## the rows mapped by SCANFI code come first; those after them are matched by species
  expect_identical(
    as.data.frame(landweb_sppEquiv(lr_sppEquiv())[seq_len(nrow(expected))]),
    as.data.frame(expected)
  )
})

test_that("landweb_sppEquiv() keeps one row per SCANFI species, and that species' other names", {
  out <- landweb_sppEquiv(lr_sppEquiv())
  coded <- out[nzchar(SCANFI)]
  expect_setequal(coded$SCANFI, names(landweb_species_map()))
  expect_identical(coded$LandWeb, unname(landweb_species_map()[coded$SCANFI]))
  expect_identical(
    as.data.frame(out[!nzchar(SCANFI), .(Latin_full, LandR, LandWeb)]),
    data.frame(
      Latin_full = c("Pinus contorta", "Populus balsamifera v. balsamifera", "Populus trichocarpa"),
      LandR = c("Pinu_con", "Popu_bal", "Popu_bal"),
      LandWeb = c("Pinu_spp", "Popu_spp", "Popu_spp")
    )
  )
})

test_that("landweb_sppEquiv() maps the NFI names for lodgepole pine and black cottonwood", {
  ## Biomass_speciesParameters drops PSP trees whose `Latin_full` is not in sppEquiv
  out <- landweb_sppEquiv(lr_sppEquiv())
  expect_identical(out[Latin_full == "Pinus contorta", LandWeb], "Pinu_spp")
  expect_identical(out[Latin_full == "Populus trichocarpa", LandWeb], "Popu_spp")
  expect_identical(out[Latin_full == "Pinus contorta", Leading], "Pine leading")
  ## not a LandWeb species; a SCANFI layer LandWeb does not use
  expect_false(any(startsWith(out$Latin_full, "Fraxinus")))
  expect_false("Pseudotsuga menziesii var. menziesii" %in% out$Latin_full)
})

test_that("rows matched by species leave each group's first row and the colours as they were", {
  out <- landweb_sppEquiv(lr_sppEquiv())
  expect_true(all(nzchar(out[!duplicated(LandWeb), SCANFI])))
  ## LandR::sppColors() uses the table's colours only when no row lacks one
  expect_false(any(is.na(out$colorHex) | !nzchar(out$colorHex)))
  expect_identical(
    out[Latin_full == "Pinus contorta", colorHex],
    out[SCANFI == "PINU_CON_LAT", colorHex]
  )
})

test_that("landweb_sppEquiv() treats an NA SCANFI code or colour like an empty one", {
  input <- lr_sppEquiv()
  for (col in names(input)) {
    data.table::set(input, i = which(input[[col]] == ""), j = col, value = NA_character_)
  }
  out <- landweb_sppEquiv(input)
  ref <- landweb_sppEquiv(lr_sppEquiv())
  expect_identical(out$Latin_full, ref$Latin_full)
  expect_identical(out$LandWeb, ref$LandWeb)
  expect_identical(out$colorHex, ref$colorHex)
})

test_that("landweb_sppEquiv() gives each merged group a single label", {
  out <- landweb_sppEquiv(lr_sppEquiv())
  merged <- out[LandWeb %in% c("Abie_spp", "Lari_spp", "Pice_gla", "Pinu_spp", "Popu_spp", "Pseu_men")]
  nLabels <- merged[, .(
    short = data.table::uniqueN(EN_generic_short),
    full = data.table::uniqueN(EN_generic_full),
    leading = data.table::uniqueN(Leading)
  ), by = LandWeb]
  expect_identical(unique(c(nLabels$short, nLabels$full, nLabels$leading)), 1L)
  expect_identical(out[LandWeb == "Pinu_spp", unique(Leading)], "Pine leading")
})

test_that("landweb_sppEquiv() labels Abie_spp as fir, cedar and hemlock included", {
  out <- landweb_sppEquiv(lr_sppEquiv())
  expect_setequal(out[LandWeb == "Abie_spp", SCANFI], c("ABIE_BAL", "ABIE_LAS", "THUJ_PLI", "TSUG_HET"))
  expect_identical(out[LandWeb == "Abie_spp", unique(Leading)], "Fir leading")
})

test_that("every group label maps back to exactly one LandWeb group", {
  out <- landweb_sppEquiv(lr_sppEquiv())
  for (col in c("EN_generic_short", "EN_generic_full", "Leading")) {
    groupsPerLabel <- out[, .(n = data.table::uniqueN(LandWeb)), by = col]$n
    expect_identical(unique(groupsPerLabel), 1L)
  }
})

test_that("landweb_sppEquiv() does not modify its input", {
  input <- lr_sppEquiv()
  before <- data.table::copy(input)
  landweb_sppEquiv(input)
  expect_identical(as.data.frame(input), as.data.frame(before))
  expect_false("LandWeb" %in% names(input))
})

test_that("landweb_sppEquiv() overwrites an existing LandWeb column", {
  input <- lr_sppEquiv()[, LandWeb := "stale"]
  expect_identical(
    as.data.frame(landweb_sppEquiv(input)),
    as.data.frame(landweb_sppEquiv(lr_sppEquiv()))
  )
})

test_that("landweb_sppEquiv() needs a LandR column to match species by", {
  expect_error(landweb_sppEquiv(lr_sppEquiv()[, !"LandR"]), "no `LandR` column")
})

test_that("landweb_sppEquiv() rejects tables it cannot map", {
  expect_snapshot(error = TRUE, {
    landweb_sppEquiv(as.data.frame(lr_sppEquiv()))
    landweb_sppEquiv(lr_sppEquiv()[, !"SCANFI"])
    landweb_sppEquiv(lr_sppEquiv()[!SCANFI %in% names(landweb_species_map())])
  })
})

test_that("generic Pinus contorta takes var. latifolia's row, not the first variety", {
  input <- lr_sppEquiv()
  ## shore pine comes first in LandR's table, so a first-match donor would be var. contorta
  expect_lt(which(input$SCANFI == "PINU_CON_CON"), which(input$SCANFI == "PINU_CON_LAT"))
  input[SCANFI == "PINU_CON_CON", colorHex := "#111111"]
  input[SCANFI == "PINU_CON_LAT", colorHex := "#222222"]
  input[Latin_full == "Pinus contorta", colorHex := ""]
  out <- landweb_sppEquiv(input)
  expect_identical(out[Latin_full == "Pinus contorta", colorHex], "#222222")
  expect_identical(out[Latin_full == "Pinus contorta", LandWeb], "Pinu_spp")
})

## landweb_dominant_sppEquiv / landweb_member_cover --------------------------------------------------

## cover by SCANFI code, shaped like WesternAlbertaUpland's (percent of each group's cover)
wauCover <- c(ABIE_LAS = 77.5, ABIE_BAL = 22.1, THUJ_PLI = 0.3, TSUG_HET = 0.1, LARI_LAR = 99.7,
              LARI_OCC = 0.3, PICE_GLA = 83, PICE_ENG = 14.7, PICE_ENG_GLA = 2.3, PICE_MAR = 100,
              PINU_CON_LAT = 92.5, PINU_BAN = 7.5, PINU_CON_CON = 0, POPU_TRE = 82.7, POPU_BAL = 11.7,
              BETU_PAP = 5.7, PSEU_MEN = 74, PSEU_MEN_GLA = 26)

test_that("landweb_dominant_sppEquiv() keeps only the dominant member's trait code in each group", {
  se <- landweb_sppEquiv(lr_sppEquiv())
  out <- landweb_dominant_sppEquiv(se, wauCover)
  kept <- out[!is.na(LANDIS_traits) & nzchar(LANDIS_traits), unique(LANDIS_traits), by = LandWeb]
  expect_identical(kept[LandWeb == "Pice_gla", V1], "PICE.GLA")
  expect_identical(kept[LandWeb == "Pinu_spp", V1], "PINU.CON.LAT")
  expect_identical(kept[LandWeb == "Abie_spp", V1], "ABIE.LAS")
  expect_identical(kept[LandWeb == "Popu_spp", V1], "POPU.TRE")
  expect_identical(kept[LandWeb == "Lari_spp", V1], "LARI.LAR")
  ## a single-code group is untouched
  expect_identical(out[LandWeb == "Pice_mar", LANDIS_traits], se[LandWeb == "Pice_mar", LANDIS_traits])
  ## only LANDIS_traits changes, and the input is not modified
  rest <- out[, !"LANDIS_traits"]
  data.table::setattr(rest, "dominant", NULL)
  expect_identical(rest, se[, !"LANDIS_traits"])
  expect_false(anyNA(se$LANDIS_traits[nzchar(se$SCANFI)]))
  dom <- attr(out, "dominant")
  expect_equal(dom[LandWeb == "Pice_gla", share], 0.83)
})

test_that("rows sharing the dominant trait code keep it (hybrid spruce shares Engelmann's)", {
  se <- landweb_sppEquiv(lr_sppEquiv())
  cover <- c(PICE_ENG = 40, PICE_ENG_GLA = 30, PICE_GLA = 50)
  out <- landweb_dominant_sppEquiv(se, cover)
  ## Engelmann + hybrid (PICE.ENG) = 70 > white spruce (PICE.GLA) = 50
  expect_identical(out[SCANFI %in% c("PICE_ENG", "PICE_ENG_GLA"), LANDIS_traits], c("PICE.ENG", "PICE.ENG"))
  expect_true(is.na(out[SCANFI == "PICE_GLA", LANDIS_traits]))
})

test_that("a group with no cover in the study area is left unchanged", {
  se <- landweb_sppEquiv(lr_sppEquiv())
  out <- landweb_dominant_sppEquiv(se, c(PICE_MAR = 10))
  expect_identical(out$LANDIS_traits, se$LANDIS_traits)
  expect_identical(nrow(attr(out, "dominant")), 0L)
})

test_that("landweb_member_cover() sums each member's cover layer, keyed by SCANFI code", {
  dir <- withr::local_tempdir()
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(10, 20, NA, 30))
  terra::writeRaster(r, file.path(dir, "SCANFI_spsCC_PICE_GLA_2020_v2_20260119_4px_SA.tif"))
  terra::writeRaster(r * 2, file.path(dir, "SCANFI_spsCC_PICE_ENG_GLA_2020_v2_20260119_4px_SA.tif"))
  terra::writeRaster(r, file.path(dir, "SCANFI_Pice_gla_4px.tif")) ## merged layer: not a member
  cov <- landweb_member_cover(dir)
  expect_identical(sort(names(cov)), c("PICE_ENG_GLA", "PICE_GLA"))
  expect_equal(cov[["PICE_GLA"]], 60)
  expect_equal(cov[["PICE_ENG_GLA"]], 120)
  expect_error(landweb_member_cover(withr::local_tempdir()), "No SCANFI_spsCC")
})

## landweb_require_species ---------------------------------------------------------------------------

test_that("landweb_require_species() returns non-empty inputs unchanged", {
  dt <- data.table::data.table(species = c("Pice_mar", "Popu_spp"))
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 2, vals = 1:8)
  lst <- list(speciesLayers = r, other = NULL)
  expect_identical(landweb_require_species(dt, "cohortData"), dt)
  expect_identical(landweb_require_species(r, "speciesLayers"), r)
  expect_identical(landweb_require_species(lst, "speciesLayers"), lst)
})

test_that("landweb_require_species() stops on an empty species input", {
  expect_snapshot(error = TRUE, {
    landweb_require_species(data.table::data.table(species = character(0)), "cohortData")
    landweb_require_species(NULL, "sppEquiv")
    landweb_require_species(list(speciesLayers = NULL), "speciesLayers")
    landweb_require_species(list(other = 1), "speciesLayers")
  })
})
