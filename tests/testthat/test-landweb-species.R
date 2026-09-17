## Fixtures -----------------------------------------------------------------------------------------

## `LandR::sppEquivalencies_CA` (LandR 1.2.0.9007): every row, with the columns the LandWeb groups
## read or relabel.
lr_sppEquiv <- function() {
  data.table::fread(
    test_path("fixtures", "sppEquivalencies_CA_LandR-1.2.0.9007.csv"),
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

## landweb_sppEquiv ---------------------------------------------------------------------------------

test_that("landweb_sppEquiv() reproduces the preamble's inline table exactly", {
  expect_identical(
    as.data.frame(landweb_sppEquiv(lr_sppEquiv())),
    as.data.frame(sppEquiv_inline_preamble_1.0.7(lr_sppEquiv()))
  )
})

test_that("landweb_sppEquiv() keeps one row per mapped species and drops the rest", {
  out <- landweb_sppEquiv(lr_sppEquiv())
  expect_setequal(out$SCANFI, names(landweb_species_map()))
  expect_identical(out$LandWeb, unname(landweb_species_map()[out$SCANFI]))
})

test_that("landweb_sppEquiv() gives each merged group a single label", {
  out <- landweb_sppEquiv(lr_sppEquiv())
  merged <- out[LandWeb %in% c("Lari_spp", "Pice_gla", "Pinu_spp", "Popu_spp", "Pseu_men")]
  nLabels <- merged[, .(
    short = data.table::uniqueN(EN_generic_short),
    full = data.table::uniqueN(EN_generic_full),
    leading = data.table::uniqueN(Leading)
  ), by = LandWeb]
  expect_identical(unique(c(nLabels$short, nLabels$full, nLabels$leading)), 1L)
  expect_identical(out[LandWeb == "Pinu_spp", unique(Leading)], "Pine leading")
})

test_that("landweb_sppEquiv() leaves Abie_spp with its species' own labels (pinned, not fixed)", {
  out <- landweb_sppEquiv(lr_sppEquiv())
  expect_setequal(
    out[LandWeb == "Abie_spp", EN_generic_full],
    c("Balsam fir", "Subalpine fir", "Western redcedar", "Western hemlock")
  )
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

test_that("landweb_sppEquiv() rejects tables it cannot map", {
  expect_snapshot(error = TRUE, {
    landweb_sppEquiv(as.data.frame(lr_sppEquiv()))
    landweb_sppEquiv(lr_sppEquiv()[, !"SCANFI"])
  })
})
