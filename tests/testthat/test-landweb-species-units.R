## `LandR::sppEquivalencies_CA` (LandR 1.2.0.9047), as in test-landweb-species.R.
lr_table_units <- function() {
  data.table::fread(
    test_path("fixtures", "sppEquivalencies_CA_LandR-1.2.0.9047.csv"),
    colClasses = "character",
    na.strings = NULL
  )
}

## A 10 x 11 landscape: cells 1-100 hold the compositions below (percent cover, each summing to 100),
## cells 101-110 are NA in every layer (outside the study area).
##   1-50   white spruce 60, lodgepole pine 40
##   51-70  aspen 70, black spruce 30
##   71-80  jack pine 50, lodgepole pine 50
##   81-90  hybrid spruce 50, white spruce 50
##   91     Engelmann spruce 40, white spruce 60
##   92-93  balsam fir 30, subalpine fir 30, white spruce 40
##   94     western redcedar 20, balsam fir 20, white spruce 60
##   95     tamarack 4, black spruce 96         (tamarack is under the 5% floor everywhere)
##   96     Douglas-fir 50, white spruce 50
##   97-100 balsam poplar 50, aspen 50
units_layers <- function() {
  spp <- c(
    "Pice_gla",
    "Pinu_con",
    "Popu_tre",
    "Pice_mar",
    "Pinu_ban",
    "Pice_eng_gla",
    "Pice_eng",
    "Abie_bal",
    "Abie_las",
    "Thuj_pli",
    "Lari_lar",
    "Pseu_men",
    "Popu_bal"
  )
  v <- matrix(0, nrow = 110L, ncol = length(spp), dimnames = list(NULL, spp))
  put <- function(cells, ...) {
    x <- list(...)
    for (s in names(x)) {
      v[cells, s] <<- x[[s]]
    }
  }
  put(1:50, Pice_gla = 60, Pinu_con = 40)
  put(51:70, Popu_tre = 70, Pice_mar = 30)
  put(71:80, Pinu_ban = 50, Pinu_con = 50)
  put(81:90, Pice_eng_gla = 50, Pice_gla = 50)
  put(91, Pice_eng = 40, Pice_gla = 60)
  put(92:93, Abie_bal = 30, Abie_las = 30, Pice_gla = 40)
  put(94, Thuj_pli = 20, Abie_bal = 20, Pice_gla = 60)
  put(95, Lari_lar = 4, Pice_mar = 96)
  put(96, Pseu_men = 50, Pice_gla = 50)
  put(97:100, Popu_bal = 50, Popu_tre = 50)
  v[101:110, ] <- NA
  r <- terra::rast(
    nrows = 10,
    ncols = 11,
    xmin = 0,
    xmax = 11,
    ymin = 0,
    ymax = 10,
    crs = "EPSG:3857",
    nlyrs = length(spp)
  )
  terra::values(r) <- v
  names(r) <- spp
  r
}

test_that("landweb_species_sppEquiv() keeps landweb_sppEquiv()'s rows, one code per species", {
  tab <- lr_table_units()
  before <- data.table::copy(tab)
  sp <- landweb_species_sppEquiv(tab)
  grp <- landweb_sppEquiv(tab)

  expect_identical(tab, before)
  expect_identical(sp[["SCANFI"]], grp[["SCANFI"]])
  expect_identical(sp[["LandR"]], grp[["LandR"]])
  expect_identical(sp[["LandWeb"]], sp[["LandR"]])
  expect_identical(sp[["LandWebGroup"]], grp[["LandWeb"]])
  expect_identical(sp[["LandWebReport"]], unname(landweb_report_map()[sp[["LandR"]]]))
})

test_that("every Pinus contorta row speaks for interior lodgepole pine", {
  sp <- landweb_species_sppEquiv(lr_table_units())
  pc <- sp[LandWeb == "Pinu_con"]
  expect_gte(nrow(pc), 3L) ## shore pine, var. latifolia, generic
  expect_identical(unique(pc[["EN_generic_short"]]), "Lg pine")
  expect_identical(unique(pc[["Leading"]]), "Lodgepole pine leading")
  expect_identical(unique(pc[["LANDIS_traits"]]), "PINU.CON.LAT")
})

test_that("every row of a species carries the same labels and colour", {
  sp <- landweb_species_sppEquiv(lr_table_units())
  perSpecies <- unique(sp[, list(LandWeb, EN_generic_short, EN_generic_full, Leading, colorHex)])
  expect_identical(anyDuplicated(perSpecies[["LandWeb"]]), 0L)
  expect_identical(perSpecies[LandWeb == "Popu_bal", EN_generic_short], "Ba poplar") ## not "Bl ctnwood"
  expect_identical(perSpecies[LandWeb == "Pseu_men", EN_generic_short], "Doug-fir")
})

test_that("every species has its own colour", {
  sp <- landweb_species_sppEquiv(lr_table_units())
  perSpecies <- unique(sp[, list(LandWeb, colorHex)])
  expect_identical(anyDuplicated(perSpecies[["LandWeb"]]), 0L) ## one colour per species
  expect_identical(anyDuplicated(perSpecies[["colorHex"]]), 0L) ## and no two species share one
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", perSpecies[["colorHex"]])))
})

test_that("species clearing the threshold run alone; the rest join, group or drop", {
  res <- landweb_species_units(
    landweb_species_sppEquiv(lr_table_units()),
    units_layers(),
    threshold = 1
  )

  expect_identical(res$hybridParent, "Pice_gla")
  expect_identical(
    res$layerMap,
    c(
      Pice_gla = "Pice_gla",
      Pinu_con = "Pinu_con",
      Popu_tre = "Popu_tre",
      Pice_mar = "Pice_mar",
      Pinu_ban = "Pinu_ban",
      Pice_eng_gla = "Pice_gla",
      Pice_eng = "Pice_gla",
      Abie_bal = "Abie_spp",
      Abie_las = "Abie_spp",
      Thuj_pli = "Abie_spp",
      Lari_lar = NA,
      Pseu_men = NA,
      Popu_bal = "Popu_bal"
    )
  )
  shares <- stats::setNames(res$units[["share"]], res$units[["unit"]])
  expect_equal(
    shares[order(names(shares))],
    c(
      Abie_spp = 1.6,
      Pice_gla = 43.4,
      Pice_mar = 7,
      Pinu_ban = 5,
      Pinu_con = 25,
      Popu_bal = 2,
      Popu_tre = 16
    )
  )
  u <- res$units[unit == "Abie_spp"]
  expect_identical(u[["kind"]], "group")
  expect_identical(u[["dominant"]], "Abie_bal") ## cedar never leads a unit
})

test_that("a unit takes its dominant member's traits, labels and colour", {
  res <- landweb_species_units(
    landweb_species_sppEquiv(lr_table_units()),
    units_layers(),
    threshold = 1
  )
  eq <- res$sppEquiv

  expect_setequal(
    unique(eq[["LandWeb"]]),
    c("Abie_spp", "Pice_gla", "Pice_mar", "Pinu_ban", "Pinu_con", "Popu_bal", "Popu_tre")
  )
  traits <- unique(eq[!is.na(LANDIS_traits) & nzchar(LANDIS_traits), list(LandWeb, LANDIS_traits)])
  expect_identical(anyDuplicated(traits[["LandWeb"]]), 0L) ## one trait code per unit
  expect_identical(traits[LandWeb == "Abie_spp", LANDIS_traits], "ABIE.BAL")
  expect_identical(traits[LandWeb == "Pice_gla", LANDIS_traits], "PICE.GLA") ## hybrid rows too
  expect_identical(traits[LandWeb == "Pinu_con", LANDIS_traits], "PINU.CON.LAT")

  ## every row of a unit carries the unit's labels and colour
  perUnit <- unique(eq[, list(LandWeb, EN_generic_short, Leading, colorHex)])
  expect_identical(anyDuplicated(perUnit[["LandWeb"]]), 0L)
  expect_identical(perUnit[LandWeb == "Pice_gla", EN_generic_short], "Whi spr")
  expect_identical(perUnit[LandWeb == "Abie_spp", EN_generic_short], "Fir")
  ## one reporting group per unit
  expect_true(all(eq[, data.table::uniqueN(LandWebReport), by = "LandWeb"][["V1"]] == 1L))

  expect_identical(names(res$sppColorVect), c(sort(unique(eq[["LandWeb"]])), "Mixed"))
  expect_identical(anyDuplicated(unname(res$sppColorVect)), 0L)
  expect_identical(
    names(res$sppColorVectReport),
    c("Wh_Spruce", "Bl_Spruce", "Pine", "Fir", "Decid", "Mixed")
  )
})

test_that("a lower threshold splits more species; cedar joins the largest fir", {
  res <- landweb_species_units(
    landweb_species_sppEquiv(lr_table_units()),
    units_layers(),
    threshold = 0.5
  )
  expect_identical(
    unname(res$layerMap[c("Abie_bal", "Abie_las", "Thuj_pli", "Pseu_men")]),
    c("Abie_bal", "Abie_las", "Abie_bal", "Pseu_men")
  )
  expect_identical(res$units[unit == "Abie_las", kind], "species")
  eq <- res$sppEquiv
  expect_identical(unique(eq[LandWeb == "Abie_las", EN_generic_short]), "Subalp fir")
})

test_that("hybrid spruce joins Engelmann spruce where Engelmann spruce has more cover", {
  r <- units_layers()
  v <- terra::values(r)
  v[1:50, "Pice_eng"] <- 60 ## Engelmann now outweighs white spruce
  v[1:50, "Pice_gla"] <- 0
  terra::values(r) <- v
  res <- landweb_species_units(landweb_species_sppEquiv(lr_table_units()), r, threshold = 1)
  expect_identical(res$hybridParent, "Pice_eng")
  expect_identical(unname(res$layerMap["Pice_eng_gla"]), "Pice_eng")
  expect_identical(
    unique(res$sppEquiv[LandWeb == "Pice_eng" & SCANFI == "PICE_ENG_GLA", LANDIS_traits]),
    "PICE.ENG"
  )
})

test_that("shares are measured over the study area only", {
  sa <- terra::vect(terra::ext(0, 11, 5, 10), crs = "EPSG:3857") ## top half: cells 1-55
  res <- landweb_species_units(
    landweb_species_sppEquiv(lr_table_units()),
    units_layers(),
    studyArea = sa
  )
  expect_setequal(
    stats::na.omit(unique(res$layerMap)),
    c("Pice_gla", "Pinu_con", "Popu_tre", "Pice_mar")
  )
})

test_that("landweb_sum_layers() sums each unit's layers and keeps NA cells", {
  layers <- units_layers()
  res <- landweb_species_units(landweb_species_sppEquiv(lr_table_units()), layers, threshold = 1)
  summed <- landweb_sum_layers(layers, res$layerMap)

  expect_setequal(names(summed), unique(stats::na.omit(res$layerMap)))
  v <- terra::values(summed)
  expect_identical(unname(v[c(1, 85, 91, 92), "Pice_gla"]), c(60, 100, 100, 40))
  expect_identical(unname(v[c(92, 94), "Abie_spp"]), c(60, 40))
  expect_true(all(is.na(v[101:110, ])))
  expect_error(landweb_sum_layers(layers, res$layerMap[-1]), "No unit given")
  expect_error(
    landweb_sum_layers(
      layers,
      stats::setNames(rep(NA_character_, terra::nlyr(layers)), names(layers))
    ),
    "drops every layer"
  )
})
