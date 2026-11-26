## a small species table: two spruce-group species, two pines, two broadleaves, white spruce
groups_equiv <- function() {
  data.table::data.table(
    LandWeb = c("Pice_mar", "Lari_lar", "Pinu_ban", "Pinu_con", "Popu_tre", "Betu_pap", "Pice_gla"),
    LandWebReport = c("Bl_Spruce", "Bl_Spruce", "Pine", "Pine", "Decid", "Decid", "Wh_Spruce"),
    Type = c("Conifer", "Conifer", "Conifer", "Conifer", "Deciduous", "Deciduous", "Conifer")
  )
}

test_that("each simulated code gets one group, and each group one Type", {
  grp <- sppEquiv_groups(groups_equiv(), "LandWeb", "LandWebReport")
  expect_identical(grp$col, "LandWebReport")
  expect_identical(
    unname(grp$map[c("Lari_lar", "Betu_pap", "Pinu_con")]),
    c("Bl_Spruce", "Decid", "Pine")
  )
  expect_identical(grp$groups[["LandWebReport"]], c("Bl_Spruce", "Decid", "Pine", "Wh_Spruce"))
  expect_identical(grp$groups[["Type"]], c("Conifer", "Deciduous", "Conifer", "Conifer"))

  ## varieties sharing a code collapse to one row
  twice <- rbind(
    groups_equiv(),
    data.table::data.table(LandWeb = "Pinu_con", LandWebReport = "Pine", Type = "Conifer")
  )
  expect_identical(sppEquiv_groups(twice, "LandWeb", "LandWebReport")$map, grp$map)
})

test_that("a code in two groups, a group of conifers and broadleaves, or a missing column stop", {
  twoGroups <- rbind(
    groups_equiv(),
    data.table::data.table(LandWeb = "Pinu_con", LandWebReport = "Fir", Type = "Conifer")
  )
  expect_error(
    sppEquiv_groups(twoGroups, "LandWeb", "LandWebReport"),
    "more than one group: Pinu_con"
  )

  mixed <- groups_equiv()
  data.table::set(mixed, i = 2L, j = "Type", value = "Deciduous") ## tamarack as a broadleaf
  expect_error(
    sppEquiv_groups(mixed, "LandWeb", "LandWebReport"),
    "both conifers and broadleaves: Bl_Spruce"
  )

  noGroup <- groups_equiv()
  data.table::set(noGroup, i = 7L, j = "LandWebReport", value = NA_character_)
  expect_error(sppEquiv_groups(noGroup, "LandWeb", "LandWebReport"), "No group .* Pice_gla")

  expect_error(sppEquiv_groups(groups_equiv(), "LandWeb", "Report"), "lacks column")
})

test_that("recode_cohorts() returns a recoded copy and stops on a species with no group", {
  grp <- sppEquiv_groups(groups_equiv(), "LandWeb", "LandWebReport")
  cd <- data.table::data.table(
    pixelGroup = c(1L, 1L, 2L),
    speciesCode = factor(c("Pice_mar", "Lari_lar", "Popu_tre")),
    age = 50L,
    B = c(400L, 400L, 200L)
  )
  before <- data.table::copy(cd)
  out <- recode_cohorts(cd, grp$map)

  expect_identical(cd, before)
  expect_identical(names(out), c("pixelGroup", "speciesCode", "B"))
  expect_identical(out[["speciesCode"]], c("Bl_Spruce", "Bl_Spruce", "Decid"))
  expect_identical(out[["B"]], cd[["B"]])

  cd2 <- data.table::data.table(pixelGroup = 1L, speciesCode = "Abie_bal", B = 100L)
  expect_error(recode_cohorts(cd2, grp$map), "No group for species: Abie_bal")
})

test_that("top_species keeps the n most abundant species, ties by name, none that are absent", {
  s <- c(Pice_mar = 40, Popu_tre = 25, Pinu_con = 25, Betu_pap = 0, Abie_bal = NA, Pice_gla = 10)
  expect_identical(top_species(s, 3), c("Pice_mar", "Pinu_con", "Popu_tre")) ## tie: Pinu_con before Popu_tre
  expect_identical(top_species(rev(s), 3), top_species(s, 3)) ## input order does not matter
  expect_identical(top_species(s, 10), c("Pice_mar", "Pinu_con", "Popu_tre", "Pice_gla")) ## 0 and NA never
  expect_identical(top_species(s, 0), character(0))
  expect_error(top_species(c(1, 2), 1), "named numeric")
  expect_error(top_species(s, NA), "non-negative")
  expect_error(top_species(s, -1), "non-negative")
  expect_error(top_species(s, c(1, 2)), "non-negative")
  expect_error(top_species(s, "3"), "non-negative")
  expect_error(top_species(s, 2.5), "whole number")
  expect_error(top_species(c(Pice_mar = 1, 5), 2), "unique names")
  expect_error(top_species(c(Pice_mar = 1, NA_real_)[1:2] |> stats::setNames(c("a", NA)), 2), "unique names")
  expect_error(top_species(c(a = 3, a = 2), 2), "unique names")
  expect_identical(top_species(c(a = -5, b = NaN, c = 2, d = -1), 4), "c") ## negative and NaN dropped
  expect_identical(top_species(c(abie = 1, Pice = 1), 2), c("Pice", "abie")) ## C-locale byte order on ties
})
