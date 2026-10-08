## `LandR::sppEquivalencies_CA` (LandR 1.2.0.9047): every row, as in test-landweb-species.R, plus
## the `Type` column that `LandR::vegTypeMapGenerator()` reads to find mixedwood stands.
lr_table <- function() {
  data.table::fread(
    test_path("fixtures", "sppEquivalencies_CA_LandR-1.2.0.9047.csv"),
    colClasses = "character",
    na.strings = NULL
  )
}

test_that("each reporting group has one code, label, type and colour", {
  g <- landweb_report_groups()
  expect_identical(nrow(g), 6L)
  for (col in c("code", "label", "colour")) {
    expect_identical(anyDuplicated(g[[col]]), 0L, label = col)
  }
  expect_true(all(g[["Type"]] %in% c("Conifer", "Deciduous")))
  expect_true(all(grepl("^#[0-9A-F]{6}$", g[["colour"]])))
})

test_that("every code LandWeb simulates has a reporting group", {
  map <- landweb_report_map()
  expect_true(all(map %in% landweb_report_groups()[["code"]]))

  ## today's merged groups, and every species in them by its LandR code
  tab <- lr_table()
  spp <- unique(tab[["LandR"]][tab[["SCANFI"]] %in% names(landweb_species_map())])
  codes <- c(unique(unname(landweb_species_map())), spp)
  expect_identical(setdiff(codes, names(map)), character(0))
})

test_that("species sit in the data dictionary's groups, with Douglas-fir on its own", {
  map <- landweb_report_map()
  expect_identical(unname(map[c("Pice_gla", "Pice_eng", "Pice_eng_gla")]), rep("Wh_Spruce", 3L))
  expect_identical(unname(map[c("Pice_mar", "Lari_lar", "Lari_occ")]), rep("Bl_Spruce", 3L))
  expect_identical(unname(map[c("Pinu_ban", "Pinu_con")]), rep("Pine", 2L))
  expect_identical(unname(map[c("Abie_bal", "Abie_las", "Thuj_pli", "Tsug_het")]), rep("Fir", 4L))
  expect_identical(unname(map[c("Popu_tre", "Popu_bal", "Betu_pap")]), rep("Decid", 3L))
  expect_identical(unname(map["Pseu_men"]), "Doug_fir")
})

test_that("a reporting group's members all have the group's Type", {
  ## LandR::vegTypeMapGenerator() merges one Type per code; a group spanning conifers and
  ## broadleaves would get two rows there and be counted twice.
  tab <- lr_table()
  map <- landweb_report_map()
  groups <- landweb_report_groups()
  spp <- unique(tab[LandR %in% names(map) & nzchar(Type), list(LandR, Type)])
  expect_identical(anyDuplicated(spp[["LandR"]]), 0L) ## one Type per species code
  groupType <- groups[["Type"]][match(map[spp[["LandR"]]], groups[["code"]])]
  expect_identical(spp[["LandR"]][spp[["Type"]] != groupType], character(0))
  ## and the table covers every species code of the map that LandR types
  expect_true(all(c("Lari_lar", "Betu_pap", "Pseu_men", "Thuj_pli") %in% spp[["LandR"]]))
})

test_that("landweb_add_report_column() adds each row's group and leaves its input alone", {
  eq <- landweb_sppEquiv(lr_table())
  before <- data.table::copy(eq)
  out <- landweb_add_report_column(eq)

  expect_identical(eq, before)
  expect_identical(out[["LandWebReport"]], unname(landweb_report_map()[eq[["LandWeb"]]]))
  expect_setequal(unique(out[["LandWebReport"]]), landweb_report_groups()[["code"]])
  ## one group per simulated code
  perCode <- out[, list(n = data.table::uniqueN(LandWebReport)), by = "LandWeb"]
  expect_true(all(perCode[["n"]] == 1L))
  ## the reporting groups regroup today's simulated groups as the data dictionary does
  byGroup <- unique(out[, list(LandWeb, LandWebReport)])
  expect_identical(
    byGroup[order(LandWeb), LandWebReport],
    c("Fir", "Bl_Spruce", "Wh_Spruce", "Bl_Spruce", "Pine", "Decid", "Doug_fir")
  )
})

test_that("landweb_add_report_column() errors on a code with no group and on an existing column", {
  expect_error(
    landweb_add_report_column(data.table::data.table(LandWeb = c("Pice_mar", "Acer_rub"))),
    "Acer_rub"
  )
  expect_error(
    landweb_add_report_column(data.table::data.table(LandWeb = "Pice_mar", LandWebReport = "x")),
    "already has"
  )
  expect_error(landweb_add_report_column(data.frame(LandWeb = "Pice_mar")), "data.table")
})

test_that("landweb_report_colours() names the requested groups plus Mixed, in table order", {
  expect_identical(names(landweb_report_colours(c("Decid", "Pine"))), c("Pine", "Decid", "Mixed"))
  all <- landweb_report_colours()
  expect_identical(names(all), c(landweb_report_groups()[["code"]], "Mixed"))
  expect_identical(anyDuplicated(unname(all)), 0L)
  expect_error(landweb_report_colours("Larch"), "Larch")
})
