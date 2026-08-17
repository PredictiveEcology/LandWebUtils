test_that("caribouRangeLayers returns the expected schema and the six jurisdictions", {
  lyrs <- caribouRangeLayers()
  expect_named(lyrs, c("juris", "key", "source", "id", "labelCols", "statusCols", "extirpated"))
  expect_setequal(lyrs$juris, c("AB", "BC", "SK", "MB", "NWT", "ON"))

  ## Ontario is easy to drop from the list and its loss is silent: `RANGE_NAME` was populated for
  ## exactly the ON+MB features of the old combined layer, so omitting ON costs every ON tenure x
  ## Caribou reporting unit.
  expect_contains(lyrs$juris, "ON")

  expect_true(all(lyrs$source %in% LandWebUtils:::.LAYER_SOURCES))
  expect_false(anyDuplicated(lyrs$juris) > 0L)
  expect_false(anyDuplicated(lyrs$key) > 0L)
})

test_that("caribouRangeLayers labels AB on SUBUNIT before LOCALRANGE (the v2 granularity)", {
  ab <- caribouRangeLayers()
  ab <- ab[ab$juris == "AB", ]
  ## A La Peche / Narraway / Redrock-Prairie Creek are three separate v2 units but share the single
  ## LOCALRANGE "West Central", so labelling on LOCALRANGE alone would collapse them into one.
  expect_identical(LandWebUtils:::.labelColList(ab$labelCols)[[1L]], "SUBUNIT")
})

test_that("caribouRangeLayers labels NWT on the ASCII-safe column", {
  nwt <- caribouRangeLayers()
  nwt <- nwt[nwt$juris == "NWT", ]
  ## `NAME` carries diacritics (Sahtu', Wek'eezhii) and this label is slugged into the refCode that
  ## keys parquet aggregates and figure filenames.
  expect_identical(nwt$labelCols, "REGION")
  expect_match(nwt$id, "/97/query\\?")
})

test_that("only AB and BC carry an extirpation rule; the others ship no status field", {
  lyrs <- caribouRangeLayers()
  withRule <- lyrs$juris[!is.na(lyrs$statusCols) & !is.na(lyrs$extirpated)]
  expect_setequal(withRule, c("AB", "BC"))

  ## BC's field name differs between its WFS (HERD_STATUS) and a 10-char-truncated shapefile export
  ## (HERD_STAT), so both must be candidates.
  bc <- lyrs[lyrs$juris == "BC", ]
  expect_setequal(LandWebUtils:::.labelColList(bc$statusCols), c("HERD_STATUS", "HERD_STAT"))
})

test_that(".dropExtirpated drops matching rows, tolerates a missing/renamed column, and no-ops on NA", {
  ## one row PER status value -- a fixed row count would recycle a shorter vector and could pass for
  ## the wrong reason (it did: the 2-value AB case was silently recycled to 3 rows)
  mk <- function(status, col = "HERD_STATUS") {
    v <- terra::vect(rep("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))", length(status)))
    v[[col]] <- status
    v
  }
  v <- mk(c("Herd", "Extirpated", "Trace occurrences"))
  expect_identical(as.integer(nrow(LandWebUtils:::.dropExtirpated(v, "HERD_STATUS,HERD_STAT", "^Extirpated$"))), 2L)

  ## the truncated shapefile spelling is found via the candidate list
  v2 <- mk(c("Herd", "Extirpated", "Herd"), col = "HERD_STAT")
  expect_identical(as.integer(nrow(LandWebUtils:::.dropExtirpated(v2, "HERD_STATUS,HERD_STAT", "^Extirpated$"))), 2L)

  ## AB's vocabulary is a prefix ("Extirp"), not the whole word
  vab <- mk(c("Active", "Extirp"), col = "STATUS")
  expect_identical(as.integer(nrow(LandWebUtils:::.dropExtirpated(vab, "STATUS", "^Extirp"))), 1L)

  ## a source with no status field, or no rule, is passed through untouched
  expect_identical(as.integer(nrow(LandWebUtils:::.dropExtirpated(v, NA_character_, NA_character_))), 3L)
  expect_identical(as.integer(nrow(LandWebUtils:::.dropExtirpated(v, "NO_SUCH_COL", "^Extirpated$"))), 3L)
})

test_that(".caribouEcotype splits boreal from mountain per source", {
  mk <- function(nm, vals) {
    v <- terra::vect("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
    v <- do.call(rbind, replicate(length(vals), v, simplify = FALSE))
    v[[nm]] <- vals
    v
  }
  ## BC: "Northern" AND "Mountain" are both mountain caribou designatable units
  bc <- mk("ECOTYPE", c("Boreal", "Northern", "Mountain"))
  expect_identical(LandWebUtils:::.caribouEcotype(bc, "BC"), c("Boreal", "Mountain", "Mountain"))

  ## AB: Banff / Jasper / West Central are the mountain ranges
  ab <- mk("LOCALRANGE", c("Chinchaga", "Jasper", "West Central", "Banff", "Cold Lake"))
  expect_identical(
    LandWebUtils:::.caribouEcotype(ab, "AB"),
    c("Boreal", "Mountain", "Mountain", "Mountain", "Boreal")
  )

  ## everything else is boreal-only within this AOI
  sk <- mk("RNGEUNIT", c("SK1", "SK2"))
  expect_identical(LandWebUtils:::.caribouEcotype(sk, "SK"), c("Boreal", "Boreal"))
})

test_that(".caribouNameFixes normalises the cross-jurisdiction spelling disagreements", {
  fx <- LandWebUtils:::.caribouNameFixes()
  ## Left unfixed these reach partner-facing figures AND block the name-grouping that stitches a
  ## cross-border range into one reporting unit.
  expect_identical(fx[["Bischto"]], "Bistcho")            ## Alberta's own typo
  expect_identical(fx[["Snake-Sahtaneh"]], "Snake-Sahtahneh") ## BC drops the second 'h'
  expect_false(any(names(fx) == unname(fx))) ## a no-op entry would be a mistake
})

test_that("reportingPolygonLayers routes caribou through the assembler", {
  lyrs <- reportingPolygonLayers()
  cb <- lyrs[lyrs$NAME_SHORT == "Caribou", ]

  ## ONE row now: the six jurisdictional sources are assembled by buildCaribouRanges(), replacing
  ## both Julie's hand-combined v10 layer (which had dropped the NWT ranges) and the interim
  ## GNWT-only supplement that patched the hole.
  expect_identical(nrow(cb), 1L)
  expect_identical(cb$source, "assembly")
  expect_identical(cb$id, "caribou")
  expect_true(cb$cross)
  expect_false(cb$isTenure)

  ## the assembler must exist, or buildReportingPolygons() would fail at run time
  expect_contains(names(LandWebUtils:::.ASSEMBLERS), cb$id)
  expect_true(is.function(LandWebUtils:::.ASSEMBLERS[[cb$id]]))
})

test_that("the ECCC national layer is NOT a reporting layer", {
  ## Reporting goes to jurisdictional partners, who each want their own management-unit names. ECCC
  ## is kept only as a documented comparison in scripts/make_caribou_reference.R.
  lyrs <- reportingPolygonLayers()
  expect_false(any(grepl("ECCC|Priority Species|All51", lyrs$key)))
  expect_false(any(grepl("data-donnees", lyrs$id[!is.na(lyrs$id)])))
})
