test_that("reportingPolygonLayers returns the expected schema", {
  lyrs <- reportingPolygonLayers()
  expect_named(lyrs, c("key", "NAME_SHORT", "isTenure", "cross", "source", "id", "labelCols", "where"))
  expect_setequal(unique(lyrs$source), c("drive", "url", "assembly"))
  expect_contains(lyrs$key, "FMA Boundaries Updated")

  ## exactly one tenure LAYER (rows sharing a NAME_SHORT are one layer), never crossed with itself
  byShort <- split(lyrs, lyrs$NAME_SHORT)
  expect_identical(sum(vapply(byShort, function(d) d$isTenure[[1L]], logical(1))), 1L)
  expect_false(any(lyrs$isTenure & lyrs$cross))

  ## caribou is now a single ASSEMBLY row -- the six jurisdictional sources live in
  ## caribouRangeLayers() and are merged by buildCaribouRanges(). See test-caribou.R.
  expect_identical(sum(lyrs$NAME_SHORT == "Caribou"), 1L)
  expect_identical(lyrs$source[lyrs$NAME_SHORT == "Caribou"], "assembly")
})

test_that("reportingPolygonLayers rejects an unknown source", {
  ## a typo would otherwise fall through to the zipped-URL path and fail later as an opaque
  ## "no vector file found"
  lyrs <- reportingPolygonLayers()
  lyrs$source[[1L]] <- "geojsn"
  expect_error(
    LandWebUtils:::.validateLayerSources(lyrs$source),
    "unknown source.*geojsn"
  )
  expect_silent(LandWebUtils:::.validateLayerSources(c("drive", "url", "geojson")))
})

test_that(".mergeLayerSources merges rows sharing a NAME_SHORT, and leaves singles alone", {
  mk <- function(nm, extra = NULL) {
    v <- terra::vect("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
    v$Name <- nm
    if (!is.null(extra)) v$extra <- extra
    terra::crs(v) <- "EPSG:3857"
    v
  }
  ## single source: untouched, keeps its other attributes
  one <- list(ANSR = mk("Montane", extra = "keep"))
  expect_identical(names(LandWebUtils:::.mergeLayerSources(one)), "ANSR")
  expect_contains(names(LandWebUtils:::.mergeLayerSources(one)[["ANSR"]]), "extra")

  ## two sources for one layer: merged, reduced to the common `Name`
  two <- list(Caribou = mk("Chinchaga", extra = "a"), Caribou = mk("Northwest Territories (NWT)"))
  m <- LandWebUtils:::.mergeLayerSources(two)
  expect_named(m, "Caribou")
  expect_setequal(as.character(m[["Caribou"]]$Name), c("Chinchaga", "Northwest Territories (NWT)"))
  expect_identical(names(m[["Caribou"]]), "Name")
})

test_that("landbaseLayers returns the expected schema", {
  lyrs <- landbaseLayers()
  expect_named(lyrs, c("key", "NAME_SHORT", "tenure", "source", "id", "layer", "status_col"))
  expect_contains(lyrs$key, "Spray Lake C5 Landbase")
  expect_false(any(is.na(lyrs$status_col)))

  ## every landbase belongs to a curated tenure -- a typo here silently disables that landbase
  expect_contains(tenureShortNames()$NAME_SHORT, lyrs$tenure)
})

test_that(".slug uses underscores, not the dots make.names() would produce", {
  expect_identical(LandWebUtils:::.slug("National Ecozones"), "National_Ecozones")
  expect_identical(LandWebUtils:::.slug("BC Biogeoclimatic zones"), "BC_Biogeoclimatic_zones")
  expect_identical(LandWebUtils:::.slug("Parks"), "Parks")
  expect_no_match(LandWebUtils:::.slug(reportingPolygonLayers()$key), "[.]", all = TRUE)
})

test_that(".slug collapses runs of non-alphanumerics and trims them", {
  expect_identical(LandWebUtils:::.slug("a  b--c"), "a_b_c")
  expect_identical(LandWebUtils:::.slug(" leading/trailing "), "leading_trailing")
})

test_that(".findVectorFile searches recursively (archives that unpack to a subdir)", {
  dir <- withr::local_tempdir()
  nested <- file.path(dir, "Ecoregions")
  dir.create(nested)
  file.create(file.path(nested, "ecoregions.shp"))

  ## a non-recursive list.files() found nothing here, so the layer was dropped silently
  expect_match(LandWebUtils:::.findVectorFile(dir), "ecoregions\\.shp$")
  expect_length(LandWebUtils:::.findVectorFile(withr::local_tempdir()), 0L)
})

test_that(".findVectorFile finds a GeoJSON (the ArcGIS REST sources are not archives)", {
  dir <- withr::local_tempdir()
  file.create(file.path(dir, "Caribou_Ranges_NWT.geojson"))
  expect_match(LandWebUtils:::.findVectorFile(dir), "\\.geojson$")
})

test_that(".findLandbaseSource prefers a File Geodatabase, then a shapefile", {
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "coverage.gdb"))
  expect_match(LandWebUtils:::.findLandbaseSource(dir), "coverage\\.gdb$")

  dir2 <- withr::local_tempdir()
  file.create(file.path(dir2, "coverage.shp"))
  expect_match(LandWebUtils:::.findLandbaseSource(dir2), "coverage\\.shp$")

  expect_null(LandWebUtils:::.findLandbaseSource(withr::local_tempdir()))
})

test_that("buildLandbasePolygons skips (no download) when no landbase applies", {
  sa <- terra::vect("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
  terra::crs(sa) <- "EPSG:3857"

  out <- buildLandbasePolygons(sa, tenures = "NoSuchTenure", destinationPath = withr::local_tempdir())
  expect_identical(out, list())

  ## regression: gating used to match a regex against the STUDY-AREA name, so once study areas
  ## became ecoregion groups ("WesternAlbertaUpland") nothing matched and no landbase was ever
  ## fetched. A study-area name must not accidentally pass the tenure gate.
  out2 <- buildLandbasePolygons(
    sa, tenures = "WesternAlbertaUpland", destinationPath = withr::local_tempdir()
  )
  expect_identical(out2, list())
})

test_that("buildLandbasePolygons gates to the landbases of the tenures present", {
  sa <- terra::vect("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
  terra::crs(sa) <- "EPSG:3857"

  ## intercept the fetch/extract so the test needs no network: record which layers pass the
  ## gate, then let each layer fall through (no source file found -> NULL) so the loop
  ## continues and every gated layer is observed.
  fetched <- character()
  testthat::local_mocked_bindings(
    drive_download_once = function(file, path, ...) {
      fetched <<- c(fetched, path)
      invisible(path)
    },
    archive_extract_once = function(...) invisible(NULL),
    .package = "workflowtools"
  )

  expect_identical(
    buildLandbasePolygons(sa, tenures = "SprayLake", destinationPath = withr::local_tempdir()),
    list()
  )
  expect_length(fetched, 1L)
  expect_match(fetched, "Spray_Lake_C5", fixed = TRUE)

  ## a tenure with several coverages pulls all of them (West Fraser's three CLS blocks)
  fetched <- character()
  buildLandbasePolygons(
    sa, tenures = "Tolko_Vand_WF_SL", destinationPath = withr::local_tempdir()
  )
  expect_length(fetched, 3L)

  ## and several tenures pull the union of theirs
  fetched <- character()
  buildLandbasePolygons(
    sa, tenures = c("SprayLake", "Manning"), destinationPath = withr::local_tempdir()
  )
  expect_length(fetched, 2L)
})

test_that("labelCols coalesces a source's name columns left to right", {
  ## Historically (LandWeb#118) the Caribou ranges arrived pre-combined and labelling them on
  ## RANGE_NAME alone left the 53 of 74 ranges named in the OTHER jurisdictions' columns unnamed --
  ## and unnamed features are dropped downstream, silently removing every tenure x Caribou unit
  ## outside Ontario/Manitoba. Caribou is now assembled per jurisdiction, so the coalescing is
  ## within a single source; these are the two real multi-column cases.
  lyrs <- caribouRangeLayers()

  ## SK: CONUNIT, falling back to RNGEUNIT -- SK1 carries the literal "Not Applicable" placeholder
  sk <- terra::vect(rep("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))", 2L))
  sk$CONUNIT <- c("SK2 West", "Not Applicable")
  sk$RNGEUNIT <- c("SK2", "SK1")
  out <- LandWebUtils:::.labelReportingLayer(sk, lyrs[lyrs$juris == "SK", ])
  expect_identical(out$Name, c("SK2 West", "SK1"))

  ## AB: SUBUNIT, falling back to LOCALRANGE
  ab <- terra::vect(rep("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))", 2L))
  ab$SUBUNIT <- c("A La Peche", NA)
  ab$LOCALRANGE <- c("West Central", "Cold Lake")
  out <- LandWebUtils:::.labelReportingLayer(ab, lyrs[lyrs$juris == "AB", ])
  expect_identical(out$Name, c("A La Peche", "Cold Lake"))
})

test_that("placeholder strings never become a reporting-unit name", {
  expect_identical(
    LandWebUtils:::.blankToNA(c("Real", "Not Applicable", "", "  ", "N/A", "none")),
    c("Real", NA, NA, NA, NA, NA)
  )
})

test_that(".labelColList splits the comma-separated priority list", {
  expect_identical(LandWebUtils:::.labelColList("A,B , C"), c("A", "B", "C"))
  expect_identical(LandWebUtils:::.labelColList("Name"), "Name")
  expect_length(LandWebUtils:::.labelColList(NA_character_), 0L)
})
