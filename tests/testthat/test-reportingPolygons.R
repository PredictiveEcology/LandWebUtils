test_that("reportingPolygonLayers returns the expected schema", {
  lyrs <- reportingPolygonLayers()
  expect_named(
    lyrs,
    c("key", "NAME_SHORT", "isTenure", "cross", "source", "id", "labelCols", "where")
  )
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
  expect_error(LandWebUtils:::.validateLayerSources(lyrs$source), "unknown source.*geojsn")
  expect_silent(LandWebUtils:::.validateLayerSources(c("drive", "url", "geojson")))
})

test_that(".mergeLayerSources merges rows sharing a NAME_SHORT, and leaves singles alone", {
  mk <- function(nm, extra = NULL) {
    v <- terra::vect("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
    v$Name <- nm
    if (!is.null(extra)) {
      v$extra <- extra
    }
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

test_that("each landbase coverage belongs to the tenure it lies in", {
  ## A coverage is crossed only with its own tenure, so one assigned to a tenure it does not
  ## overlap produces no reporting unit, without an error. S20, S21 and C5 had such tenures: S20
  ## lies in MillarWestern, S21 in WF_Tolko and C5 in Crowsnest.
  expected <- tibble::tribble(
    ~key                              , ~NAME_SHORT , ~tenure            ,
    "Spray Lake C5 Landbase"          , "LB_C5"     , "Crowsnest"        ,
    "West Fraser Blue Ridge Landbase" , "LB"        , "BlueRidge"        ,
    "West Fraser N CLS S17 Landbase"  , "LB_S17"    , "Tolko_Vand_WF_SL" ,
    "West Fraser N CLS S20 Landbase"  , "LB_S20"    , "MillarWestern"    ,
    "West Fraser N CLS S21 Landbase"  , "LB_S21"    , "WF_Tolko"         ,
    "Tolko AB North Landbase"         , "LB"        , "Tolko_Norbord_LC" ,
    "Sundre FP Landbase"              , "LB"        , "Sundre"           ,
    "Manning Landbase"                , "LB"        , "Manning"
  )
  expect_identical(landbaseLayers()[, c("key", "NAME_SHORT", "tenure")], expected)
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

  out <- buildLandbasePolygons(
    sa,
    tenures = "NoSuchTenure",
    destinationPath = withr::local_tempdir()
  )
  expect_identical(out, list())

  ## regression: gating used to match a regex against the STUDY-AREA name, so once study areas
  ## became ecoregion groups ("WesternAlbertaUpland") nothing matched and no landbase was ever
  ## fetched. A study-area name must not accidentally pass the tenure gate.
  out2 <- buildLandbasePolygons(
    sa,
    tenures = "WesternAlbertaUpland",
    destinationPath = withr::local_tempdir()
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
  gdbs <- character()
  testthat::local_mocked_bindings(.driveGdbDownloadOnce = function(id, destinationPath, ...) {
    gdbs <<- c(gdbs, id)
    NULL
  })

  expect_identical(
    buildLandbasePolygons(sa, tenures = "Crowsnest", destinationPath = withr::local_tempdir()),
    list()
  )
  expect_length(fetched, 1L)
  expect_match(fetched, "Spray_Lake_C5", fixed = TRUE)

  ## the company in a coverage's name does not decide its tenure: C5 lies in Crowsnest, and the
  ## Slave Lake consortium keeps only S17 of the three CLS coverages
  fetched <- character()
  buildLandbasePolygons(sa, tenures = "SprayLake", destinationPath = withr::local_tempdir())
  expect_length(fetched, 0L)
  expect_length(gdbs, 0L)
  buildLandbasePolygons(sa, tenures = "Tolko_Vand_WF_SL", destinationPath = withr::local_tempdir())
  expect_length(fetched, 0L)
  expect_identical(gdbs, "1oX9crrPy_PaZRwb09j9LsHB37PY3hY7s")

  ## and several tenures pull the union of theirs: S21's shapefile and S20's geodatabase
  fetched <- character()
  gdbs <- character()
  buildLandbasePolygons(
    sa,
    tenures = c("MillarWestern", "WF_Tolko"),
    destinationPath = withr::local_tempdir()
  )
  expect_length(fetched, 1L)
  expect_match(fetched, "CLS_S21", fixed = TRUE)
  expect_identical(gdbs, "1oX9crrPy_PaZRwb09j9LsHB37PY3hY7s")
})

test_that("CLS S17 and S20 are read from the cleaned geodatabase", {
  ## Their shapefiles' .dbf files are cut at the 2 GiB dBASE limit, so GDAL read only the shapes
  ## before the cut: S17 lost 1,092 km2 and S20 407 km2.
  lyrs <- landbaseLayers()
  cls <- lyrs[lyrs$NAME_SHORT %in% c("LB_S17", "LB_S20"), ]
  expect_identical(cls$source, c("drive_gdb", "drive_gdb"))
  expect_identical(cls$id, rep("1oX9crrPy_PaZRwb09j9LsHB37PY3hY7s", 2L))
  expect_identical(cls$layer, c("CLS_Clipped_S17_Clean", "CLS_Clipped_S20_Clean"))
  expect_identical(cls$status_col, c("STATUS", "STATUS"))
  expect_identical(cls$tenure, c("Tolko_Vand_WF_SL", "MillarWestern"))

  ## S21's shapefile is whole, and stays
  expect_identical(lyrs$source[lyrs$NAME_SHORT == "LB_S21"], "drive")
})

test_that("buildLandbasePolygons fetches a geodatabase once for all of its layers", {
  sa <- terra::vect("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))", crs = "EPSG:3857")
  dest <- withr::local_tempdir()
  gdb <- file.path(dest, "Landbases.gdb")

  nFetch <- 0L
  read <- list()
  testthat::local_mocked_bindings(.driveGdbDownloadOnce = function(id, destinationPath, ...) {
    nFetch <<- nFetch + 1L
    gdb
  })
  testthat::local_mocked_bindings(
    prep_landbase = function(src, aoi, status_col, layer = NULL, ...) {
      read[[length(read) + 1L]] <<- list(src = src, layer = layer, status_col = status_col)
      v <- terra::vect(c("POLYGON ((0 0, 1 0, 1 1, 0 0))", "POLYGON ((0 0, 1 1, 0 1, 0 0))"))
      v$lbstatus <- c("ACTIVE", "PASSIVE")
      v
    },
    .package = "spatialutils"
  )

  out <- buildLandbasePolygons(
    sa,
    tenures = c("Tolko_Vand_WF_SL", "MillarWestern"),
    destinationPath = dest
  )

  expect_identical(nFetch, 1L)
  expect_identical(vapply(read, `[[`, "", "src"), c(gdb, gdb))
  expect_identical(
    vapply(read, `[[`, "", "layer"),
    c("CLS_Clipped_S17_Clean", "CLS_Clipped_S20_Clean")
  )
  expect_identical(vapply(read, `[[`, "", "status_col"), c("STATUS", "STATUS"))
  expect_named(out, c("Tolko_Vand_WF_SL LB_S17", "MillarWestern LB_S20"))
  ## STATUS holds the same ACTIVE/PASSIVE values as the shapefiles' F_CONDITIO
  expect_identical(out[[1L]]$Name, c("ACTIVE", "PASSIVE"))
})

test_that(".driveGdbDownloadOnce downloads every file but lock files into <dest>/<folder name>", {
  dest <- withr::local_tempdir()
  calls <- local_fake_drive_gdb(list(
    gdb = "g",
    a00000001.gdbtable = "table one",
    a00000001.gdbtablx = "index",
    `_gdb.host.123.4.sr.lock` = "lock"
  ))

  gdb <- LandWebUtils:::.driveGdbDownloadOnce("folder", dest)

  expect_identical(gdb, file.path(dest, "Landbases.gdb"))
  ## the subfolder and the ArcGIS lock file are not downloaded
  expect_setequal(list.files(gdb), c("gdb", "a00000001.gdbtable", "a00000001.gdbtablx"))
  expect_identical(readLines(file.path(gdb, "a00000001.gdbtable"), warn = FALSE), "table one")
  ## and no temporary download is left behind
  expect_length(list.files(dest, pattern = "partial", all.files = TRUE, recursive = TRUE), 0L)
})

test_that(".driveGdbDownloadOnce skips a file present with Drive's size, and replaces one without", {
  dest <- withr::local_tempdir()
  gdb <- file.path(dest, "Landbases.gdb")
  dir.create(gdb)
  writeBin(charToRaw("g"), file.path(gdb, "gdb"))
  writeBin(charToRaw("cut"), file.path(gdb, "a00000001.gdbtable"))

  calls <- local_fake_drive_gdb(list(gdb = "g", a00000001.gdbtable = "table one"))
  LandWebUtils:::.driveGdbDownloadOnce("folder", dest)
  expect_identical(calls$downloaded, "a00000001.gdbtable")
  expect_identical(readLines(file.path(gdb, "a00000001.gdbtable"), warn = FALSE), "table one")

  ## a second call (another layer of the same gdb, or another study area) downloads nothing
  calls$downloaded <- character()
  LandWebUtils:::.driveGdbDownloadOnce("folder", dest)
  expect_length(calls$downloaded, 0L)
})

test_that(".driveGdbDownloadOnce keeps a file another worker placed during its download", {
  ## Study areas on different workers share the inputs directory, so two can fetch the gdb at
  ## once. The first worker may already be reading a file it placed; the second must not replace
  ## it with its own copy.
  dest <- withr::local_tempdir()
  gdb <- file.path(dest, "Landbases.gdb")
  local_fake_drive_gdb(list(a00000001.gdbtable = "table one"), on_download = function(name, path) {
    writeBin(charToRaw("TABLE ONE"), file.path(gdb, name))
  })

  LandWebUtils:::.driveGdbDownloadOnce("folder", dest)

  expect_identical(readLines(file.path(gdb, "a00000001.gdbtable"), warn = FALSE), "TABLE ONE")
  expect_length(list.files(dest, pattern = "partial", all.files = TRUE, recursive = TRUE), 0L)
})

test_that(".driveGdbDownloadOnce stops on a download of the wrong size, placing nothing", {
  dest <- withr::local_tempdir()
  local_fake_drive_gdb(list(a00000001.gdbtable = "table one"), on_download = function(name, path) {
    writeBin(charToRaw("table"), path)
  })

  expect_error(LandWebUtils:::.driveGdbDownloadOnce("folder", dest), "5 bytes.*Drive lists 9")
  expect_false(file.exists(file.path(dest, "Landbases.gdb", "a00000001.gdbtable")))
  expect_length(list.files(dest, pattern = "partial", all.files = TRUE, recursive = TRUE), 0L)
})

test_that(".driveGdbDownloadOnce stops on a folder holding two files of one name", {
  local_fake_drive_gdb(list(gdb = "g", gdb = "h"))
  expect_error(
    LandWebUtils:::.driveGdbDownloadOnce("folder", withr::local_tempdir()),
    "more than one file named .gdb."
  )
})

test_that(".driveGdbDownloadOnce stops on a folder whose name does not end in .gdb", {
  local_fake_drive_gdb(list(gdb = "g"), name = "Landbases")
  expect_error(
    LandWebUtils:::.driveGdbDownloadOnce("folder", withr::local_tempdir()),
    "must end in '.gdb'"
  )
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

test_that("buildReportingPolygons stops on a shapefile whose .dbf is cut short", {
  ## GDAL reads no shape past the last whole record of a cut .dbf, and only warns
  shp <- local_cut_shapefile(n = 3L, keep = 2L)
  expect_equal(nrow(suppressWarnings(terra::vect(shp))), 2L)

  testthat::local_mocked_bindings(.fetchVectorFile = function(...) shp)
  lyrs <- tibble::tibble(
    key = "Sub-regions",
    NAME_SHORT = "Sub",
    isTenure = FALSE,
    cross = TRUE,
    source = "drive",
    id = "x",
    labelCols = "Name",
    where = NA_character_
  )
  sa <- terra::vect("POLYGON ((-1 -1, 7 -1, 7 2, -1 2, -1 -1))", crs = "EPSG:3857")
  expect_error(
    buildReportingPolygons(sa, withr::local_tempdir(), layers = lyrs, min_area_km2 = 0),
    "2 complete records.*3 shapes"
  )
})

test_that("the Drive and URL layer readers stop on a shapefile whose .dbf is cut short", {
  testthat::local_mocked_bindings(
    drive_download_once = function(...) invisible(NULL),
    download_once = function(...) invisible(NULL),
    archive_extract_once = function(...) invisible(NULL),
    .package = "workflowtools"
  )
  dest <- withr::local_tempdir()
  local_cut_shapefile(keep = 2L, dir = file.path(dest, "1c5vkNPG81DB5wkAVT3CP_h2jFCUDqCoL"))
  expect_error(prepFMAs(dest), "2 complete records.*3 shapes")

  local_cut_shapefile(keep = 2L, dir = file.path(dest, "Ecoregions"))
  expect_error(prepEcoregionLayer(dest), "2 complete records.*3 shapes")
})

test_that(".findVectorFile finds a File Geodatabase first, and a .json", {
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "sub", "ranges.gdb"), recursive = TRUE)
  file.create(file.path(dir, "sub", "ranges.gdb", "a00000001.gdbtable"))
  file.create(file.path(dir, "gadm41_CAN_1.json"))
  f <- LandWebUtils:::.findVectorFile(dir)
  expect_match(f[[1]], "ranges\\.gdb$")
  expect_match(f[[2]], "gadm41_CAN_1\\.json$")
})

test_that("a download API's URL is saved as <slug>.zip, a plain archive URL under its own name", {
  dest <- character(0)
  testthat::local_mocked_bindings(
    download_once = function(url, destfile, ...) dest <<- c(dest, destfile),
    archive_extract_once = function(...) invisible(NULL),
    .package = "workflowtools"
  )
  dir <- withr::local_tempdir()
  LandWebUtils:::.fetchVectorFile(
    "url",
    "https://x.ca/api/file?path=%2Fa%2Fb.gdb.zip",
    dir,
    "Boreal_Caribou"
  )
  LandWebUtils:::.fetchVectorFile("url", "https://x.ca/zone/ecozone_shp.zip", dir, "Ecozones")
  expect_identical(basename(dest), c("Boreal_Caribou.zip", "ecozone_shp.zip"))
})

test_that("nationalReportingLayers: provinces, ecozones, protected areas, caribou, intact forest, none crossed", {
  lyrs <- nationalReportingLayers()
  expect_named(lyrs, c(names(reportingPolygonLayers()), "query", "layer", "dissolve", "labelMap"))
  expect_identical(
    lyrs$NAME_SHORT,
    c("Province", "Ecozone", "Protected", "BorealCaribou", "Intact")
  )
  expect_false(any(lyrs$isTenure | lyrs$cross))
  ## the ecozone row is reportingPolygonLayers()' own
  eco <- reportingPolygonLayers()
  expect_identical(lyrs$id[lyrs$NAME_SHORT == "Ecozone"], eco$id[eco$NAME_SHORT == "Ecozone"])
  expect_match(
    lyrs$query[lyrs$NAME_SHORT == "Protected"],
    "PA_OECM_DF = 1 AND BIOME = 'T'",
    fixed = TRUE
  )
  ## provinces from NRCan's CanVec (Open Government Licence - Canada), not GADM (non-commercial)
  prov <- lyrs[lyrs$NAME_SHORT == "Province", ]
  expect_match(prov$id, "canvec_1M_CA_Admin_shp\\.zip$")
  expect_identical(prov$layer, "geo_political_region_2")
  expect_identical(prov$labelCols, "juri_en")
  expect_false(any(grepl("gadm", lyrs$id, ignore.case = TRUE)))
  expect_identical(lyrs$dissolve, c(TRUE, FALSE, TRUE, FALSE, TRUE))
  ## named by two-letter code, as #8 had them
  codes <- prov$labelMap[[1]]
  expect_identical(
    unname(codes),
    c("AB", "BC", "MB", "NB", "NL", "NT", "NS", "NU", "ON", "PE", "QC", "SK", "YT")
  )
  expect_identical(
    codes[c("Quebec", "Newfoundland and Labrador")],
    c(Quebec = "QC", `Newfoundland and Labrador` = "NL")
  )
  expect_true(all(vapply(lyrs$labelMap[-1], is.null, logical(1))))
})

test_that("a layer's query selects features and labels them all one reporting unit", {
  ## protected areas: 2 PAs (one split in two parts) and an OECM; the query keeps the PAs as one unit
  v <- terra::vect(sf::st_sf(
    PA_OECM_DF = c(1L, 1L, 2L),
    BIOME = "T",
    geometry = sf::st_sfc(square(0, 0, 1), square(2, 0, 1), square(4, 0, 1), crs = 3857)
  ))
  ## a shapefile, read with OGR SQL as the CPCAD geodatabase is (a GeoPackage takes SQLite SQL, in
  ## which the geometry column must be selected by name)
  src <- file.path(withr::local_tempdir(), "pca.shp")
  terra::writeVector(v, src)
  testthat::local_mocked_bindings(.fetchVectorFile = function(...) src)
  lyrs <- tibble::tibble(
    key = "Protected Areas",
    NAME_SHORT = "Protected",
    isTenure = FALSE,
    cross = FALSE,
    source = "url",
    id = "x",
    labelCols = "label",
    where = NA_character_,
    query = "SELECT 'Protected' AS label FROM pca WHERE PA_OECM_DF = 1"
  )
  sa <- terra::vect("POLYGON ((-1 -1, 7 -1, 7 2, -1 2, -1 -1))", crs = "EPSG:3857")
  out <- buildReportingPolygons(
    sa,
    withr::local_tempdir(),
    targetCRS = "EPSG:3857",
    layers = lyrs,
    min_area_km2 = 0
  )
  expect_named(out, "Protected")
  expect_identical(unique(out$Protected$Name), "Protected")
  expect_equal(sum(terra::expanse(out$Protected, transform = FALSE)), 2)
})

test_that("the Intact row's query reads a GeoPackage as one dissolved unit", {
  ## two adjacent intact landscapes and one far away, in the IFL file's layer and columns
  v <- terra::vect(sf::st_sf(
    IFL_ID = c("NAM_1", "NAM_2", "NAM_3"),
    geometry = sf::st_sfc(square(0, 0, 1), square(1, 0, 1), square(10, 0, 1), crs = 3857)
  ))
  src <- file.path(withr::local_tempdir(), "IFL_2025.gpkg")
  terra::writeVector(v, src, layer = "IFL_2025")
  testthat::local_mocked_bindings(.fetchVectorFile = function(...) src)
  lyrs <- nationalReportingLayers()
  sa <- terra::vect("POLYGON ((-1 -1, 3 -1, 3 2, -1 2, -1 -1))", crs = "EPSG:3857")
  out <- buildReportingPolygons(
    sa,
    withr::local_tempdir(),
    targetCRS = "EPSG:3857",
    layers = lyrs[lyrs$NAME_SHORT == "Intact", ],
    min_area_km2 = 0
  )
  expect_named(out, "Intact")
  expect_equal(nrow(out$Intact), 1)
  expect_identical(out$Intact$Name, "Intact")
  expect_equal(sum(terra::expanse(out$Intact, transform = FALSE)), 2)
})

test_that("a GeoPackage served bare is not extracted", {
  extracted <- character(0)
  testthat::local_mocked_bindings(
    download_once = function(url, destfile, ...) file.create(destfile),
    archive_extract_once = function(archive, ...) extracted <<- c(extracted, archive),
    .package = "workflowtools"
  )
  dir <- withr::local_tempdir()
  f <- LandWebUtils:::.fetchVectorFile("url", "https://x.org/shp/IFL_2025.gpkg", dir, "Intact")
  expect_identical(basename(f), "IFL_2025.gpkg")
  expect_length(extracted, 0L)
  LandWebUtils:::.fetchVectorFile("url", "https://x.ca/zone/ecozone_shp.zip", dir, "Ecozones")
  expect_identical(basename(extracted), "ecozone_shp.zip")
})

test_that("a landbase gdb nested below the top level is not taken", {
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "sub", "x.gdb"), recursive = TRUE)
  file.create(file.path(dir, "y.shp"))
  expect_identical(basename(LandWebUtils:::.findLandbaseSource(dir)), "y.shp")
})

test_that("a shapefile-only search ignores a File Geodatabase", {
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "a.gdb"))
  file.create(file.path(dir, "b.shp"))
  expect_identical(basename(LandWebUtils:::.findVectorFile(dir, "shp")), "b.shp")
})

test_that("a vector file outranks a metadata .json", {
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("a_meta.json", "z.shp")))
  expect_identical(basename(LandWebUtils:::.findVectorFile(dir)[[1]]), "z.shp")
})

test_that("macOS resource folders are not read", {
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "__MACOSX", "R.gdb"), recursive = TRUE)
  dir.create(file.path(dir, "R.gdb"))
  f <- LandWebUtils:::.findVectorFile(dir)
  expect_length(f, 1L)
  expect_match(f, "/R\\.gdb$")
  expect_no_match(f, "__MACOSX")
})

## two shapefiles in one download, as in NRCan's CanVec Admin archive; the first is not the layer
local_two_shapefiles <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  mk <- function(name, label) {
    v <- terra::vect(sf::st_sf(juri_en = label, geometry = sf::st_sfc(square(0, 0, 1), crs = 3857)))
    f <- file.path(dir, paste0(name, ".shp"))
    terra::writeVector(v, f)
    f
  }
  c(mk("geopolitical_boundary_1", "wrong"), mk("geo_political_region_2", "Alberta"))
}

provinceRow <- function(layer) {
  tibble::tibble(
    key = "Provinces and Territories",
    NAME_SHORT = "Province",
    isTenure = FALSE,
    cross = FALSE,
    source = "url",
    id = "x",
    labelCols = "juri_en",
    where = NA_character_,
    layer = layer
  )
}

test_that("a layer's 'layer' picks one shapefile of several", {
  src <- local_two_shapefiles()
  testthat::local_mocked_bindings(.fetchVectorFile = function(...) src)
  sa <- terra::vect("POLYGON ((-1 -1, 2 -1, 2 2, -1 2, -1 -1))", crs = "EPSG:3857")
  out <- buildReportingPolygons(
    sa,
    withr::local_tempdir(),
    targetCRS = "EPSG:3857",
    layers = provinceRow("geo_political_region_2"),
    min_area_km2 = 0
  )
  expect_identical(unique(out$Province$Name), "Alberta")
})

test_that("a 'layer' matching no file warns and drops the layer", {
  src <- local_two_shapefiles()
  testthat::local_mocked_bindings(.fetchVectorFile = function(...) src)
  sa <- terra::vect("POLYGON ((-1 -1, 2 -1, 2 2, -1 2, -1 -1))", crs = "EPSG:3857")
  expect_warning(
    out <- buildReportingPolygons(
      sa,
      withr::local_tempdir(),
      targetCRS = "EPSG:3857",
      layers = provinceRow("no_such_layer"),
      min_area_km2 = 0
    ),
    "no vector file named 'no_such_layer'"
  )
  expect_length(out, 0L)
})

test_that("a dissolved layer counts overlapping features once", {
  v <- terra::vect(sf::st_sf(
    label = "Protected",
    geometry = sf::st_sfc(square(0, 0, 1), square(0.5, 0, 1), crs = 3857)
  ))
  src <- file.path(withr::local_tempdir(), "pa.shp")
  terra::writeVector(v, src)
  testthat::local_mocked_bindings(.fetchVectorFile = function(...) src)
  lyrs <- tibble::tibble(
    key = "Protected Areas",
    NAME_SHORT = "Protected",
    isTenure = FALSE,
    cross = FALSE,
    source = "url",
    id = "x",
    labelCols = "label",
    where = NA_character_,
    dissolve = TRUE
  )
  sa <- terra::vect("POLYGON ((-1 -1, 3 -1, 3 2, -1 2, -1 -1))", crs = "EPSG:3857")
  out <- buildReportingPolygons(
    sa,
    withr::local_tempdir(),
    targetCRS = "EPSG:3857",
    layers = lyrs,
    min_area_km2 = 0
  )
  expect_equal(nrow(out$Protected), 1)
  expect_equal(sum(terra::expanse(out$Protected, transform = FALSE)), 1.5)
})

test_that("a layer's labelMap renames its units", {
  v <- terra::vect(sf::st_sf(
    juri_en = c("Alberta", "British Columbia"),
    geometry = sf::st_sfc(square(0, 0, 1), square(1, 0, 1), crs = 3857)
  ))
  src <- file.path(withr::local_tempdir(), "prov.shp")
  terra::writeVector(v, src)
  testthat::local_mocked_bindings(.fetchVectorFile = function(...) src)
  lyrs <- provinceRow(NA_character_)
  lyrs$labelMap <- list(c(Alberta = "AB", `British Columbia` = "BC"))
  sa <- terra::vect("POLYGON ((-1 -1, 3 -1, 3 2, -1 2, -1 -1))", crs = "EPSG:3857")
  out <- buildReportingPolygons(
    sa,
    withr::local_tempdir(),
    targetCRS = "EPSG:3857",
    layers = lyrs,
    min_area_km2 = 0
  )
  expect_setequal(out$Province$Name, c("AB", "BC"))
})

test_that("a label missing from labelMap warns and is dropped", {
  v <- terra::vect(sf::st_sf(
    juri_en = c("Alberta", "Atlantis"),
    geometry = sf::st_sfc(square(0, 0, 1), square(1, 0, 1), crs = 3857)
  ))
  src <- file.path(withr::local_tempdir(), "prov.shp")
  terra::writeVector(v, src)
  testthat::local_mocked_bindings(.fetchVectorFile = function(...) src)
  lyrs <- provinceRow(NA_character_)
  lyrs$labelMap <- list(c(Alberta = "AB"))
  sa <- terra::vect("POLYGON ((-1 -1, 3 -1, 3 2, -1 2, -1 -1))", crs = "EPSG:3857")
  expect_warning(
    out <- buildReportingPolygons(
      sa,
      withr::local_tempdir(),
      targetCRS = "EPSG:3857",
      layers = lyrs,
      min_area_km2 = 0
    ),
    "no code in `labelMap` for .Atlantis."
  )
  expect_identical(out$Province$Name, "AB")
})

test_that("a source's CRS is read from the layer its query names", {
  gp <- file.path(withr::local_tempdir(), "two.gpkg")
  merc <- terra::vect(sf::st_sf(id = 1, geometry = sf::st_sfc(square(0, 0, 1), crs = 3857)))
  geo <- terra::vect(sf::st_sf(id = 1, geometry = sf::st_sfc(square(0, 0, 1), crs = 4326)))
  terra::writeVector(merc, gp, layer = "merc")
  terra::writeVector(geo, gp, layer = "geo", insert = TRUE)
  expect_true(sf::st_crs(LandWebUtils:::.sourceCrs(gp, "SELECT geom FROM geo")) == sf::st_crs(4326))
  expect_true(sf::st_crs(LandWebUtils:::.sourceCrs(gp, "")) == sf::st_crs(3857))
})

test_that("a reporting layer is read only near the study area, in the source's CRS", {
  near <- sf::st_polygon(list(rbind(c(-100, 50), c(-99, 50), c(-99, 51), c(-100, 51), c(-100, 50))))
  far <- sf::st_polygon(list(rbind(c(10, 50), c(11, 50), c(11, 51), c(10, 51), c(10, 50))))
  v <- terra::vect(sf::st_sf(
    label = c("near", "far"),
    geometry = sf::st_sfc(near, far, crs = 4326)
  ))
  src <- file.path(withr::local_tempdir(), "lyr.gpkg")
  terra::writeVector(v, src, layer = "lyr")
  nRead <- NULL
  testthat::local_mocked_bindings(.fetchVectorFile = function(...) src)
  testthat::local_mocked_bindings(
    read_vector = function(src, ...) {
      x <- terra::vect(src, ...)
      nRead <<- nrow(x)
      x
    },
    .package = "spatialutils"
  )
  sa <- terra::project(
    terra::vect(
      "POLYGON ((-100.5 49.5, -98.5 49.5, -98.5 51.5, -100.5 51.5, -100.5 49.5))",
      crs = "EPSG:4326"
    ),
    LandWebCRS
  )
  lyrs <- tibble::tibble(
    key = "Test",
    NAME_SHORT = "Test",
    isTenure = FALSE,
    cross = FALSE,
    source = "url",
    id = "x",
    labelCols = "label",
    where = NA_character_,
    query = "SELECT geom, label FROM lyr"
  )
  out <- buildReportingPolygons(
    sa,
    withr::local_tempdir(),
    targetCRS = LandWebCRS,
    layers = lyrs,
    min_area_km2 = 0
  )
  expect_equal(nRead, 1)
  expect_identical(out$Test$Name, "near")
})

test_that("the read filter is the study area, widened, in the source's CRS", {
  sa <- terra::project(
    terra::vect(
      "POLYGON ((-100.5 49.5, -98.5 49.5, -98.5 51.5, -100.5 51.5, -100.5 49.5))",
      crs = "EPSG:4326"
    ),
    LandWebCRS
  )
  src <- file.path(withr::local_tempdir(), "pt.gpkg")
  pt <- terra::vect(sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_point(c(-99, 50)), crs = 4326)))
  terra::writeVector(pt, src, layer = "lyr")
  for (x in list(sa, sf::st_as_sf(sa))) {
    f <- LandWebUtils:::.readFilter(x, src, "SELECT geom FROM lyr")
    expect_s4_class(f, "SpatVector")
    expect_true(terra::same.crs(f, "EPSG:4326"))
    expect_true(terra::relate(f, terra::project(sa, "EPSG:4326"), "covers")[1, 1])
  }
})
