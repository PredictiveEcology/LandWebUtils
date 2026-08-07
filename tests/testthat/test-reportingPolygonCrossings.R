## Toy fixture: two side-by-side tenures (T1 supplied as TWO features, so the dissolve is
## exercised) under two horizontal sub-region bands that span BOTH tenures -- the arrangement
## that makes the "cross before grouping" requirement observable.
##
##   y=4 +-----------+-----------+
##       |    S2     |    S2     |
##   y=2 +-----+-----+-----------+
##       |    S1     |    S1     |
##   y=0 +-----+-----+-----------+
##      x=0   2     4           8
##       <-- T1a T1b -->  <- T2 ->
.toyTenures <- function() {
  v <- terra::vect(c(
    "POLYGON ((0 0, 2 0, 2 4, 0 4, 0 0))", ## T1, part a
    "POLYGON ((2 0, 4 0, 4 4, 2 4, 2 0))", ## T1, part b
    "POLYGON ((4 0, 8 0, 8 4, 4 4, 4 0))"  ## T2
  ))
  v$Name <- c("T1", "T1", "T2")
  terra::crs(v) <- "EPSG:3857"
  v
}

.toySubregions <- function() {
  v <- terra::vect(c(
    "POLYGON ((0 0, 8 0, 8 2, 0 2, 0 0))", ## S1: spans both tenures
    "POLYGON ((0 2, 8 2, 8 4, 0 4, 0 2))"  ## S2: spans both tenures
  ))
  v$Name <- c("S1", "S2")
  terra::crs(v) <- "EPSG:3857"
  v
}

.toyLayers <- function() {
  tibble::tribble(
    ~key         , ~NAME_SHORT , ~isTenure , ~cross , ~source , ~id  , ~labelCols ,
    "tenures"    , "FMA"       , TRUE      , FALSE  , "drive" , "a"  , "Name"    ,
    "subregions" , "SUB"       , FALSE     , TRUE   , "drive" , "b"  , "Name"    ,
    "caribou"    , "Caribou"   , FALSE     , TRUE   , "drive" , "c"  , "Name"
  )
}

.areaByName <- function(v) {
  ## transform = FALSE: planar area in CRS units, so the toy coordinates give exact areas
  ## (expanse() otherwise computes geodesic area on the ellipsoid).
  d <- data.frame(
    Name = as.character(v$Name),
    area = as.numeric(terra::expanse(v, transform = FALSE))
  )
  stats::setNames(
    stats::aggregate(area ~ Name, d, sum)$area,
    stats::aggregate(area ~ Name, d, sum)$Name
  )
}

test_that("joinReportingPolygons concatenates names and keeps the class it was given", {
  tenure <- .toyTenures()[3, ] ## T2
  sub <- .toySubregions()

  z <- joinReportingPolygons(sub, tenure)
  expect_s4_class(z, "SpatVector")
  expect_setequal(as.character(z$Name), c("S1 T2", "S2 T2"))
  expect_false("Name.1" %in% names(z))

  ## sf in -> sf out
  zsf <- joinReportingPolygons(sf::st_as_sf(sub), sf::st_as_sf(tenure))
  expect_s3_class(zsf, "sf")
  expect_setequal(zsf$Name, c("S1 T2", "S2 T2"))
})

test_that("joinReportingPolygons reduces an ACTIVE/PASSIVE side so the tenure is not repeated", {
  ## both operands are already tenure-crossed, as in the v2 triple crossings
  lb <- terra::vect(c(
    "POLYGON ((0 0, 4 0, 4 2, 0 2, 0 0))",
    "POLYGON ((0 2, 4 2, 4 4, 0 4, 0 2))"
  ))
  lb$Name <- c("ACTIVE T1", "PASSIVE T1")
  terra::crs(lb) <- "EPSG:3857"

  car <- terra::vect("POLYGON ((0 0, 4 0, 4 3, 0 3, 0 0))")
  car$Name <- "RangeA T1"
  terra::crs(car) <- "EPSG:3857"

  z <- joinReportingPolygons(lb, car)
  expect_setequal(as.character(z$Name), c("ACTIVE RangeA T1", "PASSIVE RangeA T1"))
})

test_that("joinReportingPolygons re-concatenates a Name.1/Name.2 input without geometry work", {
  x <- sf::st_as_sf(.toySubregions())
  x$Name.1 <- c("a1", "a2")
  x$Name.2 <- c("b1", "b2")
  x$Name <- NULL

  z <- joinReportingPolygons(x, NULL)
  expect_identical(z$Name, c("b1 a1", "b2 a2"))
  expect_false(any(c("Name.1", "Name.2") %in% names(z)))
})

test_that("buildCrossedReportingPolygons produces one unit per tenure x layer", {
  polys <- list(FMA = .toyTenures(), SUB = .toySubregions())

  out <- buildCrossedReportingPolygons(polys, layers = .toyLayers(), min_area_km2 = 0)
  expect_setequal(names(out), c("T1 SUB", "T2 SUB"))
  expect_setequal(as.character(out[["T1 SUB"]]$Name), c("S1 T1", "S2 T1"))
  expect_setequal(as.character(out[["T2 SUB"]]$Name), c("S1 T2", "S2 T2"))
})

test_that("crossing happens BEFORE grouping: a sub-region is not pooled across tenures", {
  polys <- list(FMA = .toyTenures(), SUB = .toySubregions())
  out <- buildCrossedReportingPolygons(polys, layers = .toyLayers(), min_area_km2 = 0)

  ## S1 spans both tenures (8 x 2 = 16). Grouping the study-area-clipped layer by name gives
  ## ONE 16-unit "S1"; crossing first splits it into each tenure's own share (8 each).
  expect_equal(.areaByName(.toySubregions())[["S1"]], 16, tolerance = 1e-6)
  expect_equal(.areaByName(out[["T1 SUB"]])[["S1 T1"]], 8, tolerance = 1e-6)
  expect_equal(.areaByName(out[["T2 SUB"]])[["S1 T2"]], 8, tolerance = 1e-6)
})

test_that("crossed areas sum back to the tenure total (nothing dropped or double-counted)", {
  tenures <- .toyTenures()
  polys <- list(FMA = tenures, SUB = .toySubregions())
  out <- buildCrossedReportingPolygons(polys, layers = .toyLayers(), min_area_km2 = 0)

  tenureArea <- .areaByName(tenures)
  for (tn in names(tenureArea)) {
    crossed <- sum(.areaByName(out[[paste(tn, "SUB")]]))
    expect_equal(crossed, tenureArea[[tn]], tolerance = 1e-6)
  }
})

test_that("a tenure's disjoint parts are dissolved before crossing", {
  ## T1 is two features; without the dissolve the crossing would report only one part
  polys <- list(FMA = .toyTenures(), SUB = .toySubregions())
  out <- buildCrossedReportingPolygons(polys, layers = .toyLayers(), min_area_km2 = 0)
  expect_equal(sum(.areaByName(out[["T1 SUB"]])), 16, tolerance = 1e-6) ## 4 x 4, not 2 x 4
})

test_that("non-intersecting and unlabelled features do not mint reporting units", {
  far <- terra::vect("POLYGON ((100 100, 101 100, 101 101, 100 101, 100 100))")
  far$Name <- "Elsewhere"
  terra::crs(far) <- "EPSG:3857"

  sub <- .toySubregions()
  sub$Name <- c("S1", NA) ## an unlabelled feature must not become "NA T1"

  polys <- list(FMA = .toyTenures(), SUB = sub, Caribou = far)
  out <- buildCrossedReportingPolygons(polys, layers = .toyLayers(), min_area_km2 = 0)

  expect_false(any(grepl("Caribou", names(out)))) ## no overlap -> dropped
  expect_setequal(as.character(out[["T1 SUB"]]$Name), "S1 T1")
  expect_false(any(grepl("NA", unlist(lapply(out, function(x) as.character(x$Name))))))
})

test_that("landbases are crossed with their tenure, and the triple is built where both exist", {
  lb <- terra::vect(c(
    "POLYGON ((0 0, 4 0, 4 2, 0 2, 0 0))",
    "POLYGON ((0 2, 4 2, 4 4, 0 4, 0 2))"
  ))
  lb$Name <- c("ACTIVE", "PASSIVE")
  terra::crs(lb) <- "EPSG:3857"

  car <- terra::vect("POLYGON ((0 0, 8 0, 8 3, 0 3, 0 0))")
  car$Name <- "RangeA"
  terra::crs(car) <- "EPSG:3857"

  polys <- list(FMA = .toyTenures(), SUB = .toySubregions(), Caribou = car)
  out <- buildCrossedReportingPolygons(polys, landbases = list("T1 LB" = lb), layers = .toyLayers(), min_area_km2 = 0)

  expect_contains(names(out), c("T1 LB", "T1 LB Caribou"))
  expect_setequal(as.character(out[["T1 LB"]]$Name), c("ACTIVE T1", "PASSIVE T1"))

  ## the triple names the status, the range, and the tenure -- each exactly once
  expect_setequal(as.character(out[["T1 LB Caribou"]]$Name), c("ACTIVE RangeA T1", "PASSIVE RangeA T1"))

  ## T2 has no landbase, so no landbase units for it
  expect_false(any(startsWith(names(out), "T2 LB")))
})

test_that(".dropRepeatedToken collapses a doubled tenure name to its last occurrence", {
  ## the ACTIVE|PASSIVE reduction only fires for those literal statuses; a source coverage
  ## labelling status differently would otherwise repeat the tenure in a partner-facing label
  expect_identical(
    LandWebUtils:::.dropRepeatedToken("Contributing T1 RangeA T1", "T1"),
    "Contributing RangeA T1"
  )
  expect_identical(LandWebUtils:::.dropRepeatedToken("ACTIVE RangeA T1", "T1"), "ACTIVE RangeA T1")
})

test_that("buildCrossedReportingPolygons warns and returns empty without a tenure layer", {
  expect_warning(
    out <- buildCrossedReportingPolygons(list(SUB = .toySubregions()), layers = .toyLayers(), min_area_km2 = 0),
    "no tenure layer"
  )
  expect_identical(out, list())
})

test_that("sliver tenures that merely graze the study area are not crossed", {
  ## a third tenure with a negligible footprint: crossing it would mint a full set of
  ## reporting units (own refCode, aggregate dir, figures) over a handful of pixels
  tenures <- .toyTenures()
  sliver <- terra::vect("POLYGON ((8 0, 8.0001 0, 8.0001 0.0001, 8 0.0001, 8 0))")
  sliver$Name <- "T3"
  terra::crs(sliver) <- "EPSG:3857"
  tenures <- rbind(tenures, sliver)

  polys <- list(FMA = tenures, SUB = .toySubregions())

  ## the toy tenures are metres-scale, so any positive threshold drops all of them; use a
  ## threshold between the sliver and the real tenures instead
  expect_message(
    out <- buildCrossedReportingPolygons(polys, layers = .toyLayers(), min_area_km2 = 1e-11),
    "sliver tenure"
  )
  expect_setequal(names(out), c("T1 SUB", "T2 SUB"))
})

test_that("only the study area's OWN tenures are crossed when `members` is given", {
  ## v10 tenure polygons overlap, so a group boundary clips in neighbouring tenures assigned
  ## to a different group; crossing them splits one partner's numbers across two reports
  polys <- list(FMA = .toyTenures(), SUB = .toySubregions())

  expect_message(
    out <- buildCrossedReportingPolygons(
      polys, layers = .toyLayers(), members = "T1", min_area_km2 = 0
    ),
    "non-member tenure"
  )
  expect_setequal(names(out), "T1 SUB")

  ## NULL (the default) keeps the previous behaviour: cross every tenure present
  out_all <- buildCrossedReportingPolygons(polys, layers = .toyLayers(), min_area_km2 = 0)
  expect_setequal(names(out_all), c("T1 SUB", "T2 SUB"))

  ## a member that is not present is simply absent, not an error
  expect_silent(suppressMessages(
    out2 <- buildCrossedReportingPolygons(
      polys, layers = .toyLayers(), members = c("T1", "T2", "T99"), min_area_km2 = 0
    )
  ))
  expect_setequal(names(out2), c("T1 SUB", "T2 SUB"))
})
