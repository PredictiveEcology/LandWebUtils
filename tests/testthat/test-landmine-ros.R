## Fixture: a vegetation type map whose attribute table is 1..n, plus matching inputs.
## `labels` become the RAT's labels; "Mixed" (or anything absent from `sppEquiv`) becomes the
## `mixed` fuel type.
.ros_fixture <- function(labels = c("Pice_mar", "Pinu_spp", "Popu_spp", "Abie_spp", "Mixed"),
                         vegType = rep(seq_along(labels), each = 4L),
                         tsf = rep(c(10, 60, 200, 900), times = length(labels)),
                         flammable = 1L) {
  n <- length(vegType)
  vegTypeMap <- terra::rast(nrows = n, ncols = 1L, vals = vegType)
  levels(vegTypeMap) <- data.frame(ID = seq_along(labels), Species = labels)

  list(
    vegTypeMap = vegTypeMap,
    tsf = rep_len(tsf, n),
    flammable = rep_len(as.integer(flammable), n),
    sppEquiv = .ros_sppEquiv(),
    ROSTable = .ros_table()
  )
}

.ros_sppEquiv <- function() {
  spp <- c("Pice_mar", "Pice_gla", "Pinu_spp", "Popu_spp", "Abie_spp", "Thuj_pli")
  dt <- data.table::data.table(LandWeb = spp, LandR = spp)
  dt[, LandMine := landmine_known_species()[LandWeb]]
  dt[]
}

.ros_table <- function() {
  data.table::rbindlist(list(
    list("immature_young", "decid", 6L),
    list("mature", "decid", 9L),
    list("immature_young", "mixed", 12L),
    list("immature", "pine", 14L),
    list("mature", "mixed", 17L),
    list("immature_young", "softwood", 18L),
    list("immature_young", "spruce", 20L),
    list("mature", "pine", 21L),
    list("young", "pine", 22L),
    list("mature", "softwood", 27L),
    list("mature", "spruce", 30L)
  )) |>
    data.table::setnames(old = 1:3, new = c("age", "leading", "ros"))
}

## raster values whose attribute-table label matches no species, i.e. the `mixed` fuel type
.ros_mixed_pixels <- function(f) {
  rat <- terra::levels(f$vegTypeMap)[[1]]
  rat[[1]][!(as.character(rat[[2]]) %in% f$sppEquiv$LandWeb)]
}

.ros_call <- function(f, ROSother = 30L, ...) {
  landmine_fire_ros(
    vegTypeMap = f$vegTypeMap, rstTimeSinceFire = f$tsf, flammableMap = f$flammable,
    ROSTable = f$ROSTable, sppEquiv = f$sppEquiv, sppEquivCol = "LandWeb",
    ROSother = ROSother, ...
  )
}

## The LandMine module's `fireROS()` body, kept verbatim as the equivalence oracle. `sim`/`mod`/
## `P(sim)` reads are replaced by arguments; nothing else is changed.
.fireros_inline <- function(vegTypeMap, tsf, flammable, ROSTable, sppEquiv, sppEquivCol,
                            ROSother, knownSpecies, ROStype = "default") {
  ROS <- rep(NA_integer_, terra::ncell(vegTypeMap))

  vegType <- terra::values(vegTypeMap, mat = FALSE)
  vegTypes <- data.table::data.table(terra::levels(vegTypeMap)[[1]])

  sppNames <- LandWebUtils:::equivalentName(as.character(vegTypes[[2]]), sppEquiv, sppEquivCol)
  suppressWarnings({
    onRaster <- data.table::rbindlist(list(
      list("mixed", which(is.na(sppNames))),
      list("spruce", grep(sppNames, pattern = "Pice")),
      list("pine", grep(sppNames, pattern = "Pinu")),
      list("decid", grep(sppNames, pattern = "Popu")),
      list("softwood", grep(sppNames, pattern = "Pice|Pinu|Popu", invert = TRUE))
    ))
  })
  onRaster <- stats::na.omit(unique(onRaster, by = "V2")) |>
    data.table::setnames(old = 1:2, new = c("leading", "pixelValue")) |>
    data.table::setkeyv("pixelValue")
  onRaster[, species := sppNames]

  sppEquiv <- sppEquiv[, c("LandMine", "LandWeb")][, leading := knownSpecies[LandWeb]] |>
    stats::na.omit(on = "LandMine")
  sppEquiv <- sppEquiv[onRaster, on = c("LandMine" = "leading", "LandWeb" = "species")] |>
    unique()

  sppEquivHere <- unique(stats::na.omit(sppEquiv$LandWeb))
  haveAllKnown <- sppEquivHere %in% names(knownSpecies)
  if (!all(haveAllKnown)) {
    stop("LandMine only has rate of spread burn rates for\n",
         paste(names(knownSpecies), collapse = ", "),
         "\nMissing rate of spread for ", paste(sppEquivHere[!haveAllKnown], collapse = ", "))
  }

  sppEquiv <- unique(sppEquiv, by = c("LandMine", "leading", "pixelValue"))
  sppEquiv <- sppEquiv[ROSTable, on = "leading", allow.cartesian = TRUE, nomatch = NULL]
  sppEquiv <- sppEquiv[, c("leading", "age", "ros", "pixelValue")]
  sppEquiv <- unique(sppEquiv, by = c("age", "leading", "pixelValue"))

  sppEquiv[, used := "no"]
  sppEquiv[(used == "no") & grepl("(^|_)mature", age), used := "mature"]
  sppEquiv[(used == "no") & grepl("(^|_)immature", age), used := "immature"]
  sppEquiv[(used == "no") & grepl("(^|_)young", age), used := "young"]
  data.table::setkeyv(sppEquiv, "used")

  cuts <- list()
  if (!any(grepl("_mature$|^mature_|_mature_", sppEquiv$age))) {
    cuts[[1]] <- tsf > 120
  } else {
    cuts[[1]] <- !is.na(tsf)
  }

  if (!any(grepl("_immature$|^immature_|_immature_", sppEquiv$age))) {
    cuts[[2]] <- tsf > 40 & tsf <= 120
  } else {
    cuts[[2]] <- tsf <= 120
  }

  cuts[[3]] <- tsf <= 40

  if (!all(sppEquiv["mature"]$pixelValue %in% vegTypes[[1]])) {
    cuts[[1]] <- cuts[[1]] & vegType %in% sppEquiv["mature"]$pixelValue
  }

  if (!all(sppEquiv["immature"]$pixelValue %in% vegTypes[[1]])) {
    cuts[[2]] <- cuts[[2]] & vegType %in% sppEquiv["immature"]$pixelValue
  }

  if (all(sppEquiv["young"]$pixelValue %in% vegTypes[[1]])) {
    cuts[[3]] <- cuts[[3]] & vegType %in% sppEquiv["young"]$pixelValue
  }

  mature <- which(cuts[[1]])
  immature <- which(cuts[[2]])
  young <- which(cuts[[3]])

  if (length(mature)) {
    ROS[mature] <- sppEquiv["mature"]$ros[match(vegType[mature], sppEquiv["mature"]$pixelValue)]
  }
  if (length(immature)) {
    ROS[immature] <- sppEquiv["immature"]$ros[match(vegType[immature], sppEquiv["immature"]$pixelValue)]
  }
  if (length(young)) {
    ROS[young] <- sppEquiv["young"]$ros[match(vegType[young], sppEquiv["young"]$pixelValue)]
  }

  ROSnonflam <- switch(ROStype,
    burny = ROSTable[leading == "decid" & age == "immature_young", ros],
    NA_integer_
  )

  ROS[flammable == 1L & is.na(ROS)] <- as.integer(ROSother)
  ROS[flammable == 0L | is.na(flammable)] <- as.integer(ROSnonflam)

  ROS
}

## ---- the species/fuel-type mapping -------------------------------------------------------------

test_that("knownSpecies maps every species to one of the four table fuel types", {
  ks <- landmine_known_species()
  expect_named(ks, c(
    "Abie_spp", "Lari_spp", "Pice_gla", "Pice_mar", "Pinu_spp",
    "Popu_spp", "Pseu_men", "Thuj_pli", "Tsug_het"
  ))
  expect_setequal(unique(unname(ks)), c("softwood", "decid", "spruce", "pine"))
  ## the fifth fuel type in the ROSTable is unreachable from here -- see the mixedwood test
  expect_false("mixed" %in% ks)
})

test_that("both larch and poplar are 'decid', and all three firs are 'softwood'", {
  ks <- landmine_known_species()
  expect_equal(unname(ks[c("Lari_spp", "Popu_spp")]), c("decid", "decid"))
  expect_equal(unname(ks[c("Abie_spp", "Pseu_men", "Tsug_het")]), rep("softwood", 3L))
})

## ---- the core assignment ------------------------------------------------------------------------

test_that("each vegetation type gets its own table rate at each age class", {
  f <- .ros_fixture(labels = c("Pice_mar", "Pinu_spp", "Popu_spp", "Abie_spp"),
                    vegType = rep(1:4, each = 3L),
                    tsf = rep(c(10, 60, 200), times = 4L))
  ros <- .ros_call(f)

  ## young, immature, mature for spruce / pine / decid / softwood
  expect_equal(ros, c(
    20L, 20L, 30L, ## spruce: immature_young, immature_young, mature
    22L, 14L, 21L, ## pine: young, immature, mature
    6L, 6L, 9L, ## decid: immature_young, immature_young, mature
    18L, 18L, 27L ## softwood: immature_young, immature_young, mature
  ))
})

test_that("MIXEDWOOD receives the table's mixed rates", {
  f <- .ros_fixture(labels = c("Pice_mar", "Mixed"), vegType = c(1L, 1L, 2L, 2L),
                    tsf = c(60, 200, 60, 200))
  ros <- .ros_call(f)

  expect_equal(ros[1:2], c(20L, 30L)) ## spruce
  expect_equal(ros[3:4], c(12L, 17L)) ## mixed: immature_young, mature
})

test_that("the mixed rows of the ROSTable are live entries", {
  f <- .ros_fixture(labels = c("Pice_mar", "Mixed"), vegType = c(2L, 2L), tsf = c(60, 200))
  withMixed <- .ros_call(f)
  f$ROSTable <- f$ROSTable[leading != "mixed"]
  ## without them, mixedwood has no rate and falls through to ROSother
  expect_equal(withMixed, c(12L, 17L))
  expect_equal(.ros_call(f), c(30L, 30L))
})

test_that("LARCH still falls through to ROSother (a related defect, deliberately NOT fixed)", {
  ## `Lari_spp` matches none of Pice/Pinu/Popu so the pattern-based typing calls it "softwood",
  ## while `knownSpecies` calls it "decid"; the two disagree and the join finds no match.
  ## Filling `leading` from the pattern-derived value would override the authoritative mapping.
  f <- .ros_fixture(labels = c("Lari_spp", "Popu_spp"), vegType = c(1L, 2L), tsf = c(60, 60))
  ## larch must be IN sppEquiv, or it is simply an unmatched entry and becomes `mixed`
  f$sppEquiv <- rbind(f$sppEquiv, data.table::data.table(
    LandWeb = "Lari_spp", LandR = "Lari_spp",
    LandMine = landmine_known_species()[["Lari_spp"]]
  ))
  ros <- .ros_call(f)
  expect_equal(ros[[1]], 30L) ## larch -> ROSother, NOT decid's 6L
  expect_equal(ros[[2]], 6L) ## poplar, the same fuel type per knownSpecies, is assigned
})

## ---- age class boundaries -----------------------------------------------------------------------

test_that("the age cutoffs are inclusive upper bounds", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = rep(1L, 6L),
                    tsf = c(0, 39, 40, 41, 120, 121))
  ## pine is the only fuel type with three distinct rates, so the class is legible from the value
  expect_equal(.ros_call(f), c(22L, 22L, 22L, 14L, 14L, 21L))
})

test_that("ageCutoffs moves the boundaries", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = rep(1L, 3L), tsf = c(30, 80, 200))
  expect_equal(.ros_call(f), c(22L, 14L, 21L))
  expect_equal(.ros_call(f, ageCutoffs = c(20, 100)), c(14L, 14L, 21L))
  expect_equal(.ros_call(f, ageCutoffs = c(100, 500)), c(22L, 22L, 14L))
})

test_that("ageCutoffs must be length 2 and increasing", {
  f <- .ros_fixture()
  expect_snapshot(error = TRUE, .ros_call(f, ageCutoffs = 40))
  expect_snapshot(error = TRUE, .ros_call(f, ageCutoffs = c(120, 40)))
})

test_that("stands older than the self-test's 999-year bin are still assigned (previously unchecked)", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = rep(1L, 2L), tsf = c(998, 5000))
  expect_equal(.ros_call(f), c(21L, 21L))
})

test_that("time since fire of exactly zero is young (previously unchecked)", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = 1L, tsf = 0)
  expect_equal(.ros_call(f), 22L)
})

test_that("NA time since fire yields the flammable default, not a table rate", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = rep(1L, 2L), tsf = c(NA, 10))
  expect_equal(.ros_call(f), c(30L, 22L))
})

## ---- compound age labels ------------------------------------------------------------------------

test_that("a compound label is claimed by the OLDEST class it names", {
  f <- .ros_fixture(labels = "Pice_mar", vegType = rep(1L, 3L), tsf = c(60, 100, 200))
  ## spruce is "immature_young" + "mature": both immature ages share 20L
  expect_equal(.ros_call(f), c(20L, 20L, 30L))
})

test_that("a table with 'mature_immature' drops the upper age cut entirely", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = rep(1L, 3L), tsf = c(10, 60, 200))
  f$ROSTable <- data.table::rbindlist(list(
    list("young", "pine", 22L),
    list("mature_immature", "pine", 15L),
    list("mature", "spruce", 30L)
  )) |> data.table::setnames(c("age", "leading", "ros"))

  ## every non-NA pixel is now mature-eligible; the young class still overrides where it applies
  expect_equal(.ros_call(f), c(22L, 15L, 15L))
})

## ---- flammability -------------------------------------------------------------------------------

test_that("non-flammable pixels override any forest rate already assigned", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = rep(1L, 3L), tsf = c(10, 60, 200),
                    flammable = c(1L, 0L, NA))
  expect_equal(.ros_call(f), c(22L, NA_integer_, NA_integer_))
})

test_that("ROStype = 'burny' gives non-flammable pixels the young deciduous rate", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = rep(1L, 2L), tsf = c(10, 10),
                    flammable = c(1L, 0L))
  expect_equal(.ros_call(f, ROStype = "burny"), c(22L, 6L))
})

test_that("flammable pixels with no table rate get ROSother", {
  ## drop softwood from the table, so Abie_spp has no rate at any age
  f <- .ros_fixture(labels = c("Pinu_spp", "Abie_spp"), vegType = c(1L, 2L), tsf = c(10, 10))
  f$ROSTable <- f$ROSTable[leading != "softwood"]
  expect_equal(.ros_call(f), c(22L, 30L))
})

## ---- the young-class species filter ---------------------------------------------------------

test_that("legacy and always agree whenever some type present has a young rate", {
  ## pine is the only fuel type with a "young" row in the default table
  f <- .ros_fixture(labels = c("Pinu_spp", "Pice_mar", "Popu_spp"),
                    vegType = rep(1:3, each = 3L),
                    tsf = rep(c(10, 60, 200), times = 3L))
  expect_equal(.ros_call(f, youngFilter = "legacy"), .ros_call(f, youngFilter = "always"))
})

test_that("with NO young rate in play, legacy wipes young stands to ROSother", {
  ## no pine, so nothing in the lookup is classified "young"; legacy then skips the filter and
  ## assigns NA over the immature rate these pixels already had
  f <- .ros_fixture(labels = c("Pice_mar", "Popu_spp"), vegType = c(1L, 1L, 2L, 2L),
                    tsf = c(10, 60, 10, 60))

  expect_equal(.ros_call(f, youngFilter = "legacy"), c(30L, 20L, 30L, 6L))
  expect_equal(.ros_call(f, youngFilter = "always"), c(20L, 20L, 6L, 6L))
})

test_that("the filter never changes the mature or immature classes", {
  f <- .ros_fixture(labels = c("Pice_mar", "Popu_spp"), vegType = c(1L, 2L, 1L, 2L),
                    tsf = c(60, 60, 200, 200))
  expect_equal(.ros_call(f, youngFilter = "legacy"), .ros_call(f, youngFilter = "always"))
})

## ---- input validation ---------------------------------------------------------------------------

test_that("a vegetation map with no attribute table is a clear error, not a silent NA map", {
  f <- .ros_fixture()
  levels(f$vegTypeMap) <- NULL
  expect_snapshot(error = TRUE, .ros_call(f))
})

test_that("an attribute table with duplicate or missing values is refused", {
  ## terra refuses duplicate IDs itself, so this guard is defensive; test it directly
  ok <- data.table::data.table(ID = 1:3, Species = c("a", "b", "c"))
  expect_true(LandWebUtils:::.ros_check_rat(ok))
  expect_snapshot(error = TRUE, LandWebUtils:::.ros_check_rat(
    data.table::data.table(ID = c(1L, 1L, 3L), Species = c("a", "b", "c"))
  ))
  expect_snapshot(error = TRUE, LandWebUtils:::.ros_check_rat(
    data.table::data.table(ID = c(1L, NA, 3L), Species = c("a", "b", "c"))
  ))
})

test_that("THE ROW-ORDER CONTRACT: values come from the table, not from row position", {
  ## `LandR::vegTypeMapGenerator` writes the table in value order, but a
  ## writeRaster()/rast() round trip can return the rows permuted while keeping each ID
  ## with its own label. Indexing by position then gives whole species the wrong fuel type.
  labels <- c("Pice_mar", "Pinu_spp", "Popu_spp", "Abie_spp")
  n <- 16L
  set.seed(7)
  vegType <- sample(seq_along(labels), n, replace = TRUE)
  tsf <- sample(0:250, n, replace = TRUE)

  f <- .ros_fixture(labels = labels, vegType = vegType, tsf = tsf)
  ordered <- .ros_call(f)

  ## same map, same ID-to-label pairing, rows permuted
  perm <- c(4L, 2L, 1L, 3L)
  levels(f$vegTypeMap) <- data.frame(ID = perm, Species = labels[perm])
  expect_equal(.ros_call(f), ordered)

  ## the module's inline version, which indexed by position, does NOT survive this
  shuffled <- .fireros_inline(
    f$vegTypeMap, f$tsf, f$flammable, f$ROSTable, data.table::copy(f$sppEquiv),
    "LandWeb", 30L, landmine_known_species()
  )
  expect_false(identical(shuffled, ordered))
})

test_that("rasters on a different grid are refused", {
  f <- .ros_fixture()
  f$tsf <- f$tsf[-1L]
  expect_snapshot(error = TRUE, .ros_call(f))
})

test_that("a species with no rate of spread is named in the error", {
  f <- .ros_fixture(labels = "Pice_mar")
  f$sppEquiv <- rbind(f$sppEquiv, data.table::data.table(
    LandWeb = "Quer_spp", LandR = "Quer_spp", LandMine = NA_character_
  ))
  f$vegTypeMap <- terra::rast(nrows = 2L, ncols = 1L, vals = c(1L, 2L))
  levels(f$vegTypeMap) <- data.frame(ID = 1:2, Species = c("Pice_mar", "Quer_spp"))
  f$tsf <- c(10, 10)
  f$flammable <- c(1L, 1L)
  expect_error(.ros_call(f), "Missing rate of spread for")
})

test_that("ROSother is validated against the table's range and against mature spruce", {
  f <- .ros_fixture()
  expect_snapshot(error = TRUE, .ros_call(f, ROSother = 500L))
  expect_snapshot(error = TRUE, .ros_call(f, ROSother = 20L))
  expect_no_error(.ros_call(f, ROSother = 29L)) ## within 5% of 30
})

test_that("a table with no mature spruce row cannot validate ROSother", {
  f <- .ros_fixture()
  f$ROSTable <- f$ROSTable[!(leading == "spruce" & age == "mature")]
  expect_snapshot(error = TRUE, .ros_call(f))
})

test_that("'burny' needs an immature_young deciduous row", {
  f <- .ros_fixture()
  f$ROSTable <- f$ROSTable[!(leading == "decid" & age == "immature_young")]
  expect_snapshot(error = TRUE, .ros_call(f, ROStype = "burny"))
})

test_that("a SpatRaster and a bare vector give the same answer", {
  f <- .ros_fixture()
  fromVector <- .ros_call(f)
  tsfMap <- terra::rast(f$vegTypeMap)
  terra::values(tsfMap) <- f$tsf
  flamMap <- terra::rast(f$vegTypeMap)
  terra::values(flamMap) <- f$flammable
  fromRaster <- landmine_fire_ros(
    vegTypeMap = f$vegTypeMap, rstTimeSinceFire = tsfMap, flammableMap = flamMap,
    ROSTable = f$ROSTable, sppEquiv = f$sppEquiv, sppEquivCol = "LandWeb", ROSother = 30L
  )
  expect_equal(fromRaster, fromVector)
})

## ---- the self-test ------------------------------------------------------------------------------

test_that("the self-test catches a rate of spread that disagrees with the lookup", {
  f <- .ros_fixture(labels = "Pinu_spp", vegType = rep(1L, 3L), tsf = c(10, 60, 200))
  lookup <- LandWebUtils:::.ros_lookup(
    data.table::data.table(terra::levels(f$vegTypeMap)[[1]]),
    f$sppEquiv, "LandWeb", f$ROSTable, landmine_known_species()
  )
  cuts <- LandWebUtils:::.ros_age_cuts(lookup, f$tsf, rep(1L, 3L), 1L, c(40, 120), "legacy")

  good <- c(22L, 14L, 21L)
  expect_true(LandWebUtils:::.ros_self_test(good, rep(1L, 3L), f$tsf, cuts, lookup, c(40, 120)))
  expect_error(
    LandWebUtils:::.ros_self_test(c(99L, 14L, 21L), rep(1L, 3L), f$tsf, cuts, lookup, c(40, 120)),
    "failed its test"
  )
})

test_that("assertions = FALSE skips the self-test", {
  f <- .ros_fixture()
  expect_equal(.ros_call(f, assertions = FALSE), .ros_call(f, assertions = TRUE))
})

## ---- equivalence with the module's inline implementation -----------------------------------------

test_that("EQUIVALENCE with the module's inline fireROS over randomised landscapes", {
  set.seed(4321)
  pool <- c("Pice_mar", "Pice_gla", "Pinu_spp", "Popu_spp", "Abie_spp", "Thuj_pli", "Mixed")

  for (i in seq_len(200L)) {
    labels <- sample(pool, size = sample(2:7, 1L))
    n <- sample(10:60, 1L)
    vegType <- sample(seq_along(labels), n, replace = TRUE)
    tsf <- sample(c(0:250, NA), n, replace = TRUE)
    flammable <- sample(c(0L, 1L, NA), n, replace = TRUE, prob = c(0.2, 0.7, 0.1))
    ROStype <- sample(c("default", "burny"), 1L)

    f <- .ros_fixture(labels = labels, vegType = vegType, tsf = tsf, flammable = flammable)

    want <- .fireros_inline(
      f$vegTypeMap, f$tsf, f$flammable, f$ROSTable, data.table::copy(f$sppEquiv),
      "LandWeb", 30L, landmine_known_species(), ROStype
    )
    got <- .ros_call(f, ROStype = ROStype, youngFilter = "legacy")

    ## The two diverge by design on mixedwood -- the oracle drops it, this assigns the table's
    ## rates -- so the claim under test is that mixedwood is the ONLY place they differ.
    ## The mixed rates themselves are checked deterministically above.
    mixedPx <- .ros_mixed_pixels(f)
    keep <- !(terra::values(f$vegTypeMap, mat = FALSE) %in% mixedPx)
    expect_identical(got[keep], want[keep], info = paste("iteration", i))
  }
})

test_that("EQUIVALENCE holds when a fuel type is missing from the landscape", {
  set.seed(99)
  for (drop in c("Pice_mar", "Pinu_spp", "Popu_spp", "Abie_spp")) {
    labels <- setdiff(c("Pice_mar", "Pinu_spp", "Popu_spp", "Abie_spp"), drop)
    n <- 24L
    f <- .ros_fixture(labels = labels, vegType = sample(seq_along(labels), n, replace = TRUE),
                      tsf = sample(0:250, n, replace = TRUE), flammable = 1L)
    want <- .fireros_inline(
      f$vegTypeMap, f$tsf, f$flammable, f$ROSTable, data.table::copy(f$sppEquiv),
      "LandWeb", 30L, landmine_known_species()
    )
    expect_identical(.ros_call(f, youngFilter = "legacy"), want, info = drop)
  }
})

test_that("the lookup does not modify the caller's sppEquiv by reference", {
  f <- .ros_fixture()
  before <- data.table::copy(f$sppEquiv)
  .ros_call(f)
  expect_identical(f$sppEquiv, before)
})
