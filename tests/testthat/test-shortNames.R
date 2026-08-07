test_that("validateShortNames accepts a well-formed set", {
  expect_identical(
    validateShortNames(c("ANSR", "Caribou", "WeyCo_GrandePr")),
    c("ANSR", "Caribou", "WeyCo_GrandePr")
  )
})

test_that("validateShortNames rejects duplicates -- they silently overwrite outputs", {
  expect_error(validateShortNames(c("ANSR", "ANSR")), "duplicated short name")
  ## ...but only within a group when `by` is supplied (a tenure's landbase codes)
  expect_silent(validateShortNames(c("LB", "LB"), by = c("Sundre", "Manning")))
  expect_error(validateShortNames(c("LB", "LB"), by = c("Sundre", "Sundre")), "within a group")
})

test_that("validateShortNames rejects over-long and non-filesystem-safe names", {
  expect_error(validateShortNames(strrep("a", 17L)), "exceed 16 characters")
  expect_silent(validateShortNames(strrep("a", 16L)))

  ## anything `.slug()` would rewrite: the on-disk refCode and the figure label would drift
  expect_error(validateShortNames("has space"), "not filesystem-safe")
  expect_error(validateShortNames("Ltd."), "not filesystem-safe")
  expect_error(validateShortNames("_leading"), "not filesystem-safe")
})

test_that("validateShortNames rejects missing/empty names", {
  expect_error(validateShortNames(c("ANSR", NA)), "missing/empty")
  expect_error(validateShortNames(c("ANSR", "")), "missing/empty")
})

test_that("tenureShortNames is curated, unique, and within the length cap", {
  tn <- tenureShortNames()
  expect_named(tn, c("fma_name", "NAME_SHORT"))
  expect_gt(nrow(tn), 0L)
  expect_false(anyDuplicated(tn$fma_name) > 0L)
  ## the table validates itself on load; assert the invariants explicitly too
  expect_silent(validateShortNames(tn$NAME_SHORT))

  ## the long v10 names this exists to replace really are long
  expect_gt(max(nchar(tn$fma_name)), 60L)
})

test_that("shortNameFor looks up codes and returns NA for uncurated names", {
  expect_identical(shortNameFor("Spray Lake Sawmills (1980) Ltd."), "SprayLake")
  expect_identical(shortNameFor("Fort Providence"), "NT_FtProvidence")
  expect_identical(shortNameFor("No Such Tenure Ltd."), NA_character_)
  expect_length(shortNameFor(c("Golden TSA", "FML-2")), 2L)
})

test_that("refCodeFor slugs the layer name, and is collision-free where abbreviate() was not", {
  expect_identical(refCodeFor("lm", "ANSR"), "lm_ANSR")
  expect_identical(refCodeFor("lm", "SprayLake ANSR"), "lm_SprayLake_ANSR")
  expect_identical(
    refCodeFor("lw", "Tolko_Vand_WF_SL LB_S17 Caribou"),
    "lw_Tolko_Vand_WF_SL_LB_S17_Caribou"
  )

  ## regression (LandWeb#118): `abbreviate(x, minlength = 8)`, called one name at a time,
  ## mapped these two distinct reporting units onto the same key, so their parquet aggregates
  ## and figures overwrote each other.
  collided <- c("Canadian Forest Products Ltd. Alpine", "Crowsnest Forest Products Ltd. Alpine")
  expect_identical(
    abbreviate(collided[1], minlength = 8L),
    abbreviate(collided[2], minlength = 8L),
    ignore_attr = TRUE
  )
  expect_false(
    identical(refCodeFor("lm", "Canfor Alpine"), refCodeFor("lm", "Crowsnest Alpine"))
  )
})

test_that("every curated short name round-trips through .slug (refCode == NAME_SHORT)", {
  for (x in c(tenureShortNames()$NAME_SHORT, reportingPolygonLayers()$NAME_SHORT,
              landbaseLayers()$NAME_SHORT)) {
    expect_identical(LandWebUtils:::.slug(x), x)
  }
})
