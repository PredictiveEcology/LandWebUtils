test_that("extractFMA matches the canonical Name column when present", {
  fmas <- data.frame(
    Name = c("SprayLake", "AlPac"),
    FMA_NAME = c("Spray Lake Sawmills (1980) Ltd.", "Alberta-Pacific Forest Industries Inc.")
  )
  out <- extractFMA(fmas, "SprayLake")
  expect_equal(nrow(out), 1L)
  expect_equal(out$Name, "SprayLake")
})

test_that("extractFMA falls back to LandWebStudyAreas$Description vs FMA_NAME (v10 shapefile)", {
  ## v10 shapefile: canonical short names absent from `Name`; full name in `FMA_NAME`.
  fmas <- data.frame(
    Name = c(NA_character_, NA_character_),
    FMA_NAME = c("Spray Lake Sawmills (1980) Ltd.", "Sundre Forest Products Inc.")
  )
  out <- extractFMA(fmas, "SprayLake")
  expect_equal(nrow(out), 1L)
  expect_equal(out$FMA_NAME, "Spray Lake Sawmills (1980) Ltd.")
})

test_that("extractFMA returns no rows for a name absent from both Name and Description", {
  fmas <- data.frame(
    Name = c(NA_character_, NA_character_),
    FMA_NAME = c("Spray Lake Sawmills (1980) Ltd.", "Sundre Forest Products Inc.")
  )
  expect_equal(nrow(extractFMA(fmas, "Mistik")), 0L)
})
