test_that("equivalentName maps a value to the requested column", {
  df <- data.frame(
    KNN = c("Pice_gla", "Pice_mar", "Popu_tre"),
    Latin = c("Picea glauca", "Picea mariana", "Populus tremuloides")
  )
  expect_equal(equivalentName("Pice_mar", df, "Latin"), "Picea mariana")
  expect_equal(equivalentName("Picea glauca", df, "KNN"), "Pice_gla")
})

test_that("equivalentName returns NA for an unknown value", {
  df <- data.frame(KNN = c("Pice_gla", "Pice_mar"), Latin = c("Picea glauca", "Picea mariana"))
  expect_equal(equivalentName("nope", df, "Latin"), NA_character_)
})
