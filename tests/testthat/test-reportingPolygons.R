test_that("reportingPolygonLayers returns the expected schema", {
  lyrs <- reportingPolygonLayers()
  expect_named(lyrs, c("key", "source", "id", "labelCol"))
  expect_setequal(unique(lyrs$source), c("drive", "url"))
  expect_contains(lyrs$key, "FMA Boundaries Updated")
})
