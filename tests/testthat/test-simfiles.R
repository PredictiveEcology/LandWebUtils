test_that("findOldSimFiles() searches only replicate folders", {
  ## the output path contains "rep1", and so does a folder that is not a replicate
  outDir <- file.path(withr::local_tempdir(), "prep1", "sim")
  dirs <- file.path(outDir, c("rep01", "rep10", "rep1_old", "figures"))
  fs::dir_create(dirs)
  f <- file.path(dirs, "mySimOut.qs")
  fs::file_create(f)
  Sys.setFileTime(f, as.POSIXct("2020-01-01", tz = "UTC"))

  old <- findOldSimFiles(outDir)
  expect_length(old, 2L)
  expect_setequal(basename(dirname(old)), c("rep01", "rep10"))

  ## a relative path: fs lists the folders without a leading "./"
  withr::local_dir(outDir)
  old <- findOldSimFiles(".")
  expect_length(old, 2L)
  expect_setequal(basename(dirname(old)), c("rep01", "rep10"))
})
