## Mock the googledrive calls `.driveGdbDownloadOnce()` makes, for a Drive folder named `name`
## holding `files` (file name -> contents) plus one subfolder. `on_download(name, path)` runs inside
## each fake download, after the file is written, to stand in for another worker. Returns an
## environment whose `downloaded` records the names of the files downloaded, in order.
local_fake_drive_gdb <- function(
  files,
  name = "Landbases.gdb",
  on_download = function(name, path) NULL,
  env = parent.frame()
) {
  folderMime <- "application/vnd.google-apps.folder"
  ids <- paste0("id", seq_along(files))
  listing <- tibble::tibble(
    name = c(names(files), "subfolder"),
    id = c(ids, "idSub"),
    drive_resource = c(
      lapply(unname(files), function(x) {
        list(mimeType = "application/octet-stream", size = as.character(nchar(x, "bytes")))
      }),
      list(list(mimeType = folderMime))
    )
  )

  calls <- new.env()
  calls$downloaded <- character()
  testthat::local_mocked_bindings(
    drive_get = function(file, ...) {
      tibble::tibble(name = name, id = "folder", drive_resource = list(list(mimeType = folderMime)))
    },
    drive_ls = function(path, ...) listing,
    drive_download = function(file, path, ...) {
      i <- match(as.character(file), ids)
      calls$downloaded <- c(calls$downloaded, names(files)[[i]])
      writeBin(charToRaw(files[[i]]), path)
      on_download(names(files)[[i]], path)
      invisible(tibble::tibble(name = names(files)[[i]], local_path = path))
    },
    .package = "googledrive",
    .env = env
  )
  calls
}
