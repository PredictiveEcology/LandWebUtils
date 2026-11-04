## Write `n` unit squares to `<dir>/<name>.shp` with the attribute columns `cols` (values recycled
## to `n`), then cut its .dbf to `keep` complete records: the state a writer leaves a .dbf in when
## it stops at the 2 GiB dBASE size limit, with the header still counting all `n`. The CLS S17 and
## S20 landbase shapefiles are cut like this. Returns the .shp path.
local_cut_shapefile <- function(
  n = 3L,
  keep = n,
  cols = list(Name = c("A", "B", "C")),
  dir = withr::local_tempdir(.local_envir = env),
  name = "layer",
  env = parent.frame()
) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  shp <- file.path(dir, paste0(name, ".shp"))

  v <- do.call(
    rbind,
    lapply(seq_len(n) - 1L, function(i) terra::as.polygons(terra::ext(2 * i, 2 * i + 1, 0, 1)))
  )
  terra::crs(v) <- "EPSG:3857"
  for (k in names(cols)) {
    v[[k]] <- rep_len(cols[[k]], n)
  }
  terra::writeVector(v, shp)

  if (keep < n) {
    dbf <- sub("\\.shp$", ".dbf", shp)
    bytes <- readBin(dbf, "raw", file.size(dbf))
    hdrLen <- readBin(bytes[9:10], "integer", size = 2L, signed = FALSE, endian = "little")
    recLen <- readBin(bytes[11:12], "integer", size = 2L, signed = FALSE, endian = "little")
    writeBin(bytes[seq_len(hdrLen + keep * recLen)], dbf)
  }

  shp
}
