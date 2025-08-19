utils::globalVariables(c(
  ".", ":=", ".N", ".SD", "NPixels", "proportion"
))

#' Calculate proportion of landscape occupied by each vegetation class
#'
#' This function is recursive.
#' If `poly` is a `SpatialPolygon`, then the function
#' will enter once, and convert this to a fasterized version, and pass that into
#' the function replacing `poly`.
#' It is also recursive of passed a vector of filenames for `tsf` and `vtm`.
#'
#' @param tsf A single filename, relative or absolute, pointing to a Time Since Fire raster.
#'            Can be any format that `terra` can use.
#'
#' @param vtm A single filename, relative or absolute, pointing to a Vegetation Type Map raster.
#'            Can be any format that `terra` can use.
#'
#' @param poly A single `sf` object or a factor `SpatRaster`.
#'
#' @param ageClasses A character vector with labels for age classes to bin the `tsf` times,
#'                   e.g., `c("Young", "Immature", "Mature", "Old")`. See `.ageClasses`.
#'
#' @param ageClassCutOffs A numeric vector with the endpoints for the `ageClasses`.
#'                        Should be `length(ageClasses) + 1`. See `.ageClassCutOffs`.
#'
#' @param sppEquivCol Character giving the column name to use in `sppEquiv`.
#'
#' @param sppEquiv Species equivalency table, e.g., derived from `LandR::sppEquivalencies_CA`.
#'
#' @return A `data.table` with proportion of the pixels in each vegetation class,
#'         for each given age class within each polygon.
#'
#' @export
LeadingVegTypeByAgeClass <- function(tsf, vtm, poly, ageClassCutOffs, ageClasses,
                                     sppEquivCol, sppEquiv) {
  ## main function code
  startTime <- Sys.time()
  if (utils::tail(ageClassCutOffs, 1) != Inf) {
    ageClassCutOffs <- c(ageClassCutOffs, Inf)
  }

  ## TODO: LandWeb workaround
  if (basename(vtm[1]) %in% c("CurrentConditionVTM.grd", "CurrentConditionVTM.tif")) {
    tsf <- file.path(dirname(vtm), "CurrentConditionTSF.tif")
  }

  ## prepare tsf rasters
  timeSinceFireFilesRast <- terra::rast(tsf[1])
  timeSinceFireFilesRast[] <- timeSinceFireFilesRast[]

  ## Use this when NOT in parallel
  # timeSinceFireFilesRast <- Cache(rasterToMemory, tsf[1])

  rasTsf <- terra::classify(
    timeSinceFireFilesRast,
    cbind(
      from = ageClassCutOffs[-length(ageClassCutOffs)] - 0.1,
      to = ageClassCutOffs[-1],
      seq_along(ageClasses)
    )
  )

  levels(rasTsf) <- data.frame(ID = seq_along(ageClasses), Factor = ageClasses)

  ## prepare vtm rasters
  rasVeg <- terra::rast(vtm[1])
  rasVeg[] <- rasVeg[] # 3 seconds

  splitVal <- paste0("_", 75757575, "_") # unlikely to occur for any other reason

  ## Individual species
  nas3 <- is.na(rasVeg[]) | rasVeg[] == 0
  nas1 <- is.na(rasTsf[]) | rasTsf[] == 0
  nas <- nas3 | nas1
  name1 <- as.character(pemisc::factorValues2(rasTsf, rasTsf[][!nas])[, 1]) ## TODO: test / verify use of terra w/ factorValues2
  # as.character(terra::levels(rasTsf)[[1]]$Factor)[rasTsf[][!nas]]
  name3 <- as.character(pemisc::factorValues2(rasVeg, rasVeg[][!nas])[, 1]) ## TODO: test / verify use of terra w/ factorValues2
  ## as.character(terra::levels(rasVeg)[[1]]$Factor)[rasVeg[][!nas]]
  ff <- paste(name1, name3, sep = splitVal) # 4 seconds

  ras <- terra::rast(rasVeg)
  ffFactor <- factor(ff)
  ras[!nas] <- ffFactor # 2 seconds

  eTable <- data.frame(ID = seq_along(levels(ffFactor)), VALUE = levels(ffFactor))
  types <- strsplit(as.character(eTable$VALUE), split = splitVal)
  types <- do.call(rbind, types)

  ## ensure species names all consistent (TODO: ensure this propagates)
  whMixed <- which(types[, 2] == "Mixed")
  types[, 2] <- LandR::equivalentName(types[, 2], sppEquiv, sppEquivCol)
  types[whMixed, 2] <- "Mixed"

  levels(ras) <- data.frame(eTable, ageClass = types[, 1], vegCover = types[, 2])

  ## prepare poly factor raster
  if (is(poly, "SpatialPolygons")) {
    if (!"shinyLabel" %in% colnames(poly@data)) {
      stop("poly must have a column 'shinyLabel'")
    }

    poly <- reproducible::Cache(fasterize2, rasTsf, poly, field = "polygonNum")
  }
  levs <- terra::levels(poly)[[1]]

  ## this is same, if all values present: e.g., 1, 2, 3, 4, 5 ...,
  ## but not if missing: e.g., 1, 2, 3, 5
  levs <- pemisc::factorValues2(poly, levs$ID)  ## TODO: test / verify use factorValues2
  facVals <- pemisc::factorValues2(poly, poly[], att = c("shinyLabel", "polygonNum"))

  bb <- data.table(
    zone = facVals$shinyLabel,
    polygonID = facVals$polygonNum,
    cell = seq_len(ncell(ras))
  )

  bb[, c("ageClass", "vegCover") := pemisc::factorValues2(ras, ras[][bb$cell], att = c("ageClass", "vegCover"))]
  bb <- na.omit(bb)

  ## One species at a time -- collapse polygons with same 'zone' name
  tabulated <- bb[, list(NPixels = .N), by = c("zone", "ageClass", "vegCover")] ## dedupes the zones
  tabulated[, proportion := round(NPixels / sum(NPixels), 4), by = c("zone", "vegCover")]

  ## All species -- collapse polygons with same 'zone' name
  tabulated2 <- bb[, list(NPixels = .N), by = c("zone", "ageClass")] ## dedupes the zones
  tabulated2[, proportion := round(NPixels / sum(NPixels), 4), by = c("zone")]
  set(tabulated2, NULL, "vegCover", "All species")

  tabulated <- rbindlist(list(tabulated, tabulated2), use.names = TRUE, fill = TRUE)

  ## column containing the factor names varies, so we need to search for the right one
  colID <- which(colnames(terra::levels(rasVeg)[[1]]) %in% c("category", "Factor", "VALUE"))
  coverClasses <- terra::levels(rasVeg)[[1]][[colID]]
  if (is.factor(coverClasses)) {
    coverClasses <- levels(coverClasses)
  }

  coverClasses <- as.character(coverClasses)

  emptyID <- which(coverClasses == "")
  if (length(emptyID)) {
    coverClasses <- coverClasses[-emptyID]
  }

  if (!("All species" %in% levels(coverClasses))) {
    coverClasses <- c(coverClasses, "All species")
  }

  ## ensure species names all consistent (TODO: ensure this propagates)
  whAll <- which(coverClasses == "All species")
  whMixed <- which(coverClasses == "Mixed")
  coverClasses <- LandR::equivalentName(coverClasses, sppEquiv, sppEquivCol)
  coverClasses[c(whAll, whMixed)] <- c("All species", "Mixed")

  allCombos <- expand.grid(
    ageClass = ageClasses,
    vegCover = coverClasses,
    zone = unique(levs$shinyLabel),
    stringsAsFactors = FALSE
  )
  # allCombos$polygonID <- match(allCombos$zone, levs$shinyLabel)
  data.table::setDT(allCombos)

  tabulated <- merge(
    tabulated,
    allCombos,
    # by = c("zone", "vegCover", "ageClass", "polygonID"),
    by = c("zone", "vegCover", "ageClass"),
    all.y = TRUE
  )
  ## fill in zeros where there is no value
  tabulated[is.na(proportion), proportion := 0]
  set(
    tabulated,
    NULL,
    "label",
    paste(
      tabulated$ageClass,
      paste(gsub(basename(dirname(tsf[1])), pattern = "\\.", replacement = ""),
        basename(tsf[1]),
        sep = "_"
      ),
      sep = "."
    )
  )

  endTime <- Sys.time()
  message("    Leading cover calculation took ", format(endTime - startTime, digits = 2))

  return(tabulated)
}
