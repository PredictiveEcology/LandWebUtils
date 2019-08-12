if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ":=", ".N", ".SD", "NPixels", "proportion"))
}

#' Calculate proportion of landscape occupied by each vegetation class
#'
#' This function is recursive.
#' If \code{poly} is a \code{SpatialPolygon}, then the function
#' will enter once, and convert this to a fasterized version, and pass that into
#' the function replacing \code{poly}.
#' It is also recursive of passed a vector of filenames for \code{tsf} and \code{vtm}.
#'
#' @param tsf A single filename, relative or absolute, pointing to a Time Since Fire raster.
#'            Can be any format that \code{raster} can use.
#' @param vtm A single filename, relative or absolute, pointing to a Vegetation Type Map raster.
#'            Can be any format that \code{raster} can use.
#' @param poly A single \code{SpatialPolygonsDataFrame} object or a factor \code{RasterLayer}.
#'             This layer MUST have a column labelled \code{shinyLabel}
#' @param ageClasses A character vector with labels for age classes to bin the \code{tsf} times,
#'                   e.g., \code{c("Young", "Immature", "Mature", "Old")}
#' @param ageClassCutOffs A numeric vector with the endpoints for the \code{ageClasses}.
#'                        Should be \code{length(ageClasses) + 1}
#'
#' @param sppEquivCol TODO: description needed
#'
#' @param sppEquiv TODO: description needed
#'
#' @return A \code{data.table} with proportion of the pixels in each vegetation class,
#'         for each given age class within each polygon.
#'
#' @export
#' @importFrom data.table set setDT rbindlist
#' @importFrom raster factorValues ncell
#' @importFrom stats na.omit
#' @importFrom utils tail
LeadingVegTypeByAgeClass <- function(tsf, vtm, poly, ageClassCutOffs, ageClasses,
                                     sppEquivCol, sppEquiv) {
  # main function code
  startTime <- Sys.time()
  if (tail(ageClassCutOffs, 1) != Inf)
    ageClassCutOffs <- c(ageClassCutOffs, Inf)

  # prepare tsf rasters
  if (basename(vtm[1]) == "CurrentConditionVTM.tif") ## TODO: LandWeb workaround
    tsf <- file.path(dirname(vtm), "CurrentConditionTSF.tif")

  timeSinceFireFilesRast <- raster(tsf[1])
  timeSinceFireFilesRast[] <- timeSinceFireFilesRast[]

  # Use this when NOT in parallel
  #timeSinceFireFilesRast <- Cache(rasterToMemory, tsf[1])

  rasTsf <- reclassify(
    timeSinceFireFilesRast,
    cbind(
      from = ageClassCutOffs[-length(ageClassCutOffs)] - 0.1,
      to = ageClassCutOffs[-1],
      seq_along(ageClasses)
    )
  )

  levels(rasTsf) <- data.frame(ID = seq_along(ageClasses), Factor = ageClasses)

  # prepare vtm rasters
  rasVeg <- raster(vtm[1])
  rasVeg[] <- rasVeg[] # 3 seconds

  splitVal <- paste0("_", 75757575, "_") # unlikely to occur for any other reason

  # Individual species
  nas3 <- is.na(rasVeg[]) | rasVeg[] == 0
  nas1 <- is.na(rasTsf[]) | rasTsf[] == 0
  nas <- nas3 | nas1
  name1 <- as.character(factorValues(rasTsf, rasTsf[][!nas])[, 1])
  #as.character(raster::levels(rasTsf)[[1]]$Factor)[rasTsf[][!nas]]
  name3 <- as.character(factorValues(rasVeg, rasVeg[][!nas])[, 1])
  #as.character(raster::levels(rasVeg)[[1]]$Factor)[rasVeg[][!nas]]
  ff <- paste(name1, name3, sep = splitVal) # 4 seconds

  ras <- raster(rasVeg)
  ffFactor <- factor(ff)
  ras[!nas] <- ffFactor # 2 seconds

  eTable <- data.frame(ID = seq_along(levels(ffFactor)), VALUE = levels(ffFactor))
  types <- strsplit(as.character(eTable$VALUE), split = splitVal)
  types <- do.call(rbind, types)

  ## ensure species names all consistent (TODO: ensure this propagates)
  whMixed <- which(types[, 2] == "Mixed")
  types[, 2] <- equivalentName(types[, 2], sppEquiv, sppEquivCol)
  types[whMixed, 2] <- "Mixed"

  levels(ras) <- data.frame(eTable, ageClass = types[, 1], vegCover = types[, 2])

  # prepare poly factor raster
  if (is(poly, "SpatialPolygons")) {
    if (!"shinyLabel" %in% colnames(poly@data))
      stop("poly must have a column 'shinyLabel'")

    poly <- Cache(fasterize2, rasTsf, poly, field = "polygonNum")
  }
  levs <- raster::levels(poly)[[1]]

  # this is same, if all values present: e.g., 1, 2, 3, 4, 5 ...,
  # but not if missing: e.g., 1, 2, 3, 5
  levs <- factorValues(poly, levs$ID)
  facVals <- factorValues(
    poly,
    poly[],
    att = c("shinyLabel", "polygonNum")
  )

  bb <- data.table(
    zone = facVals$shinyLabel,
    polygonID = facVals$polygonNum,
    cell = seq_len(ncell(ras))
  )

  # add age and vegCover by reference
  bb[, c("ageClass", "vegCover") := factorValues(ras, ras[][bb$cell], att = c("ageClass", "vegCover"))]
  bb <- na.omit(bb)

  # One species at a time -- collapse polygons with same 'zone' name
  #tabulated <- bb[, list(NPixels = .N), by = c("zone", "polygonID", "ageClass", "vegCover")] ## keeps polyID
  tabulated <- bb[, list(NPixels = .N), by = c("zone", "ageClass", "vegCover")] ## dedupes the zones
  tabulated[, proportion := round(NPixels / sum(NPixels), 4), by = c("zone", "vegCover")]

  # All species -- collapse polygons with same 'zone' name
  #tabulated2 <- bb[, list(NPixels = .N), by = c("zone", "polygonID", "ageClass")] ## keeps polyID
  tabulated2 <- bb[, list(NPixels = .N), by = c("zone", "ageClass")] ## dedupes the zones
  tabulated2[, proportion := round(NPixels / sum(NPixels), 4), by = c("zone")]
  set(tabulated2, NULL, "vegCover", "All species")

  tabulated <- rbindlist(list(tabulated, tabulated2), use.names = TRUE, fill = TRUE)

  ## column containing the factor names varies, so we need to search for the right one
  colID <- which(colnames(raster::levels(rasVeg)[[1]]) %in% c("category", "Factor", "VALUE"))
  coverClasses <- raster::levels(rasVeg)[[1]][[colID]]
  if (is.factor(coverClasses))
    coverClasses <- levels(coverClasses)

  coverClasses <- as.character(coverClasses)
  emptyID <- which(coverClasses == "")
  if (length(emptyID))
    coverClasses <- coverClasses[-emptyID]

  if (!("All species" %in% levels(coverClasses)))
    coverClasses <- c(coverClasses, "All species")

  allCombos <- expand.grid(
    ageClass = ageClasses,
    vegCover = coverClasses,
    zone = levs$shinyLabel,
    stringsAsFactors = FALSE
  )
  #allCombos$polygonID <- match(allCombos$zone, levs$shinyLabel)
  data.table::setDT(allCombos)

  tabulated <- merge(
    tabulated,
    allCombos,
    #by = c("zone", "vegCover", "ageClass", "polygonID"),
    by = c("zone", "vegCover", "ageClass"),
    all.y = TRUE
  )
  # fill in zeros where there is no value
  tabulated[is.na(proportion), proportion := 0]
  set(tabulated,
      NULL,
      "label",
      paste(
        tabulated$ageClass,
        paste(gsub(basename(dirname(tsf[1])), pattern = "\\.", replacement = ""),
              basename(tsf[1]), sep = "_"),
        sep = "."
      ))

  endTime <- Sys.time()
  message("    Leading cover calculation took ",
          format(endTime - startTime, digits = 2))

  return(tabulated)
}
