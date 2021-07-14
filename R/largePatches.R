utils::globalVariables(c(":=", "sizeInHa"))

#' Calculate proportion of large patches in NRV
#'
#' TODO: needs description
#'
#' @param tsf TODO: description needed
#' @param vtm TODO: description needed
#' @param poly TODO: description needed
#' @param labelColumn TODO: description needed
#' @param id TODO: description needed
#' @param ageClassCutOffs TODO: description needed
#' @param ageClasses TODO: description needed
#' @param sppEquivCol TODO: description needed
#' @param sppEquiv TODO: description needed
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom LandR equivalentName
#' @importFrom map areaAndPolyValue fasterize2 .rasterToMemory
#' @importFrom raster levels raster reclassify
#' @importFrom reproducible Cache
LargePatches <- function(tsf, vtm, poly, labelColumn, id, ageClassCutOffs, ageClasses,
                         sppEquivCol, sppEquiv) {
  vtm <- vtm[1]

  if (basename(vtm) == "CurrentConditionVTM.tif") ## TODO: LandWeb workaround
    tsf <- file.path(dirname(vtm), "CurrentConditionTSF.tif")

  timeSinceFireFilesRast <- Cache(.rasterToMemory, tsf[1])

  tsf <- reclassify(timeSinceFireFilesRast,
                    cbind(from = ageClassCutOffs - 0.1,
                          to = c(ageClassCutOffs[-1], Inf),
                          seq_along(ageClasses)))
  levels(tsf) <- data.frame(ID = seq_along(ageClasses), Factor = ageClasses)

  poly$tmp <- factor(poly[[labelColumn]])
  rasRepPoly <- Cache(
    fasterize2,
    poly,
    emptyRaster = raster(timeSinceFireFilesRast), # doesn't need to the data -- makes Caching more effective
    field = "tmp"
  )

  # 3rd raster
  rasVeg <- Cache(.rasterToMemory, vtm)

  splitVal <- paste0("_", 75757575, "_") # unlikely to occur for any other reason

  # Individual species
  nas3 <- is.na(rasRepPoly[])
  nas2 <- is.na(rasVeg[]) | is.na(factorValues2(rasVeg, rasVeg[], att = 1))
  nas1 <- is.na(tsf[])
  nas <- nas3 | nas2 | nas1

  if (!isTRUE(all(nas))) {
    #name1a <- as.character(raster::levels(tsf)[[1]]$Factor)[tsf[][!nas]]
    name1 <- as.character(factorValues2(tsf, tsf[], att = 2)[!nas])

    colID <- which(colnames(raster::levels(rasVeg)[[1]]) %in% c("category", "Factor", "VALUE"))

    #name2a <- as.character(raster::levels(rasVeg)[[1]][[colID]])[rasVeg[][!nas]]
    name2 <- as.character(factorValues2(rasVeg, rasVeg[], att = colID)[!nas])

    # rasRepPoly will have the numeric values of the *factor* in poly$tmp, NOT
    #   the raster::levels(rasRepPoly)[[1]])
    name3 <- raster::levels(poly$tmp)[rasRepPoly[][!nas]] ## fixed 2021-05-05

    if (!identical(length(name1), length(name2)) || !identical(length(name1), length(name3)))
      stop("There is something wrong with tsf or rasVeg or rasRepPoly inside LargePatches")

    ff <- paste(name1, name2, name3, sep = splitVal) # 4 seconds
    ras <- raster(rasVeg)
    ffFactor <- factor(ff)
    ras[!nas] <- ffFactor # 2 seconds ## note: sum(!nas, na.rm = TRUE) should equal length(ffFactor)

    areaAndPolyOut <- Cache(areaAndPolyValue, ras, length = Inf) # maybe lots of NAs on edge
    eTable <- data.frame(ID = seq_along(levels(ffFactor)), VALUE = levels(ffFactor))
    types <- strsplit(as.character(eTable$VALUE), split = splitVal)
    types <- do.call(rbind, types)

    facPolygonID <- factor(types[areaAndPolyOut$polyID, 3])
    outBySpecies <- data.table(polygonID = as.numeric(facPolygonID),
                               sizeInHa = areaAndPolyOut$sizeInHa,
                               vegCover = types[areaAndPolyOut$polyID, 2],
                               rep = id,
                               ageClass = types[areaAndPolyOut$polyID, 1],
                               polygonName = as.character(facPolygonID))

    # All species combined # remove name2
    ff <- paste(name1, name3, sep = splitVal)
    ff[grepl("NA", ff)] <- NA
    ras <- raster(rasVeg)
    ffFactor <- factor(ff)
    ras[!nas] <- ffFactor

    rm(areaAndPolyOut)
    areaAndPolyOut2 <- Cache(areaAndPolyValue, ras, length = Inf) # maybe lots of NAs on edge

    eTable <- data.frame(ID = seq_along(levels(ffFactor)), VALUE = levels(ffFactor))
    types <- strsplit(as.character(eTable$VALUE), split = splitVal)
    types <- do.call(rbind, types)

    facPolygonID <- factor(types[areaAndPolyOut2$polyID, 2])

    outAllSpecies <- data.table(polygonID = as.numeric(facPolygonID),
                                sizeInHa = areaAndPolyOut2$sizeInHa,
                                vegCover = "All species",
                                rep = id,
                                ageClass = types[areaAndPolyOut2$polyID, 1],
                                polygonName = as.character(facPolygonID))

    out <- rbindlist(list(outBySpecies, outAllSpecies))
    out <- out[sizeInHa >= 100] # never will need patches smaller than 100 ha
  } else {
    out <- data.table(polygonID = character(), sizeInHa = numeric(), vegCover = character(),
                      rep = numeric(), ageClass = numeric(), polygonName = numeric())
  }

  out[!is.na(equivalentName(out$vegCover, sppEquiv, sppEquivCol)),
      vegCover := equivalentName(vegCover, sppEquiv, sppEquivCol)]
  out
}
