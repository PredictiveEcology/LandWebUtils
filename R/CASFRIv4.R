utils::globalVariables(c(
  "AGE", "GID", "keepSpecies", "pct", "value"
))

#' Load CASFRI data for LandWeb
#'
#' TODO: description needed
#'
#' @param CASFRIRas TODO: description needed
#' @param attrFile TODO: description needed
#' @param headerFile TODO: description needed
#' @template sppEquiv
#' @template sppEquivCol
#' @param type Character string. Either `"cover"` or `"age"`.
#'
#' @return TODO: description needed
#'
#' @export
#' @importFrom data.table data.table fread melt set setkey
#' @importFrom reproducible asPath Cache
loadCASFRI <- function(CASFRIRas, attrFile, headerFile, sppEquiv, sppEquivCol,
                       type = c("cover", "age")) {
  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  sppNameVectorCASFRI <- equivalentName(sppNameVector, sppEquiv,  column = "CASFRI", multi = TRUE)

  # CASFRI stuff
  CASFRIheader <- fread(headerFile, skip = 14, nrows = 49, header = FALSE, sep = "", fill = TRUE)
  header <- apply(CASFRIheader, 1, function(x) sub(pattern = "(\t+| ).*$", "", x))
  CASFRIheader <- header[nchar(header) != 0]

  wantedColumns <- grep(CASFRIheader, pattern = "^SPECIES|^GID|^AGE")

  CASFRIattr <- fread(asPath(attrFile), select = wantedColumns)

  setnames(CASFRIattr, CASFRIheader[wantedColumns])
  setkey(CASFRIattr, "GID")

  NAVals <- c("XXXX MISS", "UNDEF", "XXXX ERRC")
  numSpeciesColumns <- length(grep("SPECIES_PER", names(CASFRIattr), value = TRUE))
  if (type[1] == "cover") {
    for (i in seq(numSpeciesColumns)) {
      set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_", i)]] %in% NAVals),
          paste0("SPECIES_", i), NA_character_)
      set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_PER_", i)]] %in% NAVals),
          paste0("SPECIES_", i), NA_character_)
    }
    for (i in 1:1) {
      message("remove CASFRI entries with <15 cover as dominant species,",
              " i.e., these pixels are deemed untreed")
      CASFRIattr <- CASFRIattr[which(CASFRIattr[[paste0("SPECIES_PER_", i)]] > 15), ]
    }
    message("set CASFRI entries with <15 cover in 2nd-5th dominance class to NA")
    for (i in 2:5) {
      set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_PER_", i)]] <= 15),
          paste0("SPECIES_", i), NA_character_)
    }

    CASFRIattrLong <- melt(CASFRIattr, id.vars = c("GID"),
                           measure.vars = paste0("SPECIES_", 1:5))
    CA2 <- melt(CASFRIattr, id.vars = c("GID"),
                measure.vars = c(paste0("SPECIES_PER_", 1:5)))
    CASFRIattrLong[, pct := CA2$value]
    rm(CA2)
    CASFRIattrLong <- na.omit(CASFRIattrLong)
    CASFRIattrLong <- CASFRIattrLong[value %in% sppNameVectorCASFRI]
  } else {
    CASFRIattrLong <- CASFRIattr[, .(GID, AGE)]
    CASFRIattrLong <- CASFRIattrLong[!is.na(AGE) & AGE > -1]
  }

  CASFRIdt <- data.table(GID = CASFRIRas[], rastInd = 1:ncell(CASFRIRas))
  CASFRIdt <- CASFRIdt[!is.na(GID)]
  #CASFRIdt <- CASFRIdt[isNA == FALSE]
  setkey(CASFRIdt, GID)
  #set(CASFRIdt, NULL, "isNA", NULL)

  return(list(CASFRIattrLong = CASFRIattrLong, CASFRIdt = CASFRIdt))
}

#' `CASFRItoSpRasts`
#'
#' TODO: description and title needed
#'
#' @param CASFRIRas TODO: description needed
#' @param CASFRIattrLong TODO: description needed
#' @param CASFRIdt TODO: description needed
#' @template sppEquiv
#' @template sppEquivCol
#' @template destinationPath
#'
#' @return TODO: description needed
#'
#' @export
#' @importFrom data.table setkey
#' @importFrom magrittr %>%
#' @importFrom reproducible asPath Cache
#' @importFrom raster crs crs<- NAvalue<- raster setValues stack writeRaster
CASFRItoSpRasts <- function(CASFRIRas, CASFRIattrLong, CASFRIdt,
                            sppEquiv, sppEquivCol, destinationPath) {
  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]), ]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  # This
  sppListMergesCASFRI <- lapply(sppNameVector, function(x)
    equivalentName(x, sppEquiv,  column = "CASFRI", multi = TRUE)
  )

  ## create list and template raster
  spRasts <- list()
  spRas <- raster(CASFRIRas) %>% setValues(., NA_integer_)

  ## NOT SURE IF THESE LINES ABOUT NA are relevant -- Eliot Dec 7
  ## selected spp absent from CASFRI data
  NA_Sp <- which(is.na(sppListMergesCASFRI))#setdiff(speciesLandR, unique(keepSpecies$spGroup))

  ## All NA_Sp species codes should be in CASFRI spp list
  if (length(NA_Sp))
    warning("Not all selected species are in loadedCASFRI. Check if this is correct:\n",
            paste(paste0(keepSpecies$CASFRI[NA_Sp], collapse = ", "), "absent\n"))

  ## empty rasters for NA_sp
  for (sp in NA_Sp) {
    message("  running ", sp, ". Assigning NA, because absent from CASFRI")
    spRasts[[sp]] <- spRas
    NAval <- 65535L
    spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                           filename = asPath(file.path(destinationPath,
                                                       paste0("CASFRI_", sp, ".tif"))),
                           overwrite = TRUE, datatype = "INT2U", NAflag = NAval)
    ## NAvals need to be converted back to NAs
    NAvalue(spRasts[[sp]]) <- NAval
  }

  sppTODO <- unique(names(sppListMergesCASFRI))

  for (sp in sppTODO) {
    spCASFRI <- sppListMergesCASFRI[[sp]]
    spRasts[[sp]] <- spRas
    message("starting ", sp)
    if (length(spCASFRI) > 1)
      message("  Merging ", paste(spCASFRI, collapse = ", "), "; becoming: ", sp)
    aa2 <- CASFRIattrLong[value %in% spCASFRI][, min(100L, sum(pct)), by = GID]
    setkey(aa2, GID)
    cc <- aa2[CASFRIdt] %>% na.omit()
    rm(aa2)
    spRasts[[sp]][cc$rastInd] <- cc$V1
    message("  ", sp, " writing to disk")

    startCRS <- crs(spRasts[[sp]])
    NAval <- 255L
    spRasts[[sp]] <- writeRaster(spRasts[[sp]],
                                 filename = asPath(file.path(destinationPath,
                                                             paste0("CASFRI_", sp, ".tif"))),
                                 datatype = "INT1U", overwrite = TRUE, NAflag = NAval)
    ## NAvals need to be converted back to NAs
    NAvalue(spRasts[[sp]]) <- NAval

    if (is(spRasts[[sp]], "Raster")) {
      # Rasters need to have their disk-backed value assigned, but not shapefiles
      # This is a bug in writeRaster was spotted with crs of rastTmp became
      # +proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
      # should have stayed at
      # +proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0
      if (!identical(startCRS, crs(spRasts[[sp]])))
        crs(spRasts[[sp]]) <- startCRS
    }
    message("  ", sp, " done")
  }

  raster::stack(spRasts)
}

#' Prepare species layers from CASFRI v4
#'
#' @inheritParams LandR::prepSpeciesLayers_KNN
#'
#' @export
#' @importFrom reproducible asPath Cache prepInputs
prepSpeciesLayers_CASFRI <- function(destinationPath, outputPath,
                                     url = NULL,
                                     studyArea, rasterToMatch,
                                     sppEquiv,
                                     sppEquivCol, ...) {
  if (is.null(url))
    url <- "https://drive.google.com/file/d/1y0ofr2H0c_IEMIpx19xf3_VTBheY0C9h"

  CASFRItiffFile <- asPath(file.path(destinationPath, "Landweb_CASFRI_GIDs.tif"))
  CASFRIattrFile <- asPath(file.path(destinationPath, "Landweb_CASFRI_GIDs_attributes3.csv"))
  CASFRIheaderFile <- asPath(file.path(destinationPath, "Landweb_CASFRI_GIDs_README.txt"))

  message("  Loading CASFRI layers...")
  CASFRIRas <- Cache(prepInputs,
                     #targetFile = asPath("Landweb_CASFRI_GIDs.tif"),
                     targetFile = basename(CASFRItiffFile),
                     archive = asPath("CASFRI for Landweb.zip"),
                     url = url,
                     alsoExtract = c(CASFRItiffFile, CASFRIattrFile, CASFRIheaderFile),
                     destinationPath = destinationPath,
                     fun = "raster::raster",
                     studyArea = studyArea,
                     rasterToMatch = rasterToMatch,
                     method = "bilinear", ## ignore warning re: ngb (#5)
                     datatype = "INT4U",
                     filename2 = NULL,
                     overwrite = TRUE,
                     userTags =  c("CASFRIRas", "stable"))

  message("Load CASFRI data and headers, and convert to long format, and define species groups")

  loadedCASFRI <- Cache(loadCASFRI,
                        CASFRIRas = CASFRIRas,
                        attrFile = CASFRIattrFile,
                        headerFile = CASFRIheaderFile, ## TODO: this isn't used internally
                        sppEquiv = sppEquiv,
                        sppEquivCol = sppEquivCol,
                        type = "cover")

  message("Make stack from CASFRI data and headers")
  CASFRISpStack <- CASFRItoSpRasts(CASFRIRas = CASFRIRas,
                                   sppEquiv = sppEquiv,
                                   sppEquivCol = sppEquivCol,
                                   CASFRIattrLong = loadedCASFRI$CASFRIattrLong,
                                   CASFRIdt = loadedCASFRI$CASFRIdt,
                                   destinationPath = outputPath)

  return(CASFRISpStack)
}
