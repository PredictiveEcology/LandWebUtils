if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(":=", ".N", ".SD", "ageClass", "group", "N", "NCC",
                           "polygonName", "vegCover"))
}

#' @importFrom graphics abline axis barplot hist
#' @importFrom pemisc factorValues2
#' @keywords internal
.doPlotHistogram <- function(data, colName, colNameCC, xlim, force.min.n = FALSE, fname = NULL, ...) {
  minNumBars <- 6
  maxNumBars <- 30
  rangeNClusters <- if (isTRUE(force.min.n)) {
    range(c(0, xlim, minNumBars)) ## TODO: verify for largePatches
  } else {
    xlim
  }
  attemptedNumBars <- max(minNumBars, min(maxNumBars, diff(rangeNClusters)))
  breaksRaw <- seq(rangeNClusters[1], rangeNClusters[2], length.out = attemptedNumBars)
  prettyBreaks <- if (isTRUE(force.min.n)) {
    pretty(breaksRaw, n = attemptedNumBars, min.n = min(attemptedNumBars, minNumBars))
  } else {
    pretty(breaksRaw, n = attemptedNumBars)
  }

  dataForBreaks <- dataForHistogram <- if (NROW(data) == 0) {
    # add a bar at zero if there are no patches
    hist(0, plot = FALSE, breaks = prettyBreaks)
  } else {
    hist(data[[colName]], plot = FALSE, breaks = prettyBreaks)
  }

  breaksLabels <- dataForBreaks$breaks
  breaksInterval <- diff(breaksLabels)[1]

  histogramData <- dataForHistogram$counts / sum(dataForHistogram$counts) ## use proportion
  histogramData[is.na(histogramData)] <- 0 # NA means that there were no large patches in dt

  barplotBreaks <- seq_along(breaksLabels) - 0.5
  addAxisParams <- list(side = 1, labels = breaksLabels, at = barplotBreaks)# - min(breaksLabels))
  verticalLineAtX <- unique(data[[colNameCC]])[1] / breaksInterval + 0.5 # The barplot xaxis is 1/2 a barwidth off

  if (!is.null(fname)) png(fname, width = 800, height = 600, units = "px")
  barplot(histogramData, ...)
  if (!is.null(addAxisParams)) do.call(axis, addAxisParams)
  if (!is.null(verticalLineAtX)) abline(v = verticalLineAtX, col = "red", lwd = 3)
  if (!is.null(fname)) dev.off()
}

#' Generate histograms for large patches
#'
#' TODO: description needed
#'
#' @param map A \code{map} object.
#' @param functionName TODO: description needed
#' @param analysisGroups TODO: description needed
#' @param dPath Destination path for the resulting PNG files.
#'
#' @export
#' @importFrom data.table setnames
#' @importFrom magrittr %>%
#' @importFrom raster res
#' @importFrom reproducible checkPath
#' @importFrom tools toTitleCase
#' @importFrom utils write.csv
runHistsLargePatches <- function(map, functionName, analysisGroups, dPath) {
  dPath <- checkPath(dPath, create = TRUE)

  allRepPolys <- na.omit(map@metadata[[analysisGroups]])
  names(allRepPolys) <- allRepPolys

  lapply(allRepPolys, function(poly) {
    allData <- map@analysesData[[functionName]][["LargePatches"]][[poly]]
    if (is.null(allData))
      allData <- map@analysesData[[functionName]][[poly]] ## TODO: fix upstream

    data <- allData[!grepl("CC", group)]
    dataCC <- allData[grepl("CC", group)]

    data$rep <- as.numeric(factor(data$group)) ## TODO: workaround for incorrect rep values from LargePatches

    slices <- c("ageClass", "polygonName", "vegCover", "rep")
    slicesNoRep <- slices[slices != "rep"]

    data <- data[!ageClass == "NA" | !vegCover == "NA"]
    emptyDT <- data.table(expand.grid(ageClass = unique(data$ageClass),
                                      vegCover = unique(data$vegCover),
                                      polygonName = unique(data$polygonName),
                                      rep = unique(data$rep)))

    patchSizes <- c(100, 500, 1000, 5000) ## minPatchSize <- 100

    fout <- lapply(patchSizes, function(minPatchSize) {
      nClustersDT <- data[sizeInHa >= minPatchSize, c(N = .N, .(totalArea = sum(sizeInHa))), by = slices]
      nClustersDT <- nClustersDT[emptyDT, on = slices, nomatch = NA]
      nClustersDT[is.na(N), ':='(N = 0, totalArea = 0)]

      nClustersDT_CC <- dataCC[sizeInHa >= minPatchSize, c(NCC = .N, .(totalAreaCC = sum(sizeInHa))), by = slicesNoRep]
      nClustersDT_CC <- nClustersDT_CC[emptyDT, on = slicesNoRep, nomatch = NA]
      nClustersDT_CC[is.na(NCC), ':='(NCC = 0, totalAreaCC = 0.0)]

      nClustersDT <- nClustersDT[nClustersDT_CC, on = slices]

      try(write.csv(nClustersDT, file.path(dPath, paste0("largePatches_",
                                                         gsub(" ", "_", poly),
                                                         "_", minPatchSize, ".csv"))))

      saveDir <- checkPath(file.path(dPath, "largePatches", poly, minPatchSize), create = TRUE)
      savePng <- quote(file.path(saveDir, paste0(paste(unique(polygonName),
                                                       unique(vegCover),
                                                       unique(ageClass),
                                                       collapse = " "), ".png")))

      xlim <- c(0, max(nClustersDT[["N"]], nClustersDT[["NCC"]]))
      nClustersDT[, tryCatch(.doPlotHistogram(data = .SD,
                                              colName = "N",
                                              colNameCC = "NCC",
                                              fname = eval(savePng),
                                              # ccLine = NCC,
                                              border = "grey",
                                              col = "darkgrey",
                                              main = paste(unique(polygonName),
                                                           unique(vegCover),
                                                           unique(ageClass),
                                                           collapse = " "),
                                              space = 0,
                                              xlab = paste0("Number of patches greater than ",
                                                            minPatchSize, " ha"),
                                              ylab = "Proportion in NRV",
                                              xlim = xlim,
                                              ylim = c(0, 1),
                                              force.min.n = TRUE),
                             error = function(e) warning(e)),
                  .SDcols = c(slices, "N", "NCC"), by = slicesNoRep]

      nClustersDT[, list(filename = eval(savePng)), by = slicesNoRep]
    })
    unlist(fout)
  })
}

#' Generate histograms for leading vegetation cover
#'
#' TODO: description needed
#'
#' @param map A \code{map} object.
#' @param functionName TODO: description needed
#' @param analysisGroups TODO: description needed
#' @param dPath Destination path for the resulting PNG files.
#'
#' @export
#' @importFrom data.table setnames
#' @importFrom magrittr %>%
#' @importFrom raster res
#' @importFrom reproducible checkPath
#' @importFrom tools toTitleCase
#' @importFrom utils write.csv
#' @include misc.R
runHistsVegCover <- function(map, functionName, analysisGroups, dPath) {
  allRepPolys <- na.omit(map@metadata[[analysisGroups]])
  names(allRepPolys) <- allRepPolys

  lapply(allRepPolys, function(poly) {
    allData <- map@analysesData[[functionName]][["LeadingVegTypeByAgeClass"]][[poly]]
    if (is.null(allData))
      allData <- map@analysesData[[functionName]][[poly]] ## TODO: fix upstream
    allData <- unique(allData) ## remove duplicates; with LandWeb#89

    ## WORKAROUND inconsistent species names
    allData[["vegCover"]] <- gsub(" leading", "", allData[["vegCover"]]) %>% ## no longer needed?
      tools::toTitleCase() %>%
      as.factor() ## match CC raster names
    allData[vegCover == "Fir", vegCover := "Abie_sp"] ## so far, rep01 is the only one needing fixing
    allData[vegCover == "Wh Spruce", vegCover := "Pice_gla"]
    allData[vegCover == "Bl Spruce", vegCover := "Pice_mar"]
    allData[vegCover == "Pine", vegCover := "Pinu_sp"]
    allData[vegCover == "Decid", vegCover := "Popu_sp"]

    allData$ageClass <- factor(allData$ageClass, .ageClasses)

    data <- allData[!grepl("CC", group)]
    data$rep <- as.numeric(factor(data$group))

    dataCC <- allData[grepl("CC", group)]
    setnames(dataCC, "proportion", "proportionCC") ## rename the column to proportionCC
    dataCC <- dataCC[, c("group", "label", "NPixels") := list(NULL, NULL, NULL)]

    data2 <- dataCC[data, on = .(zone, vegCover, ageClass)]
    data2[is.na(NPixels), NPixels := 0]

    ## sum = all species + each indiv species = 2 * totalPixels
    ## NOTE: this is number of TREED pixels, which is likely smaller than the polygon area
    data2[, totalPixels := as.double(base::sum(NPixels, na.rm = TRUE)),
          by = c("group", "vegCover", "zone")]
    data2[, totalPixels2 := as.double(base::mean(totalPixels, na.rm = TRUE)),
          by = c("vegCover", "zone")] ## use mean for plot labels below
    #try(write.csv(data2, file.path(dPath, paste0("leading_", gsub(" ", "_", poly), ".csv"))))

    saveDir <- checkPath(file.path(dPath, "leading", poly), create = TRUE)
    savePng <- quote(file.path(saveDir, paste0(paste(unique(zone),
                                                     unique(vegCover),
                                                     unique(ageClass),
                                                     collapse = " "), ".png")))
    slices <- c("zone", "vegCover", "ageClass")

    data2[, tryCatch(.doPlotHistogram(data = .SD,
                                      colName = "proportion",
                                      colNameCC = "proportionCC",
                                      fname = eval(savePng),
                                      border = "grey",
                                      col = "darkgrey",
                                      main = paste(unique(zone),
                                                   unique(vegCover),
                                                   unique(ageClass),
                                                   collapse = "_"),
                                      space = 0,
                                      xlab = paste0("Proportion of forest area (total ",
                                                    format(unique(totalPixels2) *
                                                             prod(res(rasterToMatch(map))) / 1e4,
                                                           big.mark = ","),
                                                    " ha)"),
                                      ylab = "Proportion in NRV",
                                      xlim = c(0, 1),
                                      ylim = c(0, 1),
                                      force.min.n = FALSE),
                     error = function(e) warning(e)),
          .SDcols = c(slices, "rep", "proportion", "proportionCC"), by = slices]
    data2[, list(filename = eval(savePng)), by = slices]
  })
}
