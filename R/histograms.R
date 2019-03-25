if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(":=", ".N", ".SD", "ageClass", "group", "N", "NCC",
                           "polygonName", "vegCover"))
}

#' @importFrom graphics abline axis barplot hist
#' @importFrom pemisc factorValues2
.doPlotHistogram <- function(data, xlim, fname = NULL, ...) {
  minNumBars <- 6
  maxNumBars <- 30
  rangeNClusters <- range(c(0, xlim, minNumBars)) ## TODO: verify
  attemptedNumBars <- max(minNumBars, min(maxNumBars, diff(rangeNClusters)))

  breaksRaw <- seq(rangeNClusters[1], rangeNClusters[2], length.out = attemptedNumBars)
  prettyBreaks <- pretty(breaksRaw, n = attemptedNumBars, min.n = min(attemptedNumBars, minNumBars))

  dataForBreaks <- dataForHistogram <- if (NROW(data) == 0) {
    # add a bar at zero if there are no patches
    hist(0, plot = FALSE, breaks = prettyBreaks)
  } else {
    hist(data$N, plot = FALSE, breaks = prettyBreaks)
  }

  breaksLabels <- dataForBreaks$breaks
  breaksInterval <- diff(breaksLabels)[1]

  histogramData <- dataForHistogram$counts / sum(dataForHistogram$counts) ## use proportion
  histogramData[is.na(histogramData)] <- 0 # NA means that there were no large patches in dt
#browser()
  barplotBreaks <- seq_along(breaksLabels) - 0.5
  addAxisParams <- list(side = 1, labels = breaksLabels, at = barplotBreaks - min(breaksLabels))
  verticalLineAtX <- unique(data$NCC)[1] / breaksInterval + 0.5 # The barplot xaxis is 1/2 a barwidth off

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
#' @param dPath Destination path for the resulting png files.
#'
#' @export
#' @importFrom data.table setnames
#' @importFrom magrittr %>%
#' @importFrom raster res
#' @importFrom reproducible checkPath
#' @importFrom tools toTitleCase
#' @importFrom utils write.csv
runHistsLargePatches <- function(map, functionName, analysisGroups, dPath) {
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
      nClustersDT <- data[sizeInHa >= minPatchSize, .N, by = slices]
      nClustersDT <- nClustersDT[emptyDT, on = slices, nomatch = NA]
      nClustersDT[is.na(N), N := 0]

      nClustersDT_CC <- dataCC[sizeInHa >= minPatchSize, .N, by = slicesNoRep]
      #nClustersDT_CC$rep <- as.numeric(nClustersDT_CC$rep)
      setnames(nClustersDT_CC, "N", "NCC")
      nClustersDT_CC <- nClustersDT_CC[emptyDT, on = slicesNoRep, nomatch = NA]
      nClustersDT_CC[is.na(NCC), NCC := 0]

      nClustersDT <- nClustersDT[nClustersDT_CC, on = slices]

      try(write.csv(nClustersDT, file.path(dPath, paste0("largePatches_",
                                                         gsub(" ", "_", poly),
                                                         "_", minPatchSize, ".csv"))))

      saveDir <- checkPath(file.path(dPath, poly, "largePatches", minPatchSize), create = TRUE)
      savePng <- quote(file.path(saveDir, paste0(paste(unique(polygonName),
                                                       unique(vegCover),
                                                       unique(ageClass),
                                                       collapse = " "), ".png")))

      #browser()
      xlim <- c(0, max(nClustersDT$N, nClustersDT$NCC))
      nClustersDT[, tryCatch(.doPlotHistogram(data = .SD,
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
                                              ylim = c(0, 1)),
                             error = function(e) warning(e)),
                  .SDcols = c(slices, "N", "NCC"), by = slicesNoRep]

      nClustersDT[, list(filename = eval(savePng)), by = slicesNoRep]
    })
    unlist(fout)
  })
}
