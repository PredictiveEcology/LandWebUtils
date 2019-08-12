if (getRversion() >= "3.1.0") {
  utils::globalVariables(c(".", ":=", ".SD", "ageClass", "group", "proportionCC",
                           "totalPixels", "totalPixels2", "vegCover", "zone"))
}

#' @export
#' @importFrom data.table data.table
#' @importFrom graphics boxplot points
#' @importFrom grDevices dev.off png
#' @importFrom utils write.table
.doPlotBoxplot <- function(data, authStatus, fname = NULL, ageClasses,
                           fout = NULL, vegCover, zone, ...) {
  if (!is.null(fname)) png(fname, height = 600, width = 800, units = "px")
  a <- boxplot(proportion ~ as.factor(ageClass), data, ...)

  boxplotData <- data.table(zone = rep(zone, 4),
                            vegCover = rep(vegCover, 4),
                            ageClass = ageClasses,
                            proportionCC = data$proportionCC[1:4],
                            MIN = a$stats[1, ],
                            Q1 = a$stats[2, ],
                            MED = a$stats[3, ],
                            Q3 = a$stats[4, ],
                            MAX = a$stats[5, ])

  if (!is.null(fout))
    try(write.table(boxplotData, file = fout, append = TRUE, col.names = FALSE,
                    row.names = FALSE, sep = ","))

  if (isTRUE(authStatus)) {
    ids <- match(factor(data$ageClass[1:4]), data[1:4, ]$ageClass)
    points(data$proportionCC[ids], factor(data$ageClass[1:4]), col = "red", pch = 20, cex = 3)
  }
  if (!is.null(fname)) dev.off()
}

#' Generate box and whisker plots for leading vegetation cover
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
#' @importFrom map rasterToMatch
#' @importFrom raster res
#' @importFrom reproducible checkPath
#' @importFrom tools toTitleCase
#' @importFrom utils write.csv
runBoxPlotsVegCover <- function(map, functionName, analysisGroups, dPath) {
  ageClasses <- c("Young", "Immature", "Mature", "Old")
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

    allData$ageClass <- factor(allData$ageClass, ageClasses)

    data <- allData[!grepl("CC", group)]

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
    try(write.csv(data2, file.path(dPath, paste0("leading_", gsub(" ", "_", poly), ".csv"))))

    ## output the box and whisker plot ranges (quartiles, etc.)
    empty <- data.table(zone = character(0),  vegCover = character(0),
                        ageClass = character(0), proportionCC = numeric(0),
                        MIN = numeric(0), Q1 = numeric(0), MED = numeric(0),
                        Q3 = numeric(0), MAX = numeric(0))
    fout <- file.path(dPath, paste0("leading_boxplots_", gsub(" ", "_", poly), ".csv"))
    try(write.csv(empty, fout, row.names = FALSE))

    saveDir <- checkPath(file.path(dPath, poly), create = TRUE)
    savePng <- quote(file.path(saveDir, paste0(unique(paste(zone, vegCover, collapse = " ")), ".png")))
    slices <- c("zone", "vegCover")
    data2[, tryCatch(.doPlotBoxplot(data = .SD,
                                    authStatus = TRUE,
                                    col = "limegreen",
                                    fname = eval(savePng),
                                    ageClasses = ageClasses,
                                    fout = fout,
                                    vegCover = vegCover,
                                    zone = zone,
                                    horizontal = TRUE,
                                    main = unique(paste(zone, vegCover, collapse = "_")),
                                    xlab = paste0("Proportion of forest area (total ",
                                                  format(unique(totalPixels2) *
                                                           prod(res(rasterToMatch(map))) / 1e4,
                                                         big.mark = ","),
                                                  " ha)"),
                                    ylab = "Age class",
                                    ylim = c(0, 1)),
                     error = function(e) warning(e)),
          .SDcols = c("ageClass", "proportion", "proportionCC", "NPixels"), by = slices]
    data2[, list(filename = eval(savePng)), by = slices]
  })
}
