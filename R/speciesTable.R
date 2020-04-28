if (getRversion() >= "3.1.0") {
  utils::globalVariables(c("growthcurve", "longevity", "mortalityshape", "postfireregen",
                           "resproutage_max", "resproutage_min", "resproutprob",
                           "shadetolerance", "species"))
}

#' Customize species trait table values for LandWeb
#'
#' @param speciesTable The species traits table
#' @param runName      Optional. Character string indicating the LandWeb run name (scenario).
#' @param params       A named list (of parameters) of named lists (by species), with species
#'                     traits overrides (e.g., \code{list(seeddistance_eff = list(Abie_sp = 25))}).
#'
#' @export
updateSpeciesTable <- function(speciesTable, runName = NULL, params) {
  ## checks:
  traits <- names(params)
  missingTraits <- traits[!traits %in% names(speciesTable)]
  if (length(missingTraits))
    stop("The traits: ", paste(missingTraits, collapse = ", "),
         "\ndo not exist in `speciesTable`")

  for (trt in traits) {
    subParams <- params[[trt]]

    ## checks:
    spp <- names(subParams)
    missingSpp <- spp[!spp %in% speciesTable$species]
    if (length(missingSpp))
      stop("The species: ", paste(missingSpp, collapse = ", "),
           "\ndo not exist in `speciesTable$species`")

    ## this is sub ideal to convert classes, but the only solution I found.
    ## neither class(...) <- class(..) nor as(..., class(...)) were changing integer to numeric.
    ## only as.numeric(...) worked
    if (class(unlist(subParams)) != class(speciesTable[[trt]])) {
      asFUN <- get(paste0("as.", class(unlist(subParams))))
      if (!is.function(asFUN))
        stop("Can't find an `as.` method to convert class to '", class(unlist(subParams)), "'")

      speciesTable <- speciesTable[, (trt) := lapply(.SD, asFUN), .SDcols = trt]
    }

    for (sp in spp) {
      set(speciesTable, which(speciesTable$species == sp), trt, subParams[[sp]])
    }
  }

  if (!is.null(runName)) {
    if (grepl("noDispersal|aspenDispersal|highDispersal", runName)) {
      speciesTable[, postfireregen := "resprout"] ## force all species to resprout
      speciesTable[, resproutprob := 1.0]     # default 0.5
      speciesTable[, resproutage_min := 0]    # defaults vary by species
      speciesTable[, resproutage_max := 400]  # defaults vary by species

      # decrease shade tolerance
      speciesTable[species == "Abie_sp", shadetolerance := 3]     # original default 4
      speciesTable[species == "Pice_gla", shadetolerance := 2]    # original default 3
      speciesTable[species == "Pice_mar", shadetolerance := 3]    # original default 4
      speciesTable[species == "Pinu_ban", shadetolerance := 1]    # original default 1
      speciesTable[species == "Pinu_con", shadetolerance := 1]    # original default 1
      speciesTable[species == "Pinu_sp", shadetolerance := 1]     # original default 1
      speciesTable[species == "Popu_sp", shadetolerance := 1]     # original default 1
    }

    if (grepl("aspen80", runName)) {
      speciesTable[species == "Popu_sp", longevity := 80] # default 150
    }

    ## TODO: THESE SHOULD BE PULLED OUT OF THE FUNCTION AND PASSED VIA
    ## THE params ARGUMENT.
    ## growth curves:
    #   Biomass Succession User Guide p17, 0 is faster growth, 1 was the prev assumption
    speciesTable[species == "Abie_sp", growthcurve := 0]  # original default 0
    speciesTable[species == "Pice_gla", growthcurve := 1] # original default 1
    speciesTable[species == "Pice_mar", growthcurve := 1] # original default 1
    speciesTable[species == "Pinu_ban", growthcurve := 0] # original default 0
    speciesTable[species == "Pinu_con", growthcurve := 0] # original default 0
    speciesTable[species == "Pinu_sp", growthcurve := 0]  # original default 0
    speciesTable[species == "Popu_sp", growthcurve := 0]  # original default 0

    ## mortality
    speciesTable[species == "Abie_sp", mortalityshape := 15]  # default 15
    speciesTable[species == "Pice_gla", mortalityshape := 15] # default 15
    speciesTable[species == "Pice_mar", mortalityshape := 15] # default 15
    speciesTable[species == "Pinu_ban", mortalityshape := 15] # default 15
    speciesTable[species == "Pinu_con", mortalityshape := 15] # default 15
    speciesTable[species == "Pinu_sp", mortalityshape := 15]  # default 15
    speciesTable[species == "Popu_sp", mortalityshape := 25]  # default 25

    ## resprouting (normally, only aspen resprouts)
    speciesTable[species == "Popu_sp", resproutage_min := 25] # default 10
    #speciesTable[species == "Popu_sp", resproutprob := 0.1]  # default 0.5
  }

  return(speciesTable)
}
