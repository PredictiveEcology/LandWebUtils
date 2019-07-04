if (getRversion() >= "3.1.0") {
  utils::globalVariables(c("growthcurve", "longevity", "mortalityshape", "postfireregen",
                           "resproutage_min", "resproutprob", "shadetolerance", "species"))
}

#' Customize species trait table values for LandWeb
#'
#' @param speciesTable The species traits table
#' @param runName      Character string indicating the LandWeb run name (scenario)
#' @param params       A named list (of parameters) of named lists (by species), providing species
#'                     table overrides (e.g., \code{list(seeddistance_eff = list(Abie_sp = 25))}).
#'
#' @export
updateSpeciesTable <- function(speciesTable, runName, params) {
  speciesTable[species == "Abie_sp", `:=`(seeddistance_eff = params$seeddistance_eff$Abie_sp,
                                          seeddistance_max = params$seeddistance_max$Abie_sp)]
  speciesTable[species == "Pice_gla", `:=`(seeddistance_eff = params$seeddistance_eff$Pice_gla,
                                           seeddistance_max = params$seeddistance_max$Pice_gla)]
  speciesTable[species == "Pice_mar", `:=`(seeddistance_eff = params$seeddistance_eff$Pice_mar,
                                           seeddistance_max = params$seeddistance_max$Pice_mar)]
  speciesTable[species == "Pinu_ban", `:=`(seeddistance_eff = params$seeddistance_eff$Pinu_ban,
                                           seeddistance_max = params$seeddistance_max$Pinu_ban)]
  speciesTable[species == "Pinu_con", `:=`(seeddistance_eff = params$seeddistance_eff$Pinu_con,
                                           seeddistance_max = params$seeddistance_max$Pinu_con)]
  speciesTable[species == "Pinu_sp", `:=`(seeddistance_eff = params$seeddistance_eff$Pinu_sp,
                                          seeddistance_max = params$seeddistance_max$Pinu_sp)]
  speciesTable[species == "Popu_sp", `:=`(seeddistance_eff = params$seeddistance_eff$Popu_sp,
                                          seeddistance_max = params$seeddistance_max$Popu_sp)]

  ## resprouting (normally, only aspen resprouts)
  speciesTable[species == "Popu_sp", resproutage_min := 25] # default 10
  #speciesTable[species == "Popu_sp", resproutprob := 0.1]  # default 0.5
  if (grepl("noDispersal|aspenDispersal|highDispersal", runName)) {
    speciesTable[, postfireregen := "resprout"] ## force all species to resprout
    speciesTable[, resproutprob := 1.0]     # default 0.5
    speciesTable[, resproutage_min := 0]    # defaults vary by species
    speciesTable[, resproutage_max := 400]  # defaults vary by species
    #speciesTable[, shadetolerance := 5]     # defaults vary by species
  }

  ## growth curves:
  #   Biomass Succession User Guide p17, 0 is faster growth, 1 was the prev assumption
  speciesTable[species == "Abie_sp", growthcurve := 1]  # original default 0
  speciesTable[species == "Pice_gla", growthcurve := 1] # original default 1
  speciesTable[species == "Pice_mar", growthcurve := 1] # original default 1
  speciesTable[species == "Pinu_ban", growthcurve := 1] # original default 0
  speciesTable[species == "Pinu_con", growthcurve := 1] # original default 0
  speciesTable[species == "Pinu_sp", growthcurve := 1]  # original default 0
  speciesTable[species == "Popu_sp", growthcurve := 1]  # original default 0

  ## mortality
  speciesTable[species == "Abie_sp", mortalityshape := 15]  # default 15
  speciesTable[species == "Pice_gla", mortalityshape := 15] # default 15
  speciesTable[species == "Pice_mar", mortalityshape := 15] # default 15
  speciesTable[species == "Pinu_ban", mortalityshape := 15] # default 15
  speciesTable[species == "Pinu_con", mortalityshape := 15] # default 15
  speciesTable[species == "Pinu_sp", mortalityshape := 15]  # default 15
  speciesTable[species == "Popu_sp", mortalityshape := 25]  # default 25

  if (grepl("aspen80", runName)) {
    speciesTable[species == "Popu_sp", longevity := 80] # default 150
  }

  return(speciesTable)
}
