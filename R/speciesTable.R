if (getRversion() >= "3.1.0") {
  utils::globalVariables(c("growthcurve", "longevity", "mortalityshape", "postfireregen",
                           "resproutage_min", "resproutprob", "shadetolerance", "species"))
}

#' Customize species trait table values for LandWeb
#'
#' @param speciesTable The species traits table
#' @param runName      The LandWeb run name
#'
#' @export
updateSpeciesTable <- function(speciesTable, runName) {
  if (grepl("aspenDispersal", runName)) {
    ## seed dispersal (see LandWeb#96, LandWeb#112)
    speciesTable[species == "Abie_sp", `:=`(seeddistance_eff = 0, seeddistance_max = 125)]     # defaults 25, 160
    speciesTable[species == "Pice_gla", `:=`(seeddistance_eff = 0, seeddistance_max = 125)]    # defaults 100, 303
    speciesTable[species == "Pice_mar", `:=`(seeddistance_eff = 0, seeddistance_max = 125)]    # defaults 80, 200
    speciesTable[species == "Pinu_ban", `:=`(seeddistance_eff = 0, seeddistance_max = 125)]    # defaults 30, 100
    speciesTable[species == "Pinu_con", `:=`(seeddistance_eff = 0, seeddistance_max = 125)]    # defaults 30, 100
    speciesTable[species == "Pinu_sp", `:=`(seeddistance_eff = 0, seeddistance_max = 125)]     # defaults 30, 100
    speciesTable[species == "Popu_sp", `:=`(seeddistance_eff = 100, seeddistance_max = 235)]   # defaults 200, 5000
  } else if (grepl("highDispersal", runName)) {
    ## seed dispersal (see LandWeb#96, LandWeb#112)
    speciesTable[species == "Abie_sp", `:=`(seeddistance_eff = 500, seeddistance_max = 1250)]  # defaults 25, 160
    speciesTable[species == "Pice_gla", `:=`(seeddistance_eff = 500, seeddistance_max = 1250)] # defaults 100, 303
    speciesTable[species == "Pice_mar", `:=`(seeddistance_eff = 500, seeddistance_max = 1250)] # defaults 80, 200
    speciesTable[species == "Pinu_ban", `:=`(seeddistance_eff = 500, seeddistance_max = 1250)] # defaults 30, 100
    speciesTable[species == "Pinu_con", `:=`(seeddistance_eff = 500, seeddistance_max = 1250)] # defaults 30, 100
    speciesTable[species == "Pinu_sp", `:=`(seeddistance_eff = 300, seeddistance_max = 3000)]  # defaults 30, 100
    speciesTable[species == "Popu_sp", `:=`(seeddistance_eff = 300, seeddistance_max = 3000)]  # defaults 200, 500
  }

  ## resprouting (normally, only aspen resprouts)
  if (grepl("noDispersal|aspenDispersal", runName)) {
    speciesTable[, postfireregen := "resprout"] ## force all species to resprout
    speciesTable[, resproutprob := 1.0]  # default 0.5
    speciesTable[, shadetolerance := 5]  # defaults vary by species
  }
  speciesTable[species == "Popu_sp", resproutage_min := 25] # default 10
  #speciesTable[species == "Popu_sp", resproutprob := 0.1]  # default 0.5

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
