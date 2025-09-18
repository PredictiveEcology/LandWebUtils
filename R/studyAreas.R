utils::globalVariables(c("FMU_NAME", "Name"))

#' Study areas available for LandWeb simulations
#'
#' @return `LandWebStudyAreas` is a data.frame` with 4 columns:
#' - `Name` specifies the shortname of the study area for use with simulation setup;
#' - `Province` identifies the location of the study area;
#' - `ID` corresponds to the polygon ID of the Updated FMA Boundaries polygons;
#' - `Description` provides the full study area name and additional context (if any);
#'
#' @export
LandWebStudyAreas <- tibble::tribble(
  ~Name,              ~Province,  ~ID,                 ~Description,

  ## by company / FMA
  "AlPac",            "AB",       "38",                "ALPAC Forest Products Incorporated",
  "ANC",              "AB",       "39",                "ANC Timber Ltd.",
  "BlueRidge",        "AB",       "40",                "Blueridge Lumber Inc.",
  "DMI",              "AB",       "42, 43",            "Mercer Peace River Pulp Ltd. (formerly DMI)",
  "Edson",            "AB",       "53",                "West Fraser Mills Ltd. (Edson)",
  "FMANWT",           "NWT",      "59",                "Fort Resolution",
  "FMANWT2",          "NWT",      "58",                "Fort Providence",
  "LP_BC",            "BC",       "5, 6",              "Lousiana Pacific (British Columbia)",
  "LP_MB",            "MB",       "32",                "Lousiana Pacific (Manitoba)",
  "Manning",          "AB",       "45",                "Manning Diversified Forest Products Ltd.",
  "MillarWestern",    "AB",       "46",                "Millar Western Forest Products Ltd.",
  "Mistik",           "SK",       "15",                "Mistik",
  "MPR",              "AB",       "42, 43",            "Mercer Peace River Pulp Ltd. (formerly DMI)",
  "SprayLake",        "AB",       "47",                "Spray Lake Sawmills (1980) Ltd.",
  "Sundre",           "AB",       "48",                "Sundre Forest Products Inc.",
  "Tolko_AB_N",       "AB",       "50",                "Tolko Industries Ltd. (Alberta North)",
  "Tolko_AB_S",       "AB",       "44, 49, 51",        "Tolko Industries Ltd. (Alberta South)",
  "Tolko_SK",         "SK",       "22",                "Tolko (Saskatchewan)",
  "Vanderwell",       "AB",       "51, 52",            "Vanderwell Contractors (1971) Ltd.",
  "WestFraser_N",     "AB",       "44, 51, 55",        "West Fraser Mills Ltd. (Slave Lake)",
  "WestFraser_S",     "AB",       "53, 54",            "West Fraser Mills Ltd. (Edson + Hinton)",
  "WeyCo_GP",         "AB",       "56",                "Weyerhauser Company Ltd. (Grand Prairie)",
  "WeyCo_PT",         "AB",       "57",                "Weyerhauser Company Ltd. (Pembina Timberland)",
  "WeyCo_SK",         "SK",       "21",                "Weyerhauser Company Ltd. (Pasquia-Porcupine)",

  ## Alberta subregions / by FMU
  "FMU",              "AB",       NA_character_,       "Alberta Forest Management Units",
  "NWAB",             "AB",       "42, 43, 45, 50",    "Northwestern Alberta",

  ## Full LandWeb study area
  "LandWeb",          NA_character_,  NA_character_,   "full LandWeb study area",

  ## by province
  "provAB",           "AB",       "38--57",            "Alberta",
  "provMB",           "MB",       "31, 32, 33",        "Manitoba",
  "provNWT",          "NWT",      "58, 59",            "Northwest Territories",
  "provSK",           "SK",       "7--30",             "Saskatchewan",

  ## small test area
  "random",           "AB",       NA_character_,  "Small random study area used for testing"
) |>
  data.frame()

#' @param name character specifying the study area name. Must be one of `LandWebStudyAreas$Name`.
#'
#' @return `studyAreaProv` returns a character indicating the province the study area is in.
#'
#' @export
#' @rdname LandWebStudyAreas
studyAreaProv <- function(name) {
  ## use startsWith() to allow name with a runName suffix (e.g., '_v3')
  dplyr::filter(LandWebStudyAreas, startsWith(name, Name))[["Province"]]
}

#' @param prov character specifying the province abbreviation (i.e., "AB", "BC", "MB", "NWT", or "SK").
#'
#' @return `studyAreaIn` returns a logical indicating whether study area `name` is in `prov`.
#'
#' @export
#' @rdname LandWebStudyAreas
studyAreaIn <- function(name, prov) {
  Province <- studyAreaProv(name)
  is.na(Province) || (Province == prov)
}

#' Extract boundary polygon(s) for LandWeb forest management area(s)
#'
#' @param fmas  `sf` object with LandWeb FMA boundaries
#'
#' @param name  A character (regex) string to match.
#'
#' @return `sf` polygons object
#'
#' @export
extractFMA <- function(fmas, name) {
  fmas[grepl(name, fmas[["Name"]]), ]
}

#' Extract boundary polygon(s) for Alberta forest management unit(s)
#'
#' @param fmus  `sf` object with Alberta FMU boundaries
#'
#' @param name  A character (regex) string to match.
#'
#' @return `sf` polygons object
#'
#' @export
extractFMU <- function(fmus, name) {
  dplyr::filter(fmus, FMU_NAME == name)
}

#' Extract boundary polygon(s) for LandWeb study areas
#'
#' @param name  A character (regex) string to match.
#'
#' @inheritParams prepFMAs
#'
#' @return `sf` polygons object
#'
#' @export
prepStudyArea <- function(name, destinationPath, targetCRS = LandWebCRS) {
  if (!grepl(paste(LandWebStudyAreas$Name, collapse = "|"), name)) {
    stop("name ", name, ", does not contain valid study area name.\n",
         "Study area name must be one of:\n", paste(LandWebStudyAreas$Name, collapse = ", "), ".")
  }

  if (grepl("random", name)) {
    studyArea <- prepRandomStudyArea(destinationPath, targetCRS, .seed = 867)
  } else if (grepl("FMU", name)) {
    studyArea <- prepFMUs(destinationPath, targetCRS) |> extractFMU(name)
  } else {
    studyArea <- prepFMAs(destinationPath, targetCRS) |> extractFMA(name)
  }
}

#' @export
#' @rdname prepStudyArea
prepRandomStudyArea <- function(destinationPath, targetCRS = LandWebCRS, .seed = 867) {
  message(crayon::red("Using random study area."))

  ansrs <- prepANSRs(destinationPath, targetCRS) ## TODO: why is this needed/used?

  withr::with_seed(.seed, {
    ## random area in Central-East AB
    SpaDES.tools::randomPolygon(ansrs, area = 4e10)
  })
}
