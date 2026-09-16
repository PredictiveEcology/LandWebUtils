## ---- curated short names ---------------------------------------------------------------------
## Reporting units are keyed on their name: `refCode` (= `.slug(NAME_SHORT)`) names the parquet
## aggregate directory AND the figure filenames, so two units that key to the same string silently
## overwrite each other's outputs. The former keying -- `abbreviate(name, minlength = 8)` called one
## name at a time -- gives NO cross-unit uniqueness guarantee (R only makes abbreviations unique
## within a single call's vector). With the tenure x sub-region crossings that matters a lot: the
## v10 `FMA_NAME`s run to 97 characters and, crossed with the 19 Alberta Natural Subregions, only
## 327 of the 399 combinations abbreviate uniquely -- 68 collide (e.g. "Canadian Forest Products
## Ltd. Alpine" and "Crowsnest Forest Products Ltd. Alpine" both -> "CFPLtd.A").
##
## Vectorizing `abbreviate()` would fix the collisions but introduce a worse problem: the code for
## one unit would depend on which OTHER units happen to be present, so adding a reporting layer
## could silently rename another unit's outputs. Hence curated, pre-approved short names, validated
## for uniqueness and length when the tables load.

## Maximum length of a curated short name. Crossed reporting-unit names are
## "<tenure NAME_SHORT> <sub-region name>", so this bounds the crossed name (and therefore the
## figure title / facet strip label) at a width that still renders legibly.
.maxShortNameLength <- 16L

#' Validate a set of curated short names
#'
#' Enforces the invariants the output keying depends on: every short name is present,
#' filesystem-safe (unchanged by `.slug()`, so `refCode` round-trips), no longer than
#' `max_len`, and unique --- globally, or within each level of `by` (used for the per-tenure
#' landbase codes, which need only be unique inside their tenure).
#'
#' @param x character vector of candidate short names.
#' @param what character used in error messages to name the table being validated.
#' @param max_len integer maximum permitted length (default `16`).
#' @param by optional grouping vector the same length as `x`; uniqueness is then required
#'   only within each group.
#'
#' @return `x`, invisibly, if every check passes; otherwise an error.
#' @export
validateShortNames <- function(x, what = "NAME_SHORT", max_len = .maxShortNameLength, by = NULL) {
  x <- as.character(x)

  bad <- is.na(x) | !nzchar(x)
  if (any(bad)) {
    stop(what, ": ", sum(bad), " missing/empty short name(s).", call. = FALSE)
  }

  ## `.slug()`-stability: `refCode` is `.slug(NAME_SHORT)`, so a short name that the slug would
  ## rewrite (spaces, punctuation, leading/trailing separators) would not round-trip -- the
  ## on-disk key and the label in the figure title would drift apart.
  unstable <- x[.slug(x) != x]
  if (length(unstable)) {
    stop(
      what, ": ", length(unstable), " short name(s) are not filesystem-safe ",
      "(use only letters, digits, and single underscores): ",
      paste(sQuote(unstable), collapse = ", "), call. = FALSE
    )
  }

  long <- x[nchar(x) > max_len]
  if (length(long)) {
    stop(
      what, ": ", length(long), " short name(s) exceed ", max_len, " characters: ",
      paste(sprintf("%s (%d)", long, nchar(long)), collapse = ", "), call. = FALSE
    )
  }

  key <- if (is.null(by)) x else paste(by, x, sep = "\r")
  dup <- unique(x[duplicated(key)])
  if (length(dup)) {
    stop(
      what, ": duplicated short name(s)",
      if (!is.null(by)) " within a group" else "", ": ",
      paste(sQuote(dup), collapse = ", "),
      ". Short names key the parquet aggregates and figure files; duplicates overwrite.",
      call. = FALSE
    )
  }

  invisible(x)
}

#' Curated short names for the LandWeb tenures
#'
#' One row per member tenure of the v10 combined boundary layer
#' (`FMAs_LandWebltfc_map_v10_LCC1`), mapping its **coalesced identity** --- the same
#' `FMA_NAME`/`FMU_NAME`/`TSA_NUMB_1`/`FOREST_NAM`/`FML_NAME`/`Name` coalesce that
#' [build_studyarea_crosswalk()] uses as `fma_name` --- to a curated `NAME_SHORT`.
#'
#' `NAME_SHORT` is the partner-facing label: it names the tenure in every crossed
#' reporting unit ("`<NAME_SHORT> <sub-region>`"), in the figure titles and facet strips,
#' and (via `.slug()`) in the `refCode` that keys the parquet aggregates and figure files.
#' Codes are derived from the **v10 identity**, not from the legacy v2 study-area registry
#' ([LandWebStudyAreas]), because tenures change hands and the v10 layer is the authority.
#'
#' Codes are hand-curated rather than generated so that adding a reporting layer, or a new
#' member to the boundary layer, can never silently rename an existing unit's outputs.
#' Adding a member to the v10 layer therefore requires adding its code here --- an
#' uncurated member is an error, not a silent fallback (see [build_studyarea_crosswalk()]).
#'
#' @return A `tibble` with columns `fma_name` (the coalesced v10 identity) and `NAME_SHORT`.
#' @seealso [validateShortNames()], [build_studyarea_crosswalk()]
#' @export
tenureShortNames <- function() {
  out <- tibble::tribble(
    ~fma_name                                                                                         , ~NAME_SHORT        ,
    ## ---- Alberta (FMA_NAME) --------------------------------------------------------------------
    "ANC Timber Ltd."                                                                                 , "ANC"              ,
    "Alberta-Pacific Forest Industries Inc."                                                          , "AlPac"            ,
    "Blue Ridge Lumber Inc."                                                                          , "BlueRidge"        ,
    "Canadian Forest Products Ltd."                                                                   , "Canfor"           ,
    "Canfor (Whitecourt) Forest Products Ltd."                                                        , "CanforWhitecourt" ,
    "Crowsnest Forest Products Ltd."                                                                  , "Crowsnest"        ,
    "Manning Forest Products Ltd."                                                                    , "Manning"          ,
    "Mercer Peace River Pulp Ltd. (East)"                                                             , "Mercer_E"         ,
    "Mercer Peace River Pulp Ltd. (West)"                                                             , "Mercer_W"         ,
    "Millar Western Forest Products Ltd."                                                             , "MillarWestern"    ,
    "Spray Lake Sawmills (1980) Ltd."                                                                 , "SprayLake"        ,
    "Sundre Forest Products Inc."                                                                     , "Sundre"           ,
    "Tolko Industries Ltd. (High Prairie)"                                                            , "Tolko_HP"         ,
    ## the two Tolko-led consortia: coded for their partners, since "Tolko" alone is ambiguous
    "Tolko Industries Ltd., Norbord Inc. and La Crete Sawmills Ltd."                                  , "Tolko_Norbord_LC" ,
    "Tolko Industries Ltd., Vanderwell Contractors (1971) Ltd. and West Fraser Mills Ltd. (Slave Lake)", "Tolko_Vand_WF_SL",
    "Vanderwell Contractors (1971) Ltd."                                                              , "Vanderwell"       ,
    "West Fraser Mills Ltd. (Edson)"                                                                  , "WF_Edson"         ,
    "West Fraser Mills Ltd. (Hinton)"                                                                 , "WF_Hinton"        ,
    "West Fraser Mills Ltd. and Tolko Industries Ltd."                                                , "WF_Tolko"         ,
    "Weyerhaeuser Company Limited (Grande Prairie)"                                                   , "WeyCo_GrandePr"   ,
    "Weyerhaeuser Company Limited (Pembina Timberland)"                                               , "WeyCo_Pembina"    ,

    ## ---- British Columbia (TSA_NUMB_1) ---------------------------------------------------------
    "Cranbrook TSA"                                                                                   , "Cranbrook_TSA"    ,
    "Dawson Creek TSA"                                                                                , "DawsonCreek_TSA"  ,
    "Fort Nelson TSA"                                                                                 , "FortNelson_TSA"   ,
    "Fort St. John TSA"                                                                               , "FortStJohn_TSA"   ,
    "Golden TSA"                                                                                      , "Golden_TSA"       ,
    "Invermere TSA"                                                                                   , "Invermere_TSA"    ,
    "MacKenzie TSA"                                                                                   , "Mackenzie_TSA"    ,
    "Prince George TSA"                                                                               , "PrinceGeo_TSA"    ,
    "Robson Valley TSA"                                                                               , "RobsonValley_TSA" ,

    ## ---- Saskatchewan (FOREST_NAM) -------------------------------------------------------------
    ## province-prefixed: several SK forest names ("North Central", "Prince Albert") are generic
    ## enough to be confusing once they appear beside BC/ON units in a cross-study-area report.
    "Bronson-Green Lake"                                                                              , "SK_BronsonGreen"  ,
    "Glaslyn"                                                                                         , "SK_Glaslyn"       ,
    "Meadow Lake"                                                                                     , "SK_MeadowLake"    ,
    "North Central"                                                                                   , "SK_NorthCentral"  ,
    "North East"                                                                                      , "SK_NorthEast"     ,
    "Pasquia Porcupine"                                                                               , "SK_PasquiaPorc"   ,
    "Prince Albert"                                                                                   , "SK_PrinceAlbert"  ,

    ## ---- Ontario (FMU_NAME) --------------------------------------------------------------------
    ## "<X> Forest" -> "ON_<X>": the "Forest" suffix carries no information (every ON unit has it),
    ## and the prefix keeps e.g. "Caribou Forest" distinct from the Caribou Ranges reporting layer.
    "Boundary Waters Forest"                                                                          , "ON_BoundaryWtrs"  ,
    "Caribou Forest"                                                                                  , "ON_Caribou"       ,
    "Dryden Forest"                                                                                   , "ON_Dryden"        ,
    "English River Forest"                                                                            , "ON_EnglishRiver"  ,
    "Kenogami Forest"                                                                                 , "ON_Kenogami"      ,
    "Kenora Forest"                                                                                   , "ON_Kenora"        ,
    "Lac Seul Forest"                                                                                 , "ON_LacSeul"       ,
    "Ogoki Forest"                                                                                    , "ON_Ogoki"         ,
    "Red Lake Forest"                                                                                 , "ON_RedLake"       ,
    "Trout Lake Forest"                                                                               , "ON_TroutLake"     ,
    "Wabigoon Forest"                                                                                 , "ON_Wabigoon"      ,
    "Whiskey Jack Forest"                                                                             , "ON_WhiskeyJack"   ,
    "Whitefeather Forest"                                                                             , "ON_Whitefeather"  ,

    ## ---- Manitoba (FML_NAME) -------------------------------------------------------------------
    "FML-2"                                                                                           , "MB_FML2"          ,
    "FML-3"                                                                                           , "MB_FML3"          ,

    ## ---- Northwest Territories (Name; these rows carry no FMA_NAME) -----------------------------
    "Fort Providence"                                                                                 , "NT_FtProvidence"  ,
    "Fort Resolution"                                                                                 , "NT_FtResolution"
  )

  validateShortNames(out$NAME_SHORT, what = "tenureShortNames()$NAME_SHORT")
  out
}

#' Output key (`refCode`) for one analysis kind and reporting layer
#'
#' The filesystem key an analysis writes under: it names the parquet aggregate directory
#' (`_aggregates/<refCode>/`) and the figure/CSV filenames. Built by slugging the layer's
#' name --- which is its curated `NAME_SHORT`, or `"<tenure> <NAME_SHORT>"` for a crossed
#' unit --- so the on-disk key is exactly the label shown in the figures, and two units can
#' only collide if their curated short names do (which [validateShortNames()] forbids).
#'
#' This replaces the former `abbreviate(layer, minlength = 8)` keying, which gave no
#' cross-unit uniqueness guarantee: 68 of the 399 crossed FMA x ANSR names abbreviated onto
#' a name already taken, silently overwriting each other's aggregates and figures.
#'
#' @param kind character analysis kind (`"lm"`, `"pm"`, `"lw"`, `"sspm"`, ...).
#' @param layer character reporting-layer name (a `reportingPolygons` list name).
#'
#' @return a character `refCode`.
#' @export
refCodeFor <- function(kind, layer) {
  paste0(kind, "_", .slug(layer))
}

#' Look up curated short names
#'
#' Vectorised lookup of `NAME_SHORT` for a set of names, against a two-column curated
#' table (e.g. [tenureShortNames()]).
#'
#' @param x character vector of long names to look up.
#' @param table a two-column curated table; the first column holds the long names and
#'   `NAME_SHORT` the short ones (default [tenureShortNames()]).
#'
#' @return character vector the same length as `x`; `NA` where `x` has no curated code.
#' @export
shortNameFor <- function(x, table = tenureShortNames()) {
  table[["NAME_SHORT"]][match(as.character(x), table[[1L]])]
}
