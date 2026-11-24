#' @param ccAgeLayer Which age layer to read for the current condition, when `vtm` is the
#'   current-condition map (`CurrentConditionVTM.tif`): `"TSF"`, time since fire
#'   (`CurrentConditionTSF.tif`), or `"SAM"`, stand age (`CurrentConditionSAM.tif`). Ignored for
#'   the maps of simulated years, whose age layer is `tsf`. Defaults to the option
#'   `LandWebUtils.ccAgeLayer`, else `"TSF"`, as since 1.0.1.
