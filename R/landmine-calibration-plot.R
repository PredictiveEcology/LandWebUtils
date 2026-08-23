utils::globalVariables(c(
  "areaHa", "bestval", "empirical", "fitted", "group", "hi", "iteration", "lo", "model1996",
  "objective", "pctIslandsExp", "pctIslandsObs", "parameter", "shapeExp", "shapeObs", "value",
  "what"
))

#' LandMine calibration diagnostic plots
#'
#' Figures for reporting what a fire-spread calibration achieved, and whether the objective can
#' actually tell good parameters from bad. They take plain data frames so the same functions
#' serve the module vignette, the project reports, and ad-hoc checks.
#'
#' - `landmine_plot_calibShape()` -- simulated shape index against Andison's fitted curve;
#' - `landmine_plot_calibIslands()` -- simulated remnant-island fraction against Andison Table 3.5;
#' - `landmine_plot_calibDiscrimination()` -- objective values for fitted vs random parameters;
#' - `landmine_plot_calibParams()` -- spread of repeat fits within the search bounds;
#' - `landmine_plot_calibConvergence()` -- DEoptim's best value against iteration, showing
#'   whether the run had converged or was still improving when it stopped.
#'
#' @return a `ggplot` object.
#'
#' @name landmine_calibration_plots
#' @rdname landmine_calibration_plots
NULL

#' @param components The `"components"` attribute of a [landmine_optim_fitAndison()] result: one
#'                   row per fire size (and replicate) with `areaHa`, `shapeObs`, `shapeExp`,
#'                   `pctIslandsObs`, `pctIslandsExp`.
#' @param shapeRangeHa Range Andison fitted the shape relationship over; shaded on the plot,
#'                     because the `log10(AREA)^4` term extrapolates sharply outside it.
#'
#' @export
#' @rdname landmine_calibration_plots
landmine_plot_calibShape <- function(components, shapeRangeHa = c(10, 3000)) {
  rng <- range(c(components$areaHa, shapeRangeHa))
  curve <- data.frame(areaHa = 10^seq(log10(rng[1]), log10(rng[2]), length.out = 200))
  curve$empirical <- landmine_optim_shapeTarget(curve$areaHa, intercept = 1.770)
  curve$model1996 <- landmine_optim_shapeTarget(curve$areaHa, intercept = 1.881)

  ggplot() +
    geom_rect(
      data = data.frame(lo = shapeRangeHa[1], hi = shapeRangeHa[2]),
      aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
      fill = "grey90", inherit.aes = FALSE
    ) +
    geom_line(data = curve, aes(x = areaHa, y = empirical, colour = "Andison empirical")) +
    geom_line(data = curve, aes(x = areaHa, y = model1996, colour = "Andison 1996 model"),
              linetype = "dashed") +
    geom_point(data = components, aes(x = areaHa, y = shapeObs, colour = "simulated"),
               size = 2, alpha = 0.8) +
    scale_x_log10() +
    scale_colour_manual(values = c("Andison empirical" = "black",
                                   "Andison 1996 model" = "grey40",
                                   "simulated" = "firebrick")) +
    labs(
      x = "Fire area (ha, log scale)", y = "Shape index (perimeter / circle of equal area)",
      colour = NULL,
      title = "Simulated fire shape against Andison's fitted relationship",
      subtitle = paste0("Shaded band is the ", shapeRangeHa[1], "-", shapeRangeHa[2],
                        " ha range the relationship was fitted over")
    ) +
    theme_bw()
}

#' @export
#' @rdname landmine_calibration_plots
landmine_plot_calibIslands <- function(components) {
  brk <- data.frame(
    lo = c(min(c(components$areaHa, 1)), 500, 1000),
    hi = c(500, 1000, max(c(components$areaHa, 1001))),
    value = c(4.0, 6.0, 9.0)
  )
  ggplot() +
    geom_segment(data = brk, aes(x = lo, xend = hi, y = value, yend = value,
                                 colour = "Andison Table 3.5 (empirical)"),
                 linewidth = 1) +
    geom_point(data = components, aes(x = areaHa, y = pctIslandsObs, colour = "simulated"),
               size = 2, alpha = 0.8) +
    scale_x_log10() +
    scale_colour_manual(values = c("Andison Table 3.5 (empirical)" = "black",
                                   "simulated" = "firebrick")) +
    labs(
      x = "Fire area (ha, log scale)", y = "Percent of fire area in remnant islands",
      colour = NULL,
      title = "Simulated remnant islands against Andison's empirical targets",
      subtitle = "Targets step at 500 and 1,000 ha; islands under 2 ha are excluded"
    ) +
    theme_bw()
}

#' @param scores A `data.frame` with `group` (e.g. `"fitted"` / `"random"`) and `objective`.
#'
#' @export
#' @rdname landmine_calibration_plots
landmine_plot_calibDiscrimination <- function(scores) {
  missing <- setdiff(c("group", "objective"), names(scores))
  if (length(missing)) {
    stop("`scores` needs column(s): ", paste(missing, collapse = ", "),
         ". Got: ", paste(names(scores), collapse = ", "), ".", call. = FALSE)
  }
  agg <- stats::aggregate(objective ~ group, data = scores, FUN = mean)
  sep <- if (nrow(agg) == 2L) {
    noise <- stats::sd(scores$objective[scores$group == agg$group[1]])
    sprintf("separation: %.2f sd", abs(diff(agg$objective)) / noise)
  } else {
    NULL
  }
  ggplot(scores, aes(x = group, y = objective, colour = group)) +
    geom_boxplot(outlier.shape = NA, colour = "grey40") +
    geom_jitter(position = position_jitter(width = 0.12, height = 0), size = 2, alpha = 0.85) +
    labs(
      x = NULL, y = "Objective value (lower is better)", colour = NULL,
      title = "Can the objective tell good parameters from bad?",
      subtitle = sep
    ) +
    theme_bw() +
    theme(legend.position = "none")
}

#' @param params A matrix or `data.frame` of fitted parameter sets, one row per calibration run
#'               (as returned by [landmine_optim_params_read()]).
#' @param lower,upper Search bounds used for the fit, in the same parameter order.
#'
#' @export
#' @rdname landmine_calibration_plots
landmine_plot_calibParams <- function(params, lower, upper) {
  params <- as.matrix(params)
  stopifnot(ncol(params) == length(lower), length(lower) == length(upper))
  nm <- colnames(params)
  if (is.null(nm)) nm <- paste0("par", seq_len(ncol(params)))

  ## express each fit as its position within its own search interval, so parameters on
  ## different scales are comparable: 0 = lower bound, 1 = upper bound
  rel <- sweep(params, 2, lower, "-")
  rel <- sweep(rel, 2, upper - lower, "/")
  long <- data.frame(
    parameter = factor(rep(nm, each = nrow(rel)), levels = nm),
    value = as.vector(rel)
  )
  ggplot(long, aes(x = parameter, y = value)) +
    geom_jitter(position = position_jitter(width = 0.12, height = 0),
                size = 2, alpha = 0.7, colour = "firebrick") +
    geom_hline(yintercept = c(0, 1), linetype = "dotted") +
    labs(
      x = NULL, y = "Position within search interval (0 = lower bound, 1 = upper)",
      title = "Spread of repeat calibration fits",
      subtitle = paste0(nrow(params), " fits; a parameter spanning the full interval is not ",
                        "identified by the objective")
    ) +
    theme_bw()
}

#' @param deoptim A `DEoptim` result object (or its `member$bestvalit` vector).
#'
#' @export
#' @rdname landmine_calibration_plots
landmine_plot_calibConvergence <- function(deoptim) {
  bv <- if (is.numeric(deoptim)) deoptim else deoptim$member$bestvalit
  if (is.null(bv)) {
    stop("could not find `member$bestvalit`; pass a DEoptim result or the trace vector.",
         call. = FALSE)
  }
  d <- data.frame(iteration = seq_along(bv), bestval = as.numeric(bv))
  lastImp <- suppressWarnings(max(which(diff(d$bestval) < 0) + 1))
  sub <- if (is.finite(lastImp)) {
    sprintf("last improvement at iteration %d of %d; %d further iterations gained nothing",
            lastImp, nrow(d), nrow(d) - lastImp)
  } else {
    "no improvement recorded"
  }
  p <- ggplot(d, aes(x = iteration, y = bestval)) +
    geom_line(colour = "firebrick", linewidth = 0.8) +
    labs(x = "DEoptim iteration", y = "Best objective value",
         title = "Calibration convergence", subtitle = sub) +
    theme_bw()
  if (is.finite(lastImp)) {
    p <- p + geom_vline(xintercept = lastImp, linetype = "dotted")
  }
  p
}
