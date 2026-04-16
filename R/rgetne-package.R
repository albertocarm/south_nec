#' rgetne: Bayesian Promotion Time Cure Model Toolkit and Interactive App
#'
#' Fits Bayesian Promotion Time Cure Models (Yakovlev & Tsodikov bounded
#' cumulative hazard) via Stan, supports interactions and per-variable
#' priors, ships identifiability diagnostics and clinical plots, and
#' bundles an interactive Shiny app (see [cure()]).
#'
#' @section Quick start:
#' \preformatted{
#'   library(rgetne)
#'   cure()                       # launches the Shiny app
#'
#'   dat  <- prep_data(df_datos_app)
#'   mod  <- cure_model(dat,
#'     clonogenic = ~ ki67_percent * periop_therapy + stage_III,
#'     kinetic    = ~ ki67_percent + periop_therapy)
#'   cfit <- fit_model(mod)
#'   summary_coefs(cfit)
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#'
#' @importFrom rstan stan rstan_options check_divergences check_treedepth
#' @importFrom bayesplot mcmc_trace
#' @importFrom bslib page_sidebar
#' @importFrom stats as.formula complete.cases cor median model.matrix
#'   quantile sd setNames
#' @importFrom utils head capture.output
#' @importFrom grDevices rgb
#' @importFrom survival Surv survfit
#' @importFrom graphics abline legend lines par plot polygon text
## usethis namespace: end
NULL

# ggplot2 aes() uses non-standard evaluation on column names. Register those
# names so R CMD check does not warn about "no visible binding for global
# variable".
utils::globalVariables(c(
  # plot_continuous_effect
  "x", "theta_lo", "theta_hi", "theta_med",
  "cure_lo", "cure_hi", "cure_med",
  "mtime_lo", "mtime_hi", "mtime_med",
  # plot_contrast_cure_or
  "at", "or_med", "or_lo", "or_hi",
  "sig", "group", "lo", "hi", "cure",
  # plot_km / data_summary_table
  "variable", "value", "n", "pct"
))
