# =============================================================================
# Standalone launcher for shinyapps.io / rsconnect / any Shiny Server.
#
# This file is NOT part of the rgetne package source (it is listed in
# .Rbuildignore). It is a deployment script that makes the Shiny app run
# from a plain repo checkout without having to install rgetne first.
#
# Only the three libraries that the Shiny app truly needs attached to the
# search path are loaded here:
#   * shiny / bslib  -> the app UI uses their functions by bare name
#   * rstan          -> so `rstan_options()` can be called below
# Everything else (loo, bayesplot, ggplot2, gridExtra, ...) is accessed via
# fully qualified `pkg::fn()` calls from R/toolkit.R and requires no
# attachment. rsconnect still picks those packages up from the parse tree
# and installs them during deployment.
#
# Deployment:
#     rsconnect::deployApp(".", appPrimaryDoc = "app.R")
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(rstan)
})

rstan_options(auto_write = TRUE)
if (is.null(getOption("mc.cores"))) {
  options(mc.cores = max(1L, parallel::detectCores(logical = FALSE) - 1L))
}

APP_ROOT <- normalizePath(".", mustWork = TRUE)

# --- Load toolkit functions directly (no package install required) ----------
toolkit_path <- file.path(APP_ROOT, "R", "toolkit.R")
if (!file.exists(toolkit_path))
  stop("Could not locate R/toolkit.R next to app.R.")
source(toolkit_path, local = FALSE)

# --- Resolve repo-local paths to bundled assets -----------------------------
STAN_FILE <- file.path(APP_ROOT, "inst", "stan", "promotion_time_cure_v2.stan")
DATA_FILE <- file.path(APP_ROOT, "inst", "extdata", "df_datos_app.rds")
if (!file.exists(STAN_FILE)) stop("Missing Stan model: ", STAN_FILE)
if (!file.exists(DATA_FILE)) stop("Missing example data: ", DATA_FILE)

stan_model_file <- function() STAN_FILE
example_data    <- function() DATA_FILE

# Precomputed fit (avoids running Stan on the server)
PRECOMPUTED_FILE <- file.path(APP_ROOT, "inst", "extdata", "cfit_default.rds")
precomputed_cfit_path <- function() {
  if (file.exists(PRECOMPUTED_FILE)) PRECOMPUTED_FILE else NULL
}

stan_file_cached <- function() {
  src <- stan_model_file()
  dst <- file.path(tempdir(), basename(src))
  if (!file.exists(dst) || file.info(src)$mtime > file.info(dst)$mtime) {
    file.copy(src, dst, overwrite = TRUE)
  }
  dst
}

# --- Hand over to the Shiny app ---------------------------------------------
source(file.path(APP_ROOT, "inst", "shiny", "app.R"), local = FALSE)
