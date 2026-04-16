# =============================================================================
# Shiny app launchers
# =============================================================================

#' Launch the interactive Bayesian Cure Modelling application
#'
#' Opens the Shiny app shipped with the package. The app lets the user
#' choose covariates for the clonogenic and kinetic sub-models, toggle an
#' interaction, adjust MCMC settings, run the Stan sampler, and inspect
#' results, diagnostics, identifiability tests and clinical plots.
#'
#' `cure()` is a convenience alias for `cure_app()`.
#'
#' @param ... Passed to [shiny::runApp()].
#' @return Invisibly `NULL`; called for side-effect.
#' @examples
#' \dontrun{
#'   rgetne::cure()
#' }
#' @export
cure_app <- function(...) {
  # The Shiny app uses `shiny` and `bslib` functions by bare name, so both
  # namespaces must be attached before `runApp()` is called. We attach them
  # here (not via `library()` at file top, which is disallowed inside R/).
  for (pkg in c("shiny", "bslib")) {
    if (!paste0("package:", pkg) %in% search()) {
      tryCatch(attachNamespace(pkg),
               error = function(e) {
                 stop(sprintf("Package '%s' is required but not installed.", pkg),
                      call. = FALSE)
               })
    }
  }

  app_dir <- system.file("shiny", package = "rgetne")
  if (!nzchar(app_dir))
    stop("Shiny app not found inside the package. Please reinstall rgetne.")
  shiny::runApp(app_dir, ...)
}

#' @rdname cure_app
#' @export
cure <- function(...) cure_app(...)
