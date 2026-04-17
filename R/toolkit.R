# =============================================================================
# Promotion Time Cure Model -- Toolkit (package port of v4)
# =============================================================================
# Supports:
#   - Interactions (A*B, A:B)
#   - Splines (ns(), bs(), rcs via rms::rcs -> user-supplied dummies)
#   - Granular per-term priors
#   - rms-style contrasts for interactions
#   - Configurable MCMC (chains, iter, warmup, adapt_delta, ...)
# =============================================================================

.onAttach <- function(libname, pkgname) {
  # Safe defaults for Stan
  try(rstan::rstan_options(auto_write = TRUE), silent = TRUE)
  if (is.null(getOption("mc.cores"))) {
    options(mc.cores = max(1L, parallel::detectCores(logical = FALSE) - 1L))
  }
}

# =============================================================================
#' Path to the bundled Stan model
#'
#' @return Absolute path to `promotion_time_cure_v2.stan` shipped with the
#'   package.
#' @export
stan_model_file <- function() {
  p <- system.file("stan", "promotion_time_cure_v2.stan", package = "rgetne")
  if (!nzchar(p)) stop("Stan model not found; please reinstall rgetne.")
  p
}

#' Path to the bundled example dataset
#'
#' @return Absolute path to `df_datos_app.rds` shipped with the package.
#' @export
example_data <- function() {
  p <- system.file("extdata", "df_datos_app.rds", package = "rgetne")
  if (!nzchar(p)) stop("Example data not found; please reinstall rgetne.")
  p
}

# =============================================================================
#' Path to the publication dataset (English factor labels)
#'
#' Returns the path to `dataset_publicacion_english.rds` bundled under
#' `inst/extdata/`. This dataset contains the original factor-coded
#' variables (`EFS_Time`, `EFS_Status`, `OS_Time`, `OS_Status`, `Sex`,
#' `TNM_Stage`, `Primary_Surgery`, `Site_4_Groups`, `Ki67_Index`,
#' `Perioperative_Chemo`, `ECOG_PS`) used by the Cox multivariable tab
#' and the Kaplan-Meier exploration tab of the Shiny app.
#'
#' @return Character path. Errors if the file is not installed.
#' @export
publication_data <- function() {
  p <- system.file("extdata", "dataset_publicacion_english.rds",
                   package = "rgetne")
  if (!nzchar(p)) {
    stop("Publication dataset not found. Place ",
         "'dataset_publicacion_english.rds' in inst/extdata/ and reinstall.")
  }
  p
}

# =============================================================================
#' Load an RDS dataset with an EFS_time column
#'
#' Reads the file and runs it through [prep_data()] so that raw column names
#' (`time_efs`, `event_status`, factor-encoded covariates, ...) are mapped to
#' the modelling layout (`EFS_time`, `.status_bin`, binary dummies). Rows with
#' non-positive `EFS_time` are dropped by `prep_data()`.
#'
#' @param path Path to an `.rds` file.
#' @return The data.frame ready for [cure_model()] / [fit_cure_bayes()].
#' @export
load_data <- function(path) {
  df <- readRDS(path)
  df <- prep_data(df)
  message(sprintf("Data loaded: %d observations (EFS_time > 0)", nrow(df)))
  df
}

# =============================================================================
#' Pre-process the shipped example dataset for modelling
#'
#' Converts factors to integer dummies (so that `model.matrix` column names
#' equal the input column names) and sets up `EFS_time` / `.status_bin`.
#' This is the transformation expected by [cure_model()] and the Shiny app.
#'
#' @param df A data.frame with columns `time_efs`, `event_status`, and any of
#'   `primary_surgery`, `disease_stage`, `primary_site`, `ki67_percent`,
#'   `periop_therapy`.
#' @return A data.frame with standardised numeric columns:
#'   `primary_surgery_yes`, `stage_II`, `stage_III`, `site_colorectal`,
#'   `site_pancreas`, `ki67_percent`, `periop_therapy`, plus `EFS_time` and
#'   `.status_bin`.
#' @export
prep_data <- function(df) {
  stopifnot(is.data.frame(df))

  if (!"EFS_time" %in% names(df)) {
    if ("time_efs" %in% names(df))      df$EFS_time <- df$time_efs
    else if ("time" %in% names(df))     df$EFS_time <- df$time
    else stop("Input must contain a column named 'EFS_time', 'time_efs' or 'time'.")
  }
  if (!".status_bin" %in% names(df)) {
    if ("event_status" %in% names(df))  df$.status_bin <- as.integer(df$event_status)
    else if ("status" %in% names(df))   df$.status_bin <- as.integer(df$status)
    else if ("event" %in% names(df))    df$.status_bin <- as.integer(df$event)
    else stop("Input must contain a column named '.status_bin', 'event_status', 'status' or 'event'.")
  }

  df <- df[!is.na(df$EFS_time) & df$EFS_time > 0, , drop = FALSE]

  if ("primary_surgery" %in% names(df)) {
    df$primary_surgery_yes <- as.integer(df$primary_surgery == "Yes")
  }
  if ("disease_stage" %in% names(df)) {
    df$stage_II  <- as.integer(df$disease_stage == "II")
    df$stage_III <- as.integer(df$disease_stage == "III")
  }
  if ("primary_site" %in% names(df)) {
    df$site_colorectal <- as.integer(df$primary_site == "Colorectal")
    df$site_pancreas   <- as.integer(df$primary_site %in% c("P\u00e1ncreas", "Pancreas"))
  }
  if ("periop_therapy" %in% names(df)) {
    df$periop_therapy <- as.integer(df$periop_therapy)
  }

  df
}

# =============================================================================
#' Define a Promotion Time Cure Model specification
#'
#' @param dat A data.frame (already pre-processed, see [prep_data()]) with
#'   `EFS_time` and `.status_bin`.
#' @param clonogenic A one-sided formula for the incidence (log-theta) linear
#'   predictor, e.g. `~ A * B + C`.
#' @param kinetic A one-sided formula for the latency (log-lambda) linear
#'   predictor.
#' @param prior_clonogenic Scalar, unnamed vector, or named list of prior SDs
#'   for the clonogenic coefficients (see details).
#' @param prior_kinetic Same structure as `prior_clonogenic`, for latency.
#' @param prior_intercept Scalar prior SD for both intercepts.
#' @param stan_file Path to the `.stan` file. Defaults to the bundled model.
#' @param spline_vars Optional named list of pre-computed spline bases.
#' @param scale Logical. If `TRUE` (default, back-compatible) the design
#'   matrices are internally centred and divided by their SD before entering
#'   Stan. If `FALSE`, matrices are used as-is and `X_center` / `Z_center` are
#'   set to 0 and `X_scale` / `Z_scale` to 1, so that downstream helpers
#'   treating the linear predictor are identity transforms. Use `FALSE` when
#'   predictors have already been standardised upstream via
#'   [preprocess_predictors()], which lets every per-variable prior be
#'   interpreted on the same standardised scale.
#'
#' @details
#' `prior_*` accepts:
#' * a single scalar (same SD for every coefficient),
#' * a numeric vector of length K (positional),
#' * a named list: every name is grep-matched against the `model.matrix`
#'   columns; matching columns receive the supplied SD, the rest default
#'   to 1.
#'
#' The default `prior_clonogenic` is a shared SD of `2.5` across every
#' coefficient (weakly informative). The default `prior_kinetic` is a
#' named list that sets a tight SD of `0.25` on `ki67_percent` and
#' `periop_therapy` (the covariates that tend to dominate the kinetic
#' ridge in this cure model) and leaves every other coefficient at the
#' fallback SD of `1.0`. Pass a scalar, a vector or a custom list to
#' override.
#'
#' @return An object of class `cure_spec`.
#' @export
cure_model <- function(dat, clonogenic, kinetic,
                       prior_clonogenic = 2.5,
                       prior_kinetic    = list(ki67_percent   = 0.25,
                                               periop_therapy = 0.25),
                       prior_intercept  = 2.5,
                       stan_file = stan_model_file(),
                       spline_vars = NULL,
                       scale = TRUE) {

  all_terms <- unique(c(all.vars(clonogenic), all.vars(kinetic)))
  needed    <- unique(c("EFS_time", ".status_bin", all_terms))
  missing_cols <- setdiff(needed, names(dat))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in `dat`: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  dat_cc <- dat[complete.cases(dat[, needed, drop = FALSE]), , drop = FALSE]
  message(sprintf("Observations: %d -> %d complete", nrow(dat), nrow(dat_cc)))

  X_raw <- model.matrix(clonogenic, data = dat_cc)[, -1, drop = FALSE]
  Z_raw <- model.matrix(kinetic,    data = dat_cc)[, -1, drop = FALSE]

  inject <- function(M, lst) {
    if (is.null(lst)) return(M)
    for (nm in names(lst)) {
      sp_mat <- lst[[nm]]
      if (is.matrix(sp_mat)) {
        M <- M[, !colnames(M) %in% nm, drop = FALSE]
        colnames(sp_mat) <- paste0(nm, "_sp", seq_len(ncol(sp_mat)))
        M <- cbind(M, sp_mat[as.numeric(rownames(dat_cc)), , drop = FALSE])
      }
    }
    M
  }
  X_raw <- inject(X_raw, spline_vars$clonogenic)
  Z_raw <- inject(Z_raw, spline_vars$kinetic)

  stopifnot(nrow(X_raw) == nrow(Z_raw))

  X_names <- colnames(X_raw)
  Z_names <- colnames(Z_raw)

  if (ncol(X_raw) == 0 || ncol(Z_raw) == 0) {
    stop("Both clonogenic and kinetic formulas must supply at least one covariate.")
  }

  message(sprintf("N = %d, events = %d (%.0f%%), K_clono = %d, K_kinet = %d",
                  nrow(dat_cc), sum(dat_cc$.status_bin),
                  100 * mean(dat_cc$.status_bin), ncol(X_raw), ncol(Z_raw)))

  if (isTRUE(scale)) {
    X <- scale(X_raw); Z <- scale(Z_raw)
    X_center <- attr(X, "scaled:center"); X_sc <- attr(X, "scaled:scale")
    Z_center <- attr(Z, "scaled:center"); Z_sc <- attr(Z, "scaled:scale")
    X[is.na(X)] <- X_raw[is.na(X)]
    Z[is.na(Z)] <- Z_raw[is.na(Z)]
    X_sc[is.na(X_sc) | X_sc == 0] <- 1
    Z_sc[is.na(Z_sc) | Z_sc == 0] <- 1
  } else {
    X <- X_raw; Z <- Z_raw
    X_center <- setNames(rep(0, ncol(X_raw)), X_names)
    Z_center <- setNames(rep(0, ncol(Z_raw)), Z_names)
    X_sc     <- setNames(rep(1, ncol(X_raw)), X_names)
    Z_sc     <- setNames(rep(1, ncol(Z_raw)), Z_names)
  }

  resolve_prior <- function(spec, var_names, label) {
    k <- length(var_names)
    if (is.list(spec)) {
      sds <- rep(1.0, k); names(sds) <- var_names
      for (nm in names(spec)) {
        idx <- grep(nm, var_names, fixed = TRUE)
        if (length(idx) > 0) sds[idx] <- spec[[nm]]
      }
      sds
    } else if (length(spec) == 1) {
      rep(spec, k)
    } else if (length(spec) == k) {
      spec
    } else {
      stop(sprintf("%s prior: length %d doesn't match %d covariates",
                   label, length(spec), k))
    }
  }

  pr_beta  <- resolve_prior(prior_clonogenic, X_names, "clonogenic")
  pr_gamma <- resolve_prior(prior_kinetic,    Z_names, "kinetic")

  stan_data <- list(
    N = nrow(dat_cc), K_inc = ncol(X), K_lat = ncol(Z),
    X = as.matrix(X), Z = as.matrix(Z),
    t = dat_cc$EFS_time, d = as.integer(dat_cc$.status_bin),
    prior_sd_beta       = as.array(pr_beta),
    prior_sd_gamma      = as.array(pr_gamma),
    prior_sd_intercepts = prior_intercept
  )

  structure(list(
    stan_data = stan_data, stan_file = stan_file,
    dat = dat_cc,
    X_raw = X_raw, Z_raw = Z_raw,
    X_names = X_names, Z_names = Z_names,
    X_center = X_center, X_scale = X_sc,
    Z_center = Z_center, Z_scale = Z_sc,
    prior_beta = pr_beta, prior_gamma = pr_gamma,
    prior_intercept = prior_intercept,
    formula_clono = clonogenic, formula_kinet = kinetic
  ), class = "cure_spec")
}

# =============================================================================
#' Fit a Promotion Time Cure Model specification via Stan
#'
#' @param mod           A `cure_spec` object from [cure_model()].
#' @param chains,iter,warmup  MCMC settings.
#' @param adapt_delta   NUTS target acceptance.
#' @param max_treedepth NUTS max treedepth.
#' @param seed          RNG seed.
#' @param ...           Passed through to [rstan::stan()].
#'
#' @return A `cure_fit` list containing `fit` (stanfit) and `spec` (cure_spec).
#' @export
fit_model <- function(mod, chains = 4, iter = 4000, warmup = floor(iter / 2),
                      adapt_delta = 0.95, max_treedepth = 12, seed = 42, ...) {
  fit <- rstan::stan(
    file   = mod$stan_file,
    data   = mod$stan_data,
    chains = chains, iter = iter, warmup = warmup,
    seed   = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    ...
  )
  structure(list(fit = fit, spec = mod), class = "cure_fit")
}

# =============================================================================
#' MCMC diagnostic summary (R-hat, n_eff, divergences, tree-depth)
#' @param cfit A `cure_fit` object.
#' @export
mcmc_checks <- function(cfit) {
  if (!inherits(cfit, "cure_fit"))
    stop("'cfit' must be an object returned by fit_model().")
  fit <- cfit$fit; sp <- cfit$spec
  if (!inherits(fit, "stanfit"))
    stop("cfit$fit is not a stanfit object; rerun fit_model().")

  cat("\n=== MCMC CONVERGENCE SUMMARY ===\n")
  cat("Target: R-hat < 1.01 and n_eff > 400 for every monitored parameter.\n\n")

  pars <- c("beta0",
            paste0("beta[",  seq_len(sp$stan_data$K_inc), "]"),
            "gamma0",
            paste0("gamma[", seq_len(sp$stan_data$K_lat), "]"),
            "alpha")

  s <- tryCatch(
    rstan::summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary,
    error = function(e) {
      message("Could not build parameter summary: ", conditionMessage(e))
      NULL
    }
  )

  if (is.matrix(s) || is.data.frame(s)) {
    print(round(as.matrix(s), 3))

    rh  <- s[, "Rhat"]
    neff <- s[, "n_eff"]
    bad_rh  <- which(rh   > 1.01 & !is.na(rh))
    bad_neff <- which(neff < 400 & !is.na(neff))

    cat("\n-- R-hat check --\n")
    if (length(bad_rh) == 0) {
      cat("OK: all R-hat <= 1.01 (chains have mixed).\n")
    } else {
      cat(sprintf("[!] %d parameter(s) with R-hat > 1.01:\n",
                  length(bad_rh)))
      print(round(rh[bad_rh], 4))
    }

    cat("\n-- Effective sample size check --\n")
    if (length(bad_neff) == 0) {
      cat("OK: all n_eff >= 400.\n")
    } else {
      cat(sprintf("[!] %d parameter(s) with n_eff < 400:\n",
                  length(bad_neff)))
      print(round(neff[bad_neff]))
    }
  } else {
    cat("(Summary matrix not available.)\n")
  }

  cat("\n-- Sampler diagnostics --\n")
  rstan::check_divergences(fit)
  rstan::check_treedepth(fit)
  invisible(NULL)
}

# =============================================================================
#' Coefficient summary on the original (un-scaled) covariate scale
#' @param cfit A `cure_fit` object.
#' @export
summary_coefs <- function(cfit) {
  sp <- cfit$spec; post <- as.data.frame(cfit$fit)

  cat("\n===================================================================\n")
  cat("  COEFFICIENTS (original scale)\n")
  cat("===================================================================\n")

  results <- list()

  cat("\n--- CLONOGENIC (incidence, log-theta) ---\n")
  cat("    positive = more clonogenic cells = LOWER cure rate\n\n")
  for (j in seq_along(sp$X_names)) {
    col  <- paste0("beta[", j, "]")
    orig <- post[[col]] / sp$X_scale[j]
    med  <- median(orig); lo <- quantile(orig, 0.025); hi <- quantile(orig, 0.975)
    sig  <- ifelse(lo > 0, " *+", ifelse(hi < 0, " *-", ""))
    cat(sprintf("  %-45s: %6.3f [%6.3f, %6.3f]%s\n",
                sp$X_names[j], med, lo, hi, sig))
    results[[paste0("clono_", sp$X_names[j])]] <- c(median = med, lo = lo, hi = hi)
  }

  cat("\n--- KINETIC (latency, log-lambda Weibull) ---\n")
  cat("    negative = FASTER relapse\n\n")
  for (j in seq_along(sp$Z_names)) {
    col  <- paste0("gamma[", j, "]")
    orig <- post[[col]] / sp$Z_scale[j]
    med  <- median(orig); lo <- quantile(orig, 0.025); hi <- quantile(orig, 0.975)
    sig  <- ifelse(lo > 0, " *+", ifelse(hi < 0, " *-", ""))
    cat(sprintf("  %-45s: %6.3f [%6.3f, %6.3f]%s\n",
                sp$Z_names[j], med, lo, hi, sig))
    results[[paste0("kinet_", sp$Z_names[j])]] <- c(median = med, lo = lo, hi = hi)
  }

  cat(sprintf("\n  %-45s: %6.3f [%6.3f, %6.3f]\n", "alpha (Weibull shape)",
              median(post$alpha), quantile(post$alpha, 0.025),
              quantile(post$alpha, 0.975)))

  cure <- exp(-exp(post$beta0))
  cat(sprintf("\n  Cure fraction (covariates at mean): %.1f%% [%.1f%%, %.1f%%]\n",
              median(cure) * 100, quantile(cure, 0.025) * 100,
              quantile(cure, 0.975) * 100))

  invisible(results)
}

# =============================================================================
#' Compute a scenario contrast for clonogenic or kinetic components
#'
#' @param cfit       A `cure_fit` object.
#' @param component  Either `"clonogenic"` or `"kinetic"`.
#' @param var        Name of the variable whose effect we want.
#' @param levels     Named list: covariate settings in the "exposed" scenario.
#' @param ref        Named list: covariate settings in the "reference" scenario.
#' @param at         Named list: levels of interacting variables to evaluate at.
#' @param prob       Credible interval width (default 0.95).
#' @export
cure_contrast <- function(cfit,
                          component = c("clonogenic", "kinetic"),
                          var, levels, ref, at = list(), prob = 0.95) {

  component <- match.arg(component)
  sp   <- cfit$spec
  post <- as.data.frame(cfit$fit)

  if (component == "clonogenic") {
    var_names   <- sp$X_names
    center      <- sp$X_center; sc <- sp$X_scale
    coef_cols   <- paste0("beta[", seq_along(var_names), "]")
  } else {
    var_names   <- sp$Z_names
    center      <- sp$Z_center; sc <- sp$Z_scale
    coef_cols   <- paste0("gamma[", seq_along(var_names), "]")
  }

  grid <- if (length(at) == 0) data.frame(row = 1) else expand.grid(at)

  alpha_lo <- (1 - prob) / 2
  alpha_hi <- 1 - alpha_lo

  results <- data.frame()

  for (g in seq_len(nrow(grid))) {
    x_lev <- rep(0, length(var_names))
    x_ref <- rep(0, length(var_names))
    names(x_lev) <- names(x_ref) <- var_names

    for (nm in names(at)) {
      val <- grid[g, nm]
      for (j in seq_along(var_names)) {
        vn <- var_names[j]
        if (nm == vn || grepl(nm, vn, fixed = TRUE)) {
          x_lev[j] <- (val - center[j]) / sc[j]
          x_ref[j] <- (val - center[j]) / sc[j]
        }
      }
    }
    for (nm in names(levels)) {
      for (j in seq_along(var_names)) {
        if (var_names[j] == nm) x_lev[j] <- (levels[[nm]] - center[j]) / sc[j]
      }
    }
    for (nm in names(ref)) {
      for (j in seq_along(var_names)) {
        if (var_names[j] == nm) x_ref[j] <- (ref[[nm]] - center[j]) / sc[j]
      }
    }

    for (j in seq_along(var_names)) {
      vn <- var_names[j]
      if (grepl(":", vn, fixed = TRUE)) {
        parts <- strsplit(vn, ":")[[1]]
        val_a_lev <- NA; val_b_lev <- NA
        val_a_ref <- NA; val_b_ref <- NA
        for (p in parts) {
          raw_val_lev <- NULL; raw_val_ref <- NULL
          if (p %in% names(levels)) raw_val_lev <- levels[[p]]
          if (p %in% names(ref))    raw_val_ref <- ref[[p]]
          if (p %in% names(at)) {
            raw_val_lev <- grid[g, p]; raw_val_ref <- grid[g, p]
          }
          if (is.null(raw_val_lev)) raw_val_lev <- center[p]
          if (is.null(raw_val_ref)) raw_val_ref <- center[p]
          if (is.na(val_a_lev)) {
            val_a_lev <- raw_val_lev; val_a_ref <- raw_val_ref
          } else {
            val_b_lev <- raw_val_lev; val_b_ref <- raw_val_ref
          }
        }
        raw_int_lev <- val_a_lev * val_b_lev
        raw_int_ref <- val_a_ref * val_b_ref
        x_lev[j] <- (raw_int_lev - center[j]) / sc[j]
        x_ref[j] <- (raw_int_ref - center[j]) / sc[j]
      }
    }

    coef_mat <- as.matrix(post[, coef_cols])
    lp_lev <- as.vector(coef_mat %*% x_lev)
    lp_ref <- as.vector(coef_mat %*% x_ref)
    delta  <- lp_lev - lp_ref

    row <- data.frame(
      contrast = sprintf("%s: %s vs %s", var,
                         paste(names(levels), levels, sep = "=", collapse = ","),
                         paste(names(ref),    ref,    sep = "=", collapse = ",")),
      delta_median     = median(delta),
      delta_lo         = quantile(delta, alpha_lo),
      delta_hi         = quantile(delta, alpha_hi),
      exp_delta_median = median(exp(delta)),
      exp_delta_lo     = quantile(exp(delta), alpha_lo),
      exp_delta_hi     = quantile(exp(delta), alpha_hi),
      prob_pos         = mean(delta > 0),
      prob_neg         = mean(delta < 0),
      row.names        = NULL
    )
    for (nm in names(at)) row[[nm]] <- grid[g, nm]
    results <- rbind(results, row)
  }

  cat(sprintf("\n=== CONTRAST: %s (%s component) ===\n", var, component))
  cat(sprintf("  %s vs %s\n",
              paste(names(levels), levels, sep = "=", collapse = ", "),
              paste(names(ref),    ref,    sep = "=", collapse = ", ")))
  if (length(at) > 0) {
    cat(sprintf("  Evaluated at levels of: %s\n", paste(names(at), collapse = ", ")))
  }

  cat(sprintf("\n  %-12s %8s %18s %8s %18s %8s %8s\n",
              "", "delta", paste0(prob * 100, "% CI"),
              "exp(d)", paste0(prob * 100, "% CI"), "P(d>0)", "P(d<0)"))
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
  for (i in seq_len(nrow(results))) {
    at_str <- ""
    if (length(at) > 0) {
      at_vals <- sapply(names(at), function(nm) sprintf("%s=%s", nm, results[i, nm]))
      at_str  <- paste(at_vals, collapse = ", ")
    }
    cat(sprintf("  %-12s %8.3f [%6.3f, %6.3f]  %6.3f [%5.3f, %5.3f]  %.3f     %.3f\n",
                at_str,
                results$delta_median[i], results$delta_lo[i], results$delta_hi[i],
                results$exp_delta_median[i], results$exp_delta_lo[i], results$exp_delta_hi[i],
                results$prob_pos[i], results$prob_neg[i]))
  }

  invisible(results)
}

# =============================================================================
#' Posterior cross-correlation between clonogenic and kinetic coefficients
#' @param cfit A `cure_fit` object.
#' @export
test_correlation <- function(cfit) {
  fit <- cfit$fit; sp <- cfit$spec

  cat("\n=== TEST: Posterior cross-correlation (clonogenic vs kinetic) ===\n")
  cat("    |r| > 0.6 = caution | |r| > 0.8 = severe\n")
  cat("    r^2 is the fraction of posterior variance shared between the two\n")
  cat("    coefficients (overlap). |r| = 0.6 => 36% overlap;",
      "|r| = 0.8 => 64%.\n\n")

  inc_pars <- c("beta0",  paste0("beta[",  seq_len(sp$stan_data$K_inc), "]"))
  lat_pars <- c("gamma0", paste0("gamma[", seq_len(sp$stan_data$K_lat), "]"), "alpha")

  inc_labs <- c("Intercept (theta)",  sp$X_names)
  lat_labs <- c("Intercept (lambda)", sp$Z_names, "alpha")

  post_mat <- as.matrix(fit)[, c(inc_pars, lat_pars), drop = FALSE]
  cor_mat  <- stats::cor(post_mat)

  # Positional indexing avoids collisions when the same variable name
  # appears in both sub-models (e.g. ki67_percent in clonogenic and kinetic).
  n_inc <- length(inc_pars)
  n_lat <- length(lat_pars)
  cross <- cor_mat[seq_len(n_inc), n_inc + seq_len(n_lat), drop = FALSE]
  rownames(cross) <- inc_labs
  colnames(cross) <- lat_labs
  print(round(cross, 3))

  shared <- intersect(sp$X_names, sp$Z_names)
  if (length(shared) > 0) {
    cat("\n--- Shared variables ---\n")
    for (sv in shared) {
      ix <- which(sp$X_names == sv); iz <- which(sp$Z_names == sv)
      r <- as.numeric(stats::cor(
        as.matrix(fit)[, paste0("beta[",  ix, "]")],
        as.matrix(fit)[, paste0("gamma[", iz, "]")]))
      cat(sprintf("  %s (clono vs kinet): r = %.3f | r^2 = %.3f (%.0f%% shared)\n",
                  sv, r, r^2, 100 * r^2))
    }
  }

  problem <- which(abs(cross) > 0.6, arr.ind = TRUE)
  # Drop self-pairs (same row and column label), which can only arise from
  # the Intercept/alpha diagonals in pathological cases; cross is rectangular
  # so true self-overlap is impossible after the positional fix, but we keep
  # the guard to make the output robust to future refactors.
  if (nrow(problem) > 0) {
    keep <- inc_labs[problem[, 1]] != lat_labs[problem[, 2]]
    problem <- problem[keep, , drop = FALSE]
  }
  if (nrow(problem) > 0) {
    cat("\n!! PAIRS WITH |r| > 0.6:\n")
    for (k in seq_len(nrow(problem))) {
      r_val <- cross[problem[k, 1], problem[k, 2]]
      cat(sprintf("   %s <-> %s : r = %.3f | r^2 = %.3f (%.0f%% shared)%s\n",
                  inc_labs[problem[k, 1]], lat_labs[problem[k, 2]],
                  r_val, r_val^2, 100 * r_val^2,
                  ifelse(abs(r_val) > 0.8, " *** SEVERE ***", "")))
    }
  } else {
    cat("\n OK: No cross-correlations |r| > 0.6\n")
  }

  max_r <- max(abs(cross))
  verdict <- if (max_r > 0.8) "SEVERE" else if (max_r > 0.6) "CAUTION" else "OK"
  cat(sprintf("\nMax |cross-r|: %.3f | max r^2: %.3f (%.0f%% shared) -> %s\n",
              max_r, max_r^2, 100 * max_r^2, verdict))

  invisible(cross)
}

# =============================================================================
#' Incremental censoring sensitivity of the cure fraction
#' @param cfit      A `cure_fit` object.
#' @param quantiles Follow-up quantiles to refit at.
#' @param chains,iter MCMC settings used for the refits.
#' @export
test_censoring <- function(cfit, quantiles = c(1.0, 0.90, 0.80, 0.70),
                           chains = 2, iter = 2000) {
  sp <- cfit$spec; sd <- sp$stan_data

  cat("\n=== TEST: Incremental censoring ===\n")

  t_max <- max(sd$t)
  results <- list()

  for (q in quantiles) {
    t_cut <- t_max * q
    cat(sprintf("\n--- Follow-up at %.0f%% (t = %.1f) ---\n", q * 100, t_cut))

    t_c <- pmin(sd$t, t_cut)
    d_c <- as.integer(ifelse(sd$t <= t_cut, sd$d, 0L))
    cat(sprintf("  Events: %d -> %d\n", sum(sd$d), sum(d_c)))

    if (q == 1.0) {
      fc <- cfit$fit
    } else {
      sdc <- sd; sdc$t <- t_c; sdc$d <- d_c
      fc <- rstan::stan(file = sp$stan_file, data = sdc,
                        chains = chains, iter = iter,
                        warmup = floor(iter / 2),
                        seed = 42, refresh = 0,
                        control = list(adapt_delta = 0.95, max_treedepth = 12))
    }

    pc   <- as.data.frame(fc)
    cure <- exp(-exp(pc$beta0))

    results[[as.character(q)]] <- data.frame(
      quantile   = q, t_cut = t_cut, n_events = sum(d_c),
      cure_med   = median(cure),
      cure_lo    = quantile(cure, 0.025),
      cure_hi    = quantile(cure, 0.975),
      cure_width = quantile(cure, 0.975) - quantile(cure, 0.025),
      row.names  = NULL)

    cat(sprintf("  Cure: %.1f%% [%.1f%%, %.1f%%]\n",
                median(cure) * 100, quantile(cure, 0.025) * 100,
                quantile(cure, 0.975) * 100))
  }

  tab <- do.call(rbind, results); rownames(tab) <- NULL
  cat("\n=== SUMMARY ===\n")
  print(tab, digits = 3, row.names = FALSE)

  rng <- diff(range(tab$cure_med))
  wr  <- max(tab$cure_width) / min(tab$cure_width)
  cat(sprintf("\nCure range: %.1f pp | CI inflation: %.2fx\n", rng * 100, wr))
  if      (rng > 0.10)           cat("VERDICT: UNSTABLE\n")
  else if (rng > 0.05 || wr > 2) cat("VERDICT: CAUTION\n")
  else                           cat("VERDICT: STABLE\n")

  invisible(tab)
}

# =============================================================================
#' Approximate leave-one-out cross-validation (PSIS-LOO)
#' @param cfit A `cure_fit` object.
#' @export
loo_fit <- function(cfit) {
  ll <- loo::extract_log_lik(cfit$fit, parameter_name = "log_lik")
  result <- loo::loo(ll)
  print(result)
  invisible(result)
}

# =============================================================================
#' Triptych: clonogenic load, cure rate, and median failure time
#'   along a continuous covariate
#'
#' @param cfit         A `cure_fit` object.
#' @param x_var        Name of the continuous covariate to scan (must appear
#'   in the clonogenic and/or kinetic design).
#' @param fix          Named list of covariate values to hold fixed; unnamed
#'   covariates are held at their sample mean.
#' @param n_grid       Number of evaluation points on x_var.
#' @param save_path    Optional path (without extension) to save PDF+PNG.
#' @param xlab         X-axis label.
#' @export
plot_continuous_effect <- function(cfit, x_var = "ki67_percent",
                                   fix = list(), n_grid = 100,
                                   save_path = NULL, xlab = x_var) {

  sp   <- cfit$spec
  post <- as.data.frame(cfit$fit)

  if (!x_var %in% names(sp$dat))
    stop(sprintf("x_var '%s' is not present in the fitted data.", x_var))

  x_raw <- sp$dat[[x_var]]
  x_seq <- seq(min(x_raw, na.rm = TRUE), max(x_raw, na.rm = TRUE),
               length.out = n_grid)

  K_inc <- sp$stan_data$K_inc
  K_lat <- sp$stan_data$K_lat
  beta_cols  <- paste0("beta[",  seq_len(K_inc), "]")
  gamma_cols <- paste0("gamma[", seq_len(K_lat), "]")
  n_draws <- nrow(post)

  theta_mat <- cure_mat <- lambda_mat <- median_t_mat <-
    matrix(NA, n_draws, n_grid)

  for (i in seq_along(x_seq)) {
    raw_X <- setNames(sp$X_center, sp$X_names)
    raw_Z <- setNames(sp$Z_center, sp$Z_names)

    for (j in seq_along(sp$X_names)) if (sp$X_names[j] == x_var) raw_X[j] <- x_seq[i]
    for (j in seq_along(sp$Z_names)) if (sp$Z_names[j] == x_var) raw_Z[j] <- x_seq[i]

    for (nm in names(fix)) {
      for (j in seq_along(sp$X_names)) if (sp$X_names[j] == nm) raw_X[j] <- fix[[nm]]
      for (j in seq_along(sp$Z_names)) if (sp$Z_names[j] == nm) raw_Z[j] <- fix[[nm]]
    }

    for (j in seq_along(sp$X_names)) {
      if (grepl(":", sp$X_names[j], fixed = TRUE)) {
        parts <- strsplit(sp$X_names[j], ":")[[1]]
        v1 <- raw_X[parts[1]]; v2 <- raw_X[parts[2]]
        if (is.na(v1)) v1 <- sp$X_center[parts[1]]
        if (is.na(v2)) v2 <- sp$X_center[parts[2]]
        raw_X[j] <- v1 * v2
      }
    }
    for (j in seq_along(sp$Z_names)) {
      if (grepl(":", sp$Z_names[j], fixed = TRUE)) {
        parts <- strsplit(sp$Z_names[j], ":")[[1]]
        v1 <- raw_Z[parts[1]]; v2 <- raw_Z[parts[2]]
        if (is.na(v1)) v1 <- sp$Z_center[parts[1]]
        if (is.na(v2)) v2 <- sp$Z_center[parts[2]]
        raw_Z[j] <- v1 * v2
      }
    }

    x_std <- (raw_X - sp$X_center) / sp$X_scale
    z_std <- (raw_Z - sp$Z_center) / sp$Z_scale

    log_th <- post$beta0  + as.matrix(post[, beta_cols])  %*% x_std
    log_la <- post$gamma0 + as.matrix(post[, gamma_cols]) %*% z_std

    theta_mat[, i]    <- exp(log_th)
    lambda_mat[, i]   <- exp(log_la)
    cure_mat[, i]     <- exp(-theta_mat[, i])
    median_t_mat[, i] <- lambda_mat[, i] * (log(2))^(1 / post$alpha)
  }

  summ <- function(m) data.frame(
    med = apply(m, 2, median),
    lo  = apply(m, 2, quantile, 0.025),
    hi  = apply(m, 2, quantile, 0.975))

  df_plot <- data.frame(
    x         = x_seq,
    theta_med = summ(theta_mat)$med,  theta_lo = summ(theta_mat)$lo,  theta_hi = summ(theta_mat)$hi,
    cure_med  = summ(cure_mat)$med,   cure_lo  = summ(cure_mat)$lo,   cure_hi  = summ(cure_mat)$hi,
    mtime_med = summ(median_t_mat)$med, mtime_lo = summ(median_t_mat)$lo, mtime_hi = summ(median_t_mat)$hi
  )

  p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = theta_lo, ymax = theta_hi),
                         fill = "#3B82F6", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = theta_med),
                       colour = "#1E40AF", linewidth = 1.2) +
    ggplot2::labs(x = xlab, y = expression(theta),
                  title = sprintf("Clonogenic load vs %s", xlab)) +
    ggplot2::theme_minimal(base_size = 13)

  p2 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = cure_lo * 100, ymax = cure_hi * 100),
                         fill = "#10B981", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = cure_med * 100),
                       colour = "#047857", linewidth = 1.2) +
    ggplot2::labs(x = xlab, y = "P(cure) %",
                  title = sprintf("Cure rate vs %s", xlab)) +
    ggplot2::theme_minimal(base_size = 13) + ggplot2::ylim(0, NA)

  p3 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mtime_lo, ymax = mtime_hi),
                         fill = "#F59E0B", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = mtime_med),
                       colour = "#B45309", linewidth = 1.2) +
    ggplot2::labs(x = xlab, y = "Median failure time",
                  title = sprintf("Relapse speed vs %s", xlab)) +
    ggplot2::theme_minimal(base_size = 13)

  g <- gridExtra::arrangeGrob(p1, p2, p3, ncol = 3)
  gridExtra::grid.arrange(g)

  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".pdf"), g, width = 18, height = 6)
    ggplot2::ggsave(paste0(save_path, ".png"), g, width = 18, height = 6, dpi = 300)
    message(sprintf("Saved: %s.pdf/.png", save_path))
  }

  invisible(df_plot)
}

#' Back-compatibility alias
#' @inheritParams plot_continuous_effect
#' @param ki67_var deprecated, use `x_var`.
#' @export
plot_ki67_effect <- function(cfit, ki67_var = "ki67_percent",
                             fix = list(), n_grid = 100, save_path = NULL) {
  plot_continuous_effect(cfit, x_var = ki67_var, fix = fix,
                         n_grid = n_grid, save_path = save_path,
                         xlab = ki67_var)
}

# =============================================================================
#' Population survival curve implied by the cure model
#' @param cfit A `cure_fit` object.
#' @param t_grid Numeric grid of times at which to evaluate S(t).
#' @param main_title Plot title.
#' @export
plot_surv_curve <- function(cfit, t_grid = seq(0.01, 15, by = 0.1),
                            main_title = "Population survival (promotion time cure)") {
  sp <- cfit$spec; post <- as.data.frame(cfit$fit)
  K_inc <- sp$stan_data$K_inc; K_lat <- sp$stan_data$K_lat

  newX <- rep(0, K_inc); newZ <- rep(0, K_lat)
  beta_cols  <- paste0("beta[",  seq_len(K_inc), "]")
  gamma_cols <- paste0("gamma[", seq_len(K_lat), "]")
  n_draws <- nrow(post)
  S_mat <- matrix(NA, n_draws, length(t_grid))

  for (s in seq_len(n_draws)) {
    th <- exp(post$beta0[s]  + sum(newX * as.numeric(post[s, beta_cols])))
    la <- exp(post$gamma0[s] + sum(newZ * as.numeric(post[s, gamma_cols])))
    S0 <- exp(-(t_grid / la)^post$alpha[s])
    S_mat[s, ] <- exp(-th * (1 - S0))
  }

  S_med <- apply(S_mat, 2, median)
  S_lo  <- apply(S_mat, 2, quantile, 0.025)
  S_hi  <- apply(S_mat, 2, quantile, 0.975)

  opar <- par(mar = c(4, 4, 2.5, 1))
  on.exit(par(opar), add = TRUE)

  plot(t_grid, S_med, type = "l", lwd = 2, ylim = c(0, 1),
       xlab = "Time", ylab = "S(t)", main = main_title)
  polygon(c(t_grid, rev(t_grid)), c(S_lo, rev(S_hi)),
          col = rgb(0.2, 0.4, 0.8, 0.2), border = NA)
  cure <- exp(-exp(post$beta0 + as.matrix(post[, beta_cols]) %*% newX))
  abline(h = median(cure), lty = 2, col = "red")
  text(max(t_grid) * 0.75, median(cure) + 0.03,
       sprintf("Cure: %.1f%%", median(cure) * 100), col = "red")
}

# =============================================================================
#' Cure-probability odds ratio (exposed vs reference) along a modifier
#'
#' @param cfit       A `cure_fit` object.
#' @param var        Name of the variable that defines the contrast.
#' @param levels     Named list: exposed-scenario values.
#' @param ref        Named list: reference-scenario values.
#' @param at_var     Name of the continuous modifier plotted on the X axis.
#' @param at_seq     Numeric grid on `at_var` (default: observed range).
#' @param n_grid     Grid length when `at_seq = NULL`.
#' @param prob       Credible-interval width (default 0.95).
#' @param xlab,ylab,title Plot labels (English defaults).
#' @param save_path  Optional path without extension to save PDF+PNG.
#' @param log_scale  Log-scale for the OR axis (default TRUE).
#' @export
plot_contrast_cure_or <- function(cfit, var, levels, ref,
                                  at_var, at_seq = NULL, n_grid = 50,
                                  prob = 0.95,
                                  xlab  = NULL,
                                  ylab  = "Cure Odds Ratio",
                                  title = NULL, save_path = NULL,
                                  log_scale = TRUE) {

  sp   <- cfit$spec
  post <- as.data.frame(cfit$fit)

  var_names <- sp$X_names
  center    <- sp$X_center
  sc        <- sp$X_scale
  coef_cols <- paste0("beta[", seq_along(var_names), "]")

  if (is.null(at_seq)) {
    raw_vals <- sp$dat[[at_var]]
    at_seq <- seq(min(raw_vals, na.rm = TRUE),
                  max(raw_vals, na.rm = TRUE), length.out = n_grid)
  }

  alpha_lo <- (1 - prob) / 2
  alpha_hi <- 1 - alpha_lo
  n_draws  <- nrow(post)
  coef_mat <- as.matrix(post[, coef_cols])

  build_x <- function(scenario, at_val) {
    x <- rep(0, length(var_names)); names(x) <- var_names
    for (j in seq_along(var_names)) {
      if (var_names[j] == at_var) x[j] <- (at_val - center[j]) / sc[j]
    }
    for (nm in names(scenario)) {
      for (j in seq_along(var_names)) {
        if (var_names[j] == nm) x[j] <- (scenario[[nm]] - center[j]) / sc[j]
      }
    }
    for (j in seq_along(var_names)) {
      vn <- var_names[j]
      if (grepl(":", vn, fixed = TRUE)) {
        parts <- strsplit(vn, ":")[[1]]
        vals <- sapply(parts, function(p) {
          if (p %in% names(scenario)) return(scenario[[p]])
          if (p == at_var) return(at_val)
          center[p]
        })
        raw_int <- prod(vals)
        x[j] <- (raw_int - center[j]) / sc[j]
      }
    }
    x
  }

  or_mat       <- matrix(NA, n_draws, length(at_seq))
  cure_trt_mat <- matrix(NA, n_draws, length(at_seq))
  cure_ref_mat <- matrix(NA, n_draws, length(at_seq))

  for (i in seq_along(at_seq)) {
    x_trt <- build_x(levels, at_seq[i])
    x_ref <- build_x(ref,    at_seq[i])
    log_theta_trt <- post$beta0 + as.vector(coef_mat %*% x_trt)
    log_theta_ref <- post$beta0 + as.vector(coef_mat %*% x_ref)
    cure_trt <- exp(-exp(log_theta_trt))
    cure_ref <- exp(-exp(log_theta_ref))
    cure_trt_mat[, i] <- cure_trt
    cure_ref_mat[, i] <- cure_ref
    odds_trt <- cure_trt / (1 - cure_trt)
    odds_ref <- cure_ref / (1 - cure_ref)
    or_mat[, i] <- odds_trt / odds_ref
  }

  or_med <- apply(or_mat, 2, median)
  or_lo  <- apply(or_mat, 2, quantile, alpha_lo)
  or_hi  <- apply(or_mat, 2, quantile, alpha_hi)
  cure_trt_med <- apply(cure_trt_mat, 2, median) * 100
  cure_ref_med <- apply(cure_ref_mat, 2, median) * 100
  p_or_gt1     <- apply(or_mat, 2, function(x) mean(x > 1))

  df_plot <- data.frame(
    at     = at_seq,
    or_med = or_med, or_lo = or_lo, or_hi = or_hi,
    cure_trt = cure_trt_med, cure_ref = cure_ref_med,
    p_or_gt1 = p_or_gt1
  )

  if (is.null(xlab))  xlab  <- at_var
  if (is.null(title)) title <- sprintf("Cure OR: %s=%s vs %s=%s",
                                       var, paste(levels[[var]], collapse = ","),
                                       var, paste(ref[[var]],    collapse = ","))

  df_plot$sig <- ifelse(df_plot$or_lo > 1 | df_plot$or_hi < 1, "sig", "ns")

  p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = at, y = or_med)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = or_lo, ymax = or_hi),
                            colour = "#4A90D9", linewidth = 0.6, alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(fill = sig), shape = 21, size = 2.8,
                        colour = "white", stroke = 0.3) +
    ggplot2::scale_fill_manual(values = c("sig" = "#1B4F8A", "ns" = "#4A90D9"),
                               guide = "none") +
    ggplot2::labs(x = xlab, y = ylab, title = title,
                  subtitle = sprintf("%s%% credibility interval. Dark points: CI excludes OR = 1",
                                     prob * 100)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold", size = 13),
                   plot.subtitle = ggplot2::element_text(colour = "grey40", size = 10),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(size = 11))

  if (log_scale) p1 <- p1 + ggplot2::scale_y_log10()

  cure_trt_lo <- apply(cure_trt_mat, 2, quantile, alpha_lo) * 100
  cure_trt_hi <- apply(cure_trt_mat, 2, quantile, alpha_hi) * 100
  cure_ref_lo <- apply(cure_ref_mat, 2, quantile, alpha_lo) * 100
  cure_ref_hi <- apply(cure_ref_mat, 2, quantile, alpha_hi) * 100

  trt_label <- paste0(var, " = ", levels[[var]])
  ref_label <- paste0(var, " = ", ref[[var]])

  df_cure <- data.frame(
    at    = rep(at_seq, 2),
    cure  = c(cure_trt_med, cure_ref_med),
    lo    = c(cure_trt_lo,  cure_ref_lo),
    hi    = c(cure_trt_hi,  cure_ref_hi),
    group = factor(rep(c(trt_label, ref_label), each = length(at_seq)),
                   levels = c(trt_label, ref_label))
  )

  dodge_w <- diff(range(at_seq)) * 0.012

  p2 <- ggplot2::ggplot(df_cure, ggplot2::aes(x = at, y = cure, colour = group)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi),
                            linewidth = 0.5, alpha = 0.5,
                            position = ggplot2::position_dodge(width = dodge_w)) +
    ggplot2::geom_point(size = 2.2,
                        position = ggplot2::position_dodge(width = dodge_w)) +
    ggplot2::scale_colour_manual(values = c("#1B4F8A", "#B8372E")) +
    ggplot2::labs(x = xlab, y = "P(cure) %",
                  title = "Cure fraction by group", colour = "") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 10),
                   plot.title = ggplot2::element_text(face = "bold", size = 13),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(size = 11)) +
    ggplot2::ylim(0, NA)

  g <- gridExtra::arrangeGrob(p1, p2, ncol = 2)
  gridExtra::grid.arrange(g)

  if (!is.null(save_path)) {
    ggplot2::ggsave(paste0(save_path, ".pdf"), g, width = 14, height = 6)
    ggplot2::ggsave(paste0(save_path, ".png"), g, width = 14, height = 6, dpi = 300)
    message(sprintf("Saved: %s.pdf/.png", save_path))
  }

  cat("\n=== Cure OR summary ===\n")
  cat(sprintf("  %-8s %8s %18s %8s %10s %10s\n",
              at_var, "OR", sprintf("%s%% CI", prob * 100),
              "P(OR>1)", "Cure(trt)", "Cure(ref)"))
  cat("  ", paste(rep("-", 70), collapse = ""), "\n")
  for (i in seq_along(at_seq)) {
    cat(sprintf("  %-8.1f %8.2f [%5.2f, %5.2f]  %7.3f  %8.1f%%  %8.1f%%\n",
                at_seq[i], or_med[i], or_lo[i], or_hi[i],
                p_or_gt1[i], cure_trt_med[i], cure_ref_med[i]))
  }

  invisible(df_plot)
}

# =============================================================================
#' Kaplan-Meier survival curve (ggplot2, with stratification and faceting)
#'
#' Draws a Kaplan-Meier curve from a data frame containing a time and an
#' event-status column. Supports (i) stratification by one variable (colour
#' groups within each panel), (ii) faceting by a second variable (one panel
#' per level), and (iii) an optional subset filter so a panel can be built
#' for a single tumour site or stage. Log-rank p-values are added per panel
#' when at least two strata are present.
#'
#' @param dat A `data.frame`. Works with the classical preprocessed layout
#'   (`EFS_time`, `.status_bin`) or with arbitrary column names via
#'   `time_col` / `status_col`.
#' @param time_col,status_col Column names for follow-up time and event
#'   indicator. Defaults preserve backwards compatibility with [prep_data()].
#' @param strata   Optional character: name of a column to stratify the
#'   curves (colour legend within each panel).
#' @param facet_by Optional character: name of a column whose levels define
#'   the facet grid (one panel per level). May be the same as `strata`.
#' @param filter_var,filter_values Optional column name and set of values to
#'   restrict the dataset before plotting. Useful for "tumour site = Colon"
#'   type explorations.
#' @param time_max Maximum time to display on the x-axis (default: observed).
#' @param xlab,ylab,title Axis and title labels.
#' @param palette Optional character vector of colours for the strata.
#' @return A `ggplot` object (invisible in interactive sessions).
#' @export
plot_km <- function(dat,
                    strata        = NULL,
                    facet_by      = NULL,
                    filter_var    = NULL,
                    filter_values = NULL,
                    time_col      = NULL,
                    status_col    = NULL,
                    time_max      = NULL,
                    xlab          = "Time (years)",
                    ylab          = "Event-free survival",
                    title         = "Kaplan-Meier estimate",
                    palette       = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (is.null(time_col)) {
    time_col <- if ("EFS_time" %in% names(dat)) "EFS_time"
                else if ("EFS_Time" %in% names(dat)) "EFS_Time"
                else stop("Unable to find a time column; pass time_col.")
  }
  if (is.null(status_col)) {
    status_col <- if (".status_bin" %in% names(dat)) ".status_bin"
                  else if ("EFS_Status" %in% names(dat)) "EFS_Status"
                  else stop("Unable to find a status column; pass status_col.")
  }

  d <- dat
  if (!is.null(filter_var) && filter_var %in% names(d) &&
      length(filter_values) > 0) {
    d <- d[as.character(d[[filter_var]]) %in% as.character(filter_values),
           , drop = FALSE]
  }
  d <- d[!is.na(d[[time_col]]) & d[[time_col]] >= 0, , drop = FALSE]
  if (nrow(d) == 0) stop("No observations left after filtering.")

  d$..time..   <- as.numeric(d[[time_col]])
  d$..status.. <- as.integer(d[[status_col]])

  has_strata <- !is.null(strata) && strata %in% names(d)
  has_facet  <- !is.null(facet_by) && facet_by %in% names(d)

  build_curve <- function(df, group_col) {
    so <- survival::Surv(df$..time.., df$..status..)
    if (is.null(group_col)) {
      sf <- survival::survfit(so ~ 1)
      data.frame(
        time  = c(0, sf$time),
        surv  = c(1, sf$surv),
        cens  = c(FALSE, sf$n.censor > 0),
        group = "Overall",
        stringsAsFactors = FALSE
      )
    } else {
      g <- factor(df[[group_col]])
      sf <- survival::survfit(so ~ g)
      strata_lengths <- sf$strata
      grp_labels <- sub("^g=", "", names(strata_lengths))
      grp_vec <- rep(grp_labels, times = unname(strata_lengths))
      out <- data.frame(
        time  = sf$time,
        surv  = sf$surv,
        cens  = sf$n.censor > 0,
        group = grp_vec,
        stringsAsFactors = FALSE
      )
      # Prepend (0, 1) for each stratum so curves start at the top-left.
      starts <- do.call(rbind, lapply(split(out, out$group), function(z) {
        data.frame(time = 0, surv = 1, cens = FALSE,
                   group = z$group[1], stringsAsFactors = FALSE)
      }))
      rbind(starts, out)
    }
  }

  logrank_p <- function(df, group_col) {
    g <- df[[group_col]]
    if (length(unique(g[!is.na(g)])) < 2L) return(NA_real_)
    fit <- tryCatch(
      survival::survdiff(survival::Surv(df$..time.., df$..status..) ~ g),
      error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    dfr <- length(fit$n) - 1L
    stats::pchisq(fit$chisq, dfr, lower.tail = FALSE)
  }

  if (has_facet) {
    panels <- split(d, d[[facet_by]], drop = TRUE)
    curves <- do.call(rbind, lapply(names(panels), function(lv) {
      cur <- build_curve(panels[[lv]], if (has_strata) strata else NULL)
      cur$facet <- lv
      cur
    }))
    p_lbl_df <- data.frame(
      facet = names(panels),
      label = vapply(names(panels), function(lv) {
        if (!has_strata) return("")
        pv <- logrank_p(panels[[lv]], strata)
        if (is.na(pv)) "" else paste0("log-rank p ", .fmt_p(pv))
      }, character(1)),
      stringsAsFactors = FALSE
    )
  } else {
    curves <- build_curve(d, if (has_strata) strata else NULL)
    curves$facet <- "all"
    p_lbl_df <- data.frame(
      facet = "all",
      label = if (has_strata) {
        pv <- logrank_p(d, strata)
        if (is.na(pv)) "" else paste0("log-rank p ", .fmt_p(pv))
      } else "",
      stringsAsFactors = FALSE
    )
  }

  xlim_r <- if (!is.null(time_max) && time_max > 0) c(0, time_max)
            else c(0, max(curves$time, na.rm = TRUE))

  g <- ggplot2::ggplot(
        curves,
        ggplot2::aes(x = .data$time, y = .data$surv, color = .data$group)) +
    ggplot2::geom_step(linewidth = 0.9) +
    ggplot2::geom_point(
      data = curves[curves$cens, , drop = FALSE],
      shape = 3, size = 2, show.legend = FALSE) +
    ggplot2::coord_cartesian(xlim = xlim_r, ylim = c(0, 1)) +
    ggplot2::labs(x = xlab, y = ylab, title = title, color = strata) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = if (has_strata) "bottom" else "none")

  if (!is.null(palette)) {
    g <- g + ggplot2::scale_color_manual(values = palette)
  }
  if (has_facet) {
    g <- g + ggplot2::facet_wrap(~ facet)
  }
  if (any(nzchar(p_lbl_df$label))) {
    g <- g + ggplot2::geom_text(
      data = p_lbl_df[nzchar(p_lbl_df$label), , drop = FALSE],
      ggplot2::aes(x = xlim_r[1] + diff(xlim_r) * 0.02,
                   y = 0.08,
                   label = .data$label),
      inherit.aes = FALSE, hjust = 0, size = 3.3, color = "grey25")
  }

  g
}

# =============================================================================
#' Descriptive summary table of the analysis dataset
#'
#' Returns a `data.frame` with counts and percentages for binary columns and
#' mean/SD/median/IQR for continuous columns, suitable for rendering in Shiny
#' or printing to the console.
#'
#' @param dat A `data.frame` produced by [prep_data()].
#' @param vars Character vector of column names to summarize (default: all
#'   model variables plus `EFS_time` and `.status_bin`).
#' @return A `data.frame` with columns `variable`, `stat`, and `value`.
#' @export
data_summary_table <- function(dat, vars = NULL) {
  if (is.null(vars)) {
    vars <- intersect(
      c("EFS_time", ".status_bin",
        "primary_surgery_yes", "stage_II", "stage_III",
        "site_colorectal", "site_pancreas",
        "ki67_percent", "periop_therapy"),
      names(dat)
    )
  }

  rows <- list()
  for (v in vars) {
    x <- dat[[v]]
    x <- x[!is.na(x)]
    if (all(x %in% c(0, 1))) {
      n1 <- sum(x == 1)
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, stat = "n (%) = 1",
        value = sprintf("%d (%.1f%%)", n1, 100 * n1 / length(x)),
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, stat = "mean (SD)",
        value = sprintf("%.1f (%.1f)", mean(x), sd(x)),
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- data.frame(
        variable = v, stat = "median [IQR]",
        value = sprintf("%.1f [%.1f -- %.1f]",
                        median(x), quantile(x, 0.25), quantile(x, 0.75)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
