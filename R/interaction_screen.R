# =============================================================================
# Systematic exploration of pairwise interactions in the clonogenic submodel.
# Reproduces the structure of Supplementary Table A5: each interaction is fitted
# as a separate model (additive base + single interaction term) and compared
# against the additive main-effects model via LOO-CV.
# =============================================================================

# -----------------------------------------------------------------------------
#' Map a probability of direction to a qualitative evidence label
#'
#' Implements the convention used in Supplementary Table A5:
#' `**` for `pd > 0.90`, `*` for `pd > 0.80`, `Weak` for `pd > 0.60`,
#' otherwise `None`.
#'
#' @param pd Probability of direction in `[0.5, 1]`.
#' @return A character label.
#' @keywords internal
evidence_label <- function(pd) {
  if (is.na(pd))         return(NA_character_)
  if (pd >  0.90)        return("Suggestive **")
  if (pd >  0.80)        return("Suggestive *")
  if (pd >  0.60)        return("Weak")
  "None"
}

# -----------------------------------------------------------------------------
#' Resolve a coefficient column for an interaction term in a fitted spec
#'
#' Returns the index `j` such that `paste0("beta[", j, "]")` is the posterior
#' draws column for the interaction. Tries both orderings (`a:b` and `b:a`)
#' because `model.matrix()` may reorder factors.
#'
#' @keywords internal
.find_interaction_index <- function(spec, pair) {
  cand <- c(paste(pair,        collapse = ":"),
            paste(rev(pair),   collapse = ":"))
  hit  <- which(spec$X_names %in% cand)
  if (length(hit) == 0) {
    stop(sprintf("Could not locate interaction term '%s' in fitted design (X_names: %s)",
                 cand[1], paste(spec$X_names, collapse = ", ")))
  }
  hit[1]
}

# -----------------------------------------------------------------------------
#' Systematically screen pairwise interactions in the clonogenic submodel
#'
#' Fits the additive base model and one model per supplied pair (base + single
#' interaction term). For every interaction model, the function reports:
#' the posterior median and 95% credible interval of the interaction
#' coefficient on the original (un-standardised) covariate scale,
#' the probability of direction (posterior probability that the coefficient
#' has the sign of its median), and the LOO-CV ELPD difference relative to
#' the additive base model. Results are returned sorted by strength of
#' evidence so the table can be rendered directly as Supplementary Table A5.
#'
#' Conventions
#' \itemize{
#'   \item `pd > 0.90` is reported as "Suggestive **", `pd > 0.80` as
#'         "Suggestive *", `pd > 0.60` as "Weak", otherwise "None".
#'   \item ELPD differences are reported as
#'         `delta_elpd = elpd_loo(interaction) - elpd_loo(base)`.
#'         Negative values favour the additive (base) model.
#'   \item For continuous \eqn{\times} binary interactions the coefficient
#'         is the change in the slope of the continuous variable when the
#'         binary variable is switched on. For binary \eqn{\times} binary
#'         interactions it is the departure from additivity on the
#'         log-\eqn{\theta} scale.
#' }
#'
#' @param data A `data.frame` already pre-processed by [prep_data()].
#' @param base_clonogenic A one-sided formula for the additive clonogenic
#'   linear predictor, e.g.
#'   `~ ki67_percent + periop_therapy + primary_surgery_yes + stage_III + site_pancreas`.
#' @param kinetic A one-sided formula for the latency linear predictor.
#' @param pairs A list of length-2 character vectors naming the variables to
#'   interact one at a time. Defaults to all unordered pairs of the variables
#'   appearing in `base_clonogenic`.
#' @param labels Optional named character vector mapping `"a:b"` term names to
#'   pretty display labels (e.g. \code{c("ki67_percent:periop_therapy" =
#'   "Ki-67 x Periop. chemotherapy")}). Both orderings are recognised.
#' @param prior_clonogenic,prior_kinetic,prior_intercept Forwarded to
#'   [fit_cure_bayes()]. Defaults to `NULL` so weakly-informative per-variable
#'   priors are built automatically.
#' @param chains,iter,warmup,seed,adapt_delta,max_treedepth Forwarded to
#'   [fit_model()] via [fit_cure_bayes()].
#' @param prob Credible interval width. Default `0.95`.
#' @param meaningful_se Multiplier of `se_diff` above which an ELPD difference
#'   is flagged as meaningful. Default `2`, matching the table footnote.
#' @param verbose If `TRUE`, prints progress for every fitted model.
#' @param ... Additional arguments forwarded to [fit_cure_bayes()].
#'
#' @return A `data.frame` of class `c("interaction_screen", "data.frame")`
#'   with one row per supplied pair, sorted by strength of evidence and
#'   then by `pd` (descending). Columns:
#'   \describe{
#'     \item{term}{Pretty label for the interaction.}
#'     \item{term_raw}{`"a:b"` interaction column name in the design matrix.}
#'     \item{coef}{Posterior median of the interaction coefficient
#'                 on the original covariate scale.}
#'     \item{lo, hi}{Lower and upper bounds of the `prob` credible interval.}
#'     \item{pd}{Probability of direction.}
#'     \item{p_neg}{Posterior probability that the coefficient is negative.}
#'     \item{delta_elpd}{`elpd_loo(interaction) - elpd_loo(base)`.}
#'     \item{se_diff}{Standard error of `delta_elpd`.}
#'     \item{meaningful}{`TRUE` if `|delta_elpd| > meaningful_se * se_diff`
#'                       AND favouring the interaction model.}
#'     \item{evidence}{Qualitative label (Suggestive **, Suggestive *,
#'                     Weak, None).}
#'   }
#'   The returned object carries `attr(., "loo_base")`, `attr(., "fit_base")`
#'   and `attr(., "fits_int")` so the underlying objects can be inspected.
#'
#' @export
screen_interactions <- function(data,
                                base_clonogenic,
                                kinetic,
                                pairs            = NULL,
                                labels           = NULL,
                                prior_clonogenic = NULL,
                                prior_kinetic    = NULL,
                                prior_intercept  = 2.5,
                                chains           = 4,
                                iter             = 4000,
                                warmup           = floor(iter / 2),
                                seed             = 42,
                                adapt_delta      = 0.95,
                                max_treedepth    = 12,
                                prob             = 0.95,
                                meaningful_se    = 2,
                                verbose          = TRUE,
                                ...) {

  stopifnot(inherits(base_clonogenic, "formula"),
            inherits(kinetic,         "formula"))

  base_vars <- all.vars(base_clonogenic)
  if (length(base_vars) < 2)
    stop("`base_clonogenic` must contain at least two variables.")

  if (is.null(pairs)) {
    pairs <- utils::combn(base_vars, 2, simplify = FALSE)
  } else {
    pairs <- lapply(pairs, function(p) {
      if (length(p) != 2) stop("Each entry of `pairs` must be a length-2 character vector.")
      missing <- setdiff(p, base_vars)
      if (length(missing) > 0)
        stop(sprintf("Variables not in base_clonogenic: %s",
                     paste(missing, collapse = ", ")))
      as.character(p)
    })
  }

  if (verbose)
    message(sprintf("Fitting additive base model (K = %d covariates)...",
                    length(base_vars)))

  fit_base <- fit_cure_bayes(
    data             = data,
    clonogenic       = base_clonogenic,
    kinetic          = kinetic,
    prior_clonogenic = prior_clonogenic,
    prior_kinetic    = prior_kinetic,
    prior_intercept  = prior_intercept,
    chains           = chains,
    iter             = iter,
    warmup           = warmup,
    seed             = seed,
    adapt_delta      = adapt_delta,
    max_treedepth    = max_treedepth,
    refresh          = 0,
    ...
  )

  ll_base  <- loo::extract_log_lik(fit_base$fit, parameter_name = "log_lik")
  loo_base <- loo::loo(ll_base)

  alpha_lo <- (1 - prob) / 2
  alpha_hi <- 1 - alpha_lo

  rows     <- vector("list", length(pairs))
  fits_int <- vector("list", length(pairs))

  for (k in seq_along(pairs)) {
    pair      <- pairs[[k]]
    int_term  <- paste(pair, collapse = ":")
    f_int     <- stats::update(base_clonogenic,
                               stats::as.formula(paste("~ . +", int_term)))

    if (verbose)
      message(sprintf("[%d/%d] Fitting %s ...", k, length(pairs), int_term))

    fit_int <- fit_cure_bayes(
      data             = data,
      clonogenic       = f_int,
      kinetic          = kinetic,
      prior_clonogenic = prior_clonogenic,
      prior_kinetic    = prior_kinetic,
      prior_intercept  = prior_intercept,
      chains           = chains,
      iter             = iter,
      warmup           = warmup,
      seed             = seed,
      adapt_delta      = adapt_delta,
      max_treedepth    = max_treedepth,
      refresh          = 0,
      ...
    )
    fits_int[[k]] <- fit_int

    sp       <- fit_int$spec
    j        <- .find_interaction_index(sp, pair)
    post     <- as.data.frame(fit_int$fit)
    coef_col <- paste0("beta[", j, "]")

    # Recover the coefficient on the original covariate scale. Each component
    # of the interaction column was either left untouched (binary / factor) or
    # divided by its scale attribute by `preprocess_predictors()`. Rescaling
    # the coefficient by the product of those divisors undoes the standardising
    # transformation for the interaction term.
    pp        <- sp$scaling
    divisor   <- 1
    for (v in pair) {
      sv <- if (!is.null(pp$scale) && v %in% names(pp$scale)) pp$scale[[v]] else NA_real_
      if (!is.na(sv) && sv != 0) divisor <- divisor * sv
    }
    coef_orig <- post[[coef_col]] / divisor

    med   <- stats::median(coef_orig)
    lo    <- unname(stats::quantile(coef_orig, alpha_lo))
    hi    <- unname(stats::quantile(coef_orig, alpha_hi))
    p_neg <- mean(coef_orig < 0)
    pd    <- if (med >= 0) mean(coef_orig > 0) else p_neg

    ll_int  <- loo::extract_log_lik(fit_int$fit, parameter_name = "log_lik")
    loo_int <- loo::loo(ll_int)

    elpd_int  <- loo_int$estimates ["elpd_loo", "Estimate"]
    elpd_base <- loo_base$estimates["elpd_loo", "Estimate"]
    delta     <- elpd_int - elpd_base

    cmp     <- loo::loo_compare(loo_base, loo_int)
    se_diff <- if (nrow(cmp) >= 2) cmp[2, "se_diff"] else NA_real_

    meaningful <- isTRUE(delta > meaningful_se * se_diff)

    label <- int_term
    if (!is.null(labels)) {
      cands <- c(int_term, paste(rev(pair), collapse = ":"))
      hit   <- intersect(cands, names(labels))
      if (length(hit) > 0) label <- unname(labels[hit[1]])
    }

    rows[[k]] <- data.frame(
      term       = label,
      term_raw   = int_term,
      coef       = med,
      lo         = lo,
      hi         = hi,
      pd         = pd,
      p_neg      = p_neg,
      delta_elpd = delta,
      se_diff    = se_diff,
      meaningful = meaningful,
      evidence   = evidence_label(pd),
      row.names  = NULL,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)

  ev_order <- c("Suggestive **" = 1L,
                "Suggestive *"  = 2L,
                "Weak"          = 3L,
                "None"          = 4L)
  out <- out[order(ev_order[out$evidence], -out$pd), , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "loo_base") <- loo_base
  attr(out, "fit_base") <- fit_base
  attr(out, "fits_int") <- stats::setNames(fits_int,
                                           vapply(pairs, paste, character(1),
                                                  collapse = ":"))
  class(out) <- c("interaction_screen", class(out))
  out
}

# -----------------------------------------------------------------------------
#' Format an `interaction_screen` result as a publication-style table
#'
#' Returns a `data.frame` with the same column layout as Supplementary
#' Table A5: \code{Interaction term}, \code{Coefficient}, \code{95\% CrI},
#' \code{P(direction)}, \eqn{\Delta}\code{ELPD}, \code{Evidence}.
#'
#' @param x An object returned by [screen_interactions()].
#' @param digits_coef Decimals for the coefficient and CrI bounds.
#' @param digits_elpd Decimals for the ELPD difference.
#' @return A `data.frame` ready to pass to `knitr::kable()`.
#' @export
format_interaction_table <- function(x,
                                     digits_coef = 3,
                                     digits_elpd = 1) {
  if (!inherits(x, "interaction_screen"))
    stop("`x` must be an object returned by screen_interactions().")

  fmt_pct <- function(p) sprintf("%.1f%%", 100 * p)

  out <- data.frame(
    interaction_term = x$term,
    coefficient      = formatC(x$coef, format = "f", digits = digits_coef),
    cri_95           = sprintf("%s, %s",
                               formatC(x$lo, format = "f", digits = digits_coef),
                               formatC(x$hi, format = "f", digits = digits_coef)),
    p_direction      = fmt_pct(x$pd),
    delta_elpd       = formatC(x$delta_elpd, format = "f", digits = digits_elpd),
    evidence         = x$evidence,
    stringsAsFactors = FALSE
  )
  names(out) <- c("Interaction term", "Coefficient", "95% CrI",
                  "P(direction)", "\u0394ELPD", "Evidence")
  out
}

# -----------------------------------------------------------------------------
#' Print method for `interaction_screen` objects
#' @param x An object returned by [screen_interactions()].
#' @param ... Ignored.
#' @export
print.interaction_screen <- function(x, ...) {
  cat("\n=== Pairwise interaction screen (clonogenic submodel) ===\n")
  cat("ELPD differences are: elpd(interaction) - elpd(additive base).\n")
  cat("Negative values favour the additive model.\n\n")
  print(format_interaction_table(x), row.names = FALSE)
  invisible(x)
}
