# =============================================================================
# Multivariable Cox models, forest plots and combined gt tables
# =============================================================================
#
# These helpers complement the Bayesian cure-model workflow with a classical
# multivariable Cox analysis, producing publication-ready forest plots and a
# side-by-side gt table for two endpoints (typically Event-Free Survival and
# Overall Survival). The Shiny app exposes them as the "Cox models" tab.
#
# The API is dataset-agnostic: the caller supplies the time / status / variable
# columns by name, so the same functions work for `EFS_Time/EFS_Status` or
# `OS_Time/OS_Status` or any other survival pair in the input data.
# =============================================================================

# -----------------------------------------------------------------------------
# Small utilities (internal)
# -----------------------------------------------------------------------------

.fmt_hr <- function(eff, lo, hi) {
  paste0(sprintf("%.2f", eff),
         " (", sprintf("%.2f", lo), "\u2013", sprintf("%.2f", hi), ")")
}

.fmt_p <- function(p) {
  ifelse(is.na(p), "-",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

.pretty_label <- function(raw) {
  # "Var - Level:Reference"  ->  "Var (Level vs Reference)"
  m <- regmatches(raw, regexec("^(.*?) - (.*?):(.*)$", raw))[[1]]
  if (length(m) == 4L) {
    sprintf("%s (%s vs %s)", m[2], m[3], m[4])
  } else {
    raw
  }
}

# =============================================================================
#' Fit a multivariable Cox model and return a tidy forest-ready data frame
#'
#' Uses `rms::cph()` so that `summary()` produces already-transformed hazard
#' ratios with 95% confidence intervals, and Wald p-values computed from the
#' point estimate / standard error pair. Continuous covariates are contrasted
#' over their IQR by `rms::datadist`, matching the output of a standard
#' publication Cox analysis.
#'
#' @param dat A `data.frame` containing the survival columns and covariates.
#' @param time_col,status_col Character names of the time and event columns.
#' @param vars Character vector of covariate column names to include in the
#'   model.
#' @param ref_levels Optional named list `list(var = "reference_level")` that
#'   re-levels factors so the reference level is explicit. Any variable not
#'   present keeps its current reference level.
#' @return A list with:
#'   * `fit` — the fitted `cph` object;
#'   * `fp` — data frame (`Variable`, `Effect`, `CI_Lower`, `CI_Upper`,
#'     `p`, `p_lbl`, `hr_lbl`, `sig`) ready for `cox_forest_plot()` and
#'     `cox_combined_table()`;
#'   * `n`, `events` — observation and event counts used by the fit.
#' @export
cox_multivar_fit <- function(dat, time_col, status_col, vars,
                             ref_levels = NULL) {
  if (!requireNamespace("rms", quietly = TRUE)) {
    stop("Package 'rms' is required: install.packages('rms').")
  }
  stopifnot(is.data.frame(dat),
            time_col %in% names(dat),
            status_col %in% names(dat),
            all(vars %in% names(dat)))

  keep <- !is.na(dat[[time_col]]) & dat[[time_col]] >= 0
  d <- dat[keep, c(time_col, status_col, vars), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d <- droplevels(d)

  if (!is.null(ref_levels)) {
    for (nm in names(ref_levels)) {
      if (nm %in% names(d) && is.factor(d[[nm]]) &&
          ref_levels[[nm]] %in% levels(d[[nm]])) {
        d[[nm]] <- stats::relevel(d[[nm]], ref = ref_levels[[nm]])
      }
    }
  }

  dd <- rms::datadist(d[, vars, drop = FALSE])
  old_dd <- getOption("datadist")
  assign("..cox_dd..", dd, envir = .GlobalEnv)
  options(datadist = "..cox_dd..")
  on.exit({
    options(datadist = old_dd)
    if (exists("..cox_dd..", envir = .GlobalEnv)) {
      rm("..cox_dd..", envir = .GlobalEnv)
    }
  }, add = TRUE)

  rhs <- paste(vars, collapse = " + ")
  f <- stats::as.formula(
    sprintf("survival::Surv(%s, %s) ~ %s", time_col, status_col, rhs))

  fit <- rms::cph(f, data = d, x = TRUE, y = TRUE, surv = TRUE)
  s <- summary(fit)

  # rms::summary() reports each contrast as two rows: the raw beta difference
  # and its exponentiated HR. The HR rows are every second row starting from 2.
  idx_hr  <- seq(2, nrow(s), by = 2)
  idx_raw <- idx_hr - 1

  eff <- s[idx_hr,  "Effect"]
  lo  <- s[idx_hr,  "Lower 0.95"]
  hi  <- s[idx_hr,  "Upper 0.95"]
  se  <- s[idx_raw, "S.E."]
  bet <- s[idx_raw, "Effect"]
  p   <- 2 * (1 - stats::pnorm(abs(bet / se)))

  raw <- rownames(s)[idx_raw]
  labs <- vapply(raw, .pretty_label, character(1))

  fp <- data.frame(
    Variable = labs,
    Effect   = eff,
    CI_Lower = lo,
    CI_Upper = hi,
    p        = p,
    hr_lbl   = .fmt_hr(eff, lo, hi),
    p_lbl    = .fmt_p(p),
    sig      = p < 0.05,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(
    fit    = fit,
    fp     = fp,
    n      = unname(fit$stats["Obs"]),
    events = unname(fit$stats["Events"])
  )
}

# =============================================================================
#' Forest plot from a `cox_multivar_fit()` result
#'
#' @param fp Data frame from `cox_multivar_fit()$fp`.
#' @param title,subtitle Plot text.
#' @return A `ggplot` object.
#' @export
cox_forest_plot <- function(fp, title = "Forest plot",
                            subtitle = "Multivariable Cox model") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  jco <- tryCatch(ggsci::pal_jco()(2),
                  error = function(e) c("#0073C2FF", "#EFC000FF"))

  fp$sig_f <- factor(fp$sig, levels = c(FALSE, TRUE))
  fp$Variable <- factor(fp$Variable, levels = rev(fp$Variable))

  xmax_ci <- max(fp$CI_Upper, na.rm = TRUE)

  ggplot2::ggplot(fp, ggplot2::aes(x = .data$Effect, y = .data$Variable)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        color = "grey50") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$CI_Lower, xmax = .data$CI_Upper,
                   color = .data$sig_f),
      height = 0.25, linewidth = 1, show.legend = FALSE) +
    ggplot2::geom_point(
      ggplot2::aes(shape = .data$sig_f, color = .data$sig_f),
      size = 4.5, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(x = xmax_ci * 1.1, label = .data$p_lbl,
                   color = .data$sig_f),
      hjust = 0, size = 4.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = jco[2], "TRUE" = jco[1])) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 18)) +
    ggplot2::scale_x_continuous(
      "Hazard Ratio (95% CI)",
      expand = ggplot2::expansion(mult = c(0.05, 0.3))) +
    ggplot2::labs(title = title, subtitle = subtitle, y = NULL) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                            size = 17),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40",
                                            size = 13),
      axis.title.x  = ggplot2::element_text(size = 14),
      axis.text.x   = ggplot2::element_text(size = 12, color = "grey20"),
      axis.text.y   = ggplot2::element_text(size = 13, color = "grey10"))
}

# =============================================================================
#' Combined gt table for two Cox endpoints (e.g. EFS and OS)
#'
#' Produces a publication-ready gt table with two spanners, bolded p-values
#' for rows that reach statistical significance, alternating row shading and
#' a source note with counts.
#'
#' @param efs,os Results from `cox_multivar_fit()` (either list or the `fp`
#'   data frame; if a data frame is provided, `n_*` / `events_*` must be
#'   supplied explicitly).
#' @param n_efs,events_efs,n_os,events_os Integer counts. Only required when
#'   `efs` / `os` are data frames rather than full fit results.
#' @param endpoint_labels Character vector of length 2 with the column
#'   spanner labels.
#' @param extra_note Character appended to the source note.
#' @return A `gt` table object.
#' @export
cox_combined_table <- function(efs, os,
                               n_efs = NULL, events_efs = NULL,
                               n_os  = NULL, events_os  = NULL,
                               endpoint_labels = c("Event-Free Survival",
                                                   "Overall Survival"),
                               extra_note = NULL) {
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required: install.packages('gt').")
  }

  coerce <- function(x) {
    if (is.data.frame(x)) list(fp = x, n = NA_integer_, events = NA_integer_)
    else x
  }
  efs <- coerce(efs); os <- coerce(os)
  if (!is.null(n_efs))      efs$n      <- n_efs
  if (!is.null(events_efs)) efs$events <- events_efs
  if (!is.null(n_os))       os$n       <- n_os
  if (!is.null(events_os))  os$events  <- events_os

  de <- data.frame(
    Variable          = efs$fp$Variable,
    `EFS HR (95% CI)` = efs$fp$hr_lbl,
    `EFS P`           = efs$fp$p_lbl,
    efs_sig           = efs$fp$sig,
    check.names       = FALSE,
    stringsAsFactors  = FALSE
  )
  do <- data.frame(
    Variable         = os$fp$Variable,
    `OS HR (95% CI)` = os$fp$hr_lbl,
    `OS P`           = os$fp$p_lbl,
    os_sig           = os$fp$sig,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  comb <- merge(de, do, by = "Variable", all = TRUE, sort = FALSE)
  comb[is.na(comb$`EFS HR (95% CI)`), "EFS HR (95% CI)"] <- "-"
  comb[is.na(comb$`EFS P`),            "EFS P"]           <- "-"
  comb[is.na(comb$`OS HR (95% CI)`),  "OS HR (95% CI)"]  <- "-"
  comb[is.na(comb$`OS P`),             "OS P"]            <- "-"
  comb$efs_sig[is.na(comb$efs_sig)] <- FALSE
  comb$os_sig[is.na(comb$os_sig)]   <- FALSE

  note <- sprintf(
    "%s: N = %s, Events = %s. %s: N = %s, Events = %s. %s",
    endpoint_labels[1],
    ifelse(is.na(efs$n),      "-", as.character(efs$n)),
    ifelse(is.na(efs$events), "-", as.character(efs$events)),
    endpoint_labels[2],
    ifelse(is.na(os$n),       "-", as.character(os$n)),
    ifelse(is.na(os$events),  "-", as.character(os$events)),
    "HR > 1 indicates higher hazard. Bold indicates p < 0.05."
  )
  if (!is.null(extra_note) && nzchar(extra_note)) {
    note <- paste(note, extra_note)
  }

  tbl <- gt::gt(comb[, c("Variable",
                         "EFS HR (95% CI)", "EFS P",
                         "OS HR (95% CI)",  "OS P")])
  tbl <- gt::tab_header(
    tbl,
    title = gt::md("**Multivariable Cox Models**"),
    subtitle = gt::md(sprintf("*%s vs %s*",
                              endpoint_labels[1], endpoint_labels[2])))
  tbl <- gt::tab_spanner(tbl, label = endpoint_labels[1],
                         columns = c("EFS HR (95% CI)", "EFS P"))
  tbl <- gt::tab_spanner(tbl, label = endpoint_labels[2],
                         columns = c("OS HR (95% CI)",  "OS P"))
  tbl <- gt::cols_align(tbl, align = "center",
                        columns = c("EFS HR (95% CI)", "EFS P",
                                    "OS HR (95% CI)",  "OS P"))
  tbl <- gt::cols_label(tbl,
                        `EFS HR (95% CI)` = "HR (95% CI)",
                        `EFS P`           = "P value",
                        `OS HR (95% CI)`  = "HR (95% CI)",
                        `OS P`            = "P value")
  if (any(comb$efs_sig)) {
    tbl <- gt::tab_style(tbl,
                         style = gt::cell_text(weight = "bold"),
                         locations = gt::cells_body(columns = "EFS P",
                                                    rows = comb$efs_sig))
  }
  if (any(comb$os_sig)) {
    tbl <- gt::tab_style(tbl,
                         style = gt::cell_text(weight = "bold"),
                         locations = gt::cells_body(columns = "OS P",
                                                    rows = comb$os_sig))
  }
  if (nrow(comb) >= 2) {
    tbl <- gt::tab_style(
      tbl,
      style = gt::cell_fill(color = "#f0f4fb"),
      locations = gt::cells_body(rows = seq(1, nrow(comb), by = 2)))
  }
  gt::tab_source_note(tbl, source_note = note)
}

# =============================================================================
#' Convenience wrapper: fit EFS + OS and build both plots plus combined table
#'
#' @param dat,vars,ref_levels See `cox_multivar_fit()`.
#' @param efs_time,efs_status,os_time,os_status Column names for the two
#'   endpoints.
#' @return A list with `efs`, `os` (each a `cox_multivar_fit` result),
#'   `plot_efs`, `plot_os` (ggplot) and `table` (gt).
#' @export
cox_dual_endpoint <- function(dat, vars,
                              efs_time   = "EFS_Time",
                              efs_status = "EFS_Status",
                              os_time    = "OS_Time",
                              os_status  = "OS_Status",
                              ref_levels = NULL) {
  efs <- cox_multivar_fit(dat, efs_time, efs_status, vars, ref_levels)
  os  <- cox_multivar_fit(dat, os_time,  os_status,  vars, ref_levels)
  list(
    efs      = efs,
    os       = os,
    plot_efs = cox_forest_plot(efs$fp,
                               title = "Forest Plot \u2014 Event-Free Survival"),
    plot_os  = cox_forest_plot(os$fp,
                               title = "Forest Plot \u2014 Overall Survival"),
    table    = cox_combined_table(efs, os)
  )
}
