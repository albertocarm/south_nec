# =============================================================================
# Per-variable scaling, weakly-informative priors and marginal predictive
# effects for the Promotion Time Cure Model.
# =============================================================================
# This pipeline decouples predictor standardisation from model fitting so that
# per-variable priors can be specified on a common, interpretable scale. The
# response variable is never transformed. Means and scale divisors of every
# continuous predictor are stored and reused to project new data into the
# same space at prediction time.
# =============================================================================

# -----------------------------------------------------------------------------
#' Detect whether a vector is continuous
#'
#' A vector is treated as continuous when it is numeric and takes more than
#' two distinct non-missing values. Binary numeric vectors (0/1), logicals,
#' factors and character vectors are considered non-continuous.
#'
#' @param x A vector.
#' @return `TRUE` if `x` is continuous, `FALSE` otherwise.
#' @keywords internal
is_continuous <- function(x) {
  if (!is.numeric(x)) return(FALSE)
  length(unique(stats::na.omit(x))) > 2
}

# -----------------------------------------------------------------------------
#' Scale a continuous vector; leave binary or categorical vectors untouched
#'
#' Applies either a Gelman-style transformation (subtract the mean and divide
#' by twice the standard deviation) or a Z-score transformation (subtract the
#' mean and divide by one standard deviation) to continuous inputs. Binary
#' and categorical inputs are returned unchanged.
#'
#' The returned vector carries the attributes `"center"`, `"scale"` and
#' `"metodo"` describing the transformation, so callers can reproduce it on
#' new data without re-estimating the mean and SD.
#'
#' @param x A vector.
#' @param metodo One of `"gelman"` (default) or `"zscore"`.
#' @return A numeric vector of the same length as `x` with attributes
#'   `center`, `scale` and `metodo`. For non-continuous inputs the attributes
#'   are `NA` and the input is returned as-is.
#' @export
escalar_continua <- function(x, metodo = c("gelman", "zscore")) {
  metodo <- match.arg(metodo)

  if (!is_continuous(x)) {
    attr(x, "center") <- NA_real_
    attr(x, "scale")  <- NA_real_
    attr(x, "metodo") <- "none"
    return(x)
  }

  mu   <- mean(x, na.rm = TRUE)
  sd_x <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sd_x) || sd_x == 0) sd_x <- 1

  divisor <- if (metodo == "gelman") 2 * sd_x else sd_x
  z       <- (x - mu) / divisor

  attr(z, "center") <- mu
  attr(z, "scale")  <- divisor
  attr(z, "metodo") <- metodo
  z
}

# -----------------------------------------------------------------------------
#' Pre-process predictor columns, leaving the response untouched
#'
#' Applies [escalar_continua()] to every column listed in `predictors`. The
#' response column (and any other column not in `predictors`) is never
#' touched. The original mean and scale divisor of every continuous
#' predictor are stored so that [posterior_epred_cure()] can project new
#' data into the same standardised space.
#'
#' @param data A `data.frame` containing the predictors and the response.
#' @param predictors Character vector of predictor column names. The
#'   response must NOT be included.
#' @param metodo Scaling method, passed to [escalar_continua()].
#' @return A list with components
#'   * `data` - the input `data.frame` with predictor columns replaced by
#'     their scaled versions (other columns untouched),
#'   * `center` - named vector with the mean of each continuous predictor,
#'   * `scale`  - named vector with the scale divisor of each continuous
#'     predictor,
#'   * `metodo` - the scaling method used,
#'   * `continuous` - names of columns that were actually scaled,
#'   * `predictors` - the input `predictors` argument.
#' @export
preprocess_predictors <- function(data, predictors,
                                  metodo = c("gelman", "zscore")) {
  metodo <- match.arg(metodo)
  stopifnot(is.data.frame(data))

  missing_cols <- setdiff(predictors, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Predictors not found in data: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  center <- stats::setNames(rep(NA_real_, length(predictors)), predictors)
  divs   <- stats::setNames(rep(NA_real_, length(predictors)), predictors)
  continuous <- character(0)

  for (p in predictors) {
    z <- escalar_continua(data[[p]], metodo = metodo)
    c_p <- attr(z, "center")
    s_p <- attr(z, "scale")
    data[[p]]  <- as.numeric(z)
    center[p]  <- c_p
    divs[p]    <- s_p
    if (!is.na(c_p)) continuous <- c(continuous, p)
  }

  list(data = data,
       center = center, scale = divs,
       metodo = metodo,
       continuous = continuous,
       predictors = predictors)
}

# -----------------------------------------------------------------------------
#' Default weakly-informative per-variable prior SDs
#'
#' Builds a named vector of prior SDs covering every column produced by
#' `model.matrix(formula, data)`. Three tiers are applied in order:
#'
#' 1. every coefficient starts at `sd_binary` (binary / categorical default),
#' 2. coefficients whose column name matches a scaled continuous predictor
#'    are tightened to `sd_continuous`,
#' 3. coefficients whose column name matches any variable in `tight_vars`
#'    are tightened further to `sd_tight` (including interaction columns
#'    that contain the variable name).
#'
#' The default `tight_vars = c("ki67_percent", "periop_therapy")` and
#' `sd_tight = 0.25` mirror the hand-tuned priors used in the GETNE /
#' SOUTH-NEC analysis, where those two covariates drive the clonogenic /
#' kinetic posterior ridge. Pass `tight_vars = character(0)` to disable
#' the third tier.
#'
#' @param design_names Character vector with the column names of the design
#'   matrix.
#' @param continuous Character vector of names of scaled continuous
#'   predictors (typically `preprocess_predictors()$continuous`).
#' @param sd_continuous Prior SD for coefficients of standardised continuous
#'   predictors. Default `1` is weakly informative on the Gelman scale.
#' @param sd_binary Prior SD for coefficients of binary / categorical
#'   predictors. Default `2.5` matches the common `rstanarm` convention.
#' @param sd_tight Prior SD applied to every coefficient that matches a
#'   name in `tight_vars`. Default `0.25`.
#' @param tight_vars Character vector of variable names that should receive
#'   the tight prior. Default `c("ki67_percent", "periop_therapy")`.
#' @return Named numeric vector of length `length(design_names)`.
#' @keywords internal
default_prior_sds <- function(design_names,
                              continuous,
                              sd_continuous = 1.0,
                              sd_binary     = 2.5,
                              sd_tight      = 0.25,
                              tight_vars    = c("ki67_percent",
                                                "periop_therapy")) {
  sds <- rep(sd_binary, length(design_names))
  names(sds) <- design_names
  for (nm in continuous) {
    hit <- grepl(nm, design_names, fixed = TRUE)
    sds[hit] <- sd_continuous
  }
  for (nm in tight_vars) {
    hit <- grepl(nm, design_names, fixed = TRUE)
    sds[hit] <- sd_tight
  }
  sds
}

# -----------------------------------------------------------------------------
#' Fit a Promotion Time Cure Model with explicit per-variable scaling
#'
#' End-to-end wrapper that
#' (1) standardises continuous predictors via [preprocess_predictors()],
#' (2) builds a `cure_spec` with `scale = FALSE` so the prior SDs are
#' applied directly on the already-standardised design, and
#' (3) fits the Stan model. Binary and categorical predictors are not
#' rescaled; their coefficients therefore retain their natural
#' one-unit-change interpretation.
#'
#' The returned object stores the scaling metadata under `spec$scaling`,
#' which is required by [posterior_epred_cure()] to project new data.
#'
#' @param data A `data.frame` containing the response (`EFS_time`,
#'   `.status_bin`) and every column used in `clonogenic` and `kinetic`.
#' @param clonogenic,kinetic One-sided model formulas.
#' @param predictors Character vector of predictor names. Defaults to all
#'   variables appearing on the right-hand side of either formula.
#' @param metodo Scaling method for continuous predictors
#'   (`"gelman"` or `"zscore"`). See [escalar_continua()].
#' @param prior_clonogenic,prior_kinetic,prior_intercept Prior
#'   specifications forwarded to [cure_model()]. When `NULL` (default), a
#'   weakly-informative per-variable prior is built via
#'   [default_prior_sds()]: Normal(0, 1) for standardised continuous
#'   coefficients and Normal(0, 2.5) for binary / categorical ones. Users
#'   can still pass a scalar, a positional vector or a named list / named
#'   vector to override individual coefficients.
#' @param ... Passed to [fit_model()] (e.g. `chains`, `iter`, `seed`).
#' @return A `cure_fit` object whose `spec$scaling` slot carries the
#'   preprocessing metadata.
#' @export
fit_cure_bayes <- function(data,
                           clonogenic,
                           kinetic,
                           predictors       = NULL,
                           metodo           = c("gelman", "zscore"),
                           prior_clonogenic = NULL,
                           prior_kinetic    = NULL,
                           prior_intercept  = 2.5,
                           ...) {

  metodo <- match.arg(metodo)

  if (is.null(predictors)) {
    predictors <- unique(c(all.vars(clonogenic), all.vars(kinetic)))
  }

  pp <- preprocess_predictors(data, predictors, metodo = metodo)

  dat_cc <- pp$data[stats::complete.cases(
    pp$data[, unique(c("EFS_time", ".status_bin", predictors)), drop = FALSE]
  ), , drop = FALSE]

  X_names <- colnames(stats::model.matrix(clonogenic, data = dat_cc))[-1]
  Z_names <- colnames(stats::model.matrix(kinetic,    data = dat_cc))[-1]

  if (is.null(prior_clonogenic)) {
    prior_clonogenic <- default_prior_sds(X_names, pp$continuous)
  }
  if (is.null(prior_kinetic)) {
    prior_kinetic <- default_prior_sds(Z_names, pp$continuous)
  }

  mod <- cure_model(dat = pp$data,
                    clonogenic       = clonogenic,
                    kinetic          = kinetic,
                    prior_clonogenic = prior_clonogenic,
                    prior_kinetic    = prior_kinetic,
                    prior_intercept  = prior_intercept,
                    scale            = FALSE)

  mod$scaling <- pp

  cfit <- fit_model(mod, ...)
  cfit$spec$scaling <- pp
  cfit
}

# -----------------------------------------------------------------------------
#' Posterior expected predictions for a promotion time cure fit
#'
#' Projects a `newdata` `data.frame` (on the ORIGINAL, un-standardised scale
#' of the predictors) into the posterior space of the model and returns one
#' of the following quantities as a `draws x rows` matrix:
#' * `"cure_prob"`  - individual cure probability `exp(-theta)` (default),
#' * `"theta"`      - expected clonogen count `exp(log_theta)`,
#' * `"lambda"`     - latency scale `exp(log_lambda)`,
#' * `"median_time"`- implied median event time `lambda * log(2)^(1/alpha)`.
#'
#' Continuous columns of `newdata` are scaled with the same formula
#' (Gelman or Z-score) and the same mean / divisor estimated on the
#' training data and stored in `cfit$spec$scaling`. Binary and categorical
#' columns are passed through unchanged.
#'
#' @param cfit A `cure_fit` object produced by [fit_cure_bayes()].
#' @param newdata A `data.frame` with the original-scale predictors used
#'   in `cfit$spec$formula_clono` and `cfit$spec$formula_kinet`.
#' @param type One of `"cure_prob"`, `"theta"`, `"lambda"`, `"median_time"`.
#' @return A numeric matrix with `n_draws` rows and `nrow(newdata)` columns,
#'   carrying attributes `type` and `newdata` for downstream summarisation.
#' @export
posterior_epred_cure <- function(cfit, newdata,
                                 type = c("cure_prob", "theta",
                                          "lambda",    "median_time")) {
  type <- match.arg(type)
  if (!inherits(cfit, "cure_fit"))
    stop("'cfit' must be an object returned by fit_cure_bayes() / fit_model().")

  sp <- cfit$spec
  pp <- sp$scaling
  if (is.null(pp))
    stop("Model spec has no scaling metadata; refit with fit_cure_bayes().")

  nd <- newdata
  for (v in pp$continuous) {
    if (v %in% names(nd)) {
      nd[[v]] <- (nd[[v]] - pp$center[v]) / pp$scale[v]
    }
  }

  X_new <- stats::model.matrix(sp$formula_clono, data = nd)[, -1, drop = FALSE]
  Z_new <- stats::model.matrix(sp$formula_kinet, data = nd)[, -1, drop = FALSE]

  miss_X <- setdiff(sp$X_names, colnames(X_new))
  miss_Z <- setdiff(sp$Z_names, colnames(Z_new))
  if (length(miss_X) > 0 || length(miss_Z) > 0) {
    stop(sprintf("newdata does not reproduce the training design. Missing: %s",
                 paste(c(miss_X, miss_Z), collapse = ", ")))
  }
  X_new <- X_new[, sp$X_names, drop = FALSE]
  Z_new <- Z_new[, sp$Z_names, drop = FALSE]

  post <- as.data.frame(cfit$fit)
  beta_cols  <- paste0("beta[",  seq_along(sp$X_names), "]")
  gamma_cols <- paste0("gamma[", seq_along(sp$Z_names), "]")

  Beta  <- as.matrix(post[, beta_cols,  drop = FALSE])
  Gamma <- as.matrix(post[, gamma_cols, drop = FALSE])
  beta0  <- post$beta0
  gamma0 <- post$gamma0
  alpha  <- post$alpha

  log_theta  <- sweep(Beta  %*% t(X_new), 1, beta0,  "+")
  log_lambda <- sweep(Gamma %*% t(Z_new), 1, gamma0, "+")

  theta  <- exp(log_theta)
  lambda <- exp(log_lambda)

  out <- switch(type,
                cure_prob   = exp(-theta),
                theta       = theta,
                lambda      = lambda,
                median_time = lambda * (log(2))^(1 / alpha))

  attr(out, "type")    <- type
  attr(out, "newdata") <- newdata
  out
}

# -----------------------------------------------------------------------------
#' Summarise a posterior expected-prediction matrix
#'
#' Collapses the `draws x rows` matrix returned by [posterior_epred_cure()]
#' into a tidy `data.frame` with posterior mean, median and credible
#' interval for every row of the original `newdata`. The original-scale
#' `newdata` is recycled from the matrix attributes when available so the
#' caller can plot the expected prediction against the clinical scale of
#' the covariates directly.
#'
#' @param epred Matrix returned by [posterior_epred_cure()].
#' @param newdata Optional `data.frame` used for the predictions. Defaults
#'   to `attr(epred, "newdata")`.
#' @param probs Two-element vector with the lower and upper probability
#'   defining the credible interval. Default `c(0.025, 0.975)`.
#' @return A `data.frame` with columns `median`, `mean`, `lo`, `hi` and
#'   the columns of `newdata`.
#' @export
summarize_epred <- function(epred,
                            newdata = NULL,
                            probs   = c(0.025, 0.975)) {
  if (is.null(newdata)) newdata <- attr(epred, "newdata")
  stopifnot(length(probs) == 2, probs[1] < probs[2])

  qmat <- apply(epred, 2, stats::quantile, probs = probs, na.rm = TRUE)

  out <- data.frame(
    median = apply(epred, 2, stats::median, na.rm = TRUE),
    mean   = colMeans(epred, na.rm = TRUE),
    lo     = qmat[1, ],
    hi     = qmat[2, ]
  )
  type <- attr(epred, "type")
  if (!is.null(type)) out$type <- type
  if (!is.null(newdata)) out <- cbind(newdata, out, row.names = NULL)
  out
}
