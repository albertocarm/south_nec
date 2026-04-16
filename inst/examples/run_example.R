# =============================================================================
# End-to-end example: Promotion Time Cure Model with an interaction
# =============================================================================
# This script reproduces the full workflow -- data, model, fit, diagnostics,
# contrasts, plots -- using only resources bundled with rgetne. Copy-paste,
# adjust priors/formulas, and go.
# =============================================================================

library(rgetne)

# --- 1. Data ---------------------------------------------------------------
raw <- readRDS(example_data())      # or: data(df_datos_app); raw <- df_datos_app
dat <- prep_data(raw)

# --- 2. Model specification + per-variable priors --------------------------
mod <- cure_model(
  dat,
  clonogenic = ~ primary_surgery_yes + stage_II + stage_III +
                 site_colorectal + ki67_percent * periop_therapy,
  kinetic    = ~ ki67_percent + periop_therapy,
  prior_clonogenic = 2.5,
  prior_kinetic    = 0.25,
  prior_intercept  = 2.5
)

# --- 3. Fit (NUTS) ---------------------------------------------------------
cfit <- fit_model(mod,
                  chains = 4, iter = 4000, warmup = 1000,
                  adapt_delta = 0.95)

# --- 4. MCMC diagnostics ---------------------------------------------------
mcmc_checks(cfit)

# --- 5. Coefficient summary on original scale ------------------------------
summary_coefs(cfit)

# --- 6. Leave-one-out CV ---------------------------------------------------
loo_fit(cfit)

# --- 7. Identifiability tests ----------------------------------------------
test_correlation(cfit)
test_censoring(cfit)

# --- 8. Contrasts ----------------------------------------------------------
cure_contrast(cfit, "clonogenic",
              var    = "periop_therapy",
              levels = list(periop_therapy = 1),
              ref    = list(periop_therapy = 0),
              at     = list(ki67_percent = c(30, 50, 70, 90)))

cure_contrast(cfit, "clonogenic",
              var    = "ki67_percent",
              levels = list(ki67_percent = 70),
              ref    = list(ki67_percent = 35),
              at     = list(periop_therapy = c(0, 1)))

# --- 9. Plots --------------------------------------------------------------
plot_continuous_effect(cfit, x_var = "ki67_percent",
                       fix = list(periop_therapy = 0))
plot_continuous_effect(cfit, x_var = "ki67_percent",
                       fix = list(periop_therapy = 1))

plot_contrast_cure_or(
  cfit,
  var    = "periop_therapy",
  levels = list(periop_therapy = 1),
  ref    = list(periop_therapy = 0),
  at_var = "ki67_percent",
  at_seq = seq(25, 95, by = 5),
  xlab   = "Ki67 (%)",
  title  = "Cure OR: peri-operative therapy vs none"
)

plot_surv_curve(cfit)
