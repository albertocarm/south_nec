# =============================================================================
# Example: Supplementary Table A5
# Systematic exploration of pairwise interactions in the clonogenic submodel
# =============================================================================
# Each interaction is fitted as a separate model (additive base + single
# interaction term) and compared against the additive main-effects model via
# LOO-CV. Posterior medians, 95% CrI, P(direction) and ELPD differences are
# returned in a single tidy table sorted by strength of evidence.
# =============================================================================

library(rgetne)

raw <- readRDS(example_data())
dat <- prep_data(raw)

base_clono <- ~ ki67_percent + periop_therapy + primary_surgery_yes +
                 stage_III   + site_pancreas  + site_colorectal
kinet      <- ~ ki67_percent + periop_therapy

pairs <- list(
  c("ki67_percent",        "periop_therapy"),
  c("ki67_percent",        "primary_surgery_yes"),
  c("ki67_percent",        "site_pancreas"),
  c("periop_therapy",      "site_pancreas"),
  c("primary_surgery_yes", "stage_III"),
  c("primary_surgery_yes", "site_colorectal"),    # "Other site" proxy
  c("periop_therapy",      "primary_surgery_yes"),
  c("periop_therapy",      "stage_III"),
  c("stage_III",           "site_pancreas"),
  c("primary_surgery_yes", "site_pancreas")
)

labels <- c(
  "ki67_percent:periop_therapy"          = "Ki-67 \u00d7 Periop. chemotherapy",
  "ki67_percent:primary_surgery_yes"     = "Ki-67 \u00d7 Primary surgery",
  "ki67_percent:site_pancreas"           = "Ki-67 \u00d7 Pancreatic site",
  "periop_therapy:site_pancreas"         = "Periop. chemo \u00d7 Pancreatic site",
  "primary_surgery_yes:stage_III"        = "Primary surgery \u00d7 Stage III",
  "primary_surgery_yes:site_colorectal"  = "Primary surgery \u00d7 Other site",
  "periop_therapy:primary_surgery_yes"   = "Periop. chemo \u00d7 Primary surgery",
  "periop_therapy:stage_III"             = "Periop. chemo \u00d7 Stage III",
  "stage_III:site_pancreas"              = "Stage III \u00d7 Pancreatic site",
  "primary_surgery_yes:site_pancreas"    = "Primary surgery \u00d7 Pancreatic site"
)

a5 <- screen_interactions(
  data            = dat,
  base_clonogenic = base_clono,
  kinetic         = kinet,
  pairs           = pairs,
  labels          = labels,
  chains          = 4, iter = 4000, warmup = 1000, seed = 42
)

print(a5)

# Render as a publication-style table
tab <- format_interaction_table(a5)
print(tab, row.names = FALSE)
# knitr::kable(tab, align = "lrlrrl",
#              caption = "Supplementary Table A5. Pairwise interaction screen.")
