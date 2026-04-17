# =============================================================================
# rgetne :: Interactive Bayesian Promotion Time Cure Model app
# =============================================================================
# This file is invoked by `cure()` via `shiny::runApp(system.file(
# "shiny", package = "rgetne"))`. It assumes the `rgetne` package is attached.
# =============================================================================

# -----------------------------------------------------------------------------
# This file assumes it is sourced in a session where the packages `shiny`
# and `bslib` are already attached (so that bare names like `fluidPage` or
# `page_sidebar` resolve). Two entry points guarantee this:
#
#   * `rgetne::cure()` in R/cure_app.R attaches shiny + bslib before calling
#     shiny::runApp() on this directory.
#   * The standalone `app.R` at the repo root (used for shinyapps.io)
#     attaches shiny + bslib and sources R/toolkit.R before sourcing us.
#
# No library() calls here: package Imports are declared in DESCRIPTION and
# toolkit symbols are exposed as part of rgetne.
# -----------------------------------------------------------------------------
if (requireNamespace("rgetne", quietly = TRUE) &&
    !"package:rgetne" %in% search()) {
  attachNamespace("rgetne")
}

# -----------------------------------------------------------------------------
# Resolve bundled assets
# -----------------------------------------------------------------------------
# `stan_file_cached` is defined by app.R in standalone mode; provide a default
# that relies on stan_model_file() from the toolkit.
if (!exists("stan_file_cached", mode = "function")) {
  stan_file_cached <- function() {
    src <- stan_model_file()
    dst <- file.path(tempdir(), basename(src))
    if (!file.exists(dst) || file.info(src)$mtime > file.info(dst)$mtime) {
      file.copy(src, dst, overwrite = TRUE)
    }
    dst
  }
}

default_data_file <- function() example_data()

# Path to precomputed model fit (if available)
if (!exists("precomputed_cfit_path", mode = "function")) {
  precomputed_cfit_path <- function() {
    p <- system.file("extdata", "cfit_default.rds", package = "rgetne")
    if (nzchar(p) && file.exists(p)) p else NULL
  }
}

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
MODEL_VARS <- c(
  "primary_surgery_yes", "stage_II", "stage_III",
  "site_colorectal", "site_pancreas",
  "ki67_percent", "periop_therapy"
)

CONTINUOUS_VARS <- c("ki67_percent")

BINARY_VARS <- setdiff(MODEL_VARS, CONTINUOUS_VARS)

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- page_sidebar(
  title = tagList(
    span("Bayesian Promotion Time Cure Model"),
    span(
      style = paste(
        "font-size: 0.7em; margin-left: 12px; padding: 2px 10px;",
        "background:#27ae60; color:white; border-radius: 10px;",
        "vertical-align: middle;"),
      sprintf("build %s",
              tryCatch(as.character(utils::packageVersion("rgetne")),
                       error = function(e) "dev")),
      " \u00b7 JCO theme")),
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C3E50"),
  window_title = "rgetne :: Cure app",

  sidebar = sidebar(
    width = 380,
    h4("Model Specification"),

    accordion(
      open = c("Clonogenic Model (Cure Fraction)",
               "Kinetic Model (Latency)"),

      accordion_panel(
        "Data",
        helpText("A preprocessed example dataset ships with the package."),
        fileInput("user_file",
                  "Optional: upload your own .rds",
                  accept = ".rds"),
        verbatimTextOutput("data_info")
      ),

      accordion_panel(
        "Clonogenic Model (Cure Fraction)",
        helpText("Covariates affecting the probability of ultimate cure",
                 "(incidence, log-theta)."),
        checkboxGroupInput("clono_vars", "Main effects:",
                           choices  = MODEL_VARS,
                           selected = MODEL_VARS),
        hr(),
        h6("Interaction (optional)"),
        selectInput("int_var1", "Variable 1:",
                    choices = c("None", MODEL_VARS),
                    selected = "ki67_percent"),
        selectInput("int_var2", "Variable 2:",
                    choices = c("None", MODEL_VARS),
                    selected = "periop_therapy"),
        hr(),
        h6("Prior SD per variable (clonogenic)"),
        helpText(style = "font-size:0.82em;",
                 "Continuous coefficients live on the standardised scale",
                 "(see scaling method below); binary/categorical",
                 "coefficients live on the raw 0/1 scale."),
        uiOutput("prior_clono_ui")
      ),

      accordion_panel(
        "Kinetic Model (Latency)",
        helpText("Covariates affecting the time to recurrence for non-cured",
                 "patients (log-lambda)."),
        checkboxGroupInput("kinet_vars", "Main effects:",
                           choices  = MODEL_VARS,
                           selected = c("ki67_percent", "periop_therapy")),
        hr(),
        h6("Prior SD per variable (kinetic)"),
        uiOutput("prior_kinet_ui")
      ),

      accordion_panel(
        "MCMC & computation",
        selectInput("scaling_method",
                    "Scaling method for continuous predictors:",
                    choices  = c("Gelman (mean, 2*SD)" = "gelman",
                                 "Z-score (mean, SD)"  = "zscore"),
                    selected = "gelman"),
        numericInput("mcmc_chains", "Chains:",     value = 4,    min = 1, max = 8),
        numericInput("mcmc_iter",   "Iterations:", value = 4000, step = 500),
        numericInput("mcmc_warmup", "Warmup:",     value = 1000, step = 500),
        numericInput("mcmc_adapt",  "Adapt delta:", value = 0.95,
                     step = 0.01, min = 0.8, max = 0.999),
        numericInput("prior_int", "Intercepts prior SD:",
                     value = 2.5, step = 0.5, min = 0.1)
      )
    ),

    hr(),
    actionButton("run_model", "Run Bayesian Model",
                 class = "btn-primary", width = "100%"),
    helpText(
      style = "color:#e67e22; font-size:0.85em; margin-top:6px;",
      icon("clock"),
      "Model fitting typically takes 3\u201310 minutes depending on data",
      "size and MCMC settings. Please be patient while Stan compiles",
      "and samples."
    ),
    uiOutput("precomputed_badge"),
    br(),
    uiOutput("formula_preview")
  ),

  navset_card_underline(
    id = "main_tabs",

    # -------------------------------------------------------------------------
    nav_panel(
      "Overview",
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Clinical interpretation"),
          card_body(
            p("The Promotion Time Cure Model separates survival into two",
              "biological layers:"),
            tags$ul(
              tags$li(strong("Clonogenic (cure fraction):"),
                      " probability that a patient is completely cured",
                      " (no residual clonogenic cells). Driven by the",
                      " covariates in the clonogenic model."),
              tags$li(strong("Kinetic (time-to-event):"),
                      " for patients who are not cured, this is the speed",
                      " at which disease recurs. Driven by the kinetic model.")
            ),
            p("Positive clonogenic coefficient = more clonogenic load =",
              "lower cure probability. Positive kinetic coefficient =",
              "longer time-to-recurrence (slower progression).")
          )
        ),
        card(
          card_header("Technical notes"),
          card_body(
            p("Bayesian parametric survival model using a bounded cumulative",
              "hazard (Yakovlev & Tsodikov)."),
            tags$ul(
              tags$li(strong("Sampler:"), " Hamiltonian Monte Carlo (NUTS)",
                      " via rstan."),
              tags$li(strong("Baseline:"), " Weibull latency."),
              tags$li(strong("Priors:"), " weakly informative Normal priors",
                      " on every coefficient, configurable on the sidebar.")
            ),
            p("Assess convergence via R-hat (< 1.01) and ESS (> 400) in the",
              "Diagnostics tab.")
          )
        )
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Data",
      layout_columns(
        col_widths = c(8, 4),
        uiOutput("data_source_info"),
        uiOutput("data_source_ui")
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header("Preview (first 10 rows)"),
          tableOutput("data_preview")
        ),
        card(
          card_header("Descriptive statistics"),
          tableOutput("data_desc_table")
        )
      ),
      layout_columns(
        col_widths = c(4, 8),
        card(
          card_header("Explore a variable"),
          uiOutput("explore_var_ui"),
          radioButtons("explore_type", "Plot type:",
                       choices = c("Auto"      = "auto",
                                   "Density"   = "density",
                                   "Histogram" = "histogram",
                                   "Bar"       = "bar",
                                   "Pie"       = "pie"),
                       selected = "auto", inline = TRUE)
        ),
        card(
          card_header("Distribution"),
          plotOutput("explore_plot", height = "380px")
        )
      ),
      card(
        card_header("Ki-67 density and categorical overview"),
        layout_columns(
          col_widths = c(6, 6),
          plotOutput("overview_ki67", height = "320px"),
          plotOutput("overview_efs",  height = "320px")
        ),
        plotOutput("overview_cats", height = "460px")
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Kaplan-Meier",
      layout_sidebar(
        sidebar = sidebar(
          width = 300,
          h5("KM Settings"),
          helpText(
            style = "font-size:0.82em; color:#555;",
            "Uses the publication dataset (factor-coded).",
            "Stratify, facet and filter to explore interactions."),
          selectInput("km_outcome", "Outcome:",
                      choices = c("Event-Free Survival" = "EFS",
                                  "Overall Survival"    = "OS"),
                      selected = "EFS"),
          selectInput("km_strata", "Stratify by (colour):",
                      choices = "None", selected = "None"),
          selectInput("km_facet", "Facet by (panels):",
                      choices = "None", selected = "None"),
          selectInput("km_filter_var", "Restrict to a subset of:",
                      choices = "None", selected = "None"),
          uiOutput("km_filter_levels_ui"),
          numericInput("km_tmax", "Max time (blank = full range):",
                       value = NA, min = 0, step = 1),
          actionButton("km_btn", "Plot Kaplan-Meier",
                       class = "btn-info", width = "100%")
        ),
        card(
          card_header(textOutput("km_card_title")),
          plotOutput("km_plot", height = "520px")
        )
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Cox models",
      layout_sidebar(
        sidebar = sidebar(
          width = 320,
          h5("Cox multivariable"),
          helpText(
            style = "font-size:0.82em; color:#555;",
            "Fits a multivariable Cox model on the publication dataset for",
            "both EFS and OS, builds forest plots and a combined gt table."),
          uiOutput("cox_vars_ui"),
          actionButton("cox_btn", "Fit EFS + OS",
                       class = "btn-success", width = "100%")
        ),
        navset_card_underline(
          nav_panel(
            "Forest plots",
            layout_columns(
              col_widths = c(6, 6),
              card(card_header("Event-Free Survival"),
                   plotOutput("cox_plot_efs", height = "440px")),
              card(card_header("Overall Survival"),
                   plotOutput("cox_plot_os",  height = "440px"))
            )
          ),
          nav_panel(
            "Combined table",
            card(card_header("Multivariable Cox Models (EFS + OS)"),
                 gt::gt_output("cox_table"))
          )
        )
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Cure model \u2014 results",
      verbatimTextOutput("summary_out")
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Cure model \u2014 diagnostics",
      h5("MCMC trace plots"),
      plotOutput("mcmc_trace_plot", height = "320px"),
      hr(),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("MCMC convergence"),
             verbatimTextOutput("mcmc_out")),
        card(card_header("Posterior cross-correlation"),
             verbatimTextOutput("cor_out"))
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Cure model \u2014 identifiability",
      helpText("The incremental censoring test refits the model at decreasing",
               "follow-up quantiles to check whether the cure fraction is",
               "identifiable from the available events. This can be slow."),
      actionButton("run_cens", "Run censoring sensitivity",
                   class = "btn-warning"),
      helpText(
        style = "color:#e67e22; font-size:0.85em; margin-top:6px;",
        icon("clock"),
        "This test refits the model multiple times at different",
        "follow-up horizons. It can take 5\u201315 minutes to complete."
      ),
      hr(),
      verbatimTextOutput("cens_out")
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Cure model \u2014 visualisations",
      layout_sidebar(
        sidebar = sidebar(
          width = 300,
          h5("Plot configuration"),
          selectInput("var_mod", "Continuous X-axis:",
                      choices  = CONTINUOUS_VARS,
                      selected = "ki67_percent"),
          selectInput("var_contraste", "Stratification variable:",
                      choices  = BINARY_VARS,
                      selected = "periop_therapy"),
          numericInput("val_expuesto", "Exposed level:",     value = 1),
          numericInput("val_referencia", "Reference level:", value = 0),
          actionButton("plot_btn", "Generate plots", class = "btn-success")
        ),
        card(
          card_header("Cure probability by continuous variable (reference group)"),
          plotOutput("plot_triptico_ref", height = "360px")
        ),
        card(
          card_header("Cure probability by continuous variable (exposed group)"),
          plotOutput("plot_triptico_exp", height = "360px")
        ),
        card(
          card_header("Cure odds ratio (exposed vs reference)"),
          plotOutput("plot_or_cure", height = "500px")
        ),
        card(
          card_header("Baseline population survival curve"),
          plotOutput("plot_surv", height = "360px")
        )
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Cure model \u2014 marginal effects",
      layout_sidebar(
        sidebar = sidebar(
          width = 320,
          h5("Scenario grid"),
          helpText(style = "font-size:0.85em;",
                   "Pick a continuous X-axis (scanned over its observed",
                   "range) and an optional binary variable to split the",
                   "lines. All other covariates are held fixed at their",
                   "sample mean (continuous) or 0 (binary)."),
          selectInput("me_target", "Quantity:",
                      choices = c("Cure probability"  = "cure_prob",
                                  "Clonogen count"    = "theta",
                                  "Latency scale"     = "lambda",
                                  "Median event time" = "median_time"),
                      selected = "cure_prob"),
          selectInput("me_xvar", "Continuous X-axis:",
                      choices  = CONTINUOUS_VARS,
                      selected = "ki67_percent"),
          selectInput("me_groupvar", "Group by (binary):",
                      choices  = c("None", BINARY_VARS),
                      selected = "periop_therapy"),
          numericInput("me_ngrid", "Grid resolution:",
                       value = 40, min = 5, max = 200, step = 5),
          actionButton("me_btn", "Compute marginal effects",
                       class = "btn-success", width = "100%")
        ),
        card(
          card_header("Posterior expected prediction"),
          plotOutput("me_plot", height = "480px")
        ),
        card(
          card_header("Summary table"),
          tableOutput("me_table")
        )
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "About",
      layout_columns(
        col_widths = c(12),
        card(
          card_header("About rgetne"),
          card_body(
            p(strong("rgetne"), "is the reproducibility companion for the",
              "GETNE/SOUTH-NEC study. It fits Bayesian Promotion Time Cure",
              "Models via Stan and ships this interactive Shiny front-end."),
            h6("Reference paper"),
            p(tags$em("Perioperative Chemotherapy and Long-Term Cure in",
                      "Resected Grade 3 Gastroenteropancreatic Neuroendocrine",
                      "Carcinomas: The GETNE/SOUTH-NEC Study."),
              " Manuscript under review."),
            h6("What this app does"),
            tags$ul(
              tags$li("Fits a Bayesian promotion-time cure model (Yakovlev &",
                      "Tsodikov) that separates treatment effects into a",
                      tags$strong("clonogenic"), "(cure fraction) and",
                      tags$strong("kinetic"), "(time-to-recurrence) component."),
              tags$li("Uses data from 176 patients with stage I-III GEP-NEC",
                      "from the Spanish GETNE registry (R-GETNE)."),
              tags$li("Supports interactive covariate selection, interactions,",
                      "per-variable priors, MCMC diagnostics, Kaplan-Meier",
                      "plots, cure probability curves, and identifiability",
                      "tests.")
            ),
            h6("Links"),
            p("Source code: ",
              tags$a(href = "https://github.com/albertocarm/rgetne",
                     "github.com/albertocarm/rgetne"))
          )
        )
      )
    )
  )
)

# Poor-man's %||% if not attached yet
`%||%` <- function(a, b) if (is.null(a)) b else a

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  rval <- reactiveValues(data = NULL, cfit = NULL, cens_out = NULL,
                         precomputed = FALSE, me_res = NULL,
                         pub = NULL, cox_res = NULL)

  # --- Publication dataset (factor-coded) ------------------------------------
  # Used by the Kaplan-Meier exploration tab and the Cox multivariable tab.
  # Loaded lazily from inst/extdata/dataset_publicacion_english.rds.
  observe({
    if (!is.null(rval$pub)) return()
    pub_path <- tryCatch(publication_data(), error = function(e) NULL)
    if (is.null(pub_path) || !file.exists(pub_path)) return()
    tryCatch({
      rval$pub <- readRDS(pub_path)
    }, error = function(e) {
      showNotification(paste("Could not load publication dataset:", e$message),
                       type = "warning")
    })
  })

  # Publication dataset augmented with a Ki-67 quartile factor so that it can
  # be selected as strata, facet or filter in the KM and Cox tabs the same
  # way as any other categorical variable.
  pub_aug <- reactive({
    df <- rval$pub; req(df)
    if ("Ki67_Index" %in% names(df)) {
      q <- stats::quantile(df$Ki67_Index, c(.25, .5, .75), na.rm = TRUE)
      df$Ki67_Q4 <- cut(df$Ki67_Index,
                        breaks = c(-Inf, q, Inf),
                        labels = c(sprintf("Q1 [<=%.0f%%]",    q[1]),
                                   sprintf("Q2 (%.0f-%.0f%%]", q[1], q[2]),
                                   sprintf("Q3 (%.0f-%.0f%%]", q[2], q[3]),
                                   sprintf("Q4 (>%.0f%%)",    q[3])),
                        include.lowest = TRUE)
    }
    df
  })

  pub_factor_vars <- reactive({
    df <- pub_aug(); req(df)
    keep <- vapply(df, function(x) {
      is.factor(x) || is.character(x) ||
        (is.numeric(x) && length(unique(stats::na.omit(x))) <= 6)
    }, logical(1))
    drop <- c("EFS_Time", "EFS_Status", "OS_Time", "OS_Status",
              "Ki67_Index", "Ki67_Percent")
    setdiff(names(df)[keep], drop)
  })

  pub_model_vars <- reactive({
    df <- rval$pub; req(df)
    candidates <- c("Sex", "TNM_Stage", "Primary_Surgery", "Site_4_Groups",
                    "Ki67_Index", "Perioperative_Chemo", "ECOG_PS")
    intersect(candidates, names(df))
  })

  # Dataset driving the Data tab: prefer the publication dataset when
  # available (richer set of covariates), fall back to the preprocessed one.
  data_explorer_df <- reactive({
    if (!is.null(rval$pub)) pub_aug() else rval$data
  })
  data_explorer_src <- reactive({
    if (!is.null(rval$pub)) "publication" else "preprocessed"
  })

  # --- Per-variable prior UI -------------------------------------------------
  # Shiny rejects input IDs containing ":" (interpreted as a type separator),
  # so interaction terms like "ki67_percent:periop_therapy" are sanitised to
  # "ki67_percent__x__periop_therapy" for the ID while the original name is
  # preserved as label and as the key of the prior list passed to Stan.
  safe_id  <- function(v) gsub(":", "__x__", v, fixed = TRUE)
  TIGHT_VARS <- c("ki67_percent", "periop_therapy")
  default_prior_for_var <- function(v, component = c("clonogenic", "kinetic")) {
    component <- match.arg(component)
    if (component == "clonogenic") return(2.5)
    parts <- strsplit(v, ":", fixed = TRUE)[[1]]
    if (any(parts %in% TIGHT_VARS))      return(0.25)
    if (any(parts %in% CONTINUOUS_VARS)) return(1.0)
    2.5
  }

  build_prior_inputs <- function(id_prefix, selected_vars,
                                 interaction_term = NULL,
                                 component = c("clonogenic", "kinetic")) {
    component <- match.arg(component)
    if (length(selected_vars) == 0) {
      return(helpText("Select at least one main effect."))
    }
    rows <- lapply(selected_vars, function(v) {
      numericInput(paste0(id_prefix, safe_id(v)),
                   label = v,
                   value = default_prior_for_var(v, component),
                   step  = 0.1, min = 0.01)
    })
    if (!is.null(interaction_term)) {
      rows <- c(rows,
                list(numericInput(paste0(id_prefix, safe_id(interaction_term)),
                                  label = interaction_term,
                                  value = default_prior_for_var(interaction_term,
                                                                component),
                                  step  = 0.1, min = 0.01)))
    }
    do.call(tagList, rows)
  }

  interaction_name <- reactive({
    v1 <- input$int_var1; v2 <- input$int_var2
    if (is.null(v1) || is.null(v2) || v1 == "None" || v2 == "None" || v1 == v2)
      return(NULL)
    paste(v1, v2, sep = ":")
  })

  output$prior_clono_ui <- renderUI({
    build_prior_inputs("prior_clono_",
                       input$clono_vars,
                       interaction_name(),
                       component = "clonogenic")
  })

  output$prior_kinet_ui <- renderUI({
    build_prior_inputs("prior_kinet_",
                       input$kinet_vars,
                       interaction_term = NULL,
                       component = "kinetic")
  })

  collect_prior_list <- function(id_prefix, selected_vars,
                                 interaction_term = NULL) {
    terms <- c(selected_vars, interaction_term)
    out <- list()
    for (v in terms) {
      val <- input[[paste0(id_prefix, safe_id(v))]]
      if (!is.null(val) && is.numeric(val) && val > 0) out[[v]] <- val
    }
    out
  }

  # --- Load data (default or uploaded) ---------------------------------------
  observe({
    uploaded <- !is.null(input$user_file)
    path <- if (uploaded) input$user_file$datapath else default_data_file()
    tryCatch({
      raw <- readRDS(path)
      dat <- prep_data(raw)
      rval$data <- dat
      rval$cfit <- NULL
      rval$cens_out <- NULL
      rval$precomputed <- FALSE

      # Auto-load precomputed model if using default data
      if (!uploaded) {
        pre_path <- precomputed_cfit_path()
        if (!is.null(pre_path)) {
          tryCatch({
            rval$cfit <- readRDS(pre_path)
            rval$precomputed <- TRUE
            showNotification(
              sprintf("Data loaded: %d obs. Precomputed model results ready.",
                      nrow(dat)),
              type = "message", duration = 8)
          }, error = function(e) {
            showNotification(sprintf("Data loaded: %d observations.", nrow(dat)),
                             type = "message")
          })
        } else {
          showNotification(sprintf("Data loaded: %d observations.", nrow(dat)),
                           type = "message")
        }
      } else {
        showNotification(sprintf("Data loaded: %d observations.", nrow(dat)),
                         type = "message")
      }
    }, error = function(e) {
      showNotification(paste("Failed to load data:", e$message),
                       type = "error", duration = NULL)
    })
  })

  output$precomputed_badge <- renderUI({
    if (isTRUE(rval$precomputed)) {
      helpText(
        style = "color:#27ae60; font-size:0.85em; margin-top:4px;",
        icon("check-circle"),
        "Precomputed results loaded. You can browse all tabs immediately.",
        "Press 'Run Bayesian Model' to refit with different settings."
      )
    }
  })

  output$data_info <- renderPrint({
    req(rval$data)
    cat(sprintf("n = %d rows, %d columns\n",
                nrow(rval$data), ncol(rval$data)))
    cat(sprintf("events = %d (%.0f%%)\n",
                sum(rval$data$.status_bin),
                100 * mean(rval$data$.status_bin)))
  })

  output$data_source_info <- renderUI({
    src <- data_explorer_src()
    n   <- nrow(data_explorer_df() %||% data.frame())
    k   <- ncol(data_explorer_df() %||% data.frame())
    badge <- if (src == "publication") "Publication dataset (factor-coded)"
             else "Preprocessed cure-model dataset"
    helpText(
      style = "font-size:0.95em;",
      icon("database"), " Showing: ", strong(badge),
      sprintf(" (n = %d, %d columns).", n, k))
  })

  output$data_source_ui <- renderUI({
    NULL
  })

  output$data_preview <- renderTable({
    df <- data_explorer_df(); req(df)
    utils::head(df, 10)
  })

  output$data_desc_table <- renderTable({
    df <- data_explorer_df(); req(df)
    if (data_explorer_src() == "publication") {
      rows <- lapply(names(df), function(v) {
        x <- df[[v]]
        x <- x[!is.na(x)]
        if (length(x) == 0) return(NULL)
        if (is.numeric(x) && length(unique(x)) > 6) {
          data.frame(variable = v, stat = "median [IQR]",
                     value = sprintf("%.1f [%.1f-%.1f]",
                                     stats::median(x),
                                     stats::quantile(x, .25),
                                     stats::quantile(x, .75)),
                     stringsAsFactors = FALSE)
        } else {
          tab <- table(factor(x))
          lbl <- paste(sprintf("%s: %d (%.0f%%)",
                               names(tab), as.integer(tab),
                               100 * as.integer(tab) / sum(tab)),
                       collapse = "; ")
          data.frame(variable = v, stat = "n (%)",
                     value = lbl, stringsAsFactors = FALSE)
        }
      })
      do.call(rbind, rows)
    } else {
      data_summary_table(df)
    }
  }, striped = TRUE, hover = TRUE)

  output$explore_var_ui <- renderUI({
    df <- data_explorer_df(); req(df)
    ch <- setdiff(names(df), c(".status_bin"))
    default <- if ("Ki67_Index" %in% ch) "Ki67_Index"
               else if ("ki67_percent" %in% ch) "ki67_percent"
               else ch[1]
    selectInput("explore_var", "Variable:", choices = ch, selected = default)
  })

  output$explore_plot <- renderPlot({
    df <- data_explorer_df(); req(df)
    v  <- input$explore_var; req(v, v %in% names(df))
    plot_variable(df, v, type = input$explore_type %||% "auto")
  })

  output$overview_ki67 <- renderPlot({
    df <- data_explorer_df(); req(df)
    v <- if ("Ki67_Index" %in% names(df)) "Ki67_Index"
         else if ("ki67_percent" %in% names(df)) "ki67_percent"
         else return(NULL)
    plot_variable(df, v, type = "density",
                  title = "Ki-67 proliferation index")
  })

  output$overview_efs <- renderPlot({
    df <- data_explorer_df(); req(df)
    v <- if ("EFS_Time" %in% names(df)) "EFS_Time"
         else if ("EFS_time" %in% names(df)) "EFS_time"
         else return(NULL)
    plot_variable(df, v, type = "density",
                  title = "Event-free survival time")
  })

  output$overview_cats <- renderPlot({
    df <- data_explorer_df(); req(df)
    cat_vars <- setdiff(pub_factor_vars(), "Ki67_Q4")
    if (length(cat_vars) == 0) return(NULL)
    plots <- lapply(cat_vars, function(v) plot_variable(df, v, type = "bar"))
    gridExtra::grid.arrange(grobs = plots,
                            ncol = min(3, length(plots)))
  })

  # --- Kaplan-Meier (publication dataset) ------------------------------------
  km_state <- reactiveValues(triggered = FALSE)

  observe({
    req(rval$pub)
    choices <- c("None", pub_factor_vars())
    updateSelectInput(session, "km_strata",     choices = choices,
                      selected = isolate(input$km_strata)     %||% "None")
    updateSelectInput(session, "km_facet",      choices = choices,
                      selected = isolate(input$km_facet)      %||% "None")
    updateSelectInput(session, "km_filter_var", choices = choices,
                      selected = isolate(input$km_filter_var) %||% "None")
  })

  output$km_filter_levels_ui <- renderUI({
    v <- input$km_filter_var
    if (is.null(v) || v == "None" || is.null(rval$pub)) return(NULL)
    lv <- sort(unique(as.character(rval$pub[[v]])))
    selectInput("km_filter_levels", paste("Levels of", v, "to keep:"),
                choices = lv, selected = lv, multiple = TRUE)
  })

  observeEvent(input$km_btn, { km_state$triggered <- TRUE })

  output$km_card_title <- renderText({
    if (isTRUE(km_state$triggered) && !is.null(input$km_outcome)) {
      if (input$km_outcome == "OS") "Kaplan-Meier \u2014 Overall Survival"
      else "Kaplan-Meier \u2014 Event-Free Survival"
    } else "Kaplan-Meier estimate"
  })

  # Prevent picking the same variable for strata and facet: when the user
  # chooses one dropdown, collapse the other to "None" if they collide.
  observeEvent(input$km_strata, {
    if (!is.null(input$km_facet) && !is.null(input$km_strata) &&
        input$km_strata != "None" && input$km_strata == input$km_facet) {
      updateSelectInput(session, "km_facet", selected = "None")
      showNotification(
        "Stratify and facet cannot be the same variable; facet reset.",
        type = "warning", duration = 4)
    }
  }, ignoreInit = TRUE)
  observeEvent(input$km_facet, {
    if (!is.null(input$km_facet) && !is.null(input$km_strata) &&
        input$km_facet != "None" && input$km_facet == input$km_strata) {
      updateSelectInput(session, "km_strata", selected = "None")
      showNotification(
        "Stratify and facet cannot be the same variable; stratify reset.",
        type = "warning", duration = 4)
    }
  }, ignoreInit = TRUE)

  output$km_plot <- renderPlot({
    req(km_state$triggered, rval$pub)
    outcome <- input$km_outcome %||% "EFS"
    time_col   <- if (outcome == "OS") "OS_Time"   else "EFS_Time"
    status_col <- if (outcome == "OS") "OS_Status" else "EFS_Status"
    strat <- if (!is.null(input$km_strata) && input$km_strata != "None")
               input$km_strata else NULL
    fct   <- if (!is.null(input$km_facet) && input$km_facet != "None")
               input$km_facet else NULL
    if (!is.null(strat) && !is.null(fct) && identical(strat, fct)) {
      fct <- NULL
    }
    fvar  <- if (!is.null(input$km_filter_var) && input$km_filter_var != "None")
               input$km_filter_var else NULL
    flev  <- if (!is.null(fvar)) input$km_filter_levels else NULL
    tmax  <- if (!is.na(input$km_tmax) && input$km_tmax > 0)
               input$km_tmax else NULL
    plot_km(pub_aug(),
            time_col = time_col, status_col = status_col,
            strata = strat, facet_by = fct,
            filter_var = fvar, filter_values = flev,
            time_max = tmax,
            ylab = if (outcome == "OS") "Overall survival"
                   else "Event-free survival",
            title = if (outcome == "OS") "Kaplan-Meier \u2014 OS"
                    else "Kaplan-Meier \u2014 EFS")
  })

  # --- Cox multivariable (publication dataset) -------------------------------
  output$cox_vars_ui <- renderUI({
    req(rval$pub)
    choices <- pub_factor_vars()
    if ("Ki67_Index" %in% names(rval$pub)) choices <- c("Ki67_Index", choices)
    defaults <- intersect(pub_model_vars(), choices)
    checkboxGroupInput("cox_vars", "Covariates:",
                       choices  = choices,
                       selected = defaults)
  })

  observeEvent(input$cox_btn, {
    req(rval$pub, input$cox_vars)
    vars <- input$cox_vars
    if (length(vars) < 1L) {
      showNotification("Pick at least one covariate.", type = "warning")
      return()
    }
    df <- pub_aug()
    has_efs <- all(c("EFS_Time", "EFS_Status") %in% names(df))
    has_os  <- all(c("OS_Time",  "OS_Status")  %in% names(df))
    if (!has_efs && !has_os) {
      showNotification("Publication dataset lacks EFS_* and OS_* columns.",
                       type = "error")
      return()
    }
    withProgress(message = "Fitting Cox models\u2026", value = 0.2, {
      res <- tryCatch({
        if (has_efs && has_os) {
          cox_dual_endpoint(df, vars)
        } else {
          tc <- if (has_efs) "EFS_Time"   else "OS_Time"
          sc <- if (has_efs) "EFS_Status" else "OS_Status"
          fit <- cox_multivar_fit(df, tc, sc, vars)
          list(efs = fit, os = fit,
               plot_efs = cox_forest_plot(fit$fp),
               plot_os  = cox_forest_plot(fit$fp),
               table    = cox_combined_table(fit, fit))
        }
      }, error = function(e) {
        showNotification(paste("Cox fit failed:", e$message),
                         type = "error", duration = NULL)
        NULL
      })
      incProgress(0.8)
      rval$cox_res <- res
    })
  })

  output$cox_plot_efs <- renderPlot({
    req(rval$cox_res); rval$cox_res$plot_efs
  })
  output$cox_plot_os <- renderPlot({
    req(rval$cox_res); rval$cox_res$plot_os
  })
  output$cox_table <- gt::render_gt({
    req(rval$cox_res); rval$cox_res$table
  })

  # --- Formula builder -------------------------------------------------------
  build_formulas <- reactive({
    clono_terms <- input$clono_vars
    if (length(clono_terms) == 0) clono_terms <- "1"
    f_clono_str <- paste("~", paste(clono_terms, collapse = " + "))

    if (!is.null(input$int_var1) && !is.null(input$int_var2) &&
        input$int_var1 != "None" && input$int_var2 != "None" &&
        input$int_var1 != input$int_var2) {
      f_clono_str <- paste(f_clono_str, "+",
                           paste(input$int_var1, "*", input$int_var2))
    }

    kinet_terms <- input$kinet_vars
    if (length(kinet_terms) == 0) kinet_terms <- "1"
    f_kinet_str <- paste("~", paste(kinet_terms, collapse = " + "))

    list(
      clono     = as.formula(f_clono_str),
      kinet     = as.formula(f_kinet_str),
      clono_str = f_clono_str,
      kinet_str = f_kinet_str
    )
  })

  output$formula_preview <- renderUI({
    f <- build_formulas()
    tagList(
      h6("Active formulas:"),
      tags$code(style = "color:#2C3E50; font-size:0.85em;", paste("Clono:", f$clono_str)),
      br(),
      tags$code(style = "color:#7F8C8D; font-size:0.85em;", paste("Kinet:", f$kinet_str))
    )
  })

  # --- Run model -------------------------------------------------------------
  observeEvent(input$run_model, {
    req(rval$data)

    if (length(input$clono_vars) == 0 || length(input$kinet_vars) == 0) {
      showNotification(
        "Both the clonogenic and kinetic models require at least one covariate.",
        type = "error", duration = NULL)
      return(NULL)
    }

    withProgress(message = "Compiling & running Stan model...", value = 0, {
      tryCatch({
        forms <- build_formulas()

        incProgress(0.1, detail = "Standardising predictors")

        clono_priors <- collect_prior_list("prior_clono_",
                                           input$clono_vars,
                                           interaction_name())
        kinet_priors <- collect_prior_list("prior_kinet_",
                                           input$kinet_vars)

        predictors <- unique(c(input$clono_vars, input$kinet_vars,
                               if (!is.null(interaction_name()))
                                 c(input$int_var1, input$int_var2)))

        incProgress(0.15, detail = "NUTS sampling (this can take a while)")

        cfit <- fit_cure_bayes(
          data             = rval$data,
          clonogenic       = forms$clono,
          kinetic          = forms$kinet,
          predictors       = predictors,
          metodo           = input$scaling_method,
          prior_clonogenic = if (length(clono_priors) > 0) clono_priors else NULL,
          prior_kinetic    = if (length(kinet_priors) > 0) kinet_priors else NULL,
          prior_intercept  = input$prior_int,
          chains           = input$mcmc_chains,
          iter             = input$mcmc_iter,
          warmup           = input$mcmc_warmup,
          adapt_delta      = input$mcmc_adapt
        )

        rval$cfit     <- cfit
        rval$cens_out <- NULL
        rval$me_res   <- NULL
        rval$precomputed <- FALSE
        showNotification("Model successfully fitted.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error during model run:", e$message),
                         type = "error", duration = NULL)
      })
    })
  })

  # --- Summary / diagnostics -------------------------------------------------
  output$summary_out <- renderPrint({
    req(rval$cfit)
    summary_coefs(rval$cfit)
  })

  output$mcmc_out <- renderPrint({
    req(rval$cfit)
    mcmc_checks(rval$cfit)
  })

  output$mcmc_trace_plot <- renderPlot({
    req(rval$cfit)
    bayesplot::mcmc_trace(rval$cfit$fit,
                          pars = c("beta0", "gamma0", "alpha"))
  })

  output$cor_out <- renderPrint({
    req(rval$cfit)
    test_correlation(rval$cfit)
  })

  # --- Censoring sensitivity (on demand) -------------------------------------
  observeEvent(input$run_cens, {
    req(rval$cfit)
    withProgress(message = "Refitting at censored horizons...", value = 0, {
      tryCatch({
        tab <- utils::capture.output(
          test_censoring(rval$cfit, chains = 2, iter = 1000)
        )
        rval$cens_out <- paste(tab, collapse = "\n")
      }, error = function(e) {
        rval$cens_out <- paste("Error:", e$message)
      })
    })
  })

  output$cens_out <- renderPrint({
    if (is.null(rval$cens_out))
      cat("Press 'Run censoring sensitivity' to start.\n")
    else
      cat(rval$cens_out)
  })

  # --- Plots (on demand) -----------------------------------------------------
  plot_state <- reactiveValues(ref = NULL, exp = NULL, or = NULL, surv = NULL)

  observeEvent(input$plot_btn, {
    req(rval$cfit)

    v_mod <- input$var_mod
    v_con <- input$var_contraste
    v_exp <- input$val_expuesto
    v_ref <- input$val_referencia

    # Validate continuous var is in clonogenic design
    if (!(v_mod %in% rval$cfit$spec$X_names ||
          v_mod %in% rval$cfit$spec$Z_names)) {
      showNotification(
        sprintf("'%s' is not in the current model. Add it as a main effect.",
                v_mod),
        type = "warning", duration = 8)
      return(NULL)
    }
    if (!(v_con %in% rval$cfit$spec$X_names)) {
      showNotification(
        sprintf("Stratification variable '%s' is not in the clonogenic model.",
                v_con),
        type = "warning", duration = 8)
      return(NULL)
    }

    list_ref <- setNames(list(v_ref), v_con)
    list_exp <- setNames(list(v_exp), v_con)

    plot_state$ref  <- list(cfit = rval$cfit, x_var = v_mod, fix = list_ref)
    plot_state$exp  <- list(cfit = rval$cfit, x_var = v_mod, fix = list_exp)
    plot_state$or   <- list(cfit = rval$cfit, var = v_con,
                            levels = list_exp, ref = list_ref, at_var = v_mod)
    plot_state$surv <- list(cfit = rval$cfit)
  })

  output$plot_triptico_ref <- renderPlot({
    req(plot_state$ref)
    plot_continuous_effect(
      plot_state$ref$cfit,
      x_var = plot_state$ref$x_var,
      fix   = plot_state$ref$fix,
      xlab  = plot_state$ref$x_var)
  })
  output$plot_triptico_exp <- renderPlot({
    req(plot_state$exp)
    plot_continuous_effect(
      plot_state$exp$cfit,
      x_var = plot_state$exp$x_var,
      fix   = plot_state$exp$fix,
      xlab  = plot_state$exp$x_var)
  })
  output$plot_or_cure <- renderPlot({
    req(plot_state$or)
    plot_contrast_cure_or(
      plot_state$or$cfit,
      var     = plot_state$or$var,
      levels  = plot_state$or$levels,
      ref     = plot_state$or$ref,
      at_var  = plot_state$or$at_var,
      save_path = NULL)
  })
  output$plot_surv <- renderPlot({
    req(plot_state$surv)
    plot_surv_curve(plot_state$surv$cfit)
  })

  # --- Marginal effects (posterior_epred_cure) -------------------------------
  me_target_label <- function(type) {
    switch(type,
           cure_prob   = "Cure probability",
           theta       = "Clonogen count (theta)",
           lambda      = "Latency scale (lambda)",
           median_time = "Median event time")
  }

  observeEvent(input$me_btn, {
    req(rval$cfit, rval$data)

    sp <- rval$cfit$spec
    if (is.null(sp$scaling)) {
      showNotification(
        paste("This model has no scaling metadata (it was fit with the",
              "legacy pipeline or loaded from the precomputed cache). Click",
              "'Run Bayesian Model' to refit with the new per-variable prior",
              "pipeline and enable marginal effects."),
        type = "warning", duration = 12)
      return(NULL)
    }

    x_var <- input$me_xvar
    g_var <- input$me_groupvar
    n_grid <- max(5L, as.integer(input$me_ngrid))

    # Covariates required by the training formulas
    req_vars <- unique(c(all.vars(sp$formula_clono),
                         all.vars(sp$formula_kinet)))

    x_raw <- rval$data[[x_var]]
    if (is.null(x_raw)) {
      showNotification(sprintf("Variable '%s' is not in the data.", x_var),
                       type = "error")
      return(NULL)
    }
    x_seq <- seq(min(x_raw, na.rm = TRUE),
                 max(x_raw, na.rm = TRUE),
                 length.out = n_grid)

    use_group <- !is.null(g_var) && g_var != "None" && g_var %in% req_vars
    group_levels <- if (use_group) c(0L, 1L) else NA_integer_

    newdata <- expand.grid(x = x_seq, g = group_levels,
                           stringsAsFactors = FALSE)
    names(newdata) <- c(x_var, if (use_group) g_var else "._drop")
    if (!use_group) newdata[["._drop"]] <- NULL

    # Hold remaining covariates at their mean (continuous) or 0 (binary)
    for (v in setdiff(req_vars, names(newdata))) {
      col <- rval$data[[v]]
      if (is.null(col)) {
        newdata[[v]] <- 0
      } else if (v %in% CONTINUOUS_VARS) {
        newdata[[v]] <- mean(col, na.rm = TRUE)
      } else {
        newdata[[v]] <- 0L
      }
    }

    withProgress(message = "Computing posterior predictions...", value = 0.3, {
      tryCatch({
        ep <- posterior_epred_cure(rval$cfit, newdata, type = input$me_target)
        res <- summarize_epred(ep, newdata = newdata)
        res$x_var <- x_var
        res$group <- if (use_group) factor(newdata[[g_var]]) else factor("all")
        res$target <- input$me_target
        rval$me_res <- list(res = res, x_var = x_var,
                            use_group = use_group, g_var = g_var,
                            target = input$me_target)
      }, error = function(e) {
        showNotification(paste("Marginal effects failed:", e$message),
                         type = "error", duration = 10)
      })
    })
  })

  output$me_plot <- renderPlot({
    me <- rval$me_res
    req(me)

    df <- me$res
    ylab <- me_target_label(me$target)
    x_col <- me$x_var
    if (me$target == "cure_prob") {
      df$median <- df$median * 100
      df$lo     <- df$lo * 100
      df$hi     <- df$hi * 100
      ylab <- paste0(ylab, " (%)")
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = median))
    if (me$use_group) {
      df$group <- factor(df[[me$g_var]])
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = median,
                                            colour = group, fill = group)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                             alpha = 0.18, colour = NA) +
        ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::scale_colour_manual(
          values = c("0" = "#1B4F8A", "1" = "#B8372E"),
          name = me$g_var) +
        ggplot2::scale_fill_manual(
          values = c("0" = "#1B4F8A", "1" = "#B8372E"),
          name = me$g_var)
    } else {
      p <- p +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                             alpha = 0.2, fill = "#1B4F8A") +
        ggplot2::geom_line(linewidth = 1.1, colour = "#1B4F8A")
    }

    p +
      ggplot2::labs(x = x_col, y = ylab,
                    title = sprintf("%s vs %s", ylab, x_col),
                    subtitle = "Posterior median with 95% credible interval") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                     plot.subtitle = ggplot2::element_text(colour = "grey40"),
                     panel.grid.minor = ggplot2::element_blank())
  })

  output$me_table <- renderTable({
    me <- rval$me_res
    req(me)
    df <- me$res
    keep <- c(me$x_var,
              if (me$use_group) me$g_var,
              "median", "lo", "hi")
    out <- df[, intersect(keep, names(df)), drop = FALSE]
    head(out, 20)
  }, striped = TRUE, hover = TRUE, digits = 4)
}

shinyApp(ui, server)
