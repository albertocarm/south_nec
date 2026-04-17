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
  title = "Bayesian Promotion Time Cure Model",
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
      card(
        card_header("Variable summary"),
        verbatimTextOutput("data_summary")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Ki-67 distribution"),
          plotOutput("hist_ki67", height = "260px")
        ),
        card(
          card_header("Event-free survival time distribution"),
          plotOutput("hist_efs", height = "260px")
        )
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Kaplan-Meier",
      layout_sidebar(
        sidebar = sidebar(
          width = 280,
          h5("KM Settings"),
          selectInput("km_strata", "Stratify by:",
                      choices  = c("None", BINARY_VARS),
                      selected = "periop_therapy"),
          numericInput("km_tmax", "Max time (leave blank for full range):",
                       value = NA, min = 0, step = 1),
          actionButton("km_btn", "Plot Kaplan-Meier",
                       class = "btn-info", width = "100%")
        ),
        card(
          card_header("Kaplan-Meier estimate of event-free survival"),
          plotOutput("km_plot", height = "480px")
        )
      )
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Results summary",
      verbatimTextOutput("summary_out")
    ),

    # -------------------------------------------------------------------------
    nav_panel(
      "Diagnostics",
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
      "Identifiability",
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
      "Visualisations",
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
      "Marginal effects",
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
            p("Carmona-Bayonas A, Garcia-Carbonero R, Jimenez-Fonseca P,",
              "et al.",
              tags$em("Perioperative Chemotherapy and Long-Term Cure in",
                      "Resected Grade 3 Gastroenteropancreatic Neuroendocrine",
                      "Carcinomas: The GETNE/SOUTH-NEC Study."), "2026."),
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
                     "github.com/albertocarm/rgetne")),
            p("Citation: Carmona-Bayonas A et al. rgetne: Bayesian Promotion",
              "Time Cure Model Toolkit, 2026.",
              tags$a(href = "https://github.com/albertocarm/rgetne",
                     "https://github.com/albertocarm/rgetne"))
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
                         precomputed = FALSE, me_res = NULL)

  # --- Per-variable prior UI -------------------------------------------------
  default_prior_for_var <- function(v) if (v %in% CONTINUOUS_VARS) 1.0 else 2.5

  build_prior_inputs <- function(id_prefix, selected_vars,
                                 interaction_term = NULL) {
    if (length(selected_vars) == 0) {
      return(helpText("Select at least one main effect."))
    }
    rows <- lapply(selected_vars, function(v) {
      numericInput(paste0(id_prefix, v),
                   label = v,
                   value = default_prior_for_var(v),
                   step  = 0.1, min = 0.01)
    })
    if (!is.null(interaction_term)) {
      rows <- c(rows,
                list(numericInput(paste0(id_prefix, interaction_term),
                                  label = interaction_term,
                                  value = 1.0,
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
                       interaction_name())
  })

  output$prior_kinet_ui <- renderUI({
    build_prior_inputs("prior_kinet_",
                       input$kinet_vars,
                       interaction_term = NULL)
  })

  collect_prior_list <- function(id_prefix, selected_vars,
                                 interaction_term = NULL) {
    terms <- c(selected_vars, interaction_term)
    out <- list()
    for (v in terms) {
      val <- input[[paste0(id_prefix, v)]]
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

  output$data_preview <- renderTable({
    req(rval$data)
    head(rval$data, 10)
  })

  output$data_summary <- renderPrint({
    req(rval$data)
    summary(rval$data[, intersect(MODEL_VARS, names(rval$data)), drop = FALSE])
  })

  output$data_desc_table <- renderTable({
    req(rval$data)
    data_summary_table(rval$data)
  }, striped = TRUE, hover = TRUE)

  output$hist_ki67 <- renderPlot({
    req(rval$data)
    opar <- par(mar = c(4, 4, 2, 1))
    on.exit(par(opar), add = TRUE)
    hist(rval$data$ki67_percent, breaks = 20,
         col = "#4A90D9", border = "white",
         main = "Ki-67 (%)", xlab = "Ki-67 (%)", ylab = "Frequency")
  })

  output$hist_efs <- renderPlot({
    req(rval$data)
    opar <- par(mar = c(4, 4, 2, 1))
    on.exit(par(opar), add = TRUE)
    hist(rval$data$EFS_time, breaks = 25,
         col = "#1B4F8A", border = "white",
         main = "EFS time", xlab = "Time (years)", ylab = "Frequency")
  })

  # --- Kaplan-Meier -----------------------------------------------------------
  km_state <- reactiveValues(sf = NULL)

  observeEvent(input$km_btn, {
    req(rval$data)
    km_state$sf <- TRUE
  })

  output$km_plot <- renderPlot({
    req(km_state$sf)
    strat <- if (!is.null(input$km_strata) && input$km_strata != "None")
               input$km_strata else NULL
    tmax  <- if (!is.na(input$km_tmax) && input$km_tmax > 0)
               input$km_tmax else NULL
    plot_km(rval$data, strata = strat, time_max = tmax)
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
