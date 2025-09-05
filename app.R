# app.R — PROMIS Change Predictor for Mental Health Clinic
# --------------------------------------------------------
# What this app does
# - Takes baseline and follow-up PHQ-9 or GAD-7 scores
# - Lets the user choose a treatment program (models are program-specific)
# - Applies the appropriate prediction model to estimate PROMIS score change
# - Shows a 95% prediction interval (PI) as a plausible range for the change
#
# How to plug in your real models
# - Replace the example `model_store` below with your fitted models (e.g.,
#   saved as RDS files or a data frame of coefficients). The app expects, for
#   each (program, measure) pair, an intercept and coefficients for the baseline
#   score and the change in the score, plus an outcome SD (sigma) for PIs.
# - If you have full model objects (e.g., lm/glmnet/xgboost), modify the
#   `predict_promis_change()` function to call `predict()` accordingly and supply
#   the predictive SD for the PI (or use bootstrap/quantile PIs).
#
# NOTE ON DIRECTION: By convention here, PROMIS change > 0 means improved QoL.
# If PHQ-9 or GAD-7 decreases (improvement), we expect PROMIS change to increase.
# So in the example models below, the coefficient on (followup - baseline) is
# negative. Adjust to your data as needed.

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(glue)


# 1) Load models ---------------------------------------------------------------

# Source the data build (creates the RDS). In production, you'd bake this step and
# only do readRDS(); for dev this is handy.
# source("app_synthetic_inputs.R", local = TRUE)

# Pick which registry to expose in the app:
model_registry <- readRDS("models/model_registry_group.rds")

program_choices <- sort(unique(model_registry$program))
measure_choices <- c("Auto", "PHQ-9", "GAD-7")


# 2) Core prediction helper (now uses predict()) -------------------------------

predict_promis_change <- function(program, measure, baseline, followup) {
  stopifnot(length(program) == 1, length(measure) == 1)
  if (!measure %in% c("PHQ-9", "GAD-7")) stop("measure must be 'PHQ-9' or 'GAD-7'")
  
  row <- model_registry %>% dplyr::filter(program == !!program, measure == !!measure)
  if (nrow(row) != 1L) stop(glue("No model found for program '{program}' and measure '{measure}'."))
  
  change <- followup - baseline
  
  # newdata with reference covariates
  nd <- as.list(row$cov_ref[[1]])
  nd[[row$base_var]] <- baseline
  nd[[row$dvar]]     <- change
  newdata <- as.data.frame(nd, stringsAsFactors = TRUE)
  
  # compute 80/90/95 prediction intervals
  levels <- c(0.80, 0.90, 0.95)
  intervals <- purrr::map_dfr(levels, function(lv) {
    pr <- predict(row$model[[1]], newdata = newdata, interval = "prediction", level = lv)
    tibble::tibble(level = lv, fit = pr[1,"fit"], lwr = pr[1,"lwr"], upr = pr[1,"upr"])
  }) %>% dplyr::arrange(level)
  
  # keep 95% as the headline; also keep lo/hi for backward compat with your UI
  list(
    change     = change,
    yhat       = intervals$fit[intervals$level == 0.95],
    lo         = intervals$lwr[intervals$level == 0.95],
    hi         = intervals$upr[intervals$level == 0.95],
    intervals  = intervals,
    model_row  = row
  )
}

# Utility to infer which measure to use
infer_measure <- function(sel_measure, phq_has, gad_has) {
  if (sel_measure %in% c("PHQ-9", "GAD-7")) return(sel_measure)
  # Auto: prefer the one with complete inputs; if both, default to PHQ-9
  if (phq_has && !gad_has) return("PHQ-9")
  if (!phq_has && gad_has) return("GAD-7")
  if (phq_has && gad_has) return("PHQ-9")
  return(NA_character_)
}

# Input validators
is_complete_pair <- function(baseline, followup) {
  is.finite(baseline) && is.finite(followup)
}

# 3) UI ------------------------------------------------------------
ui <- page_fillable(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  title = "PROMIS Change Predictor",
  layout_column_wrap(
    width = 1/3,
    gap = "20px",
    card(
      header = "Program & Measure",
      card_body(
        selectInput("program", "Treatment Program", choices = program_choices, selected = program_choices[1]),
        radioButtons("measure_choice", "Which measure drives the prediction?", choices = measure_choices, selected = "Auto"),
        helpText("Choose 'Auto' to use whichever score has complete baseline & follow-up inputs (PHQ-9 preferred if both).")
      )
    ),
    card(
      header = "PHQ-9 Inputs (0–27)",
      card_body(
        numericInput("phq_base", "Baseline PHQ-9", value = NA, min = 0, max = 27, step = 1),
        numericInput("phq_follow", "Follow-up PHQ-9", value = NA, min = 0, max = 27, step = 1)
      )
    ),
    card(
      header = "GAD-7 Inputs (0–21)",
      card_body(
        numericInput("gad_base", "Baseline GAD-7", value = NA, min = 0, max = 21, step = 1),
        numericInput("gad_follow", "Follow-up GAD-7", value = NA, min = 0, max = 21, step = 1)
      )
    )
  ),
  layout_column_wrap(
    width = 1/2,
    gap = "20px",
    card(
      header = "Prediction",
      card_body(
        uiOutput("prediction_text"),
        plotOutput("pi_plot", height = "260px"),
        uiOutput("equation_text")
      )
    ),
    card(
      header = "Inputs & Model",
      card_body(
        tableOutput("inputs_table"),
        tableOutput("model_table")
      )
    )
  )
)


# 4) Server updates ------------------------------------------------------------

server <- function(input, output, session) {
  
  measure_active <- reactive({
    phq_ok <- is_complete_pair(input$phq_base, input$phq_follow)
    gad_ok <- is_complete_pair(input$gad_base, input$gad_follow)
    infer_measure(input$measure_choice, phq_ok, gad_ok)
  })
  
  pred <- reactive({
    req(input$program)
    m <- measure_active()
    validate(need(!is.na(m), "Enter both baseline and follow-up for PHQ-9 or GAD-7 (or select the other measure)."))
    
    if (m == "PHQ-9") {
      baseline <- input$phq_base; followup <- input$phq_follow
    } else {
      baseline <- input$gad_base; followup <- input$gad_follow
    }
    
    validate(
      need(isTRUE(baseline >= 0), "Baseline must be non-negative."),
      need(isTRUE(followup >= 0), "Follow-up must be non-negative."),
      need(any(model_registry$program == input$program & model_registry$measure == m),
           glue("No model available for {input$program} using {m}."))
    )
    
    res <- predict_promis_change(input$program, m, baseline, followup)
    list(
      measure   = m,
      baseline  = baseline,
      followup  = followup,
      change    = res$change,
      yhat      = res$yhat,
      lo        = res$lo,
      hi        = res$hi,
      intervals = res$intervals,   # <-- add this line
      model_row = res$model_row
    )
  })
  
  output$prediction_text <- renderUI({
    p  <- pred()
    iv <- p$intervals
    fmt <- function(x) sprintf("%.1f", x)
    HTML(glue(
      "<h4 style='margin-top:0;'>Predicted PROMIS change: <b>{fmt(p$yhat)}</b></h4>
     <div>Prediction intervals:</div>
     <ul style='margin-top:4px;'>
       <li>80%: <b>{fmt(iv$lwr[iv$level==0.80])} to {fmt(iv$upr[iv$level==0.80])}</b></li>
       <li>90%: <b>{fmt(iv$lwr[iv$level==0.90])} to {fmt(iv$upr[iv$level==0.90])}</b></li>
       <li>95%: <b>{fmt(iv$lwr[iv$level==0.95])} to {fmt(iv$upr[iv$level==0.95])}</b></li>
     </ul>
     <div style='color:#6c757d;'>Positive values indicate improvement.</div>"
    ))
  })
  
  
  # (optionally bump height in UI: plotOutput("pi_plot", height = "300px"))
  output$pi_plot <- renderPlot({
    p  <- pred()
    iv <- p$intervals %>%
      dplyr::mutate(
        level_lab = factor(paste0(level*100, "% PI"),
                           levels = c("95% PI","90% PI","80% PI")) # top→bottom
      )
    
    point_df <- tibble::tibble(
      level_lab = factor("Point", levels = c("95% PI","90% PI","80% PI","Point")),
      fit = as.numeric(p$yhat)
    )
    
    ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
      geom_errorbarh(
        data = iv,
        aes(y = level_lab, xmin = lwr, xmax = upr),
        height = 0.20, linewidth = 1
      ) +
      geom_point(
        data = point_df,
        aes(x = fit, y = level_lab),
        size = 3
      ) +
      labs(x = "Predicted ΔPROMIS (Follow-up − Baseline)", y = NULL,
           title = "Prediction Intervals (80/90/95%)") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.title.y = element_blank())
  })
  
  
  output$inputs_table <- renderTable({
    p <- pred()
    tibble::tibble(
      Program = input$program,
      Measure = p$measure,
      Baseline = p$baseline,
      Follow_up = p$followup,
      Change = p$change
    )
  })
  
  # Show the actual coefficients used (from the lm object)
  output$model_table <- renderTable({
    m <- measure_active(); req(m)
    row <- model_registry %>% filter(program == input$program, measure == m)
    req(nrow(row) == 1)
    co <- coef(row$model[[1]])
    tibble::tibble(
      Term = names(co),
      Coefficient = unname(co)
    )
  })
  
  # Equation text (pulls coefficients directly; uses your selected baseline/change names)
  output$equation_text <- renderUI({
    m <- measure_active(); req(m)
    row <- model_registry %>% filter(program == input$program, measure == m)
    req(nrow(row) == 1)
    co <- coef(row$model[[1]])
    b0 <- round(unname(co["(Intercept)"]), 3)
    bB <- round(unname(co[row$base_var]), 3)
    bD <- round(unname(co[row$dvar]), 3)
    HTML(glue(
      "<div style='margin-top:10px;color:#6c757d;'>
        Model: ΔPROMIS = {b0} + {bB}×{row$base_var} + {bD}×(Follow-up − Baseline)
        <br/>PI computed via lm prediction interval.
      </div>"
    ))
  })
}

shinyApp(ui, server)