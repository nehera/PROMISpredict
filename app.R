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

# --------------------------------------------------------
# 1) Example model registry (REPLACE with your real models)
# --------------------------------------------------------
# Each row defines a linear prediction of PROMIS change:
#   ΔPROMIS = intercept + b_baseline * baseline + b_change * (followup - baseline)
#   PI: yhat ± 1.96 * sigma  (simple approximation)

source("app_synthetic_inputs.R", local = TRUE)

# choose which synthetic registry to expose in the app:
model_store <- model_store_type   # or model_store_overall / model_store_group

program_choices <- sort(unique(model_store$program))
measure_choices <- c("Auto", "PHQ-9", "GAD-7")


# --------------------------------------------------------
# 2) Core prediction helper
# --------------------------------------------------------
predict_promis_change <- function(program, measure, baseline, followup) {
  stopifnot(length(program) == 1, length(measure) == 1)
  if (!measure %in% c("PHQ-9", "GAD-7")) {
    stop("measure must be 'PHQ-9' or 'GAD-7'")
  }
  row <- model_store %>% filter(program == !!program, measure == !!measure)
  if (nrow(row) != 1L) {
    stop(glue("No model found for program '{program}' and measure '{measure}'."))
  }
  change <- followup - baseline
  yhat <- row$intercept + row$b_baseline * baseline + row$b_change * change
  # Simple normal-approx PI; replace with model-based PI if available
  lo <- as.numeric(yhat - 1.96 * row$sigma)
  hi <- as.numeric(yhat + 1.96 * row$sigma)
  list(
    yhat = as.numeric(yhat),
    lo = lo,
    hi = hi,
    change = change
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

# --------------------------------------------------------
# 3) UI
# --------------------------------------------------------
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

# --------------------------------------------------------
# 4) Server
# --------------------------------------------------------
server <- function(input, output, session) {
  
  # Which measure will we actually use?
  measure_active <- reactive({
    phq_ok <- is_complete_pair(input$phq_base, input$phq_follow)
    gad_ok <- is_complete_pair(input$gad_base, input$gad_follow)
    infer_measure(input$measure_choice, phq_ok, gad_ok)
  })
  
  # Reactive prediction
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
      need(isTRUE(followup >= 0), "Follow-up must be non-negative.")
    )
    
    # Ensure the selected program has a model for this measure
    validate(need(any(model_store$program == input$program & model_store$measure == m),
                  glue("No model available for {input$program} using {m}.")))
    
    res <- predict_promis_change(input$program, m, baseline, followup)
    list(
      measure = m,
      baseline = baseline,
      followup = followup,
      change = res$change,
      yhat = res$yhat,
      lo = res$lo,
      hi = res$hi
    )
  })
  
  # Prediction text
  output$prediction_text <- renderUI({
    p <- pred()
    HTML(glue(
      "<h4 style='margin-top:0;'>Predicted PROMIS change: <b>{sprintf('%.1f', p$yhat)}</b></h4>\n",
      "<div>95% PI: <b>", sprintf("%.1f to %.1f", p$lo, p$hi), "</b></div>",
      "<div style='color:#6c757d;'>Positive values indicate improvement.</div>"
    ))
  })
  
  # Plot with point & PI
  output$pi_plot <- renderPlot({
    p <- pred()
    df <- tibble::tibble(yhat = p$yhat, lo = p$lo, hi = p$hi)
    ggplot(df, aes(x = 1, y = yhat)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.05, linewidth = 1) +
      geom_point(size = 3) +
      scale_x_continuous(limits = c(0.8, 1.2)) +
      labs(x = NULL, y = "PROMIS Change", title = "Predicted Change with 95% PI") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
  })
  
  # Show inputs
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
  
  # Show model row used
  output$model_table <- renderTable({
    m <- measure_active()
    req(m)
    model_store %>%
      filter(program == input$program, measure == m) %>%
      rename(Intercept = intercept, Baseline_coef = b_baseline, Change_coef = b_change, Sigma = sigma)
  })
  
  # Equation text
  output$equation_text <- renderUI({
    m <- measure_active(); req(m)
    row <- model_store %>% filter(program == input$program, measure == m)
    req(nrow(row) == 1)
    HTML(glue(
      "<div style='margin-top:10px;color:#6c757d;'>",
      "Model: ΔPROMIS = {round(row$intercept,3)} + {round(row$b_baseline,3)}×Baseline + {round(row$b_change,3)}×(Follow-up − Baseline)",
      "<br/>PI ≈ ŷ ± 1.96×{round(row$sigma,2)}",
      "</div>"
    ))
  })
}

shinyApp(ui, server)
