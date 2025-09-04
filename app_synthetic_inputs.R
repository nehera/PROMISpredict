library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(broom)

# Synthetic data generator
generate_longitudinal_data <- function(n_individuals = 200, n_visits_range = c(2, 5), seed = 123) {
  set.seed(seed)
  MRN <- seq_len(n_individuals)
  df_list <- lapply(MRN, function(id) {
    n_visits <- sample(n_visits_range[1]:n_visits_range[2], 1)
    PHQ9 <- pmax(0, round(rnorm(n_visits, 10, 5)))
    GAD7 <- pmax(0, round(rnorm(n_visits, 8, 4)))
    base_PROMIS <- rnorm(1, 50, 10)
    delta_PROMIS <- rnorm(1, 0, 5)
    PROMIS_score <- base_PROMIS + 0.7 * PHQ9 + 0.5 * GAD7 + rnorm(n_visits, 0, 4) + seq_len(n_visits) * delta_PROMIS / n_visits
    time <- sort(sample(seq(0, 365, by = 3), n_visits))  # days
    
    data.frame(
      MRN = id,
      time = time,
      PHQ9 = PHQ9,
      GAD7 = GAD7,
      PROMIS = PROMIS_score,
      SVC_DEPT_TYPE = sample(c("Program", "Clinical"), 1),
      SVC_DEPT_GROUP = sample(LETTERS[1:5], 1),
      SEX = sample(c("Male", "Female"), 1),
      RACE = sample(c("White", "Black", "Asian", "Other"), 1),
      ETHNICITY = sample(c("Hispanic", "Non-Hispanic"), 1),
      AGE = rnorm(1, 45, 10),
      PEDIATRIC = sample(c(0, 1), 1),
      PROMIS_quarter = sample(1:4, n_visits, replace = TRUE)
    )
  })
  bind_rows(df_list)
}

set.seed(1)
synthetic_long <- generate_longitudinal_data(n_individuals = 1000)

# Pick each MRN's first and last visit as baseline/follow-up
panels <- synthetic_long %>%
  arrange(MRN, time) %>%
  group_by(MRN) %>%
  slice(c(1, n())) %>%
  mutate(visit = c("baseline","followup")) %>%
  ungroup() %>%
  select(MRN, visit, time, PROMIS, PHQ9, GAD7, SVC_DEPT_TYPE, SVC_DEPT_GROUP,
         SEX, RACE, ETHNICITY, AGE, PEDIATRIC, PROMIS_quarter)

wide <- panels %>%
  pivot_wider(
    id_cols = c(MRN, SVC_DEPT_TYPE, SVC_DEPT_GROUP, SEX, RACE, ETHNICITY, AGE, PEDIATRIC),
    names_from = visit,
    values_from = c(PROMIS, PHQ9, GAD7, time, PROMIS_quarter),
    names_sep = "_"
  ) %>%
  filter(!is.na(PROMIS_baseline), !is.na(PROMIS_followup)) %>%
  mutate(
    dPROMIS = PROMIS_followup - PROMIS_baseline,
    dPHQ9   = PHQ9_followup   - PHQ9_baseline,
    dGAD7   = GAD7_followup   - GAD7_baseline
  )

# Cross-sectional dataset = baseline row
cs_df <- panels %>%
  filter(visit == "baseline") %>%
  rename(
    PROMIS_cs = PROMIS, PHQ9_cs = PHQ9, GAD7_cs = GAD7,
    PROMIS_quarter_cs = PROMIS_quarter
  )

# Longitudinal dataset = change with baseline covariates joined
long_df <- wide %>%
  left_join(
    cs_df %>% select(MRN, PROMIS_quarter_cs),
    by = "MRN"
  )

# Covariates used in both models (baseline values)
covars <- c("SEX", "RACE", "ETHNICITY", "AGE", "PEDIATRIC", "PROMIS_quarter_cs")

# Utility: RMSE as predictive sigma
rmse <- function(model) sqrt(mean(residuals(model)^2, na.rm = TRUE))

# Cross-sectional (baseline PROMIS ~ reverse-coded PRO + covariates)
fit_cs_adjusted <- function(data, measure = c("PHQ9","GAD7"), strata = NULL) {
  measure <- match.arg(measure)
  pro_var <- paste0(measure, "_cs")
  # reverse-code so higher = better health
  d <- data %>%
    mutate(PRO_rev = max(.data[[pro_var]], na.rm = TRUE) - .data[[pro_var]])
  
  base_fml <- paste0("PROMIS_cs ~ PRO_rev + ", paste(covars, collapse = " + "))
  if (!is.null(strata)) {
    fml <- as.formula(paste0(base_fml, " + ", strata, " + PRO_rev:", strata))
  } else {
    fml <- as.formula(base_fml)
  }
  fit <- lm(fml, data = d)
  list(fit = fit, sigma = rmse(fit))
}

# Longitudinal (ΔPROMIS ~ baseline PRO + ΔPRO + covariates)
# NOTE: We DO NOT reverse-code here; instead we flip signs of coefficients on return
fit_long_adjusted <- function(data, measure = c("PHQ9","GAD7"), strata = NULL) {
  measure <- match.arg(measure)
  base_var <- paste0(measure, "_baseline")
  dvar     <- paste0("d", measure)
  
  # Build formula
  rhs <- c(base_var, dvar, covars)
  base_fml <- paste("dPROMIS ~", paste(rhs, collapse = " + "))
  if (!is.null(strata)) {
    fml <- as.formula(paste0(base_fml, " + ", strata, " + ", dvar, ":", strata))
  } else {
    fml <- as.formula(base_fml)
  }
  fit <- lm(fml, data = data)
  list(fit = fit, sigma = rmse(fit))
}

# Build per-program rows from a linear model
coef_to_row <- function(fit_obj, sigma, program_label, measure, base_var, dvar) {
  cf <- coef(fit_obj)
  intercept  <- unname(cf["(Intercept)"])
  b_base     <- unname(ifelse(base_var %in% names(cf), cf[base_var], 0))
  b_change   <- unname(ifelse(dvar     %in% names(cf), cf[dvar],     0))
  
  # Flip signs so that negative b_change => ΔPRO<0 (improvement) -> ΔPROMIS>0
  tibble::tibble(
    program    = program_label,
    measure    = measure,
    intercept  = intercept,
    b_baseline = +b_base,   # baseline PRO as-is (higher worse => usually positive; keep as fitted)
    b_change   = -b_change, # sign-flip for the app convention
    sigma      = sigma
  )
}

# Fit and assemble model_store
make_model_store <- function(cs_data, long_data,
                             stratify = c("overall","type","group"),
                             measures = c("PHQ-9","GAD-7")) {
  stratify <- match.arg(stratify)
  out <- list()
  
  # determine groups/program labels
  groups <- switch(
    stratify,
    overall = "Overall",
    type    = sort(unique(long_data$SVC_DEPT_TYPE)),
    group   = sort(unique(long_data$SVC_DEPT_GROUP))
  )
  
  for (grp in groups) {
    cs_df_use   <- cs_data
    long_df_use <- long_data
    lbl <- grp
    if (stratify == "type") {
      cs_df_use   <- cs_df_use   %>% filter(SVC_DEPT_TYPE == grp)
      long_df_use <- long_df_use %>% filter(SVC_DEPT_TYPE == grp)
    } else if (stratify == "group") {
      cs_df_use   <- cs_df_use   %>% filter(SVC_DEPT_GROUP == grp)
      long_df_use <- long_df_use %>% filter(SVC_DEPT_GROUP == grp)
    }
    
    for (m in measures) {
      # cross-sectional (adjusted) — not needed by the app, but returned invisibly if you want to inspect
      cs_fit <- fit_cs_adjusted(cs_df_use, measure = ifelse(m == "PHQ-9","PHQ9","GAD7"),
                                strata = NULL)
      # longitudinal (adjusted) — the one your app consumes
      long_fit <- fit_long_adjusted(
        long_df_use,
        measure = ifelse(m == "PHQ-9","PHQ9","GAD7"),
        strata = NULL
      )
      base_var <- paste0(ifelse(m=="PHQ-9","PHQ9","GAD7"), "_baseline")
      dvar     <- paste0("d", ifelse(m=="PHQ-9","PHQ9","GAD7"))
      
      row <- coef_to_row(
        fit_obj = long_fit$fit,
        sigma   = long_fit$sigma,
        program_label = lbl,
        measure = m,
        base_var = base_var,
        dvar     = dvar
      )
      out[[length(out) + 1]] <- row
    }
  }
  dplyr::bind_rows(out)
}

# Build three flavors (choose one for your app)
model_store_overall <- make_model_store(cs_df, long_df, stratify = "overall")
model_store_type    <- make_model_store(cs_df, long_df, stratify = "type")
model_store_group   <- make_model_store(cs_df, long_df, stratify = "group")

# Example: use type-level models in your app
model_store <- model_store_type %>%
  rename(program = program) %>%
  arrange(program, measure)

model_store
