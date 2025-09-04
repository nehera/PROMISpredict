library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(broom)

dir.create("models", showWarnings = FALSE, recursive = TRUE)

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

# helper to compute reference covariates for prediction
covariate_reference <- function(df, covars) {
  refs <- lapply(covars, function(v) {
    x <- df[[v]]
    if (is.numeric(x)) {
      mean(x, na.rm = TRUE)
    } else if (is.logical(x)) {
      as.logical(round(mean(as.integer(x), na.rm = TRUE)))
    } else {
      # modal level
      tab <- sort(table(x), decreasing = TRUE)
      nm <- names(tab)[1]
      if (is.factor(x)) factor(nm, levels = levels(x)) else nm
    }
  })
  setNames(refs, covars)
}

fit_long_adjusted <- function(data, measure = c("PHQ9","GAD7"), strata = NULL) {
  measure <- match.arg(measure)
  base_var <- paste0(measure, "_baseline")
  dvar     <- paste0("d", measure)
  
  rhs <- c(base_var, dvar, covars)
  base_fml <- paste("dPROMIS ~", paste(rhs, collapse = " + "))
  if (!is.null(strata)) {
    fml <- as.formula(paste0(base_fml, " + ", strata, " + ", dvar, ":", strata))
  } else {
    fml <- as.formula(base_fml)
  }
  fit <- lm(fml, data = data)
  
  list(
    fit    = fit,
    sigma  = rmse(fit),
    base_var = base_var,
    dvar     = dvar,
    cov_ref  = covariate_reference(data, covars)
  )
}

make_model_registry <- function(cs_data, long_data,
                                stratify = c("overall","type","group"),
                                measures = c("PHQ-9","GAD-7")) {
  stratify <- match.arg(stratify)
  groups <- switch(
    stratify,
    overall = "Overall",
    type    = sort(unique(long_data$SVC_DEPT_TYPE)),
    group   = sort(unique(long_data$SVC_DEPT_GROUP))
  )
  
  out <- list()
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
      long_fit <- fit_long_adjusted(
        long_df_use,
        measure = ifelse(m == "PHQ-9","PHQ9","GAD7"),
        strata  = NULL
      )
      
      out[[length(out)+1]] <- tibble::tibble(
        program  = lbl,
        measure  = m,
        model    = list(long_fit$fit),
        sigma    = long_fit$sigma,
        base_var = long_fit$base_var,
        dvar     = long_fit$dvar,
        cov_ref  = list(long_fit$cov_ref)
      )
    }
  }
  bind_rows(out)
}

# create registries (you can keep all three if useful)
model_registry_overall <- make_model_registry(cs_df, long_df, stratify = "overall")
model_registry_type    <- make_model_registry(cs_df, long_df, stratify = "type")
model_registry_group   <- make_model_registry(cs_df, long_df, stratify = "group")

# choose which one your app will use by default; save all for flexibility
saveRDS(model_registry_overall, file = "models/model_registry_overall.rds")
saveRDS(model_registry_type,    file = "models/model_registry_type.rds")
saveRDS(model_registry_group,   file = "models/model_registry_group.rds")