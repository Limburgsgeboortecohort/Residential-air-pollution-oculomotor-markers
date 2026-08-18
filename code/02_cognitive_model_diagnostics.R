# ==============================================================================
# 02_cognitive_model_diagnostics.R
#
# Model diagnostics and linear-versus-natural-spline sensitivity analyses.
#
# This script:
#   1. Generates residuals-versus-fitted plots
#   2. Generates normal Q-Q plots
#   3. Compares linear models with natural spline models (3 df) for selected
#      primary exposure-outcome associations
# ==============================================================================

source(file.path("code", "00_setup_cognitive_analysis.R"))


# ==============================================================================
# 1. Re-run primary results to identify selected associations
# ==============================================================================

primary_results <- purrr::map_dfr(
  all_exposures,
  function(exposure) {
    purrr::map_dfr(
      outcomes,
      function(outcome) {
        fit_linear_model(
          df = dat,
          outcome = outcome,
          exposure = exposure
        )
      }
    )
  }
)


# ==============================================================================
# 2. Associations shown in Supplemental Figure 1 and Supplemental Table 1
# ==============================================================================

# These eight exposure-outcome combinations correspond exactly to the
# associations presented in Supplemental Figure 1 and Supplemental Table 1.

selected_associations <- tibble::tribble(
  ~exposure,             ~outcome,
  "BC_trim1",            "ERc_SacV",
  "BC_first1000days",    "ERc_SacV",
  "BC_Lastm",            "ERc_SacV",
  "NO2_Lastm",           "ERc_SacV",
  "NO2_Lastm",           "ERc_SacN",
  "PM25_trim1",          "ERc_SacV",
  "PM25_Lastm",          "ERc_SacV",
  "PM25_Lastm",          "ERc_SacN"
)

diagnostic_dir <- file.path(
  results_dir,
  "model_diagnostics"
)

dir.create(
  diagnostic_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==============================================================================
# 3. Residuals-versus-fitted and normal Q-Q plots
# ==============================================================================

make_diagnostic_plots <- function(exposure, outcome) {

  fitted <- fit_lm_for_diagnostics(
    df = dat,
    outcome = outcome,
    exposure = exposure
  )

  model <- fitted$model

  residual_data <- tibble::tibble(
    fitted = stats::fitted(model),
    residual = stats::residuals(model),
    standardized_residual = stats::rstandard(model)
  )

  # Residuals versus fitted
  p_resid <- ggplot2::ggplot(
    residual_data,
    ggplot2::aes(
      x = fitted,
      y = residual
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = FALSE
    ) +
    ggplot2::labs(
      x = "Fitted values",
      y = "Residuals"
    ) +
    ggplot2::theme_minimal(base_size = 14)

  ggplot2::ggsave(
    filename = file.path(
      diagnostic_dir,
      paste0(
        "residuals_vs_fitted__",
        outcome,
        "__",
        exposure,
        ".png"
      )
    ),
    plot = p_resid,
    width = 7,
    height = 5,
    dpi = 300
  )

  # Normal Q-Q plot
  p_qq <- ggplot2::ggplot(
    residual_data,
    ggplot2::aes(
      sample = standardized_residual
    )
  ) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line() +
    ggplot2::labs(
      x = "Theoretical quantiles",
      y = "Standardized residuals"
    ) +
    ggplot2::theme_minimal(base_size = 14)

  ggplot2::ggsave(
    filename = file.path(
      diagnostic_dir,
      paste0(
        "normal_qq__",
        outcome,
        "__",
        exposure,
        ".png"
      )
    ),
    plot = p_qq,
    width = 7,
    height = 5,
    dpi = 300
  )

  invisible(NULL)
}

if (nrow(selected_associations) > 0) {
  purrr::pwalk(
    selected_associations,
    make_diagnostic_plots
  )
}


# ==============================================================================
# 4. Linear versus natural spline models
# ==============================================================================

compare_linear_spline <- function(exposure, outcome) {

  fitted <- fit_lm_for_diagnostics(
    df = dat,
    outcome = outcome,
    exposure = exposure
  )

  analysis_data <- fitted$data
  linear_formula <- fitted$formula
  linear_model <- fitted$model

  numeric_use <- unique(c(
    numeric_covariates,
    if (
      outcome == "ERc_PupD" &&
      pupil_baseline_covariate %in% names(analysis_data)
    ) {
      pupil_baseline_covariate
    }
  ))

  rhs_covariates <- c(
    intersect(
      numeric_use,
      names(analysis_data)
    ),
    intersect(
      factor_covariates,
      names(analysis_data)
    )
  )

  spline_rhs <- c(
    paste0(
      "splines::ns(",
      exposure,
      ", df = 3)"
    ),
    rhs_covariates
  )

  spline_formula <- stats::as.formula(
    paste(
      outcome,
      "~",
      paste(
        spline_rhs,
        collapse = " + "
      )
    )
  )

  spline_model <- stats::lm(
    spline_formula,
    data = analysis_data
  )

  comparison <- stats::anova(
    linear_model,
    spline_model,
    test = "F"
  )

  tibble::tibble(
    exposure = exposure,
    outcome = outcome,
    n = nrow(analysis_data),
    p_nonlinearity = comparison$`Pr(>F)`[2],
    aic_linear = stats::AIC(linear_model),
    aic_spline = stats::AIC(spline_model),
    linear_formula = paste(
      deparse(linear_formula),
      collapse = " "
    ),
    spline_formula = paste(
      deparse(spline_formula),
      collapse = " "
    )
  )
}

spline_results <- tibble::tibble()

if (nrow(selected_associations) > 0) {
  spline_results <- purrr::map2_dfr(
    selected_associations$exposure,
    selected_associations$outcome,
    compare_linear_spline
  )
}


# ==============================================================================
# 5. Export
# ==============================================================================

output_file <- file.path(
  results_dir,
  "02_cognitive_model_diagnostics.xlsx"
)

workbook <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  workbook,
  "Selected_associations"
)
openxlsx::writeData(
  workbook,
  "Selected_associations",
  selected_associations
)

openxlsx::addWorksheet(
  workbook,
  "Spline_comparison"
)
openxlsx::writeData(
  workbook,
  "Spline_comparison",
  spline_results
)

openxlsx::saveWorkbook(
  workbook,
  output_file,
  overwrite = TRUE
)

message("Cognitive model diagnostics completed.")
message("Results written to: ", output_file)
message("Diagnostic plots written to: ", diagnostic_dir)
