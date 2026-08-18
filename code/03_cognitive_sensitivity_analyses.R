# ==============================================================================
# 03_cognitive_sensitivity_analyses.R
#
# Sensitivity analyses for cognitive eye-tracking outcomes.
#
# This script reproduces the cognitive section of Supplemental Table 4.
# Sensitivity analyses are restricted to exposure-outcome associations with
# nominal p < 0.05 in the primary cognitive analysis.
#
# Analyses included:
#   1. Primary model
#   2. Additional adjustment for time of day of follow-up assessment
#   3. Additional adjustment for reported average sleep duration
#   4. Additional adjustment of early-life models for cumulative lifetime
#      exposure to the corresponding pollutant
#   5. Mutual adjustment of trimester-specific models for the two remaining
#      trimester-specific exposure variables, entered as separate covariates
# ==============================================================================

source(file.path("code", "00_setup_cognitive_analysis.R"))


# ==============================================================================
# 1. Identify primary cognitive associations with nominal p < 0.05
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

selected_pairs <- primary_results |>
  dplyr::filter(
    !is.na(p_value),
    p_value < 0.05
  ) |>
  dplyr::select(
    exposure,
    outcome
  ) |>
  dplyr::distinct()


# ==============================================================================
# 2. Sensitivity-model function
# ==============================================================================

run_sensitivity_model <- function(
    exposure,
    outcome,
    analysis_type
) {

  if (analysis_type == "primary") {

    result <- fit_linear_model(
      df = dat,
      outcome = outcome,
      exposure = exposure
    )

  } else if (analysis_type == "time_of_day") {

    result <- fit_linear_model(
      df = dat,
      outcome = outcome,
      exposure = exposure,
      extra_factor = time_of_day_covariate
    )

  } else if (analysis_type == "sleep_duration") {

    result <- fit_linear_model(
      df = dat,
      outcome = outcome,
      exposure = exposure,
      extra_numeric = sleep_covariate
    )

  } else if (analysis_type == "lifetime_exposure") {

    # Only applicable to early-life exposures.
    if (!exposure %in% early_life_exposures) {
      return(tibble::tibble())
    }

    adjuster <- lifetime_adjuster(exposure)

    if (length(adjuster) == 0) {
      return(tibble::tibble())
    }

    result <- fit_linear_model(
      df = dat,
      outcome = outcome,
      exposure = exposure,
      extra_numeric = adjuster
    )

  } else if (analysis_type == "trimester_specific") {

    # Only applicable to trimester-specific exposure windows.
    # The other two trimesters for the same pollutant are entered separately.
    adjusters <- trimester_adjusters(exposure)

    if (length(adjusters) == 0) {
      return(tibble::tibble())
    }

    result <- fit_linear_model(
      df = dat,
      outcome = outcome,
      exposure = exposure,
      extra_numeric = adjusters
    )

  } else {

    stop(
      "Unknown sensitivity-analysis type: ",
      analysis_type
    )
  }

  result |>
    dplyr::mutate(
      analysis_type = analysis_type
    )
}


# ==============================================================================
# 3. Run sensitivity analyses for selected primary associations
# ==============================================================================

analysis_types <- c(
  "primary",
  "time_of_day",
  "sleep_duration",
  "lifetime_exposure",
  "trimester_specific"
)

sensitivity_long <- purrr::pmap_dfr(
  selected_pairs,
  function(exposure, outcome) {

    purrr::map_dfr(
      analysis_types,
      function(analysis_type) {

        run_sensitivity_model(
          exposure = exposure,
          outcome = outcome,
          analysis_type = analysis_type
        )
      }
    )
  }
) |>
  dplyr::mutate(
    hypothesis_family = vapply(
      exposure,
      exposure_family,
      FUN.VALUE = character(1)
    )
  ) |>
  dplyr::arrange(
    outcome,
    exposure,
    factor(
      analysis_type,
      levels = analysis_types
    )
  )


# ==============================================================================
# 4. Wide table for comparison with Supplemental Table 4
# ==============================================================================

sensitivity_wide <- sensitivity_long |>
  dplyr::select(
    outcome,
    exposure,
    analysis_type,
    n,
    estimate_per_iqr,
    ci_low_per_iqr,
    ci_high_per_iqr,
    p_value
  ) |>
  tidyr::pivot_wider(
    names_from = analysis_type,
    values_from = c(
      n,
      estimate_per_iqr,
      ci_low_per_iqr,
      ci_high_per_iqr,
      p_value
    ),
    names_sep = "__"
  )


# ==============================================================================
# 5. Export
# ==============================================================================

output_file <- file.path(
  results_dir,
  "03_cognitive_sensitivity_analyses.xlsx"
)

workbook <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  workbook,
  "Sensitivity_long"
)
openxlsx::writeData(
  workbook,
  "Sensitivity_long",
  sensitivity_long
)

openxlsx::addWorksheet(
  workbook,
  "Sensitivity_wide"
)
openxlsx::writeData(
  workbook,
  "Sensitivity_wide",
  sensitivity_wide
)

openxlsx::saveWorkbook(
  workbook,
  output_file,
  overwrite = TRUE
)

message("Cognitive sensitivity analyses completed.")
message("Results written to: ", output_file)
