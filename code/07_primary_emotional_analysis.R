# ==============================================================================
# 07_primary_emotional_analysis.R
#
# Primary linear regression analyses for emotional attention-bias outcomes.
#
# This script reproduces the analyses reported in Supplemental Tables 8 and 10.
# ==============================================================================

source(
  file.path(
    "code",
    "06_setup_emotional_analysis.R"
  )
)


# ==============================================================================
# 1. Primary models
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
) |>
  dplyr::mutate(
    hypothesis_family = vapply(
      exposure,
      exposure_family,
      FUN.VALUE = character(1)
    ),
    meff = vapply(
      hypothesis_family,
      get_meff,
      FUN.VALUE = numeric(1)
    ),
    corrected_alpha = alpha / meff,
    p_value_meff = pmin(
      p_value * meff,
      1
    ),
    significant_meff = p_value <
      corrected_alpha
  ) |>
  dplyr::arrange(
    hypothesis_family,
    outcome,
    exposure
  )


# ==============================================================================
# 2. Meff summary
# ==============================================================================

meff_summary <- tibble::tibble(
  hypothesis_family = c(
    "early_life",
    "recent"
  ),
  n_exposures = c(
    length(
      meff_early$variables
    ),
    length(
      meff_recent$variables
    )
  ),
  meff = c(
    meff_early$meff,
    meff_recent$meff
  ),
  corrected_alpha = c(
    alpha /
      meff_early$meff,
    alpha /
      meff_recent$meff
  )
)


# ==============================================================================
# 3. Export
# ==============================================================================

output_file <- file.path(
  results_dir,
  "07_primary_emotional_analysis.xlsx"
)

workbook <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  workbook,
  "Primary_results"
)

openxlsx::writeData(
  workbook,
  "Primary_results",
  primary_results
)

openxlsx::addWorksheet(
  workbook,
  "Meff_summary"
)

openxlsx::writeData(
  workbook,
  "Meff_summary",
  meff_summary
)

if (
  !is.null(
    meff_early$correlation_matrix
  )
) {

  openxlsx::addWorksheet(
    workbook,
    "Meff_corr_early"
  )

  openxlsx::writeData(
    workbook,
    "Meff_corr_early",
    as.data.frame(
      meff_early$correlation_matrix
    ),
    rowNames = TRUE
  )
}

if (
  !is.null(
    meff_recent$correlation_matrix
  )
) {

  openxlsx::addWorksheet(
    workbook,
    "Meff_corr_recent"
  )

  openxlsx::writeData(
    workbook,
    "Meff_corr_recent",
    as.data.frame(
      meff_recent$correlation_matrix
    ),
    rowNames = TRUE
  )
}

openxlsx::saveWorkbook(
  workbook,
  output_file,
  overwrite = TRUE
)

message(
  "Primary emotional attention-bias analysis completed."
)

message(
  "Results written to: ",
  output_file
)
