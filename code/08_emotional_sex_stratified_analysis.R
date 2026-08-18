# ==============================================================================
# 08_emotional_sex_stratified_analysis.R
#
# Sex-stratified emotional attention-bias analyses.
#
# This script reproduces the analyses reported in Supplemental Tables 9 and 11.
# Models use the same specification as the primary analysis, except that child
# sex is omitted from the covariate set within each sex stratum.
# ==============================================================================

source(
  file.path(
    "code",
    "06_setup_emotional_analysis.R"
  )
)


# ==============================================================================
# 1. Check sex variable
# ==============================================================================

if (
  !"Seks" %in% names(dat)
) {
  stop(
    "The sex variable 'Seks' is not present in the analysis dataset."
  )
}

sex_levels <- levels(
  droplevels(
    dat$Seks
  )
)

factor_covariates_stratified <- setdiff(
  factor_covariates,
  "Seks"
)


# ==============================================================================
# 2. Run sex-stratified analyses
# ==============================================================================

sex_stratified_results <- purrr::map_dfr(
  sex_levels,
  function(sex_level) {

    dat_sex <- dat |>
      dplyr::filter(
        Seks == sex_level
      ) |>
      droplevels()

    purrr::map_dfr(
      all_exposures,
      function(exposure) {

        purrr::map_dfr(
          outcomes,
          function(outcome) {

            fit_linear_model(
              df = dat_sex,
              outcome = outcome,
              exposure = exposure,
              factor_covars =
                factor_covariates_stratified
            ) |>
              dplyr::mutate(
                sex_stratum = sex_level
              )
          }
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
    sex_stratum,
    hypothesis_family,
    outcome,
    exposure
  )


# ==============================================================================
# 3. Export
# ==============================================================================

output_file <- file.path(
  results_dir,
  "08_emotional_sex_stratified_analysis.xlsx"
)

openxlsx::write.xlsx(
  sex_stratified_results,
  file = output_file,
  overwrite = TRUE
)

message(
  "Sex-stratified emotional attention-bias analysis completed."
)

message(
  "Results written to: ",
  output_file
)
