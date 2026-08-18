# ==============================================================================
# 06_setup_emotional_analysis.R
#
# Shared setup and helper functions for the emotional attention-bias analyses.
#
# Manuscript:
# Residential Early-Life Air Pollution Exposure and Oculomotor Markers of
# Cognitive Workload in Children: A Prospective Birth Cohort Study
#
# Cohort:
# ENVIRONAGE (ENVIRonmental influence ON early AGEing)
#
# Software:
# R version 4.4.1
#
# Data availability:
# Individual-level ENVIRONAGE data are not included in this repository.
# ==============================================================================


# ==============================================================================
# 1. Settings
# ==============================================================================

data_path  <- file.path("data", "analysis_dataset.xlsx")
sheet_name <- 1

results_dir <- file.path("results", "emotional_attention_bias")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05


# ==============================================================================
# 2. Required packages
# ==============================================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "purrr",
  "tibble",
  "tidyr",
  "poolr",
  "openxlsx"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  stop(
    "The following required packages are not installed: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))


# ==============================================================================
# 3. Load data
# ==============================================================================

if (!file.exists(data_path)) {
  stop(
    "Analysis dataset not found at: ", data_path,
    "\nIndividual-level ENVIRONAGE data are not distributed with this repository."
  )
}

dat <- readxl::read_excel(
  data_path,
  sheet = sheet_name
)

# Standardize whole-childhood variable names if spaces are present.
names(dat)[names(dat) == "PM25_Whole childhood"] <- "PM25_Whole_childhood"
names(dat)[names(dat) == "NO2_Whole childhood"]  <- "NO2_Whole_childhood"
names(dat)[names(dat) == "BC_Whole childhood"]   <- "BC_Whole_childhood"


# ==============================================================================
# 4. Variable definitions
# ==============================================================================

# Emotional attention-bias outcomes
outcomes <- c(
  "FixTD_biasHappy",
  "FixTD_biasAnger",
  "FixTD_biasFear"
)

# Early-life exposure hypothesis family
early_life_exposures <- c(
  "BC_trim1", "BC_trim2", "BC_trim3", "BC_whole_preg",
  "PM25_trim1", "PM25_trim2", "PM25_trim3", "PM25_whole_preg",
  "NO2_trim1", "NO2_trim2", "NO2_trim3", "NO2_whole_preg",
  "BC_first1000days", "PM25_first1000days", "NO2_first1000days"
)

# Recent exposure hypothesis family
recent_exposures <- c(
  "BC_Lastw", "BC_Lastm", "BC_last1y",
  "PM25_Lastw", "PM25_Lastm", "PM25_last1y",
  "NO2_Lastw", "NO2_Lastm", "NO2_last1y",
  "BC_urine"
)

all_exposures <- c(
  early_life_exposures,
  recent_exposures
)

# Primary model covariates
numeric_covariates <- c(
  "AgeM",
  "AgeChild"
)

factor_covariates <- c(
  "DiplH",
  "Smoke",
  "SeasonDel",
  "Parity",
  "Seks",
  "Ethn"
)

# ES_WsTD was evaluated during development of the analysis but was not included
# in the final manuscript models reproduced here.


# ==============================================================================
# 5. Data preparation
# ==============================================================================

numize <- function(x) {
  if (is.numeric(x)) return(x)

  x <- as.character(x)
  x[x == ""] <- NA_character_

  suppressWarnings(
    as.numeric(
      gsub(",", ".", x, fixed = TRUE)
    )
  )
}

numeric_variables <- intersect(
  c(
    all_exposures,
    outcomes,
    numeric_covariates
  ),
  names(dat)
)

dat[numeric_variables] <- lapply(
  dat[numeric_variables],
  numize
)

# Preserve source coding while treating categorical variables as factors.
for (v in intersect(factor_covariates, names(dat))) {
  dat[[v]] <- factor(dat[[v]])
}

# Urinary BC is already log10-transformed in the manuscript analysis dataset
# and is therefore not transformed again here.


# ==============================================================================
# 6. Helper functions
# ==============================================================================

exposure_family <- function(exposure) {

  if (exposure %in% early_life_exposures) {
    return("early_life")
  }

  if (exposure %in% recent_exposures) {
    return("recent")
  }

  NA_character_
}


# ==============================================================================
# 7. Li & Ji effective number of independent exposure tests
# ==============================================================================

calculate_meff <- function(df, exposures) {

  exposures <- exposures[
    exposures %in% names(df)
  ]

  if (length(exposures) < 2) {
    return(
      list(
        meff = length(exposures),
        variables = exposures,
        correlation_matrix = NULL
      )
    )
  }

  exposure_data <- df[
    ,
    exposures,
    drop = FALSE
  ]

  keep <- vapply(
    exposure_data,
    function(x) {

      x <- suppressWarnings(
        as.numeric(x)
      )

      sum(!is.na(x)) >= 10 &&
        is.finite(
          stats::sd(
            x,
            na.rm = TRUE
          )
        ) &&
        stats::sd(
          x,
          na.rm = TRUE
        ) > 0
    },
    logical(1)
  )

  exposure_data <- exposure_data[
    ,
    keep,
    drop = FALSE
  ]

  variables_used <- names(
    exposure_data
  )

  if (length(variables_used) < 2) {
    return(
      list(
        meff = length(variables_used),
        variables = variables_used,
        correlation_matrix = NULL
      )
    )
  }

  correlation_matrix <- stats::cor(
    exposure_data,
    use = "pairwise.complete.obs",
    method = "spearman"
  )

  correlation_matrix[
    !is.finite(correlation_matrix)
  ] <- 0

  diag(correlation_matrix) <- 1

  meff <- poolr::meff(
    correlation_matrix,
    method = "liji"
  )

  list(
    meff = as.numeric(meff),
    variables = variables_used,
    correlation_matrix = correlation_matrix
  )
}

meff_early <- calculate_meff(
  dat,
  early_life_exposures
)

meff_recent <- calculate_meff(
  dat,
  recent_exposures
)

get_meff <- function(family) {

  if (family == "early_life") {
    return(meff_early$meff)
  }

  if (family == "recent") {
    return(meff_recent$meff)
  }

  NA_real_
}


# ==============================================================================
# 8. Primary linear regression function
# ==============================================================================

fit_linear_model <- function(
    df,
    outcome,
    exposure,
    numeric_covars = numeric_covariates,
    factor_covars = factor_covariates
) {

  if (
    !outcome %in% names(df) ||
    !exposure %in% names(df)
  ) {
    return(tibble::tibble())
  }

  variables_needed <- unique(
    c(
      outcome,
      exposure,
      numeric_covars,
      factor_covars
    )
  )

  variables_present <- intersect(
    variables_needed,
    names(df)
  )

  analysis_data <- df |>
    dplyr::select(
      dplyr::all_of(
        variables_present
      )
    ) |>
    dplyr::filter(
      stats::complete.cases(.)
    )

  # Drop categorical covariates with only one observed level.
  for (
    v in intersect(
      factor_covars,
      names(analysis_data)
    )
  ) {

    analysis_data[[v]] <- droplevels(
      as.factor(
        analysis_data[[v]]
      )
    )

    if (
      nlevels(
        analysis_data[[v]]
      ) < 2
    ) {
      analysis_data[[v]] <- NULL
    }
  }

  n <- nrow(
    analysis_data
  )

  if (
    n < 10 ||
    !exposure %in% names(
      analysis_data
    )
  ) {

    return(
      tibble::tibble(
        exposure = exposure,
        outcome = outcome,
        n = n,
        estimate = NA_real_,
        standard_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        exposure_iqr = NA_real_,
        estimate_per_iqr = NA_real_,
        ci_low_per_iqr = NA_real_,
        ci_high_per_iqr = NA_real_,
        model_formula = NA_character_
      )
    )
  }

  rhs <- c(
    exposure,
    intersect(
      numeric_covars,
      names(analysis_data)
    ),
    intersect(
      factor_covars,
      names(analysis_data)
    )
  )

  formula <- stats::as.formula(
    paste(
      outcome,
      "~",
      paste(
        rhs,
        collapse = " + "
      )
    )
  )

  model <- stats::lm(
    formula,
    data = analysis_data
  )

  if (
    !exposure %in% names(
      stats::coef(model)
    )
  ) {

    return(
      tibble::tibble(
        exposure = exposure,
        outcome = outcome,
        n = n,
        estimate = NA_real_,
        standard_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        exposure_iqr = NA_real_,
        estimate_per_iqr = NA_real_,
        ci_low_per_iqr = NA_real_,
        ci_high_per_iqr = NA_real_,
        model_formula = paste(
          deparse(formula),
          collapse = " "
        )
      )
    )
  }

  coefficient_table <- summary(model)$coefficients

  beta <- unname(
    coefficient_table[
      exposure,
      "Estimate"
    ]
  )

  standard_error <- unname(
    coefficient_table[
      exposure,
      "Std. Error"
    ]
  )

  statistic <- unname(
    coefficient_table[
      exposure,
      "t value"
    ]
  )

  p_value <- unname(
    coefficient_table[
      exposure,
      "Pr(>|t|)"
    ]
  )

  critical_value <- stats::qt(
    0.975,
    df = model$df.residual
  )

  ci_low <- beta -
    critical_value *
    standard_error

  ci_high <- beta +
    critical_value *
    standard_error

  exposure_iqr <- stats::IQR(
    analysis_data[[exposure]],
    na.rm = TRUE
  )

  tibble::tibble(
    exposure = exposure,
    outcome = outcome,
    n = n,
    estimate = beta,
    standard_error = standard_error,
    statistic = statistic,
    p_value = p_value,
    ci_low = ci_low,
    ci_high = ci_high,
    exposure_iqr = exposure_iqr,
    estimate_per_iqr = beta * exposure_iqr,
    ci_low_per_iqr = ci_low * exposure_iqr,
    ci_high_per_iqr = ci_high * exposure_iqr,
    model_formula = paste(
      deparse(formula),
      collapse = " "
    )
  )
}
