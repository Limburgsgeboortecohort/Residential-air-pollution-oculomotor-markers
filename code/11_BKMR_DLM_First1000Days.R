# ==============================================================================
# BKMR-DLM analysis: first 1000 days from conception
# ==============================================================================
#
# Purpose
# -------
# R code used for the first-1000-days BKMR-DLM analyses reported in the
# manuscript.
#
# Exposure time scale
# -------------------
# Exposure histories are ordered from conception onward:
# day0   = conception
# day1   = 1 day after conception
# ...
# day999 = approximately 24 months after birth
#
# Developmental exposure windows
# ------------------------------
# trimester1  = days 0-89
# trimester2  = days 90-179
# trimester3  = days 180-269
# months0_6   = days 270-449
# months6_12  = days 450-629
# months12_18 = days 630-809
# months18_24 = days 810-989
#
# Main model specification
# ------------------------
# Lag basis: natural spline ("ns"), df = 3
# Kernel: polynomial
# Polynomial degree: 2
# MCMC iterations: 10,000
# Burn-in: 5,000
#
# Covariate adjustment
# --------------------
# All outcomes: maternal age, child age at assessment, conception-season
# sine/cosine terms, maternal education, smoking during pregnancy, parity,
# child sex, and ethnicity.
# ERc_PupD only: additionally adjusted for initial pupil size.
#
# Random-number generation
# ------------------------
# A fixed RNG seed can optionally be specified below for reproducible reruns.
#
# Data requirements
# -----------------
# This script assumes that the following protected study objects have already
# been created by the study's data-preparation workflow:
#   - pm25_complete
#   - no2_complete
#   - bc_complete
#   - Final_dataset_eye_tracking_study
#   - mosos_sub (must contain ID and Concep14d)
#
# Individual-level study data are not included in the public repository because
# of privacy and ethical restrictions.
#
# Outputs are written to the current working directory. Some generated output
# files contain participant-level information and must not be uploaded publicly.
# ==============================================================================

library(regimes)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)

# ------------------------------------------------------------------------------
# REQUIRED INPUT OBJECTS
# ------------------------------------------------------------------------------

required_objects <- c(
  "pm25_complete",
  "no2_complete",
  "bc_complete",
  "Final_dataset_eye_tracking_study",
  "mosos_sub"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Required input object(s) not found: ",
    paste(missing_objects, collapse = ", "),
    ". These protected study inputs are not included in the public repository."
  )
}

if (!all(c("ID", "Concep14d") %in% names(mosos_sub))) {
  stop("mosos_sub must contain ID and Concep14d.")
}

# ------------------------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------------------------

outcomes <- c(
  "ERc_FixN",
  "ERc_SacN",
  "ERc_PupD",
  "ERc_SacV",
  "FixTD_biasHappy",
  "FixTD_biasAnger",
  "FixTD_biasFear"
)

base_num_covs <- c(
  "AgeM",
  "AgeChild",
  "sin_conception_season",
  "cos_conception_season"
)

fac_covs <- c(
  "DiplM",
  "Smoke",
  "Parity",
  "Seks",
  "Ethn"
)

# Pupil-size correction used for ERc_PupD only.
# The first-1000-days pupil analysis used PupD_Cor. PupD_Corr is accepted as an
# alternative spelling if present in the protected analytic dataset.
if ("PupD_Cor" %in% names(Final_dataset_eye_tracking_study)) {
  pupil_covariate <- "PupD_Cor"
} else if ("PupD_Corr" %in% names(Final_dataset_eye_tracking_study)) {
  pupil_covariate <- "PupD_Corr"
} else {
  stop(
    "Neither PupD_Cor nor PupD_Corr exists in Final_dataset_eye_tracking_study. ",
    "ERc_PupD requires adjustment for initial pupil size."
  )
}

all_num_covs <- unique(c(base_num_covs, pupil_covariate))

use_gaussian_kernel <- FALSE
poly_degree <- 2

# Optional RNG seed for reproducible reruns
rng_seed <- NULL   # e.g., 12345

if (!is.null(rng_seed)) {
  set.seed(rng_seed)
}

# ------------------------------------------------------------
# Helper: mean if at least one daily value exists
# ------------------------------------------------------------

mean_if_any_observed <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# ------------------------------------------------------------
# Create first-1000-days exposure windows
# ------------------------------------------------------------

make_windows <- function(dat, n_days = 1000) {
  
  dat$ID <- as.character(dat$ID)
  
  day_cols <- grep("^X_day[0-9]+$", names(dat), value = TRUE)
  
  if (length(day_cols) == 0) {
    stop("No exposure columns found. Expected names like X_day0, X_day1, ...")
  }
  
  day_num <- as.numeric(sub("^X_day", "", day_cols))
  day_cols <- day_cols[order(day_num)]
  
  xdat <- dat[, day_cols, drop = FALSE]
  
  daily <- t(apply(xdat, 1, function(x) {
    
    x <- as.numeric(x)
    valid_idx <- which(!is.na(x))
    
    out <- rep(NA_real_, n_days)
    
    if (length(valid_idx) > 0) {
      
      conception_idx <- max(valid_idx)
      first_idx <- max(1, conception_idx - n_days + 1)
      
      vals <- rev(x[first_idx:conception_idx])
      
      out[seq_along(vals)] <- vals
    }
    
    out
  }))
  
  colnames(daily) <- paste0("day", 0:(n_days - 1))
  rownames(daily) <- dat$ID
  
  windows <- list(
    trimester1   = 0:89,
    trimester2   = 90:179,
    trimester3   = 180:269,
    months0_6    = 270:449,
    months6_12   = 450:629,
    months12_18  = 630:809,
    months18_24  = 810:989
  )
  
  window_mat <- data.frame(ID = dat$ID)
  
  for (w in names(windows)) {
    cols <- paste0("day", windows[[w]])
    
    window_mat[[w]] <- apply(
      daily[, cols, drop = FALSE],
      1,
      mean_if_any_observed
    )
  }
  
  list(
    daily_matrix = daily,
    window_matrix = window_mat
  )
}

# ------------------------------------------------------------
# Create exposure windows
# ------------------------------------------------------------

pm25_exp <- make_windows(pm25_complete)
no2_exp  <- make_windows(no2_complete)
bc_exp   <- make_windows(bc_complete)

pm25_windows <- pm25_exp$window_matrix
no2_windows  <- no2_exp$window_matrix
bc_windows   <- bc_exp$window_matrix

cat("Missing windows PM25:\n")
print(colSums(is.na(pm25_windows)))

cat("Missing windows NO2:\n")
print(colSums(is.na(no2_windows)))

cat("Missing windows BC:\n")
print(colSums(is.na(bc_windows)))

# ------------------------------------------------------------
# Align exposures to analytic dataset
# ------------------------------------------------------------

Final_dataset_eye_tracking_study$ID <- as.character(Final_dataset_eye_tracking_study$ID)

pm25_windows$ID <- as.character(pm25_windows$ID)
no2_windows$ID  <- as.character(no2_windows$ID)
bc_windows$ID   <- as.character(bc_windows$ID)

ID_vector <- Final_dataset_eye_tracking_study$ID

pm25_sub <- pm25_windows[
  match(ID_vector, pm25_windows$ID),
  -1,
  drop = FALSE
]

no2_sub <- no2_windows[
  match(ID_vector, no2_windows$ID),
  -1,
  drop = FALSE
]

bc_sub <- bc_windows[
  match(ID_vector, bc_windows$ID),
  -1,
  drop = FALSE
]

rownames(pm25_sub) <- ID_vector
rownames(no2_sub)  <- ID_vector
rownames(bc_sub)   <- ID_vector

pm25_sub <- as.matrix(pm25_sub)
no2_sub  <- as.matrix(no2_sub)
bc_sub   <- as.matrix(bc_sub)

# ------------------------------------------------------------
# Add conception-season covariates
# ------------------------------------------------------------

if (!all(c("ID", "Concep14d") %in% names(mosos_sub))) {
  stop("mosos_sub must contain ID and Concep14d.")
}

conception_dates_sub <- mosos_sub[
  match(ID_vector, mosos_sub$ID),
]

stopifnot(all(ID_vector == conception_dates_sub$ID))

conception_dates_sub$Concep14d <- as.Date(conception_dates_sub$Concep14d)

conception_doy <- as.numeric(
  format(conception_dates_sub$Concep14d, "%j")
)

Final_dataset_eye_tracking_study$sin_conception_season <- sin(
  2 * pi * conception_doy / 365.25
)

Final_dataset_eye_tracking_study$cos_conception_season <- cos(
  2 * pi * conception_doy / 365.25
)

# ------------------------------------------------------------
# Base complete-case restriction:
# exposures + covariates only, not all outcomes
# ------------------------------------------------------------

dat0 <- Final_dataset_eye_tracking_study

dat0[all_num_covs] <- lapply(dat0[all_num_covs], as.numeric)
dat0[fac_covs] <- lapply(dat0[fac_covs], as.factor)

complete_base_idx <-
  complete.cases(pm25_sub) &
  complete.cases(no2_sub) &
  complete.cases(bc_sub) &
  complete.cases(dat0[, c(base_num_covs, fac_covs)])

cat("Total N before restriction:", nrow(dat0), "\n")
cat("Excluded due to exposure/covariate missingness:", sum(!complete_base_idx), "\n")
cat("N available before outcome-specific exclusion:", sum(complete_base_idx), "\n")

PM25_base <- pm25_sub[complete_base_idx, , drop = FALSE]
NO2_base  <- no2_sub[complete_base_idx, , drop = FALSE]
BC_base   <- bc_sub[complete_base_idx, , drop = FALSE]

dat_base <- dat0[complete_base_idx, ]

stopifnot(all(rownames(PM25_base) == dat_base$ID))
stopifnot(all(rownames(NO2_base) == dat_base$ID))
stopifnot(all(rownames(BC_base) == dat_base$ID))

exp.dat_base <- list(
  pm25 = PM25_base,
  no2  = NO2_base,
  bc   = BC_base
)

lapply(exp.dat_base, dim)
sapply(exp.dat_base, function(x) sum(is.na(x)))
colSums(is.na(dat_base[, c(base_num_covs, fac_covs)]))

# ------------------------------------------------------------
# Diagnostics: excluded children
# ------------------------------------------------------------

id_check <- data.frame(
  ID = ID_vector,
  in_pm25 = ID_vector %in% pm25_windows$ID,
  in_no2  = ID_vector %in% no2_windows$ID,
  in_bc   = ID_vector %in% bc_windows$ID
)

missing_exposure_by_window <- data.frame(
  variable = c(
    paste0("pm25_", colnames(pm25_sub)),
    paste0("no2_", colnames(no2_sub)),
    paste0("bc_", colnames(bc_sub))
  ),
  n_missing = c(
    colSums(is.na(pm25_sub)),
    colSums(is.na(no2_sub)),
    colSums(is.na(bc_sub))
  )
)

missing_detail <- data.frame(ID = ID_vector)

for (v in colnames(pm25_sub)) {
  missing_detail[[paste0("pm25_", v)]] <- is.na(pm25_sub[, v])
}

for (v in colnames(no2_sub)) {
  missing_detail[[paste0("no2_", v)]] <- is.na(no2_sub[, v])
}

for (v in colnames(bc_sub)) {
  missing_detail[[paste0("bc_", v)]] <- is.na(bc_sub[, v])
}

missing_detail$n_missing_exposure_windows <- rowSums(missing_detail[, -1])

missing_detail_excluded <- missing_detail %>%
  filter(n_missing_exposure_windows > 0)

missing_detail_long <- missing_detail_excluded %>%
  pivot_longer(
    cols = -c(ID, n_missing_exposure_windows),
    names_to = "window",
    values_to = "missing"
  ) %>%
  filter(missing) %>%
  arrange(ID, window)

write.xlsx(
  list(
    id_check = id_check,
    missing_exposure_by_window = missing_exposure_by_window,
    missing_detail_excluded = missing_detail_excluded,
    missing_detail_long = missing_detail_long
  ),
  file = "diagnose_first1000_exposure_exclusions.xlsx",
  overwrite = TRUE
)

# ------------------------------------------------------------
# Fit BKMR-DLM models with outcome-specific N
# Polynomial kernel
# ------------------------------------------------------------

fits <- list()
summaries <- list()
model_N <- data.frame()

for (outcome in outcomes) {

  cat("\nRunning outcome:", outcome, "\n")

  num_covs_out <- base_num_covs

  if (outcome == "ERc_PupD") {
    num_covs_out <- c(base_num_covs, pupil_covariate)
    cat("Additional adjustment for ERc_PupD:", pupil_covariate, "\n")
  }

  y_raw <- as.numeric(dat_base[[outcome]])

  keep <- complete.cases(
    y_raw,
    dat_base[, c(num_covs_out, fac_covs)]
  )

  cat("N used:", sum(keep), "\n")
  
  dat_out <- dat_base[keep, ]
  y_out <- y_raw[keep]
  
  exp.dat_out <- lapply(exp.dat_base, function(x) {
    x[keep, , drop = FALSE]
  })
  
  cov.dat_out <- model.matrix(
    as.formula(
      paste("~", paste(c(num_covs_out, fac_covs), collapse = " + "))
    ),
    data = dat_out
  )[, -1, drop = FALSE]
  
  stopifnot(length(y_out) == nrow(cov.dat_out))
  stopifnot(all(sapply(exp.dat_out, nrow) == length(y_out)))
  stopifnot(sum(is.na(y_out)) == 0)
  stopifnot(sum(is.na(cov.dat_out)) == 0)
  stopifnot(all(sapply(exp.dat_out, function(x) sum(is.na(x))) == 0))
  
  model_N <- rbind(
    model_N,
    data.frame(
      outcome = outcome,
      N_used = length(y_out)
    )
  )
  
  fit <- bkmrdlm(
    y = y_out,
    x = exp.dat_out,
    z = cov.dat_out,
    niter = 10000,
    nburn = 5000,
    basis.opts = list(type = "ns", df = 3),
    gaussian = use_gaussian_kernel,
    polydegree = poly_degree
  )
  
  fits[[outcome]] <- fit
  summaries[[outcome]] <- summary(fit)
  
  save(
    fit,
    file = paste0(
      "bkmrdlm_first1000days_looseWindows_df3_poly",
      poly_degree,
      "_",
      outcome,
      ".RData"
    )
  )
}

print(model_N)

fits_first1000 <- fits
exp.dat_first1000_base <- exp.dat_base
dat_base_first1000 <- dat_base

# ------------------------------------------------------------
# Extract lag-window weights
# ------------------------------------------------------------

extract_weights <- function(fit, outcome_name) {
  
  weights <- as.data.frame(summary(fit)$weights)
  
  window_labels <- colnames(fit$x[[1]])
  
  weights$pollutant <- names(fit$x)[weights$m]
  weights$window <- window_labels[weights$t]
  weights$outcome <- outcome_name
  
  weights %>%
    select(
      outcome,
      pollutant,
      m,
      window,
      t,
      mean,
      lower,
      upper
    )
}

weights_all <- bind_rows(
  lapply(names(fits), function(outcome) {
    extract_weights(
      fit = fits[[outcome]],
      outcome_name = outcome
    )
  })
)

weights_all <- weights_all %>%
  mutate(
    window = factor(
      window,
      levels = c(
        "trimester1",
        "trimester2",
        "trimester3",
        "months0_6",
        "months6_12",
        "months12_18",
        "months18_24"
      ),
      labels = c(
        "Trimester 1",
        "Trimester 2",
        "Trimester 3",
        "0–6 months",
        "6–12 months",
        "12–18 months",
        "18–24 months"
      )
    ),
    excludes_zero = lower > 0 | upper < 0
  )

weights_summary <- weights_all %>%
  group_by(outcome, pollutant) %>%
  summarise(
    mean_abs_weight = mean(abs(mean), na.rm = TRUE),
    max_abs_weight = max(abs(mean), na.rm = TRUE),
    window_of_max_abs_weight = as.character(window[which.max(abs(mean))]),
    .groups = "drop"
  ) %>%
  arrange(outcome, desc(mean_abs_weight))

critical_windows <- weights_all %>%
  filter(excludes_zero) %>%
  arrange(outcome, pollutant, t)

write.xlsx(
  list(
    model_N = model_N,
    weights_all = weights_all,
    weights_summary = weights_summary,
    critical_windows = critical_windows
  ),
  file = paste0(
    "bkmrdlm_first1000days_looseWindows_weights_outcomeSpecificN_poly",
    poly_degree,
    ".xlsx"
  ),
  overwrite = TRUE
)

write.xlsx(
  list(
    PM25_first1000_base = as.data.frame(PM25_base),
    NO2_first1000_base  = as.data.frame(NO2_base),
    BC_first1000_base   = as.data.frame(BC_base),
    analytic_dataset_base = dat_base
  ),
  file = "first1000days_exposure_windows_base_sample.xlsx",
  rowNames = TRUE,
  overwrite = TRUE
)

# ------------------------------------------------------------
# Plot weight functions
# ------------------------------------------------------------

for (outcome in names(fits)) {
  
  plot_dat <- weights_all %>%
    filter(outcome == !!outcome)
  
  p <- ggplot(plot_dat, aes(x = window, y = mean, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.20
    ) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ pollutant, scales = "free_y") +
    labs(
      x = "Developmental exposure window",
      y = "Estimated lag weight"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = paste0(
      "bkmrdlm_first1000days_looseWindows_weights_poly",
      poly_degree,
      "_",
      outcome,
      ".png"
    ),
    plot = p,
    width = 9,
    height = 5,
    dpi = 600
  )
}

# ============================================================
# FIRST 1000 DAYS: MIXTURE EFFECTS
# Uses exp.dat_base, not old exp.dat
# Polynomial kernel prediction is consistent with model fitting
# ============================================================

stopifnot(exists("fits"))
stopifnot(exists("exp.dat_base"))
stopifnot(length(fits) > 0)
stopifnot(length(exp.dat_base) > 0)

# ------------------------------------------------------------
# Prediction helper for h(x)
# ------------------------------------------------------------
# Custom posterior prediction helper retained from the accepted analysis.

predict_bkmrdlm_h_new <- function(
    fit,
    xnew,
    gaussian_kernel = FALSE,
    polydegree = 2
) {
  
  M <- length(fit$x)
  n <- length(fit$y)
  n_iter <- nrow(fit$beta)
  n_new <- nrow(xnew[[1]])
  
  xnew_std <- vector("list", M)
  
  for (m in seq_len(M)) {
    xnew_std[[m]] <- as.matrix(
      (xnew[[m]] - fit$xscale$mean[m]) / fit$xscale$sd[m]
    )
  }
  
  weights <- vector("list", M)
  Xtheta_old <- matrix(NA_real_, n, M)
  Xtheta_new <- matrix(NA_real_, n_new, M)
  h_post <- matrix(NA_real_, n_iter, n_new)
  
  for (m in seq_len(M)) {
    weights[[m]] <-
      fit$theta[[m]] %*%
      t(fit$basis[[m]]$psi) /
      sqrt(nrow(fit$basis[[m]]$psi)) /
      sqrt(fit$rho[, m])
  }
  
  for (s in seq_len(n_iter)) {
    
    r <- as.numeric(fit$y - fit$z %*% fit$beta[s, ])
    
    for (m in seq_len(M)) {
      Xtheta_old[, m] <- fit$x[[m]] %*% weights[[m]][s, ]
      Xtheta_new[, m] <- xnew_std[[m]] %*% weights[[m]][s, ]
    }
    
    K_old <- matrix(NA_real_, n, n)
    K_newold <- matrix(NA_real_, n_new, n)
    
    if (gaussian_kernel) {
      
      for (i in seq_len(n)) {
        K_old[i, ] <- fit$tau2[s] *
          exp(-rowSums(t(t(Xtheta_old) - Xtheta_old[i, ]))^2)
      }
      
      for (i in seq_len(n_new)) {
        K_newold[i, ] <- fit$tau2[s] *
          exp(-rowSums(t(t(Xtheta_old) - Xtheta_new[i, ]))^2)
      }
      
    } else {
      
      for (i in seq_len(n)) {
        K_old[i, ] <- fit$tau2[s] *
          (1 + Xtheta_old %*% Xtheta_old[i, ])^polydegree
      }
      
      for (i in seq_len(n_new)) {
        K_newold[i, ] <- fit$tau2[s] *
          (1 + Xtheta_old %*% Xtheta_new[i, ])^polydegree
      }
    }
    
    diag(K_old) <- diag(K_old) + 1 + 1e-8
    
    K_old_inv <- tryCatch(
      chol2inv(chol(K_old)),
      error = function(e) {
        chol2inv(chol(K_old + diag(1e-6, nrow(K_old))))
      }
    )
    
    h_post[s, ] <- as.numeric(K_newold %*% K_old_inv %*% r)
  }
  
  h_post
}

summarise_contrast <- function(post_diff) {
  data.frame(
    mean = mean(post_diff, na.rm = TRUE),
    sd = sd(post_diff, na.rm = TRUE),
    lower = as.numeric(quantile(post_diff, 0.025, na.rm = TRUE)),
    upper = as.numeric(quantile(post_diff, 0.975, na.rm = TRUE)),
    pr_gt_0 = mean(post_diff > 0, na.rm = TRUE),
    significant = !(
      quantile(post_diff, 0.025, na.rm = TRUE) <= 0 &
        quantile(post_diff, 0.975, na.rm = TRUE) >= 0
    )
  )
}

make_profile <- function(exp.dat, q = 0.5) {
  lapply(exp.dat, function(x) {
    matrix(
      apply(x, 2, quantile, probs = q, na.rm = TRUE),
      nrow = 1
    )
  })
}

window_names <- colnames(exp.dat_base[[1]])

window_label_lookup <- c(
  trimester1 = "Trimester 1",
  trimester2 = "Trimester 2",
  trimester3 = "Trimester 3",
  months0_6 = "0–6 months",
  months6_12 = "6–12 months",
  months12_18 = "12–18 months",
  months18_24 = "18–24 months"
)

# ------------------------------------------------------------
# A. Window-specific joint mixture effects
# ------------------------------------------------------------

window_mixture_effects <- list()

for (outcome in names(fits)) {
  
  cat("\nWindow-specific joint mixture effects:", outcome, "\n")
  
  fit <- fits[[outcome]]
  x_med <- make_profile(exp.dat_base, q = 0.50)
  
  out_list <- list()
  
  for (w in seq_along(window_names)) {
    
    x_low <- x_med
    x_high <- x_med
    
    for (pollutant in names(exp.dat_base)) {
      
      x_low[[pollutant]][1, w] <-
        quantile(exp.dat_base[[pollutant]][, w], 0.25, na.rm = TRUE)
      
      x_high[[pollutant]][1, w] <-
        quantile(exp.dat_base[[pollutant]][, w], 0.75, na.rm = TRUE)
    }
    
    h_low <- predict_bkmrdlm_h_new(
      fit = fit,
      xnew = x_low,
      gaussian_kernel = use_gaussian_kernel,
      polydegree = poly_degree
    )
    
    h_high <- predict_bkmrdlm_h_new(
      fit = fit,
      xnew = x_high,
      gaussian_kernel = use_gaussian_kernel,
      polydegree = poly_degree
    )
    
    diff_post <- as.numeric(h_high - h_low)
    
    window_label <- ifelse(
      window_names[w] %in% names(window_label_lookup),
      window_label_lookup[[window_names[w]]],
      window_names[w]
    )
    
    out_list[[window_names[w]]] <-
      summarise_contrast(diff_post) %>%
      mutate(
        outcome = outcome,
        contrast = "window_joint_mixture_75_vs_25",
        window = window_names[w],
        window_label = window_label,
        comparison = paste0(
          "All pollutants shifted jointly from 25th to 75th percentile in ",
          window_label,
          "; other exposure windows fixed at medians"
        )
      )
  }
  
  window_mixture_effects[[outcome]] <- bind_rows(out_list)
}

window_mixture_effects_all <- bind_rows(window_mixture_effects) %>%
  select(
    outcome,
    contrast,
    window,
    window_label,
    mean,
    sd,
    lower,
    upper,
    pr_gt_0,
    significant,
    comparison
  )

# ------------------------------------------------------------
# B. Overall cumulative joint mixture effects
# ------------------------------------------------------------

overall_mixture_effects <- list()

for (outcome in names(fits)) {
  
  cat("\nOverall joint mixture effect:", outcome, "\n")
  
  fit <- fits[[outcome]]
  
  x_low_all <- make_profile(exp.dat_base, q = 0.25)
  x_high_all <- make_profile(exp.dat_base, q = 0.75)
  
  h_low_all <- predict_bkmrdlm_h_new(
    fit = fit,
    xnew = x_low_all,
    gaussian_kernel = use_gaussian_kernel,
    polydegree = poly_degree
  )
  
  h_high_all <- predict_bkmrdlm_h_new(
    fit = fit,
    xnew = x_high_all,
    gaussian_kernel = use_gaussian_kernel,
    polydegree = poly_degree
  )
  
  diff_post <- as.numeric(h_high_all - h_low_all)
  
  overall_mixture_effects[[outcome]] <-
    summarise_contrast(diff_post) %>%
    mutate(
      outcome = outcome,
      contrast = "overall_joint_mixture_75_vs_25_all_windows",
      comparison = "All pollutants shifted jointly from 25th to 75th percentile across all first-1000-days exposure windows"
    )
}

overall_mixture_effects_all <- bind_rows(overall_mixture_effects) %>%
  select(
    outcome,
    contrast,
    mean,
    sd,
    lower,
    upper,
    pr_gt_0,
    significant,
    comparison
  )

write.xlsx(
  list(
    window_joint_mixture_effects = window_mixture_effects_all,
    overall_joint_mixture_effects = overall_mixture_effects_all
  ),
  file = paste0(
    "bkmrdlm_first1000days_joint_mixture_effects_poly",
    poly_degree,
    ".xlsx"
  ),
  overwrite = TRUE
)

window_mixture_effects_all
overall_mixture_effects_all

# ==============================================================================
# COMBINED POSTERIOR PLOTS
# ==============================================================================
# The fitted objects and weight estimates created above are used directly.

# Ensure predictUnivariate() sees the evaluated kernel settings in the stored call.
fits <- lapply(fits, function(fit) {
  fit$call$gaussian <- FALSE
  fit$call$polydegree <- poly_degree
  fit
})

weights_all <- weights_all %>%
  mutate(
    outcome = factor(outcome, levels = outcomes),
    window = factor(
      as.character(window),
      levels = c(
        "Trimester 1",
        "Trimester 2",
        "Trimester 3",
        "0–6 months",
        "6–12 months",
        "12–18 months",
        "18–24 months"
      )
    )
  )

# ============================================================
# 1. COMBINED WEIGHT PLOT
# ============================================================

combined_plot_dat <- weights_all %>%
  mutate(
    outcome = factor(outcome, levels = outcomes)
  )

p_all <- ggplot(
  combined_plot_dat,
  aes(x = window, y = mean, group = 1)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.20
  ) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    outcome ~ pollutant,
    scales = "free_y"
  ) +
  labs(
    title = paste0(
      "Lag-window weights BKMR-DLM, polynomial kernel degree ",
      poly_degree
    ),
    x = "Developmental exposure window",
    y = "Estimated lag weight"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0)
  )

print(p_all)

ggsave(
  filename = paste0(
    "bkmrdlm_first1000days_looseWindows_weights_poly",
    poly_degree,
    "_ALL_OUTCOMES.png"
  ),
  plot = p_all,
  width = 12,
  height = 20,
  dpi = 600
)

# ============================================================
# 2. DOSE-RESPONSE PLOTS LIKE REGIMES VIGNETTE
# ============================================================

dose_response_all <- list()

for (outcome in names(fits)) {
  
  cat("\nPredicting dose-response:", outcome, "\n")
  
  univar_pred <- predictUnivariate(
    fits[[outcome]],
    points = 30
  )
  
  univar_pred$outcome <- outcome
  
  dose_response_all[[outcome]] <- univar_pred
}

dose_response_all <- bind_rows(dose_response_all) %>%
  mutate(
    outcome = factor(outcome, levels = outcomes)
  )

write.xlsx(
  list(
    dose_response = dose_response_all
  ),
  file = paste0(
    "bkmrdlm_first1000days_dose_response_poly",
    poly_degree,
    "_ALL_OUTCOMES.xlsx"
  ),
  overwrite = TRUE
)

p_dose_all <- ggplot(
  dose_response_all,
  aes(x = E, y = mean, ymin = lower, ymax = upper)
) +
  geom_ribbon(alpha = 0.30) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    outcome ~ reorder(name, m),
    scales = "free"
  ) +
  labs(
    title = paste0(
      "Dose-response estimates BKMR-DLM, polynomial kernel degree ",
      poly_degree
    ),
    x = "Exposure level",
    y = "Estimated h(x)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_dose_all)

ggsave(
  filename = paste0(
    "bkmrdlm_first1000days_dose_response_poly",
    poly_degree,
    "_ALL_OUTCOMES.png"
  ),
  plot = p_dose_all,
  width = 12,
  height = 20,
  dpi = 600
)

# ============================================================
# 3. INTERACTION PLOTS LIKE REGIMES VIGNETTE
# ============================================================

interaction_all <- list()

for (outcome in names(fits)) {
  
  cat("\nPredicting interactions:", outcome, "\n")
  
  fit <- fits[[outcome]]
  
  combodat <- NULL
  
  for (q in c(0.10, 0.25, 0.75, 0.90)) {
    
    for (m in seq_along(fit$x)) {
      
      tmp <- predictUnivariate(
        fit,
        points = 30,
        crossM = m,
        qtl = q
      )
      
      tmp$outcome <- outcome
      tmp$qtl_value <- q
      
      combodat <- rbind(combodat, tmp)
    }
  }
  
  # Remove self-interactions
  combodat <- combodat[
    combodat$m != combodat$cross_m,
  ]
  
  interaction_all[[outcome]] <- combodat
}

interaction_all <- bind_rows(interaction_all) %>%
  mutate(
    outcome = factor(outcome, levels = outcomes),
    qtl_label = factor(
      qtl_value,
      levels = c(0.10, 0.25, 0.75, 0.90),
      labels = c("10th", "25th", "75th", "90th")
    )
  )

write.xlsx(
  list(
    interactions = interaction_all
  ),
  file = paste0(
    "bkmrdlm_first1000days_interactions_poly",
    poly_degree,
    "_ALL_OUTCOMES.xlsx"
  ),
  overwrite = TRUE
)

p_interaction_all <- ggplot(
  interaction_all,
  aes(
    x = E,
    y = mean,
    color = qtl_label,
    group = qtl_label
  )
) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    outcome + reorder(cross, cross_m) ~ reorder(name, m),
    scales = "free"
  ) +
  labs(
    title = paste0(
      "Interaction estimates BKMR-DLM, polynomial kernel degree ",
      poly_degree
    ),
    x = "Exposure level",
    y = "Estimated h(x)",
    color = "Other exposure\nquantile"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p_interaction_all)

ggsave(
  filename = paste0(
    "bkmrdlm_first1000days_interactions_poly",
    poly_degree,
    "_ALL_OUTCOMES.png"
  ),
  plot = p_interaction_all,
  width = 12,
  height = 32,
  dpi = 600
)

cat("\nDone. Polynomial kernel plots exported.\n")

# ==============================================================================
# SOFTWARE INFORMATION
# ==============================================================================
capture.output(
  sessionInfo(),
  file = "sessionInfo.txt"
)

cat("\nSession information written to sessionInfo.txt\n")
