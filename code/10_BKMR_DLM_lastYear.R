# ==============================================================================
# BKMR-DLM analysis: monthly air-pollution exposures during the year before FU2
# ==============================================================================
# Purpose
# -------
# R code used for the BKMR-DLM analyses reported in the manuscript.
#
# Exposure time scale
# -------------------
# X_day0  = day of FU2 / closest available day
# X_day1  = 1 day before FU2
# month0  = days 0-29 before FU2
# month11 = days 330-359 before FU2
#
# Main model specification
# ------------------------
# Lag basis: natural spline ("ns"), df = 3
# Kernel: polynomial
# Polynomial degree: 2
# MCMC iterations: 10,000
# Burn-in: 5,000
#
# Random-number generation
# ------------------------
# A fixed RNG seed can optionally be specified below for reproducible reruns.
#
# Data requirements
# -----------------
# This script assumes that the following R objects have already been created by
# the study's data-preparation workflow:
#   - PM25_fu2_withMeanExp
#   - NO2_fu2_withMeanExp
#   - BC_fu2_withMeanExp
#   - Final_dataset_eye_tracking_study
#
# Individual-level study data are not included in the public repository because
# of privacy and ethical restrictions.
#
# Outputs are written to the current working directory.
# ==============================================================================

library(regimes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)

# ------------------------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------------------------

poly_degree <- 2

# Optional RNG seed for reproducible reruns
rng_seed <- NULL   # e.g., 12345

if (!is.null(rng_seed)) {
  set.seed(rng_seed)
}

required_objects <- c(
  "PM25_fu2_withMeanExp",
  "NO2_fu2_withMeanExp",
  "BC_fu2_withMeanExp",
  "Final_dataset_eye_tracking_study"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Required input object(s) not found: ",
    paste(missing_objects, collapse = ", "),
    ". See README.md for required inputs."
  )
}

# ------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------

poly_degree <- 2
outcomes <- c(
  "ERc_FixN",
  "ERc_SacN",
  "ERc_PupD",
  "ERc_SacV",
  "FixTD_biasHappy",
  "FixTD_biasAnger",
  "FixTD_biasFear"
)

num_covs <- c("AgeM", "AgeChild", "sin_season", "cos_season")
fac_covs <- c("DiplM", "Smoke", "Parity", "Seks", "Ethn")

# ------------------------------------------------------------
# Function: daily wide data -> monthly exposure matrix
# ------------------------------------------------------------

make_first360_monthly <- function(df, id_col = "ID") {
  
  day_cols <- grep("^X_day[0-9]+$", names(df), value = TRUE)
  day_nums <- as.integer(sub("X_day", "", day_cols))
  
  day_order <- order(day_nums)
  day_cols <- day_cols[day_order]
  day_nums <- day_nums[day_order]
  
  needed_days <- 0:359
  
  if (!all(needed_days %in% day_nums)) {
    stop("Not all X_day0 to X_day359 columns are present.")
  }
  
  day_cols_360 <- day_cols[match(needed_days, day_nums)]
  daily_mat <- as.matrix(df[, day_cols_360])
  
  monthly_mat <- t(apply(daily_mat, 1, function(x) {
    
    sapply(0:11, function(m) {
      
      idx <- (m * 30 + 1):((m + 1) * 30)
      vals <- x[idx]
      
      if (all(is.na(vals))) {
        NA_real_
      } else {
        mean(vals, na.rm = TRUE)
      }
    })
  }))
  
  colnames(monthly_mat) <- paste0("month", 0:11)
  rownames(monthly_mat) <- as.character(df[[id_col]])
  
  monthly_mat
}

# ------------------------------------------------------------
# Clean exposure datasets
# ------------------------------------------------------------

PM25_fu2_clean <- PM25_fu2_withMeanExp[, -c(2:5)]
NO2_fu2_clean  <- NO2_fu2_withMeanExp[, -c(2:5)]
BC_fu2_clean   <- BC_fu2_withMeanExp[, -c(2:5)]

# ------------------------------------------------------------
# Create monthly exposure matrices
# ------------------------------------------------------------

PM25_fu2_monthly <- make_first360_monthly(PM25_fu2_clean, id_col = "ID")
NO2_fu2_monthly  <- make_first360_monthly(NO2_fu2_clean,  id_col = "ID")
BC_fu2_monthly   <- make_first360_monthly(BC_fu2_clean,   id_col = "ID")

cat("\nExposure dimensions:\n")
print(dim(PM25_fu2_monthly))
print(dim(NO2_fu2_monthly))
print(dim(BC_fu2_monthly))

cat("\nMissing exposure values:\n")
print(sum(is.na(PM25_fu2_monthly)))
print(sum(is.na(NO2_fu2_monthly)))
print(sum(is.na(BC_fu2_monthly)))

# ------------------------------------------------------------
# Restrict to eye-tracking study IDs and align order
# ------------------------------------------------------------

Final_dataset_eye_tracking_study$ID <- as.character(Final_dataset_eye_tracking_study$ID)
ID_vector <- Final_dataset_eye_tracking_study$ID

PM25_sub <- PM25_fu2_monthly[
  match(ID_vector, rownames(PM25_fu2_monthly)),
]

NO2_sub <- NO2_fu2_monthly[
  match(ID_vector, rownames(NO2_fu2_monthly)),
]

BC_sub <- BC_fu2_monthly[
  match(ID_vector, rownames(BC_fu2_monthly)),
]

rownames(PM25_sub) <- ID_vector
rownames(NO2_sub)  <- ID_vector
rownames(BC_sub)   <- ID_vector

stopifnot(all(rownames(PM25_sub) == ID_vector))
stopifnot(all(rownames(NO2_sub) == ID_vector))
stopifnot(all(rownames(BC_sub) == ID_vector))

# ------------------------------------------------------------
# Joint complete-case exposure exclusion
# ------------------------------------------------------------

complete_exp_idx <-
  complete.cases(PM25_sub) &
  complete.cases(NO2_sub) &
  complete.cases(BC_sub)

cat("\nExcluded due to exposure missingness:", sum(!complete_exp_idx), "\n")

PM25_final <- PM25_sub[complete_exp_idx, , drop = FALSE]
NO2_final  <- NO2_sub[complete_exp_idx, , drop = FALSE]
BC_final   <- BC_sub[complete_exp_idx, , drop = FALSE]

exp.dat <- list(
  pm25 = PM25_final,
  no2  = NO2_final,
  bc   = BC_final
)

cat("\nFinal exposure dimensions:\n")
print(lapply(exp.dat, dim))

# ------------------------------------------------------------
# Retrieve FU2 date and create seasonal terms
# ------------------------------------------------------------

FU2_dates <- PM25_fu2_withMeanExp[, c("ID", "date_FU2")]
FU2_dates$ID <- as.character(FU2_dates$ID)

final_IDs <- rownames(exp.dat$pm25)

FU2_dates_sub <- FU2_dates[
  match(final_IDs, FU2_dates$ID),
]

stopifnot(all(final_IDs == FU2_dates_sub$ID))

FU2_dates_sub$date_FU2 <- as.Date(FU2_dates_sub$date_FU2)
FU2_dates_sub$doy <- as.numeric(format(FU2_dates_sub$date_FU2, "%j"))

FU2_dates_sub$sin_season <- sin(2 * pi * FU2_dates_sub$doy / 365.25)
FU2_dates_sub$cos_season <- cos(2 * pi * FU2_dates_sub$doy / 365.25)

# ------------------------------------------------------------
# Create analytic dataset
# ------------------------------------------------------------

dat_model <- Final_dataset_eye_tracking_study[
  match(final_IDs, Final_dataset_eye_tracking_study$ID),
]

dat_model$sin_season <- FU2_dates_sub$sin_season
dat_model$cos_season <- FU2_dates_sub$cos_season

dat_model[num_covs] <- lapply(dat_model[num_covs], as.numeric)
dat_model[fac_covs] <- lapply(dat_model[fac_covs], as.factor)

stopifnot(all(dat_model$ID == final_IDs))

# ============================================================
# FIT BKMR-DLM MODELS
# Polynomial kernel, degree 2
# ============================================================

fits <- list()
summaries <- list()
model_N <- data.frame()

for (outcome in outcomes) {
  
  cat("\nRunning outcome:", outcome, "\n")
  
  out.var <- as.numeric(dat_model[[outcome]])
  
  keep <- complete.cases(
    out.var,
    dat_model[, c(num_covs, fac_covs)]
  )
  
  exp.dat_out <- lapply(exp.dat, function(x) {
    x[keep, , drop = FALSE]
  })
  
  dat_out <- dat_model[keep, ]
  out.var_out <- out.var[keep]
  
  cov.dat_out <- model.matrix(
    as.formula(
      paste("~", paste(c(num_covs, fac_covs), collapse = " + "))
    ),
    data = dat_out
  )[, -1, drop = FALSE]
  
  stopifnot(length(out.var_out) == nrow(cov.dat_out))
  stopifnot(all(sapply(exp.dat_out, nrow) == length(out.var_out)))
  stopifnot(sum(is.na(out.var_out)) == 0)
  stopifnot(sum(is.na(cov.dat_out)) == 0)
  stopifnot(all(sapply(exp.dat_out, function(x) sum(is.na(x))) == 0))
  
  cat("N used:", length(out.var_out), "\n")
  
  model_N <- rbind(
    model_N,
    data.frame(
      outcome = outcome,
      N_used = length(out.var_out)
    )
  )
  
  fit <- bkmrdlm(
    y = out.var_out,
    x = exp.dat_out,
    z = cov.dat_out,
    niter = 10000,
    nburn = 5000,
    basis.opts = list(type = "ns", df = 3),
    gaussian = FALSE,
    polydegree = poly_degree
  )
  
  fit$call$gaussian <- FALSE
  fit$call$polydegree <- poly_degree
  
  fits[[outcome]] <- fit
  summaries[[outcome]] <- summary(fit)
  
  save(
    fit,
    file = paste0(
      "bkmrdlm_fu2_first360_monthly_df3_poly",
      poly_degree,
      "_seasonAdj_",
      outcome,
      ".RData"
    )
  )
}

print(model_N)

# ------------------------------------------------------------
# Save final monthly exposure matrices
# ------------------------------------------------------------

write.xlsx(
  list(
    PM25_monthly = as.data.frame(PM25_final),
    NO2_monthly  = as.data.frame(NO2_final),
    BC_monthly   = as.data.frame(BC_final),
    analytic_dataset = dat_model,
    model_N = model_N
  ),
  file = paste0(
    "fu2_first360_monthly_exposures_complete_poly",
    poly_degree,
    ".xlsx"
  ),
  rowNames = TRUE,
  overwrite = TRUE
)

# ============================================================
# EXTRACT MONTHLY WEIGHTS
# ============================================================

extract_weights <- function(fit, outcome_name, exp.dat) {
  
  weights <- as.data.frame(summary(fit)$weights)
  
  weights$pollutant <- names(exp.dat)[weights$m]
  weights$month <- weights$t - 1
  weights$outcome <- outcome_name
  
  weights %>%
    select(
      outcome,
      pollutant,
      m,
      month,
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
      outcome_name = outcome,
      exp.dat = exp.dat
    )
  })
)

weights_all <- weights_all %>%
  mutate(
    excludes_zero = lower > 0 | upper < 0,
    outcome = factor(outcome, levels = outcomes)
  )

weights_summary <- weights_all %>%
  group_by(outcome, pollutant) %>%
  summarise(
    mean_abs_weight = mean(abs(mean), na.rm = TRUE),
    max_abs_weight  = max(abs(mean), na.rm = TRUE),
    month_of_max_abs_weight = month[which.max(abs(mean))],
    .groups = "drop"
  ) %>%
  arrange(outcome, desc(mean_abs_weight))

critical_windows <- weights_all %>%
  filter(excludes_zero) %>%
  arrange(outcome, pollutant, month)

write.xlsx(
  list(
    model_N = model_N,
    weights_all = weights_all,
    weights_summary = weights_summary,
    critical_windows = critical_windows
  ),
  file = paste0(
    "bkmrdlm_fu2_first360_monthly_df3_weights_poly",
    poly_degree,
    "_all_outcomes.xlsx"
  ),
  overwrite = TRUE
)
# ============================================================
# INDIVIDUAL WEIGHT PLOTS
# ============================================================

for (outcome in names(fits)) {
  
  plot_dat <- weights_all %>%
    filter(outcome == !!outcome)
  
  p <- ggplot(plot_dat, aes(x = month, y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ pollutant, scales = "free_y") +
    scale_x_continuous(breaks = 0:11, limits = c(0, 11)) +
    labs(
      x = "Months before FU2",
      y = "Estimated lag weight",
      title = paste0("Lag weight functions: ", outcome),
      subtitle = paste0(
        "Polynomial kernel degree ",
        poly_degree,
        "; month0 = days 0–29 before FU2"
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = paste0(
      "bkmrdlm_fu2_first360_monthly_df3_weights_poly",
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
# PANEL FIGURE: ALL WEIGHT PLOTS
# ============================================================

p_weights_all <- ggplot(
  weights_all,
  aes(x = month, y = mean)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    outcome ~ pollutant,
    scales = "free_y"
  ) +
  scale_x_continuous(breaks = 0:11, limits = c(0, 11)) +
  labs(
    title = paste0(
      "Lag-window weights BKMR-DLM, FU2 monthly exposure, polynomial kernel degree ",
      poly_degree
    ),
    x = "Months before FU2",
    y = "Estimated lag weight"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_weights_all)

ggsave(
  filename = paste0(
    "bkmrdlm_fu2_first360_monthly_weights_poly",
    poly_degree,
    "_ALL_OUTCOMES.png"
  ),
  plot = p_weights_all,
  width = 12,
  height = 20,
  dpi = 600
)

# ============================================================
# DOSE-RESPONSE PLOTS LIKE REGIMES VIGNETTE
# ============================================================

fits <- lapply(fits, function(fit) {
  fit$call$gaussian <- FALSE
  fit$call$polydegree <- poly_degree
  fit
})

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
    "bkmrdlm_fu2_first360_monthly_dose_response_poly",
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
      "Dose-response estimates BKMR-DLM, FU2 monthly exposure, polynomial kernel degree ",
      poly_degree
    ),
    x = "Standardized exposure level",
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
    "bkmrdlm_fu2_first360_monthly_dose_response_poly",
    poly_degree,
    "_ALL_OUTCOMES.png"
  ),
  plot = p_dose_all,
  width = 12,
  height = 20,
  dpi = 600
)

# ============================================================
# INTERACTION PLOTS LIKE REGIMES VIGNETTE
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
    "bkmrdlm_fu2_first360_monthly_interactions_poly",
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
      "Interaction estimates BKMR-DLM, FU2 monthly exposure, polynomial kernel degree ",
      poly_degree
    ),
    x = "Standardized exposure level",
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
    "bkmrdlm_fu2_first360_monthly_interactions_poly",
    poly_degree,
    "_ALL_OUTCOMES.png"
  ),
  plot = p_interaction_all,
  width = 12,
  height = 32,
  dpi = 600
)

# ============================================================
# PREDICTION CONTRASTS:
# 1. Single-pollutant exposure-response contrasts
# 2. Monthly joint mixture effects
# 3. Overall joint mixture effect
# ============================================================

get_univar_contrast <- function(fit, outcome_name, points = 30) {
  
  fit$call$gaussian <- FALSE
  fit$call$polydegree <- poly_degree
  
  univar <- predictUnivariate(fit, points = points)
  
  contrast <- univar %>%
    group_by(name) %>%
    summarise(
      idx_low  = which.min(abs(E - quantile(E, 0.25, na.rm = TRUE))),
      idx_high = which.min(abs(E - quantile(E, 0.75, na.rm = TRUE))),
      E_low  = E[idx_low],
      E_high = E[idx_high],
      pred_low  = mean[idx_low],
      pred_high = mean[idx_high],
      lower_low  = lower[idx_low],
      upper_low  = upper[idx_low],
      lower_high = lower[idx_high],
      upper_high = upper[idx_high],
      effect_75_vs_25 = pred_high - pred_low,
      effect_lower = lower_high - upper_low,
      effect_upper = upper_high - lower_low,
      significant = !(effect_lower <= 0 & effect_upper >= 0),
      .groups = "drop"
    ) %>%
    mutate(outcome = outcome_name) %>%
    select(
      outcome,
      pollutant = name,
      E_low,
      E_high,
      pred_low,
      pred_high,
      effect_75_vs_25,
      effect_lower,
      effect_upper,
      significant,
      lower_low,
      upper_low,
      lower_high,
      upper_high
    )
  
  contrast
}

single_pollutant_contrasts_all <- bind_rows(
  lapply(names(fits), function(outcome) {
    get_univar_contrast(
      fit = fits[[outcome]],
      outcome_name = outcome,
      points = 30
    )
  })
)

# Custom posterior prediction helper used for the joint-mixture contrasts below.
# This section is retained unchanged computationally from the accepted analysis.
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

# In the exported tables, `significant` means that the 95% posterior interval
# for the contrast excludes zero; it is not a frequentist p-value.
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

# ------------------------------------------------------------
# Monthly joint mixture effects
# NOTE: percentile profiles are based on exp.dat (the exposure-complete sample),
# while each fitted outcome model may use a smaller outcome-specific complete-case
# subset. This is retained as in the accepted analysis and should match manuscript
# wording about the reference distribution.
# ------------------------------------------------------------

monthly_mixture_effects <- list()

for (outcome in names(fits)) {
  
  cat("\nMonthly joint mixture effects:", outcome, "\n")
  
  fit <- fits[[outcome]]
  x_med <- make_profile(exp.dat, q = 0.50)
  out_month <- list()
  
  for (month in 0:11) {
    
    x_low <- x_med
    x_high <- x_med
    
    for (pollutant in names(exp.dat)) {
      
      x_low[[pollutant]][1, month + 1] <-
        quantile(exp.dat[[pollutant]][, month + 1], 0.25, na.rm = TRUE)
      
      x_high[[pollutant]][1, month + 1] <-
        quantile(exp.dat[[pollutant]][, month + 1], 0.75, na.rm = TRUE)
    }
    
    h_low <- predict_bkmrdlm_h_new(
      fit = fit,
      xnew = x_low,
      gaussian_kernel = FALSE,
      polydegree = poly_degree
    )
    
    h_high <- predict_bkmrdlm_h_new(
      fit = fit,
      xnew = x_high,
      gaussian_kernel = FALSE,
      polydegree = poly_degree
    )
    
    diff_post <- as.numeric(h_high - h_low)
    
    out_month[[paste0("month", month)]] <-
      summarise_contrast(diff_post) %>%
      mutate(
        outcome = outcome,
        contrast = "monthly_joint_mixture_75_vs_25",
        month = month,
        window = paste0("month", month),
        comparison = paste0(
          "All pollutants shifted jointly from 25th to 75th percentile at month ",
          month,
          "; other months fixed at medians"
        )
      )
  }
  
  monthly_mixture_effects[[outcome]] <- bind_rows(out_month)
}

monthly_mixture_effects_all <- bind_rows(monthly_mixture_effects) %>%
  select(
    outcome,
    contrast,
    month,
    window,
    mean,
    sd,
    lower,
    upper,
    pr_gt_0,
    significant,
    comparison
  )

# ------------------------------------------------------------
# Overall joint mixture effect over full year
# ------------------------------------------------------------

overall_mixture_effects <- list()

for (outcome in names(fits)) {
  
  cat("\nOverall joint mixture effect:", outcome, "\n")
  
  fit <- fits[[outcome]]
  
  x_low_all <- make_profile(exp.dat, q = 0.25)
  x_high_all <- make_profile(exp.dat, q = 0.75)
  
  h_low_all <- predict_bkmrdlm_h_new(
    fit = fit,
    xnew = x_low_all,
    gaussian_kernel = FALSE,
    polydegree = poly_degree
  )
  
  h_high_all <- predict_bkmrdlm_h_new(
    fit = fit,
    xnew = x_high_all,
    gaussian_kernel = FALSE,
    polydegree = poly_degree
  )
  
  diff_post <- as.numeric(h_high_all - h_low_all)
  
  overall_mixture_effects[[outcome]] <-
    summarise_contrast(diff_post) %>%
    mutate(
      outcome = outcome,
      contrast = "overall_joint_mixture_75_vs_25_all_months",
      comparison = "All pollutants shifted jointly from 25th to 75th percentile across all 12 months"
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
    single_pollutant_within_mixture = single_pollutant_contrasts_all,
    monthly_joint_mixture_effects = monthly_mixture_effects_all,
    overall_joint_mixture_effects = overall_mixture_effects_all
  ),
  file = paste0(
    "bkmrdlm_fu2_first360_monthly_prediction_effects_poly",
    poly_degree,
    ".xlsx"
  ),
  overwrite = TRUE
)

cat("\nDone. Polynomial kernel FU2 monthly BKMR-DLM script completed.\n")

###################################################################################
# ============================================================
# EXPORT JOINT MIXTURE EFFECTS TO EXCEL
# Standalone export section retained from the analysis workflow
# ============================================================

# ------------------------------------------------------------
# Check objects exist
# ------------------------------------------------------------

stopifnot(exists("monthly_mixture_effects_all"))
stopifnot(exists("overall_mixture_effects_all"))
stopifnot(exists("poly_degree"))

# ------------------------------------------------------------
# Export monthly + overall joint mixture effects
# ------------------------------------------------------------

write.xlsx(
  list(
    monthly_joint_mixture_effects = monthly_mixture_effects_all,
    overall_joint_mixture_effects = overall_mixture_effects_all
  ),
  file = paste0(
    "bkmrdlm_fu2_first360_monthly_joint_mixture_effects_poly",
    poly_degree,
    ".xlsx"
  ),
  overwrite = TRUE
)

cat(
  "\nExported:\n",
  paste0(
    "bkmrdlm_fu2_first360_monthly_joint_mixture_effects_poly",
    poly_degree,
    ".xlsx"
  ),
  "\n"
)

# ============================================================
# OVERALL CUMULATIVE MIXTURE DOSE-RESPONSE CURVES
# Entire 12-month period jointly shifted
# ============================================================

overall_mixture_dose_response <- list()

mix_q_seq <- seq(0.05, 0.95, by = 0.05)

for (outcome in names(fits)) {
  
  cat("\nOverall mixture dose-response:", outcome, "\n")
  
  fit <- fits[[outcome]]
  
  h_list <- list()
  
  for (q in mix_q_seq) {
    
    x_q <- make_profile(exp.dat, q = q)
    
    h_q <- predict_bkmrdlm_h_new(
      fit = fit,
      xnew = x_q,
      gaussian_kernel = FALSE,
      polydegree = poly_degree
    )
    
    h_list[[as.character(q)]] <- data.frame(
      q = q,
      mean = mean(h_q),
      lower = quantile(h_q, 0.025),
      upper = quantile(h_q, 0.975),
      outcome = outcome
    )
  }
  
  overall_mixture_dose_response[[outcome]] <-
    bind_rows(h_list)
}

overall_mixture_dose_response_all <-
  bind_rows(overall_mixture_dose_response)

# ------------------------------------------------------------
# Export Excel
# ------------------------------------------------------------

write.xlsx(
  list(
    overall_mixture_dose_response =
      overall_mixture_dose_response_all
  ),
  file = paste0(
    "bkmrdlm_fu2_first360_monthly_overall_mixture_dose_response_poly",
    poly_degree,
    ".xlsx"
  ),
  overwrite = TRUE
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_mixdose <- ggplot(
  overall_mixture_dose_response_all,
  aes(
    x = q,
    y = mean,
    ymin = lower,
    ymax = upper
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_ribbon(alpha = 0.25) +
  geom_line(linewidth = 0.9) +
  facet_wrap(
    ~ outcome,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = paste0(
      "Overall cumulative mixture dose-response\n",
      "Polynomial kernel degree ",
      poly_degree
    ),
    x = "Joint mixture percentile across all 12 months",
    y = "Estimated h(x)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank()
  )

print(p_mixdose)

ggsave(
  filename = paste0(
    "bkmrdlm_fu2_first360_monthly_overall_mixture_dose_response_poly",
    poly_degree,
    ".png"
  ),
  plot = p_mixdose,
  width = 10,
  height = 12,
  dpi = 600
)

# ==============================================================================
# REPRODUCIBILITY INFORMATION
# ==============================================================================
# Run this after the analysis to record the R and package versions used for the
# public rerun. Note that this records the current environment, not necessarily
# the historical environment of the originally accepted run.
capture.output(
  sessionInfo(),
  file = "sessionInfo.txt"
)

cat("\nSession information written to sessionInfo.txt\n")
