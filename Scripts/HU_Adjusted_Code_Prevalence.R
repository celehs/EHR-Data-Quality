#!/usr/bin/env Rscript

# ================================================================
# Synthetic EHR code prevalence by institution with HU adjustment
# ================================================================
#
# What this script does
# ---------------------
# 1. Reads one or more long-format synthetic EHR files from a folder.
#    Expected columns: patient_num, feature_id, start_date
# 2. Treats each input file as one "institution" (site).
# 3. Calculates healthcare utilization (HU) per patient within site
#    as the number of rows / encounters / observations in that file.
# 4. Estimates code prevalence in two ways:
#      - Unadjusted: crude prevalence by site
#      - HU-adjusted: pooled logistic regression for each code using
#        log(1 + HU), then standardization over a site-specific
#        synthetic HU distribution generated from a negative binomial.
# 5. Writes result tables and saves clean plots.
#
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

# -----------------------------
# User settings
# -----------------------------
input_dir <- "../Data"
output_dir <- file.path(input_dir, "HU_adjusted_results")

# Files to include. By default, all four synthetic files are used.
input_files <- c(
  "fake_codified.csv",
  "fake_nlp.csv",
  "fake_codified_v2.csv",
  "fake_nlp_v2.csv"
)

# Optional: restrict to a specific calendar year.
# Set to NULL to use all dates.
filter_year <- NULL

# Number of synthetic HU draws per site for standardization.
# Larger = smoother estimates, slightly slower runtime.
n_sim_hu <- 5000

# Minimum number of patients with a code across all sites required
# before fitting the adjusted model.
min_positive_patients <- 5

# Number of top codes to show in the summary bar plots.
top_n_codes <- 20

# Prefix filter for plotting. Set to NULL to keep all codes.
# Example: c("PheCode", "LOINC", "RXNORM")
plot_prefix_filter <- NULL

# -----------------------------
# Helper functions
# -----------------------------

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

clean_site_name <- function(filename) {
  tools::file_path_sans_ext(basename(filename))
}

# Estimate negative binomial parameters from utilization counts.
# We use a method-of-moments estimate with a stable fallback when the
# observed variance is not greater than the mean.
estimate_nb_params <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) {
    stop("No valid utilization values available.")
  }
  
  mu <- mean(x)
  v <- stats::var(x)
  
  if (is.na(v) || v <= mu || mu == 0) {
    # Near-Poisson fallback: use a large size to approximate Poisson.
    size <- 1e6
  } else {
    size <- mu^2 / (v - mu)
  }
  
  list(mu = mu, size = size)
}

simulate_site_hu <- function(x, n_sim = 5000) {
  pars <- estimate_nb_params(x)
  sim <- rnbinom(n = n_sim, mu = pars$mu, size = pars$size)
  list(sim_hu = sim, mu = pars$mu, size = pars$size)
}

# Fit pooled logistic regression for one code and compute site-specific
# standardized prevalence using the site's own synthetic HU distribution.
fit_one_code <- function(y, log_hu, site, sim_hu_by_site, code_id) {
  out_unadj <- tapply(y, site, mean)
  out_n <- tapply(y, site, length)
  out_cases <- tapply(y, site, sum)
  
  # Skip unstable fits for very rare or invariant outcomes.
  if (sum(y) < min_positive_patients || length(unique(y)) < 2) {
    return(data.frame(
      site = names(out_unadj),
      feature_id = code_id,
      prevalence_unadjusted = as.numeric(out_unadj),
      prevalence_adjusted = NA_real_,
      n_patients = as.integer(out_n),
      n_positive = as.integer(out_cases),
      model_converged = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  
  fit <- tryCatch(
    suppressWarnings(glm(y ~ log_hu, family = binomial())),
    error = function(e) NULL
  )
  
  if (is.null(fit) || any(!is.finite(stats::coef(fit)))) {
    return(data.frame(
      site = names(out_unadj),
      feature_id = code_id,
      prevalence_unadjusted = as.numeric(out_unadj),
      prevalence_adjusted = NA_real_,
      n_patients = as.integer(out_n),
      n_positive = as.integer(out_cases),
      model_converged = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  
  beta <- stats::coef(fit)
  adjusted <- vapply(names(sim_hu_by_site), function(s) {
    hu_sim <- sim_hu_by_site[[s]]
    mean(stats::plogis(beta[1] + beta[2] * log1p(hu_sim)))
  }, numeric(1))
  
  data.frame(
    site = names(out_unadj),
    feature_id = code_id,
    prevalence_unadjusted = as.numeric(out_unadj),
    prevalence_adjusted = as.numeric(adjusted[names(out_unadj)]),
    n_patients = as.integer(out_n),
    n_positive = as.integer(out_cases),
    model_converged = TRUE,
    intercept = unname(beta[1]),
    beta_log_hu = unname(beta[2]),
    stringsAsFactors = FALSE
  )
}

make_top_code_plot <- function(df, value_col, title, outfile, top_n = 20, prefix_filter = NULL) {
  plot_df <- df
  
  if (!is.null(prefix_filter)) {
    pattern <- paste0("^(", paste(prefix_filter, collapse = "|"), ")")
    plot_df <- plot_df %>% filter(grepl(pattern, feature_id))
  }
  
  keep_codes <- plot_df %>%
    group_by(feature_id) %>%
    summarize(all_sites_present = all(!is.na(.data[[value_col]])),
              max_value = max(.data[[value_col]], na.rm = TRUE),
              .groups = "drop") %>%
    filter(all_sites_present) %>%
    slice_max(order_by = max_value, n = top_n, with_ties = FALSE) %>%
    pull(feature_id)
  
  plot_df <- plot_df %>%
    filter(feature_id %in% keep_codes) %>%
    mutate(feature_id = reorder(feature_id, .data[[value_col]]))
  
  p <- ggplot(plot_df, aes(x = feature_id, y = .data[[value_col]], fill = site)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75) +
    coord_flip() +
    labs(
      title = title,
      x = "Feature ID",
      y = "Prevalence",
      fill = "Site"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9)
    )
  
  ggsave(outfile, p, width = 11, height = 8.5, dpi = 300)
  invisible(p)
}

make_hu_distribution_plot <- function(patient_hu, outfile) {
  p <- ggplot(patient_hu, aes(x = hu)) +
    geom_histogram(bins = 30, alpha = 0.7) +
    facet_wrap(~ site, scales = "free_y") +
    labs(
      title = "Healthcare utilization (HU) distribution by site",
      x = "HU = number of records per patient",
      y = "Number of patients"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(outfile, p, width = 10, height = 7, dpi = 300)
  invisible(p)
}

make_adjustment_comparison_plot <- function(df, outfile) {
  delta_df <- df %>%
    mutate(delta_adjustment = prevalence_adjusted - prevalence_unadjusted) %>%
    filter(!is.na(delta_adjustment)) %>%
    group_by(feature_id) %>%
    summarize(mean_abs_delta = mean(abs(delta_adjustment), na.rm = TRUE), .groups = "drop") %>%
    { dplyr::slice_max(., order_by = mean_abs_delta, n = min(top_n_codes, nrow(.)), with_ties = FALSE) } %>%
    inner_join(df, by = "feature_id") %>%
    mutate(delta_adjustment = prevalence_adjusted - prevalence_unadjusted,
           feature_id = reorder(feature_id, delta_adjustment))
  
  p <- ggplot(delta_df, aes(x = feature_id, y = delta_adjustment, fill = site)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75) +
    coord_flip() +
    labs(
      title = "Change after HU adjustment (adjusted - unadjusted)",
      x = "Feature ID",
      y = "Change in prevalence",
      fill = "Site"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9)
    )
  
  ggsave(outfile, p, width = 11, height = 8.5, dpi = 300)
  invisible(p)
}

# -----------------------------
# 1. Read and stack all files
# -----------------------------
safe_dir_create(output_dir)

full_paths <- file.path(input_dir, input_files)
missing_files <- full_paths[!file.exists(full_paths)]
if (length(missing_files) > 0) {
  stop(
    "These input files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Reading input files...")
long_list <- lapply(full_paths, function(fp) {
  dt <- fread(fp)
  
  required_cols <- c("patient_num", "feature_id", "start_date")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("File ", basename(fp), " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dt <- dt[, ..required_cols]
  dt[, start_date := as.Date(start_date)]
  dt[, site := clean_site_name(fp)]
  
  if (!is.null(filter_year)) {
    dt <- dt[format(start_date, "%Y") == as.character(filter_year)]
  }
  
  dt
})

long_data <- rbindlist(long_list, use.names = TRUE, fill = TRUE)

if (nrow(long_data) == 0) {
  stop("No rows available after reading files and applying filters.")
}

fwrite(long_data, file.path(output_dir, "stacked_long_data.csv"))

# ------------------------------------------------------
# 2. Build patient-level utilization and code indicators
# ------------------------------------------------------
message("Building patient-level utilization table...")

patient_hu <- long_data %>%
  count(site, patient_num, name = "hu") %>%
  mutate(log_hu = log1p(hu))

fwrite(as.data.table(patient_hu), file.path(output_dir, "patient_level_hu.csv"))

# Distinct patient-code pairs: one row means that patient had the code.
patient_code <- long_data %>%
  distinct(site, patient_num, feature_id)

# Create a stable patient key for matrix construction.
patient_index <- patient_hu %>%
  arrange(site, patient_num) %>%
  mutate(patient_id = row_number())

code_index <- patient_code %>%
  distinct(feature_id) %>%
  arrange(feature_id) %>%
  mutate(code_id = row_number())

patient_code_indexed <- patient_code %>%
  inner_join(patient_index, by = c("site", "patient_num")) %>%
  inner_join(code_index, by = "feature_id")

message("Constructing sparse patient-by-code matrix...")
X <- sparseMatrix(
  i = patient_code_indexed$patient_id,
  j = patient_code_indexed$code_id,
  x = 1L,
  dims = c(nrow(patient_index), nrow(code_index)),
  dimnames = list(NULL, code_index$feature_id)
)

# --------------------------------------------------
# 3. Fit site-specific HU distributions (NB by site)
# --------------------------------------------------
message("Fitting site-specific HU distributions...")

site_hu_info <- patient_hu %>%
  group_by(site) %>%
  summarize(
    mean_hu = mean(hu),
    var_hu = var(hu),
    n_patients = n(),
    .groups = "drop"
  )

sim_hu_by_site <- lapply(split(patient_hu$hu, patient_hu$site), simulate_site_hu, n_sim = n_sim_hu)
sim_hu_draws <- lapply(names(sim_hu_by_site), function(s) {
  data.frame(site = s, sim_hu = sim_hu_by_site[[s]]$sim_hu)
}) %>% bind_rows()

nb_params <- data.frame(
  site = names(sim_hu_by_site),
  mu = sapply(sim_hu_by_site, function(x) x$mu),
  size = sapply(sim_hu_by_site, function(x) x$size),
  stringsAsFactors = FALSE
)

fwrite(as.data.table(site_hu_info), file.path(output_dir, "site_hu_summary.csv"))
fwrite(as.data.table(nb_params), file.path(output_dir, "site_hu_negative_binomial_parameters.csv"))
fwrite(as.data.table(sim_hu_draws), file.path(output_dir, "site_simulated_hu_draws.csv"))

# -------------------------------------------------------
# 4. Estimate code prevalence with and without adjustment
# -------------------------------------------------------
message("Estimating code prevalence for each feature...")

site_vec <- patient_index$site
log_hu_vec <- patient_index$log_hu

results_list <- vector("list", length = ncol(X))
for (j in seq_len(ncol(X))) {
  y <- as.integer(X[, j] > 0)
  results_list[[j]] <- fit_one_code(
    y = y,
    log_hu = log_hu_vec,
    site = site_vec,
    sim_hu_by_site = lapply(sim_hu_by_site, `[[`, "sim_hu"),
    code_id = colnames(X)[j]
  )
  
  if (j %% 100 == 0 || j == ncol(X)) {
    message(sprintf("  Processed %s / %s codes", j, ncol(X)))
  }
}

results <- bind_rows(results_list) %>%
  arrange(feature_id, site)

results <- results %>%
  mutate(
    adjustment_delta = prevalence_adjusted - prevalence_unadjusted,
    abs_adjustment_delta = abs(adjustment_delta)
  )

fwrite(as.data.table(results), file.path(output_dir, "code_prevalence_by_site_adjusted_and_unadjusted.csv"))

# Summary table at the code level.
code_summary <- results %>%
  group_by(feature_id) %>%
  summarize(
    mean_unadjusted = mean(prevalence_unadjusted, na.rm = TRUE),
    mean_adjusted = mean(prevalence_adjusted, na.rm = TRUE),
    max_adjusted = max(prevalence_adjusted, na.rm = TRUE),
    mean_abs_adjustment = mean(abs_adjustment_delta, na.rm = TRUE),
    total_positive = sum(n_positive, na.rm = TRUE),
    all_models_converged = all(model_converged),
    .groups = "drop"
  ) %>%
  arrange(desc(max_adjusted))

fwrite(as.data.table(code_summary), file.path(output_dir, "code_summary.csv"))

# -----------------------------
# 5. Make plots
# -----------------------------
message("Making plots...")

make_hu_distribution_plot(
  patient_hu = patient_hu,
  outfile = file.path(output_dir, "plot_hu_distribution_by_site.png")
)

make_top_code_plot(
  df = results,
  value_col = "prevalence_unadjusted",
  title = sprintf("Top %s features by unadjusted prevalence", top_n_codes),
  outfile = file.path(output_dir, "plot_top_codes_unadjusted.png"),
  top_n = top_n_codes,
  prefix_filter = plot_prefix_filter
)

make_top_code_plot(
  df = results,
  value_col = "prevalence_adjusted",
  title = sprintf("Top %s features by HU-adjusted prevalence", top_n_codes),
  outfile = file.path(output_dir, "plot_top_codes_hu_adjusted.png"),
  top_n = top_n_codes,
  prefix_filter = plot_prefix_filter
)

make_adjustment_comparison_plot(
  df = results,
  outfile = file.path(output_dir, "plot_adjustment_delta.png")
)

# Optional: save session info for reproducibility.
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("Done. Outputs written to: ", output_dir)

# -----------------------------
# 6. Brief console summary
# -----------------------------
message("\nTop features by adjusted prevalence:")
print(utils::head(code_summary, 10))
