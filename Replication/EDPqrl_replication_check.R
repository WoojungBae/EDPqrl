# ==============================================================================
# PSQC simulation tables
#
# Manuscript table order:
#   1. Main-text Table 1
#   2. Supplementary Table S1
#   3. Supplementary Table S2
#   4. Supplementary Table S3
#
# This script generates only the simulation tables reported in the manuscript.
# ==============================================================================

# ------------------------------------------------------------------------------
# R packages
# ------------------------------------------------------------------------------
library(Rcpp)
library(RcppArmadillo)

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

library(knitr)
library(kableExtra)
library(htmltools)
library(xml2)

# ------------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
source_dir = ""
results_dir = "Results"
truth_dir = ""

source_path = if (nzchar(source_dir)) source_dir else "."
results_path = if (nzchar(results_dir)) results_dir else "."
truth_path = if (nzchar(truth_dir)) truth_dir else "."

# ------------------------------------------------------------------------------
# Source files
# ------------------------------------------------------------------------------
source(file.path(source_path, "EDPqrl_r.R"))
Rcpp::sourceCpp(file.path(source_path, "EDPqrl_cpp.cpp"))

# ------------------------------------------------------------------------------
# General settings
# ------------------------------------------------------------------------------
CE = TRUE
esttype = "mean"

nu_vec = c(0, 1, 2, 3)

# Quantile levels reported in the manuscript:
#   Scenarios 1 and 2: approximately 20% and 40% censoring
#   Scenarios 3 and 4: approximately 60% and 80% censoring
rho_by_scenario = list(
  `1` = c(0.30, 0.60),
  `2` = c(0.30, 0.60),
  `3` = c(0.10, 0.20),
  `4` = c(0.10, 0.20)
)

# The simulation replication script saves:
#   nominal 95% CrI results with "_level5"
#   nominal 99% CrI results with "_level1"
main_level = 0.05
supp_99_level = 0.01

# ==============================================================================
# Optional preprocessing: calculate the common true PSQC
#
# This calculation is not a manuscript table. It creates True_PSQC.txt, which
# is used as the population truth when calculating bias, RMSE, and coverage.
#
# The original calculation is computationally intensive:
#   10,000 Monte Carlo datasets, each with 10,000 observations.
#
# Leave regenerate_true_psqc = FALSE when True_PSQC.txt is already supplied.
# ==============================================================================

regenerate_true_psqc = FALSE
true_psqc_file = file.path(truth_path, "True_PSQC.txt")

if (regenerate_true_psqc) {

  set.seed(1)

  truth_rho_vec = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60)
  truth_nu_vec = c(0, 1, 2, 3)

  n_sim = 1e4
  n_sample = 1e4

  param_grid = expand.grid(
    rho = truth_rho_vec,
    nu = truth_nu_vec
  )

  param_grid = param_grid[, c("nu", "rho")]
  n_params = nrow(param_grid)

  # The true PSQC depends only on Y0 and Y1 and is common across scenarios.
  scn_id = 1

  sum_results = matrix(
    0,
    nrow = n_params,
    ncol = 3
  )

  valid_results = integer(n_params)

  cat("Starting true PSQC calculation...\n")

  for (ii in seq_len(n_sim)) {

    dat = generate_data(
      Scn = scn_id,
      n = n_sample
    )

    for (jj in seq_len(n_params)) {

      curr_nu = param_grid$nu[jj]
      curr_rho = param_grid$rho[jj]

      # Time-specific always-survivor (TAS) principal stratum
      idx_tas = (
        dat$Y0 > curr_nu &
          dat$Y1 > curr_nu
      )

      Y0_tas = dat$Y0[idx_tas] - curr_nu
      Y1_tas = dat$Y1[idx_tas] - curr_nu

      Q0 = if (length(Y0_tas) > 0) {
        as.numeric(
          quantile(
            Y0_tas,
            probs = curr_rho,
            names = FALSE
          )
        )
      } else {
        NA_real_
      }

      Q1 = if (length(Y1_tas) > 0) {
        as.numeric(
          quantile(
            Y1_tas,
            probs = curr_rho,
            names = FALSE
          )
        )
      } else {
        NA_real_
      }

      diff = Q1 - Q0

      if (all(is.finite(c(Q0, Q1, diff)))) {
        sum_results[jj, 1] = sum_results[jj, 1] + Q0
        sum_results[jj, 2] = sum_results[jj, 2] + Q1
        sum_results[jj, 3] = sum_results[jj, 3] + diff

        valid_results[jj] = valid_results[jj] + 1L
      }
    }

    if (ii %% 1000 == 0) {
      cat(sprintf("Completed %d / %d\n", ii, n_sim))
    }
  }

  if (any(valid_results == 0)) {
    stop(
      "No valid true PSQC values were obtained for at least one ",
      "(nu, rho) combination."
    )
  }

  truth_base = cbind(
    param_grid,
    sum_results / valid_results
  )

  colnames(truth_base) = c(
    "nu",
    "rho",
    "True_Res_Y0",
    "True_Res_Y1",
    "True_Diff"
  )

  truth_base = as.data.frame(
    truth_base,
    stringsAsFactors = FALSE
  )

  truth_base = truth_base[
    order(truth_base$nu, truth_base$rho),
    ,
    drop = FALSE
  ]

  rownames(truth_base) = NULL

  write.table(
    truth_base,
    file = true_psqc_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  cat("Saved: ", true_psqc_file, "\n", sep = "")
}

if (!file.exists(true_psqc_file)) {
  stop(
    "Missing true PSQC file: ", true_psqc_file, "\n",
    "Supply True_PSQC.txt or set regenerate_true_psqc = TRUE."
  )
}

true_psqc = read.table(
  true_psqc_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

required_truth_cols = c("nu", "rho", "True_Diff")

if (!all(required_truth_cols %in% names(true_psqc))) {
  stop(
    "True_PSQC.txt must contain columns: ",
    paste(required_truth_cols, collapse = ", ")
  )
}

# ==============================================================================
# Formatting helpers
# ==============================================================================

fmt2 = function(x) {
  x = as.numeric(x)
  x[is.finite(x) & abs(x) < 0.005] = 0
  sprintf("%.2f", x)
}

fmt_cp = function(x) {
  sprintf("%.1f", 100 * as.numeric(x))
}

level_tag = function(level) {
  format(
    100 * level,
    trim = TRUE,
    scientific = FALSE
  )
}

# ==============================================================================
# Read simulation output and calculate operating characteristics
# ==============================================================================

get_res_matrix = function(
    model_name,
    scn,
    N,
    level,
    results_path = results_path,
    true_psqc = true_psqc
) {

  file_model = file.path(
    results_path,
    paste0(
      model_name,
      "_Scn", scn,
      "_N", N,
      "_level", level_tag(level),
      ".txt"
    )
  )

  if (!file.exists(file_model)) {
    stop("Missing simulation result file: ", file_model)
  }

  df_model = read.table(
    file_model,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # The first simulation replicate writes the column header.
  if (
    nrow(df_model) > 0 &&
      as.character(df_model[1, 1]) == "run_ID"
  ) {
    df_model = df_model[-1, , drop = FALSE]
  }

  df_model[] = lapply(
    df_model,
    function(x) {
      suppressWarnings(as.numeric(as.character(x)))
    }
  )

  # E_rho_nu_diff_result contains nine entries for each (nu, rho) row.
  block_vars = c(
    "psi_z0",
    "psi_z1",
    "nu",
    "rho",
    "est",
    "sd",
    "q_low",
    "q_up",
    "ci_length"
  )

  K = (ncol(df_model) - 1) / length(block_vars)

  if (!is.finite(K) || K != floor(K)) {
    stop(
      "Unexpected number of columns in: ",
      file_model
    )
  }

  K = as.integer(K)

  colnames(df_model) = c(
    "run_ID",
    paste0(
      rep(block_vars, K),
      "_row",
      rep(seq_len(K), each = length(block_vars))
    )
  )

  results_mat = tibble::tibble(
    psi_z0 = numeric(K),
    psi_z1 = numeric(K),
    nu = numeric(K),
    rho = numeric(K),
    True = numeric(K),
    Est = numeric(K),
    Bias = numeric(K),
    RMSE = numeric(K),
    Cover = numeric(K),
    n_est_used = integer(K),
    n_cover_used = integer(K)
  )

  for (kk in seq_len(K)) {

    results_mat$psi_z0[kk] =
      df_model[[paste0("psi_z0_row", kk)]][1]

    results_mat$psi_z1[kk] =
      df_model[[paste0("psi_z1_row", kk)]][1]

    results_mat$nu[kk] =
      df_model[[paste0("nu_row", kk)]][1]

    results_mat$rho[kk] =
      df_model[[paste0("rho_row", kk)]][1]

    truth_index = which(
      abs(true_psqc$nu - results_mat$nu[kk]) < 1e-8 &
        abs(true_psqc$rho - results_mat$rho[kk]) < 1e-8
    )

    if (length(truth_index) == 0) {
      true_value = NA_real_
    } else {
      true_value = true_psqc$True_Diff[truth_index[1]]
    }

    vals_est =
      df_model[[paste0("est_row", kk)]]

    vals_low =
      df_model[[paste0("q_low_row", kk)]]

    vals_up =
      df_model[[paste0("q_up_row", kk)]]

    vals_length =
      df_model[[paste0("ci_length_row", kk)]]

    # Point-estimate metrics use the posterior estimate from this result file.
    mask_est =
      is.finite(vals_est) &
      abs(vals_est) < 1000

    # Coverage additionally requires a valid interval.
    mask_cover =
      mask_est &
      is.finite(vals_low) &
      is.finite(vals_up) &
      is.finite(vals_length) &
      abs(vals_low) < 1000 &
      abs(vals_up) < 1000 &
      abs(vals_length) < 1000 &
      vals_low <= vals_up

    results_mat$n_est_used[kk] =
      sum(mask_est)

    results_mat$n_cover_used[kk] =
      sum(mask_cover)

    results_mat$True[kk] =
      true_value

    if (any(mask_est) && is.finite(true_value)) {

      est_clean = vals_est[mask_est]

      results_mat$Est[kk] =
        mean(est_clean)

      results_mat$Bias[kk] =
        mean(est_clean - true_value)

      results_mat$RMSE[kk] =
        sqrt(mean((est_clean - true_value)^2))

    } else {

      results_mat$Est[kk] =
        NA_real_

      results_mat$Bias[kk] =
        NA_real_

      results_mat$RMSE[kk] =
        NA_real_
    }

    if (any(mask_cover) && is.finite(true_value)) {

      low_clean = vals_low[mask_cover]
      up_clean = vals_up[mask_cover]

      results_mat$Cover[kk] =
        mean(
          true_value >= low_clean &
            true_value <= up_clean
        )

    } else {

      results_mat$Cover[kk] =
        NA_real_
    }
  }

  rho_keep = rho_by_scenario[[as.character(scn)]]

  if (is.null(rho_keep)) {
    stop("No rho specification found for scenario ", scn)
  }

  results_mat = results_mat %>%
    filter(
      vapply(
        rho,
        function(x) {
          any(abs(x - rho_keep) < 1e-8)
        },
        logical(1)
      )
    ) %>%
    arrange(nu, rho)

  if (nrow(results_mat) != length(nu_vec) * length(rho_keep)) {
    stop(
      "The expected (nu, rho) combinations were not found in ",
      file_model
    )
  }

  results_mat
}

# ==============================================================================
# Common performance-table helpers
# ==============================================================================

prepare_result_block = function(df) {
  df %>%
    transmute(
      nu,
      rho,
      True,
      Bias,
      RMSE,
      Cover
    )
}

make_four_setting_block = function(
    model_name,
    left_500,
    left_1500,
    right_500,
    right_1500
) {

  left_500 = prepare_result_block(left_500)
  left_1500 = prepare_result_block(left_1500)
  right_500 = prepare_result_block(right_500)
  right_1500 = prepare_result_block(right_1500)

  stopifnot(
    all(left_500$nu == left_1500$nu),
    all(left_500$rho == left_1500$rho),
    all(left_500$nu == right_500$nu),
    all(left_500$rho == right_500$rho),
    all(left_500$nu == right_1500$nu),
    all(left_500$rho == right_1500$rho)
  )

  tibble::tibble(
    Model = model_name,
    nu = left_500$nu,
    rho = left_500$rho,
    True = left_500$True,

    Bias1 = left_500$Bias,
    RMSE1 = left_500$RMSE,
    CP1 = left_500$Cover,

    Bias2 = left_1500$Bias,
    RMSE2 = left_1500$RMSE,
    CP2 = left_1500$Cover,

    Bias3 = right_500$Bias,
    RMSE3 = right_500$RMSE,
    CP3 = right_500$Cover,

    Bias4 = right_1500$Bias,
    RMSE4 = right_1500$RMSE,
    CP4 = right_1500$Cover
  )
}

format_four_setting_display = function(df) {
  df %>%
    mutate(
      True = fmt2(True),

      Bias1 = fmt2(Bias1),
      RMSE1 = fmt2(RMSE1),
      CP1 = fmt_cp(CP1),

      Bias2 = fmt2(Bias2),
      RMSE2 = fmt2(RMSE2),
      CP2 = fmt_cp(CP2),

      Bias3 = fmt2(Bias3),
      RMSE3 = fmt2(RMSE3),
      CP3 = fmt_cp(CP3),

      Bias4 = fmt2(Bias4),
      RMSE4 = fmt2(RMSE4),
      CP4 = fmt_cp(CP4)
    )
}

blank_model_and_nu = function(df) {

  df %>%
    mutate(
      Model = as.character(Model),
      nu = as.character(nu)
    ) %>%
    group_by(Model) %>%
    mutate(
      nu = ifelse(
        duplicated(nu),
        "",
        nu
      )
    ) %>%
    ungroup() %>%
    mutate(
      Model = ifelse(
        duplicated(Model),
        "",
        Model
      )
    )
}

make_performance_table_html = function(
    df,
    caption,
    censoring_labels,
    sample_size_labels,
    collapse_nu = TRUE
) {

  display_df = format_four_setting_display(df)

  clean_header_names = c(
    "Model", "\u03bd", "\u03c1", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  )

  header_n = stats::setNames(
    c(4, 3, 3, 3, 3),
    c(" ", sample_size_labels)
  )

  header_censoring = stats::setNames(
    c(4, 6, 6),
    c(" ", censoring_labels)
  )

  tbl = kableExtra::kbl(
    display_df,
    format = "html",
    escape = FALSE,
    align = "c",
    col.names = clean_header_names,
    caption = caption
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "condensed"),
      full_width = FALSE,
      font_size = 11
    ) %>%
    kableExtra::add_header_above(header_n) %>%
    kableExtra::add_header_above(header_censoring) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7, background = "grey90") %>%
    kableExtra::column_spec(8:10, background = "grey80") %>%
    kableExtra::column_spec(11:13, background = "grey90") %>%
    kableExtra::column_spec(14:16, background = "grey80")

  if (collapse_nu) {
    tbl = tbl %>%
      kableExtra::collapse_rows(
        columns = c(1, 2),
        valign = "top"
      )
  } else {
    tbl = tbl %>%
      kableExtra::collapse_rows(
        columns = 1,
        valign = "top"
      )
  }

  htmltools::html_print(
    htmltools::HTML(as.character(tbl)),
    viewer = getOption("viewer")
  )

  invisible(tbl)
}

make_performance_table_latex = function(
    df,
    caption,
    label,
    censoring_labels,
    sample_size_labels,
    font_size = 7
) {

  display_df = format_four_setting_display(df) %>%
    blank_model_and_nu()

  clean_header_names = c(
    "Model", "$\\nu$", "$\\rho$", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  )

  header_n = stats::setNames(
    c(4, 3, 3, 3, 3),
    c(" ", sample_size_labels)
  )

  header_censoring = stats::setNames(
    c(4, 6, 6),
    c(" ", censoring_labels)
  )

  first_model_rows = which(display_df$Model != "")

  out = kableExtra::kbl(
    display_df,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", rep("c", ncol(display_df) - 1)),
    col.names = clean_header_names,
    caption = caption,
    label = label,
    linesep = ""
  ) %>%
    kableExtra::kable_styling(
      latex_options = c("scale_down"),
      full_width = FALSE,
      font_size = font_size
    ) %>%
    kableExtra::add_header_above(
      header_n,
      escape = FALSE
    ) %>%
    kableExtra::add_header_above(
      header_censoring,
      escape = FALSE
    ) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7, background = "#f2f2f2") %>%
    kableExtra::column_spec(8:10, background = "#d9d9d9") %>%
    kableExtra::column_spec(11:13, background = "#f2f2f2") %>%
    kableExtra::column_spec(14:16, background = "#d9d9d9")

  if (length(first_model_rows) > 1) {
    out = out %>%
      kableExtra::row_spec(
        first_model_rows[-1] - 1,
        extra_latex_after =
          "\n\\cmidrule(l{3pt}r{3pt}){1-16}\n"
      )
  }

  out = as.character(out)

  out = sub(
    "\\\\begin\\{table\\}",
    "\\\\begin{table}[tbp]",
    out
  )

  out
}

# ==============================================================================
# Load simulation results used in Tables 1, S2, and S3
# ==============================================================================

# ------------------------------------------------------------------------------
# Nominal 95% CrIs: N = 500 and N = 1500
# Used in Main-text Table 1
# ------------------------------------------------------------------------------

edpmm_s1_500 =
  get_res_matrix("EDPMM", 1, 500, main_level)

edpmm_s1_1500 =
  get_res_matrix("EDPMM", 1, 1500, main_level)

edpmm_s2_500 =
  get_res_matrix("EDPMM", 2, 500, main_level)

edpmm_s2_1500 =
  get_res_matrix("EDPMM", 2, 1500, main_level)

edpmm_s3_500 =
  get_res_matrix("EDPMM", 3, 500, main_level)

edpmm_s3_1500 =
  get_res_matrix("EDPMM", 3, 1500, main_level)

edpmm_s4_500 =
  get_res_matrix("EDPMM", 4, 500, main_level)

edpmm_s4_1500 =
  get_res_matrix("EDPMM", 4, 1500, main_level)


dpmm_s1_500 =
  get_res_matrix("DPMM", 1, 500, main_level)

dpmm_s1_1500 =
  get_res_matrix("DPMM", 1, 1500, main_level)

dpmm_s2_500 =
  get_res_matrix("DPMM", 2, 500, main_level)

dpmm_s2_1500 =
  get_res_matrix("DPMM", 2, 1500, main_level)

dpmm_s3_500 =
  get_res_matrix("DPMM", 3, 500, main_level)

dpmm_s3_1500 =
  get_res_matrix("DPMM", 3, 1500, main_level)

dpmm_s4_500 =
  get_res_matrix("DPMM", 4, 500, main_level)

dpmm_s4_1500 =
  get_res_matrix("DPMM", 4, 1500, main_level)


# ------------------------------------------------------------------------------
# Nominal 95% CrIs: N = 5000
# Used in Supplementary Table S2
# ------------------------------------------------------------------------------

edpmm_s3_5000 =
  get_res_matrix("EDPMM", 3, 5000, main_level)

edpmm_s4_5000 =
  get_res_matrix("EDPMM", 4, 5000, main_level)

dpmm_s3_5000 =
  get_res_matrix("DPMM", 3, 5000, main_level)

dpmm_s4_5000 =
  get_res_matrix("DPMM", 4, 5000, main_level)


# ------------------------------------------------------------------------------
# Nominal 99% CrIs: N = 500 and N = 1500
# Used in Supplementary Table S3
# ------------------------------------------------------------------------------

edpmm_s3_500_99 =
  get_res_matrix("EDPMM", 3, 500, supp_99_level)

edpmm_s3_1500_99 =
  get_res_matrix("EDPMM", 3, 1500, supp_99_level)

edpmm_s4_500_99 =
  get_res_matrix("EDPMM", 4, 500, supp_99_level)

edpmm_s4_1500_99 =
  get_res_matrix("EDPMM", 4, 1500, supp_99_level)


dpmm_s3_500_99 =
  get_res_matrix("DPMM", 3, 500, supp_99_level)

dpmm_s3_1500_99 =
  get_res_matrix("DPMM", 3, 1500, supp_99_level)

dpmm_s4_500_99 =
  get_res_matrix("DPMM", 4, 500, supp_99_level)

dpmm_s4_1500_99 =
  get_res_matrix("DPMM", 4, 1500, supp_99_level)


# ==============================================================================
# Main-text Table 1
#
# Simulation performance of the EDPMM and DPMM for estimating the marginal
# PSQC. The table reports N = 500 and N = 1500 under approximately 20%, 40%,
# 60%, and 80% right censoring using nominal 95% credible intervals.
#
# Quantile levels:
#   20% and 40% censoring: rho = 0.3 and 0.6
#   60% and 80% censoring: rho = 0.1 and 0.2
# ==============================================================================

table1_top = bind_rows(
  make_four_setting_block(
    "EDPMM",
    edpmm_s1_500,
    edpmm_s1_1500,
    edpmm_s2_500,
    edpmm_s2_1500
  ),
  make_four_setting_block(
    "DPMM",
    dpmm_s1_500,
    dpmm_s1_1500,
    dpmm_s2_500,
    dpmm_s2_1500
  )
)

table1_bottom = bind_rows(
  make_four_setting_block(
    "EDPMM",
    edpmm_s3_500,
    edpmm_s3_1500,
    edpmm_s4_500,
    edpmm_s4_1500
  ),
  make_four_setting_block(
    "DPMM",
    dpmm_s3_500,
    dpmm_s3_1500,
    dpmm_s4_500,
    dpmm_s4_1500
  )
)

table1_caption = paste0(
  "Simulation performance of the EDPMM and DPMM for estimating the marginal ",
  "PSQC at landmark time $\\nu$ and quantile level $\\rho$. Results are ",
  "reported for $N\\in\\{500,1500\\}$ and approximately 20\\%, 40\\%, ",
  "60\\%, and 80\\% right censoring. For each setting, the table gives the ",
  "true PSQC, empirical bias, RMSE, and empirical CP of the nominal 95\\% ",
  "credible interval across 1000 simulated datasets."
)

# HTML viewer: combine the 20%/40% and 60%/80% panels into one table.
table1_top_display = format_four_setting_display(table1_top)
table1_bottom_display = format_four_setting_display(table1_bottom)

table1_col_names = c(
  "Model", "\u03bd", "\u03c1", "True",
  "Bias", "RMSE", "CP",
  "Bias", "RMSE", "CP",
  "Bias", "RMSE", "CP",
  "Bias", "RMSE", "CP"
)

table1_top_html = kableExtra::kbl(
  table1_top_display,
  format = "html",
  escape = FALSE,
  align = "c",
  col.names = table1_col_names,
  caption = "Main-text Table 1"
) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE,
    font_size = 11
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "N = 500" = 3,
      "N = 1500" = 3,
      "N = 500" = 3,
      "N = 1500" = 3
    )
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "20%" = 6,
      "40%" = 6
    )
  ) %>%
  kableExtra::column_spec(5:7, background = "grey90") %>%
  kableExtra::column_spec(8:10, background = "grey80") %>%
  kableExtra::column_spec(11:13, background = "grey90") %>%
  kableExtra::column_spec(14:16, background = "grey80") %>%
  kableExtra::collapse_rows(
    columns = c(1, 2),
    valign = "top"
  )

table1_bottom_html = kableExtra::kbl(
  table1_bottom_display,
  format = "html",
  escape = FALSE,
  align = "c",
  col.names = table1_col_names
) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE,
    font_size = 11
  ) %>%
  kableExtra::column_spec(5:7, background = "grey90") %>%
  kableExtra::column_spec(8:10, background = "grey80") %>%
  kableExtra::column_spec(11:13, background = "grey90") %>%
  kableExtra::column_spec(14:16, background = "grey80") %>%
  kableExtra::collapse_rows(
    columns = c(1, 2),
    valign = "top"
  )

extract_tbody_rows = function(tbl_html) {

  doc = xml2::read_html(
    paste0(
      "<html><body>",
      as.character(tbl_html),
      "</body></html>"
    )
  )

  tbody = xml2::xml_find_first(
    doc,
    ".//tbody"
  )

  paste(
    vapply(
      xml2::xml_children(tbody),
      as.character,
      character(1)
    ),
    collapse = "\n"
  )
}

table1_top_html_text = as.character(table1_top_html)
table1_bottom_rows = extract_tbody_rows(table1_bottom_html)

table1_middle_rows = paste0(
  "\n<tr>",
  '<td colspan="16" style="padding:0;height:8px;',
  'border-top:3px solid #444;background:white;"></td>',
  "</tr>",
  "\n<tr>",
  '<td colspan="4" style="border:none;background:white;"></td>',
  '<td colspan="6" style="text-align:center;font-weight:bold;',
  'font-size:120%;border:none;border-bottom:1px solid #bbb;',
  'background:white;">60%</td>',
  '<td colspan="6" style="text-align:center;font-weight:bold;',
  'font-size:120%;border:none;border-bottom:1px solid #bbb;',
  'background:white;">80%</td>',
  "</tr>\n"
)

table1_viewer_html = sub(
  "</tbody>",
  paste0(
    table1_middle_rows,
    table1_bottom_rows,
    "\n</tbody>"
  ),
  table1_top_html_text,
  fixed = TRUE
)

htmltools::html_print(
  htmltools::HTML(table1_viewer_html),
  viewer = getOption("viewer")
)

# The manuscript presents the two censoring panels as one table.
table1_latex_df = bind_rows(
  table1_top %>% mutate(.panel = 1L),
  table1_bottom %>% mutate(.panel = 2L)
) %>%
  format_four_setting_display() %>%
  mutate(
    Model = as.character(Model),
    nu = as.character(nu)
  ) %>%
  group_by(.panel, Model) %>%
  mutate(
    nu = ifelse(
      duplicated(nu),
      "",
      nu
    )
  ) %>%
  ungroup() %>%
  group_by(.panel) %>%
  mutate(
    Model = ifelse(
      duplicated(Model),
      "",
      Model
    )
  ) %>%
  ungroup() %>%
  select(-.panel)

table1_mid_header = paste0(
  "\n\\midrule\n",
  "\\multicolumn{4}{c}{} & ",
  "\\multicolumn{6}{c}{60\\%} & ",
  "\\multicolumn{6}{c}{80\\%} \\\\\n",
  "\\cmidrule(lr){5-10}\\cmidrule(lr){11-16}\n"
)

table1_tex = kableExtra::kbl(
  table1_latex_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = c("l", rep("c", 15)),
  col.names = c(
    "Model", "$\\nu$", "$\\rho$", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  ),
  caption = table1_caption,
  label = "tableEDPqrl_sim_main",
  linesep = ""
) %>%
  kableExtra::kable_styling(
    latex_options = c("scale_down"),
    full_width = FALSE,
    font_size = 7
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "$N=500$" = 3,
      "$N=1500$" = 3,
      "$N=500$" = 3,
      "$N=1500$" = 3
    ),
    escape = FALSE
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "20\\%" = 6,
      "40\\%" = 6
    ),
    escape = FALSE
  ) %>%
  kableExtra::column_spec(4, background = "white") %>%
  kableExtra::column_spec(5:7, background = "#f2f2f2") %>%
  kableExtra::column_spec(8:10, background = "#d9d9d9") %>%
  kableExtra::column_spec(11:13, background = "#f2f2f2") %>%
  kableExtra::column_spec(14:16, background = "#d9d9d9") %>%
  kableExtra::row_spec(
    8,
    extra_latex_after =
      "\n\\cmidrule(l{3pt}r{3pt}){1-16}\n"
  ) %>%
  kableExtra::row_spec(
    16,
    extra_latex_after = table1_mid_header
  ) %>%
  kableExtra::row_spec(
    24,
    extra_latex_after =
      "\n\\cmidrule(l{3pt}r{3pt}){1-16}\n"
  ) %>%
  as.character()

table1_tex = sub(
  "\\\\begin\\{table\\}",
  "\\\\begin{table}[tbp]",
  table1_tex
)

cat(
  "\n\n================ MAIN-TEXT TABLE 1 ================\n"
)
cat(table1_tex, sep = "\n")


# ==============================================================================
# Supplementary Table S1
#
# Observed failure-time information in the TAS principal stratum across
# censoring settings.
#
# Each table cell reports:
#   percentage of TAS members with an observed failure
#   (percentage of the full population that is in the TAS stratum and has
#    an observed failure)
#
# The same Monte Carlo population is used for all censoring settings so that
# the TAS principal-stratum probabilities are identical across settings.
# ==============================================================================

tas_nu_vec = c(0, 1, 2, 3)

censoring_settings = tibble::tibble(
  censoring = c(20, 40, 60, 80),
  cstar = c(3.20, 1.79, 0.53, -0.95)
)

# Values used for Supplementary Table S1 in the manuscript
n_mc = 1e6
chunk_size = 5e4
trt_intercept = 0.2
mc_seed = 20260711

estimate_tas_information = function(
    censoring,
    cstar,
    nu_vec = c(0, 1, 2, 3),
    n_mc = 1e6,
    chunk_size = 5e4,
    trt_intercept = 0.2,
    seed = 1
) {

  set.seed(seed)

  n_principal = numeric(length(nu_vec))
  n_observed_failure = numeric(length(nu_vec))

  n_total = 0
  n_event = 0

  n_chunks = ceiling(n_mc / chunk_size)

  for (bb in seq_len(n_chunks)) {

    n_current = min(
      chunk_size,
      n_mc - n_total
    )

    dat = generate_data(
      Scn = 1,
      n = n_current,
      cstar = cstar,
      trt_intercept = trt_intercept
    )

    Dobs = as.logical(dat$Dobs)

    n_event =
      n_event + sum(Dobs)

    for (jj in seq_along(nu_vec)) {

      nu = nu_vec[jj]

      tas_indicator = (
        dat$Y0 > nu &
          dat$Y1 > nu
      )

      observed_failure_indicator = (
        tas_indicator &
          Dobs
      )

      n_principal[jj] =
        n_principal[jj] + sum(tas_indicator)

      n_observed_failure[jj] =
        n_observed_failure[jj] +
        sum(observed_failure_indicator)
    }

    n_total =
      n_total + n_current
  }

  tibble::tibble(
    censoring = censoring,
    cstar = cstar,
    nu = nu_vec,
    n_total = n_total,
    n_event = n_event,
    n_principal = n_principal,
    n_observed_failure = n_observed_failure,

    actual_censoring_percent =
      100 * (1 - n_event / n_total),

    principal_percent =
      100 * n_principal / n_total,

    observed_failure_percent =
      100 * n_observed_failure / n_total,

    conditional_failure_percent =
      100 * n_observed_failure / n_principal
  )
}

tas_information = purrr::pmap_dfr(
  censoring_settings,
  function(censoring, cstar) {

    estimate_tas_information(
      censoring = censoring,
      cstar = cstar,
      nu_vec = tas_nu_vec,
      n_mc = n_mc,
      chunk_size = chunk_size,
      trt_intercept = trt_intercept,
      seed = mc_seed
    )
  }
)

# Check that the same TAS population was used for every censoring setting.
tas_count_check = tas_information %>%
  group_by(nu) %>%
  summarise(
    n_distinct_principal =
      n_distinct(n_principal),
    .groups = "drop"
  )

if (any(tas_count_check$n_distinct_principal != 1)) {
  stop(
    "TAS counts differ across censoring settings. ",
    "Check that the same Monte Carlo seed is used."
  )
}

table_s1_df = tas_information %>%
  mutate(
    cell = sprintf(
      "%.1f (%.1f)",
      conditional_failure_percent,
      observed_failure_percent
    )
  ) %>%
  select(
    censoring,
    nu,
    cell
  ) %>%
  pivot_wider(
    names_from = nu,
    values_from = cell,
    names_prefix = "nu_"
  ) %>%
  arrange(censoring) %>%
  select(
    censoring,
    nu_0,
    nu_1,
    nu_2,
    nu_3
  )

table_s1_caption = paste0(
  "Observed failure-time information in the TAS principal stratum across ",
  "censoring settings. Each entry gives the percentage of TAS members with ",
  "an observed failure, $100\\Pr(D=1\\mid Y^{0}>\\nu,Y^{1}>\\nu)$. The ",
  "percentage of the full Monte Carlo population that both belongs to the TAS ",
  "principal stratum and has an observed failure, ",
  "$100\\Pr(Y^{0}>\\nu,Y^{1}>\\nu,D=1)$, is shown in parentheses. These two ",
  "quantities are related by ",
  "$\\Pr(Y^{0}>\\nu,Y^{1}>\\nu,D=1)=",
  "\\Pr(D=1\\mid Y^{0}>\\nu,Y^{1}>\\nu)",
  "\\Pr(Y^{0}>\\nu,Y^{1}>\\nu)$. The TAS principal-stratum probabilities, ",
  "$100\\Pr(Y^{0}>\\nu,Y^{1}>\\nu)$, are 100.0\\%, 66.3\\%, 47.5\\%, and ",
  "36.3\\% at $\\nu=0,1,2,$ and 3, respectively. The Monte Carlo population ",
  "size is 1,000,000 for each censoring setting."
)

table_s1_viewer = kableExtra::kbl(
  table_s1_df,
  format = "html",
  escape = FALSE,
  align = "c",
  col.names = c(
    "Censoring (%)",
    "\u03bd = 0",
    "\u03bd = 1",
    "\u03bd = 2",
    "\u03bd = 3"
  ),
  caption = "Supplementary Table S1"
) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE,
    font_size = 11
  )

htmltools::html_print(
  htmltools::HTML(as.character(table_s1_viewer)),
  viewer = getOption("viewer")
)

table_s1_tex = kableExtra::kbl(
  table_s1_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "ccccc",
  col.names = c(
    "Censoring (\\%)",
    "$\\nu=0$",
    "$\\nu=1$",
    "$\\nu=2$",
    "$\\nu=3$"
  ),
  caption = table_s1_caption,
  label = "tableEDPqrl_SIM_pop_supp",
  linesep = ""
) %>%
  kableExtra::kable_styling(
    full_width = FALSE,
    font_size = 8
  ) %>%
  as.character()

table_s1_tex = sub(
  "\\\\begin\\{table\\}",
  "\\\\begin{table}[tbp]",
  table_s1_tex
)

cat(
  "\n\n================ SUPPLEMENTARY TABLE S1 ================\n"
)
cat(table_s1_tex, sep = "\n")


# ==============================================================================
# Supplementary Table S2
#
# Additional simulation performance under approximately 60% and 80% censoring
# with N = 5000 using nominal 95% credible intervals.
# ==============================================================================

make_two_setting_block = function(
    model_name,
    left,
    right
) {

  left = prepare_result_block(left)
  right = prepare_result_block(right)

  stopifnot(
    all(left$nu == right$nu),
    all(left$rho == right$rho)
  )

  tibble::tibble(
    Model = model_name,
    nu = left$nu,
    rho = left$rho,
    True = left$True,

    Bias1 = left$Bias,
    RMSE1 = left$RMSE,
    CP1 = left$Cover,

    Bias2 = right$Bias,
    RMSE2 = right$RMSE,
    CP2 = right$Cover
  )
}

format_two_setting_display = function(df) {
  df %>%
    mutate(
      True = fmt2(True),

      Bias1 = fmt2(Bias1),
      RMSE1 = fmt2(RMSE1),
      CP1 = fmt_cp(CP1),

      Bias2 = fmt2(Bias2),
      RMSE2 = fmt2(RMSE2),
      CP2 = fmt_cp(CP2)
    )
}

blank_model_and_nu_two = function(df) {

  blank_model_and_nu(df)
}

table_s2_df = bind_rows(
  make_two_setting_block(
    "EDPMM",
    edpmm_s3_5000,
    edpmm_s4_5000
  ),
  make_two_setting_block(
    "DPMM",
    dpmm_s3_5000,
    dpmm_s4_5000
  )
)

table_s2_caption = paste0(
  "Additional simulation performance of the EDPMM and DPMM for estimating ",
  "the marginal PSQC under approximately 60\\% and 80\\% censoring with ",
  "$N=5000$. The table reports the true PSQC, empirical bias, RMSE, and CP ",
  "of nominal 95\\% credible intervals at landmark time $\\nu$ and quantile ",
  "level $\\rho$."
)

table_s2_display = format_two_setting_display(
  table_s2_df
)

table_s2_viewer = kableExtra::kbl(
  table_s2_display,
  format = "html",
  escape = FALSE,
  align = "c",
  col.names = c(
    "Model", "\u03bd", "\u03c1", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  ),
  caption = "Supplementary Table S2"
) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE,
    font_size = 11
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "60%" = 3,
      "80%" = 3
    )
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "N = 5000" = 6
    )
  ) %>%
  kableExtra::column_spec(5:7, background = "grey90") %>%
  kableExtra::column_spec(8:10, background = "grey80") %>%
  kableExtra::collapse_rows(
    columns = c(1, 2),
    valign = "top"
  )

htmltools::html_print(
  htmltools::HTML(as.character(table_s2_viewer)),
  viewer = getOption("viewer")
)

table_s2_latex_df = table_s2_display %>%
  blank_model_and_nu_two()

table_s2_tex = kableExtra::kbl(
  table_s2_latex_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = c("l", rep("c", 9)),
  col.names = c(
    "Model", "$\\nu$", "$\\rho$", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  ),
  caption = table_s2_caption,
  label = "tableEDPqrl_sim34_suppN5000",
  linesep = ""
) %>%
  kableExtra::kable_styling(
    latex_options = c("scale_down"),
    full_width = FALSE,
    font_size = 8
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "60\\%" = 3,
      "80\\%" = 3
    ),
    escape = FALSE
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "$N=5000$" = 6
    ),
    escape = FALSE
  ) %>%
  kableExtra::column_spec(5:7, background = "#f2f2f2") %>%
  kableExtra::column_spec(8:10, background = "#d9d9d9") %>%
  kableExtra::row_spec(
    nrow(table_s2_df) / 2,
    extra_latex_after =
      "\n\\cmidrule(l{3pt}r{3pt}){1-10}\n"
  ) %>%
  as.character()

table_s2_tex = sub(
  "\\\\begin\\{table\\}",
  "\\\\begin{table}[tbp]",
  table_s2_tex
)

cat(
  "\n\n================ SUPPLEMENTARY TABLE S2 ================\n"
)
cat(table_s2_tex, sep = "\n")


# ==============================================================================
# Supplementary Table S3
#
# Simulation performance under approximately 60% and 80% censoring using
# nominal 99% credible intervals for N = 500 and N = 1500.
# ==============================================================================

table_s3_df = bind_rows(
  make_four_setting_block(
    "EDPMM",
    edpmm_s3_500_99,
    edpmm_s3_1500_99,
    edpmm_s4_500_99,
    edpmm_s4_1500_99
  ),
  make_four_setting_block(
    "DPMM",
    dpmm_s3_500_99,
    dpmm_s3_1500_99,
    dpmm_s4_500_99,
    dpmm_s4_1500_99
  )
)

table_s3_caption = paste0(
  "Simulation performance under approximately 60\\% and 80\\% censoring ",
  "using nominal 99\\% credible intervals for ",
  "$N\\in\\{500,1500\\}$. The table reports the true marginal PSQC, ",
  "empirical bias, RMSE, and CP at landmark time $\\nu$ and quantile level ",
  "$\\rho$."
)

table_s3_viewer = make_performance_table_html(
  df = table_s3_df,
  caption = "Supplementary Table S3",
  censoring_labels = c("60%", "80%"),
  sample_size_labels = c(
    "N = 500", "N = 1500",
    "N = 500", "N = 1500"
  )
)

table_s3_tex = make_performance_table_latex(
  df = table_s3_df,
  caption = table_s3_caption,
  label = "tableEDPqrl_sim34_supp",
  censoring_labels = c("60\\%", "80\\%"),
  sample_size_labels = c(
    "$N=500$", "$N=1500$",
    "$N=500$", "$N=1500$"
  ),
  font_size = 7
)

cat(
  "\n\n================ SUPPLEMENTARY TABLE S3 ================\n"
)
cat(table_s3_tex, sep = "\n")


# ==============================================================================
# End
# ==============================================================================

cat(
  "\n\nCreated manuscript table objects:\n",
  "  table1_tex\n",
  "  table_s1_tex\n",
  "  table_s2_tex\n",
  "  table_s3_tex\n",
  sep = ""
)
Library
/
EDPqrl
/
PSQC_simulation_tables_replication.R


# ==============================================================================
# PSQC simulation tables
#
# Manuscript table order:
#   1. Main-text Table 1
#   2. Supplementary Table S1
#   3. Supplementary Table S2
#   4. Supplementary Table S3
#
# This script generates only the simulation tables reported in the manuscript.
# ==============================================================================

# ------------------------------------------------------------------------------
# R packages
# ------------------------------------------------------------------------------
library(Rcpp)
library(RcppArmadillo)

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

library(knitr)
library(kableExtra)
library(htmltools)
library(xml2)

# ------------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
source_dir = ""
results_dir = "Results"
truth_dir = ""

source_path = if (nzchar(source_dir)) source_dir else "."
results_path = if (nzchar(results_dir)) results_dir else "."
truth_path = if (nzchar(truth_dir)) truth_dir else "."

# ------------------------------------------------------------------------------
# Source files
# ------------------------------------------------------------------------------
source(file.path(source_path, "EDPqrl_r.R"))
Rcpp::sourceCpp(file.path(source_path, "EDPqrl_cpp.cpp"))

# ------------------------------------------------------------------------------
# General settings
# ------------------------------------------------------------------------------
CE = TRUE
esttype = "mean"

nu_vec = c(0, 1, 2, 3)

# Quantile levels reported in the manuscript:
#   Scenarios 1 and 2: approximately 20% and 40% censoring
#   Scenarios 3 and 4: approximately 60% and 80% censoring
rho_by_scenario = list(
  `1` = c(0.30, 0.60),
  `2` = c(0.30, 0.60),
  `3` = c(0.10, 0.20),
  `4` = c(0.10, 0.20)
)

# The simulation replication script saves:
#   nominal 95% CrI results with "_level5"
#   nominal 99% CrI results with "_level1"
main_level = 0.05
supp_99_level = 0.01

# ==============================================================================
# Optional preprocessing: calculate the common true PSQC
#
# This calculation is not a manuscript table. It creates True_PSQC.txt, which
# is used as the population truth when calculating bias, RMSE, and coverage.
#
# The original calculation is computationally intensive:
#   10,000 Monte Carlo datasets, each with 10,000 observations.
#
# Leave regenerate_true_psqc = FALSE when True_PSQC.txt is already supplied.
# ==============================================================================

regenerate_true_psqc = FALSE
true_psqc_file = file.path(truth_path, "True_PSQC.txt")

if (regenerate_true_psqc) {

  set.seed(1)

  truth_rho_vec = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60)
  truth_nu_vec = c(0, 1, 2, 3)

  n_sim = 1e4
  n_sample = 1e4

  param_grid = expand.grid(
    rho = truth_rho_vec,
    nu = truth_nu_vec
  )

  param_grid = param_grid[, c("nu", "rho")]
  n_params = nrow(param_grid)

  # The true PSQC depends only on Y0 and Y1 and is common across scenarios.
  scn_id = 1

  sum_results = matrix(
    0,
    nrow = n_params,
    ncol = 3
  )

  valid_results = integer(n_params)

  cat("Starting true PSQC calculation...\n")

  for (ii in seq_len(n_sim)) {

    dat = generate_data(
      Scn = scn_id,
      n = n_sample
    )

    for (jj in seq_len(n_params)) {

      curr_nu = param_grid$nu[jj]
      curr_rho = param_grid$rho[jj]

      # Time-specific always-survivor (TAS) principal stratum
      idx_tas = (
        dat$Y0 > curr_nu &
          dat$Y1 > curr_nu
      )

      Y0_tas = dat$Y0[idx_tas] - curr_nu
      Y1_tas = dat$Y1[idx_tas] - curr_nu

      Q0 = if (length(Y0_tas) > 0) {
        as.numeric(
          quantile(
            Y0_tas,
            probs = curr_rho,
            names = FALSE
          )
        )
      } else {
        NA_real_
      }

      Q1 = if (length(Y1_tas) > 0) {
        as.numeric(
          quantile(
            Y1_tas,
            probs = curr_rho,
            names = FALSE
          )
        )
      } else {
        NA_real_
      }

      diff = Q1 - Q0

      if (all(is.finite(c(Q0, Q1, diff)))) {
        sum_results[jj, 1] = sum_results[jj, 1] + Q0
        sum_results[jj, 2] = sum_results[jj, 2] + Q1
        sum_results[jj, 3] = sum_results[jj, 3] + diff

        valid_results[jj] = valid_results[jj] + 1L
      }
    }

    if (ii %% 1000 == 0) {
      cat(sprintf("Completed %d / %d\n", ii, n_sim))
    }
  }

  if (any(valid_results == 0)) {
    stop(
      "No valid true PSQC values were obtained for at least one ",
      "(nu, rho) combination."
    )
  }

  truth_base = cbind(
    param_grid,
    sum_results / valid_results
  )

  colnames(truth_base) = c(
    "nu",
    "rho",
    "True_Res_Y0",
    "True_Res_Y1",
    "True_Diff"
  )

  truth_base = as.data.frame(
    truth_base,
    stringsAsFactors = FALSE
  )

  truth_base = truth_base[
    order(truth_base$nu, truth_base$rho),
    ,
    drop = FALSE
  ]

  rownames(truth_base) = NULL

  write.table(
    truth_base,
    file = true_psqc_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  cat("Saved: ", true_psqc_file, "\n", sep = "")
}

if (!file.exists(true_psqc_file)) {
  stop(
    "Missing true PSQC file: ", true_psqc_file, "\n",
    "Supply True_PSQC.txt or set regenerate_true_psqc = TRUE."
  )
}

true_psqc = read.table(
  true_psqc_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

required_truth_cols = c("nu", "rho", "True_Diff")

if (!all(required_truth_cols %in% names(true_psqc))) {
  stop(
    "True_PSQC.txt must contain columns: ",
    paste(required_truth_cols, collapse = ", ")
  )
}

# ==============================================================================
# Formatting helpers
# ==============================================================================

fmt2 = function(x) {
  x = as.numeric(x)
  x[is.finite(x) & abs(x) < 0.005] = 0
  sprintf("%.2f", x)
}

fmt_cp = function(x) {
  sprintf("%.1f", 100 * as.numeric(x))
}

level_tag = function(level) {
  format(
    100 * level,
    trim = TRUE,
    scientific = FALSE
  )
}

# ==============================================================================
# Read simulation output and calculate operating characteristics
# ==============================================================================

get_res_matrix = function(
    model_name,
    scn,
    N,
    level,
    results_path = results_path,
    true_psqc = true_psqc
) {

  file_model = file.path(
    results_path,
    paste0(
      model_name,
      "_Scn", scn,
      "_N", N,
      "_level", level_tag(level),
      ".txt"
    )
  )

  if (!file.exists(file_model)) {
    stop("Missing simulation result file: ", file_model)
  }

  df_model = read.table(
    file_model,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # The first simulation replicate writes the column header.
  if (
    nrow(df_model) > 0 &&
      as.character(df_model[1, 1]) == "run_ID"
  ) {
    df_model = df_model[-1, , drop = FALSE]
  }

  df_model[] = lapply(
    df_model,
    function(x) {
      suppressWarnings(as.numeric(as.character(x)))
    }
  )

  # E_rho_nu_diff_result contains nine entries for each (nu, rho) row.
  block_vars = c(
    "psi_z0",
    "psi_z1",
    "nu",
    "rho",
    "est",
    "sd",
    "q_low",
    "q_up",
    "ci_length"
  )

  K = (ncol(df_model) - 1) / length(block_vars)

  if (!is.finite(K) || K != floor(K)) {
    stop(
      "Unexpected number of columns in: ",
      file_model
    )
  }

  K = as.integer(K)

  colnames(df_model) = c(
    "run_ID",
    paste0(
      rep(block_vars, K),
      "_row",
      rep(seq_len(K), each = length(block_vars))
    )
  )

  results_mat = tibble::tibble(
    psi_z0 = numeric(K),
    psi_z1 = numeric(K),
    nu = numeric(K),
    rho = numeric(K),
    True = numeric(K),
    Est = numeric(K),
    Bias = numeric(K),
    RMSE = numeric(K),
    Cover = numeric(K),
    n_est_used = integer(K),
    n_cover_used = integer(K)
  )

  for (kk in seq_len(K)) {

    results_mat$psi_z0[kk] =
      df_model[[paste0("psi_z0_row", kk)]][1]

    results_mat$psi_z1[kk] =
      df_model[[paste0("psi_z1_row", kk)]][1]

    results_mat$nu[kk] =
      df_model[[paste0("nu_row", kk)]][1]

    results_mat$rho[kk] =
      df_model[[paste0("rho_row", kk)]][1]

    truth_index = which(
      abs(true_psqc$nu - results_mat$nu[kk]) < 1e-8 &
        abs(true_psqc$rho - results_mat$rho[kk]) < 1e-8
    )

    if (length(truth_index) == 0) {
      true_value = NA_real_
    } else {
      true_value = true_psqc$True_Diff[truth_index[1]]
    }

    vals_est =
      df_model[[paste0("est_row", kk)]]

    vals_low =
      df_model[[paste0("q_low_row", kk)]]

    vals_up =
      df_model[[paste0("q_up_row", kk)]]

    vals_length =
      df_model[[paste0("ci_length_row", kk)]]

    # Point-estimate metrics use the posterior estimate from this result file.
    mask_est =
      is.finite(vals_est) &
      abs(vals_est) < 1000

    # Coverage additionally requires a valid interval.
    mask_cover =
      mask_est &
      is.finite(vals_low) &
      is.finite(vals_up) &
      is.finite(vals_length) &
      abs(vals_low) < 1000 &
      abs(vals_up) < 1000 &
      abs(vals_length) < 1000 &
      vals_low <= vals_up

    results_mat$n_est_used[kk] =
      sum(mask_est)

    results_mat$n_cover_used[kk] =
      sum(mask_cover)

    results_mat$True[kk] =
      true_value

    if (any(mask_est) && is.finite(true_value)) {

      est_clean = vals_est[mask_est]

      results_mat$Est[kk] =
        mean(est_clean)

      results_mat$Bias[kk] =
        mean(est_clean - true_value)

      results_mat$RMSE[kk] =
        sqrt(mean((est_clean - true_value)^2))

    } else {

      results_mat$Est[kk] =
        NA_real_

      results_mat$Bias[kk] =
        NA_real_

      results_mat$RMSE[kk] =
        NA_real_
    }

    if (any(mask_cover) && is.finite(true_value)) {

      low_clean = vals_low[mask_cover]
      up_clean = vals_up[mask_cover]

      results_mat$Cover[kk] =
        mean(
          true_value >= low_clean &
            true_value <= up_clean
        )

    } else {

      results_mat$Cover[kk] =
        NA_real_
    }
  }

  rho_keep = rho_by_scenario[[as.character(scn)]]

  if (is.null(rho_keep)) {
    stop("No rho specification found for scenario ", scn)
  }

  results_mat = results_mat %>%
    filter(
      vapply(
        rho,
        function(x) {
          any(abs(x - rho_keep) < 1e-8)
        },
        logical(1)
      )
    ) %>%
    arrange(nu, rho)

  if (nrow(results_mat) != length(nu_vec) * length(rho_keep)) {
    stop(
      "The expected (nu, rho) combinations were not found in ",
      file_model
    )
  }

  results_mat
}

# ==============================================================================
# Common performance-table helpers
# ==============================================================================

prepare_result_block = function(df) {
  df %>%
    transmute(
      nu,
      rho,
      True,
      Bias,
      RMSE,
      Cover
    )
}

make_four_setting_block = function(
    model_name,
    left_500,
    left_1500,
    right_500,
    right_1500
) {

  left_500 = prepare_result_block(left_500)
  left_1500 = prepare_result_block(left_1500)
  right_500 = prepare_result_block(right_500)
  right_1500 = prepare_result_block(right_1500)

  stopifnot(
    all(left_500$nu == left_1500$nu),
    all(left_500$rho == left_1500$rho),
    all(left_500$nu == right_500$nu),
    all(left_500$rho == right_500$rho),
    all(left_500$nu == right_1500$nu),
    all(left_500$rho == right_1500$rho)
  )

  tibble::tibble(
    Model = model_name,
    nu = left_500$nu,
    rho = left_500$rho,
    True = left_500$True,

    Bias1 = left_500$Bias,
    RMSE1 = left_500$RMSE,
    CP1 = left_500$Cover,

    Bias2 = left_1500$Bias,
    RMSE2 = left_1500$RMSE,
    CP2 = left_1500$Cover,

    Bias3 = right_500$Bias,
    RMSE3 = right_500$RMSE,
    CP3 = right_500$Cover,

    Bias4 = right_1500$Bias,
    RMSE4 = right_1500$RMSE,
    CP4 = right_1500$Cover
  )
}

format_four_setting_display = function(df) {
  df %>%
    mutate(
      True = fmt2(True),

      Bias1 = fmt2(Bias1),
      RMSE1 = fmt2(RMSE1),
      CP1 = fmt_cp(CP1),

      Bias2 = fmt2(Bias2),
      RMSE2 = fmt2(RMSE2),
      CP2 = fmt_cp(CP2),

      Bias3 = fmt2(Bias3),
      RMSE3 = fmt2(RMSE3),
      CP3 = fmt_cp(CP3),

      Bias4 = fmt2(Bias4),
      RMSE4 = fmt2(RMSE4),
      CP4 = fmt_cp(CP4)
    )
}

blank_model_and_nu = function(df) {

  df %>%
    mutate(
      Model = as.character(Model),
      nu = as.character(nu)
    ) %>%
    group_by(Model) %>%
    mutate(
      nu = ifelse(
        duplicated(nu),
        "",
        nu
      )
    ) %>%
    ungroup() %>%
    mutate(
      Model = ifelse(
        duplicated(Model),
        "",
        Model
      )
    )
}

make_performance_table_html = function(
    df,
    caption,
    censoring_labels,
    sample_size_labels,
    collapse_nu = TRUE
) {

  display_df = format_four_setting_display(df)

  clean_header_names = c(
    "Model", "\u03bd", "\u03c1", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  )

  header_n = stats::setNames(
    c(4, 3, 3, 3, 3),
    c(" ", sample_size_labels)
  )

  header_censoring = stats::setNames(
    c(4, 6, 6),
    c(" ", censoring_labels)
  )

  tbl = kableExtra::kbl(
    display_df,
    format = "html",
    escape = FALSE,
    align = "c",
    col.names = clean_header_names,
    caption = caption
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "condensed"),
      full_width = FALSE,
      font_size = 11
    ) %>%
    kableExtra::add_header_above(header_n) %>%
    kableExtra::add_header_above(header_censoring) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7, background = "grey90") %>%
    kableExtra::column_spec(8:10, background = "grey80") %>%
    kableExtra::column_spec(11:13, background = "grey90") %>%
    kableExtra::column_spec(14:16, background = "grey80")

  if (collapse_nu) {
    tbl = tbl %>%
      kableExtra::collapse_rows(
        columns = c(1, 2),
        valign = "top"
      )
  } else {
    tbl = tbl %>%
      kableExtra::collapse_rows(
        columns = 1,
        valign = "top"
      )
  }

  htmltools::html_print(
    htmltools::HTML(as.character(tbl)),
    viewer = getOption("viewer")
  )

  invisible(tbl)
}

make_performance_table_latex = function(
    df,
    caption,
    label,
    censoring_labels,
    sample_size_labels,
    font_size = 7
) {

  display_df = format_four_setting_display(df) %>%
    blank_model_and_nu()

  clean_header_names = c(
    "Model", "$\\nu$", "$\\rho$", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  )

  header_n = stats::setNames(
    c(4, 3, 3, 3, 3),
    c(" ", sample_size_labels)
  )

  header_censoring = stats::setNames(
    c(4, 6, 6),
    c(" ", censoring_labels)
  )

  first_model_rows = which(display_df$Model != "")

  out = kableExtra::kbl(
    display_df,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", rep("c", ncol(display_df) - 1)),
    col.names = clean_header_names,
    caption = caption,
    label = label,
    linesep = ""
  ) %>%
    kableExtra::kable_styling(
      latex_options = c("scale_down"),
      full_width = FALSE,
      font_size = font_size
    ) %>%
    kableExtra::add_header_above(
      header_n,
      escape = FALSE
    ) %>%
    kableExtra::add_header_above(
      header_censoring,
      escape = FALSE
    ) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7, background = "#f2f2f2") %>%
    kableExtra::column_spec(8:10, background = "#d9d9d9") %>%
    kableExtra::column_spec(11:13, background = "#f2f2f2") %>%
    kableExtra::column_spec(14:16, background = "#d9d9d9")

  if (length(first_model_rows) > 1) {
    out = out %>%
      kableExtra::row_spec(
        first_model_rows[-1] - 1,
        extra_latex_after =
          "\n\\cmidrule(l{3pt}r{3pt}){1-16}\n"
      )
  }

  out = as.character(out)

  out = sub(
    "\\\\begin\\{table\\}",
    "\\\\begin{table}[tbp]",
    out
  )

  out
}

# ==============================================================================
# Load simulation results used in Tables 1, S2, and S3
# ==============================================================================

# ------------------------------------------------------------------------------
# Nominal 95% CrIs: N = 500 and N = 1500
# Used in Main-text Table 1
# ------------------------------------------------------------------------------

edpmm_s1_500 =
  get_res_matrix("EDPMM", 1, 500, main_level)

edpmm_s1_1500 =
  get_res_matrix("EDPMM", 1, 1500, main_level)

edpmm_s2_500 =
  get_res_matrix("EDPMM", 2, 500, main_level)

edpmm_s2_1500 =
  get_res_matrix("EDPMM", 2, 1500, main_level)

edpmm_s3_500 =
  get_res_matrix("EDPMM", 3, 500, main_level)

edpmm_s3_1500 =
  get_res_matrix("EDPMM", 3, 1500, main_level)

edpmm_s4_500 =
  get_res_matrix("EDPMM", 4, 500, main_level)

edpmm_s4_1500 =
  get_res_matrix("EDPMM", 4, 1500, main_level)


dpmm_s1_500 =
  get_res_matrix("DPMM", 1, 500, main_level)

dpmm_s1_1500 =
  get_res_matrix("DPMM", 1, 1500, main_level)

dpmm_s2_500 =
  get_res_matrix("DPMM", 2, 500, main_level)

dpmm_s2_1500 =
  get_res_matrix("DPMM", 2, 1500, main_level)

dpmm_s3_500 =
  get_res_matrix("DPMM", 3, 500, main_level)

dpmm_s3_1500 =
  get_res_matrix("DPMM", 3, 1500, main_level)

dpmm_s4_500 =
  get_res_matrix("DPMM", 4, 500, main_level)

dpmm_s4_1500 =
  get_res_matrix("DPMM", 4, 1500, main_level)


# ------------------------------------------------------------------------------
# Nominal 95% CrIs: N = 5000
# Used in Supplementary Table S2
# ------------------------------------------------------------------------------

edpmm_s3_5000 =
  get_res_matrix("EDPMM", 3, 5000, main_level)

edpmm_s4_5000 =
  get_res_matrix("EDPMM", 4, 5000, main_level)

dpmm_s3_5000 =
  get_res_matrix("DPMM", 3, 5000, main_level)

dpmm_s4_5000 =
  get_res_matrix("DPMM", 4, 5000, main_level)


# ------------------------------------------------------------------------------
# Nominal 99% CrIs: N = 500 and N = 1500
# Used in Supplementary Table S3
# ------------------------------------------------------------------------------

edpmm_s3_500_99 =
  get_res_matrix("EDPMM", 3, 500, supp_99_level)

edpmm_s3_1500_99 =
  get_res_matrix("EDPMM", 3, 1500, supp_99_level)

edpmm_s4_500_99 =
  get_res_matrix("EDPMM", 4, 500, supp_99_level)

edpmm_s4_1500_99 =
  get_res_matrix("EDPMM", 4, 1500, supp_99_level)


dpmm_s3_500_99 =
  get_res_matrix("DPMM", 3, 500, supp_99_level)

dpmm_s3_1500_99 =
  get_res_matrix("DPMM", 3, 1500, supp_99_level)

dpmm_s4_500_99 =
  get_res_matrix("DPMM", 4, 500, supp_99_level)

dpmm_s4_1500_99 =
  get_res_matrix("DPMM", 4, 1500, supp_99_level)


# ==============================================================================
# Main-text Table 1
#
# Simulation performance of the EDPMM and DPMM for estimating the marginal
# PSQC. The table reports N = 500 and N = 1500 under approximately 20%, 40%,
# 60%, and 80% right censoring using nominal 95% credible intervals.
#
# Quantile levels:
#   20% and 40% censoring: rho = 0.3 and 0.6
#   60% and 80% censoring: rho = 0.1 and 0.2
# ==============================================================================

table1_top = bind_rows(
  make_four_setting_block(
    "EDPMM",
    edpmm_s1_500,
    edpmm_s1_1500,
    edpmm_s2_500,
    edpmm_s2_1500
  ),
  make_four_setting_block(
    "DPMM",
    dpmm_s1_500,
    dpmm_s1_1500,
    dpmm_s2_500,
    dpmm_s2_1500
  )
)

table1_bottom = bind_rows(
  make_four_setting_block(
    "EDPMM",
    edpmm_s3_500,
    edpmm_s3_1500,
    edpmm_s4_500,
    edpmm_s4_1500
  ),
  make_four_setting_block(
    "DPMM",
    dpmm_s3_500,
    dpmm_s3_1500,
    dpmm_s4_500,
    dpmm_s4_1500
  )
)

table1_caption = paste0(
  "Simulation performance of the EDPMM and DPMM for estimating the marginal ",
  "PSQC at landmark time $\\nu$ and quantile level $\\rho$. Results are ",
  "reported for $N\\in\\{500,1500\\}$ and approximately 20\\%, 40\\%, ",
  "60\\%, and 80\\% right censoring. For each setting, the table gives the ",
  "true PSQC, empirical bias, RMSE, and empirical CP of the nominal 95\\% ",
  "credible interval across 1000 simulated datasets."
)

# HTML viewer: combine the 20%/40% and 60%/80% panels into one table.
table1_top_display = format_four_setting_display(table1_top)
table1_bottom_display = format_four_setting_display(table1_bottom)

table1_col_names = c(
  "Model", "\u03bd", "\u03c1", "True",
  "Bias", "RMSE", "CP",
  "Bias", "RMSE", "CP",
  "Bias", "RMSE", "CP",
  "Bias", "RMSE", "CP"
)

table1_top_html = kableExtra::kbl(
  table1_top_display,
  format = "html",
  escape = FALSE,
  align = "c",
  col.names = table1_col_names,
  caption = "Main-text Table 1"
) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE,
    font_size = 11
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "N = 500" = 3,
      "N = 1500" = 3,
      "N = 500" = 3,
      "N = 1500" = 3
    )
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "20%" = 6,
      "40%" = 6
    )
  ) %>%
  kableExtra::column_spec(5:7, background = "grey90") %>%
  kableExtra::column_spec(8:10, background = "grey80") %>%
  kableExtra::column_spec(11:13, background = "grey90") %>%
  kableExtra::column_spec(14:16, background = "grey80") %>%
  kableExtra::collapse_rows(
    columns = c(1, 2),
    valign = "top"
  )

table1_bottom_html = kableExtra::kbl(
  table1_bottom_display,
  format = "html",
  escape = FALSE,
  align = "c",
  col.names = table1_col_names
) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE,
    font_size = 11
  ) %>%
  kableExtra::column_spec(5:7, background = "grey90") %>%
  kableExtra::column_spec(8:10, background = "grey80") %>%
  kableExtra::column_spec(11:13, background = "grey90") %>%
  kableExtra::column_spec(14:16, background = "grey80") %>%
  kableExtra::collapse_rows(
    columns = c(1, 2),
    valign = "top"
  )

extract_tbody_rows = function(tbl_html) {

  doc = xml2::read_html(
    paste0(
      "<html><body>",
      as.character(tbl_html),
      "</body></html>"
    )
  )

  tbody = xml2::xml_find_first(
    doc,
    ".//tbody"
  )

  paste(
    vapply(
      xml2::xml_children(tbody),
      as.character,
      character(1)
    ),
    collapse = "\n"
  )
}

table1_top_html_text = as.character(table1_top_html)
table1_bottom_rows = extract_tbody_rows(table1_bottom_html)

table1_middle_rows = paste0(
  "\n<tr>",
  '<td colspan="16" style="padding:0;height:8px;',
  'border-top:3px solid #444;background:white;"></td>',
  "</tr>",
  "\n<tr>",
  '<td colspan="4" style="border:none;background:white;"></td>',
  '<td colspan="6" style="text-align:center;font-weight:bold;',
  'font-size:120%;border:none;border-bottom:1px solid #bbb;',
  'background:white;">60%</td>',
  '<td colspan="6" style="text-align:center;font-weight:bold;',
  'font-size:120%;border:none;border-bottom:1px solid #bbb;',
  'background:white;">80%</td>',
  "</tr>\n"
)

table1_viewer_html = sub(
  "</tbody>",
  paste0(
    table1_middle_rows,
    table1_bottom_rows,
    "\n</tbody>"
  ),
  table1_top_html_text,
  fixed = TRUE
)

htmltools::html_print(
  htmltools::HTML(table1_viewer_html),
  viewer = getOption("viewer")
)

# The manuscript presents the two censoring panels as one table.
table1_latex_df = bind_rows(
  table1_top %>% mutate(.panel = 1L),
  table1_bottom %>% mutate(.panel = 2L)
) %>%
  format_four_setting_display() %>%
  mutate(
    Model = as.character(Model),
    nu = as.character(nu)
  ) %>%
  group_by(.panel, Model) %>%
  mutate(
    nu = ifelse(
      duplicated(nu),
      "",
      nu
    )
  ) %>%
  ungroup() %>%
  group_by(.panel) %>%
  mutate(
    Model = ifelse(
      duplicated(Model),
      "",
      Model
    )
  ) %>%
  ungroup() %>%
  select(-.panel)

table1_mid_header = paste0(
  "\n\\midrule\n",
  "\\multicolumn{4}{c}{} & ",
  "\\multicolumn{6}{c}{60\\%} & ",
  "\\multicolumn{6}{c}{80\\%} \\\\\n",
  "\\cmidrule(lr){5-10}\\cmidrule(lr){11-16}\n"
)

table1_tex = kableExtra::kbl(
  table1_latex_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = c("l", rep("c", 15)),
  col.names = c(
    "Model", "$\\nu$", "$\\rho$", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  ),
  caption = table1_caption,
  label = "tableEDPqrl_sim_main",
  linesep = ""
) %>%
  kableExtra::kable_styling(
    latex_options = c("scale_down"),
    full_width = FALSE,
    font_size = 7
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "$N=500$" = 3,
      "$N=1500$" = 3,
      "$N=500$" = 3,
      "$N=1500$" = 3
    ),
    escape = FALSE
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "20\\%" = 6,
      "40\\%" = 6
    ),
    escape = FALSE
  ) %>%
  kableExtra::column_spec(4, background = "white") %>%
  kableExtra::column_spec(5:7, background = "#f2f2f2") %>%
  kableExtra::column_spec(8:10, background = "#d9d9d9") %>%
  kableExtra::column_spec(11:13, background = "#f2f2f2") %>%
  kableExtra::column_spec(14:16, background = "#d9d9d9") %>%
  kableExtra::row_spec(
    8,
    extra_latex_after =
      "\n\\cmidrule(l{3pt}r{3pt}){1-16}\n"
  ) %>%
  kableExtra::row_spec(
    16,
    extra_latex_after = table1_mid_header
  ) %>%
  kableExtra::row_spec(
    24,
    extra_latex_after =
      "\n\\cmidrule(l{3pt}r{3pt}){1-16}\n"
  ) %>%
  as.character()

table1_tex = sub(
  "\\\\begin\\{table\\}",
  "\\\\begin{table}[tbp]",
  table1_tex
)

cat(
  "\n\n================ MAIN-TEXT TABLE 1 ================\n"
)
cat(table1_tex, sep = "\n")


# ==============================================================================
# Supplementary Table S1
#
# Observed failure-time information in the TAS principal stratum across
# censoring settings.
#
# Each table cell reports:
#   percentage of TAS members with an observed failure
#   (percentage of the full population that is in the TAS stratum and has
#    an observed failure)
#
# The same Monte Carlo population is used for all censoring settings so that
# the TAS principal-stratum probabilities are identical across settings.
# ==============================================================================

tas_nu_vec = c(0, 1, 2, 3)

censoring_settings = tibble::tibble(
  censoring = c(20, 40, 60, 80),
  cstar = c(3.20, 1.79, 0.53, -0.95)
)

# Values used for Supplementary Table S1 in the manuscript
n_mc = 1e6
chunk_size = 5e4
trt_intercept = 0.2
mc_seed = 20260711

estimate_tas_information = function(
    censoring,
    cstar,
    nu_vec = c(0, 1, 2, 3),
    n_mc = 1e6,
    chunk_size = 5e4,
    trt_intercept = 0.2,
    seed = 1
) {

  set.seed(seed)

  n_principal = numeric(length(nu_vec))
  n_observed_failure = numeric(length(nu_vec))

  n_total = 0
  n_event = 0

  n_chunks = ceiling(n_mc / chunk_size)

  for (bb in seq_len(n_chunks)) {

    n_current = min(
      chunk_size,
      n_mc - n_total
    )

    dat = generate_data(
      Scn = 1,
      n = n_current,
      cstar = cstar,
      trt_intercept = trt_intercept
    )

    Dobs = as.logical(dat$Dobs)

    n_event =
      n_event + sum(Dobs)

    for (jj in seq_along(nu_vec)) {

      nu = nu_vec[jj]

      tas_indicator = (
        dat$Y0 > nu &
          dat$Y1 > nu
      )

      observed_failure_indicator = (
        tas_indicator &
          Dobs
      )

      n_principal[jj] =
        n_principal[jj] + sum(tas_indicator)

      n_observed_failure[jj] =
        n_observed_failure[jj] +
        sum(observed_failure_indicator)
    }

    n_total =
      n_total + n_current
  }

  tibble::tibble(
    censoring = censoring,
    cstar = cstar,
    nu = nu_vec,
    n_total = n_total,
    n_event = n_event,
    n_principal = n_principal,
    n_observed_failure = n_observed_failure,

    actual_censoring_percent =
      100 * (1 - n_event / n_total),

    principal_percent =
      100 * n_principal / n_total,

    observed_failure_percent =
      100 * n_observed_failure / n_total,

    conditional_failure_percent =
      100 * n_observed_failure / n_principal
  )
}

tas_information = purrr::pmap_dfr(
  censoring_settings,
  function(censoring, cstar) {

    estimate_tas_information(
      censoring = censoring,
      cstar = cstar,
      nu_vec = tas_nu_vec,
      n_mc = n_mc,
      chunk_size = chunk_size,
      trt_intercept = trt_intercept,
      seed = mc_seed
    )
  }
)

# Check that the same TAS population was used for every censoring setting.
tas_count_check = tas_information %>%
  group_by(nu) %>%
  summarise(
    n_distinct_principal =
      n_distinct(n_principal),
    .groups = "drop"
  )

if (any(tas_count_check$n_distinct_principal != 1)) {
  stop(
    "TAS counts differ across censoring settings. ",
    "Check that the same Monte Carlo seed is used."
  )
}

table_s1_df = tas_information %>%
  mutate(
    cell = sprintf(
      "%.1f (%.1f)",
      conditional_failure_percent,
      observed_failure_percent
    )
  ) %>%
  select(
    censoring,
    nu,
    cell
  ) %>%
  pivot_wider(
    names_from = nu,
    values_from = cell,
    names_prefix = "nu_"
  ) %>%
  arrange(censoring) %>%
  select(
    censoring,
    nu_0,
    nu_1,
    nu_2,
    nu_3
  )

table_s1_caption = paste0(
  "Observed failure-time information in the TAS principal stratum across ",
  "censoring settings. Each entry gives the percentage of TAS members with ",
  "an observed failure, $100\\Pr(D=1\\mid Y^{0}>\\nu,Y^{1}>\\nu)$. The ",
  "percentage of the full Monte Carlo population that both belongs to the TAS ",
  "principal stratum and has an observed failure, ",
  "$100\\Pr(Y^{0}>\\nu,Y^{1}>\\nu,D=1)$, is shown in parentheses. These two ",
  "quantities are related by ",
  "$\\Pr(Y^{0}>\\nu,Y^{1}>\\nu,D=1)=",
  "\\Pr(D=1\\mid Y^{0}>\\nu,Y^{1}>\\nu)",
  "\\Pr(Y^{0}>\\nu,Y^{1}>\\nu)$. The TAS principal-stratum probabilities, ",
  "$100\\Pr(Y^{0}>\\nu,Y^{1}>\\nu)$, are 100.0\\%, 66.3\\%, 47.5\\%, and ",
  "36.3\\% at $\\nu=0,1,2,$ and 3, respectively. The Monte Carlo population ",
  "size is 1,000,000 for each censoring setting."
)

table_s1_viewer = kableExtra::kbl(
  table_s1_df,
  format = "html",
  escape = FALSE,
  align = "c",
  col.names = c(
    "Censoring (%)",
    "\u03bd = 0",
    "\u03bd = 1",
    "\u03bd = 2",
    "\u03bd = 3"
  ),
  caption = "Supplementary Table S1"
) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE,
    font_size = 11
  )

htmltools::html_print(
  htmltools::HTML(as.character(table_s1_viewer)),
  viewer = getOption("viewer")
)

table_s1_tex = kableExtra::kbl(
  table_s1_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "ccccc",
  col.names = c(
    "Censoring (\\%)",
    "$\\nu=0$",
    "$\\nu=1$",
    "$\\nu=2$",
    "$\\nu=3$"
  ),
  caption = table_s1_caption,
  label = "tableEDPqrl_SIM_pop_supp",
  linesep = ""
) %>%
  kableExtra::kable_styling(
    full_width = FALSE,
    font_size = 8
  ) %>%
  as.character()

table_s1_tex = sub(
  "\\\\begin\\{table\\}",
  "\\\\begin{table}[tbp]",
  table_s1_tex
)

cat(
  "\n\n================ SUPPLEMENTARY TABLE S1 ================\n"
)
cat(table_s1_tex, sep = "\n")


# ==============================================================================
# Supplementary Table S2
#
# Additional simulation performance under approximately 60% and 80% censoring
# with N = 5000 using nominal 95% credible intervals.
# ==============================================================================

make_two_setting_block = function(
    model_name,
    left,
    right
) {

  left = prepare_result_block(left)
  right = prepare_result_block(right)

  stopifnot(
    all(left$nu == right$nu),
    all(left$rho == right$rho)
  )

  tibble::tibble(
    Model = model_name,
    nu = left$nu,
    rho = left$rho,
    True = left$True,

    Bias1 = left$Bias,
    RMSE1 = left$RMSE,
    CP1 = left$Cover,

    Bias2 = right$Bias,
    RMSE2 = right$RMSE,
    CP2 = right$Cover
  )
}

format_two_setting_display = function(df) {
  df %>%
    mutate(
      True = fmt2(True),

      Bias1 = fmt2(Bias1),
      RMSE1 = fmt2(RMSE1),
      CP1 = fmt_cp(CP1),

      Bias2 = fmt2(Bias2),
      RMSE2 = fmt2(RMSE2),
      CP2 = fmt_cp(CP2)
    )
}

blank_model_and_nu_two = function(df) {

  blank_model_and_nu(df)
}

table_s2_df = bind_rows(
  make_two_setting_block(
    "EDPMM",
    edpmm_s3_5000,
    edpmm_s4_5000
  ),
  make_two_setting_block(
    "DPMM",
    dpmm_s3_5000,
    dpmm_s4_5000
  )
)

table_s2_caption = paste0(
  "Additional simulation performance of the EDPMM and DPMM for estimating ",
  "the marginal PSQC under approximately 60\\% and 80\\% censoring with ",
  "$N=5000$. The table reports the true PSQC, empirical bias, RMSE, and CP ",
  "of nominal 95\\% credible intervals at landmark time $\\nu$ and quantile ",
  "level $\\rho$."
)

table_s2_display = format_two_setting_display(
  table_s2_df
)

table_s2_viewer = kableExtra::kbl(
  table_s2_display,
  format = "html",
  escape = FALSE,
  align = "c",
  col.names = c(
    "Model", "\u03bd", "\u03c1", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  ),
  caption = "Supplementary Table S2"
) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("hover", "condensed"),
    full_width = FALSE,
    font_size = 11
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "60%" = 3,
      "80%" = 3
    )
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "N = 5000" = 6
    )
  ) %>%
  kableExtra::column_spec(5:7, background = "grey90") %>%
  kableExtra::column_spec(8:10, background = "grey80") %>%
  kableExtra::collapse_rows(
    columns = c(1, 2),
    valign = "top"
  )

htmltools::html_print(
  htmltools::HTML(as.character(table_s2_viewer)),
  viewer = getOption("viewer")
)

table_s2_latex_df = table_s2_display %>%
  blank_model_and_nu_two()

table_s2_tex = kableExtra::kbl(
  table_s2_latex_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = c("l", rep("c", 9)),
  col.names = c(
    "Model", "$\\nu$", "$\\rho$", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  ),
  caption = table_s2_caption,
  label = "tableEDPqrl_sim34_suppN5000",
  linesep = ""
) %>%
  kableExtra::kable_styling(
    latex_options = c("scale_down"),
    full_width = FALSE,
    font_size = 8
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "60\\%" = 3,
      "80\\%" = 3
    ),
    escape = FALSE
  ) %>%
  kableExtra::add_header_above(
    c(
      " " = 4,
      "$N=5000$" = 6
    ),
    escape = FALSE
  ) %>%
  kableExtra::column_spec(5:7, background = "#f2f2f2") %>%
  kableExtra::column_spec(8:10, background = "#d9d9d9") %>%
  kableExtra::row_spec(
    nrow(table_s2_df) / 2,
    extra_latex_after =
      "\n\\cmidrule(l{3pt}r{3pt}){1-10}\n"
  ) %>%
  as.character()

table_s2_tex = sub(
  "\\\\begin\\{table\\}",
  "\\\\begin{table}[tbp]",
  table_s2_tex
)

cat(
  "\n\n================ SUPPLEMENTARY TABLE S2 ================\n"
)
cat(table_s2_tex, sep = "\n")


# ==============================================================================
# Supplementary Table S3
#
# Simulation performance under approximately 60% and 80% censoring using
# nominal 99% credible intervals for N = 500 and N = 1500.
# ==============================================================================

table_s3_df = bind_rows(
  make_four_setting_block(
    "EDPMM",
    edpmm_s3_500_99,
    edpmm_s3_1500_99,
    edpmm_s4_500_99,
    edpmm_s4_1500_99
  ),
  make_four_setting_block(
    "DPMM",
    dpmm_s3_500_99,
    dpmm_s3_1500_99,
    dpmm_s4_500_99,
    dpmm_s4_1500_99
  )
)

table_s3_caption = paste0(
  "Simulation performance under approximately 60\\% and 80\\% censoring ",
  "using nominal 99\\% credible intervals for ",
  "$N\\in\\{500,1500\\}$. The table reports the true marginal PSQC, ",
  "empirical bias, RMSE, and CP at landmark time $\\nu$ and quantile level ",
  "$\\rho$."
)

table_s3_viewer = make_performance_table_html(
  df = table_s3_df,
  caption = "Supplementary Table S3",
  censoring_labels = c("60%", "80%"),
  sample_size_labels = c(
    "N = 500", "N = 1500",
    "N = 500", "N = 1500"
  )
)

table_s3_tex = make_performance_table_latex(
  df = table_s3_df,
  caption = table_s3_caption,
  label = "tableEDPqrl_sim34_supp",
  censoring_labels = c("60\\%", "80\\%"),
  sample_size_labels = c(
    "$N=500$", "$N=1500$",
    "$N=500$", "$N=1500$"
  ),
  font_size = 7
)

cat(
  "\n\n================ SUPPLEMENTARY TABLE S3 ================\n"
)
cat(table_s3_tex, sep = "\n")


# ==============================================================================
# End
# ==============================================================================

cat(
  "\n\nCreated manuscript table objects:\n",
  "  table1_tex\n",
  "  table_s1_tex\n",
  "  table_s2_tex\n",
  "  table_s3_tex\n",
  sep = ""
)