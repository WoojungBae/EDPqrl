library(Rcpp)
library(RcppArmadillo)

library(MASS)
library(dplyr)
library(tidyverse)

library(knitr)
library(kableExtra)

# Set your directory (source code)
# dir = "/Users/WooJung/Documents/Rproject/EDPqrl/Replication"
setwd(dir)
source("EDPqrl_r.R")
sourceCpp("EDPqrl_cpp.cpp")

# Set your directory (source code)
setwd(paste0(dir,"/Results"))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
CE = TRUE
esttype = "mean"

# ==============================================================================
# Helper: Calculate Metrics for a Single Scenario
# ==============================================================================
get_res_matrix <- function(model_name, scn, N, extra = NULL) {
  file_true  <- paste0("True_Values.txt")
  if(!file.exists(file_true)) {
    warning(paste("Missing file:", file_true))
    return(NULL)
  }
  df_true_static <- read.table(file_true, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (is.null(extra)) {
    file_model <- paste0(model_name, "_Scn", scn, "_N", N, ".txt")
  } else {
    file_model <- paste0(model_name, "_Scn", scn, "_N", N, extra, ".txt")
  }
  if(!file.exists(file_model)) {
    warning(paste("Missing file:", file_model))
    return(NULL)
  }
  df_model <- read.table(file_model, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  
  if (as.character(df_model[1,1]) == "run_ID") df_model <- df_model[-1, ]
  df_model[] <- lapply(df_model, function(x) suppressWarnings(as.numeric(as.character(x))))
  
  # Determine number of blocks K dynamically (9 vars per block)
  block_vars <- c("eta_z0", "eta_z1", "nu", "rho", "est", "sd", "q25", "q95", "len")
  K <- (ncol(df_model) - 1) / length(block_vars)
  
  colnames(df_model) <- c("run_ID", paste0(rep(block_vars, K), "_row", rep(1:K, each=length(block_vars))))
  
  results_mat <- data.frame(
    eta_z0 = numeric(K),
    eta_z1 = numeric(K),
    nu     = numeric(K),
    rho    = numeric(K),
    True   = numeric(K),
    Est    = numeric(K),
    Bias   = numeric(K),
    RMSE   = numeric(K),
    Cover  = numeric(K)
  )
  
  for(k in 1:K) {
    # Extract metadata for this block
    results_mat$eta_z0[k] <- df_model[[paste0("eta_z0_row", k)]][1]
    results_mat$eta_z1[k] <- df_model[[paste0("eta_z1_row", k)]][1]
    results_mat$nu[k]     <- df_model[[paste0("nu_row", k)]][1]
    results_mat$rho[k]    <- df_model[[paste0("rho_row", k)]][1]
    
    # Dynamically match the True value based on nu and rho
    match_idx <- which(df_true_static$nu == results_mat$nu[k] & df_true_static$rho == results_mat$rho[k])
    if (length(match_idx) > 0) {
      val_true_scalar <- df_true_static$True_Diff[match_idx[1]]
    } else {
      val_true_scalar <- NA
    }
    
    vals_est <- df_model[[paste0("est_row", k)]]
    vals_low <- df_model[[paste0("q25_row", k)]]
    vals_up  <- df_model[[paste0("q95_row", k)]]
    vals_len <- df_model[[paste0("len_row", k)]]
    
    mask <- (is.finite(vals_est) & abs(vals_est) < 1000 &
               is.finite(vals_len) & abs(vals_len) < 1000)
    
    valid_idx <- which(mask)
    
    if(length(valid_idx) > 0) {
      est_clean <- vals_est[valid_idx]
      low_clean <- vals_low[valid_idx]
      up_clean  <- vals_up[valid_idx]
      
      results_mat$True[k]  <- val_true_scalar
      results_mat$Est[k]   <- mean(est_clean, na.rm = TRUE)
      results_mat$Bias[k]  <- mean(est_clean - val_true_scalar, na.rm = TRUE)
      results_mat$RMSE[k]  <- sqrt(mean((est_clean - val_true_scalar)^2, na.rm = TRUE))
      results_mat$Cover[k] <- mean(val_true_scalar >= low_clean & val_true_scalar <= up_clean, na.rm = TRUE)
      
    } else {
      results_mat[k, 5:9] <- NA
    }
  }
  
  formatted_res <- results_mat
  # Only apply string formatting to the metric columns, keep metadata numeric
  formatted_res[, 5:9] <- apply(results_mat[, 5:9], 2, function(x) sprintf("%.3f", x))
  
  return(formatted_res)
}

generate_full_comparison_table_viewer <- function(
    edpmm_s1_500, edpmm_s1_1500,
    edpmm_s2_500, edpmm_s2_1500,
    edpmm_s3_500, edpmm_s3_1500,
    edpmm_s4_500, edpmm_s4_1500,
    dpmm_s1_500,  dpmm_s1_1500,
    dpmm_s2_500,  dpmm_s2_1500,
    dpmm_s3_500,  dpmm_s3_1500,
    dpmm_s4_500,  dpmm_s4_1500,
    table_title
) {
  
  make_block <- function(
    ed_left_500, ed_left_1500,
    ed_right_500, ed_right_1500,
    dp_left_500, dp_left_1500,
    dp_right_500, dp_right_1500) {
    
    scale_cp <- function(df) {
      if ("Cover" %in% names(df)) {
        df$Cover <- 100 * as.numeric(as.character(df$Cover))
      }
      df
    }
    
    ed_left_500   <- scale_cp(ed_left_500)
    ed_left_1500  <- scale_cp(ed_left_1500)
    ed_right_500  <- scale_cp(ed_right_500)
    ed_right_1500 <- scale_cp(ed_right_1500)
    
    dp_left_500   <- scale_cp(dp_left_500)
    dp_left_1500  <- scale_cp(dp_left_1500)
    dp_right_500  <- scale_cp(dp_right_500)
    dp_right_1500 <- scale_cp(dp_right_1500)
    
    meta_labels <- ed_left_500[, c("nu", "rho")]
    metrics <- c("Bias", "RMSE", "Cover")
    
    ed_block <- data.frame(
      Model = "EDPMM",
      meta_labels,
      True = ed_left_500$True,
      ed_left_500[, metrics],
      ed_left_1500[, metrics],
      ed_right_500[, metrics],
      ed_right_1500[, metrics],
      stringsAsFactors = FALSE
    )
    
    dp_block <- data.frame(
      Model = "DPMM",
      meta_labels,
      True = dp_left_500$True,
      dp_left_500[, metrics],
      dp_left_1500[, metrics],
      dp_right_500[, metrics],
      dp_right_1500[, metrics],
      stringsAsFactors = FALSE
    )
    
    rbind(ed_block, dp_block)
  }
  
  final_top <- make_block(
    edpmm_s1_500, edpmm_s1_1500,
    edpmm_s2_500, edpmm_s2_1500,
    dpmm_s1_500,  dpmm_s1_1500,
    dpmm_s2_500,  dpmm_s2_1500
  )
  
  final_bottom <- make_block(
    edpmm_s3_500, edpmm_s3_1500,
    edpmm_s4_500, edpmm_s4_1500,
    dpmm_s3_500,  dpmm_s3_1500,
    dpmm_s4_500,  dpmm_s4_1500
  )
  
  clean_header_names <- c(
    "Model", "ν", "ρ", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  )
  
  header_top <- c(
    " " = 4,
    "N = 500" = 3, "N = 1500" = 3,
    "N = 500" = 3, "N = 1500" = 3
  )
  
  header_scn <- c(
    " " = 4,
    "20%" = 6,
    "40%" = 6
  )
  
  top_tbl <- kableExtra::kbl(
    final_top,
    format = "html",
    escape = FALSE,
    align = "c",
    col.names = clean_header_names,
    caption = table_title
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "condensed"),
      full_width = FALSE,
      font_size = 11
    ) %>%
    kableExtra::add_header_above(header_top) %>%
    kableExtra::add_header_above(header_scn) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7,   background = "grey90") %>%
    kableExtra::column_spec(8:10,  background = "grey80") %>%
    kableExtra::column_spec(11:13, background = "grey90") %>%
    kableExtra::column_spec(14:16, background = "grey80") %>%
    kableExtra::collapse_rows(columns = 1:2, valign = "top")
  
  bottom_tbl <- kableExtra::kbl(
    final_bottom,
    format = "html",
    escape = FALSE,
    align = "c",
    col.names = clean_header_names
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "condensed"),
      full_width = FALSE,
      font_size = 11
    ) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7,   background = "grey90") %>%
    kableExtra::column_spec(8:10,  background = "grey80") %>%
    kableExtra::column_spec(11:13, background = "grey90") %>%
    kableExtra::column_spec(14:16, background = "grey80") %>%
    kableExtra::collapse_rows(columns = 1:2, valign = "top")
  
  extract_tbody_rows <- function(tbl_html) {
    doc <- xml2::read_html(
      paste0("<html><body>", as.character(tbl_html), "</body></html>")
    )
    tbody <- xml2::xml_find_first(doc, ".//tbody")
    paste(vapply(xml2::xml_children(tbody), as.character, character(1)), collapse = "\n")
  }
  
  top_html <- as.character(top_tbl)
  bottom_rows <- extract_tbody_rows(bottom_tbl)
  
  middle_rows <- paste0(
    '\n<tr>',
    '<td colspan="16" style="padding:0;height:8px;border-top:3px solid #444;background:white;"></td>',
    '</tr>',
    '\n<tr>',
    '<td colspan="4" style="border:none;background:white;"></td>',
    '<td colspan="6" style="text-align:center;font-weight:bold;font-size:120%;border:none;border-bottom:1px solid #bbb;background:white;">60%</td>',
    '<td colspan="6" style="text-align:center;font-weight:bold;font-size:120%;border:none;border-bottom:1px solid #bbb;background:white;">80%</td>',
    '</tr>\n'
  )
  
  final_html <- sub(
    "</tbody>",
    paste0(middle_rows, bottom_rows, "\n</tbody>"),
    top_html,
    fixed = TRUE
  )
  
  htmltools::html_print(
    htmltools::HTML(final_html),
    viewer = getOption("viewer")
  )
  
  invisible(NULL)
}

generate_full_comparison_table_latex <- function(
    edpmm_s1_500, edpmm_s1_1500,
    edpmm_s2_500, edpmm_s2_1500,
    edpmm_s3_500, edpmm_s3_1500,
    edpmm_s4_500, edpmm_s4_1500,
    dpmm_s1_500,  dpmm_s1_1500,
    dpmm_s2_500,  dpmm_s2_1500,
    dpmm_s3_500,  dpmm_s3_1500,
    dpmm_s4_500,  dpmm_s4_1500,
    table_title
) {
  
  build_panel <- function(model_name,
                          sleft_500, sleft_1500,
                          sright_500, sright_1500) {
    
    scale_cp <- function(df) {
      if ("Cover" %in% names(df)) {
        df$Cover <- 100 * as.numeric(as.character(df$Cover))
      }
      df
    }
    
    sleft_500   <- scale_cp(sleft_500)
    sleft_1500  <- scale_cp(sleft_1500)
    sright_500  <- scale_cp(sright_500)
    sright_1500 <- scale_cp(sright_1500)
    
    meta <- sleft_500[, c("nu", "rho")]
    metrics <- c("Bias", "RMSE", "Cover")
    n <- nrow(meta)
    
    out <- data.frame(
      Model = c(model_name, rep("", n - 1)),
      nu    = c(as.character(meta$nu[1]),
                ifelse(meta$nu[-1] == meta$nu[-n], "", as.character(meta$nu[-1]))),
      rho   = meta$rho,
      True  = sleft_500$True,
      sleft_500[, metrics],
      sleft_1500[, metrics],
      sright_500[, metrics],
      sright_1500[, metrics],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    names(out) <- c(
      "Model", "nu", "rho", "True",
      "Bias1", "RMSE1", "CP1",
      "Bias2", "RMSE2", "CP2",
      "Bias3", "RMSE3", "CP3",
      "Bias4", "RMSE4", "CP4"
    )
    
    out
  }
  
  top_edpmm <- build_panel(
    "EDPMM",
    edpmm_s1_500, edpmm_s1_1500,
    edpmm_s2_500, edpmm_s2_1500
  )
  
  top_dpmm <- build_panel(
    "DPMM",
    dpmm_s1_500, dpmm_s1_1500,
    dpmm_s2_500, dpmm_s2_1500
  )
  
  bot_edpmm <- build_panel(
    "EDPMM",
    edpmm_s3_500, edpmm_s3_1500,
    edpmm_s4_500, edpmm_s4_1500
  )
  
  bot_dpmm <- build_panel(
    "DPMM",
    dpmm_s3_500, dpmm_s3_1500,
    dpmm_s4_500, dpmm_s4_1500
  )
  
  final_res <- rbind(top_edpmm, top_dpmm, bot_edpmm, bot_dpmm)
  
  clean_header_names <- c(
    "Model", "$\\nu$", "$\\rho$", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  )
  
  header_top <- c(
    " " = 4,
    "N = 500" = 3, "N = 1500" = 3,
    "N = 500" = 3, "N = 1500" = 3
  )
  
  header_scn <- c(
    " " = 4,
    "20\\%" = 6,
    "40\\%" = 6
  )
  
  escape_latex_caption <- function(x) {
    x <- gsub("\\\\", "\\\\textbackslash{}", x)
    x <- gsub("([&%#$])", "\\\\\\1", x, perl = TRUE)
    x <- gsub("_", "\\\\_", x, fixed = TRUE)
    x
  }
  
  mid_block <- paste0(
    "\n\\midrule\n",
    "\\multicolumn{4}{c}{} & \\multicolumn{6}{c}{60\\%} & \\multicolumn{6}{c}{80\\%} \\\\\n",
    "\n\\midrule\n"
  )
  
  out <- kableExtra::kbl(
    final_res,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", rep("c", ncol(final_res) - 1)),
    col.names = clean_header_names,
    caption = escape_latex_caption(table_title),
    linesep = ""
  ) %>%
    kableExtra::kable_styling(
      latex_options = c("scale_down"),
      full_width = FALSE,
      font_size = 7
    ) %>%
    kableExtra::add_header_above(header_top) %>%
    kableExtra::add_header_above(header_scn) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7,   background = "#f2f2f2") %>%
    kableExtra::column_spec(8:10,  background = "#d9d9d9") %>%
    kableExtra::column_spec(11:13, background = "#f2f2f2") %>%
    kableExtra::column_spec(14:16, background = "#d9d9d9") %>%
    kableExtra::row_spec(8,  extra_latex_after = "\n\\cmidrule(l{3pt}r{3pt}){1-16}\n") %>%
    kableExtra::row_spec(16, extra_latex_after = mid_block) %>%
    kableExtra::row_spec(24, extra_latex_after = "\n\\cmidrule(l{3pt}r{3pt}){1-16}\n")
  
  out <- as.character(out)
  out <- sub("\\\\begin\\{table\\}", "\\\\begin\\{table\\}[tbp]", out)
  
  return(cat(out))
}

# ==============================================================================
# Load Data Matrices
# ==============================================================================
# --- EDPMM Data ---
edpmm_s1_500  <- get_res_matrix("EDPMM", 1, 500)
edpmm_s1_1500 <- get_res_matrix("EDPMM", 1, 1500)
edpmm_s2_500  <- get_res_matrix("EDPMM", 2, 500)
edpmm_s2_1500 <- get_res_matrix("EDPMM", 2, 1500)
edpmm_s3_500  <- get_res_matrix("EDPMM", 3, 500)
edpmm_s3_1500 <- get_res_matrix("EDPMM", 3, 1500)
edpmm_s4_500  <- get_res_matrix("EDPMM", 4, 500)
edpmm_s4_1500 <- get_res_matrix("EDPMM", 4, 1500)

# --- DPMM Data ---
dpmm_s1_500  <- get_res_matrix("DPMM", 1, 500)
dpmm_s1_1500 <- get_res_matrix("DPMM", 1, 1500)
dpmm_s2_500  <- get_res_matrix("DPMM", 2, 500)
dpmm_s2_1500 <- get_res_matrix("DPMM", 2, 1500)
dpmm_s3_500  <- get_res_matrix("DPMM", 3, 500)
dpmm_s3_1500 <- get_res_matrix("DPMM", 3, 1500)
dpmm_s4_500  <- get_res_matrix("DPMM", 4, 500)
dpmm_s4_1500 <- get_res_matrix("DPMM", 4, 1500)

# ==============================================================================
# Generate Tables (N=500 vs N=1500)
# ==============================================================================
# --- Table: EDPMM & DPMM (Scenario 1 vs 2 across N) ---
generate_full_comparison_table_viewer(
  edpmm_s1_500, edpmm_s1_1500,
  edpmm_s2_500, edpmm_s2_1500,
  edpmm_s3_500, edpmm_s3_1500,
  edpmm_s4_500, edpmm_s4_1500,
  dpmm_s1_500,  dpmm_s1_1500,
  dpmm_s2_500,  dpmm_s2_1500,
  dpmm_s3_500,  dpmm_s3_1500,
  dpmm_s4_500,  dpmm_s4_1500,
  "Table: EDPMM and DPMM performance"
)

generate_full_comparison_table_latex(
  edpmm_s1_500, edpmm_s1_1500,
  edpmm_s2_500, edpmm_s2_1500,
  edpmm_s3_500, edpmm_s3_1500,
  edpmm_s4_500, edpmm_s4_1500,
  dpmm_s1_500,  dpmm_s1_1500,
  dpmm_s2_500,  dpmm_s2_1500,
  dpmm_s3_500,  dpmm_s3_1500,
  dpmm_s4_500,  dpmm_s4_1500,
  "Table: EDPMM and DPMM performance"
)

# ==============================================================================
# Generate Tables (N=500 vs N=1500) when alpha = 0.01
# ==============================================================================
generate_s4_extra_table_viewer <- function(
    edpmm_s4_500, edpmm_s4_1500,
    dpmm_s4_500,  dpmm_s4_1500,
    table_title = "Scenario 4, alpha = 0.01"
) {
  
  prep_one <- function(df, N_label) {
    out <- df[, c("nu", "rho", "True", "Bias", "RMSE", "Cover")]
    
    for (nm in c("True", "Bias", "RMSE", "Cover")) {
      out[[nm]] <- as.numeric(as.character(out[[nm]]))
    }
    
    out$Cover <- 100 * out$Cover
    
    names(out)[names(out) == "Bias"]  <- paste0("Bias_", N_label)
    names(out)[names(out) == "RMSE"]  <- paste0("RMSE_", N_label)
    names(out)[names(out) == "Cover"] <- paste0("CP_", N_label)
    
    out
  }
  
  combine_model <- function(df_500, df_1500, model_name) {
    x500  <- prep_one(df_500,  "500")
    x1500 <- prep_one(df_1500, "1500")[, c("nu", "rho", "Bias_1500", "RMSE_1500", "CP_1500")]
    
    out <- merge(x500, x1500, by = c("nu", "rho"), all = TRUE)
    
    out$Model <- model_name
    
    out <- out[, c(
      "Model", "nu", "rho", "True",
      "Bias_500", "RMSE_500", "CP_500",
      "Bias_1500", "RMSE_1500", "CP_1500"
    )]
    
    out[order(out$nu, out$rho), ]
  }
  
  ed_block <- combine_model(edpmm_s4_500, edpmm_s4_1500, "EDPMM")
  dp_block <- combine_model(dpmm_s4_500,  dpmm_s4_1500,  "DPMM")
  
  final_df <- rbind(ed_block, dp_block)
  display_df <- final_df
  
  cols_3dec <- c(
    "True",
    "Bias_500", "RMSE_500",
    "Bias_1500", "RMSE_1500"
  )
  
  cols_cp <- c("CP_500", "CP_1500")
  
  display_df[cols_3dec] <- lapply(display_df[cols_3dec], function(x) sprintf("%.3f", x))
  display_df[cols_cp]   <- lapply(display_df[cols_cp],   function(x) sprintf("%.1f", x))
  
  clean_header_names <- c(
    "Model", "ν", "ρ", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  )
  
  header_top <- c(
    " " = 4,
    "N = 500" = 3,
    "N = 1500" = 3
  )
  
  tbl <- kableExtra::kbl(
    display_df,
    format = "html",
    escape = FALSE,
    align = "c",
    col.names = clean_header_names,
    caption = table_title
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "condensed"),
      full_width = FALSE,
      font_size = 11
    ) %>%
    kableExtra::add_header_above(header_top) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7,  background = "grey90") %>%
    kableExtra::column_spec(8:10, background = "grey80") %>%
    kableExtra::collapse_rows(columns = 1, valign = "top")
  
  htmltools::html_print(
    htmltools::HTML(as.character(tbl)),
    viewer = getOption("viewer")
  )
  
  invisible(final_df)
}

generate_s4_extra_table_latex <- function(
    edpmm_s4_500, edpmm_s4_1500,
    dpmm_s4_500,  dpmm_s4_1500,
    table_title = "Scenario 4, alpha = 0.01"
) {
  
  build_panel <- function(model_name, df_500, df_1500) {
    
    scale_cp <- function(df) {
      if ("Cover" %in% names(df)) {
        df$Cover <- 100 * as.numeric(as.character(df$Cover))
      }
      df
    }
    
    df_500  <- scale_cp(df_500)
    df_1500 <- scale_cp(df_1500)
    
    for (nm in c("True", "Bias", "RMSE", "Cover")) {
      df_500[[nm]]  <- as.numeric(as.character(df_500[[nm]]))
      df_1500[[nm]] <- as.numeric(as.character(df_1500[[nm]]))
    }
    
    meta <- df_500[, c("nu", "rho")]
    n <- nrow(meta)
    
    out <- data.frame(
      Model = c(model_name, rep("", n - 1)),
      nu    = c(as.character(meta$nu[1]),
                ifelse(meta$nu[-1] == meta$nu[-n], "", as.character(meta$nu[-1]))),
      rho   = meta$rho,
      True  = sprintf("%.3f", df_500$True),
      Bias1 = sprintf("%.3f", df_500$Bias),
      RMSE1 = sprintf("%.3f", df_500$RMSE),
      CP1   = sprintf("%.1f", df_500$Cover),
      Bias2 = sprintf("%.3f", df_1500$Bias),
      RMSE2 = sprintf("%.3f", df_1500$RMSE),
      CP2   = sprintf("%.1f", df_1500$Cover),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    out
  }
  
  ed_block <- build_panel("EDPMM", edpmm_s4_500, edpmm_s4_1500)
  dp_block <- build_panel("DPMM",  dpmm_s4_500,  dpmm_s4_1500)
  
  final_res <- rbind(ed_block, dp_block)
  
  clean_header_names <- c(
    "Model", "$\\nu$", "$\\rho$", "True",
    "Bias", "RMSE", "CP",
    "Bias", "RMSE", "CP"
  )
  
  header_top <- c(
    " " = 4,
    "N = 500" = 3,
    "N = 1500" = 3
  )
  
  escape_latex_caption <- function(x) {
    x <- gsub("\\\\", "\\\\textbackslash{}", x)
    x <- gsub("([&%#$])", "\\\\\\1", x, perl = TRUE)
    x <- gsub("_", "\\\\_", x, fixed = TRUE)
    x
  }
  
  out <- kableExtra::kbl(
    final_res,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", rep("c", ncol(final_res) - 1)),
    col.names = clean_header_names,
    caption = escape_latex_caption(table_title),
    linesep = ""
  ) %>%
    kableExtra::kable_styling(
      latex_options = c("scale_down"),
      full_width = FALSE,
      font_size = 8
    ) %>%
    kableExtra::add_header_above(header_top) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7,  background = "#f2f2f2") %>%
    kableExtra::column_spec(8:10, background = "#d9d9d9") %>%
    kableExtra::row_spec(
      nrow(ed_block),
      extra_latex_after = "\n\\cmidrule(l{3pt}r{3pt}){1-10}\n"
    )
  
  out <- as.character(out)
  out <- sub("\\\\begin\\{table\\}", "\\\\begin\\{table\\}[tbp]", out)
  
  cat(out)
}

level = 0.01

edpmm_s4_500  <- get_res_matrix("EDPMM", 4, 500, paste0("level", 100*level))
edpmm_s4_1500 <- get_res_matrix("EDPMM", 4, 1500, paste0("level", 100*level))
dpmm_s4_500  <- get_res_matrix("DPMM", 4, 500, paste0("level", 100*level))
dpmm_s4_1500 <- get_res_matrix("DPMM", 4, 1500, paste0("level", 100*level))

generate_s4_extra_table_viewer(
  edpmm_s4_500, edpmm_s4_1500, dpmm_s4_500,  dpmm_s4_1500,
  table_title = "Scenario 4, alpha = 0.01"
)

generate_s4_extra_table_latex(
  edpmm_s4_500, edpmm_s4_1500, dpmm_s4_500,  dpmm_s4_1500,
  table_title = "Scenario 4, alpha = 0.01"
)

# ==============================================================================
# Generate Tables (N=5000)
# ==============================================================================
generate_s4_extra_table_viewer <- function(
    edpmm_s4_5000,
    dpmm_s4_5000,
    table_title = "Scenario 4, alpha = 0.01"
) {
  
  prep_one <- function(df, model_name) {
    out <- df[, c("nu", "rho", "True", "Bias", "RMSE", "Cover")]
    
    for (nm in c("True", "Bias", "RMSE", "Cover")) {
      out[[nm]] <- as.numeric(as.character(out[[nm]]))
    }
    
    out$Cover <- 100 * out$Cover
    out$Model <- model_name
    
    out <- out[, c("Model", "nu", "rho", "True", "Bias", "RMSE", "Cover")]
    out[order(out$nu, out$rho), ]
  }
  
  ed_block <- prep_one(edpmm_s4_5000, "EDPMM")
  dp_block <- prep_one(dpmm_s4_5000, "DPMM")
  
  final_df <- rbind(ed_block, dp_block)
  display_df <- final_df
  
  cols_3dec <- c("True", "Bias", "RMSE")
  cols_cp <- c("Cover")
  
  display_df[cols_3dec] <- lapply(display_df[cols_3dec], function(x) sprintf("%.3f", x))
  display_df[cols_cp]   <- lapply(display_df[cols_cp],   function(x) sprintf("%.1f", x))
  
  clean_header_names <- c("Model", "ν", "ρ", "True", "Bias", "RMSE", "CP")
  
  header_top <- c(
    " " = 4,
    "N = 5000" = 3
  )
  
  tbl <- kableExtra::kbl(
    display_df,
    format = "html",
    escape = FALSE,
    align = "c",
    col.names = clean_header_names,
    caption = table_title
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("hover", "condensed"),
      full_width = FALSE,
      font_size = 11
    ) %>%
    kableExtra::add_header_above(header_top) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7, background = "grey90") %>%
    kableExtra::collapse_rows(columns = 1, valign = "top")
  
  htmltools::html_print(
    htmltools::HTML(as.character(tbl)),
    viewer = getOption("viewer")
  )
  
  invisible(final_df)
}

generate_s4_extra_table_latex <- function(
    edpmm_s4_5000,
    dpmm_s4_5000,
    table_title = "Scenario 4, alpha = 0.01"
) {
  
  build_panel <- function(model_name, df_5000) {
    
    scale_cp <- function(df) {
      if ("Cover" %in% names(df)) {
        df$Cover <- 100 * as.numeric(as.character(df$Cover))
      }
      df
    }
    
    df_5000 <- scale_cp(df_5000)
    
    for (nm in c("True", "Bias", "RMSE", "Cover")) {
      df_5000[[nm]] <- as.numeric(as.character(df_5000[[nm]]))
    }
    
    meta <- df_5000[, c("nu", "rho")]
    n <- nrow(meta)
    
    out <- data.frame(
      Model = c(model_name, rep("", n - 1)),
      nu    = c(as.character(meta$nu[1]),
                ifelse(meta$nu[-1] == meta$nu[-n], "", as.character(meta$nu[-1]))),
      rho   = meta$rho,
      True  = sprintf("%.3f", df_5000$True),
      Bias  = sprintf("%.3f", df_5000$Bias),
      RMSE  = sprintf("%.3f", df_5000$RMSE),
      CP    = sprintf("%.1f", df_5000$Cover),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    out
  }
  
  ed_block <- build_panel("EDPMM", edpmm_s4_5000)
  dp_block <- build_panel("DPMM",  dpmm_s4_5000)
  
  final_res <- rbind(ed_block, dp_block)
  
  clean_header_names <- c("Model", "$\\nu$", "$\\rho$", "True", "Bias", "RMSE", "CP")
  
  header_top <- c(
    " " = 4,
    "N = 5000" = 3
  )
  
  escape_latex_caption <- function(x) {
    x <- gsub("\\\\", "\\\\textbackslash{}", x)
    x <- gsub("([&%#$])", "\\\\\\1", x, perl = TRUE)
    x <- gsub("_", "\\\\_", x, fixed = TRUE)
    x
  }
  
  out <- kableExtra::kbl(
    final_res,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", rep("c", ncol(final_res) - 1)),
    col.names = clean_header_names,
    caption = escape_latex_caption(table_title),
    linesep = ""
  ) %>%
    kableExtra::kable_styling(
      latex_options = c("scale_down"),
      full_width = FALSE,
      font_size = 8
    ) %>%
    kableExtra::add_header_above(header_top) %>%
    kableExtra::column_spec(4, background = "white") %>%
    kableExtra::column_spec(5:7, background = "#f2f2f2") %>%
    kableExtra::row_spec(
      nrow(ed_block),
      extra_latex_after = "\n\\cmidrule(l{3pt}r{3pt}){1-7}\n"
    )
  
  out <- as.character(out)
  out <- sub("\\\\begin\\{table\\}", "\\\\begin\\{table\\}[tbp]", out)
  
  cat(out)
}

edpmm_s4_5e3 <- get_res_matrix("EDPMM", 4, 5000)
dpmm_s4_5e3 <- get_res_matrix("DPMM", 4, 5000)

generate_s4_extra_table_viewer(
  edpmm_s4_5e3, dpmm_s4_5e3, 
  table_title = "Scenario 4; N = 5000"
)

generate_s4_extra_table_latex(
  edpmm_s4_5e3, dpmm_s4_5e3, 
  table_title = "Scenario 4; N = 5000"
)
