# Please set the working directory "ROOT_DIR" before running.
# To run only a specific simulation setting, identify the corresponding row(s) in all_config and replace "seq_len(nrow(all_config))" in line 159 with the desired configuration index or indices. For example, if the target setting corresponds to the first row of all_config, use 1; if several specific settings are required, specify their corresponding indices, such as c(1, 5, 10).
# The number of Monte Carlo replications can be controlled by the n_sim parameter in line 8. To change the number of replications, simply modify the value of n_sim. For example, n_sim <- 1000 runs 1,000 Monte Carlo replications for each selected simulation setting.
ROOT_DIR <- "Z:/User/Documents/paper2"
setwd(ROOT_DIR)
source("dhsic.test.R")
source("01_code_functions.R")
n_sim <- 1000
base_seed <- 123
d_ratio <- 0.2
rho_list <- c(0, 0.2, 0.5, 0.8)
scenario_list <- c(1, 2, 3, 4)
n_list <- c(200, 500, 1000)
base_config <- expand.grid(
  rho = rho_list,
  scenario = scenario_list,
  n = n_list,
  p = 20
)
extra_config_1 <- expand.grid(
  rho = rho_list,
  scenario = scenario_list,
  n = 200,
  p = 100
)
extra_config_2 <- expand.grid(
  rho = rho_list,
  scenario = scenario_list,
  n = 500,
  p = 200
)
all_config <- rbind(base_config, extra_config_1, extra_config_2)
all_config$cfg_id <- seq_len(nrow(all_config))
out_root <- file.path(ROOT_DIR, "simulation_results")
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)
run_single_rep <- function(cfg, rep_i) {
  unique_seed <- base_seed + rep_i
  set.seed(unique_seed)
  p_curr <- cfg$p
  res <- tryCatch({
    data <- generate_data_flexible(
      n = cfg$n,
      p = p_curr,
      d = d_ratio,
      rho = cfg$rho,
      scenario = as.character(cfg$scenario)
    )
    X <- data$X
    T_val <- data$T
    Y <- data$Y
    full_df <- data.frame(Trt = T_val, Y = Y, X)
    vars <- paste0("X", 1:p_curr)
    res_HSIC_ADL <- numeric(p_curr)
    res_HSIC_Y <- numeric(p_curr)
    res_OAL <- numeric(p_curr)
    res_GOAL <- numeric(p_curr)
    res_GOALDER <- numeric(p_curr)
    class_HSIC_ADL <- rep(0, p_curr)
    class_HSIC_Y <- rep(0, p_curr)
    names(class_HSIC_ADL) <- vars
    names(class_HSIC_Y) <- vars
    idx <- get_keep_idx(X, T_val, Y)
    if (length(idx) > 0) {
      sel_res <- Unified_OAL_Selection(X[, idx, drop = FALSE], T_val, Y, method = "HSIC_ADL")
      res_HSIC_ADL[idx] <- sel_res$selected
      class_HSIC_ADL[idx] <- sel_res$var_class
    }
    if (length(idx) > 0) {
      sel_resY <- Unified_OAL_Selection(X[, idx, drop = FALSE], T_val, Y, method = "HSIC_Y")
      res_HSIC_Y[idx] <- sel_resY$selected
      class_HSIC_Y[idx] <- sel_resY$var_class
    }
    sel_oal <- Unified_OAL_Selection(X, T_val, Y, method = "Linear")
    res_OAL <- sel_oal$selected
    g <- tryCatch(GOAL(full_df, var.list = vars, covar = NULL), error = function(e) NULL)
    if (!is.null(g) && !is.null(g$selectedVar)) {
      pos <- match(g$selectedVar, vars)
      pos <- pos[!is.na(pos)]
      if (length(pos) > 0) res_GOAL[pos] <- 1
    }
    gd <- tryCatch(GOALDeR(full_df, var.list = vars, covar = NULL), error = function(e) NULL)
    if (!is.null(gd) && !is.null(gd$selectedVar)) {
      pos <- match(gd$selectedVar, vars)
      pos <- pos[!is.na(pos)]
      if (length(pos) > 0) res_GOALDER[pos] <- 1
    }
    ate_unadj <- estimate_ate_npcbps(X, T_val, Y, integer(0))
    ate_HSIC_ADL <- estimate_ate_npcbps(X, T_val, Y, which(res_HSIC_ADL == 1))
    ate_HSIC_Y <- estimate_ate_npcbps(X, T_val, Y, which(res_HSIC_Y == 1))
    ate_OAL <- estimate_ate_npcbps(X, T_val, Y, which(res_OAL == 1))
    ate_GOAL <- estimate_ate_npcbps(X, T_val, Y, which(res_GOAL == 1))
    ate_GOALDER <- estimate_ate_npcbps(X, T_val, Y, which(res_GOALDER == 1))
    ate_all <- estimate_ate_npcbps(X, T_val, Y, 1:p_curr)
    ate_x1_x5 <- estimate_ate_npcbps(X, T_val, Y, 1:5)
    ate_x1_x10 <- estimate_ate_npcbps(X, T_val, Y, 1:10)
    ate_x1_x15 <- estimate_ate_npcbps(X, T_val, Y, 1:15)
    out_list <- list(
      rep = rep_i,
      scenario = cfg$scenario,
      rho = cfg$rho,
      n = cfg$n,
      p = p_curr,
      res_HSIC_ADL = res_HSIC_ADL,
      class_HSIC_ADL = class_HSIC_ADL,
      res_HSIC_Y = res_HSIC_Y,
      class_HSIC_Y = class_HSIC_Y,
      res_OAL = res_OAL,
      res_GOAL = res_GOAL,
      res_GOALDER = res_GOALDER,
      ate_estimates = c(
        Unadjusted = ate_unadj,
        HSIC_ADL = ate_HSIC_ADL,
        HSIC_Y = ate_HSIC_Y,
        OAL = ate_OAL,
        GOAL = ate_GOAL,
        GOALDeR = ate_GOALDER,
        All_vars = ate_all,
        True_Confounders = ate_x1_x5,
        True_Conf_Outcome = ate_x1_x10,
        True_Conf_Outcome_IV = ate_x1_x15
      ),
      success = TRUE
    )
    return(out_list)
  }, error = function(e) {
    out_list <- list(
      rep = rep_i,
      scenario = cfg$scenario,
      rho = cfg$rho,
      n = cfg$n,
      p = p_curr,
      res_HSIC_ADL = rep(0, p_curr),
      class_HSIC_ADL = rep(NA, p_curr),
      res_HSIC_Y = rep(0, p_curr),
      class_HSIC_Y = rep(NA, p_curr),
      res_OAL = rep(0, p_curr),
      res_GOAL = rep(0, p_curr),
      res_GOALDER = rep(0, p_curr),
      ate_estimates = c(
        Unadjusted = NA,
        HSIC_ADL = NA,
        HSIC_Y = NA,
        OAL = NA,
        GOAL = NA,
        GOALDeR = NA,
        All_vars = NA,
        True_Confounders = NA,
        True_Conf_Outcome = NA,
        True_Conf_Outcome_IV = NA
      ),
      success = FALSE,
      error_msg = e$message
    )
    return(out_list)
  })
  return(res)
}
# Run all configured replications sequentially.
for (run_cfg_idx in seq_len(nrow(all_config))) {
  cfg <- all_config[run_cfg_idx, , drop = FALSE]
  folder_name <- sprintf("scen%d_n%d_p%d_rho%.1f", cfg$scenario, cfg$n, cfg$p, cfg$rho)
  save_dir <- file.path(out_root, folder_name)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  for (rep_i in seq_len(n_sim)) {
    rds_file <- file.path(save_dir, sprintf("rep_%03d.rds", rep_i))
    if (file.exists(rds_file)) {
      existing_res <- tryCatch(readRDS(rds_file), error = function(e) NULL)
      if (!is.null(existing_res) && isTRUE(existing_res$success)) {
        next
      }
    }
    res_out <- run_single_rep(cfg, rep_i)
    saveRDS(res_out, rds_file)
  }
}
library(data.table)
# Result root directory and simulation folders.
res_root <- out_root
folder_list <- list.dirs(res_root, recursive = FALSE, full.names = TRUE)
# Define the true ATE for each scenario.
get_true_ate <- function(scenario_val) {
  true_ate_map <- c("1" = 2.0, "2" = 2.0, "3" = 2.0, "4" = 2.0)
  val <- true_ate_map[as.character(scenario_val)]
  return(as.numeric(val))
}
# Calculate the average variable-selection frequency over a variable group.
calc_freq <- function(mat, idx) {
  if (length(idx) == 0) return(NA_real_)
  mean(rowMeans(mat[, idx, drop = FALSE]))
}
# Calculate classification accuracy for a true variable class.
get_type_acc <- function(mat_class, target_type, idx_set) {
  if (length(idx_set) == 0 || all(is.na(mat_class))) return(NA_real_)
  sub_mat <- mat_class[, idx_set, drop = FALSE]
  mean(sub_mat == target_type, na.rm = TRUE)
}
all_final <- list()
# Process each simulation folder and infer the variable dimension from the folder name.
for (fdr in folder_list) {
  fn <- basename(fdr)
  rds_files <- list.files(fdr, pattern = "^rep_.*\\.rds$", full.names = TRUE)
  n_rep <- length(rds_files)
  if (n_rep == 0) {
    next
  }
  scen_val <- as.integer(sub(".*scen([0-9]+).*", "\\1", fn))
  n_val <- as.integer(sub(".*n([0-9]+).*", "\\1", fn))
  p_total <- as.integer(sub(".*p([0-9]+).*", "\\1", fn))
  rho_val <- as.numeric(sub(".*rho([0-9.]+).*", "\\1", fn))
  important_vars_idx <- 1:min(10, p_total)
  noimportant_vars_idx <- if (p_total >= 11) 11:p_total else integer(0)
  IV_vars_idx <- if (p_total >= 15) 11:15 else if (p_total >= 11) 11:p_total else integer(0)
  useless_vars_idx <- if (p_total >= 16) 16:p_total else integer(0)
  confounder_idx <- intersect(1:5, seq_len(p_total))
  predictor_idx <- intersect(6:10, seq_len(p_total))
  iv_idx <- if (p_total >= 15) 11:15 else if (p_total >= 11) 11:p_total else integer(0)
  useless_idx <- if (p_total >= 16) 16:p_total else integer(0)
  mat_HSIC_ADL <- matrix(0, n_rep, p_total)
  mat_HSIC_Y <- matrix(0, n_rep, p_total)
  mat_OAL <- matrix(0, n_rep, p_total)
  mat_GOAL <- matrix(0, n_rep, p_total)
  mat_GOALDER <- matrix(0, n_rep, p_total)
  mat_class_HSIC_ADL <- matrix(NA, n_rep, p_total)
  mat_class_HSIC_Y <- matrix(NA, n_rep, p_total)
  ate_methods <- c(
    "Unadjusted",
    "HSIC_ADL",
    "HSIC_Y",
    "OAL",
    "GOAL",
    "GOALDeR",
    "All_vars",
    "True_Confounders",
    "True_Conf_Outcome",
    "True_Conf_Outcome_IV"
  )
  mat_ate <- matrix(NA, n_rep, length(ate_methods))
  colnames(mat_ate) <- ate_methods
  for (r in seq_along(rds_files)) {
    one <- readRDS(rds_files[r])
    if (!is.null(one$res_HSIC_ADL)) mat_HSIC_ADL[r, ] <- one$res_HSIC_ADL
    if (!is.null(one$class_HSIC_ADL)) mat_class_HSIC_ADL[r, ] <- one$class_HSIC_ADL
    if (!is.null(one$res_HSIC_Y)) mat_HSIC_Y[r, ] <- one$res_HSIC_Y
    if (!is.null(one$class_HSIC_Y)) mat_class_HSIC_Y[r, ] <- one$class_HSIC_Y
    if (!is.null(one$res_OAL)) mat_OAL[r, ] <- one$res_OAL
    if (!is.null(one$res_GOAL)) mat_GOAL[r, ] <- one$res_GOAL
    if (!is.null(one$res_GOALDER)) mat_GOALDER[r, ] <- one$res_GOALDER
    if (!is.null(one$ate_estimates)) {
      est_vec <- one$ate_estimates
      for (m in ate_methods) {
        if (m %in% names(est_vec)) mat_ate[r, m] <- est_vec[m]
      }
    }
  }
  true_ate <- get_true_ate(scen_val)
  methods_full <- c(
    "Unadjusted",
    "HSIC_ADL",
    "HSIC_Y",
    "OAL",
    "GOAL",
    "GOALDeR",
    "All_vars",
    "True_Confounders",
    "True_Conf_Outcome",
    "True_Conf_Outcome_IV"
  )
  bias_vec <- numeric(length(methods_full))
  sd_vec <- numeric(length(methods_full))
  mse_vec <- numeric(length(methods_full))
  for (i in seq_along(methods_full)) {
    m_ate <- methods_full[i]
    vals <- mat_ate[, m_ate]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      bias_vec[i] <- mean(vals) - true_ate
      mse_vec[i] <- mean((vals - true_ate)^2)
    } else {
      bias_vec[i] <- NA
      mse_vec[i] <- NA
    }
  }
  v_imp_freq <- c(
    NA,
    calc_freq(mat_HSIC_ADL, important_vars_idx),
    calc_freq(mat_HSIC_Y, important_vars_idx),
    calc_freq(mat_OAL, important_vars_idx),
    calc_freq(mat_GOAL, important_vars_idx),
    calc_freq(mat_GOALDER, important_vars_idx),
    rep(NA, 4)
  )
  v_noimp_freq <- c(
    NA,
    calc_freq(mat_HSIC_ADL, noimportant_vars_idx),
    calc_freq(mat_HSIC_Y, noimportant_vars_idx),
    calc_freq(mat_OAL, noimportant_vars_idx),
    calc_freq(mat_GOAL, noimportant_vars_idx),
    calc_freq(mat_GOALDER, noimportant_vars_idx),
    rep(NA, 4)
  )
  v_iv_freq <- c(
    NA,
    calc_freq(mat_HSIC_ADL, IV_vars_idx),
    calc_freq(mat_HSIC_Y, IV_vars_idx),
    calc_freq(mat_OAL, IV_vars_idx),
    calc_freq(mat_GOAL, IV_vars_idx),
    calc_freq(mat_GOALDER, IV_vars_idx),
    rep(NA, 4)
  )
  v_useless_freq <- c(
    NA,
    calc_freq(mat_HSIC_ADL, useless_vars_idx),
    calc_freq(mat_HSIC_Y, useless_vars_idx),
    calc_freq(mat_OAL, useless_vars_idx),
    calc_freq(mat_GOAL, useless_vars_idx),
    calc_freq(mat_GOALDER, useless_vars_idx),
    rep(NA, 4)
  )
  acc_confounder <- c(
    NA,
    get_type_acc(mat_class_HSIC_ADL, 1, confounder_idx),
    get_type_acc(mat_class_HSIC_Y, 1, confounder_idx),
    rep(NA, 7)
  )
  acc_predictor <- c(
    NA,
    get_type_acc(mat_class_HSIC_ADL, 2, predictor_idx),
    get_type_acc(mat_class_HSIC_Y, 2, predictor_idx),
    rep(NA, 7)
  )
  acc_iv <- c(
    NA,
    get_type_acc(mat_class_HSIC_ADL, 3, iv_idx),
    get_type_acc(mat_class_HSIC_Y, 3, iv_idx),
    rep(NA, 7)
  )
  acc_useless <- c(
    NA,
    get_type_acc(mat_class_HSIC_ADL, 0, useless_idx),
    get_type_acc(mat_class_HSIC_Y, 0, useless_idx),
    rep(NA, 7)
  )
  df_cum <- data.frame(
    scenario = scen_val,
    n = n_val,
    p = p_total,
    rho = rho_val,
    method = methods_full,
    important_selection_frequency = v_imp_freq,
    nonimportant_selection_frequency = v_noimp_freq,
    IV_selection_frequency = v_iv_freq,
    useless_selection_frequency = v_useless_freq,
    confounder_classification_accuracy = acc_confounder,
    predictor_classification_accuracy = acc_predictor,
    IV_classification_accuracy = acc_iv,
    useless_classification_accuracy = acc_useless,
    ATE_Bias = bias_vec,
    ATE_MSE = mse_vec
  )
  num_cols <- sapply(df_cum, is.numeric)
  cols_to_round <- setdiff(names(df_cum)[num_cols], c("scenario", "n", "p"))
  df_cum[, cols_to_round] <- round(df_cum[, cols_to_round], 4)
  all_final[[fn]] <- df_cum
big_df <- rbindlist(all_final, use.names = TRUE, fill = TRUE)
write.csv(
  big_df,
  file.path(ROOT_DIR, "results_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
print(big_df)
}