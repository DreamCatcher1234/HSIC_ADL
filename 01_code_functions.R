library(MASS)
library(glmnet)
library(ranger)
library(independenceWeights)
library(lqa)
library(GOAL)
library(cdcsis)
library(Matrix)
library(extracat)
library(mclust)
library(CBPS)
library(here)
# Fit a random forest and return predictions for a test set.
fit_rf_pred <- function(x_train, y_train, x_test, ntree = 300) {
  x_train <- as.data.frame(x_train)
  x_test <- as.data.frame(x_test)
  fit <- ranger(y = y_train, x = x_train, num.trees = ntree)
  pred <- predict(fit, data = x_test)$predictions
  return(pred)
}
# Compute the maximum absolute pairwise correlation after removing constant columns.
# The function name is retained for compatibility with the existing simulation code.
get_design_condition_number <- function(X_mat) {
  X_mat <- as.matrix(X_mat)
  if (ncol(X_mat) <= 1) return(0)
  var_vals <- apply(X_mat, 2, var)
  valid_cols <- which(var_vals > 1e-8)
  if (length(valid_cols) <= 1) return(0)
  X_sub <- X_mat[, valid_cols, drop = FALSE]
  if (ncol(X_sub) <= 1) return(0)
  cor_mat <- abs(cor(X_sub))
  diag(cor_mat) <- 0
  max_cor <- max(cor_mat, na.rm = TRUE)
  return(max_cor)
}
gaussian_grammat_rcpp <- function(x, bandwidth, n = nrow(x), d = ncol(x)) {
  x <- as.matrix(x)
  sq_dist <- as.matrix(dist(x))^2
  K <- exp(-sq_dist / (2.0 * bandwidth^2))
  return(K)
}

discrete_grammat_rcpp <- function(x, n = nrow(x), d = ncol(x)) {
  x <- as.matrix(x)
  K <- as.matrix(dist(x, method = "maximum")) == 0
  mode(K) <- "logical"
  return(K)
}

median_bandwidth_rcpp <- function(x, n = nrow(x), d = ncol(x)) {
  x <- as.matrix(x)
  len <- ifelse(n > 1000, 1000, n)
  x_sub <- x[1:len, , drop = FALSE]
  
  dists <- as.numeric(dist(x_sub)^2)
  lentot <- length(dists)
  middle <- lentot %/% 2
  
  bw_sq <- sort(dists, partial = middle + 1)[middle + 1]
  bandwidth <- sqrt(bw_sq * 0.5)
  return(bandwidth)
}
# Compute the weighted distance-correlation criterion used by the custom GOALDeR implementation.
DWDC_function <- function(DataM, varlist, trt.var, wgt, beta) {
  diff_vec <- rep(NA, length(beta))
  names(diff_vec) <- varlist
  for (jj in 1:length(varlist)) {
    diff_vec[jj] <- abs(wdcor(x = DataM[, trt.var], y = DataM[, varlist[jj]], w = DataM[, wgt]))
  }
  wdiff_vec <- diff_vec * abs(beta)
  wAMD <- c(sum(wdiff_vec))
  ret <- list(diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD)
  return(ret)
}
# Custom GOALDeR implementation used as a comparison method.
GOALDeR <- function(data, var.list, covar, Trt = "Trt", out = "Y",
                    lambda_vec = c(-10, -5, -2, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
                    gamma_convergence = 2) {
  Data <- data
  n <- dim(Data)[1]
  temp.mean <- colMeans(Data[, var.list])
  Temp.mean <- matrix(temp.mean, ncol = length(var.list), nrow = nrow(Data), byrow = TRUE)
  Data[, var.list] <- Data[, var.list] - Temp.mean
  temp.sd <- apply(Data[var.list], FUN = sd, MARGIN = 2)
  Temp.sd <- matrix(temp.sd, ncol = length(var.list), nrow = nrow(Data), byrow = TRUE)
  Data[, var.list] <- Data[, var.list] / Temp.sd
  rm(list = c("temp.mean", "Temp.mean", "temp.sd", "Temp.sd"))
  betaXY_cor <- NA
  for (i in 1:length(var.list)) {
    var.name <- var.list[i]
    betaXY_cor[i] <- as.numeric(cdcor(Data[, var.name], Data[, out], Data[, Trt])$statistic)
  }
  betaXY <- betaXY_cor / max(betaXY_cor)
  names(betaXY) <- var.list
  lambda_vec <- lambda_vec
  names(lambda_vec) <- as.character(lambda_vec)
  gamma_convergence_factor <- gamma_convergence
  gamma_vals <- 2 * (gamma_convergence_factor - lambda_vec + 1)
  names(gamma_vals) <- names(lambda_vec)
  wAMD_vec <- rep(NA, length(lambda_vec))
  DR_estimated <- vector(mode = "list", length(lambda_vec))
  DR_pseudo_outcome <- matrix(NA, nrow = dim(Data)[1], ncol = length(lambda_vec))
  names(DR_estimated) <- names(wAMD_vec) <- names(lambda_vec)
  colnames(DR_pseudo_outcome) <- names(lambda_vec)
  Dcow.var <- list()
  w.full.form <- formula(paste(Trt, "~", paste(c(covar, var.list), collapse = "+")))
  for(lil in names(lambda_vec)) {
    v <- which(names(lambda_vec) == lil)
    il <- lambda_vec[lil]
    ig <- gamma_vals[lil]
    oal_pen <- adaptive.lasso(lambda = n^(il), al.weights = abs(betaXY)^(-ig))
    logit_oal <- lqa.formula(w.full.form, data = Data, penalty = oal_pen, family = gaussian())
    coeff_XA <- coef(logit_oal)[var.list]
    Dcow.var[[v]] <- names(coeff_XA)[which(round(coeff_XA, 5) != 0)]
    names(Dcow.var)[[v]] <- paste("lambda", lil, sep = "_")
    if (length(Dcow.var[[v]]) != 0) {
      w.model <- formula(paste(Trt, "~", paste(Dcow.var[[v]], collapse = "+")))
      Dcow_fit <- independence_weights(Data[, Trt], Data[, Dcow.var[[v]]])
      Data[, paste("w", lil, sep = "")] <- Dcow_fit$weights
    } else {
      Data[, paste("w", lil, sep = "")] <- 1
    }
    wAMD_vec[lil] <- DWDC_function(DataM = Data, varlist = names(betaXY), trt.var = Trt,
                                   wgt = paste("w", lil, sep = ""), beta = (betaXY_cor)^2)$wAMD
  }
  Svar <- Dcow.var[[which.min(wAMD_vec)]]
  lambda <- names(wAMD_vec)[which.min(wAMD_vec)]
  w_lil <- names(wAMD_vec)[which.min(wAMD_vec)]
  fw <- Data[, paste("w", w_lil, sep = "")]
  GOALDeR_results <- list(
    selectedVar = Svar,
    lambda = lambda,
    fw = fw
  )
  return(GOALDeR_results)
}
# Scale continuous columns using training-set moments while leaving discrete-like columns unchanged.
smart_scale <- function(X_train, X_test = NULL) {
  X_train <- as.matrix(X_train)
  p <- ncol(X_train)
  mean_vec <- rep(NA_real_, p)
  sd_vec <- rep(NA_real_, p)
  X_train_scaled <- X_train
  for (j in 1:p) {
    if (length(unique(X_train[, j])) <= 10 || all(X_train[, j] == floor(X_train[, j]))) {
      mean_vec[j] <- NA
      sd_vec[j] <- NA
    } else {
      mean_vec[j] <- mean(X_train[, j])
      sd_vec[j] <- sd(X_train[, j])
      if (sd_vec[j] > 1e-8) {
        X_train_scaled[, j] <- (X_train[, j] - mean_vec[j]) / sd_vec[j]
      }
    }
  }
  if (is.null(X_test)) {
    return(list(X_scaled = X_train_scaled, mean = mean_vec, sd = sd_vec))
  } else {
    X_test <- as.matrix(X_test)
    X_test_scaled <- X_test
    for (j in 1:p) {
      if (!is.na(mean_vec[j]) && !is.na(sd_vec[j]) && sd_vec[j] > 1e-8) {
        X_test_scaled[, j] <- (X_test[, j] - mean_vec[j]) / sd_vec[j]
      }
    }
    return(list(X_train = X_train_scaled, X_test = X_test_scaled))
  }
}
# Two-stage variable screening:
# Step 1 uses Elastic Net double selection.
# Step 2 uses cross-fitted residuals and HSIC tests on the candidate set.
get_keep_idx <- function(X, T_val, Y, alpha = 0.05, K = 5, cn_threshold = 0.5, enet_alpha = 0.5) {
  n <- nrow(X)
  p <- ncol(X)
  X_raw_mat <- as.matrix(X)
  T_raw_vec <- as.vector(T_val)
  Y_raw <- as.vector(Y)
  fit_T <- cv.glmnet(X_raw_mat, T_raw_vec, alpha = enet_alpha, family = "gaussian")
  idx_T <- which(as.vector(coef(fit_T, s = "lambda.1se"))[-1] != 0)
  XT_mat <- cbind(T = T_raw_vec, X_raw_mat)
  fit_Y <- cv.glmnet(XT_mat, Y_raw, alpha = enet_alpha, family = "gaussian")
  coef_X_in_Y <- as.vector(coef(fit_Y, s = "lambda.1se"))[-(1:2)]
  idx_Y <- which(coef_X_in_Y != 0)
  idx_step1 <- sort(union(idx_T, idx_Y))
  if (length(idx_step1) == 0) return(integer(0))
  set.seed(123)
  folds <- sample(rep(1:K, length.out = n))
  k_len <- length(idx_step1)
  p_vals_T <- numeric(k_len)
  p_vals_Y <- numeric(k_len)
  for (idx_pos in 1:k_len) {
    j <- idx_step1[idx_pos]
    Xj_raw_vec <- X_raw_mat[, j]
    other_idx <- setdiff(idx_step1, j)
    res_T <- numeric(n)
    res_Y <- numeric(n)
    if (length(other_idx) == 0) {
      res_T <- scale(T_raw_vec, scale = FALSE)
      res_Y <- scale(Y_raw, scale = FALSE)
    } else {
      otherX_full <- X_raw_mat[, other_idx, drop = FALSE]
      cn_otherX <- get_design_condition_number(otherX_full)
      use_rf <- (cn_otherX >= cn_threshold)
      for (k in 1:K) {
        train_idx <- which(folds != k)
        test_idx <- which(folds == k)
        scaled_otherX <- smart_scale(otherX_full[train_idx, , drop = FALSE],
                                     otherX_full[test_idx, , drop = FALSE])
        otherX_train <- scaled_otherX$X_train
        otherX_test <- scaled_otherX$X_test
        T_mean <- mean(T_raw_vec[train_idx])
        T_sd <- sd(T_raw_vec[train_idx])
        T_train <- (T_raw_vec[train_idx] - T_mean) / T_sd
        T_test <- (T_raw_vec[test_idx] - T_mean) / T_sd
        TY_train <- T_raw_vec[train_idx]
        TY_test <- T_raw_vec[test_idx]
        Y_train <- Y_raw[train_idx]
        Y_test <- Y_raw[test_idx]
        XT_train <- cbind(T = T_train, otherX_train)
        XT_test <- cbind(T = T_test, otherX_test)
        if (use_rf) {
          pred_T <- fit_rf_pred(otherX_train, TY_train, otherX_test)
          res_T[test_idx] <- TY_test - pred_T
        } else {
          fit_T_lin <- cv.glmnet(x = as.matrix(otherX_train), y = TY_train, alpha = 0.5)
          res_T[test_idx] <- TY_test - as.vector(predict(fit_T_lin, newx = as.matrix(otherX_test), s = "lambda.min"))
        }
        if (use_rf) {
          pred_Y <- fit_rf_pred(XT_train, Y_train, XT_test)
          res_Y[test_idx] <- Y_test - pred_Y
        } else {
          fit_Y_lin <- cv.glmnet(x = as.matrix(XT_train), y = Y_train, alpha = 0.5)
          res_Y[test_idx] <- Y_test - as.vector(predict(fit_Y_lin, newx = as.matrix(XT_test), s = "lambda.min"))
        }
      }
    }
    test_T <- dhsic.test(X = as.matrix(Xj_raw_vec), Y = as.matrix(res_T), alpha = alpha,
                         method = "eigenvalue", kernel = "gaussian", matrix.input = FALSE)
    p_vals_T[idx_pos] <- test_T$p.value
    test_Y <- dhsic.test(X = as.matrix(Xj_raw_vec), Y = as.matrix(res_Y), alpha = alpha,
                         method = "eigenvalue", kernel = "gaussian", matrix.input = FALSE)
    p_vals_Y[idx_pos] <- test_Y$p.value
  }
  p_adj_T <- p.adjust(p_vals_T, method = "BH")
  p_adj_Y <- p.adjust(p_vals_Y, method = "BH")
  keep_logic <- (p_adj_T <= alpha) | (p_adj_Y <= alpha)
  keep_idx <- idx_step1[keep_logic]
  return(keep_idx)
}
# Generate the simulation data under the scenarios defined in the uploaded source code.
generate_data_flexible <- function(n = 200, p = 20, d = 0, rho = 0.5, scenario = 1) {
  Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
  X_raw <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  n_disc <- floor(p * d)
  if (n_disc > 0) {
    disc_idx <- c()
    if (n_disc >= 1) disc_idx <- c(disc_idx, 1)
    if (n_disc >= 2) disc_idx <- c(disc_idx, 6)
    if (n_disc >= 3) disc_idx <- c(disc_idx, 11)
    remaining <- n_disc - 3
    if (remaining > 0 && p >= 16) {
      num_block4 <- min(remaining, p - 15)
      disc_idx <- c(disc_idx, 16:(16 + num_block4 - 1))
    }
    for (j in disc_idx) {
      X_raw[, j] <- ifelse(X_raw[, j] > 0, 1, 0)
    }
  }
  X <- as.data.frame(X_raw)
  colnames(X) <- paste0("X", 1:p)
  L_1_5 <- rowSums(X[, 1:5])
  L_6_10 <- rowSums(X[, 6:10])
  L_11_15 <- rowSums(X[, 11:15])
  W1 <- exp(X$X1 / 2)
  W2 <- X$X2 / (1 + exp(X$X1))
  W3 <- (X$X1 * X$X3 / 25 + 0.6)^3
  W4 <- (X$X2 + X$X4 + 20)^2
  W5 <- X$X5
  NL_1_5 <- W1 + W2 + W3 + W4 + W5
  W6 <- exp(X$X6 / 2)
  W7 <- X$X7 / (1 + exp(X$X6))
  W8 <- (X$X6 * X$X8 / 25 + 0.6)^3
  W9 <- (X$X7 + X$X9 + 20)^2
  W10 <- X$X10
  NL_6_10 <- W6 + W7 + W8 + W9 + W10
  W11 <- exp(X$X11 / 2)
  W12 <- X$X12 / (1 + exp(X$X11))
  W13 <- (X$X11 * X$X13 / 25 + 0.6)^3
  W14 <- (X$X12 + X$X14 + 20)^2
  W15 <- X$X15
  NL_11_15 <- W11 + W12 + W13 + W14 + W15
  scen_str <- as.character(scenario)
  switch(scen_str,
         "1" = {
           T_val <- L_1_5 + L_11_15 + rnorm(n, 0, 1)
           Y_val <- 2 * T_val + 0.6 * L_1_5 + 0.6 * L_6_10 + rnorm(n, 0, 1)
         },
         "2" = {
           T_val <- 0.4 * L_1_5 + L_11_15 + rnorm(n, 0, 1)
           Y_val <- 2 * T_val + 0.2 * L_1_5 + 0.6 * L_6_10 + rnorm(n, 0, 1)
         },
         "3" = {
           T_val <- L_1_5 + 1.8 * L_11_15 + rnorm(n, 0, 1)
           Y_val <- 2 * T_val + 0.6 * L_1_5 + 0.6 * L_6_10 + rnorm(n, 0, 1)
         },
         "4" = {
           T_val <- NL_1_5 + NL_11_15 + rnorm(n, 0, 1)
           Y_val <- 2 * T_val + 0.6 * NL_1_5 + 0.6 * NL_6_10 + rnorm(n, 0, 1)
         },
         stop(paste0("Unsupported scenario: ", scen_str))
  )
  return(list(X = X, T = T_val, Y = Y_val))
}
# Unified outcome-adaptive selection interface.
# method = "HSIC_ADL" uses the treatment-to-outcome dependence ratio; method = "HSIC_Y" uses the outcome-based ratio.
Unified_OAL_Selection <- function(X, T_val, Y, K = 5, cn_threshold = 0.5, epsilon = 0.1, window_size = 5, method = "HSIC_ADL") {
  n <- nrow(X)
  p <- ncol(X)
  var_names <- colnames(X)
  wts_uni <- rep(1, n)
  ratio_base <- numeric(p)
  names(ratio_base) <- var_names
  if (method %in% c("HSIC_ADL", "HSIC_Y")) {
    X_raw_mat <- as.matrix(X)
    T_raw_vec <- as.vector(T_val)
    Y_raw <- as.vector(Y)
    set.seed(123)
    folds <- sample(rep(1:K, length.out = n))
    betaXY_cor <- numeric(p)
    betaXT_dcor <- numeric(p)
    for (j in 1:p) {
      Xj_raw_vec <- X_raw_mat[, j]
      otherX_full <- X_raw_mat[, -j, drop = FALSE]
      cn_otherX <- get_design_condition_number(otherX_full)
      use_rf <- (cn_otherX >= cn_threshold)
      res_T_full <- numeric(n)
      res_Y_full <- numeric(n)
      for (k in 1:K) {
        train_idx <- which(folds != k)
        test_idx <- which(folds == k)
        scaled_otherX <- smart_scale(otherX_full[train_idx, , drop = FALSE], otherX_full[test_idx, , drop = FALSE])
        otherX_train <- scaled_otherX$X_train
        otherX_test <- scaled_otherX$X_test
        T_mean <- mean(T_raw_vec[train_idx])
        T_sd <- sd(T_raw_vec[train_idx])
        T_train <- (T_raw_vec[train_idx] - T_mean) / T_sd
        T_test <- (T_raw_vec[test_idx] - T_mean) / T_sd
        TY_train <- T_raw_vec[train_idx]
        TY_test <- T_raw_vec[test_idx]
        Y_train <- Y_raw[train_idx]
        Y_test <- Y_raw[test_idx]
        XT_train <- cbind(T = T_train, otherX_train)
        XT_test <- cbind(T = T_test, otherX_test)
        if (use_rf) {
          pred_T <- fit_rf_pred(otherX_train, TY_train, otherX_test)
          res_T_full[test_idx] <- TY_test - pred_T
          pred_Y <- fit_rf_pred(XT_train, Y_train, XT_test)
          res_Y_full[test_idx] <- Y_test - pred_Y
        } else {
          fit_T <- cv.glmnet(x = as.matrix(otherX_train), y = TY_train, alpha = 0)
          res_T_full[test_idx] <- TY_test - as.vector(predict(fit_T, newx = as.matrix(otherX_test), s = "lambda.min"))
          fit_Y <- cv.glmnet(x = as.matrix(XT_train), y = Y_raw[train_idx], alpha = 0)
          res_Y_full[test_idx] <- Y_test - as.vector(predict(fit_Y, newx = as.matrix(XT_test), s = "lambda.min"))
        }
      }
      betaXY_cor[j] <- calculate_weighted_hsic(w = wts_uni, X = as.matrix(Xj_raw_vec), A = res_Y_full)
      betaXT_dcor[j] <- calculate_weighted_hsic(w = wts_uni, X = as.matrix(Xj_raw_vec), A = res_T_full)
    }
    betaXY_cor <- as.numeric(betaXY_cor)
    betaXT_dcor <- as.numeric(betaXT_dcor)
    if (method == "HSIC_ADL") {
      ratio_base <- betaXT_dcor / betaXY_cor
    } else if (method == "HSIC_Y") {
      ratio_base <- 1 / betaXY_cor
    }
    names(ratio_base) <- var_names
  } else {
    X_mat_global <- as.matrix(X)
    T_scaled_global <- as.vector(scale(T_val))
    lm_fit <- lm(Y ~ T_scaled_global + X_mat_global)
    coef_all <- coef(lm_fit)
    coef_X <- coef_all[-(1:2)]
    coef_X[is.na(coef_X)] <- 1e-8
    metric_val <- 1 / abs(coef_X)
    ratio_base <- metric_val
    names(ratio_base) <- var_names
  }
  alpha_grid <- seq(0.1, 1000, by = 1)
  trunc_ratios <- numeric(length(alpha_grid))
  for (k in seq_along(alpha_grid)) {
    a <- alpha_grid[k]
    W_j <- ratio_base^a
    Z_j <- pmax(W_j, epsilon)
    trunc_ratios[k] <- mean(Z_j <= epsilon)
  }
  n_grid <- length(alpha_grid)
  if (window_size >= n_grid) window_size <- 3
  n_windows <- n_grid - window_size + 1
  win_vars <- numeric(n_windows)
  win_means <- numeric(n_windows)
  for (i in 1:n_windows) {
    sub_ratios <- trunc_ratios[i:(i + window_size - 1)]
    win_vars[i] <- var(sub_ratios)
    win_means[i] <- mean(sub_ratios)
  }
  valid_win_idx <- which(win_means > 0 & win_means <= 0.5)
  if (length(valid_win_idx) > 0) {
    best_win_sub_idx <- valid_win_idx[which.min(win_vars[valid_win_idx])]
    opt_alpha_idx <- best_win_sub_idx + floor(window_size / 2)
  } else {
    opt_alpha_idx <- which(trunc_ratios > 0)[1]
    if (is.na(opt_alpha_idx)) opt_alpha_idx <- 1
  }
  opt_alpha <- alpha_grid[opt_alpha_idx]
  W_opt <- ratio_base^opt_alpha
  gmm_fit <- tryCatch({
    Mclust(W_opt, G = 3, modelNames = "V", verbose = FALSE)
  }, error = function(e) NULL)
  if (!is.null(gmm_fit) && !is.null(gmm_fit$classification)) {
    clusters <- gmm_fit$classification
    means <- gmm_fit$parameters$mean
    sorted_cluster_idx <- order(means)
    vp_cluster <- sorted_cluster_idx[1]
    vc_cluster <- sorted_cluster_idx[2]
    viv_cluster <- sorted_cluster_idx[3]
    VP_idx <- which(clusters == vp_cluster)
    VC_idx <- which(clusters == vc_cluster)
    VIV_idx <- which(clusters == viv_cluster)
  } else {
    q_vals <- quantile(W_opt, probs = c(1 / 3, 2 / 3))
    VP_idx <- which(W_opt <= q_vals[1])
    VC_idx <- which(W_opt > q_vals[1] & W_opt <= q_vals[2])
    VIV_idx <- which(W_opt > q_vals[2])
  }
  var_class <- numeric(p)
  names(var_class) <- var_names
  var_class[VC_idx] <- 1
  var_class[VP_idx] <- 2
  var_class[VIV_idx] <- 3
  selected_idx <- sort(c(VC_idx, VP_idx))

  selected <- numeric(p)
  selected[selected_idx] <- 1
  names(selected) <- var_names
  res <- list(
    selected = selected,
    opt_alpha = opt_alpha,
    trunc_ratios = trunc_ratios,
    alpha_grid = alpha_grid,
    W_values = W_opt,
    var_class = var_class,
    VP_idx = VP_idx,
    VIV_idx = VIV_idx,
    selected_idx = selected_idx,
    selected_names = var_names[selected_idx]
  )
  return(res)
}
# Compute an RBF Gram matrix, using an identity-style kernel for discrete-like inputs.
gaussian_kernel <- function(X, sigma = NULL) {
  X <- as.matrix(X)
  unique_vals <- length(unique(as.vector(X)))
  if (unique_vals <= 10) {
    K <- as.matrix(dist(X, method = "manhattan")) == 0
    mode(K) <- "numeric"
    return(K)
  }
  if (is.null(sigma)) {
    dists <- as.matrix(dist(X))
    sigma <- median(dists[upper.tri(dists)]) / 2
    if (is.na(sigma) || sigma < 1e-8) sigma <- 1.0
  }
  K <- exp(-as.matrix(dist(X))^2 / (2 * sigma^2 + 1e-10))
  return(K)
}
# Compute the weighted HSIC objective for a pair of variables.
calculate_weighted_hsic <- function(w, X, A, sigma_X = NULL, sigma_A = NULL) {
  X <- as.matrix(X)
  A <- as.matrix(A)
  n <- length(w)
  if (ncol(X) == 0) return(1e9)
  K <- gaussian_kernel(X, sigma = sigma_X)
  L <- gaussian_kernel(as.matrix(A), sigma = sigma_A)
  H <- Diagonal(n) - matrix(1 / n, nrow = n, ncol = n)
  K_tilde <- H %*% K %*% H
  L_tilde <- H %*% L %*% H
  Q <- as.matrix(K_tilde * L_tilde)
  hsic_value <- as.numeric(t(w) %*% Q %*% w) / (n^2)
  return(hsic_value)
}
# Estimate the ATE using nonparametric CBPS weights followed by a weighted linear regression.
estimate_ate_npcbps <- function(X, T_val, Y, selected_idx) {
  if (length(selected_idx) == 0) {
    fit_unadj <- lm(Y ~ T_val)
    return(as.numeric(coef(fit_unadj)["T_val"]))
  }
  X_sub <- as.matrix(X[, selected_idx, drop = FALSE])
  df_cbps <- data.frame(T_val = T_val, X_sub)
  fit_weights <- tryCatch({
    cbps_mod <- CBPS::npCBPS(T_val ~ ., data = df_cbps, print.level = 0)
    cbps_mod$weights
  }, error = function(e) {
    rep(1, length(Y))
  })
  fit_lm <- lm(Y ~ T_val, weights = fit_weights)
  return(as.numeric(coef(fit_lm)["T_val"]))
}
