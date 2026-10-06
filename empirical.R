library(data.table)
library(dplyr)
library(ggplot2)
library(showtext)
ROOT_DIR <- here()
setwd(ROOT_DIR)
DATA_FILE <- "depmap_l1000_causal_data.csv"
SOURCE_FILE <- "01_code_functions.R"
setwd(ROOT_DIR)
source(SOURCE_FILE)
source("dhsic.test.R")
n_boot <- 1000
base_seed <- 123
out_dir <- file.path(ROOT_DIR, "empirical_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
df_raw <- fread(DATA_FILE, data.table = FALSE)
y_col_name <- colnames(df_raw)[1]
t_col_name <- colnames(df_raw)[2]
x_col_names <- colnames(df_raw)[3:ncol(df_raw)]
run_single_boot <- function(b_idx, df_data, x_cols, t_col, y_col) {
  curr_seed <- base_seed + b_idx
  set.seed(curr_seed)
  boot_sample_idx <- sample.int(nrow(df_data), size = nrow(df_data), replace = TRUE)
  data_b <- df_data[boot_sample_idx, ]
  Y <- data_b[[y_col]]
  T_val <- data_b[[t_col]]
  X <- as.matrix(data_b[, x_cols, drop = FALSE])
  p_curr <- length(x_cols)
  res_HSIC_ADL <- numeric(p_curr)
  res_HSIC_Y <- numeric(p_curr)
  res_OAL <- numeric(p_curr)
  names(res_HSIC_ADL) <- names(res_HSIC_Y) <- names(res_OAL) <- x_cols
  res <- tryCatch({
    idx <- get_keep_idx(X, T_val, Y)
    if (length(idx) > 0) {
      sel_res <- Unified_OAL_Selection(X[, idx, drop = FALSE], T_val, Y, method = "HSIC_ADL")
      res_HSIC_ADL[idx] <- sel_res$selected
      sel_resY <- Unified_OAL_Selection(X[, idx, drop = FALSE], T_val, Y, method = "HSIC_Y")
      res_HSIC_Y[idx] <- sel_resY$selected
    }
    sel_oal <- Unified_OAL_Selection(X, T_val, Y, method = "Linear")
    res_OAL <- sel_oal$selected
    ate_unadj <- estimate_ate_npcbps(X, T_val, Y, integer(0))
    ate_HSIC_ADL <- estimate_ate_npcbps(X, T_val, Y, which(res_HSIC_ADL == 1))
    ate_HSIC_Y <- estimate_ate_npcbps(X, T_val, Y, which(res_HSIC_Y == 1))
    ate_OAL <- estimate_ate_npcbps(X, T_val, Y, which(res_OAL == 1))
    list(
      boot_rep = b_idx,
      seed = curr_seed,
      var_names = x_cols,
      res_HSIC_ADL = res_HSIC_ADL,
      res_HSIC_Y = res_HSIC_Y,
      res_OAL = res_OAL,
      ate_estimates = c(
        Unadjusted = ate_unadj,
        HSIC_ADL = ate_HSIC_ADL,
        HSIC_Y = ate_HSIC_Y,
        OAL = ate_OAL
      ),
      success = TRUE
    )
  }, error = function(e) {
    list(
      boot_rep = b_idx,
      seed = curr_seed,
      var_names = x_cols,
      res_HSIC_ADL = rep(0, p_curr),
      res_HSIC_Y = rep(0, p_curr),
      res_OAL = rep(0, p_curr),
      ate_estimates = rep(NA_real_, 4),
      success = FALSE,
      error_msg = e$message
    )
  })
  return(res)
}
# Run bootstrap loop
for (b_i in seq_len(n_boot)) {
  curr_seed <- base_seed + b_i
  rds_file <- file.path(out_dir, sprintf("rep_%d.rds", curr_seed))
  if (file.exists(rds_file)) next
  res_out <- run_single_boot(b_i, df_raw, x_col_names, t_col_name, y_col_name)
  saveRDS(res_out, rds_file)
}
gc()
# Aggregate bootstrap results
rds_files <- list.files(out_dir, pattern = "^rep_.*\\.rds$", full.names = TRUE)
ate_list <- list()
HSIC_ADL_mat_list <- list()
var_names <- NULL
for (f in rds_files) {
  res <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.null(res) && isTRUE(res$success)) {
    ate_list[[length(ate_list) + 1]] <- res$ate_estimates
    if (is.null(var_names)) var_names <- res$var_names
    HSIC_ADL_mat_list[[length(HSIC_ADL_mat_list) + 1]] <- res$res_HSIC_ADL
  }
}
# Print mean ATE estimates
ate_matrix <- do.call(rbind, ate_list)
mean_ate <- colMeans(ate_matrix, na.rm = TRUE)
print(mean_ate)
# Calculate selection frequency for plotting
mat_HSIC_ADL <- do.call(rbind, HSIC_ADL_mat_list)
HSIC_ADL_freq_df <- data.frame(
  Gene = var_names,
  Selection_Frequency = colMeans(mat_HSIC_ADL, na.rm = TRUE),
  stringsAsFactors = FALSE
)
# Initialize plotting configurations
showtext_auto()
showtext_opts(dpi = 300)
tryCatch({ font_add("YaHei", "msyh.ttc") }, error = function(e) {})
FONT_GENE_LABEL_FIG1 <- 3.0
FONT_PCT_LABEL_FIG1 <- 1.2
FONT_AXIS_TITLE_FIG1 <- 5
FONT_AXIS_TICK_X <- 4
rich_spectrum_palette <- c(
  "#1a237e", "#0288d1", "#009688", "#8bc34a", "#ffeb3b",
  "#fbc02d", "#f57c00", "#d32f2f", "#880e4f"
)
# Figure 1: Selected core variables with frequency >= 0.60
top_selected_df <- HSIC_ADL_freq_df %>%
  filter(Selection_Frequency >= 0.60) %>%
  arrange(desc(Selection_Frequency)) %>%
  mutate(
    Gene = factor(Gene, levels = rev(Gene)),
    Pct_Label = sprintf("%.1f%%", Selection_Frequency * 100)
  )
n_top_genes <- nrow(top_selected_df)
p_fig1_bar <- ggplot(top_selected_df, aes(x = Gene, y = Selection_Frequency)) +
  geom_col(aes(fill = Selection_Frequency), width = 0.78) +
  geom_text(
    aes(label = Pct_Label),
    color = "#222222",
    hjust = -0.10,
    size = FONT_PCT_LABEL_FIG1,
    fontface = "plain"
  ) +
  coord_flip(clip = "off") +
  scale_fill_gradientn(
    colors = rich_spectrum_palette,
    limits = c(0.6, 1.0),
    labels = scales::percent,
    name = "Selection\nProb."
  ) +
  scale_y_continuous(
    limits = c(0, 1.15),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent,
    expand = c(0, 0)
  ) +
  labs(
    x = "Selected Landmark Genes",
    y = "Bootstrap Inclusion Probability"
  ) +
  theme_minimal(base_size = 6) +
  theme(
    axis.text.y = element_text(size = FONT_GENE_LABEL_FIG1, color = "#222222", margin = margin(r = 0.4)),
    axis.text.x = element_text(size = FONT_AXIS_TICK_X, color = "#222222"),
    axis.title = element_text(size = FONT_AXIS_TITLE_FIG1, color = "#111111", face = "bold"),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "#cccccc", linewidth = 0.2),
    axis.line.y = element_line(color = "#cccccc", linewidth = 0.2),
    legend.position = "right",
    legend.title = element_text(size = 4.8, face = "bold"),
    legend.text = element_text(size = 4.2),
    legend.key.width = unit(0.15, "cm"),
    legend.key.height = unit(0.35, "cm"),
    plot.margin = margin(t = 3, r = 12, b = 3, l = 3)
  )
single_col_height <- max(1.4, n_top_genes * 0.055)
ggsave("Figure1_Selected_Core_Variables_SingleColumn.png", plot = p_fig1_bar, width = 3.35, height = single_col_height, dpi = 300)
# Figure 2: All variables plot
all_genes_df <- HSIC_ADL_freq_df %>%
  arrange(desc(Selection_Frequency)) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene)))
p_fig2_appendix <- ggplot(all_genes_df, aes(x = Gene, y = Selection_Frequency)) +
  geom_col(aes(fill = Selection_Frequency), width = 0.85) +
  geom_hline(yintercept = 0.60, linetype = "dashed", color = "#d32f2f", linewidth = 0.2) +
  coord_flip() +
  scale_fill_gradientn(
    colors = rich_spectrum_palette,
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent,
    name = "Bootstrap Inclusion Rate"
  ) +
  scale_y_continuous(
    limits = c(0, 1.01),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent,
    expand = c(0, 0)
  ) +
  labs(
    x = "All 978 L1000 Landmark Genes",
    y = "Bootstrap Inclusion Probability"
  ) +
  theme_bw(base_size = 6) +
  theme(
    axis.text.y = element_text(size = 0.9, color = "#333333", margin = margin(r = 0.2)),
    axis.text.x = element_text(size = 4.5, color = "#222222"),
    axis.title = element_text(size = 5.5, face = "bold", color = "#111111"),
    legend.position = "top",
    legend.key.width = unit(0.6, "cm"),
    legend.key.height = unit(0.12, "cm"),
    legend.title = element_text(size = 4.8, face = "bold"),
    legend.text = element_text(size = 4.2),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 3, r = 5, b = 3, l = 3)
  )
ggsave(
  filename = "All_978_Genes_SingleColumn.pdf",
  plot = p_fig2_appendix,
  width = 3.35,
  height = 25,
  units = "in",
  limitsize = FALSE,
  device = cairo_pdf
)
