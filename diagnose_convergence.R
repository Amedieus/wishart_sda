library(data.table)
library(ggplot2)

# ==========
outdir <- "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q_4"

# ==========
diag_files <- list.files(
  path = outdir,
  pattern = "^mcmc_diag_t[0-9]+\\.csv$",
  full.names = TRUE
)

if (length(diag_files) == 0) {
  stop("没有找到任何 mcmc_diag_t*.csv 文件，请检查 outdir 路径是否正确。")
}

# ===== 读取并合并 =====
diag_dt <- rbindlist(lapply(diag_files, fread), fill = TRUE)

# 从文件名里提取 t（双保险；即使 csv 里已有 t_index 也可以核对）
diag_dt[, file_name := basename(diag_files[match(.I, seq_len(.N))]), by = .I]

# 更稳妥的方式：重新按文件循环读并加 file_t
diag_dt <- rbindlist(lapply(diag_files, function(f) {
  dt <- fread(f)
  dt[, diag_file := basename(f)]
  dt[, file_t := as.integer(sub("^mcmc_diag_t([0-9]+)\\.csv$", "\\1", basename(f)))]
  dt
}), fill = TRUE)

# 看看结构
print(names(diag_dt))
print(head(diag_dt))

# ===== 整体摘要 =====
summary(diag_dt$rhat_max)
summary(diag_dt$rhat_median)
summary(diag_dt$ess_min)
summary(diag_dt$ess_median)

# ===== NA 情况 =====
colSums(is.na(diag_dt))

# ===== 收敛标签 =====
diag_dt[, rhat_ok := !is.na(rhat_max) & rhat_max <= 1.05]
diag_dt[, ess_ok  := !is.na(ess_min)  & ess_min  >= 100]
diag_dt[, conv_ok := rhat_ok & ess_ok]

# 看比例
diag_dt[, .(
  n_blocks = .N,
  rhat_ok_n = sum(rhat_ok, na.rm = TRUE),
  ess_ok_n  = sum(ess_ok,  na.rm = TRUE),
  conv_ok_n = sum(conv_ok, na.rm = TRUE),
  rhat_ok_pct = mean(rhat_ok, na.rm = TRUE),
  ess_ok_pct  = mean(ess_ok,  na.rm = TRUE),
  conv_ok_pct = mean(conv_ok, na.rm = TRUE)
)]

# ===== 有问题的 block =====
bad_blocks <- diag_dt[
  is.na(rhat_max) | is.na(ess_min) | rhat_max > 1.05 | ess_min < 100
][order(file_t, block_id)]

print(bad_blocks)

# 保存一份方便你后面看
fwrite(bad_blocks, file.path(outdir, "mcmc_diag_bad_blocks.csv"))

# ===== 每个时间点汇总 =====
diag_by_time <- diag_dt[, .(
  n_blocks = .N,
  rhat_max_max = max(rhat_max, na.rm = TRUE),
  rhat_max_median = median(rhat_max, na.rm = TRUE),
  ess_min_min = min(ess_min, na.rm = TRUE),
  ess_min_median = median(ess_min, na.rm = TRUE),
  bad_block_n = sum(is.na(rhat_max) | is.na(ess_min) | rhat_max > 1.05 | ess_min < 100)
), by = .(file_t)]

diag_by_time <- diag_by_time[order(file_t)]

print(diag_by_time)

fwrite(diag_by_time, file.path(outdir, "mcmc_diag_by_time_summary.csv"))

# ===== 哪些 site 组合经常出问题 =====
diag_by_site <- diag_dt[, .(
  n_times = .N,
  worst_rhat = max(rhat_max, na.rm = TRUE),
  worst_ess  = min(ess_min, na.rm = TRUE),
  bad_n = sum(is.na(rhat_max) | is.na(ess_min) | rhat_max > 1.05 | ess_min < 100)
), by = .(site_ids)][order(-bad_n, -worst_rhat, worst_ess)]

print(head(diag_by_site, 30))

fwrite(diag_by_site, file.path(outdir, "mcmc_diag_by_site_summary.csv"))

p1 <- ggplot(diag_by_time, aes(x = file_t, y = rhat_max_max)) +
  geom_line() +
  geom_hline(yintercept = 1.05, linetype = 2) +
  labs(
    title = "Worst Rhat by time step",
    x = "Time index (t)",
    y = "Max Rhat"
  ) +
  theme_bw()

print(p1)
ggsave(file.path(outdir, "mcmc_diag_plot_rhat_by_time.png"), p1, width = 8, height = 5)

p2 <- ggplot(diag_by_time, aes(x = file_t, y = ess_min_min)) +
  geom_line() +
  geom_hline(yintercept = 100, linetype = 2) +
  labs(
    title = "Worst ESS by time step",
    x = "Time index (t)",
    y = "Min ESS"
  ) +
  theme_bw()

print(p2)
ggsave(file.path(outdir, "mcmc_diag_plot_ess_by_time.png"), p2, width = 8, height = 5)

p3 <- ggplot(diag_dt[!is.na(rhat_max)], aes(x = rhat_max)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 1.05, linetype = 2) +
  labs(
    title = "Distribution of block-wise max Rhat",
    x = "Rhat max",
    y = "Count"
  ) +
  theme_bw()

print(p3)
ggsave(file.path(outdir, "mcmc_diag_plot_rhat_hist.png"), p3, width = 8, height = 5)

p4 <- ggplot(diag_dt[!is.na(ess_min)], aes(x = ess_min)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 100, linetype = 2) +
  labs(
    title = "Distribution of block-wise min ESS",
    x = "ESS min",
    y = "Count"
  ) +
  theme_bw()

print(p4)
ggsave(file.path(outdir, "mcmc_diag_plot_ess_hist.png"), p4, width = 8, height = 5)

# ===== 自动诊断结论 =====
overall <- diag_dt[, .(
  total_blocks = .N,
  pct_rhat_ok = mean(rhat_max <= 1.05, na.rm = TRUE),
  pct_ess_ok  = mean(ess_min >= 100, na.rm = TRUE),
  pct_both_ok = mean(rhat_max <= 1.05 & ess_min >= 100, na.rm = TRUE),
  worst_rhat  = max(rhat_max, na.rm = TRUE),
  worst_ess   = min(ess_min, na.rm = TRUE)
)]

print(overall)

if (overall$pct_both_ok > 0.95 && overall$worst_rhat <= 1.05 && overall$worst_ess >= 100) {
  cat("结论：整体上 MCMC 基本收敛，绝大多数 block 可接受。\n")
} else if (overall$pct_both_ok > 0.8) {
  cat("结论：大部分 block 可接受，但仍有一部分 block 收敛不佳，需要重点检查 bad blocks。\n")
} else {
  cat("结论：MCMC 收敛情况不理想，较多 block 的 Rhat 或 ESS 不达标。\n")
}