library(data.table)
library(ggplot2)

nee_model_to_obs <- function(x) {
  x * 1e8 / 1.157407
}

plot_point_density <- function(sda.outputs, obs.mean = NULL, obs.cov = NULL, index,
                               state_order = c("AbvGrndWood", "NEE", "Qle", "LAI", "SoilMoistFrac", "TotSoilCarb"),
                               n_grid = 400) {
  
  forecast_df <- sda.outputs$forecast
  analysis_df <- sda.outputs$analysis
  
  if (is.null(forecast_df) || is.null(analysis_df)) {
    stop("forecast 或 analysis 不存在。请确认你传入的是单个 timestep 的 sda.outputs 对象。")
  }
  
  if (is.null(obs.mean)) obs.mean <- sda.outputs$obs.mean
  if (is.null(obs.cov))  obs.cov  <- sda.outputs$obs.cov
  
  if (is.null(obs.mean) || is.null(obs.cov)) {
    stop("obs.mean 或 obs.cov 不存在。")
  }
  
  forecast_mat <- as.matrix(forecast_df)
  analysis_mat <- as.matrix(analysis_df)
  
  n_sites  <- length(obs.mean)
  n_states <- length(state_order)
  
  if (ncol(forecast_mat) != n_sites * n_states) {
    stop(sprintf("forecast列数=%d, 但 n_sites * n_states = %d * %d = %d，不匹配。",
                 ncol(forecast_mat), n_sites, n_states, n_sites * n_states))
  }
  
  if (index < 1 || index > n_sites) {
    stop(sprintf("index 超出范围。当前 index 必须在 1 到 %d 之间。", n_sites))
  }
  
  col_start <- (index - 1) * n_states + 1
  col_end   <- index * n_states
  cols_this <- col_start:col_end
  
  fc_sub <- forecast_mat[, cols_this, drop = FALSE]
  an_sub <- analysis_mat[, cols_this, drop = FALSE]
  colnames(fc_sub) <- state_order
  colnames(an_sub) <- state_order
  
  # ===== 关键修改：把 analysis 里的 NEE 转成 obs 单位 =====
  if ("NEE" %in% colnames(an_sub)) {
    an_sub[, "NEE"] <- nee_model_to_obs(an_sub[, "NEE"])
  }
  
  # 如果你想连 forecast 的 NEE 也一起转，就把下面这段取消注释
  # if ("NEE" %in% colnames(fc_sub)) {
  #   fc_sub[, "NEE"] <- nee_model_to_obs(fc_sub[, "NEE"])
  # }
  
  obs_mean_i <- obs.mean[[index]]
  obs_cov_i  <- obs.cov[[index]]
  
  if (is.matrix(obs_mean_i) || is.data.frame(obs_mean_i)) {
    if (nrow(obs_mean_i) == 1) {
      obs_mu <- as.numeric(obs_mean_i[1, ])
      obs_names <- colnames(obs_mean_i)
    } else if (ncol(obs_mean_i) == 1) {
      obs_mu <- as.numeric(obs_mean_i[, 1])
      obs_names <- rownames(obs_mean_i)
    } else {
      stop("obs.mean[[index]] 格式不支持。")
    }
  } else {
    obs_mu <- as.numeric(obs_mean_i)
    obs_names <- names(obs_mean_i)
  }
  
  if (is.null(obs_names)) {
    obs_names <- paste0("obs_", seq_along(obs_mu))
  }
  
  if (length(obs_mu) == 1) {
    if (length(obs_cov_i) == 1) {
      obs_sd <- sqrt(as.numeric(obs_cov_i))
    } else if (is.matrix(obs_cov_i) && all(dim(obs_cov_i) == c(1, 1))) {
      obs_sd <- sqrt(obs_cov_i[1, 1])
    } else {
      stop("单变量 obs.cov 格式不合法。")
    }
  } else {
    if (!is.matrix(obs_cov_i)) {
      stop("多变量 obs.cov 应该是协方差矩阵。")
    }
    if (nrow(obs_cov_i) != length(obs_mu) || ncol(obs_cov_i) != length(obs_mu)) {
      stop("obs.cov 维度和 obs.mean 长度不一致。")
    }
    obs_sd <- sqrt(diag(obs_cov_i))
  }
  
  obs_dt <- data.table(
    variable = obs_names,
    obs_mean = obs_mu,
    obs_sd   = obs_sd
  )
  
  density_list <- list()
  
  for (v in state_order) {
    fc_vals <- fc_sub[, v]
    an_vals <- an_sub[, v]
    
    fc_vals <- fc_vals[is.finite(fc_vals)]
    an_vals <- an_vals[is.finite(an_vals)]
    
    obs_row <- obs_dt[variable == v]
    obs_has <- nrow(obs_row) == 1 && is.finite(obs_row$obs_sd) && obs_row$obs_sd > 0
    
    x_candidates <- c(fc_vals, an_vals)
    if (obs_has) {
      x_candidates <- c(
        x_candidates,
        obs_row$obs_mean - 4 * obs_row$obs_sd,
        obs_row$obs_mean + 4 * obs_row$obs_sd
      )
    }
    
    x_min <- min(x_candidates, na.rm = TRUE)
    x_max <- max(x_candidates, na.rm = TRUE)
    
    if (!is.finite(x_min) || !is.finite(x_max) || x_min == x_max) {
      x_min <- x_min - 1
      x_max <- x_max + 1
    }
    
    x_grid <- seq(x_min, x_max, length.out = n_grid)
    
    if (length(unique(fc_vals)) >= 2) {
      d_fc <- density(fc_vals, from = x_min, to = x_max, n = n_grid, na.rm = TRUE)
      dt_fc <- data.table(variable = v, source = "forecast", x = d_fc$x, y = d_fc$y)
    } else {
      s <- sd(fc_vals, na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1e-6
      dt_fc <- data.table(
        variable = v, source = "forecast", x = x_grid,
        y = dnorm(x_grid, mean(fc_vals, na.rm = TRUE), s)
      )
    }
    
    if (length(unique(an_vals)) >= 2) {
      d_an <- density(an_vals, from = x_min, to = x_max, n = n_grid, na.rm = TRUE)
      dt_an <- data.table(variable = v, source = "analysis", x = d_an$x, y = d_an$y)
    } else {
      s <- sd(an_vals, na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1e-6
      dt_an <- data.table(
        variable = v, source = "analysis", x = x_grid,
        y = dnorm(x_grid, mean(an_vals, na.rm = TRUE), s)
      )
    }
    
    if (obs_has) {
      dt_obs <- data.table(
        variable = v,
        source   = "obs",
        x        = x_grid,
        y        = dnorm(x_grid, mean = obs_row$obs_mean, sd = obs_row$obs_sd)
      )
      density_list[[length(density_list) + 1]] <- rbind(dt_fc, dt_an, dt_obs)
    } else {
      density_list[[length(density_list) + 1]] <- rbind(dt_fc, dt_an)
    }
  }
  
  density_dt <- rbindlist(density_list)
  
  ggplot(density_dt, aes(x = x, y = y, color = source)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ variable, scales = "free", ncol = 2) +
    theme_bw(base_size = 13) +
    labs(
      title = paste0("Index = ", index),
      x = "Value",
      y = "Density",
      color = NULL
    )
}

p <- plot_point_density(
  sda.outputs = sda.outputs,
  index = 1
)

print(p)