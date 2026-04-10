plot_bad_point_prob_distributions <- function(csv_path,
                                              output_dir = dirname(csv_path),
                                              prefix = "bad_points",
                                              bins = 35,
                                              write_per_state = TRUE,
                                              max_states = 12,
                                              verbose = TRUE) {
  if (!file.exists(csv_path)) {
    stop("Input CSV does not exist: ", normalizePath(csv_path, mustWork = FALSE))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    stop("Input CSV has no rows: ", csv_path)
  }

  required_cols <- c("obs_used", "forecast_at_obs", "analysis_at_obs")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in input CSV: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  for (cc in required_cols) {
    df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  }

  make_long <- function(d) {
    out <- rbind(
      data.frame(
        series = "Observation",
        value = d$obs_used,
        state_var = d$state_var %||% NA_character_,
        site_id = d$site_id %||% NA_character_,
        stringsAsFactors = FALSE
      ),
      data.frame(
        series = "Forecast",
        value = d$forecast_at_obs,
        state_var = d$state_var %||% NA_character_,
        site_id = d$site_id %||% NA_character_,
        stringsAsFactors = FALSE
      ),
      data.frame(
        series = "Analysis",
        value = d$analysis_at_obs,
        state_var = d$state_var %||% NA_character_,
        site_id = d$site_id %||% NA_character_,
        stringsAsFactors = FALSE
      )
    )
    out$series <- factor(out$series, levels = c("Observation", "Forecast", "Analysis"))
    out <- out[is.finite(out$value), , drop = FALSE]
    out
  }

  `%||%` <- function(x, y) if (is.null(x)) y else x
  long_df <- make_long(df)
  if (nrow(long_df) == 0) {
    stop("No finite values found in obs/forecast/analysis columns.")
  }

  summarize_one <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(c(
        n = 0,
        mean = NA_real_,
        sd = NA_real_,
        median = NA_real_,
        p05 = NA_real_,
        p25 = NA_real_,
        p75 = NA_real_,
        p95 = NA_real_
      ))
    }
    c(
      n = length(x),
      mean = mean(x),
      sd = if (length(x) > 1) stats::sd(x) else NA_real_,
      median = stats::median(x),
      p05 = stats::quantile(x, 0.05, names = FALSE),
      p25 = stats::quantile(x, 0.25, names = FALSE),
      p75 = stats::quantile(x, 0.75, names = FALSE),
      p95 = stats::quantile(x, 0.95, names = FALSE)
    )
  }

  stats_tbl <- do.call(
    rbind,
    lapply(
      split(long_df$value, long_df$series),
      function(x) as.data.frame(as.list(summarize_one(x)), stringsAsFactors = FALSE)
    )
  )
  stats_tbl$series <- rownames(stats_tbl)
  rownames(stats_tbl) <- NULL
  stats_tbl <- stats_tbl[, c("series", "n", "mean", "sd", "median", "p05", "p25", "p75", "p95")]

  stats_csv <- file.path(output_dir, paste0(prefix, "_distribution_summary.csv"))
  utils::write.csv(stats_tbl, stats_csv, row.names = FALSE)

  color_vals <- c(
    Observation = "#1f78b4",
    Forecast = "#ff7f00",
    Analysis = "#33a02c"
  )

  p_density <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = value, color = series, fill = series)
  ) +
    ggplot2::geom_density(alpha = 0.18, linewidth = 0.9, adjust = 1.05) +
    ggplot2::scale_color_manual(values = color_vals) +
    ggplot2::scale_fill_manual(values = color_vals) +
    ggplot2::labs(
      title = "Probability density at bad assimilation points",
      subtitle = sprintf("Input rows: %d bad points", nrow(df)),
      x = "Value",
      y = "Density",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p_density_file <- file.path(output_dir, paste0(prefix, "_density_overlay.png"))
  ggplot2::ggsave(p_density_file, p_density, width = 10.5, height = 6.2, dpi = 180)

  p_hist <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = value, y = ggplot2::after_stat(density), fill = series, color = series)
  ) +
    ggplot2::geom_histogram(alpha = 0.35, bins = bins, linewidth = 0.25) +
    ggplot2::geom_density(alpha = 0.06, linewidth = 0.9) +
    ggplot2::facet_wrap(~series, scales = "free_y", ncol = 1) +
    ggplot2::scale_color_manual(values = color_vals) +
    ggplot2::scale_fill_manual(values = color_vals) +
    ggplot2::labs(
      title = "Distribution by data source at bad assimilation points",
      x = "Value",
      y = "Density"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  p_hist_file <- file.path(output_dir, paste0(prefix, "_hist_density_facets.png"))
  ggplot2::ggsave(p_hist_file, p_hist, width = 10.5, height = 9.2, dpi = 180)

  p_ecdf <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = value, color = series)
  ) +
    ggplot2::stat_ecdf(linewidth = 1.0) +
    ggplot2::scale_color_manual(values = color_vals) +
    ggplot2::labs(
      title = "Empirical cumulative probability at bad assimilation points",
      x = "Value",
      y = "Cumulative probability",
      color = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p_ecdf_file <- file.path(output_dir, paste0(prefix, "_ecdf_overlay.png"))
  ggplot2::ggsave(p_ecdf_file, p_ecdf, width = 10.5, height = 6.2, dpi = 180)

  plot_files <- c(p_density_file, p_hist_file, p_ecdf_file)

  if (isTRUE(write_per_state) && "state_var" %in% names(df)) {
    state_counts <- sort(table(df$state_var), decreasing = TRUE)
    keep_states <- names(state_counts)[seq_len(min(length(state_counts), as.integer(max_states)))]
    keep_states <- keep_states[is.finite(match(keep_states, keep_states))]
    keep_states <- keep_states[nzchar(keep_states)]

    if (length(keep_states) > 0) {
      d_state <- df[df$state_var %in% keep_states, , drop = FALSE]
      long_state <- make_long(d_state)
      long_state$state_var <- factor(long_state$state_var, levels = keep_states)

      p_state <- ggplot2::ggplot(
        long_state,
        ggplot2::aes(x = value, color = series, fill = series)
      ) +
        ggplot2::geom_density(alpha = 0.12, linewidth = 0.7, adjust = 1.05) +
        ggplot2::facet_wrap(~state_var, scales = "free", ncol = 3) +
        ggplot2::scale_color_manual(values = color_vals) +
        ggplot2::scale_fill_manual(values = color_vals) +
        ggplot2::labs(
          title = sprintf(
            "Probability density by state variable (top %d by count)",
            length(keep_states)
          ),
          x = "Value",
          y = "Density",
          color = NULL,
          fill = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(size = 8)
        )

      p_state_file <- file.path(output_dir, paste0(prefix, "_density_by_state.png"))
      ggplot2::ggsave(p_state_file, p_state, width = 12, height = 8.5, dpi = 180)
      plot_files <- c(plot_files, p_state_file)
    }
  }

  if (isTRUE(verbose)) {
    message("Wrote distribution summary: ", stats_csv)
    message("Saved plot files:")
    for (pf in plot_files) {
      message("  - ", pf)
    }
  }

  invisible(list(
    input_csv = csv_path,
    summary_csv = stats_csv,
    plot_files = plot_files,
    summary = stats_tbl
  ))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- strsplit(args, "=", fixed = TRUE)
  keys <- vapply(kv, function(x) x[1], character(1))
  vals <- vapply(kv, function(x) if (length(x) >= 2) x[2] else "", character(1))
  arg <- setNames(as.list(vals), keys)

  get_arg <- function(name, default = NULL) {
    val <- arg[[paste0("--", name)]]
    if (is.null(val) || !nzchar(val)) {
      return(default)
    }
    val
  }

  parse_bool <- function(x, default = TRUE) {
    if (is.null(x) || !nzchar(x)) {
      return(default)
    }
    x <- tolower(trimws(x))
    if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
    if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
    default
  }

  csv_path <- get_arg("csv_path", NULL)
  outdir <- get_arg("outdir", NULL)
  if (is.null(csv_path) && !is.null(outdir)) {
    csv_path <- file.path(outdir, "bad_assimilation_diagnostics", "bad_assimilation_points.csv")
  }
  if (is.null(csv_path)) {
    stop("Provide --csv_path=/path/to/bad_assimilation_points.csv or --outdir=/path/to/sda_output")
  }

  plot_bad_point_prob_distributions(
    csv_path = csv_path,
    output_dir = get_arg("output_dir", dirname(csv_path)),
    prefix = get_arg("prefix", "bad_points"),
    bins = as.integer(get_arg("bins", "35")),
    write_per_state = parse_bool(get_arg("write_per_state", "true"), TRUE),
    max_states = as.integer(get_arg("max_states", "12")),
    verbose = TRUE
  )
}

# Example:
# Rscript plot_bad_point_prob_distributions.R \
#   --outdir=/path/to/output_inter_q \
#   --output_dir=/path/to/output_inter_q/bad_assimilation_diagnostics \
#   --prefix=bad_points \
#   --bins=35 \
#   --write_per_state=true \
#   --max_states=12
