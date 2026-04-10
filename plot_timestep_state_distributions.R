plot_timestep_state_distributions <- function(outdir,
                                              timestep = NULL,
                                              index = NULL,
                                              bad_points_csv = NULL,
                                              index_column = "rank",
                                              output_dir = file.path(outdir, "timestep_distribution_plots"),
                                              prefix = NULL,
                                              bins = 30,
                                              max_states = 0,
                                              verbose = TRUE) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  parse_timestep <- function(path) {
    as.integer(sub("^sda\\.output([0-9]+)\\.Rdata$", "\\1", basename(path)))
  }

  as_named_vector <- function(x, prefix = "state") {
    x <- suppressWarnings(as.numeric(x))
    nm <- names(x)
    if (is.null(nm)) {
      nm <- paste0(prefix, "_", seq_along(x))
    } else {
      bad <- is.na(nm) | !nzchar(nm)
      nm[bad] <- paste0(prefix, "_", which(bad))
    }
    names(x) <- nm
    x
  }

  make_site_map <- function(site_ids, n_state) {
    site_ids <- as.character(site_ids)
    if (length(site_ids) == 0 || n_state <= 0) {
      return(rep(NA_character_, n_state))
    }
    if (n_state %% length(site_ids) == 0) {
      return(rep(site_ids, each = n_state / length(site_ids)))
    }
    rep(site_ids, length.out = n_state)
  }

  is_block <- function(x) {
    is.list(x) && !is.null(x$data) && !is.null(x$constant)
  }

  is_block_list <- function(x) {
    is.list(x) && length(x) > 0 && is_block(x[[1]])
  }

  resolve_timestep <- function(timestep, index, bad_points_csv, index_column) {
    if (!is.null(timestep) && nzchar(as.character(timestep))) {
      ts <- suppressWarnings(as.integer(timestep))
      if (!is.finite(ts) || ts <= 0) {
        stop("Invalid --timestep value: ", timestep)
      }
      return(list(
        timestep = ts,
        selection_row = NULL,
        selection_note = "Resolved from --timestep."
      ))
    }

    idx <- suppressWarnings(as.integer(index))
    if (!is.finite(idx) || idx <= 0) {
      stop("Provide either --timestep, or --index with --bad_points_csv.")
    }
    if (is.null(bad_points_csv) || !file.exists(bad_points_csv)) {
      stop("bad_points_csv not found: ", bad_points_csv %||% "NULL")
    }

    tbl <- utils::read.csv(bad_points_csv, stringsAsFactors = FALSE)
    if (nrow(tbl) == 0) {
      stop("bad_points_csv has no rows: ", bad_points_csv)
    }
    if (!("timestep" %in% names(tbl))) {
      stop("bad_points_csv must contain column 'timestep'.")
    }

    row_id <- integer(0)
    if (!is.null(index_column) && nzchar(index_column) && (index_column %in% names(tbl))) {
      col_num <- suppressWarnings(as.integer(tbl[[index_column]]))
      row_id <- which(col_num == idx)
      if (length(row_id) == 0) {
        row_id <- which(as.character(tbl[[index_column]]) == as.character(index))
      }
    }

    if (length(row_id) == 0 && idx <= nrow(tbl)) {
      row_id <- idx
    }
    if (length(row_id) == 0) {
      stop(
        "Could not resolve index=", idx,
        " with index_column='", index_column, "' in: ", bad_points_csv
      )
    }

    sel <- tbl[row_id[1], , drop = FALSE]
    ts <- suppressWarnings(as.integer(sel$timestep[1]))
    if (!is.finite(ts) || ts <= 0) {
      stop("Resolved row has invalid timestep.")
    }

    list(
      timestep = ts,
      selection_row = sel,
      selection_note = sprintf(
        "Resolved from index=%d (%s) in bad_points_csv.",
        idx,
        if (index_column %in% names(tbl)) index_column else "row number"
      )
    )
  }

  extract_block_list <- function(enkf, timestep) {
    block_all <- enkf$block.list.all
    if (is.null(block_all)) {
      return(NULL)
    }

    # Some outputs store block list directly (not nested by timestep).
    if (is_block_list(block_all)) {
      return(block_all)
    }

    # Standard case: block.list.all[[timestep]] is a block list.
    if (length(block_all) >= timestep && is_block_list(block_all[[timestep]])) {
      return(block_all[[timestep]])
    }

    # Fallback to the last non-null block list.
    valid <- which(vapply(block_all, is_block_list, logical(1)))
    if (length(valid) > 0) {
      return(block_all[[max(valid)]])
    }

    NULL
  }

  summarize_values <- function(x) {
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

  if (!dir.exists(outdir)) {
    stop("outdir does not exist: ", normalizePath(outdir, mustWork = FALSE))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  resolved <- resolve_timestep(
    timestep = timestep,
    index = index,
    bad_points_csv = bad_points_csv,
    index_column = index_column
  )
  timestep <- resolved$timestep
  if (is.null(prefix) || !nzchar(prefix)) {
    prefix <- sprintf("timestep_%04d", timestep)
  }

  sda_files <- list.files(
    outdir,
    pattern = "^sda\\.output[0-9]+\\.Rdata$",
    full.names = TRUE
  )
  if (length(sda_files) == 0) {
    stop("No sda.output*.Rdata files found in outdir: ", outdir)
  }
  sda_files <- sda_files[order(parse_timestep(sda_files))]
  ts_all <- parse_timestep(sda_files)
  file_idx <- which(ts_all == timestep)
  if (length(file_idx) == 0) {
    stop(
      "Could not find file for timestep=", timestep,
      ". Available timesteps include: ", paste(head(ts_all, 20), collapse = ", ")
    )
  }
  sda_file <- sda_files[file_idx[1]]

  env <- new.env(parent = emptyenv())
  load(sda_file, envir = env)
  if (!exists("sda.outputs", envir = env, inherits = FALSE)) {
    stop("Missing object 'sda.outputs' in file: ", basename(sda_file))
  }
  sda_outputs <- get("sda.outputs", envir = env)
  enkf <- sda_outputs$enkf.params
  if (is.null(enkf)) {
    stop("Missing sda.outputs$enkf.params in file: ", basename(sda_file))
  }

  block_list <- extract_block_list(enkf, timestep = timestep)
  if (is.null(block_list) || length(block_list) == 0) {
    stop("No valid block list found for timestep ", timestep)
  }

  rows <- list()
  for (b in seq_along(block_list)) {
    block <- block_list[[b]]
    y_obs <- suppressWarnings(as.numeric(block$data$y.censored))
    h_idx <- suppressWarnings(as.integer(block$constant$H))
    muf <- as_named_vector(block$data$muf %||% numeric(), prefix = "state")
    mufa <- as_named_vector(block$update$mufa %||% numeric(), prefix = "state")
    n_state <- length(muf)
    if (n_state == 0 || length(y_obs) == 0 || length(h_idx) == 0) {
      next
    }

    n_obs <- min(length(y_obs), length(h_idx))
    y_obs <- y_obs[seq_len(n_obs)]
    h_idx <- h_idx[seq_len(n_obs)]
    valid <- is.finite(h_idx) & h_idx > 0 & h_idx <= n_state

    state_names <- names(muf)
    forecast <- rep(NA_real_, n_obs)
    analysis <- rep(NA_real_, n_obs)
    forecast[valid] <- suppressWarnings(as.numeric(muf[h_idx[valid]]))
    analysis[valid] <- suppressWarnings(as.numeric(mufa[h_idx[valid]]))

    site_ids <- as.character(block$site.ids %||% paste0("block_", b))
    site_map <- make_site_map(site_ids, n_state)
    state_var <- rep(NA_character_, n_obs)
    site_id <- rep(NA_character_, n_obs)
    state_var[valid] <- state_names[h_idx[valid]]
    site_id[valid] <- site_map[h_idx[valid]]

    rows[[length(rows) + 1L]] <- data.frame(
      timestep = timestep,
      block_id = b,
      obs_idx = seq_len(n_obs),
      state_idx = h_idx,
      state_var = state_var,
      site_id = site_id,
      obs = y_obs,
      forecast = forecast,
      analysis = analysis,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {
    stop("No observation rows could be extracted for timestep ", timestep)
  }

  points_df <- do.call(rbind, rows)
  keep <- is.finite(points_df$obs) | is.finite(points_df$forecast) | is.finite(points_df$analysis)
  points_df <- points_df[keep, , drop = FALSE]
  if (nrow(points_df) == 0) {
    stop("No finite obs/forecast/analysis values extracted.")
  }

  long_df <- rbind(
    data.frame(
      timestep = points_df$timestep,
      block_id = points_df$block_id,
      obs_idx = points_df$obs_idx,
      state_var = points_df$state_var,
      site_id = points_df$site_id,
      source = "Observation",
      value = points_df$obs,
      stringsAsFactors = FALSE
    ),
    data.frame(
      timestep = points_df$timestep,
      block_id = points_df$block_id,
      obs_idx = points_df$obs_idx,
      state_var = points_df$state_var,
      site_id = points_df$site_id,
      source = "Forecast",
      value = points_df$forecast,
      stringsAsFactors = FALSE
    ),
    data.frame(
      timestep = points_df$timestep,
      block_id = points_df$block_id,
      obs_idx = points_df$obs_idx,
      state_var = points_df$state_var,
      site_id = points_df$site_id,
      source = "Analysis",
      value = points_df$analysis,
      stringsAsFactors = FALSE
    )
  )
  long_df$source <- factor(long_df$source, levels = c("Observation", "Forecast", "Analysis"))
  long_df <- long_df[is.finite(long_df$value), , drop = FALSE]
  if (nrow(long_df) == 0) {
    stop("No finite values available for plotting.")
  }

  if (!("state_var" %in% names(long_df))) {
    long_df$state_var <- "unknown"
  }
  long_df$state_var[is.na(long_df$state_var) | !nzchar(long_df$state_var)] <- "unknown"

  state_levels <- names(sort(table(long_df$state_var), decreasing = TRUE))
  if (is.finite(max_states) && max_states > 0 && length(state_levels) > max_states) {
    state_levels <- state_levels[seq_len(max_states)]
    long_df <- long_df[long_df$state_var %in% state_levels, , drop = FALSE]
  }
  long_df$state_var <- factor(long_df$state_var, levels = state_levels)

  overall_split <- split(long_df$value, as.character(long_df$source))
  overall_summary <- do.call(
    rbind,
    lapply(names(overall_split), function(nm) {
      sm <- summarize_values(overall_split[[nm]])
      out <- as.data.frame(as.list(sm), stringsAsFactors = FALSE)
      out$source <- nm
      out
    })
  )
  overall_summary <- overall_summary[, c("source", "n", "mean", "sd", "median", "p05", "p25", "p75", "p95")]

  by_state_summary <- do.call(
    rbind,
    lapply(split(long_df, list(long_df$state_var, long_df$source), drop = TRUE), function(d) {
      sm <- summarize_values(d$value)
      out <- as.data.frame(as.list(sm), stringsAsFactors = FALSE)
      out$state_var <- as.character(d$state_var[1])
      out$source <- as.character(d$source[1])
      out
    })
  )
  by_state_summary <- by_state_summary[
    ,
    c("state_var", "source", "n", "mean", "sd", "median", "p05", "p25", "p75", "p95")
  ]
  by_state_summary <- by_state_summary[order(by_state_summary$state_var, by_state_summary$source), , drop = FALSE]

  points_csv <- file.path(output_dir, paste0(prefix, "_obs_forecast_analysis_points.csv"))
  summary_overall_csv <- file.path(output_dir, paste0(prefix, "_summary_overall.csv"))
  summary_state_csv <- file.path(output_dir, paste0(prefix, "_summary_by_state.csv"))
  utils::write.csv(points_df, points_csv, row.names = FALSE)
  utils::write.csv(overall_summary, summary_overall_csv, row.names = FALSE)
  utils::write.csv(by_state_summary, summary_state_csv, row.names = FALSE)

  colors <- c(
    Observation = "#1f78b4",
    Forecast = "#ff7f00",
    Analysis = "#33a02c"
  )

  p_overall <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = value, color = source, fill = source)
  ) +
    ggplot2::geom_density(alpha = 0.16, linewidth = 0.9, adjust = 1.05) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = sprintf("Timestep %d: combined distribution", timestep),
      subtitle = "Observation vs Forecast vs Analysis across all state variables",
      x = "Value",
      y = "Density",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  overall_png <- file.path(output_dir, paste0(prefix, "_density_overall.png"))
  ggplot2::ggsave(overall_png, p_overall, width = 11, height = 6.5, dpi = 180)

  n_state <- length(state_levels)
  ncol_facets <- if (n_state <= 4) 2 else if (n_state <= 9) 3 else 4
  nrow_facets <- max(1, ceiling(n_state / ncol_facets))
  facet_height <- max(6, min(30, 2.2 * nrow_facets + 1.5))

  p_by_state <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = value, y = ggplot2::after_stat(density), color = source, fill = source)
  ) +
    ggplot2::geom_histogram(alpha = 0.28, bins = bins, position = "identity", linewidth = 0.2) +
    ggplot2::geom_density(alpha = 0.08, linewidth = 0.7, adjust = 1.05) +
    ggplot2::facet_wrap(~state_var, scales = "free_y", ncol = ncol_facets) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = sprintf("Timestep %d: distribution by state variable", timestep),
      subtitle = "Each panel is one state variable",
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

  state_png <- file.path(output_dir, paste0(prefix, "_density_by_state.png"))
  ggplot2::ggsave(state_png, p_by_state, width = 14, height = facet_height, dpi = 180)

  p_ecdf <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = value, color = source)
  ) +
    ggplot2::stat_ecdf(linewidth = 0.95) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = sprintf("Timestep %d: ECDF overlay", timestep),
      x = "Value",
      y = "Cumulative probability",
      color = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  ecdf_png <- file.path(output_dir, paste0(prefix, "_ecdf_overall.png"))
  ggplot2::ggsave(ecdf_png, p_ecdf, width = 11, height = 6.5, dpi = 180)

  plot_files <- c(overall_png, state_png, ecdf_png)

  if (isTRUE(verbose)) {
    message(resolved$selection_note)
    message("Selected timestep: ", timestep)
    message("Loaded file: ", basename(sda_file))
    message("Rows extracted: ", nrow(points_df))
    message("State variables plotted: ", n_state)
    message("Output directory: ", normalizePath(output_dir, mustWork = FALSE))
    message("Outputs:")
    message("  - ", points_csv)
    message("  - ", summary_overall_csv)
    message("  - ", summary_state_csv)
    for (pf in plot_files) {
      message("  - ", pf)
    }
  }

  invisible(list(
    timestep = timestep,
    source_file = sda_file,
    selection = resolved,
    point_table = points_df,
    long_table = long_df,
    summary_overall = overall_summary,
    summary_by_state = by_state_summary,
    csv_files = c(points = points_csv, overall = summary_overall_csv, by_state = summary_state_csv),
    plot_files = plot_files
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

  parse_int <- function(x, default = NULL) {
    if (is.null(x) || !nzchar(x)) {
      return(default)
    }
    y <- suppressWarnings(as.integer(x))
    if (!is.finite(y)) default else y
  }

  outdir <- get_arg("outdir", NULL)
  if (is.null(outdir)) {
    stop("Missing required argument: --outdir=/path/to/output_dir")
  }

  timestep <- parse_int(get_arg("timestep", NULL), NULL)
  index <- parse_int(get_arg("index", NULL), NULL)
  max_states <- parse_int(get_arg("max_states", "0"), 0)

  bad_points_csv <- get_arg("bad_points_csv", NULL)
  if (is.null(bad_points_csv)) {
    bad_points_csv <- file.path(outdir, "bad_assimilation_diagnostics", "bad_assimilation_points.csv")
  }

  plot_timestep_state_distributions(
    outdir = outdir,
    timestep = timestep,
    index = index,
    bad_points_csv = bad_points_csv,
    index_column = get_arg("index_column", "rank"),
    output_dir = get_arg("output_dir", file.path(outdir, "timestep_distribution_plots")),
    prefix = get_arg("prefix", NULL),
    bins = parse_int(get_arg("bins", "30"), 30),
    max_states = max_states,
    verbose = TRUE
  )
}

# Example 1: directly specify timestep
# Rscript plot_timestep_state_distributions.R \
#   --outdir=/path/to/output_inter_q \
#   --timestep=25 \
#   --output_dir=/path/to/output_inter_q/timestep_distribution_plots
#
# Example 2: resolve timestep from bad-point rank/index
# Rscript plot_timestep_state_distributions.R \
#   --outdir=/path/to/output_inter_q \
#   --index=1 \
#   --bad_points_csv=/path/to/output_inter_q/bad_assimilation_diagnostics/bad_assimilation_points.csv \
#   --index_column=rank
