plot_timestep_state_distributions <- function(outdir,
                                              timestep = NULL,
                                              index = NULL,
                                              bad_points_csv = NULL,
                                              index_column = "rank",
                                              output_dir = file.path(outdir, "timestep_distribution_plots"),
                                              prefix = NULL,
                                              bins = 30,
                                              site_ids = NULL,
                                              facet_mode = c("site_state", "state", "site"),
                                              max_states = 0,
                                              max_facets = 0,
                                              verbose = TRUE) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  parse_timestep <- function(path) {
    as.integer(sub("^sda\\.output([0-9]+)\\.Rdata$", "\\1", basename(path)))
  }

  parse_csv_values <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    if (length(x) > 1) {
      x <- paste(x, collapse = ",")
    }
    x <- trimws(as.character(x))
    if (!nzchar(x)) {
      return(NULL)
    }
    vals <- strsplit(x, ",", fixed = TRUE)[[1]]
    vals <- trimws(vals)
    vals[nzchar(vals)]
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

    n_state <- max(length(muf), length(mufa))
    if (n_state == 0) {
      next
    }

    state_names <- rep("", n_state)
    if (length(muf) > 0) {
      idx <- seq_len(min(length(muf), n_state))
      state_names[idx] <- names(muf)[idx]
    }
    if (length(mufa) > 0) {
      idx <- seq_len(min(length(mufa), n_state))
      mufa_names <- names(mufa)[idx]
      need_fill <- !nzchar(state_names[idx]) | is.na(state_names[idx])
      state_names[idx][need_fill] <- mufa_names[need_fill]
    }
    bad_nm <- !nzchar(state_names) | is.na(state_names)
    state_names[bad_nm] <- paste0("state_", which(bad_nm))

    forecast <- rep(NA_real_, n_state)
    analysis <- rep(NA_real_, n_state)
    if (length(muf) > 0) {
      idx <- seq_len(min(length(muf), n_state))
      forecast[idx] <- suppressWarnings(as.numeric(muf[idx]))
    }
    if (length(mufa) > 0) {
      idx <- seq_len(min(length(mufa), n_state))
      analysis[idx] <- suppressWarnings(as.numeric(mufa[idx]))
    }

    obs_state <- rep(NA_real_, n_state)
    n_obs_mapped <- integer(n_state)
    if (length(y_obs) > 0 && length(h_idx) > 0) {
      n_obs <- min(length(y_obs), length(h_idx))
      y_obs <- y_obs[seq_len(n_obs)]
      h_idx <- h_idx[seq_len(n_obs)]
      valid <- is.finite(h_idx) & h_idx > 0 & h_idx <= n_state & is.finite(y_obs)
      if (any(valid)) {
        obs_agg <- tapply(y_obs[valid], h_idx[valid], function(v) mean(v, na.rm = TRUE))
        obs_idx_int <- suppressWarnings(as.integer(names(obs_agg)))
        obs_ok <- is.finite(obs_idx_int) & obs_idx_int > 0 & obs_idx_int <= n_state
        if (any(obs_ok)) {
          obs_state[obs_idx_int[obs_ok]] <- as.numeric(obs_agg[obs_ok])
        }

        obs_n <- table(h_idx[valid])
        obs_n_idx <- suppressWarnings(as.integer(names(obs_n)))
        cnt_ok <- is.finite(obs_n_idx) & obs_n_idx > 0 & obs_n_idx <= n_state
        if (any(cnt_ok)) {
          n_obs_mapped[obs_n_idx[cnt_ok]] <- as.integer(obs_n[cnt_ok])
        }
      }
    }

    site_ids <- as.character(block$site.ids %||% paste0("block_", b))
    site_map <- make_site_map(site_ids, n_state)
    state_idx <- seq_len(n_state)

    rows[[length(rows) + 1L]] <- data.frame(
      timestep = timestep,
      block_id = b,
      obs_idx = NA_integer_,
      state_idx = state_idx,
      state_var = state_names,
      site_id = site_map,
      n_obs_mapped = n_obs_mapped,
      obs = obs_state,
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

  points_df$site_id <- as.character(points_df$site_id)
  points_df$site_id[is.na(points_df$site_id) | !nzchar(points_df$site_id)] <- "unknown_site"
  points_df$state_var <- as.character(points_df$state_var)
  points_df$state_var[is.na(points_df$state_var) | !nzchar(points_df$state_var)] <- "unknown_state"

  requested_sites <- parse_csv_values(site_ids)
  if (!is.null(requested_sites) && length(requested_sites) > 0) {
    points_df <- points_df[points_df$site_id %in% requested_sites, , drop = FALSE]
    if (nrow(points_df) == 0) {
      stop(
        "No rows remained after site filter. Requested site_ids: ",
        paste(requested_sites, collapse = ", ")
      )
    }
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

  long_df$site_id <- as.character(long_df$site_id)
  long_df$site_id[is.na(long_df$site_id) | !nzchar(long_df$site_id)] <- "unknown_site"
  long_df$state_var <- as.character(long_df$state_var)
  long_df$state_var[is.na(long_df$state_var) | !nzchar(long_df$state_var)] <- "unknown_state"

  # Optionally keep only the most represented state variables.
  state_levels_all <- names(sort(table(long_df$state_var), decreasing = TRUE))
  if (is.finite(max_states) && max_states > 0 && length(state_levels_all) > max_states) {
    keep_states <- state_levels_all[seq_len(max_states)]
    long_df <- long_df[long_df$state_var %in% keep_states, , drop = FALSE]
  }
  if (nrow(long_df) == 0) {
    stop("No rows remained after max_states filtering.")
  }
  state_levels <- names(sort(table(long_df$state_var), decreasing = TRUE))
  long_df$state_var <- factor(long_df$state_var, levels = state_levels)

  long_df$facet_site <- as.character(long_df$site_id)
  long_df$facet_state <- as.character(long_df$state_var)
  long_df$facet_site_state <- paste(long_df$facet_site, long_df$facet_state, sep = " | ")

  facet_mode <- match.arg(facet_mode)
  facet_col <- switch(
    facet_mode,
    site_state = "facet_site_state",
    state = "facet_state",
    site = "facet_site"
  )

  facet_levels <- names(sort(table(long_df[[facet_col]]), decreasing = TRUE))
  if (is.finite(max_facets) && max_facets > 0 && length(facet_levels) > max_facets) {
    facet_levels <- facet_levels[seq_len(max_facets)]
    long_df <- long_df[long_df[[facet_col]] %in% facet_levels, , drop = FALSE]
  }
  if (nrow(long_df) == 0) {
    stop("No rows remained after max_facets filtering.")
  }
  facet_levels <- names(sort(table(long_df[[facet_col]]), decreasing = TRUE))
  long_df[[facet_col]] <- factor(long_df[[facet_col]], levels = facet_levels)

  # Keep point table aligned with plotted site/state combinations.
  site_state_keep <- unique(data.frame(
    site_id = as.character(long_df$site_id),
    state_var = as.character(long_df$state_var),
    stringsAsFactors = FALSE
  ))
  site_state_keep$key <- paste(site_state_keep$site_id, site_state_keep$state_var, sep = "||")
  points_key <- paste(points_df$site_id, points_df$state_var, sep = "||")
  points_df <- points_df[points_key %in% site_state_keep$key, , drop = FALSE]
  if (nrow(points_df) == 0) {
    stop("No point rows remained after facet filtering.")
  }

  lookup_tbl <- unique(
    points_df[
      ,
      c("timestep", "block_id", "site_id", "state_idx", "state_var", "n_obs_mapped"),
      drop = FALSE
    ]
  )
  lookup_tbl$obs_available <- lookup_tbl$n_obs_mapped > 0
  lookup_tbl <- lookup_tbl[order(lookup_tbl$site_id, lookup_tbl$state_var, lookup_tbl$state_idx), , drop = FALSE]

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

  by_site_state_summary <- do.call(
    rbind,
    lapply(
      split(long_df, list(long_df$site_id, long_df$state_var, long_df$source), drop = TRUE),
      function(d) {
        sm <- summarize_values(d$value)
        out <- as.data.frame(as.list(sm), stringsAsFactors = FALSE)
        out$site_id <- as.character(d$site_id[1])
        out$state_var <- as.character(d$state_var[1])
        out$source <- as.character(d$source[1])
        out
      }
    )
  )
  by_site_state_summary <- by_site_state_summary[
    ,
    c("site_id", "state_var", "source", "n", "mean", "sd", "median", "p05", "p25", "p75", "p95")
  ]
  by_site_state_summary <- by_site_state_summary[
    order(by_site_state_summary$site_id, by_site_state_summary$state_var, by_site_state_summary$source),
    ,
    drop = FALSE
  ]

  points_csv <- file.path(output_dir, paste0(prefix, "_obs_forecast_analysis_points.csv"))
  lookup_csv <- file.path(output_dir, paste0(prefix, "_site_state_lookup.csv"))
  summary_overall_csv <- file.path(output_dir, paste0(prefix, "_summary_overall.csv"))
  summary_state_csv <- file.path(output_dir, paste0(prefix, "_summary_by_state.csv"))
  summary_site_state_csv <- file.path(output_dir, paste0(prefix, "_summary_by_site_state.csv"))
  utils::write.csv(points_df, points_csv, row.names = FALSE)
  utils::write.csv(lookup_tbl, lookup_csv, row.names = FALSE)
  utils::write.csv(overall_summary, summary_overall_csv, row.names = FALSE)
  utils::write.csv(by_state_summary, summary_state_csv, row.names = FALSE)
  utils::write.csv(by_site_state_summary, summary_site_state_csv, row.names = FALSE)

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

  n_facets <- length(facet_levels)
  ncol_facets <- if (n_facets <= 4) 2 else if (n_facets <= 9) 3 else 4
  nrow_facets <- max(1, ceiling(n_facets / ncol_facets))
  facet_height <- max(6, min(30, 2.2 * nrow_facets + 1.5))

  facet_subtitle <- switch(
    facet_mode,
    site_state = "Each panel is one site | state variable",
    state = "Each panel is one state variable (all sites combined)",
    site = "Each panel is one site (all state variables combined)"
  )
  facet_formula <- stats::as.formula(paste0("~", facet_col))

  p_by_facet <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = value, y = ggplot2::after_stat(density), color = source, fill = source)
  ) +
    ggplot2::geom_histogram(alpha = 0.28, bins = bins, position = "identity", linewidth = 0.2) +
    ggplot2::geom_density(alpha = 0.08, linewidth = 0.7, adjust = 1.05) +
    ggplot2::facet_wrap(facet_formula, scales = "free_y", ncol = ncol_facets) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = sprintf("Timestep %d: distribution by %s", timestep, facet_mode),
      subtitle = facet_subtitle,
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

  facet_png <- file.path(output_dir, paste0(prefix, "_density_by_", facet_mode, ".png"))
  ggplot2::ggsave(facet_png, p_by_facet, width = 14, height = facet_height, dpi = 180)

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

  plot_files <- c(overall_png, facet_png, ecdf_png)

  if (isTRUE(verbose)) {
    message(resolved$selection_note)
    message("Selected timestep: ", timestep)
    message("Loaded file: ", basename(sda_file))
    message("Rows extracted: ", nrow(points_df))
    message("Unique state variables plotted: ", length(unique(as.character(long_df$state_var))))
    message("Unique sites plotted: ", length(unique(as.character(long_df$site_id))))
    message("Output directory: ", normalizePath(output_dir, mustWork = FALSE))
    message("Outputs:")
    message("  - ", points_csv)
    message("  - ", lookup_csv)
    message("  - ", summary_overall_csv)
    message("  - ", summary_state_csv)
    message("  - ", summary_site_state_csv)
    for (pf in plot_files) {
      message("  - ", pf)
    }
  }

  invisible(list(
    timestep = timestep,
    source_file = sda_file,
    selection = resolved,
    point_table = points_df,
    lookup_table = lookup_tbl,
    long_table = long_df,
    summary_overall = overall_summary,
    summary_by_state = by_state_summary,
    summary_by_site_state = by_site_state_summary,
    csv_files = c(
      points = points_csv,
      lookup = lookup_csv,
      overall = summary_overall_csv,
      by_state = summary_state_csv,
      by_site_state = summary_site_state_csv
    ),
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
  max_facets <- parse_int(get_arg("max_facets", "0"), 0)

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
    site_ids = get_arg("site_ids", NULL),
    facet_mode = get_arg("facet_mode", "site_state"),
    max_states = max_states,
    max_facets = max_facets,
    verbose = TRUE
  )
}

# Example 1: directly specify timestep
# Rscript plot_timestep_state_distributions.R \
#   --outdir=/path/to/output_inter_q \
#   --timestep=25 \
#   --facet_mode=site_state \
#   --site_ids=US-Ha1,US-MMS \
#   --max_facets=36 \
#   --output_dir=/path/to/output_inter_q/timestep_distribution_plots
#
# Example 2: resolve timestep from bad-point rank/index
# Rscript plot_timestep_state_distributions.R \
#   --outdir=/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q_2 \
#   --index=6 \
#   --bad_points_csv=/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q_2/bad_assimilation_diagnostics/bad_assimilation_points.csv \
#   --index_column=rank
