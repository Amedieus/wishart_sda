extract_bad_assimilation_points <- function(outdir,
                                            output_dir = file.path(outdir, "bad_assimilation_diagnostics"),
                                            n_bad = 120,
                                            rhat_threshold = 1.1,
                                            ess_threshold = 200,
                                            write_csv = TRUE,
                                            write_plots = TRUE,
                                            verbose = TRUE) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  parse_timestep <- function(path) {
    as.integer(sub("^sda\\.output([0-9]+)\\.Rdata$", "\\1", basename(path)))
  }

  as_named_vector <- function(x, prefix = "state") {
    x <- suppressWarnings(as.numeric(x))
    nm <- names(x)
    if (is.null(nm)) {
      nm <- paste0(prefix, seq_along(x))
    } else {
      bad <- is.na(nm) | !nzchar(nm)
      nm[bad] <- paste0(prefix, which(bad))
    }
    names(x) <- nm
    x
  }

  invert_precision <- function(r_inv) {
    if (is.null(r_inv)) {
      return(NULL)
    }
    if (is.null(dim(r_inv))) {
      val <- suppressWarnings(as.numeric(r_inv))
      if (length(val) == 0 || !is.finite(val[1]) || val[1] == 0) {
        return(matrix(NA_real_, nrow = 1, ncol = 1))
      }
      return(matrix(1 / val[1], nrow = 1, ncol = 1))
    }
    r_inv <- as.matrix(r_inv)
    tryCatch(
      solve(r_inv),
      error = function(e) suppressWarnings(MASS::ginv(r_inv))
    )
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

  to_date_safe <- function(x) {
    if (is.null(x) || length(x) == 0) {
      return(as.Date(NA))
    }
    x1 <- x[[1]]
    d <- suppressWarnings(tryCatch(as.Date(x1), error = function(e) as.Date(NA)))
    if (!is.na(d)) {
      return(d)
    }
    if (is.character(x1) && nchar(x1) >= 10) {
      d2 <- suppressWarnings(tryCatch(as.Date(substr(x1, 1, 10)), error = function(e) as.Date(NA)))
      if (!is.na(d2)) {
        return(d2)
      }
    }
    if (is.numeric(x1) && is.finite(x1)) {
      d3 <- suppressWarnings(tryCatch(
        as.Date(as.POSIXct(x1, origin = "1970-01-01", tz = "UTC")),
        error = function(e) as.Date(NA)
      ))
      if (!is.na(d3)) {
        return(d3)
      }
    }
    as.Date(NA)
  }

  infer_date_label <- function(sda_outputs, t_idx) {
    d <- as.Date(NA)
    restart <- sda_outputs$restart.list
    if (!is.null(restart) && length(restart) > 0) {
      for (el in restart) {
        if (!is.list(el)) next
        cand <- to_date_safe(el$stop.time %||% el$stop_time)
        if (!is.na(cand)) {
          d <- cand
          break
        }
      }
    }
    if (is.na(d)) {
      return(paste0("timestep_", t_idx))
    }
    as.character(d)
  }

  is_converged <- function(block_diag_summary, rhat_threshold, ess_threshold) {
    if (is.null(block_diag_summary)) {
      return(NA)
    }
    rhat_max <- suppressWarnings(as.numeric(block_diag_summary$rhat_max))
    ess_min <- suppressWarnings(as.numeric(block_diag_summary$ess_min))
    if (length(rhat_max) == 0 || length(ess_min) == 0 ||
        !is.finite(rhat_max[1]) || !is.finite(ess_min[1])) {
      return(NA)
    }
    (rhat_max[1] <= rhat_threshold) && (ess_min[1] >= ess_threshold)
  }

  get_q_diag <- function(block, n_obs) {
    q_type <- suppressWarnings(as.integer(block$constant$q.type))
    if (!is.finite(q_type)) {
      return(rep(NA_real_, n_obs))
    }
    if (q_type == 3) {
      aq <- suppressWarnings(as.numeric(block$data$aq))
      bq <- suppressWarnings(as.numeric(block$data$bq))
      q <- aq / bq
      if (length(q) < n_obs) {
        q <- c(q, rep(NA_real_, n_obs - length(q)))
      }
      return(q[seq_len(n_obs)])
    }

    q_obj <- block$update$q_post_mean %||% block$data$aq
    if (is.null(q_obj)) {
      return(rep(NA_real_, n_obs))
    }
    if (is.null(dim(q_obj))) {
      q <- suppressWarnings(as.numeric(q_obj))
      if (length(q) < n_obs) {
        q <- c(q, rep(NA_real_, n_obs - length(q)))
      }
      return(q[seq_len(n_obs)])
    }

    q_mat <- as.matrix(q_obj)
    q <- diag(q_mat)
    if (length(q) < n_obs) {
      q <- c(q, rep(NA_real_, n_obs - length(q)))
    }
    q[seq_len(n_obs)]
  }

  get_r_diag <- function(block, n_obs) {
    r_cov <- invert_precision(block$data$r)
    if (is.null(r_cov)) {
      return(rep(NA_real_, n_obs))
    }
    r <- suppressWarnings(as.numeric(diag(r_cov)))
    if (length(r) < n_obs) {
      r <- c(r, rep(NA_real_, n_obs - length(r)))
    }
    r[seq_len(n_obs)]
  }

  robust_scale <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    med <- stats::median(x[is.finite(x)], na.rm = TRUE)
    if (!is.finite(med) || med <= 0) {
      med <- mean(x[is.finite(x)], na.rm = TRUE)
    }
    if (!is.finite(med) || med <= 0) {
      med <- 1
    }
    x / med
  }

  format_num <- function(x, digits = 4) {
    ifelse(
      is.finite(x),
      formatC(x, format = "fg", digits = digits),
      "NA"
    )
  }

  to_long_table <- function(df, columns, labels) {
    out <- vector("list", length(columns))
    for (i in seq_along(columns)) {
      col <- columns[i]
      out[[i]] <- data.frame(
        rank = df$rank,
        metric = labels[i],
        value = as.character(df[[col]]),
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, out)
  }

  sda_files <- list.files(
    outdir,
    pattern = "^sda\\.output[0-9]+\\.Rdata$",
    full.names = TRUE
  )
  if (length(sda_files) == 0) {
    stop("No files matched '^sda.output[0-9]+.Rdata$' in: ", normalizePath(outdir, mustWork = FALSE))
  }
  sda_files <- sda_files[order(parse_timestep(sda_files))]

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  rows_all <- list()

  for (f in sda_files) {
    t_idx <- parse_timestep(f)
    env <- new.env(parent = emptyenv())
    load(f, envir = env)

    if (!exists("sda.outputs", envir = env, inherits = FALSE)) {
      warning("Skipping file without sda.outputs: ", basename(f))
      next
    }
    sda_outputs <- get("sda.outputs", envir = env)
    enkf <- sda_outputs$enkf.params
    if (is.null(enkf) || is.null(enkf$block.list.all)) {
      warning("Skipping file with missing sda.outputs$enkf.params$block.list.all: ", basename(f))
      next
    }

    block_all <- enkf$block.list.all
    block_list <- NULL
    if (length(block_all) >= t_idx && !is.null(block_all[[t_idx]])) {
      block_list <- block_all[[t_idx]]
    } else {
      non_null <- which(!vapply(block_all, is.null, logical(1)))
      if (length(non_null) > 0) {
        block_list <- block_all[[max(non_null)]]
      }
    }
    if (is.null(block_list) || length(block_list) == 0) {
      next
    }

    date_label <- infer_date_label(sda_outputs, t_idx)

    for (b in seq_along(block_list)) {
      block <- block_list[[b]]
      y_obs <- suppressWarnings(as.numeric(block$data$y.censored))
      h_idx <- suppressWarnings(as.integer(block$constant$H))
      muf <- as_named_vector(block$data$muf %||% numeric(), prefix = "state")
      mufa <- as_named_vector(block$update$mufa %||% numeric(), prefix = "state")
      state_names <- names(muf)
      n_state <- length(muf)

      if (length(y_obs) == 0 || length(h_idx) == 0 || n_state == 0) {
        next
      }

      n_obs <- min(length(y_obs), length(h_idx))
      y_obs <- y_obs[seq_len(n_obs)]
      h_idx <- h_idx[seq_len(n_obs)]

      q_diag <- get_q_diag(block, n_obs)
      r_diag <- get_r_diag(block, n_obs)

      forecast_at_obs <- rep(NA_real_, n_obs)
      analysis_at_obs <- rep(NA_real_, n_obs)
      valid_state <- is.finite(h_idx) & h_idx > 0 & h_idx <= n_state
      forecast_at_obs[valid_state] <- suppressWarnings(as.numeric(muf[h_idx[valid_state]]))
      analysis_at_obs[valid_state] <- suppressWarnings(as.numeric(mufa[h_idx[valid_state]]))

      residual_forecast <- forecast_at_obs - y_obs
      residual_analysis <- analysis_at_obs - y_obs
      improvement_abs <- abs(residual_forecast) - abs(residual_analysis)

      site_ids <- as.character(block$site.ids %||% paste0("block_", b))
      site_map <- make_site_map(site_ids, n_state)
      state_var <- rep(NA_character_, n_obs)
      site_id <- rep(NA_character_, n_obs)
      state_var[valid_state] <- state_names[h_idx[valid_state]]
      site_id[valid_state] <- site_map[h_idx[valid_state]]

      diag_summary <- block$diag$summary
      conv <- is_converged(diag_summary, rhat_threshold = rhat_threshold, ess_threshold = ess_threshold)
      rhat_max <- suppressWarnings(as.numeric(diag_summary$rhat_max))[1]
      ess_min <- suppressWarnings(as.numeric(diag_summary$ess_min))[1]

      rows_all[[length(rows_all) + 1L]] <- data.frame(
        timestep = t_idx,
        date = date_label,
        block_id = b,
        obs_idx = seq_len(n_obs),
        state_idx = h_idx,
        state_var = state_var,
        site_id = site_id,
        obs_used = y_obs,
        forecast_at_obs = forecast_at_obs,
        analysis_at_obs = analysis_at_obs,
        residual_forecast = residual_forecast,
        residual_analysis = residual_analysis,
        improvement_abs = improvement_abs,
        q_diag = q_diag,
        r_diag = r_diag,
        rhat_max = ifelse(is.finite(rhat_max), rhat_max, NA_real_),
        ess_min = ifelse(is.finite(ess_min), ess_min, NA_real_),
        converged = ifelse(is.na(conv), NA, conv),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows_all) == 0) {
    stop("No observation-level diagnostics could be extracted from the provided SDA outputs.")
  }

  all_points <- do.call(rbind, rows_all)
  finite_idx <- is.finite(all_points$residual_analysis) & is.finite(all_points$improvement_abs)
  if (!any(finite_idx)) {
    stop("No finite residual/improvement values were found.")
  }

  all_points$abs_residual_analysis <- abs(all_points$residual_analysis)
  neg_improvement <- pmax(0, -all_points$improvement_abs)
  all_points$badness_score <- robust_scale(all_points$abs_residual_analysis) + robust_scale(neg_improvement)

  ord <- order(
    -all_points$badness_score,
    -all_points$abs_residual_analysis,
    all_points$improvement_abs,
    all_points$timestep
  )
  all_points <- all_points[ord, , drop = FALSE]

  n_bad <- max(1L, min(as.integer(n_bad), nrow(all_points)))
  bad_points <- all_points[seq_len(n_bad), , drop = FALSE]
  bad_points$rank <- seq_len(nrow(bad_points))

  table_csv <- file.path(output_dir, "bad_assimilation_points.csv")
  all_csv <- file.path(output_dir, "all_assimilation_points.csv")
  if (isTRUE(write_csv)) {
    utils::write.csv(bad_points, table_csv, row.names = FALSE)
    utils::write.csv(all_points, all_csv, row.names = FALSE)
  }

  plot_files <- character(0)
  if (isTRUE(write_plots)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' is required for plotting. CSV files were still written.")
    } else {
      plot_tbl <- data.frame(
        rank = bad_points$rank,
        timestep = bad_points$timestep,
        state_var = bad_points$state_var,
        site_id = bad_points$site_id,
        residual = format_num(bad_points$residual_analysis),
        improvement = format_num(bad_points$improvement_abs),
        q = format_num(bad_points$q_diag),
        r = format_num(bad_points$r_diag),
        converged = ifelse(is.na(bad_points$converged), "NA", ifelse(bad_points$converged, "yes", "no")),
        stringsAsFactors = FALSE
      )

      cols <- c("rank", "timestep", "state_var", "site_id", "residual", "improvement", "q", "r", "converged")
      labels <- c("rank", "timestep", "state var", "site id", "residual (a-o)", "improvement", "q", "r", "converged")
      long_tbl <- to_long_table(plot_tbl, columns = cols, labels = labels)
      long_tbl$metric <- factor(long_tbl$metric, levels = labels)
      long_tbl$rank <- factor(long_tbl$rank, levels = rev(sort(unique(long_tbl$rank))))

      p_table <- ggplot2::ggplot(long_tbl, ggplot2::aes(x = metric, y = rank)) +
        ggplot2::geom_tile(fill = "grey98", color = "grey88", linewidth = 0.25) +
        ggplot2::geom_text(ggplot2::aes(label = value), size = 3) +
        ggplot2::labs(
          title = sprintf("Top %d bad assimilation points", nrow(plot_tbl)),
          subtitle = "Residual = analysis - observation; Improvement = |forecast - observation| - |analysis - observation|",
          x = NULL,
          y = "badness rank (1 = worst)"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 28, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 8)
        )

      p_table_file <- file.path(output_dir, "bad_assimilation_points_table.png")
      ggplot2::ggsave(
        filename = p_table_file,
        plot = p_table,
        width = max(11, 0.9 * length(labels)),
        height = max(6, 0.30 * nrow(plot_tbl) + 2),
        dpi = 180
      )
      plot_files <- c(plot_files, p_table_file)

      # Companion plot to show global distribution and highlight bad points.
      all_points$converged_label <- ifelse(
        is.na(all_points$converged), "unknown",
        ifelse(all_points$converged, "yes", "no")
      )
      bad_points$key <- paste(bad_points$timestep, bad_points$block_id, bad_points$obs_idx, sep = "_")
      all_points$key <- paste(all_points$timestep, all_points$block_id, all_points$obs_idx, sep = "_")
      all_points$is_top_bad <- all_points$key %in% bad_points$key

      p_scatter <- ggplot2::ggplot(
        all_points,
        ggplot2::aes(x = improvement_abs, y = residual_analysis, color = converged_label)
      ) +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
        ggplot2::geom_point(alpha = 0.45, size = 1.4) +
        ggplot2::geom_point(
          data = all_points[all_points$is_top_bad, , drop = FALSE],
          shape = 21,
          fill = "red",
          color = "black",
          size = 2.6,
          stroke = 0.25
        ) +
        ggplot2::labs(
          title = "Residual vs improvement (top bad points highlighted)",
          x = "Improvement vs forecast (higher is better)",
          y = "Residual (analysis - observation)",
          color = "Converged"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

      p_scatter_file <- file.path(output_dir, "residual_vs_improvement_scatter.png")
      ggplot2::ggsave(
        filename = p_scatter_file,
        plot = p_scatter,
        width = 10.5,
        height = 6.5,
        dpi = 180
      )
      plot_files <- c(plot_files, p_scatter_file)
    }
  }

  if (isTRUE(verbose)) {
    message("Processed timesteps: ", length(unique(all_points$timestep)))
    message("Total point diagnostics: ", nrow(all_points))
    message("Top bad points kept: ", nrow(bad_points))
    message("CSV outputs: ", normalizePath(output_dir, mustWork = FALSE))
    if (length(plot_files) > 0) {
      message("Plot files:")
      for (pf in plot_files) message("  - ", pf)
    }
  }

  invisible(list(
    all_points = all_points,
    bad_points = bad_points,
    csv_files = c(all_assimilation = all_csv, bad_points = table_csv),
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
    if (is.null(val) || !nzchar(val)) return(default)
    val
  }

  parse_bool <- function(x, default = TRUE) {
    if (is.null(x) || !nzchar(x)) return(default)
    x <- tolower(trimws(x))
    if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
    if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
    default
  }

  outdir <- get_arg("outdir", NULL)
  if (is.null(outdir)) {
    stop("Missing required argument: --outdir=/path/to/sda/output")
  }

  extract_bad_assimilation_points(
    outdir = outdir,
    output_dir = get_arg("output_dir", file.path(outdir, "bad_assimilation_diagnostics")),
    n_bad = as.integer(get_arg("n_bad", "120")),
    rhat_threshold = as.numeric(get_arg("rhat_threshold", "1.1")),
    ess_threshold = as.numeric(get_arg("ess_threshold", "200")),
    write_csv = parse_bool(get_arg("write_csv", "true"), TRUE),
    write_plots = parse_bool(get_arg("write_plots", "true"), TRUE),
    verbose = TRUE
  )
}

# Example:
# Rscript plot_bad_assimilation_points.R \
#   --outdir=/path/to/output_inter_q \
#   --output_dir=/path/to/output_inter_q/bad_assimilation_diagnostics \
#   --n_bad=100 \
#   --rhat_threshold=1.1 \
#   --ess_threshold=200
