extract_q_r_xmod <- function(outdir,
                             output_dir = file.path(outdir, "q_r_xmod_tables"),
                             write_csv = TRUE,
                             verbose = TRUE,
                             write_plots = TRUE,
                             plot_dir = file.path(output_dir, "plots"),
                             plot_width = 14,
                             plot_height = 8,
                             plot_dpi = 160) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  parse_timestep <- function(path) {
    as.integer(sub("^sda\\.output([0-9]+)\\.Rdata$", "\\1", basename(path)))
  }

  as_named_vector <- function(x, prefix = "state") {
    x <- as.numeric(x)
    nm <- names(x)
    if (is.null(nm)) {
      nm <- paste0(prefix, "_", seq_along(x))
    } else {
      empty <- is.na(nm) | nm == ""
      nm[empty] <- paste0(prefix, "_", which(empty))
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
      if (!is.finite(val) || val == 0) {
        return(matrix(NA_real_, nrow = 1, ncol = 1))
      }
      return(matrix(1 / val, nrow = 1, ncol = 1))
    }

    if (nrow(r_inv) == 1 && ncol(r_inv) == 1) {
      val <- suppressWarnings(as.numeric(r_inv[1, 1]))
      if (!is.finite(val) || val == 0) {
        return(matrix(NA_real_, nrow = 1, ncol = 1))
      }
      return(matrix(1 / val, nrow = 1, ncol = 1))
    }

    out <- tryCatch(
      solve(r_inv),
      error = function(e) {
        suppressWarnings(MASS::ginv(r_inv))
      }
    )
    out
  }

  make_site_map <- function(site_ids, n_state) {
    if (length(site_ids) == 0 || n_state == 0) {
      return(rep(NA_character_, n_state))
    }
    if (n_state %% length(site_ids) == 0) {
      return(rep(site_ids, each = n_state / length(site_ids)))
    }
    rep(site_ids, length.out = n_state)
  }

  build_meta_for_obs <- function(obs_idx, site_map, state_names, n_obs, fallback_site) {
    idx <- suppressWarnings(as.integer(obs_idx))
    idx <- idx[is.finite(idx) & idx > 0 & idx <= length(state_names)]

    if (length(idx) >= n_obs) {
      idx <- idx[seq_len(n_obs)]
      return(data.frame(
        obs_idx = seq_len(n_obs),
        state_local_idx = idx,
        site_id = site_map[idx],
        state_var = state_names[idx],
        stringsAsFactors = FALSE
      ))
    }

    state_idx <- seq_len(n_obs)
    state_idx[state_idx > length(state_names)] <- NA_integer_
    data.frame(
      obs_idx = seq_len(n_obs),
      state_local_idx = state_idx,
      site_id = rep(fallback_site, n_obs),
      state_var = ifelse(
        is.na(state_idx),
        paste0("obs_", seq_len(n_obs)),
        state_names[state_idx]
      ),
      stringsAsFactors = FALSE
    )
  }

  matrix_to_long <- function(mat,
                             row_meta,
                             col_meta = row_meta,
                             value_name = "value") {
    if (is.null(mat) || is.null(dim(mat))) {
      return(data.frame())
    }
    nr <- nrow(mat)
    nc <- ncol(mat)
    if (nrow(row_meta) != nr || nrow(col_meta) != nc) {
      return(data.frame())
    }

    out <- vector("list", nr * nc)
    k <- 1L
    for (i in seq_len(nr)) {
      for (j in seq_len(nc)) {
        out[[k]] <- data.frame(
          row_obs_idx = i,
          col_obs_idx = j,
          row_site_id = row_meta$site_id[i],
          col_site_id = col_meta$site_id[j],
          row_state_var = row_meta$state_var[i],
          col_state_var = col_meta$state_var[j],
          stringsAsFactors = FALSE
        )
        out[[k]][[value_name]] <- mat[i, j]
        k <- k + 1L
      }
    }
    do.call(rbind, out)
  }

  safe_bind <- function(x) {
    x <- Filter(function(obj) !is.null(obj) && nrow(obj) > 0, x)
    if (length(x) == 0) {
      return(data.frame())
    }
    do.call(rbind, x)
  }

  safe_order <- function(df, cols) {
    if (nrow(df) == 0 || length(cols) == 0) {
      return(df)
    }
    cols <- cols[cols %in% names(df)]
    if (length(cols) == 0) {
      return(df)
    }
    ord <- do.call(order, df[cols])
    df[ord, , drop = FALSE]
  }

  as_char_scalar <- function(x, fallback) {
    if (length(x) == 0) {
      return(fallback)
    }
    y <- as.character(x[1])
    if (!is.na(y) && nzchar(y)) y else fallback
  }

  to_date_safe <- function(x) {
    if (is.null(x) || length(x) == 0) {
      return(as.Date(NA))
    }
    x1 <- x[[1]]
    d <- tryCatch(as.Date(x1), error = function(e) as.Date(NA))
    if (!is.na(d)) {
      return(d)
    }
    if (is.character(x1) && nchar(x1) >= 10) {
      d2 <- tryCatch(as.Date(substr(x1, 1, 10)), error = function(e) as.Date(NA))
      if (!is.na(d2)) {
        return(d2)
      }
    }
    if (is.numeric(x1) && is.finite(x1)) {
      d3 <- tryCatch(
        as.Date(as.POSIXct(x1, origin = "1970-01-01", tz = "UTC")),
        error = function(e) as.Date(NA)
      )
      if (!is.na(d3)) {
        return(d3)
      }
    }
    as.Date(NA)
  }

  infer_date_year <- function(sda_outputs, t_idx) {
    d <- as.Date(NA)

    restart <- sda_outputs$restart.list
    if (!is.null(restart) && length(restart) > 0) {
      for (el in restart) {
        if (!is.list(el)) {
          next
        }
        cand <- to_date_safe(el$stop.time %||% el$stop_time)
        if (!is.na(cand)) {
          d <- cand
          break
        }
      }
    }

    if (is.na(d)) {
      enkf <- sda_outputs$enkf.params
      rst2 <- enkf$RestartList
      if (!is.null(rst2) && length(rst2) > 0) {
        for (el in rst2) {
          if (!is.list(el)) {
            next
          }
          cand <- to_date_safe(el$stop.time %||% el$stop_time)
          if (!is.na(cand)) {
            d <- cand
            break
          }
        }
      }
    }

    date_label <- if (!is.na(d)) as.character(d) else paste0("timestep_", t_idx)
    year_value <- if (!is.na(d)) as.integer(format(d, "%Y")) else NA_integer_
    list(date_label = date_label, year = year_value)
  }

  lookup_obs_raw <- function(obs_t, site_id, state_var) {
    if (is.null(obs_t) || is.na(site_id) || is.na(state_var)) {
      return(NA_real_)
    }
    if (!is.list(obs_t) || !(site_id %in% names(obs_t))) {
      return(NA_real_)
    }
    obs_site <- obs_t[[site_id]]
    if (is.null(obs_site) || length(obs_site) == 0) {
      return(NA_real_)
    }
    if (is.null(names(obs_site)) || !(state_var %in% names(obs_site))) {
      return(NA_real_)
    }
    val <- suppressWarnings(as.numeric(obs_site[[state_var]]))
    if (length(val) == 0) NA_real_ else val[1]
  }

  generate_residual_relation_plots <- function(df,
                                               plot_dir,
                                               plot_width,
                                               plot_height,
                                               plot_dpi,
                                               verbose = TRUE) {
    if (nrow(df) == 0) {
      return(character(0))
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' not available. Skipping plot generation.")
      return(character(0))
    }

    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    }

    safe_write_plot <- function(p, file_name) {
      out_path <- file.path(plot_dir, file_name)
      ok <- tryCatch(
        {
          ggplot2::ggsave(
            filename = out_path,
            plot = p,
            width = plot_width,
            height = plot_height,
            dpi = plot_dpi
          )
          TRUE
        },
        error = function(e) {
          warning("Failed to write plot ", out_path, ": ", conditionMessage(e))
          FALSE
        }
      )
      if (ok) out_path else NA_character_
    }

    build_pair_time_plot <- function(value_col, title_text, y_label_other) {
      keep <- c("timestep", "date", "site_id", "state_var", "residual_xmod_minus_obs_scaled", value_col)
      keep <- keep[keep %in% names(df)]
      if (!all(c("timestep", "site_id", "state_var", "residual_xmod_minus_obs_scaled", value_col) %in% keep)) {
        return(NULL)
      }
      d <- df[, keep, drop = FALSE]
      d <- d[
        is.finite(d$residual_xmod_minus_obs_scaled) & is.finite(d[[value_col]]),
        ,
        drop = FALSE
      ]
      if (nrow(d) == 0) {
        return(NULL)
      }

      long_df <- rbind(
        data.frame(
          timestep = d$timestep,
          date = d$date,
          site_id = d$site_id,
          state_var = d$state_var,
          metric = y_label_other,
          value = d[[value_col]],
          stringsAsFactors = FALSE
        ),
        data.frame(
          timestep = d$timestep,
          date = d$date,
          site_id = d$site_id,
          state_var = d$state_var,
          metric = "Residual (X.mod - Obs)",
          value = d$residual_xmod_minus_obs_scaled,
          stringsAsFactors = FALSE
        )
      )

      ggplot2::ggplot(
        long_df,
        ggplot2::aes(x = timestep, y = value, color = site_id, group = site_id)
      ) +
        ggplot2::geom_line(alpha = 0.7) +
        ggplot2::geom_point(alpha = 0.8, size = 1.4) +
        ggplot2::facet_grid(metric ~ state_var, scales = "free_y") +
        ggplot2::labs(
          title = title_text,
          x = "Timestep",
          y = "Value",
          color = "Site"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    }

    build_pair_scatter_plot <- function(value_col, x_label, title_text) {
      keep <- c("timestep", "site_id", "state_var", "residual_xmod_minus_obs_scaled", value_col)
      keep <- keep[keep %in% names(df)]
      if (!all(c("timestep", "site_id", "state_var", "residual_xmod_minus_obs_scaled", value_col) %in% keep)) {
        return(NULL)
      }
      d <- df[, keep, drop = FALSE]
      d <- d[
        is.finite(d$residual_xmod_minus_obs_scaled) & is.finite(d[[value_col]]),
        ,
        drop = FALSE
      ]
      if (nrow(d) == 0) {
        return(NULL)
      }

      p <- ggplot2::ggplot(
        d,
        ggplot2::aes_string(
          x = value_col,
          y = "residual_xmod_minus_obs_scaled",
          color = "timestep"
        )
      ) +
        ggplot2::geom_point(alpha = 0.75, size = 1.6) +
        ggplot2::facet_wrap(~state_var, scales = "free") +
        ggplot2::labs(
          title = title_text,
          x = x_label,
          y = "Residual (X.mod - Obs, scaled)",
          color = "Timestep"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

      if (nrow(d) >= 3) {
        p <- p + ggplot2::geom_smooth(
          method = "lm",
          se = FALSE,
          color = "black",
          linewidth = 0.6,
          inherit.aes = FALSE,
          mapping = ggplot2::aes_string(
            x = value_col,
            y = "residual_xmod_minus_obs_scaled"
          )
        )
      }
      p
    }

    files <- character(0)

    p_q_time <- build_pair_time_plot(
      value_col = "q_diag",
      title_text = "Q and residual over time by state variable",
      y_label_other = "Q (diag / representative)"
    )
    if (!is.null(p_q_time)) {
      files <- c(files, safe_write_plot(p_q_time, "q_residual_over_time_by_state.png"))
    }

    p_r_time <- build_pair_time_plot(
      value_col = "r_diag",
      title_text = "R and residual over time by state variable",
      y_label_other = "R (diag)"
    )
    if (!is.null(p_r_time)) {
      files <- c(files, safe_write_plot(p_r_time, "r_residual_over_time_by_state.png"))
    }

    p_q_scatter <- build_pair_scatter_plot(
      value_col = "q_diag",
      x_label = "Q (diag / representative)",
      title_text = "Residual vs Q by state variable (colored by timestep)"
    )
    if (!is.null(p_q_scatter)) {
      files <- c(files, safe_write_plot(p_q_scatter, "residual_vs_q_scatter_by_state.png"))
    }

    p_r_scatter <- build_pair_scatter_plot(
      value_col = "r_diag",
      x_label = "R (diag)",
      title_text = "Residual vs R by state variable (colored by timestep)"
    )
    if (!is.null(p_r_scatter)) {
      files <- c(files, safe_write_plot(p_r_scatter, "residual_vs_r_scatter_by_state.png"))
    }

    files <- files[!is.na(files)]
    if (verbose && length(files) > 0) {
      message("Saved plots to: ", normalizePath(plot_dir, mustWork = FALSE))
    }
    files
  }

  sda_files <- list.files(
    outdir,
    pattern = "^sda\\.output[0-9]+\\.Rdata$",
    full.names = TRUE
  )

  if (length(sda_files) == 0) {
    stop(
      paste0(
        "No files matched '^sda.output[0-9]+.Rdata$' in: ",
        normalizePath(outdir, mustWork = FALSE)
      )
    )
  }

  sda_files <- sda_files[order(parse_timestep(sda_files))]

  if (write_csv && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  q_diag_all <- list()
  q_matrix_all <- list()
  r_diag_all <- list()
  r_matrix_all <- list()
  xmod_all <- list()
  residual_all <- list()

  for (f in sda_files) {
    t_idx <- parse_timestep(f)
    env <- new.env(parent = emptyenv())
    load(f, envir = env)

    if (!exists("sda.outputs", envir = env, inherits = FALSE)) {
      warning("Skipping file without sda.outputs object: ", basename(f))
      next
    }
    sda_outputs <- get("sda.outputs", envir = env)
    enkf <- sda_outputs$enkf.params

    if (is.null(enkf)) {
      warning("Skipping file with NULL sda.outputs$enkf.params: ", basename(f))
      next
    }

    block_all <- enkf$block.list.all
    block_list <- NULL
    if (!is.null(block_all) && length(block_all) >= t_idx && !is.null(block_all[[t_idx]])) {
      block_list <- block_all[[t_idx]]
    } else if (!is.null(block_all)) {
      non_null <- which(!vapply(block_all, is.null, logical(1)))
      if (length(non_null) > 0) {
        block_list <- block_all[[max(non_null)]]
      }
    }

    if (is.null(block_list)) {
      warning("No block.list found for timestep ", t_idx, " in ", basename(f))
      next
    }

    date_info <- infer_date_year(sda_outputs, t_idx)
    date_label <- date_info$date_label
    year_value <- date_info$year
    obs_raw_t <- sda_outputs$obs.mean

    for (b in seq_along(block_list)) {
      block <- block_list[[b]]

      site_ids <- as.character(block$site.ids %||% paste0("block_", b))
      fallback_site <- as_char_scalar(site_ids, paste0("block_", b))

      muf <- as_named_vector(block$data$muf %||% numeric(), prefix = "state")
      state_names <- names(muf)
      n_state <- length(muf)
      site_map <- make_site_map(site_ids, n_state)

      xmod <- as_named_vector(block$update$mufa %||% rep(NA_real_, n_state), prefix = "state")
      if (length(xmod) != n_state) {
        xmod <- rep(NA_real_, n_state)
        names(xmod) <- state_names
      }

      xmod_df <- data.frame(
        timestep = t_idx,
        date = date_label,
        year = year_value,
        block_id = b,
        site_id = site_map,
        state_local_idx = seq_len(n_state),
        state_var = state_names,
        xmod_mean = as.numeric(xmod),
        stringsAsFactors = FALSE
      )
      xmod_all[[length(xmod_all) + 1L]] <- xmod_df

      obs_idx <- block$constant$H %||% integer(0)
      q_type <- block$constant$q.type %||% NA_integer_

      if (q_type == 3) {
        aq <- suppressWarnings(as.numeric(block$data$aq))
        bq <- suppressWarnings(as.numeric(block$data$bq))
        n_obs <- max(length(aq), length(bq), length(obs_idx), 0L)
        if (n_obs > 0) {
          meta <- build_meta_for_obs(
            obs_idx = obs_idx,
            site_map = site_map,
            state_names = state_names,
            n_obs = n_obs,
            fallback_site = fallback_site
          )
          aq_pad <- rep(NA_real_, n_obs)
          bq_pad <- rep(NA_real_, n_obs)
          aq_pad[seq_along(aq)] <- aq
          bq_pad[seq_along(bq)] <- bq

          q_mean <- aq_pad / bq_pad
          q_mode <- ifelse(aq_pad > 1 & is.finite(bq_pad) & bq_pad > 0, (aq_pad - 1) / bq_pad, NA_real_)

          q_diag_df <- data.frame(
            timestep = t_idx,
            date = date_label,
            year = year_value,
            block_id = b,
            q_type = q_type,
            obs_idx = seq_len(n_obs),
            site_id = meta$site_id,
            state_var = meta$state_var,
            state_local_idx = meta$state_local_idx,
            aq = aq_pad,
            bq = bq_pad,
            q_diag = q_mean,
            q_mode = q_mode,
            stringsAsFactors = FALSE
          )
          q_diag_all[[length(q_diag_all) + 1L]] <- q_diag_df
        }
      } else if (q_type == 4) {
        q_mat <- block$data$aq
        if (!is.null(q_mat)) {
          if (is.null(dim(q_mat))) {
            q_mat <- matrix(as.numeric(q_mat), nrow = 1, ncol = 1)
          }
          n_obs <- nrow(q_mat)
          meta <- build_meta_for_obs(
            obs_idx = obs_idx,
            site_map = site_map,
            state_names = state_names,
            n_obs = n_obs,
            fallback_site = fallback_site
          )

          q_diag_df <- data.frame(
            timestep = t_idx,
            date = date_label,
            year = year_value,
            block_id = b,
            q_type = q_type,
            obs_idx = seq_len(n_obs),
            site_id = meta$site_id,
            state_var = meta$state_var,
            state_local_idx = meta$state_local_idx,
            aq = diag(q_mat),
            bq = rep_len(as.numeric(block$data$bq), n_obs),
            q_diag = diag(q_mat),
            q_mode = NA_real_,
            stringsAsFactors = FALSE
          )
          q_diag_all[[length(q_diag_all) + 1L]] <- q_diag_df

          q_mat_long <- matrix_to_long(q_mat, meta, value_name = "q_value")
          if (nrow(q_mat_long) > 0) {
            q_mat_long$timestep <- t_idx
            q_mat_long$date <- date_label
            q_mat_long$year <- year_value
            q_mat_long$block_id <- b
            q_mat_long$q_type <- q_type
            q_mat_long$bq <- as.numeric(block$data$bq)[1]
            q_matrix_all[[length(q_matrix_all) + 1L]] <- q_mat_long
          }
        }
      }

      r_mat <- invert_precision(block$data$r)
      if (!is.null(r_mat)) {
        n_obs <- nrow(r_mat)
        meta <- build_meta_for_obs(
          obs_idx = obs_idx,
          site_map = site_map,
          state_names = state_names,
          n_obs = n_obs,
          fallback_site = fallback_site
        )

        r_diag_df <- data.frame(
          timestep = t_idx,
          date = date_label,
          year = year_value,
          block_id = b,
          obs_idx = seq_len(n_obs),
          site_id = meta$site_id,
          state_var = meta$state_var,
          state_local_idx = meta$state_local_idx,
          r_diag = diag(r_mat),
          stringsAsFactors = FALSE
        )
        r_diag_all[[length(r_diag_all) + 1L]] <- r_diag_df

        r_mat_long <- matrix_to_long(r_mat, meta, value_name = "r_value")
        if (nrow(r_mat_long) > 0) {
          r_mat_long$timestep <- t_idx
          r_mat_long$date <- date_label
          r_mat_long$year <- year_value
          r_mat_long$block_id <- b
          r_matrix_all[[length(r_matrix_all) + 1L]] <- r_mat_long
        }
      }

      y_obs <- suppressWarnings(as.numeric(block$data$y.censored))
      n_obs_res <- length(y_obs)
      if (n_obs_res > 0) {
        meta_res <- build_meta_for_obs(
          obs_idx = obs_idx,
          site_map = site_map,
          state_names = state_names,
          n_obs = n_obs_res,
          fallback_site = fallback_site
        )

        xmod_at_obs <- rep(NA_real_, n_obs_res)
        valid_idx <- !is.na(meta_res$state_local_idx) &
          meta_res$state_local_idx > 0 &
          meta_res$state_local_idx <= length(xmod)
        xmod_at_obs[valid_idx] <- as.numeric(xmod[meta_res$state_local_idx[valid_idx]])

        residual_scaled <- xmod_at_obs - y_obs

        obs_raw <- vapply(
          seq_len(n_obs_res),
          function(i) lookup_obs_raw(obs_raw_t, meta_res$site_id[i], meta_res$state_var[i]),
          numeric(1)
        )

        residual_df <- data.frame(
          timestep = t_idx,
          date = date_label,
          year = year_value,
          block_id = b,
          q_type = q_type,
          obs_idx = seq_len(n_obs_res),
          site_id = meta_res$site_id,
          state_var = meta_res$state_var,
          state_local_idx = meta_res$state_local_idx,
          xmod_mean = xmod_at_obs,
          obs_used_scaled = y_obs,
          residual_xmod_minus_obs_scaled = residual_scaled,
          obs_raw = obs_raw,
          stringsAsFactors = FALSE
        )
        residual_all[[length(residual_all) + 1L]] <- residual_df
      }
    }
  }

  q_diag <- safe_bind(q_diag_all)
  q_matrix <- safe_bind(q_matrix_all)
  r_diag <- safe_bind(r_diag_all)
  r_matrix <- safe_bind(r_matrix_all)
  xmod <- safe_bind(xmod_all)
  residual <- safe_bind(residual_all)

  residual_qr <- residual
  join_keys <- c("timestep", "block_id", "obs_idx", "site_id", "state_var")

  if (nrow(residual_qr) > 0 && nrow(q_diag) > 0 && all(join_keys %in% names(q_diag))) {
    q_keep <- unique(c(join_keys, "q_diag", "q_mode"))
    q_keep <- q_keep[q_keep %in% names(q_diag)]
    residual_qr <- merge(
      residual_qr,
      q_diag[, q_keep, drop = FALSE],
      by = join_keys,
      all.x = TRUE
    )
    if ("date.x" %in% names(residual_qr)) {
      residual_qr$date <- residual_qr$date.x
      residual_qr$date.x <- NULL
    } else if ("date.y" %in% names(residual_qr) && !("date" %in% names(residual_qr))) {
      residual_qr$date <- residual_qr$date.y
      residual_qr$date.y <- NULL
    }
    if ("year.x" %in% names(residual_qr)) {
      residual_qr$year <- residual_qr$year.x
      residual_qr$year.x <- NULL
    } else if ("year.y" %in% names(residual_qr) && !("year" %in% names(residual_qr))) {
      residual_qr$year <- residual_qr$year.y
      residual_qr$year.y <- NULL
    }
  } else if (nrow(residual_qr) > 0) {
    residual_qr$q_diag <- NA_real_
    residual_qr$q_mode <- NA_real_
  }

  if (nrow(residual_qr) > 0 && nrow(r_diag) > 0 && all(join_keys %in% names(r_diag))) {
    r_keep <- unique(c(join_keys, "r_diag"))
    r_keep <- r_keep[r_keep %in% names(r_diag)]
    residual_qr <- merge(
      residual_qr,
      r_diag[, r_keep, drop = FALSE],
      by = join_keys,
      all.x = TRUE
    )
    if ("date.x" %in% names(residual_qr)) {
      residual_qr$date <- residual_qr$date.x
      residual_qr$date.x <- NULL
    } else if ("date.y" %in% names(residual_qr) && !("date" %in% names(residual_qr))) {
      residual_qr$date <- residual_qr$date.y
      residual_qr$date.y <- NULL
    }
    if ("year.x" %in% names(residual_qr)) {
      residual_qr$year <- residual_qr$year.x
      residual_qr$year.x <- NULL
    } else if ("year.y" %in% names(residual_qr) && !("year" %in% names(residual_qr))) {
      residual_qr$year <- residual_qr$year.y
      residual_qr$year.y <- NULL
    }
  } else if (nrow(residual_qr) > 0 && !("r_diag" %in% names(residual_qr))) {
    residual_qr$r_diag <- NA_real_
  }

  q_xmod <- data.frame()
  if (all(c("timestep", "block_id", "site_id", "state_var", "q_diag") %in% names(q_diag)) &&
      all(c("timestep", "block_id", "site_id", "state_var", "xmod_mean") %in% names(xmod))) {
    q_xmod <- merge(
      q_diag[, c("timestep", "date", "year", "block_id", "site_id", "state_var", "q_diag"), drop = FALSE],
      xmod[, c("timestep", "date", "year", "block_id", "site_id", "state_var", "xmod_mean"), drop = FALSE],
      by = c("timestep", "block_id", "site_id", "state_var"),
      all = FALSE
    )
    if ("date.x" %in% names(q_xmod)) {
      q_xmod$date <- q_xmod$date.x
      q_xmod$date.x <- NULL
    } else if ("date.y" %in% names(q_xmod)) {
      q_xmod$date <- q_xmod$date.y
      q_xmod$date.y <- NULL
    }
    if ("year.x" %in% names(q_xmod)) {
      q_xmod$year <- q_xmod$year.x
      q_xmod$year.x <- NULL
    } else if ("year.y" %in% names(q_xmod)) {
      q_xmod$year <- q_xmod$year.y
      q_xmod$year.y <- NULL
    }
  }

  r_xmod <- data.frame()
  if (all(c("timestep", "block_id", "site_id", "state_var", "r_diag") %in% names(r_diag)) &&
      all(c("timestep", "block_id", "site_id", "state_var", "xmod_mean") %in% names(xmod))) {
    r_xmod <- merge(
      r_diag[, c("timestep", "date", "year", "block_id", "site_id", "state_var", "r_diag"), drop = FALSE],
      xmod[, c("timestep", "date", "year", "block_id", "site_id", "state_var", "xmod_mean"), drop = FALSE],
      by = c("timestep", "block_id", "site_id", "state_var"),
      all = FALSE
    )
    if ("date.x" %in% names(r_xmod)) {
      r_xmod$date <- r_xmod$date.x
      r_xmod$date.x <- NULL
    } else if ("date.y" %in% names(r_xmod)) {
      r_xmod$date <- r_xmod$date.y
      r_xmod$date.y <- NULL
    }
    if ("year.x" %in% names(r_xmod)) {
      r_xmod$year <- r_xmod$year.x
      r_xmod$year.x <- NULL
    } else if ("year.y" %in% names(r_xmod)) {
      r_xmod$year <- r_xmod$year.y
      r_xmod$year.y <- NULL
    }
  }

  summarize_relationship <- function(df, value_col) {
    if (nrow(df) == 0) {
      return(data.frame())
    }
    split_key <- paste(df$site_id, df$state_var, sep = "||")
    chunks <- split(df, split_key)
    out <- lapply(chunks, function(d) {
      vals <- d[[value_col]]
      x <- d$xmod_mean
      ok <- is.finite(vals) & is.finite(x)
      n <- sum(ok)
      corr <- if (n >= 2) stats::cor(vals[ok], x[ok]) else NA_real_
      slope <- if (n >= 2) {
        coef(stats::lm(x[ok] ~ vals[ok]))[2]
      } else {
        NA_real_
      }
      data.frame(
        site_id = d$site_id[1],
        state_var = d$state_var[1],
        n_timestep = n,
        correlation = corr,
        slope_xmod_on_value = slope,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, out)
  }

  q_xmod_summary <- summarize_relationship(q_xmod, "q_diag")
  r_xmod_summary <- summarize_relationship(r_xmod, "r_diag")

  summarize_residual_qr <- function(df, group_cols) {
    if (nrow(df) == 0) {
      return(data.frame())
    }
    group_cols <- group_cols[group_cols %in% names(df)]
    if (length(group_cols) == 0) {
      return(data.frame())
    }
    split_key <- apply(df[, group_cols, drop = FALSE], 1, paste, collapse = "||")
    chunks <- split(df, split_key)
    out <- lapply(chunks, function(d) {
      base <- d[1, group_cols, drop = FALSE]

      res <- d$residual_xmod_minus_obs_scaled
      qv <- d$q_diag
      rv <- d$r_diag

      ok_res <- is.finite(res)
      ok_q <- ok_res & is.finite(qv)
      ok_r <- ok_res & is.finite(rv)

      n_res <- sum(ok_res)
      n_q <- sum(ok_q)
      n_r <- sum(ok_r)

      corr_r_q <- if (n_q >= 2) stats::cor(res[ok_q], qv[ok_q]) else NA_real_
      corr_r_r <- if (n_r >= 2) stats::cor(res[ok_r], rv[ok_r]) else NA_real_

      slope_r_on_q <- if (n_q >= 2) coef(stats::lm(res[ok_q] ~ qv[ok_q]))[2] else NA_real_
      slope_r_on_r <- if (n_r >= 2) coef(stats::lm(res[ok_r] ~ rv[ok_r]))[2] else NA_real_

      data.frame(
        base,
        n_rows = nrow(d),
        n_timestep = length(unique(d$timestep)),
        n_residual = n_res,
        residual_mean = if (n_res > 0) mean(res[ok_res]) else NA_real_,
        residual_sd = if (n_res > 1) stats::sd(res[ok_res]) else NA_real_,
        q_mean = if (n_q > 0) mean(qv[ok_q]) else NA_real_,
        r_mean = if (n_r > 0) mean(rv[ok_r]) else NA_real_,
        n_pair_res_q = n_q,
        corr_residual_q = corr_r_q,
        slope_residual_on_q = slope_r_on_q,
        n_pair_res_r = n_r,
        corr_residual_r = corr_r_r,
        slope_residual_on_r = slope_r_on_r,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, out)
  }

  residual_qr_by_site_year_state <- summarize_residual_qr(
    residual_qr,
    c("site_id", "year", "state_var")
  )
  residual_qr_by_site_year <- summarize_residual_qr(
    residual_qr,
    c("site_id", "year")
  )

  results <- list(
    q_diag = safe_order(q_diag, c("timestep", "site_id", "state_var")),
    q_matrix = safe_order(q_matrix, c("timestep", "block_id", "row_obs_idx", "col_obs_idx")),
    r_diag = safe_order(r_diag, c("timestep", "site_id", "state_var")),
    r_matrix = safe_order(r_matrix, c("timestep", "block_id", "row_obs_idx", "col_obs_idx")),
    xmod = safe_order(xmod, c("timestep", "site_id", "state_var")),
    residual = safe_order(residual, c("timestep", "site_id", "state_var")),
    residual_qr = safe_order(residual_qr, c("timestep", "site_id", "state_var")),
    residual_qr_by_site_year_state = safe_order(
      residual_qr_by_site_year_state,
      c("site_id", "year", "state_var")
    ),
    residual_qr_by_site_year = safe_order(
      residual_qr_by_site_year,
      c("site_id", "year")
    ),
    q_xmod = safe_order(q_xmod, c("timestep", "site_id", "state_var")),
    r_xmod = safe_order(r_xmod, c("timestep", "site_id", "state_var")),
    q_xmod_summary = safe_order(q_xmod_summary, c("site_id", "state_var")),
    r_xmod_summary = safe_order(r_xmod_summary, c("site_id", "state_var"))
  )

  if (write_csv) {
    utils::write.csv(results$q_diag, file.path(output_dir, "q_diag_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$q_matrix, file.path(output_dir, "q_matrix_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$r_diag, file.path(output_dir, "r_diag_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$r_matrix, file.path(output_dir, "r_matrix_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$xmod, file.path(output_dir, "xmod_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$residual, file.path(output_dir, "residual_xmod_minus_obs_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$residual_qr, file.path(output_dir, "residual_with_qr_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(
      results$residual_qr_by_site_year_state,
      file.path(output_dir, "residual_qr_relationship_by_site_year_state.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      results$residual_qr_by_site_year,
      file.path(output_dir, "residual_qr_relationship_by_site_year.csv"),
      row.names = FALSE
    )
    utils::write.csv(results$q_xmod, file.path(output_dir, "q_vs_xmod_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$r_xmod, file.path(output_dir, "r_vs_xmod_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$q_xmod_summary, file.path(output_dir, "q_vs_xmod_summary.csv"), row.names = FALSE)
    utils::write.csv(results$r_xmod_summary, file.path(output_dir, "r_vs_xmod_summary.csv"), row.names = FALSE)

    if (verbose) {
      message("Saved tables to: ", normalizePath(output_dir, mustWork = FALSE))
    }
  }

  plot_files <- character(0)
  if (isTRUE(write_plots)) {
    plot_files <- generate_residual_relation_plots(
      df = results$residual_qr,
      plot_dir = plot_dir,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_dpi = plot_dpi,
      verbose = verbose
    )
  }
  results$plot_files <- plot_files

  if (verbose) {
    timesteps <- if ("timestep" %in% names(results$xmod)) {
      length(unique(results$xmod$timestep))
    } else {
      0
    }
    message("Timesteps processed: ", timesteps)
    message("Rows - Q diag: ", nrow(results$q_diag),
            ", R diag: ", nrow(results$r_diag),
            ", X.mod: ", nrow(results$xmod),
            ", Residual: ", nrow(results$residual))
    message("Plot files: ", length(results$plot_files))
  }

  results
}

extract_block_mcmc_param_diag <- function(outdir,
                                          file_id = NA_integer_,
                                          output_csv = NULL,
                                          write_csv = FALSE,
                                          verbose = TRUE) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  parse_timestep <- function(path) {
    as.integer(sub("^sda\\.output([0-9]+)\\.Rdata$", "\\1", basename(path)))
  }

  infer_year <- function(sda_outputs) {
    restart <- sda_outputs$restart.list
    if (!is.null(restart) && length(restart) > 0) {
      for (el in restart) {
        if (!is.list(el)) next
        cand <- el$stop.time %||% el$stop_time
        d <- tryCatch(as.Date(cand), error = function(e) as.Date(NA))
        if (!is.na(d)) return(as.integer(format(d, "%Y")))
      }
    }
    as.integer(NA)
  }

  parse_param_indices <- function(par_name) {
    inside <- sub("^[^\\[]*\\[", "", par_name)
    inside <- sub("\\]$", "", inside)
    nums <- regmatches(inside, gregexpr("[0-9]+", inside))[[1]]
    as.integer(nums)
  }

  sda_files <- list.files(outdir, pattern = "^sda\\.output[0-9]+\\.Rdata$", full.names = TRUE)
  if (length(sda_files) == 0) {
    stop("No files matched '^sda.output[0-9]+.Rdata$' in outdir.")
  }
  sda_files <- sda_files[order(parse_timestep(sda_files))]

  rows <- list()
  idx <- 1L
  for (f in sda_files) {
    t_idx <- parse_timestep(f)
    env <- new.env(parent = emptyenv())
    load(f, envir = env)
    if (!exists("sda.outputs", envir = env, inherits = FALSE)) next

    sda_outputs <- get("sda.outputs", envir = env)
    enkf <- sda_outputs$enkf.params
    block_all <- enkf$block.list.all
    if (is.null(block_all) || length(block_all) < t_idx || is.null(block_all[[t_idx]])) next

    year_value <- infer_year(sda_outputs)
    block_list <- block_all[[t_idx]]

    for (b in seq_along(block_list)) {
      block <- block_list[[b]]
      by_param <- block$diag$by_param
      if (is.null(by_param) || nrow(by_param) == 0) next
      by_param <- as.data.frame(by_param, stringsAsFactors = FALSE)

      if (!("par_name" %in% names(by_param))) next
      if (!("rhat" %in% names(by_param))) by_param$rhat <- NA_real_
      if (!("ess" %in% names(by_param))) by_param$ess <- NA_real_
      if (!("param_type" %in% names(by_param))) by_param$param_type <- sub("\\[.*$", "", by_param$par_name)
      if (!("state_var" %in% names(by_param))) by_param$state_var <- NA_character_

      state_names <- names(block$data$muf %||% numeric())
      if (is.null(state_names)) state_names <- character(0)
      h_idx <- as.integer(block$constant$H %||% integer(0))
      used_lookup <- rep(FALSE, length(state_names))
      valid_h <- h_idx[is.finite(h_idx) & h_idx > 0 & h_idx <= length(state_names)]
      if (length(valid_h) > 0) used_lookup[unique(valid_h)] <- TRUE

      out <- by_param
      out$file_id <- as.integer(file_id)
      out$year <- year_value
      out$block_id <- b
      out$site_id <- as.integer(if (length(block$site.ids %||% integer(0)) > 0) block$site.ids[1] else NA_integer_)
      out$state_index <- NA_integer_
      out$used_in_H <- NA

      ix <- grepl("^X(\\.mod)?\\[[0-9]+\\]$", out$par_name)
      if (any(ix)) {
        k <- vapply(out$par_name[ix], function(nm) parse_param_indices(nm)[1], integer(1))
        out$state_index[ix] <- k

        miss_var <- is.na(out$state_var[ix]) | out$state_var[ix] == ""
        if (any(miss_var)) {
          fill <- rep(NA_character_, length(k))
          ok <- k > 0 & k <= length(state_names)
          fill[ok] <- state_names[k[ok]]
          out$state_var[ix][miss_var] <- fill[miss_var]
        }

        used <- rep(NA, length(k))
        ok_h <- k > 0 & k <= length(used_lookup)
        used[ok_h] <- used_lookup[k[ok_h]]
        out$used_in_H[ix] <- used
      }

      iq <- grepl("^q\\[", out$par_name)
      if (any(iq)) {
        out$state_index[iq] <- NA_integer_
        out$state_var[iq] <- NA_character_
        out$used_in_H[iq] <- TRUE
      }

      out$used_in_H <- as.logical(out$used_in_H)

      out <- out[, c(
        "par_name", "rhat", "ess", "param_type", "file_id", "year",
        "block_id", "site_id", "state_index", "state_var", "used_in_H"
      ), drop = FALSE]

      rows[[idx]] <- out
      idx <- idx + 1L
    }
  }

  result <- if (length(rows) == 0) {
    data.frame(
      par_name = character(0),
      rhat = numeric(0),
      ess = numeric(0),
      param_type = character(0),
      file_id = integer(0),
      year = integer(0),
      block_id = integer(0),
      site_id = integer(0),
      state_index = integer(0),
      state_var = character(0),
      used_in_H = logical(0),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }

  if (isTRUE(write_csv)) {
    if (is.null(output_csv)) output_csv <- file.path(outdir, "mcmc_param_diag_long.csv")
    utils::write.csv(result, output_csv, row.names = FALSE)
    if (verbose) message("Saved: ", normalizePath(output_csv, mustWork = FALSE))
  }
  if (verbose) message("Rows: ", nrow(result))
  result
}


# ------------------------------ Usage ---------------------------------
# source("inspect_q_r_xmod.R")
# out <- extract_q_r_xmod(
#   outdir = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q",
#   output_dir = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q/q_r_xmod_tables"
# )
#
# View per-timestep/site values:
# head(out$q_diag)
# head(out$r_diag)
# head(out$xmod)
#
# View relationship with X.mod (pointwise and summary):
# head(out$q_xmod)
# head(out$r_xmod)
# out$q_xmod_summary
# out$r_xmod_summary
# out$plot_files
#
# Build a long parameter-diagnostics table (shape shown in the question):
# diag_long <- extract_block_mcmc_param_diag(
#   outdir = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q",
#   file_id = 1,
#   write_csv = TRUE
# )
# head(diag_long)
