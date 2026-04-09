extract_q_r_xmod <- function(outdir,
                             output_dir = file.path(outdir, "q_r_xmod_tables"),
                             write_csv = TRUE,
                             verbose = TRUE) {
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

    date_label <- paste0("timestep_", t_idx)

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
          r_mat_long$block_id <- b
          r_matrix_all[[length(r_matrix_all) + 1L]] <- r_mat_long
        }
      }
    }
  }

  q_diag <- safe_bind(q_diag_all)
  q_matrix <- safe_bind(q_matrix_all)
  r_diag <- safe_bind(r_diag_all)
  r_matrix <- safe_bind(r_matrix_all)
  xmod <- safe_bind(xmod_all)

  q_xmod <- data.frame()
  if (all(c("timestep", "block_id", "site_id", "state_var", "q_diag") %in% names(q_diag)) &&
      all(c("timestep", "block_id", "site_id", "state_var", "xmod_mean") %in% names(xmod))) {
    q_xmod <- merge(
      q_diag[, c("timestep", "date", "block_id", "site_id", "state_var", "q_diag"), drop = FALSE],
      xmod[, c("timestep", "date", "block_id", "site_id", "state_var", "xmod_mean"), drop = FALSE],
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
  }

  r_xmod <- data.frame()
  if (all(c("timestep", "block_id", "site_id", "state_var", "r_diag") %in% names(r_diag)) &&
      all(c("timestep", "block_id", "site_id", "state_var", "xmod_mean") %in% names(xmod))) {
    r_xmod <- merge(
      r_diag[, c("timestep", "date", "block_id", "site_id", "state_var", "r_diag"), drop = FALSE],
      xmod[, c("timestep", "date", "block_id", "site_id", "state_var", "xmod_mean"), drop = FALSE],
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

  results <- list(
    q_diag = safe_order(q_diag, c("timestep", "site_id", "state_var")),
    q_matrix = safe_order(q_matrix, c("timestep", "block_id", "row_obs_idx", "col_obs_idx")),
    r_diag = safe_order(r_diag, c("timestep", "site_id", "state_var")),
    r_matrix = safe_order(r_matrix, c("timestep", "block_id", "row_obs_idx", "col_obs_idx")),
    xmod = safe_order(xmod, c("timestep", "site_id", "state_var")),
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
    utils::write.csv(results$q_xmod, file.path(output_dir, "q_vs_xmod_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$r_xmod, file.path(output_dir, "r_vs_xmod_by_timestep_site.csv"), row.names = FALSE)
    utils::write.csv(results$q_xmod_summary, file.path(output_dir, "q_vs_xmod_summary.csv"), row.names = FALSE)
    utils::write.csv(results$r_xmod_summary, file.path(output_dir, "r_vs_xmod_summary.csv"), row.names = FALSE)

    if (verbose) {
      message("Saved tables to: ", normalizePath(output_dir, mustWork = FALSE))
    }
  }

  if (verbose) {
    timesteps <- if ("timestep" %in% names(results$xmod)) {
      length(unique(results$xmod$timestep))
    } else {
      0
    }
    message("Timesteps processed: ", timesteps)
    message("Rows - Q diag: ", nrow(results$q_diag),
            ", R diag: ", nrow(results$r_diag),
            ", X.mod: ", nrow(results$xmod))
  }

  results
}


# ------------------------------ Usage ---------------------------------
# source("inspect_q_r_xmod.R")
# out <- extract_q_r_xmod(
#   outdir = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output",
#   output_dir = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output/q_r_xmod_tables"
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
