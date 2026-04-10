compute_r0_from_forecast <- function(outdir = NULL,
                                     forecast_rdata = NULL,
                                     variables = NULL,
                                     corr_source = c("increment", "level"),
                                     shrink_lambda = 0.5,
                                     abs_threshold = 0.1,
                                     max_abs_corr = 0.6,
                                     min_eig = 1e-6,
                                     output_csv = NULL,
                                     output_rds = NULL,
                                     output_r_script = NULL,
                                     verbose = TRUE) {
  safe_as_numeric <- function(x, default = NA_real_) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 0 || !is.finite(x[1])) return(default)
    x[1]
  }

  parse_timestep <- function(path) {
    as.integer(sub("^sda\\.output([0-9]+)\\.Rdata$", "\\1", basename(path)))
  }

  as_forecast_matrix <- function(x) {
    if (is.null(x)) return(NULL)
    x_mat <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(x_mat)) return(NULL)
    storage.mode(x_mat) <- "double"
    x_mat
  }

  extract_forecast_list <- function(outdir, forecast_rdata) {
    if (!is.null(forecast_rdata)) {
      env <- new.env(parent = emptyenv())
      load(forecast_rdata, envir = env)
      if (!exists("forecast.all", envir = env, inherits = FALSE)) {
        stop("forecast_rdata does not contain object 'forecast.all'.")
      }
      out <- get("forecast.all", envir = env)
      if (!is.list(out) || length(out) == 0) {
        stop("'forecast.all' is empty or not a list.")
      }
      return(out)
    }

    if (is.null(outdir)) {
      stop("Provide either 'outdir' (with sda.output*.Rdata) or 'forecast_rdata'.")
    }

    files <- list.files(
      outdir,
      pattern = "^sda\\.output[0-9]+\\.Rdata$",
      full.names = TRUE
    )
    if (length(files) == 0) {
      stop("No sda.output*.Rdata files found in outdir.")
    }
    files <- files[order(parse_timestep(files))]

    out <- vector("list", length(files))
    for (i in seq_along(files)) {
      env <- new.env(parent = emptyenv())
      load(files[i], envir = env)
      if (!exists("sda.outputs", envir = env, inherits = FALSE)) {
        stop("Missing 'sda.outputs' in file: ", basename(files[i]))
      }
      sda_out <- get("sda.outputs", envir = env)
      if (is.null(sda_out$forecast)) {
        stop("Missing 'sda.outputs$forecast' in file: ", basename(files[i]))
      }
      out[[i]] <- sda_out$forecast
    }
    out
  }

  make_spd_correlation <- function(R, min_eig = 1e-6) {
    R <- as.matrix(R)
    R <- (R + t(R)) / 2
    diag(R) <- 1

    eig <- eigen(R, symmetric = TRUE)
    vals <- eig$values
    vecs <- eig$vectors
    vals[!is.finite(vals)] <- min_eig
    vals[vals < min_eig] <- min_eig
    R_spd <- vecs %*% diag(vals, nrow = length(vals)) %*% t(vecs)
    R_spd <- (R_spd + t(R_spd)) / 2

    # Convert to correlation form
    d <- sqrt(diag(R_spd))
    d[!is.finite(d) | d <= 0] <- 1
    D_inv <- diag(1 / d, nrow = length(d))
    R_spd <- D_inv %*% R_spd %*% D_inv
    R_spd <- (R_spd + t(R_spd)) / 2
    diag(R_spd) <- 1
    R_spd
  }

  matrix_to_r0_script <- function(R0, var_order) {
    p <- length(var_order)
    lines <- c(
      "## Auto-generated R0 (forecast-based correlation prior)",
      sprintf("R0 <- matrix(0, nrow = %d, ncol = %d)", p, p),
      sprintf("rownames(R0) <- c(%s)", paste(sprintf('"%s"', var_order), collapse = ", ")),
      sprintf("colnames(R0) <- c(%s)", paste(sprintf('"%s"', var_order), collapse = ", ")),
      "diag(R0) <- 1"
    )

    for (i in seq_len(p)) {
      if (i == p) break
      for (j in (i + 1):p) {
        val <- R0[i, j]
        if (is.finite(val) && abs(val) > 0) {
          lines <- c(
            lines,
            sprintf("R0[%d, %d] <- R0[%d, %d] <- %.6f", i, j, j, i, val)
          )
        }
      }
    }
    lines
  }

  corr_source <- match.arg(corr_source)
  shrink_lambda <- safe_as_numeric(shrink_lambda, 0.5)
  abs_threshold <- safe_as_numeric(abs_threshold, 0.1)
  max_abs_corr <- safe_as_numeric(max_abs_corr, 0.6)
  min_eig <- safe_as_numeric(min_eig, 1e-6)

  if (!is.finite(shrink_lambda)) shrink_lambda <- 0.5
  shrink_lambda <- max(0, min(1, shrink_lambda))
  if (!is.finite(abs_threshold) || abs_threshold < 0) abs_threshold <- 0
  if (!is.finite(max_abs_corr) || max_abs_corr <= 0) max_abs_corr <- 1
  max_abs_corr <- min(1, max_abs_corr)
  if (!is.finite(min_eig) || min_eig <= 0) min_eig <- 1e-6

  forecast_list <- extract_forecast_list(outdir = outdir, forecast_rdata = forecast_rdata)
  nt <- length(forecast_list)
  if (nt < 2 && identical(corr_source, "increment")) {
    stop("Need at least 2 timesteps for corr_source='increment'.")
  }

  key_to_var <- character(0)
  key_to_site <- character(0)
  state_order <- character(0)
  mean_vectors <- vector("list", nt)

  for (t in seq_len(nt)) {
    f_mat <- as_forecast_matrix(forecast_list[[t]])
    if (is.null(f_mat) || ncol(f_mat) == 0) {
      stop("Forecast at timestep ", t, " cannot be converted to a non-empty matrix.")
    }

    if (is.null(colnames(f_mat))) {
      colnames(f_mat) <- paste0("state_", seq_len(ncol(f_mat)))
    }
    var_names <- colnames(f_mat)
    vars_in_t <- unique(var_names)
    state_order <- c(state_order, vars_in_t[!vars_in_t %in% state_order])

    site_attr <- attr(forecast_list[[t]], "Site")
    if (is.null(site_attr) || length(site_attr) != ncol(f_mat)) {
      site_ids <- paste0("site_", seq_len(ncol(f_mat)))
    } else {
      site_ids <- as.character(site_attr)
    }

    col_keys <- paste(site_ids, var_names, sep = "||")
    means <- colMeans(f_mat, na.rm = TRUE)
    names(means) <- col_keys
    mean_vectors[[t]] <- means

    new_keys <- setdiff(col_keys, names(key_to_var))
    if (length(new_keys) > 0) {
      key_to_var[new_keys] <- var_names[match(new_keys, col_keys)]
      key_to_site[new_keys] <- site_ids[match(new_keys, col_keys)]
    }
  }

  var_order <- if (is.null(variables)) {
    unique(state_order)
  } else {
    vars <- as.character(variables)
    vars[vars %in% state_order]
  }
  if (length(var_order) < 2) {
    stop("Need at least 2 variables after filtering.")
  }

  all_keys <- sort(unique(unlist(lapply(mean_vectors, names), use.names = FALSE)))
  x_t <- matrix(NA_real_, nrow = nt, ncol = length(all_keys))
  colnames(x_t) <- all_keys
  rownames(x_t) <- paste0("t", seq_len(nt))

  for (t in seq_len(nt)) {
    v <- mean_vectors[[t]]
    x_t[t, names(v)] <- as.numeric(v)
  }

  keys_df <- data.frame(
    key = all_keys,
    site_id = unname(key_to_site[all_keys]),
    state_var = unname(key_to_var[all_keys]),
    stringsAsFactors = FALSE
  )

  sites <- unique(keys_df$site_id)
  sample_list <- vector("list", length(sites) * if (corr_source == "increment") max(0, nt - 1) else nt)
  row_names <- character(0)
  k <- 1L

  for (sid in sites) {
    row_idx <- which(keys_df$site_id == sid & keys_df$state_var %in% var_order)
    if (length(row_idx) == 0) next
    site_key_df <- keys_df[row_idx, , drop = FALSE]

    # one key per variable for this site
    key_by_var <- setNames(rep(NA_character_, length(var_order)), var_order)
    for (vv in var_order) {
      key_hits <- site_key_df$key[site_key_df$state_var == vv]
      if (length(key_hits) > 0) {
        key_by_var[[vv]] <- key_hits[1]
      }
    }

    if (identical(corr_source, "increment")) {
      for (tt in 2:nt) {
        row_val <- rep(NA_real_, length(var_order))
        names(row_val) <- var_order
        for (vv in var_order) {
          key <- key_by_var[[vv]]
          if (!is.na(key)) {
            x_now <- x_t[tt, key]
            x_prev <- x_t[tt - 1, key]
            if (is.finite(x_now) && is.finite(x_prev)) {
              row_val[vv] <- x_now - x_prev
            }
          }
        }
        if (sum(is.finite(row_val)) >= 2) {
          sample_list[[k]] <- row_val
          row_names <- c(row_names, paste0("site_", sid, "_t", tt))
          k <- k + 1L
        }
      }
    } else {
      for (tt in seq_len(nt)) {
        row_val <- rep(NA_real_, length(var_order))
        names(row_val) <- var_order
        for (vv in var_order) {
          key <- key_by_var[[vv]]
          if (!is.na(key)) {
            x_now <- x_t[tt, key]
            if (is.finite(x_now)) {
              row_val[vv] <- x_now
            }
          }
        }
        if (sum(is.finite(row_val)) >= 2) {
          sample_list[[k]] <- row_val
          row_names <- c(row_names, paste0("site_", sid, "_t", tt))
          k <- k + 1L
        }
      }
    }
  }

  sample_list <- Filter(Negate(is.null), sample_list)
  if (length(sample_list) < 3) {
    stop("Not enough valid rows to estimate variable correlation.")
  }

  sample_mat <- do.call(rbind, sample_list)
  sample_mat <- as.matrix(sample_mat)
  colnames(sample_mat) <- var_order
  rownames(sample_mat) <- row_names

  p <- ncol(sample_mat)
  if (p < 2) {
    stop("Need at least 2 variables for correlation.")
  }

  n_pair <- matrix(0L, nrow = p, ncol = p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      n_pair[i, j] <- sum(is.finite(sample_mat[, i]) & is.finite(sample_mat[, j]))
    }
  }
  colnames(n_pair) <- rownames(n_pair) <- var_order

  R_emp <- suppressWarnings(stats::cor(sample_mat, use = "pairwise.complete.obs"))
  if (is.null(dim(R_emp))) {
    stop("Correlation computation failed.")
  }
  R_emp <- as.matrix(R_emp)
  R_emp <- (R_emp + t(R_emp)) / 2
  diag(R_emp) <- 1
  R_emp[!is.finite(R_emp)] <- 0
  diag(R_emp) <- 1

  R_shrunk <- shrink_lambda * diag(p) + (1 - shrink_lambda) * R_emp
  R_shrunk <- (R_shrunk + t(R_shrunk)) / 2
  diag(R_shrunk) <- 1

  R_thresh <- R_shrunk
  offdiag <- row(R_thresh) != col(R_thresh)
  R_thresh[offdiag & abs(R_thresh) < abs_threshold] <- 0
  R_thresh[offdiag] <- pmax(pmin(R_thresh[offdiag], max_abs_corr), -max_abs_corr)
  R_thresh <- (R_thresh + t(R_thresh)) / 2
  diag(R_thresh) <- 1

  R0 <- make_spd_correlation(R_thresh, min_eig = min_eig)
  offdiag0 <- row(R0) != col(R0)
  R0[offdiag0] <- pmax(pmin(R0[offdiag0], max_abs_corr), -max_abs_corr)
  R0 <- make_spd_correlation(R0, min_eig = min_eig)

  colnames(R_emp) <- rownames(R_emp) <- var_order
  colnames(R_shrunk) <- rownames(R_shrunk) <- var_order
  colnames(R_thresh) <- rownames(R_thresh) <- var_order
  colnames(R0) <- rownames(R0) <- var_order

  if (!is.null(output_csv)) {
    out_long <- data.frame(
      var_i = rep(var_order, each = p),
      var_j = rep(var_order, times = p),
      n_pair = as.vector(n_pair),
      corr_emp = as.vector(R_emp),
      corr_shrunk = as.vector(R_shrunk),
      corr_thresholded = as.vector(R_thresh),
      corr_final = as.vector(R0),
      stringsAsFactors = FALSE
    )
    utils::write.csv(out_long, output_csv, row.names = FALSE)
  }

  if (!is.null(output_r_script)) {
    writeLines(matrix_to_r0_script(R0, var_order), con = output_r_script)
  }

  if (!is.null(output_rds)) {
    saveRDS(
      list(
        corr_source = corr_source,
        settings = list(
          shrink_lambda = shrink_lambda,
          abs_threshold = abs_threshold,
          max_abs_corr = max_abs_corr,
          min_eig = min_eig
        ),
        var_order = var_order,
        n_sample_rows = nrow(sample_mat),
        n_pair = n_pair,
        corr_emp = R_emp,
        corr_shrunk = R_shrunk,
        corr_thresholded = R_thresh,
        R0 = R0
      ),
      file = output_rds
    )
  }

  if (isTRUE(verbose)) {
    message("Timesteps: ", nt)
    message("corr_source: ", corr_source)
    message("Variables (order): ", paste(var_order, collapse = ", "))
    message("Sample rows used: ", nrow(sample_mat))
    message("R0 off-diagonal entries:")
    for (i in seq_len(p)) {
      if (i == p) break
      for (j in (i + 1):p) {
        if (abs(R0[i, j]) > 0) {
          message("  ", var_order[i], " ~ ", var_order[j], " : ", signif(R0[i, j], 6))
        }
      }
    }
  }

  invisible(list(
    corr_source = corr_source,
    var_order = var_order,
    sample_matrix = sample_mat,
    n_pair = n_pair,
    corr_emp = R_emp,
    corr_shrunk = R_shrunk,
    corr_thresholded = R_thresh,
    R0 = R0
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

  vars_raw <- get_arg("variables", NULL)
  vars <- if (is.null(vars_raw)) NULL else strsplit(vars_raw, ",", fixed = TRUE)[[1]]

  compute_r0_from_forecast(
    outdir = get_arg("outdir", NULL),
    forecast_rdata = get_arg("forecast_rdata", NULL),
    variables = vars,
    corr_source = get_arg("corr_source", "increment"),
    shrink_lambda = get_arg("shrink_lambda", "0.5"),
    abs_threshold = get_arg("abs_threshold", "0.1"),
    max_abs_corr = get_arg("max_abs_corr", "0.6"),
    min_eig = get_arg("min_eig", "1e-6"),
    output_csv = get_arg("output_csv", NULL),
    output_rds = get_arg("output_rds", NULL),
    output_r_script = get_arg("output_r_script", NULL),
    verbose = TRUE
  )
}

res <- compute_r0_from_forecast(
  outdir = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_forecast"
)