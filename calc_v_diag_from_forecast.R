compute_v_diag_from_forecast <- function(outdir = NULL,
                                         forecast_rdata = NULL,
                                         variables = NULL,
                                         normalize = TRUE,
                                         normalize_mode = c("variable", "site_variable"),
                                         output_csv = NULL,
                                         output_rds = NULL,
                                         verbose = TRUE) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  safe_var <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) {
      return(NA_real_)
    }
    stats::var(x)
  }

  parse_timestep <- function(path) {
    as.integer(sub("^sda\\.output([0-9]+)\\.Rdata$", "\\1", basename(path)))
  }

  as_forecast_matrix <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    x_mat <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(x_mat)) {
      return(NULL)
    }
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

  forecast_list <- extract_forecast_list(outdir = outdir, forecast_rdata = forecast_rdata)
  nt <- length(forecast_list)
  if (nt < 2) {
    stop("Need at least 2 timesteps to compute (x_t - x_{t-1}).")
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

  all_keys <- sort(unique(unlist(lapply(mean_vectors, names), use.names = FALSE)))
  x_t <- matrix(NA_real_, nrow = nt, ncol = length(all_keys))
  colnames(x_t) <- all_keys
  rownames(x_t) <- paste0("t", seq_len(nt))

  for (t in seq_len(nt)) {
    v <- mean_vectors[[t]]
    x_t[t, names(v)] <- as.numeric(v)
  }

  # Keep raw-difference variances for diagnostics/reference.
  dx_raw <- x_t[2:nt, , drop = FALSE] - x_t[1:(nt - 1), , drop = FALSE]
  var_dx_raw_by_key <- apply(dx_raw, 2, safe_var)

  # Unit harmonization helpers.
  zscore_col <- function(x) {
    x_ok <- x[is.finite(x)]
    mu <- if (length(x_ok) > 0) mean(x_ok) else 0
    sdv <- if (length(x_ok) > 1) stats::sd(x_ok) else 1
    if (!is.finite(mu)) mu <- 0
    if (!is.finite(sdv) || sdv <= 0) sdv <- 1
    (x - mu) / sdv
  }

  normalize_mode <- match.arg(normalize_mode)
  normalize <- isTRUE(normalize)
  x_used <- x_t
  if (normalize) {
    if (normalize_mode == "site_variable") {
      # z-score each site-variable trajectory over time.
      x_used <- apply(x_t, 2, zscore_col)
    } else if (normalize_mode == "variable") {
      # z-score by state variable across all sites and timesteps.
      x_used <- x_t
      state_by_col <- unname(key_to_var[colnames(x_t)])
      for (sv in unique(state_by_col)) {
        idx <- which(state_by_col == sv)
        vals <- as.vector(x_t[, idx, drop = FALSE])
        vals <- vals[is.finite(vals)]
        mu <- if (length(vals) > 0) mean(vals) else 0
        sdv <- if (length(vals) > 1) stats::sd(vals) else 1
        if (!is.finite(mu)) mu <- 0
        if (!is.finite(sdv) || sdv <= 0) sdv <- 1
        x_used[, idx] <- (x_t[, idx, drop = FALSE] - mu) / sdv
      }
    }
  }
  if (is.vector(x_used)) {
    x_used <- matrix(x_used, nrow = nrow(x_t), ncol = ncol(x_t))
    colnames(x_used) <- colnames(x_t)
    rownames(x_used) <- rownames(x_t)
  }

  dx_used <- x_used[2:nt, , drop = FALSE] - x_used[1:(nt - 1), , drop = FALSE]
  var_dx_norm_by_key <- apply(dx_used, 2, safe_var)

  key_tbl <- data.frame(
    key = names(var_dx_norm_by_key),
    site_id = unname(key_to_site[names(var_dx_norm_by_key)]),
    state_var = unname(key_to_var[names(var_dx_norm_by_key)]),
    var_dx_raw = as.numeric(var_dx_raw_by_key[names(var_dx_norm_by_key)]),
    var_dx_norm = as.numeric(var_dx_norm_by_key),
    stringsAsFactors = FALSE
  )

  key_tbl <- key_tbl[is.finite(key_tbl$var_dx_norm) & is.finite(key_tbl$var_dx_raw), , drop = FALSE]
  if (nrow(key_tbl) == 0) {
    stop("No finite variance estimates from (x_t - x_{t-1}).")
  }

  var_by_variable_norm <- stats::aggregate(
    var_dx_norm ~ state_var,
    data = key_tbl,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(var_by_variable_norm)[names(var_by_variable_norm) == "var_dx_norm"] <- "mean_var_dx_norm"

  var_by_variable_raw <- stats::aggregate(
    var_dx_raw ~ state_var,
    data = key_tbl,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(var_by_variable_raw)[names(var_by_variable_raw) == "var_dx_raw"] <- "mean_var_dx_raw"
  var_by_variable <- merge(var_by_variable_norm, var_by_variable_raw, by = "state_var", all = TRUE)

  if (!is.null(variables)) {
    variables <- as.character(variables)
    keep <- variables[variables %in% var_by_variable$state_var]
    if (length(keep) == 0) {
      stop("None of the provided 'variables' were found in forecast columns.")
    }
    var_by_variable <- var_by_variable[match(keep, var_by_variable$state_var), , drop = FALSE]
  } else {
    # Preserve state order from forecast columns (not alphabetical merge order).
    keep_order <- state_order[state_order %in% var_by_variable$state_var]
    var_by_variable <- var_by_variable[match(keep_order, var_by_variable$state_var), , drop = FALSE]
  }

  tau_raw <- mean(var_by_variable$mean_var_dx_raw, na.rm = TRUE)
  tau_norm <- mean(var_by_variable$mean_var_dx_norm, na.rm = TRUE)
  if (!is.finite(tau_raw) || tau_raw <= 0) {
    stop("Computed tau_raw is non-finite or <= 0; cannot normalize.")
  }
  if (!is.finite(tau_norm) || tau_norm <= 0) {
    stop("Computed tau_norm is non-finite or <= 0; cannot normalize.")
  }
  tau <- if (normalize) tau_norm else tau_raw

  var_by_variable$tau_raw <- tau_raw
  var_by_variable$tau_norm <- tau_norm
  var_by_variable$tau <- tau
  var_by_variable$v_diag_raw <- var_by_variable$mean_var_dx_raw / tau_raw
  var_by_variable$v_diag_norm <- var_by_variable$mean_var_dx_norm / tau_norm
  var_by_variable$v_diag <- if (normalize) var_by_variable$v_diag_norm else var_by_variable$v_diag_raw

  # optional outputs
  if (!is.null(output_csv)) {
    utils::write.csv(var_by_variable, output_csv, row.names = FALSE)
  }
  if (!is.null(output_rds)) {
    saveRDS(
      list(
        normalized = normalize,
        normalize_mode = normalize_mode,
        tau = tau,
        tau_raw = tau_raw,
        tau_norm = tau_norm,
        v_diag = stats::setNames(var_by_variable$v_diag, var_by_variable$state_var),
        v_diag_raw = stats::setNames(var_by_variable$v_diag_raw, var_by_variable$state_var),
        v_diag_norm = stats::setNames(var_by_variable$v_diag_norm, var_by_variable$state_var),
        by_variable = var_by_variable,
        by_site_variable = key_tbl
      ),
      file = output_rds
    )
  }

  if (isTRUE(verbose)) {
    message("Timesteps: ", nt, " (pairs: ", nt - 1, ")")
    if (normalize) {
      message("Using normalized trajectories (mode = ", normalize_mode, ").")
      message("Tau_norm (mean variance of normalized x_t - x_{t-1}): ", signif(tau_norm, 6))
      message("Tau_raw  (mean variance of raw x_t - x_{t-1}): ", signif(tau_raw, 6))
    } else {
      message("Using raw trajectories (no normalization).")
      message("Tau_raw (mean variance of raw x_t - x_{t-1}): ", signif(tau_raw, 6))
      message("Tau_norm (for reference): ", signif(tau_norm, 6))
    }
    message("v_diag:")
    for (i in seq_len(nrow(var_by_variable))) {
      message("  ", var_by_variable$state_var[i], " = ", signif(var_by_variable$v_diag[i], 6))
    }
  }

  invisible(list(
    normalized = normalize,
    normalize_mode = normalize_mode,
    tau = tau,
    tau_raw = tau_raw,
    tau_norm = tau_norm,
    v_diag = stats::setNames(var_by_variable$v_diag, var_by_variable$state_var),
    v_diag_raw = stats::setNames(var_by_variable$v_diag_raw, var_by_variable$state_var),
    v_diag_norm = stats::setNames(var_by_variable$v_diag_norm, var_by_variable$state_var),
    by_variable = var_by_variable,
    by_site_variable = key_tbl
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

  vars_raw <- get_arg("variables", NULL)
  vars <- if (is.null(vars_raw)) NULL else strsplit(vars_raw, ",", fixed = TRUE)[[1]]
  normalize_mode <- get_arg("normalize_mode", "variable")

  compute_v_diag_from_forecast(
    outdir = get_arg("outdir", NULL),
    forecast_rdata = get_arg("forecast_rdata", NULL),
    variables = vars,
    normalize = parse_bool(get_arg("normalize", "true"), TRUE),
    normalize_mode = normalize_mode,
    output_csv = get_arg("output_csv", NULL),
    output_rds = get_arg("output_rds", NULL),
    verbose = TRUE
  )
}
