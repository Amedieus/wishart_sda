MCMC_block_function <- function(block) {
  ##
  print("start debug version of MCMC")
  # ---------- helpers ----------
  .cat <- function(...) cat(..., "\n")
  .safe_dim <- function(x) {
    d <- dim(x)
    if (is.null(d)) return(NA_character_)
    paste(d, collapse = " x ")
  }
  .safe_len <- function(x) {
    if (is.null(x)) return(NA_integer_)
    length(x)
  }
  .safe_names <- function(x) {
    n <- names(x)
    if (is.null(n)) return(character(0))
    n
  }
  
  # ---------- global tryCatch ----------
  tryCatch({
    .cat("\n========================================")
    .cat("ENTER MCMC_block_function_debug")
    .cat("t =", block$t)
    .cat("q.type =", block$constant$q.type)
    .cat("YN =", block$constant$YN)
    .cat("========================================")
    
    .cat("Top-level names(block):")
    print(names(block))
    
    .cat("\n--- block$data summary ---")
    print(names(block$data))
    if (!is.null(block$data$muf)) {
      .cat("length(block$data$muf) =", length(block$data$muf))
      .cat("names(block$data$muf) =")
      print(names(block$data$muf))
    }
    if (!is.null(block$data$pf)) {
      .cat("dim(block$data$pf) =", .safe_dim(block$data$pf))
    }
    if (!is.null(block$data$y.censored)) {
      .cat("length(block$data$y.censored) =", length(block$data$y.censored))
      .cat("block$data$y.censored =")
      print(block$data$y.censored)
    }
    if (!is.null(block$data$aq)) {
      .cat("length(block$data$aq) =", .safe_len(block$data$aq))
      print(block$data$aq)
    }
    if (!is.null(block$data$bq)) {
      .cat("length(block$data$bq) =", .safe_len(block$data$bq))
      print(block$data$bq)
    }
    
    .cat("\n--- block$constant summary ---")
    print(names(block$constant))
    print(block$constant)
    
    if (!is.null(block$H)) {
      .cat("\n--- block$H summary ---")
      .cat("dim(block$H) =", .safe_dim(block$H))
      print(block$H)
    }
    
    .cat("\n--- block$Inits summary ---")
    print(names(block$Inits))
    if (!is.null(block$Inits$X)) {
      .cat("length/block$Inits$X =", .safe_len(block$Inits$X))
      print(block$Inits$X)
    }
    if (!is.null(block$Inits$X.mod)) {
      .cat("length/block$Inits$X.mod =", .safe_len(block$Inits$X.mod))
      print(block$Inits$X.mod)
    }
    if (!is.null(block$Inits$q)) {
      .cat("length/block$Inits$q =", .safe_len(block$Inits$q))
      print(block$Inits$q)
    }
    
    # ---------- basic consistency checks ----------
    .cat("\n--- preflight consistency checks ---")
    
    if (is.null(block$data$muf)) {
      stop("block$data$muf is NULL", call. = FALSE)
    }
    if (is.null(block$data$pf)) {
      stop("block$data$pf is NULL", call. = FALSE)
    }
    if (is.null(block$data$y.censored)) {
      stop("block$data$y.censored is NULL", call. = FALSE)
    }
    if (is.null(block$constant$YN)) {
      stop("block$constant$YN is NULL", call. = FALSE)
    }
    if (length(block$data$y.censored) != block$constant$YN) {
      stop(
        paste0(
          "length(block$data$y.censored) != block$constant$YN : ",
          length(block$data$y.censored), " vs ", block$constant$YN
        ),
        call. = FALSE
      )
    }
    
    state_n <- length(block$data$muf)
    .cat("state_n =", state_n)
    
    pf_dim <- dim(block$data$pf)
    if (is.null(pf_dim) || length(pf_dim) != 2) {
      stop("block$data$pf is not a 2D matrix", call. = FALSE)
    }
    if (pf_dim[1] != state_n || pf_dim[2] != state_n) {
      stop(
        paste0(
          "dim(block$data$pf) does not match state_n: ",
          paste(pf_dim, collapse = " x "), " vs state_n=", state_n
        ),
        call. = FALSE
      )
    }
    
    if (!is.null(block$H)) {
      H_dim <- dim(block$H)
      if (!is.null(H_dim) && length(H_dim) == 2) {
        if (H_dim[2] != state_n) {
          stop(
            paste0(
              "ncol(block$H) != state_n : ",
              H_dim[2], " vs ", state_n
            ),
            call. = FALSE
          )
        }
      }
    }
    
    # ---------- nimble options ----------
    nimbleOptions(
      verbose = FALSE,
      MCMCprogressBar = FALSE,
      checkNimbleFunction = FALSE,
      checkDuplicateNodeDefinitions = FALSE
    )
    
    .cat("\n--- building nimble model ---")
    model_pred <- nimble::nimbleModel(
      GEF.MultiSite.Nimble,
      data = block$data,
      inits = block$Inits,
      constants = block$constant,
      name = "base"
    )
    
    .cat("nimble model built OK")
    
    # ---------- configure MCMC ----------
    .cat("\n--- configureMCMC ---")
    conf <- nimble::configureMCMC(model_pred, print = FALSE)
    conf$setMonitors(c("X", "X.mod", "q"))
    
    samplerLists <- conf$getSamplers()
    samplerNumberOffset <- length(samplerLists)
    
    .cat("length(samplerLists) =", samplerNumberOffset)
    
    sampler_targets <- tryCatch(
      purrr::map_chr(samplerLists, ~ .x$target),
      error = function(e) {
        .cat("Could not extract sampler targets with map_chr; fallback to manual extraction")
        unlist(purrr::map(samplerLists, ~ .x$target))
      }
    )
    
    .cat("sampler targets:")
    print(sampler_targets)
    
    # ---------- q.type == 4 sampler replacement ----------
    if (block$constant$q.type == 4) {
      .cat("\n--- replacing all samplers with ESS for q.type == 4 ---")
      samplerLists <- purrr::map(samplerLists, function(l) {
        l$setName("ess")
        l
      })
    }
    
    conf$setSamplers(samplerLists)
    
    # ---------- X.mod sampler replacement ----------
    .cat("\n--- checking X.mod sampler ---")
    xmod_ind <- grep("X.mod", sampler_targets)
    .cat("xmod_ind =")
    print(xmod_ind)
    
    if (length(xmod_ind) == 1) {
      .cat("Exactly one X.mod sampler found. Replacing with ESS.")
      .cat("target =", samplerLists[[xmod_ind]]$target)
      conf$removeSampler(samplerLists[[xmod_ind]]$target)
      conf$addSampler(
        target = samplerLists[[xmod_ind]]$target,
        type = "ess",
        control = list(
          propCov = block$data$pf,
          adaptScaleOnly = TRUE,
          latents = "X",
          pfOptimizeNparticles = TRUE
        )
      )
    } else {
      .cat(
        "WARNING: expected exactly one X.mod sampler, found",
        length(xmod_ind),
        "- skipping X.mod replacement."
      )
    }
    
    .cat("\n--- adding toggle samplers for y.censored ---")
    for (i in seq_len(block$constant$YN)) {
      node_name <- paste0("y.censored[", i, "]")
      .cat("adding toggle sampler for", node_name)
      conf$addSampler(node_name, "toggle", control = list(type = "RW"))
    }
    
    # ---------- compile ----------
    .cat("\n--- buildMCMC / compileNimble ---")
    Rmcmc <- nimble::buildMCMC(conf)
    Cmodel <- nimble::compileNimble(model_pred)
    Cmcmc <- nimble::compileNimble(Rmcmc, project = model_pred, showCompilerOutput = FALSE)
    .cat("compile OK")
    
    # ---------- toggle disable ----------
    .cat("\n--- checking compiled samplerFunctions ---")
    samplerFns <- Cmcmc$samplerFunctions
    nsam <- length(samplerFns)
    .cat("length(Cmcmc$samplerFunctions) =", nsam)
    .cat("samplerNumberOffset =", samplerNumberOffset)
    .cat("YN =", block$constant$YN)
    
    if (!any(is.na(block$data$y.censored))) {
      .cat("No NA in y.censored; attempt to disable toggle samplers.")
      for (i in seq_len(block$constant$YN)) {
        idx <- samplerNumberOffset + i
        .cat("toggle sampler idx =", idx)
        if (idx > nsam) {
          stop(
            paste0(
              "toggle sampler index out of bounds: idx=",
              idx,
              ", nsam=",
              nsam,
              ", samplerNumberOffset=",
              samplerNumberOffset,
              ", YN=",
              block$constant$YN
            ),
            call. = FALSE
          )
        }
        valueInCompiledNimbleFunction(samplerFns[[idx]], "toggle", 0)
      }
    } else {
      .cat("y.censored contains NA; keep toggle samplers active.")
    }
    
    # ---------- run MCMC ----------
    .cat("\n--- runMCMC ---")
    dat <- runMCMC(
      Cmcmc,
      niter = block$MCMC$niter,
      nburnin = block$MCMC$nburnin,
      thin = block$MCMC$nthin,
      nchains = block$MCMC$nchain
    )
    .cat("runMCMC finished OK")
    
    .cat("\n--- posterior output summary ---")
    .cat("class(dat) =")
    print(class(dat))
    .cat("dim(dat) =")
    print(dim(dat))
    .cat("colnames(dat) head =")
    print(utils::head(colnames(dat), 30))
    
    # ---------- q posterior update ----------
    M <- colMeans(dat)
    block$update$aq <- block$Inits$q
    
    if (block$constant$q.type == 3) {
      .cat("\n--- q.type == 3 posterior update ---")
      aq <- bq <- rep(NA, length(block$data$y.censored))
      for (i in seq_along(aq)) {
        q_name <- paste0("q[", i, "]")
        .cat("checking posterior column:", q_name)
        if (!(q_name %in% colnames(dat))) {
          stop(
            paste0("Missing posterior column for ", q_name),
            call. = FALSE
          )
        }
        q_var <- stats::var(dat[, q_name])
        if (!is.finite(q_var) || q_var == 0) {
          stop(
            paste0("Variance of posterior ", q_name, " is invalid: ", q_var),
            call. = FALSE
          )
        }
        aq[i] <- (mean(dat[, q_name]))^2 / q_var
        bq[i] <- mean(dat[, q_name]) / q_var
      }
      
      .cat("aq =")
      print(aq)
      .cat("bq =")
      print(bq)
      
      block$aqq[, block$t + 1] <- block$aqq[, block$t]
      block$aqq[block$constant$H, block$t + 1] <- aq
      block$bqq[, block$t + 1] <- block$bqq[, block$t]
      block$bqq[block$constant$H, block$t + 1] <- bq
      
    } else if (block$constant$q.type == 4) {
      .cat("\n--- q.type == 4 posterior update ---")
      mq_idx <- grep("q", colnames(dat))
      .cat("length(mq_idx) =", length(mq_idx))
      if (length(mq_idx) == 0) {
        stop("No q columns found in posterior output for q.type == 4", call. = FALSE)
      }
      mq <- dat[, mq_idx, drop = FALSE]
      
      Hlen <- length(block$constant$H)
      .cat("length(block$constant$H) =", Hlen)
      .cat("ncol(mq) =", ncol(mq))
      
      qbar_vec <- apply(mq, 2, mean)
      if (length(qbar_vec) != Hlen * Hlen) {
        stop(
          paste0(
            "Cannot reshape q posterior mean into matrix: length(qbar_vec)=",
            length(qbar_vec),
            ", expected=",
            Hlen * Hlen
          ),
          call. = FALSE
        )
      }
      
      q.bar <- matrix(qbar_vec, Hlen, Hlen)
      
      wish.df <- function(Om, X, i, j, col) {
        (Om[i, j]^2 + Om[i, i] * Om[j, j]) / stats::var(X[, col])
      }
      
      col_idx <- matrix(seq_len(Hlen^2), Hlen, Hlen)
      WV <- matrix(0, Hlen, Hlen)
      
      for (i in seq_len(Hlen)) {
        for (j in seq_len(Hlen)) {
          WV[i, j] <- wish.df(q.bar, X = mq, i = i, j = j, col = col_idx[i, j])
        }
      }
      
      bq <- mean(WV)
      if (bq < block$constant$YN) bq <- block$constant$YN
      aq <- solve(q.bar) * bq
      
      block$aqq[, , block$t + 1] <- GrabFillMatrix(block$aqq[, , block$t], block$constant$H, aq)
      block$bqq[block$t + 1] <- bq
    }
    
    # ---------- posterior indices ----------
    .cat("\n--- extracting posterior indices for X / X.mod ---")
    iX <- grep("X[", colnames(dat), fixed = TRUE)
    iX.mod <- grep("X.mod[", colnames(dat), fixed = TRUE)
    
    .cat("length(iX) =", length(iX))
    .cat("length(iX.mod) =", length(iX.mod))
    .cat("iX =")
    print(iX)
    .cat("iX.mod =")
    print(iX.mod)
    
    if (length(iX) == 0) {
      stop("No posterior columns matched X[", call. = FALSE)
    }
    if (length(iX.mod) == 0) {
      stop("No posterior columns matched X.mod[", call. = FALSE)
    }
    
    # ---------- mua / pa ----------
    .cat("\n--- computing mua / pa ---")
    if (length(iX) == 1) {
      mua <- mean(dat[, iX])
      pa <- stats::var(dat[, iX])
    } else {
      mua <- colMeans(dat[, iX, drop = FALSE])
      pa <- stats::cov(dat[, iX, drop = FALSE])
    }
    
    # ---------- mufa / pfa ----------
    .cat("\n--- computing mufa / pfa ---")
    if (!any(is.na(block$data$y.censored))) {
      if (is.null(block$H)) {
        stop("block$H is NULL while trying to build X.all.inds", call. = FALSE)
      }
      H <- colSums(block$H)
      obs.inds <- which(H == 1)
      non.obs.inds <- which(H == 0)
      
      .cat("H =")
      print(H)
      .cat("obs.inds =")
      print(obs.inds)
      .cat("non.obs.inds =")
      print(non.obs.inds)
      
      X.all.inds <- H
      if (length(obs.inds) > 0) {
        if (length(iX) < length(obs.inds)) {
          stop(
            paste0(
              "length(iX) < length(obs.inds): ",
              length(iX), " < ", length(obs.inds)
            ),
            call. = FALSE
          )
        }
        X.all.inds[obs.inds] <- iX[seq_along(obs.inds)]
      }
      if (length(non.obs.inds) > 0) {
        if (max(non.obs.inds) > length(iX.mod)) {
          stop(
            paste0(
              "non.obs.inds exceed length(iX.mod): max(non.obs.inds)=",
              max(non.obs.inds),
              ", length(iX.mod)=",
              length(iX.mod)
            ),
            call. = FALSE
          )
        }
        X.all.inds[non.obs.inds] <- iX.mod[non.obs.inds]
      }
      
      .cat("X.all.inds =")
      print(X.all.inds)
      
      mufa <- colMeans(dat[, X.all.inds, drop = FALSE])
      pfa <- stats::cov(dat[, X.all.inds, drop = FALSE])
    } else {
      mufa <- colMeans(dat[, iX.mod, drop = FALSE])
      pfa <- stats::cov(dat[, iX.mod, drop = FALSE])
    }
    
    .cat("\n--- SUCCESS: leaving MCMC_block_function_debug ---")
    
    block$update <- list(
      aq = if (exists("aq")) aq else NULL,
      bq = if (exists("bq")) bq else NULL,
      mua = mua,
      pa = pa,
      mufa = mufa,
      pfa = pfa
    )
    
    return(block)
    
  }, error = function(e) {
    .cat("\n########################################")
    .cat("MCMC_block_function_debug FAILED")
    .cat("t =", if (!is.null(block$t)) block$t else NA)
    .cat("q.type =", if (!is.null(block$constant$q.type)) block$constant$q.type else NA)
    .cat("YN =", if (!is.null(block$constant$YN)) block$constant$YN else NA)
    .cat("ERROR MESSAGE:")
    .cat(conditionMessage(e))
    .cat("########################################\n")
    stop(e)
  })
}