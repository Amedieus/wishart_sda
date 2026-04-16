MCMC_block_function <- function(block) {
  # disable printing out messages.
  nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE, checkNimbleFunction = FALSE, checkDuplicateNodeDefinitions = FALSE)
  #build nimble model
  #TODO: harmonize the MCMC code between block-based and general analysis functions to reduce the complexity of code.
  model_pred <- nimble::nimbleModel(GEF.MultiSite.Nimble,
                                    data = block$data,
                                    inits = block$Inits,
                                    constants = block$constant,
                                    name = 'base')
  #configure MCMC
  conf <- nimble::configureMCMC(model_pred, print=FALSE)
  conf$setMonitors(c("X", "X.mod", "q"))
  
  #Handle samplers
  #hear we change the RW_block sampler to the ess sampler 
  #because it has a better performance of MVN sampling
  samplerLists <- conf$getSamplers()
  samplerNumberOffset <- length(samplerLists)
  if (block$constant$q.type == 4) {
    #if we have wishart q
    #everything should be sampled with ess sampler.
    samplerLists %>% purrr::map(function(l){l$setName("ess")})
  }
  conf$setSamplers(samplerLists)
  
  #add Pf as propCov in the control list of the X.mod nodes.
  X.mod.ind <- which(grepl("X.mod", samplerLists %>% purrr::map(~ .x$target) %>% unlist()))
  conf$removeSampler(samplerLists[[X.mod.ind]]$target)
  conf$addSampler(target = samplerLists[[X.mod.ind]]$target, type = "ess",
                  control = list(propCov= block$data$pf, adaptScaleOnly = TRUE,
                                 latents = "X", pfOptimizeNparticles = TRUE))
  #add toggle Y sampler.
  # Revise 1: only those who needed is y.censored
  na_idx <- which(is.na(block$data$y.censored))
  
  if (length(na_idx) > 0) {
    for (i in na_idx) {
      conf$addSampler(
        target = paste0("y.censored[", i, "]"),
        type = "toggle",
        control = list(type = "RW")
      )
    }
  }
  
    # conf$printSamplers()
  #compile MCMC
  Rmcmc <- nimble::buildMCMC(conf)
  Cmodel <- nimble::compileNimble(model_pred)
  Cmcmc <- nimble::compileNimble(Rmcmc, project = model_pred, showCompilerOutput = FALSE)
  
  # Revise 2
  #if we don't have any NA in the Y.
  # if (!any(is.na(block$data$y.censored))) {
  #   #add toggle Y sampler.
  #   for(i in 1:block$constant$YN) {
  #     valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[samplerNumberOffset+i]], 'toggle', 0)
  #   }
  # }
  
  #run MCMC
  dat <- runMCMC(Cmcmc, niter = block$MCMC$niter, nburnin = block$MCMC$nburnin, thin = block$MCMC$nthin, nchains = block$MCMC$nchain)
  #update aq, bq, mua, and pa
  M <- colMeans(dat)
  block$update$aq <- block$Inits$q
  if (block$constant$q.type == 3) {
    #if it's a vector q case
    aq <- bq <- rep(NA, length(block$data$y.censored))
    for (i in seq_along(aq)) {
      CHAR <- paste0("[", i, "]")
      aq[i] <- (mean(dat[, paste0("q", CHAR)]))^2/stats::var(dat[, paste0("q", CHAR)])
      bq[i] <- mean(dat[, paste0("q", CHAR)])/stats::var(dat[, paste0("q", CHAR)])
    }
    #update aqq and bqq
    block$aqq[,block$t+1] <- block$aqq[, block$t]
    block$aqq[block$constant$H, block$t+1] <- aq
    block$bqq[,block$t+1] <- block$bqq[, block$t]
    block$bqq[block$constant$H, block$t+1] <- bq
  } else if (block$constant$q.type == 4) {
    #previous updates
    mq <- dat[,  grep("q", colnames(dat))]  # Omega, Precision
    q.bar <- matrix(apply(mq, 2, mean),
                    length(block$constant$H),
                    length(block$constant$H)
    )
    wish.df <- function(Om, X, i, j, col) {
      (Om[i, j]^2 + Om[i, i] * Om[j, j]) / stats::var(X[, col])
    }
    col <- matrix(1:length(block$constant$H) ^ 2,
                  length(block$constant$H),
                  length(block$constant$H))
    WV  <- matrix(0, length(block$constant$H), length(block$constant$H))
    for (i in seq_along(block$constant$H)) {
      for (j in seq_along(block$constant$H)) {
        WV[i, j] <- wish.df(q.bar, X = mq, i = i, j = j, col = col[i, j])
      }
    }
    bq <- mean(WV)
    if (bq < block$constant$YN) {
      bq <- block$constant$YN
    }
    aq <- solve(q.bar) * bq
    block$aqq[,,block$t+1] <- GrabFillMatrix(block$aqq[,,block$t], block$constant$H, aq)
    block$bqq[block$t+1] <- bq
  }
  #update mua and pa; mufa, and pfa
  iX <- grep("X[", colnames(dat), fixed = TRUE)
  iX.mod <- grep("X.mod[", colnames(dat), fixed = TRUE)
  if (length(iX) == 1) {
    mua <- mean(dat[, iX])
    pa <- stats::var(dat[, iX])
  } else {
    mua <- colMeans(dat[, iX])
    pa <- stats::cov(dat[, iX])
  }
  # construct X.all object.
  # NA only occurs when there is zero observation.
  if (!any(is.na(block$data$y.censored))) {
    H <- colSums(block$H)
    obs.inds <- which(H == 1)
    non.obs.inds <- which(H == 0)
    X.all.inds <- H
    X.all.inds[obs.inds] <- iX
    X.all.inds[non.obs.inds] <- iX.mod[non.obs.inds]
    mufa <- colMeans(dat[, X.all.inds])
    pfa <- stats::cov(dat[, X.all.inds])
  } else {
    mufa <- colMeans(dat[, iX.mod])
    pfa <- stats::cov(dat[, iX.mod])
  }
  #return values.
  block$update <- list(aq = aq, bq = bq, mua = mua, pa = pa, mufa = mufa, pfa = pfa)
  return(block)
}

