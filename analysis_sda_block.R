##' @title analysis_sda_block
##' @name  analysis_sda_block
##' @author Dongchen Zhang
##' 
##' @param settings  pecan standard multi-site settings list.  
##' @param block.list.all Lists of forecast and analysis outputs for each time point of each block. If t=1, we initialize those outputs of each block with NULL from the `sda.enkf.multisite` function.
##' @param X A matrix contains ensemble forecasts with the dimensions of `[ensemble number, site number * number of state variables]`. The columns are matched with the site.ids and state variable names of the inside the `FORECAST` object in the `sda.enkf.multisite` script. 
##' @param obs.mean Lists of date times named by time points, which contains lists of sites named by site ids, which contains observation means for each state variables of each site for each time point. 
##' @param obs.cov   Lists of date times named by time points, which contains lists of sites named by site ids, which contains observation covariances for all state variables of each site for each time point. 
##' @param t time point in format of YYYY-MM-DD.
##' @param nt total length of time steps, corresponding to the `nt` variable in the `sda.enkf.multisite` function.
##' @param MCMC.args arguments for the MCMC sampling, details can be found in the roxygen strucutre for control list in the `sda.enkf.multisite` function.
##' @param block.list.all.pre pre-existed block.list.all object for passing the aqq and bqq to the current SDA run, the default is NULL. Details can be found in the roxygen structure for `pre_enkf_params` of the `sda.enkf.multisite` function.
##' @details This function will add data and constants into each block that are needed for the MCMC sampling.
##'  
##' @description This function provides the block-based MCMC sampling approach.
##' 
##' @return It returns the `build.block.xy` object and the analysis results.
##' @importFrom dplyr %>%
analysis_sda_block <- function (settings, block.list.all, X, obs.mean, obs.cov, t, nt, MCMC.args, block.list.all.pre = NULL) {
  # grab cores from settings.
  cores <- as.numeric(settings$state.data.assimilation$batch.settings$general.job$cores)
  # if we didn't assign number of CPUs in the settings.
  if (length(cores) == 0 | is.null(cores)) {
    cores <- parallel::detectCores()
  }
  cores <- cores - 1
  # if we only have one CPU.
  if (cores < 1) cores <- 1
  #convert from vector values to block lists.
  if ("try-error" %in% class(try(block.results <- PEcAnAssimSequential:::build.block.xy(settings = settings, 
                                                                 block.list.all = block.list.all, 
                                                                 X = X, 
                                                                 obs.mean = obs.mean, 
                                                                 obs.cov = obs.cov, 
                                                                 t = t)))) {
    PEcAn.logger::logger.severe("Something wrong within the build.block.xy function.")
    return(0)
  }
  #grab block.list and H from the results.
  block.list.all <- block.results[[1]]
  H <- block.results[[2]]
  Y <- block.results[[3]]
  R <- block.results[[4]]
  
  #update q.
  if ("try-error" %in% class(try(block.list.all <- PEcAnAssimSequential:::update_q(block.list.all, t, nt, aqq.Init = as.numeric(settings$state.data.assimilation$aqq.Init),
                                                            bqq.Init = as.numeric(settings$state.data.assimilation$bqq.Init),
                                                            MCMC_dat = NULL,
                                                            block.list.all.pre)))) {
    PEcAn.logger::logger.severe("Something wrong within the update_q function.")
    return(0)
  }
  
  #add initial conditions for the MCMC sampling.
  if ("try-error" %in% class(try(block.list.all[[t]] <- PEcAnAssimSequential:::MCMC_Init(block.list.all[[t]], X)))) {
    PEcAn.logger::logger.severe("Something wrong within the MCMC_Init function.")
    return(0)
  }
  
  #update MCMC args.
  block.list.all[[t]] <- block.list.all[[t]] %>% 
    purrr::map(function(l){
      l$MCMC <- MCMC.args
      l
    })
  
  #parallel for loop over each block.
  PEcAn.logger::logger.info(paste0("Running MCMC ", "for ", length(block.list.all[[t]]), " blocks"))
  cl <- parallel::makeCluster(as.numeric(cores))
  doSNOW::registerDoSNOW(cl)
  l <- NULL
  if ("try-error" %in% class(try(block.list.all[[t]] <- foreach::foreach(l = block.list.all[[t]], 
                                                                         .packages = c("Kendall", 
                                                                                       "purrr", 
                                                                                       "nimble", 
                                                                                       "PEcAnAssimSequential"),
                                                                         .export = c("MCMC_block_function")) %dopar% {MCMC_block_function(l)}))) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    PEcAn.logger::logger.severe("Something wrong within the MCMC_block_function function.")
    return(0)
  }
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  PEcAn.logger::logger.info("Completed!")
  
  #convert from block lists to vector values.
  if ("try-error" %in% class(try(V <- PEcAnAssimSequential:::block.2.vector(block.list.all[[t]], X, H, settings$state.data.assimilation$adjustment)))) {
    PEcAn.logger::logger.severe("Something wrong within the block.2.vector function.")
    return(0)
  }
  
  #return values
  return(list(block.list.all = block.list.all,
              mu.f = V$mu.f,
              Pf = V$Pf,
              mu.a = V$mu.a,
              Pa = V$Pa,
              Y = Y,
              R = R,
              analysis = V$analysis))
}
