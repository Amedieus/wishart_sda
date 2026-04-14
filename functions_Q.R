######### Change the 

#' sda_matchparam
#'
#' @name sda_matchparam
#' @author Alexis Helgeson
#' 
#' @param settings settings object passed from sda.enkf_MultiSite
#' @param ensemble.samples taken from sample.Rdata object
#' @param site.ids character object passed from sda.enkf_MultiSite
#' @param nens number of ensemble members in model runs, taken from restart$runids
#'
#' @return new.params object used to 
sda_matchparam <- function(settings, ensemble.samples, site.ids, nens){
  #reformatting params
  new.params <- list()
  all.pft.names <- names(ensemble.samples)
  
  #loop over each site.
  for (i in seq_along(site.ids)) {
    #match pft name
    site.pft.name <- settings[[i]]$run$site$site.pft$pft.name
    if(is.null(site.pft.name)){
      site_pft = utils::read.csv(settings[[i]]$run$inputs$pft.site$path)
      site.pft.name = site_pft$pft[site_pft$site == settings[[i]]$run$site$id]
    }
    which.pft <- which(all.pft.names==site.pft.name)
    
    site.param <- list()
    site.samples <- ensemble.samples[which.pft]
    for (j in seq_len(nens)) {
      site.param[[j]] <- lapply(site.samples, function(x, n) {
        x[j, ] }, n = j)
    } 
    new.params[[i]] <- site.param
  }
  names(new.params) <- site.ids
  
  return(new.params)
}




write.ensemble.configs <- function (defaults, ensemble.samples, settings, model, input_design, 
                                    clean = FALSE, write.to.db = TRUE, restart = NULL, samples = NULL, 
                                    rename = TRUE) 
{
  
  # print(names(settings$run$inputs))
  
  for (input_tag in names(settings$run$inputs)) {
    input <- settings$run$inputs[[input_tag]]
    input_paths <- input$path
    if (is.null(input_paths) || length(input_paths) == 0) {
      PEcAn.logger::logger.error("Input", sQuote(input_tag), 
                                 "has no paths specified")
    }
    if (length(input_paths) > 1 && !(input_tag %in% names(settings$ensemble$samplingspace))) {
      PEcAn.logger::logger.error("Input", sQuote(input_tag), 
                                 "has", length(input_paths), "paths but no sampling method.", 
                                 "Add <samplingspace> for this input in pecan.xml")
    }
  }
  con <- NULL
  ####
  write_config.SIPNET <- write.config.SIPNET
  my.write.config <- paste("write_config.", model, sep = "")
  my.write_restart <- paste0("write_restart.", model)
  if (is.null(ensemble.samples)) {
    return(list(runs = NULL, ensemble.id = NULL))
  }
  if (!is.null(settings$database$bety$write)) {
    write.to.db <- as.logical(settings$database$bety$write)
  }
  if (write.to.db) {
    con <- try(PEcAn.DB::db.open(settings$database$bety))
    on.exit(try(PEcAn.DB::db.close(con), silent = TRUE), 
            add = TRUE)
    if (inherits(con, "try-error")) {
      con <- NULL
      PEcAn.logger::logger.warn("We were not able to successfully establish a connection with BETYdb ")
    }
  }
  if (!is.null(settings$workflow$id)) {
    workflow.id <- settings$workflow$id
  }
  else {
    workflow.id <- -1
  }
  if (is.null(restart)) {
    if (!is.null(con) && write.to.db) {
      ensemble.id <- PEcAn.DB::db.query(paste0("INSERT INTO ensembles (runtype, workflow_id) ", 
                                               "VALUES ('ensemble', ", format(workflow.id, scientific = FALSE), 
                                               ")", "RETURNING id"), con = con)[["id"]]
      for (pft in defaults) {
        PEcAn.DB::db.query(paste0("INSERT INTO posteriors_ensembles (posterior_id, ensemble_id) ", 
                                  "values (", pft$posteriorid, ", ", ensemble.id, 
                                  ")"), con = con)
      }
    }
    else {
      ensemble.id <- NA
    }
    if (!is.null(con)) {
      required_tags <- dplyr::tbl(con, "models") %>% dplyr::filter(.data$id == 
                                                                     !!as.numeric(settings$model$id)) %>% dplyr::inner_join(dplyr::tbl(con, 
                                                                                                                                       "modeltypes_formats"), by = c("modeltype_id")) %>% 
        dplyr::collect() %>% dplyr::filter(.data$required == 
                                             TRUE) %>% dplyr::pull("tag")
    }
    else {
      required_tags <- c("met", "parameters")
    }
    samp <- settings$ensemble$samplingspace
    parents <- lapply(samp, "[[", "parent")
    order <- names(samp)[lapply(parents, function(tr) which(names(samp) %in% 
                                                              tr)) %>% unlist()]
    samp.ordered <- samp[c(order, names(samp)[!(names(samp) %in% 
                                                  order)])]
    if (is.null(samples)) {
      samples <- list()
      input_tags <- names(settings$run$inputs)
      for (input_tag in input_tags) {
        if (input_tag %in% colnames(input_design)) {
          input_paths <- settings$run$inputs[[input_tag]]$path
          input_indices <- input_design[[input_tag]]
          samples[[input_tag]] <- list(samples = lapply(input_indices, 
                                                        function(idx) input_paths[[idx]]))
        }
      }
    }
    required_tags %>% purrr::walk(function(r_tag) {
      if (is.null(samples[[r_tag]]) & r_tag != "parameters") 
        samples[[r_tag]]$samples <<- rep(settings$run$inputs[[tolower(r_tag)]]$path[1], 
                                         settings$ensemble$size)
    })
    site.pfts.vec <- settings$run$site$site.pft %>% unlist %>% 
      as.character
    if (!is.null(site.pfts.vec)) {
      defined.pfts <- settings$pfts %>% purrr::map("name") %>% 
        unlist %>% as.character
      if (length(which(site.pfts.vec %in% defined.pfts)) > 
          0) 
        ensemble.samples <- ensemble.samples[site.pfts.vec[which(site.pfts.vec %in% 
                                                                   defined.pfts)]]
      if (length(which(!(site.pfts.vec %in% defined.pfts))) > 
          0) 
        PEcAn.logger::logger.warn(paste0("The following pfts are specified for the siteid ", 
                                         settings$run$site$id, " but they are not defined as a pft in pecan.xml:", 
                                         site.pfts.vec[which(!(site.pfts.vec %in% defined.pfts))], 
                                         collapse = ","))
    }
    if (is.null(samp$parameters)) 
      samples$parameters$samples <- ensemble.samples %>% 
      purrr::map(~.x[rep(1, settings$ensemble$size), 
      ])
    if (is.null(samples$parameters$samples)) 
      samples$parameters$samples <- ensemble.samples
    inputs <- names(settings$run$inputs)
    inputs <- inputs[grepl(".id$", inputs)]
    runs <- data.frame()
    for (i in seq_len(settings$ensemble$size)) {
      if (!is.null(con) && write.to.db) {
        paramlist <- paste("ensemble=", i, sep = "")
        run.id <- PEcAn.DB::db.query(paste0("INSERT INTO runs (model_id, site_id, start_time, finish_time, outdir, ensemble_id, parameter_list) ", 
                                            "values ('", settings$model$id, "', '", settings$run$site$id, 
                                            "', '", settings$run$start.date, "', '", settings$run$end.date, 
                                            "', '", settings$run$outdir, "', ", ensemble.id, 
                                            ", '", paramlist, "') ", "RETURNING id"), con = con)[["id"]]
        if (!is.null(inputs)) {
          for (x in inputs) {
            PEcAn.DB::db.query(paste0("INSERT INTO inputs_runs (input_id, run_id) ", 
                                      "values (", settings$run$inputs[[x]], ", ", 
                                      run.id, ")"), con = con)
          }
        }
      }
      else {
        run.id <- PEcAn.utils::get.run.id("ENS", PEcAn.utils::left.pad.zeros(i, 
                                                                             5), site.id = settings$run$site$id)
      }
      runs[i, "id"] <- run.id
      if (clean) {
        unlink(file.path(settings$rundir, run.id))
        unlink(file.path(settings$modeloutdir, run.id))
      }
      dir.create(file.path(settings$rundir, run.id), recursive = TRUE)
      dir.create(file.path(settings$modeloutdir, run.id), 
                 recursive = TRUE)
      input_info <- ""
      for (input_i in seq_along(settings$run$inputs)) {
        input_tag <- names(settings$run$inputs)[[input_i]]
        if (!is.null(samples[[input_tag]])) {
          settings$run$inputs[[input_tag]][["path"]] <- samples[[input_tag]][["samples"]][[i]]
          input_info <- paste0(input_info, format(input_tag, 
                                                  width = 12, justify = "left"), ": ", samples[[input_tag]]$samples[[i]], 
                               "\n")
        }
      }
      # cat("runtype     : ensemble\n", "workflow id : ", 
      #     format(workflow.id, scientific = FALSE), "\n", 
      #     "ensemble id : ", format(ensemble.id, scientific = FALSE), 
      #     "\n", "run         : ", i, "/", settings$ensemble$size, 
      #     "\n", "run id      : ", format(run.id, scientific = FALSE), 
      #     "\n", "pft names   : ", as.character(lapply(settings$pfts, 
      #                                                 function(x) x[["name"]])), "\n", "model       : ", 
      #     model, "\n", "model id    : ", format(settings$model$id, 
      #                                           scientific = FALSE), "\n", "site        : ", 
      #     settings$run$site$name, "\n", "site  id    : ", 
      #     format(settings$run$site$id, scientific = FALSE), 
      #     "\n", input_info, "start date  : ", settings$run$start.date, 
      #     "\n", "end date    : ", settings$run$end.date, 
      #     "\n", "hostname    : ", settings$host$name, "\n", 
      #     "rundir      : ", file.path(settings$host$rundir, 
      #                                 run.id), "\n", "outdir      : ", file.path(settings$host$outdir, 
      #                                                                            run.id), "\n", file = file.path(settings$rundir, 
      # run.id, "README.txt"))
      for (input_i in seq_along(settings$run$inputs)) {
        input_tag <- names(settings$run$inputs)[[input_i]]
        input <- settings$run$inputs[[input_tag]]
        if (!input_tag %in% names(samples)) {
          settings$run$inputs[[input_tag]]$path <- input$path[1]
          #### See if it works
          # print(input$path[1])
        }
        else {
          settings$run$inputs[[input_tag]]$path <- samples[[input_tag]][["samples"]][[i]]
        }
      }
      ####
      # print("reading IC files")
      # # message("IC path for this run:")
      # print(settings$run$inputs$poolinitcond$path)
      
      ##### Test code
      # —— 放在 do.call(my.write.config, ...) 之前 —— 
      # p <- settings$run$inputs$poolinitcond$path
      # 
      # # 1) 彻底“清洁”这个路径
      # if (is.list(p)) p <- unlist(p, use.names = FALSE)
      # p <- as.character(p)
      # p <- p[1]
      # p <- trimws(p)
      # 
      # cat("IC CHECK | type:", paste(class(p), collapse=","),
      #     "| len:", length(p),
      #     "| nchar:", nchar(p), "\n")
      # cat("IC CHECK | path:", p, "\n")
      # cat("IC CHECK | exists:", file.exists(p),
      #     "| access(read):", tryCatch(file.access(p, 4), error = function(e) NA), "\n")
      # cat("IC CHECK | size(B):", tryCatch(file.info(p)$size, error=function(e) NA), "\n")
      # 
      # 2) 直接在这里试开 NetCDF（和 write.config 里用的是同一条路径）
      # ok <- tryCatch({
      #   nc <- ncdf4::nc_open(p); on.exit(ncdf4::nc_close(nc), add = TRUE)
      #   cat("IC CHECK | nc_open OK; vars:", paste(names(nc$var), collapse=","), "\n")
      #   TRUE
      # }, error = function(e) {
      #   cat("IC CHECK | nc_open FAIL:", conditionMessage(e), "\n")
      #   FALSE
      # })
      # 
      # # 3) 若需要，还可以提前试 prepare_pools（和里面一致）
      # pp <- tryCatch({
      #   PEcAn.data.land::prepare_pools(p, constants = list(sla = 10))  # sla 给个随便的正数
      # }, error = function(e) e)
      # if (inherits(pp, "error") || is.null(pp)) {
      #   cat("IC CHECK | prepare_pools FAIL\n")
      # } else {
      #   cat("IC CHECK | prepare_pools OK; names:", paste(names(pp), collapse=","), "\n")
      # }
      # 
      do.call(my.write.config, args = list(defaults = defaults, 
                                           trait.values = lapply(samples$parameters$samples, 
                                                                 function(x, n) {
                                                                   x[n, , drop = FALSE]
                                                                 }, n = i), settings = settings, run.id = run.id))
      # cat(format(run.id, scientific = FALSE), file = file.path(settings$rundir, 
      #                                                          "runs.txt"), sep = "\n", append = TRUE)
    }
    return(invisible(list(runs = runs, ensemble.id = ensemble.id, 
                          samples = samples)))
  }
  else {
    inputs <- restart$inputs
    run.id <- restart$runid
    new.params <- restart$new.params
    new.state <- restart$new.state
    ensemble.id <- restart$ensemble.id
    site.pfts.vec <- settings$run$site$site.pft %>% unlist %>% 
      as.character
    if (!is.null(site.pfts.vec)) {
      defined.pfts <- settings$pfts %>% purrr::map("name") %>% 
        unlist %>% as.character
      if (length(which(site.pfts.vec %in% defined.pfts)) > 
          0) 
        new.params <- new.params %>% purrr::map(~list(.x[[which(site.pfts.vec %in% 
                                                                  defined.pfts)]], restart = .x$restart))
      if (length(which(!(site.pfts.vec %in% defined.pfts))) > 
          0) 
        PEcAn.logger::logger.warn(paste0("The following pfts are specified for the siteid ", 
                                         settings$run$site$id, " but they are not defined as a pft in pecan.xml:", 
                                         site.pfts.vec[which(!(site.pfts.vec %in% defined.pfts))]))
    }
    for (j in 1:length(run.id)) {
      if (!file.exists(file.path(settings$rundir, run.id[[j]]))) {
        dir.create(file.path(settings$rundir, run.id[[j]]))
      }
    }
    for (i in seq_len(settings$ensemble$size)) {
      input_list <- list()
      for (input_tag in names(inputs)) {
        if (!is.null(inputs[[input_tag]]$samples[[i]])) 
          input_list[[input_tag]] <- list(path = inputs[[input_tag]]$samples[[i]])
      }
      do.call(my.write_restart, args = list(outdir = settings$host$outdir, 
                                            runid = run.id[[i]], start.time = restart$start.time, 
                                            stop.time = restart$stop.time, settings = settings, 
                                            new.state = new.state[i, ], new.params = new.params[[i]], 
                                            inputs = input_list, RENAME = rename))
    }
    params <- new.params
    return(invisible(list(runs = data.frame(id = run.id), 
                          ensemble.id = ensemble.id, samples = inputs)))
  }
}

write.config.SIPNET <- function (defaults, trait.values, settings, run.id, inputs = NULL, 
                                 IC = NULL, restart = NULL, spinup = NULL) 
{
  # print("11: new WCS function")
  template.in <- system.file("sipnet.in", package = "PEcAn.SIPNET")
  # print("checkpoint1")
  config.text <- readLines(con = template.in, n = -1)
  # print("checkpoint2")
  writeLines(config.text, con = file.path(settings$rundir, 
                                          run.id, "sipnet.in"))
  # print("checkpoint3")
  template.clim <- settings$run$inputs$met$path
  if (!is.null(inputs)) {
    if ("met" %in% names(inputs)) {
      template.clim <- inputs$met$path
    }
  }
  # print("checkpoint4")
  PEcAn.logger::logger.info(paste0("Writing SIPNET configs with input ", 
                                   template.clim))
  rundir <- file.path(settings$host$rundir, as.character(run.id))
  outdir <- file.path(settings$host$outdir, as.character(run.id))
  # print("checkpoint2")
  if (is.null(settings$host$qsub) && (settings$host$name == 
                                      "localhost")) {
    rundir <- file.path(settings$rundir, as.character(run.id))
    outdir <- file.path(settings$modeloutdir, as.character(run.id))
  }
  if (!is.null(settings$model$jobtemplate) && file.exists(settings$model$jobtemplate)) {
    jobsh <- readLines(con = settings$model$jobtemplate, 
                       n = -1)
  }
  else {
    jobsh <- readLines(con = system.file("template.job", 
                                         package = "PEcAn.SIPNET"), n = -1)
  }
  hostsetup <- ""
  if (!is.null(settings$model$prerun)) {
    hostsetup <- paste(hostsetup, sep = "\n", paste(settings$model$prerun, 
                                                    collapse = "\n"))
  }
  if (!is.null(settings$host$prerun)) {
    hostsetup <- paste(hostsetup, sep = "\n", paste(settings$host$prerun, 
                                                    collapse = "\n"))
  }
  cdosetup <- ""
  if (!is.null(settings$host$cdosetup)) {
    cdosetup <- paste(cdosetup, sep = "\n", paste(settings$host$cdosetup, 
                                                  collapse = "\n"))
  }
  hostteardown <- ""
  if (!is.null(settings$model$postrun)) {
    hostteardown <- paste(hostteardown, sep = "\n", paste(settings$model$postrun, 
                                                          collapse = "\n"))
  }
  if (!is.null(settings$host$postrun)) {
    hostteardown <- paste(hostteardown, sep = "\n", paste(settings$host$postrun, 
                                                          collapse = "\n"))
  }
  cpruncmd <- cpoutcmd <- rmoutdircmd <- rmrundircmd <- ""
  if (!is.null(settings$host$rabbitmq)) {
    settings$host$rabbitmq$cpfcmd <- ifelse(is.null(settings$host$rabbitmq$cpfcmd), 
                                            "", settings$host$rabbitmq$cpfcmd)
    cpruncmd <- gsub("@OUTDIR@", settings$host$rundir, settings$host$rabbitmq$cpfcmd)
    cpruncmd <- gsub("@OUTFOLDER@", rundir, cpruncmd)
    cpoutcmd <- gsub("@OUTDIR@", settings$host$outdir, settings$host$rabbitmq$cpfcmd)
    cpoutcmd <- gsub("@OUTFOLDER@", outdir, cpoutcmd)
    rmoutdircmd <- paste("rm", file.path(outdir, "*"))
    rmrundircmd <- paste("rm", file.path(rundir, "*"))
  }
  jobsh <- gsub("@HOST_SETUP@", hostsetup, jobsh)
  jobsh <- gsub("@CDO_SETUP@", cdosetup, jobsh)
  jobsh <- gsub("@HOST_TEARDOWN@", hostteardown, jobsh)
  jobsh <- gsub("@SITE_LAT@", settings$run$site$lat, jobsh)
  jobsh <- gsub("@SITE_LON@", settings$run$site$lon, jobsh)
  jobsh <- gsub("@SITE_MET@", template.clim, jobsh)
  jobsh <- gsub("@OUTDIR@", outdir, jobsh)
  jobsh <- gsub("@RUNDIR@", rundir, jobsh)
  jobsh <- gsub("@START_DATE@", settings$run$start.date, jobsh)
  jobsh <- gsub("@END_DATE@", settings$run$end.date, jobsh)
  jobsh <- gsub("@BINARY@", settings$model$binary, jobsh)
  jobsh <- gsub("@REVISION@", settings$model$revision, jobsh)
  jobsh <- gsub("@CPRUNCMD@", cpruncmd, jobsh)
  jobsh <- gsub("@CPOUTCMD@", cpoutcmd, jobsh)
  jobsh <- gsub("@RMOUTDIRCMD@", rmoutdircmd, jobsh)
  jobsh <- gsub("@RMRUNDIRCMD@", rmrundircmd, jobsh)
  
  # print("checkpoint5")
  
  if (is.null(settings$state.data.assimilation$NC.Prefix)) {
    settings$state.data.assimilation$NC.Prefix <- "sipnet.out"
  }
  # print("checkpoint6")
  jobsh <- gsub("@PREFIX@", settings$state.data.assimilation$NC.Prefix, 
                jobsh)
  if (is.null(settings$state.data.assimilation$NC.Overwrite)) {
    settings$state.data.assimilation$NC.Overwrite <- FALSE
  }
  # print("checkpoint7")
  jobsh <- gsub("@OVERWRITE@", settings$state.data.assimilation$NC.Overwrite, 
                jobsh)
  if (is.null(settings$state.data.assimilation$FullYearNC)) {
    settings$state.data.assimilation$FullYearNC <- FALSE
  }
  # print("checkpoint8")
  jobsh <- gsub("@CONFLICT@", settings$state.data.assimilation$FullYearNC, 
                jobsh)
  if (is.null(settings$model$delete.raw)) {
    settings$model$delete.raw <- FALSE
  }
  # print("checkpoint9")
  jobsh <- gsub("@DELETE.RAW@", settings$model$delete.raw, 
                jobsh)
  writeLines(jobsh, con = file.path(settings$rundir, run.id, 
                                    "job.sh"))
  Sys.chmod(file.path(settings$rundir, run.id, "job.sh"))
  template.paramSpatial <- system.file("template.param-spatial", 
                                       package = "PEcAn.SIPNET")
  file.copy(template.paramSpatial, file.path(settings$rundir, 
                                             run.id, "sipnet.param-spatial"))
  template.param <- system.file("template.param", package = "PEcAn.SIPNET")
  if ("default.param" %in% names(settings$model)) {
    template.param <- settings$model$default.param
  }
  # print("checkpoint10")
  param <- utils::read.table(template.param)
  trait_names_all_pfts <- as.vector(sapply(trait.values, names))
  dup_traitnames <- trait_names_all_pfts[duplicated(trait_names_all_pfts)]
  if (length(dup_traitnames) > 0) {
    PEcAn.logger::logger.warn("Multiple trait values given for parameters", 
                              paste(dQuote(dup_traitnames), collapse = ", "), "write.config.SIPNET will use the value it sees last.")
  }
  # print("checkpoint11")
  for (pft in seq_along(trait.values)) {
    pft.traits <- unlist(trait.values[[pft]])
    pft.names <- names(pft.traits)
    constant.traits <- unlist(defaults[[1]]$constants)
    constant.names <- names(constant.traits)
    for (i in seq_along(constant.traits)) {
      ind <- match(constant.names[i], pft.names)
      if (is.na(ind)) {
        pft.names <- c(pft.names, constant.names[i])
        pft.traits <- c(pft.traits, constant.traits[i])
      }
      else {
        pft.traits[ind] <- constant.traits[i]
      }
    }
    # print("checkpoint12")
    pft.names <- pft.names[pft.traits != "NA" & !is.na(pft.traits)]
    pft.traits <- pft.traits[pft.traits != "NA" & !is.na(pft.traits)]
    pft.traits <- as.numeric(pft.traits)
    leafC <- 0.48
    if ("leafC" %in% pft.names) {
      leafC <- pft.traits[which(pft.names == "leafC")]
      id <- which(param[, 1] == "cFracLeaf")
      param[id, 2] <- leafC * 0.01
    }
    SLA <- NA
    id <- which(param[, 1] == "leafCSpWt")
    if ("SLA" %in% pft.names) {
      SLA <- pft.traits[which(pft.names == "SLA")]
      param[id, 2] <- 1000 * leafC * 0.01/SLA
    }
    else {
      SLA <- 1000 * leafC/param[id, 2]
    }
    Amax <- NA
    id <- which(param[, 1] == "aMax")
    if ("Amax" %in% pft.names) {
      Amax <- pft.traits[which(pft.names == "Amax")]
      param[id, 2] <- Amax * SLA
    }
    else {
      Amax <- param[id, 2] * SLA
    }
    if ("AmaxFrac" %in% pft.names) {
      param[which(param[, 1] == "aMaxFrac"), 2] <- pft.traits[which(pft.names == 
                                                                      "AmaxFrac")]
    }
    if ("extinction_coefficient" %in% pft.names) {
      param[which(param[, 1] == "attenuation"), 2] <- pft.traits[which(pft.names == 
                                                                         "extinction_coefficient")]
    }
    if ("leaf_respiration_rate_m2" %in% pft.names) {
      Rd <- pft.traits[which(pft.names == "leaf_respiration_rate_m2")]
      id <- which(param[, 1] == "baseFolRespFrac")
      param[id, 2] <- max(min(Rd/Amax, 1), 0)
    }
    if ("Vm_low_temp" %in% pft.names) {
      param[which(param[, 1] == "psnTMin"), 2] <- pft.traits[which(pft.names == 
                                                                     "Vm_low_temp")]
    }
    if ("psnTOpt" %in% pft.names) {
      param[which(param[, 1] == "psnTOpt"), 2] <- pft.traits[which(pft.names == 
                                                                     "psnTOpt")]
    }
    if ("growth_resp_factor" %in% pft.names) {
      param[which(param[, 1] == "growthRespFrac"), 2] <- pft.traits[which(pft.names == 
                                                                            "growth_resp_factor")]
    }
    if ("half_saturation_PAR" %in% pft.names) {
      param[which(param[, 1] == "halfSatPar"), 2] <- pft.traits[which(pft.names == 
                                                                        "half_saturation_PAR")]
    }
    if ("stomatal_slope.BB" %in% pft.names) {
      id <- which(param[, 1] == "m_ballBerry")
      param[id, 2] <- pft.traits[which(pft.names == "stomatal_slope.BB")]
    }
    if ("dVPDSlope" %in% pft.names) {
      param[which(param[, 1] == "dVpdSlope"), 2] <- pft.traits[which(pft.names == 
                                                                       "dVPDSlope")]
    }
    if ("dVpdExp" %in% pft.names) {
      param[which(param[, 1] == "dVpdExp"), 2] <- pft.traits[which(pft.names == 
                                                                     "dVpdExp")]
    }
    if ("leaf_turnover_rate" %in% pft.names) {
      param[which(param[, 1] == "leafTurnoverRate"), 2] <- pft.traits[which(pft.names == 
                                                                              "leaf_turnover_rate")]
    }
    if ("wueConst" %in% pft.names) {
      param[which(param[, 1] == "wueConst"), 2] <- pft.traits[which(pft.names == 
                                                                      "wueConst")]
    }
    if ("veg_respiration_Q10" %in% pft.names) {
      param[which(param[, 1] == "vegRespQ10"), 2] <- pft.traits[which(pft.names == 
                                                                        "veg_respiration_Q10")]
    }
    if ("stem_respiration_rate" %in% pft.names) {
      vegRespQ10 <- param[which(param[, 1] == "vegRespQ10"), 
                          2]
      id <- which(param[, 1] == "baseVegResp")
      stem_resp_g <- (((pft.traits[which(pft.names == "stem_respiration_rate")]) * 
                         (44.0096/1e+06) * (12.01/44.0096))/1000) * 86400
      param[id, 2] <- stem_resp_g * vegRespQ10^(-25/10)
    }
    if ("root_turnover_rate" %in% pft.names) {
      id <- which(param[, 1] == "fineRootTurnoverRate")
      param[id, 2] <- pft.traits[which(pft.names == "root_turnover_rate")]
    }
    if ("fine_root_respiration_Q10" %in% pft.names) {
      param[which(param[, 1] == "fineRootQ10"), 2] <- pft.traits[which(pft.names == 
                                                                         "fine_root_respiration_Q10")]
    }
    if ("root_respiration_rate" %in% pft.names) {
      fineRootQ10 <- param[which(param[, 1] == "fineRootQ10"), 
                           2]
      id <- which(param[, 1] == "baseFineRootResp")
      root_resp_rate_g <- (((pft.traits[which(pft.names == 
                                                "root_respiration_rate")]) * (44.0096/1e+06) * 
                              (12.01/44.0096))/1000) * 86400
      param[id, 2] <- root_resp_rate_g * fineRootQ10^(-25/10)
    }
    if ("coarse_root_respiration_Q10" %in% pft.names) {
      param[which(param[, 1] == "coarseRootQ10"), 2] <- pft.traits[which(pft.names == 
                                                                           "coarse_root_respiration_Q10")]
    }
    alloc_params <- c("root_allocation_fraction", "wood_allocation_fraction", 
                      "leaf_allocation_fraction")
    if (all(alloc_params %in% pft.names)) {
      sum_alloc <- pft.traits[which(pft.names == "root_allocation_fraction")] + 
        pft.traits[which(pft.names == "wood_allocation_fraction")] + 
        pft.traits[which(pft.names == "leaf_allocation_fraction")]
      if (sum_alloc > 1) {
        PEcAn.logger::logger.warn("Sum of allocation parameters exceeds 1 for runid = ", 
                                  run.id, "- This won't break anything since SIPNET has internal check, but notice that such combinations might not take effect in the outputs.")
      }
    }
    if ("root_allocation_fraction" %in% pft.names) {
      param[which(param[, 1] == "fineRootAllocation"), 
            2] <- pft.traits[which(pft.names == "root_allocation_fraction")]
    }
    if ("wood_allocation_fraction" %in% pft.names) {
      param[which(param[, 1] == "woodAllocation"), 2] <- pft.traits[which(pft.names == 
                                                                            "wood_allocation_fraction")]
    }
    if ("leaf_allocation_fraction" %in% pft.names) {
      param[which(param[, 1] == "leafAllocation"), 2] <- pft.traits[which(pft.names == 
                                                                            "leaf_allocation_fraction")]
    }
    if ("wood_turnover_rate" %in% pft.names) {
      param[which(param[, 1] == "woodTurnoverRate"), 2] <- pft.traits[which(pft.names == 
                                                                              "wood_turnover_rate")]
    }
    if ("soil_respiration_Q10" %in% pft.names) {
      param[which(param[, 1] == "soilRespQ10"), 2] <- pft.traits[which(pft.names == 
                                                                         "soil_respiration_Q10")]
    }
    if ("som_respiration_rate" %in% pft.names) {
      param[which(param[, 1] == "baseSoilResp"), 2] <- pft.traits[which(pft.names == 
                                                                          "som_respiration_rate")]
    }
    if ("turn_over_time" %in% pft.names) {
      id <- which(param[, 1] == "litterBreakdownRate")
      param[id, 2] <- pft.traits[which(pft.names == "turn_over_time")]
    }
    if ("frozenSoilEff" %in% pft.names) {
      param[which(param[, 1] == "frozenSoilEff"), 2] <- pft.traits[which(pft.names == 
                                                                           "frozenSoilEff")]
    }
    if ("frozenSoilFolREff" %in% pft.names) {
      param[which(param[, 1] == "frozenSoilFolREff"), 2] <- pft.traits[which(pft.names == 
                                                                               "frozenSoilFolREff")]
    }
    if ("soilWHC" %in% pft.names) {
      param[which(param[, 1] == "soilWHC"), 2] <- pft.traits[which(pft.names == 
                                                                     "soilWHC")]
    }
    if (TRUE) {
      param[which(param[, 1] == "soilRespQ10Cold"), 2] <- param[which(param[, 
                                                                            1] == "soilRespQ10"), 2]
      param[which(param[, 1] == "baseSoilRespCold"), 2] <- param[which(param[, 
                                                                             1] == "baseSoilResp"), 2] * 0.25
    }
    if ("immedEvapFrac" %in% pft.names) {
      param[which(param[, 1] == "immedEvapFrac"), 2] <- pft.traits[which(pft.names == 
                                                                           "immedEvapFrac")]
    }
    if ("leafWHC" %in% pft.names) {
      param[which(param[, 1] == "leafPoolDepth"), 2] <- pft.traits[which(pft.names == 
                                                                           "leafWHC")]
    }
    if ("waterRemoveFrac" %in% pft.names) {
      param[which(param[, 1] == "waterRemoveFrac"), 2] <- pft.traits[which(pft.names == 
                                                                             "waterRemoveFrac")]
    }
    if ("fastFlowFrac" %in% pft.names) {
      param[which(param[, 1] == "fastFlowFrac"), 2] <- pft.traits[which(pft.names == 
                                                                          "fastFlowFrac")]
    }
    if ("rdConst" %in% pft.names) {
      param[which(param[, 1] == "rdConst"), 2] <- pft.traits[which(pft.names == 
                                                                     "rdConst")]
    }
    if ("GDD" %in% pft.names) {
      param[which(param[, 1] == "gddLeafOn"), 2] <- pft.traits[which(pft.names == 
                                                                       "GDD")]
    }
    if ("fracLeafFall" %in% pft.names) {
      param[which(param[, 1] == "fracLeafFall"), 2] <- pft.traits[which(pft.names == 
                                                                          "fracLeafFall")]
    }
    if ("leafGrowth" %in% pft.names) {
      param[which(param[, 1] == "leafGrowth"), 2] <- pft.traits[which(pft.names == 
                                                                        "leafGrowth")]
    }
    if (!is.null(settings$run$inputs$leaf_phenology)) {
      obs_year_start <- lubridate::year(settings$run$start.date)
      obs_year_end <- lubridate::year(settings$run$end.date)
      # if (obs_year_start != obs_year_end) {
      #   PEcAn.logger::logger.info("Start.date and end.date are not in the same year.", 
      #                             "Using phenological data from start year only.")
      # }
      leaf_pheno_path <- settings$run$inputs$leaf_phenology$path
      if (!is.null(leaf_pheno_path)) {
        leafphdata <- utils::read.csv(leaf_pheno_path)
        leafOnDay <- leafphdata$leafonday[leafphdata$year == 
                                            obs_year_start & leafphdata$site_id == settings$run$site$id]
        leafOffDay <- leafphdata$leafoffday[leafphdata$year == 
                                              obs_year_start & leafphdata$site_id == settings$run$site$id]
        if (!is.na(leafOnDay)) {
          param[which(param[, 1] == "leafOnDay"), 2] <- leafOnDay
        }
        if (!is.na(leafOffDay)) {
          param[which(param[, 1] == "leafOffDay"), 2] <- leafOffDay
        }
      }
      else {
        PEcAn.logger::logger.info("No phenology data were found.", 
                                  "Please consider running `PEcAn.data.remote::extract_phenology_MODIS`", 
                                  "to get the parameter file.")
      }
    }
  }
  # print("checkpoint13")
  if (length(settings$run$inputs$soil_physics$path) > 0) {
    template.soil_physics <- settings$run$inputs$soil_physics$path
    if (!is.null(inputs)) {
      if ("soil_physics" %in% names(inputs)) {
        template.soil_physics <- inputs$soil_physics$path
      }
    }
    # print("breakpoint1")
    # print(length(template.soil_physics))
    if (length(template.soil_physics) != 1) {
      PEcAn.logger::logger.warn(paste0("No single soil physical parameter file was found for ", 
                                       run.id))
    }
    else {
      # print("breakpoint1.5")
      # print("template.soil_physics: ")
      # print(template.soil_physics)
      
      # check the type of template.soil_physics
      # print(class(template.soil_physics))
      # Change from list to character
      template.soil_physics <- as.character(template.soil_physics)
      soil_IC_list <- PEcAn.data.land::pool_ic_netcdf2list(template.soil_physics)
      # print("breakpoint2")
      if ("volume_fraction_of_water_in_soil_at_saturation" %in% 
          names(soil_IC_list$vals)) {
        if ("depth" %in% names(soil_IC_list$dims)) {
          thickness <- c(soil_IC_list$dims$depth[1], 
                         diff(soil_IC_list$dims$depth))
          thickness <- PEcAn.utils::ud_convert(thickness, 
                                               "m", "cm")
          soilWHC_total <- sum(unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"]) * 
                                 thickness)
          if (thickness[1] <= 10) {
            param[which(param[, 1] == "litterWHC"), 2] <- unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"])[1] * 
              thickness[1]
          }
        }
        else {
          PEcAn.logger::logger.warn("No depth info was found in the soil file. Will use the default or user-specified soil depth")
          thickness <- 100
          soilWHC_total <- soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"] * 
            thickness
        }
        param[which(param[, 1] == "soilWHC"), 2] <- soilWHC_total
      }
      if ("soil_hydraulic_conductivity_at_saturation" %in% 
          names(soil_IC_list$vals)) {
        param[which(param[, 1] == "litWaterDrainRate"), 
              2] <- PEcAn.utils::ud_convert(unlist(soil_IC_list$vals["soil_hydraulic_conductivity_at_saturation"])[1], 
                                            "m s-1", "cm day-1")
      }
    }
  }
  if (!is.null(IC)) {
    ic.names <- names(IC)
    plant_wood_vars <- c("AbvGrndWood", "abvGrndWoodFrac", 
                         "coarseRootFrac", "fineRootFrac")
    if (all(plant_wood_vars %in% ic.names)) {
      # print("breakpoint3")
      if (IC$abvGrndWoodFrac < 0.05) {
        wood_total_C <- IC$AbvGrndWood
      }
      else {
        wood_total_C <- IC$AbvGrndWood/IC$abvGrndWoodFrac
      }
      if (is.infinite(wood_total_C) | is.nan(wood_total_C) | 
          wood_total_C < 0) {
        wood_total_C <- 0
        if (round(IC$AbvGrndWood) > 0 & round(IC$abvGrndWoodFrac, 
                                              3) == 0) 
          PEcAn.logger::logger.warn(paste0("There is a major problem with ", 
                                           run.id, " in either the model's parameters or IC.", 
                                           "Because the ABG is estimated=", IC$AbvGrndWood, 
                                           " while AGB Frac is estimated=", IC$abvGrndWoodFrac))
      }
      param[which(param[, 1] == "plantWoodInit"), 2] <- wood_total_C
      param[which(param[, 1] == "coarseRootFrac"), 2] <- IC$coarseRootFrac
      param[which(param[, 1] == "fineRootFrac"), 2] <- IC$fineRootFrac
    }
    if ("lai" %in% ic.names) {
      param[which(param[, 1] == "laiInit"), 2] <- IC$lai
    }
    if ("litter_carbon_content" %in% ic.names) {
      param[which(param[, 1] == "litterInit"), 2] <- IC$litter_carbon_content
    }
    if ("soil" %in% ic.names) {
      param[which(param[, 1] == "soilInit"), 2] <- IC$soil
    }
    if ("litter_mass_content_of_water" %in% ic.names) {
      param[which(param[, 1] == "litterWFracInit"), 2] <- IC$litter_mass_content_of_water/(param[which(param[, 
                                                                                                             1] == "litterWHC"), 2] * 10)
    }
    if ("soilWater" %in% ic.names) {
      param[which(param[, 1] == "soilWFracInit"), 2] <- IC$soilWater/(param[which(param[, 
                                                                                        1] == "soilWHC"), 2] * 10)
    }
    if ("soilWFrac" %in% ic.names) {
      param[which(param[, 1] == "soilWFracInit"), 2] <- IC$soilWFrac
    }
    if ("SWE" %in% ic.names) {
      param[which(param[, 1] == "snowInit"), 2] <- IC$SWE
    }
    if ("microbe" %in% ic.names) {
      param[which(param[, 1] == "microbeInit"), 2] <- IC$microbe
    }
  }
  else if (length(settings$run$inputs$poolinitcond$path) > 
           0) {
    IC.path <- settings$run$inputs$poolinitcond$path
    if (length(IC.path) > 1) {
      PEcAn.logger::logger.error("write.config.SIPNET needs one poolinitcond path", 
                                 "got", length(IC.path))
    }
    IC.pools <- PEcAn.data.land::prepare_pools(IC.path, constants = list(sla = SLA))
    if (!is.null(IC.pools)) {
      IC.nc <- ncdf4::nc_open(IC.path)
      ic_ncvars_to_try <- c("nee", "SoilMoistFrac", "SWE", 
                            "date_of_budburst", "date_of_senescence", "Microbial Biomass C")
      ic_has_ncvars <- ic_ncvars_to_try %in% names(IC.nc$var)
      names(ic_has_ncvars) <- ic_ncvars_to_try
      if ("wood" %in% names(IC.pools)) {
        param[param[, 1] == "plantWoodInit", 2] <- PEcAn.utils::ud_convert(IC.pools$wood, 
                                                                           "kg m-2", "g m-2")
      }
      lai <- IC.pools$LAI
      if (!is.na(lai) && is.numeric(lai)) {
        param[param[, 1] == "laiInit", 2] <- lai
      }
      is_deciduous_pft <- isTRUE(param[param[, 1] == "fracLeafFall", 
                                       2] > 0.5)
      start_day <- lubridate::yday(settings$run$start.date)
      starts_with_leaves <- (start_day >= param[param[, 
                                                      1] == "leafOnDay", 2] && start_day <= param[param[, 
                                                                                                        1] == "leafOffDay", 2])
      if (is_deciduous_pft && !starts_with_leaves) {
        param[param[, 1] == "laiInit", 2] <- 0
      }
      if (ic_has_ncvars[["nee"]]) {
        nee <- ncdf4::ncvar_get(IC.nc, "nee")
        if (!is.na(nee) && is.numeric(nee)) {
          param[param[, 1] == "neeInit", 2] <- nee
        }
      }
      if ("litter" %in% names(IC.pools)) {
        param[param[, 1] == "litterInit", 2] <- PEcAn.utils::ud_convert(IC.pools$litter, 
                                                                        "g m-2", "g m-2")
      }
      if ("soil" %in% names(IC.pools)) {
        param[param[, 1] == "soilInit", 2] <- PEcAn.utils::ud_convert(sum(IC.pools$soil), 
                                                                      "kg m-2", "g m-2")
      }
      if (ic_has_ncvars[["SoilMoistFrac"]]) {
        soilWFrac <- ncdf4::ncvar_get(IC.nc, "SoilMoistFrac")
        if (!is.na(soilWFrac) && is.numeric(soilWFrac)) {
          param[param[, 1] == "soilWFracInit", 2] <- sum(soilWFrac)/100
          litterWFrac <- soilWFrac
        }
      }
      if (ic_has_ncvars[["SWE"]]) {
        snow <- ncdf4::ncvar_get(IC.nc, "SWE")
        if (!is.na(snow) && is.numeric(snow)) {
          param[param[, 1] == "snowInit", 2] <- PEcAn.utils::ud_convert(snow, 
                                                                        "kg m-2", "g cm-2")
        }
      }
      if (ic_has_ncvars[["date_of_budburst"]]) {
        leafOnDay <- ncdf4::ncvar_get(IC.nc, "date_of_budburst")
        if (!is.na(leafOnDay) && is.numeric(leafOnDay)) {
          param[param[, 1] == "leafOnDay", 2] <- leafOnDay
        }
      }
      if (ic_has_ncvars[["date_of_senescence"]]) {
        leafOffDay <- ncdf4::ncvar_get(IC.nc, "date_of_senescence")
        if (!is.na(leafOffDay) && is.numeric(leafOffDay)) {
          param[param[, 1] == "leafOffDay", 2] <- leafOffDay
        }
      }
      if (ic_has_ncvars[["Microbial Biomass C"]]) {
        microbe <- ncdf4::ncvar_get(IC.nc, "Microbial Biomass C")
        if (!is.na(microbe) && is.numeric(microbe)) {
          param[param[, 1] == "microbeInit", 2] <- PEcAn.utils::ud_convert(microbe, 
                                                                           "mg kg-1", "mg g-1")
        }
      }
      ncdf4::nc_close(IC.nc)
    }
    else {
      PEcAn.logger::logger.error("Bad initial conditions filepath; keeping defaults")
    }
  }
  else {
  }
  if (!is.null(settings$run$inputs$soilmoisture)) {
    if (!is.null(settings$run$inputs$soilmoisture$path)) {
      soil.path <- settings$run$inputs$soilmoisture$path
      soilWFrac <- ncdf4::ncvar_get(ncdf4::nc_open(soil.path), 
                                    varid = "mass_fraction_of_unfrozen_water_in_soil_moisture")
      param[which(param[, 1] == "soilWFracInit"), 2] <- soilWFrac
    }
  }
  if (file.exists(file.path(settings$rundir, run.id, "sipnet.param"))) {
    file.rename(file.path(settings$rundir, run.id, "sipnet.param"), 
                file.path(settings$rundir, run.id, paste0("sipnet_", 
                                                          lubridate::year(settings$run$start.date), "_", 
                                                          lubridate::year(settings$run$end.date), ".param")))
  }
  utils::write.table(param, file.path(settings$rundir, run.id, 
                                      "sipnet.param"), row.names = FALSE, col.names = FALSE, 
                     quote = FALSE)
  # print("end of write.ensemble.config")
}



#' build_X
#' 
#' @name build_X
#' @author Alexis Helgeson
#' 
#' @description builds X matrix for SDA
#'
#' @param new.params object created from sda_matchparam, passed from sda.enkf_MultiSite
#' @param nens number of ensemble members i.e. runs
#' @param read_restart_times passed from sda.enkf_MultiSite
#' @param settings settings object, passed from sda.enkf_MultiSite
#' @param outdir location of previous run output folder containing .nc files
#' @param out.configs object created for build_X passed from sda.enkf_MultiSite
#' @param t Default t=1, for function to work within time loop
#' @param var.names list of state variables taken from settings object
#' @param my.read_restart object that points to the model restart function i.e. read_restart.SIPNET
#' @param restart_flag flag if it's a restart stage. Default is FALSE.
#'
#' @return X ready to be passed to SDA Analysis code
build_X <- function(out.configs, settings, new.params, nens, read_restart_times, outdir, t = 1, var.names, my.read_restart, restart_flag = FALSE){
  
  # print("Start buildx")
  
  if(t == 1 & restart_flag){
    
    # print("checkpoint1")
    
    reads <-
      furrr::future_pmap(list(out.configs %>% `class<-`(c("list")), settings, new.params),function(configs,my_settings,siteparams) {
        # Loading the model package - this is required bc of the furrr
        #library(paste0("PEcAn.",settings$model$type), character.only = TRUE)
        #source("~/pecan/models/sipnet/R/read_restart.SIPNET.R")
        
        X_tmp <- vector("list", 2)
        
        # print("checkpoint2")
        
        
        for (i in seq_len(nens)) {
          X_tmp[[i]] <- do.call( my.read_restart,
                                 args = list(
                                   outdir = outdir,
                                   runid = my_settings$run$id[i] %>% as.character(),
                                   stop.time = read_restart_times[t+1],
                                   settings = my_settings,
                                   var.names = var.names,
                                   params = siteparams[[i]]
                                 )
          )
          
        }
        return(X_tmp)
      })
    
  }else{
    reads <-
      furrr::future_pmap(list(out.configs %>% `class<-`(c("list")), settings, new.params),function(configs,my_settings,siteparams) {
        
        X_tmp <- vector("list", 2)
        
        # print("checkpoint3")
        ####### Modify the output that contains nc files
        
        for (i in seq_len(nens)) {
          # print("loop begins")
          # print(configs$runs$id[i] %>% as.character())
          X_tmp[[i]] <- do.call( my.read_restart,
                                 args = list(
                                   outdir = outdir,
                                   runid = configs$runs$id[i] %>% as.character(),
                                   stop.time = read_restart_times[t+1],
                                   var.names = var.names,
                                   params = siteparams[[i]]
                                 )
          )
          # print("successful loop")
        }
        return(X_tmp)
      })
  }
  return(reads)
  # print("End of BuildX")
}



read.restart.SIPNET <- function (outdir, runid, stop.time, settings, var.names, params) 
{
  # print("Start read_restart.SIPNET")
  
  prior.sla <- params[[which(!names(params) %in% c("soil", 
                                                   "soil_SDA", "restart"))[1]]]$SLA
  forecast <- list()
  params$restart <- c()
  state.vars <- c("SWE", "SoilMoist", "SoilMoistFrac", "AbvGrndWood", "NEE","Qle",
                  "TotSoilCarb", "LAI", "litter_carbon_content", "fine_root_carbon_content", 
                  "coarse_root_carbon_content", "litter_mass_content_of_water")
  # print("breakpoint1")
  
  params$restart <- rep(NA, length(setdiff(state.vars, var.names)))
  names(params$restart) <- setdiff(state.vars, var.names)
  ens <- read.output(runid = runid, outdir = file.path(outdir, 
                                                       runid), start.year = lubridate::year(stop.time), end.year = lubridate::year(stop.time), 
                     variables = c(state.vars, "time_bounds"))
  
  #### See if ens is okay
  PEcAn.logger::logger.info(paste("read.output OK for runid:", runid))
  
  # 1) 打印变量名、每个变量长度（NULL 记为 NA）
  # cat("ENS names: ", paste(names(ens), collapse=", "), "\n")
  # len_or_na <- function(x) if (is.null(x)) NA_integer_ else length(x)
  # print(sapply(ens, len_or_na))
  
  # 2) 重点变量的前几项（如果存在）
  # peek <- function(x, n=5) if (is.null(x)) NULL else utils::head(x, n)
  # cat("Head(AbvGrndWood):\n"); print(peek(ens$AbvGrndWood))
  # cat("Head(fine_root_carbon_content):\n"); print(peek(ens$fine_root_carbon_content))
  # cat("Head(coarse_root_carbon_content):\n"); print(peek(ens$coarse_root_carbon_content))
  # cat("Head(LAI):\n"); print(peek(ens$LAI))
  
  # 3) time_bounds 维度和示例
  # if (!is.null(ens$time_bounds)) {
  #   cat("time_bounds dim: ", paste(dim(ens$time_bounds), collapse=" x "), "\n")
  #   cat("time_bounds[1, 1:6]:\n"); print(ens$time_bounds[1, 1:min(6, ncol(ens$time_bounds))])
  # } else {
  #   cat("time_bounds is NULL\n")
  # }
  
  # 4) 若关键变量缺失，直接给出可读的报错信息（避免在 ud_convert 里才崩）
  must_have <- c("AbvGrndWood","fine_root_carbon_content","coarse_root_carbon_content")
  missing <- must_have[ sapply(ens[must_have], is.null) ]
  if (length(missing) > 0) {
    stop("Missing variables in NetCDF for runid ", runid, ": ",
         paste(missing, collapse=", "),
         " | Available: ", paste(names(ens), collapse=", "),
         call. = FALSE)
  }
  
  start.time <- as.Date(paste0(lubridate::year(stop.time), 
                               "-01-01"))
  time_var <- ens$time_bounds[1, ]
  real_time <- as.POSIXct(time_var * 3600 * 24, origin = start.time)
  last <- which(as.Date(real_time) == as.Date(stop.time))[length(which(as.Date(real_time) == 
                                                                         as.Date(stop.time)))]
  
  # print("breakpoint3")
  # print(list(
  #   last = last,
  #   has_AbvGrndWood = !is.null(ens$AbvGrndWood),
  #   len_AbvGrndWood = if (is.null(ens$AbvGrndWood)) NA_integer_ else length(ens$AbvGrndWood),
  #   has_fine = !is.null(ens$fine_root_carbon_content),
  #   len_fine = if (is.null(ens$fine_root_carbon_content)) NA_integer_ else length(ens$fine_root_carbon_content),
  #   has_coarse = !is.null(ens$coarse_root_carbon_content),
  #   len_coarse = if (is.null(ens$coarse_root_carbon_content)) NA_integer_ else length(ens$coarse_root_carbon_content)
  # ))
  
  if (is.null(ens$AbvGrndWood)) {
    stop("ens$AbvGrndWood is NULL. Available vars: ", paste(names(ens), collapse=", "), call. = FALSE)
  }
  if (is.na(last) || length(ens$AbvGrndWood) < last) {
    stop("Index 'last' invalid for AbvGrndWood: last=", last,
         " len=", length(ens$AbvGrndWood), call. = FALSE)
  }
  
  if ("AbvGrndWood" %in% var.names) {
    # print("breakpoint3.1")
    forecast[[length(forecast) + 1]] <- PEcAn.utils::ud_convert(ens$AbvGrndWood[last], 
                                                                "kg/m^2", "Mg/ha")
    names(forecast[[length(forecast)]]) <- c("AbvGrndWood")
    wood_total_C <- ens$AbvGrndWood[last] + ens$fine_root_carbon_content[last] + 
      ens$coarse_root_carbon_content[last]
    if (wood_total_C <= 0) 
      wood_total_C <- 1e-04
    params$restart["abvGrndWoodFrac"] <- ens$AbvGrndWood[last]/wood_total_C
    params$restart["coarseRootFrac"] <- ens$coarse_root_carbon_content[last]/wood_total_C
    params$restart["fineRootFrac"] <- ens$fine_root_carbon_content[last]/wood_total_C
  }
  else {
    params$restart["AbvGrndWood"] <- PEcAn.utils::ud_convert(ens$AbvGrndWood[last], 
                                                             "kg/m^2", "g/m^2")
    wood_total_C <- ens$AbvGrndWood[last] + ens$fine_root_carbon_content[last] + 
      ens$coarse_root_carbon_content[last]
    if (wood_total_C <= 0) 
      wood_total_C <- 1e-04
    params$restart["abvGrndWoodFrac"] <- ens$AbvGrndWood[last]/wood_total_C
    params$restart["coarseRootFrac"] <- ens$coarse_root_carbon_content[last]/wood_total_C
    params$restart["fineRootFrac"] <- ens$fine_root_carbon_content[last]/wood_total_C
  }
  
  if ("GWBI" %in% var.names) {
    forecast[[length(forecast) + 1]] <- PEcAn.utils::ud_convert(mean(ens$GWBI), 
                                                                "kg/m^2/s", "Mg/ha/yr")
    names(forecast[[length(forecast)]]) <- c("GWBI")
  }
  if ("NEE" %in% var.names) {
    forecast[[length(forecast) + 1]] <-
      nee_model_to_obs(mean(ens$NEE))
    names(forecast[[length(forecast)]]) <- "NEE"
  }
  
  if ("Qle" %in% var.names) {
    forecast[[length(forecast) + 1]] <- mean(ens$Qle)
    names(forecast[[length(forecast)]]) <- c("Qle")
  }
  if ("leaf_carbon_content" %in% var.names) {
    forecast[[length(forecast) + 1]] <- ens$leaf_carbon_content[last]
    names(forecast[[length(forecast)]]) <- c("LeafC")
  }
  if ("LAI" %in% var.names) {
    forecast[[length(forecast) + 1]] <- ens$LAI[last]
    names(forecast[[length(forecast)]]) <- c("LAI")
  }
  else {
    params$restart["LAI"] <- ens$LAI[last]
  }
  if ("litter_carbon_content" %in% var.names) {
    forecast[[length(forecast) + 1]] <- ens$litter_carbon_content[last]
    names(forecast[[length(forecast)]]) <- c("litter_carbon_content")
  }
  else {
    params$restart["litter_carbon_content"] <- PEcAn.utils::ud_convert(ens$litter_carbon_content[last], 
                                                                       "kg m-2", "g m-2")
  }
  if ("litter_mass_content_of_water" %in% var.names) {
    forecast[[length(forecast) + 1]] <- ens$litter_mass_content_of_water[last]
    names(forecast[[length(forecast)]]) <- c("litter_mass_content_of_water")
  }
  else {
    params$restart["litter_mass_content_of_water"] <- ens$litter_mass_content_of_water[last]
  }
  if ("SoilMoist" %in% var.names) {
    forecast[[length(forecast) + 1]] <- ens$SoilMoist[last]
    names(forecast[[length(forecast)]]) <- c("SoilMoist")
  }
  else {
    params$restart["SoilMoist"] <- ens$SoilMoist[last]
  }
  if ("SoilMoistFrac" %in% var.names) {
    forecast[[length(forecast) + 1]] <- ens$SoilMoistFrac[last] * 
      100
    names(forecast[[length(forecast)]]) <- c("SoilMoistFrac")
  }
  else {
    params$restart["SoilMoistFrac"] <- ens$SoilMoistFrac[last]
  }
  if ("SWE" %in% var.names) {
    forecast[[length(forecast) + 1]] <- ens$SWE[last]
    names(forecast[[length(forecast)]]) <- c("SWE")
  }
  else {
    params$restart["SWE"] <- ens$SWE[last]/10
  }
  if ("TotLivBiom" %in% var.names) {
    forecast[[length(forecast) + 1]] <- PEcAn.utils::ud_convert(ens$TotLivBiom[last], 
                                                                "kg/m^2", "Mg/ha")
    names(forecast[[length(forecast)]]) <- c("TotLivBiom")
  }
  if ("TotSoilCarb" %in% var.names) {
    forecast[[length(forecast) + 1]] <- ens$TotSoilCarb[last]
    names(forecast[[length(forecast)]]) <- c("TotSoilCarb")
  }
  else {
    params$restart["TotSoilCarb"] <- PEcAn.utils::ud_convert(ens$TotSoilCarb[last], 
                                                             "kg m-2", "g m-2")
  }
  params$restart <- stats::na.omit(params$restart)
  # print(runid)
  X_tmp <- list(X = unlist(forecast), params = params)
  ### Check
  # print(X_tmp)
  return(X_tmp)
  # print("end of write.config.SIPNET")
}



read.output <- function (runid, outdir, start.year = NA, end.year = NA, variables = "GPP", 
                         dataframe = FALSE, pft.name = NULL, ncfiles = NULL, verbose = FALSE, 
                         print_summary = TRUE) 
  
{
  
  ##### Modify this specific outdir in read.output to read those nc files
  
  # print("Start debugging ens in read.output")
  # print(runid)
  # print(outdir)
  # print(variables)
  
  if ((missing(runid) || missing(outdir)) && is.null(ncfiles)) {
    PEcAn.logger::logger.severe("`runid` or `outdir` is missing, and `ncfiles` is NULL.", 
                                "Either provide both `runid` and `outdir`, or explicitly specify all `ncfiles`.")
  }
  if (!is.null(ncfiles)) {
    if (missing(runid)) 
      runid <- NULL
    if (missing(outdir)) 
      outdir <- NULL
  }
  if (is.null(ncfiles)) {
    ncfiles_sub <- list.files(path = outdir, pattern = "^-?[[:digit:]]{4}\\.nc$", 
                              full.names = FALSE)
    ncfiles <- file.path(outdir, ncfiles_sub)
  }
  else {
    ncfiles_sub <- basename(ncfiles)
  }
  if (!is.na(start.year)) {
    if (lubridate::is.instant(start.year)) {
      start.year <- lubridate::year(start.year)
    }
    else if (is.character(start.year)) {
      start.year <- as.numeric(start.year)
    }
    else if (is.numeric(start.year)) {
      if (start.year%%1 != 0) {
        PEcAn.logger::logger.severe("Start year `", start.year, 
                                    "` is numeric, but not an integer.")
      }
    }
    else {
      PEcAn.logger::logger.severe("`start.year` must be of type numeric, character, Date, or POSIXt", 
                                  "but given `", start.year, "` which is type `", 
                                  typeof(start.year), "`.")
    }
  }
  if (!is.na(end.year)) {
    if (lubridate::is.instant(end.year)) {
      end.year <- lubridate::year(end.year)
    }
    else if (is.character(end.year)) {
      end.year <- as.numeric(end.year)
    }
    else if (is.numeric(end.year)) {
      if (end.year%%1 != 0) {
        PEcAn.logger::logger.severe("End year `", end.year, 
                                    "` is numeric, but not an integer.")
      }
    }
    else {
      PEcAn.logger::logger.error("`end.year` must be of type numeric, character, Date, or POSIXt", 
                                 "but given `", end.year, "` which is type `", 
                                 typeof(end.year), "`.")
    }
  }
  nc_years <- suppressWarnings(as.numeric(gsub("^(-?[[:digit:]]{4})\\.nc", 
                                               "\\1", ncfiles_sub)))
  if (any(is.na(nc_years))) {
    PEcAn.logger::logger.debug("Unable to deduce NetCDF file years from their names. ", 
                               "Setting `nc_years` to length 0 numeric vector.")
    nc_years <- numeric()
  }
  if (is.na(start.year)) {
    PEcAn.logger::logger.debug("Missing start year.")
    if (length(nc_years) != 0) 
      start.year <- min(nc_years)
  }
  if (is.na(end.year)) {
    PEcAn.logger::logger.debug("Missing end year.")
    if (length(nc_years) != 0) 
      end.year <- max(nc_years)
  }
  if (!is.na(start.year) && !is.na(end.year)) {
    keep <- which(nc_years >= start.year & nc_years <= end.year)
    ncfiles <- ncfiles[keep]
  }
  else {
    PEcAn.logger::logger.info("No start or end year provided; reading output for all years")
  }
  origin_year <- start.year
  if (!is.finite(origin_year)) {
    PEcAn.logger::logger.warn("Invalid (or missing) origin year `", 
                              origin_year, "`. ", "Setting origin year to 1970.")
    origin_year <- 1970
  }
  run_origin <- paste0(origin_year, "-01-01")
  nofiles <- FALSE
  if (length(ncfiles) == 0) {
    PEcAn.logger::logger.warn("read.output: no netCDF files of model output present", 
                              "for runid = ", runid, " in ", outdir, " for years ", 
                              start.year, ":", end.year, ".", "Returning NA.")
    if (length(nc_years) > 0) {
      PEcAn.logger::logger.info("netCDF files for other years present: ", 
                                nc_years)
    }
    nofiles <- TRUE
  }
  else {
    PEcAn.logger::logger.info("Reading the following files: ", 
                              normalizePath(ncfiles))
  }
  result <- list()
  if (nofiles) {
    PEcAn.logger::logger.info("No files found. Returning all NA.")
    result <- lapply(variables, function(x) NA)
  }
  else {
    for (ncfile in ncfiles) {
      if (verbose) 
        PEcAn.logger::logger.debug("Processing file: ", 
                                   ncfile)
      nc <- ncdf4::nc_open(ncfile, verbose = verbose)
      if (is.null(variables)) {
        variables <- names(nc[["var"]])
        PEcAn.logger::logger.info("Variables not specified. Reading output for all variables, which are as follows: ", 
                                  paste(variables, collapse = ", "))
      }
      if (dataframe) {
        seconds <- ud_convert(nc$dim$time$vals, nc$dim$time$units, 
                              paste("seconds since", run_origin))
        result[["posix"]] <- abind::abind(result[["posix"]], 
                                          seconds)
      }
      for (v in variables) {
        if (verbose) 
          PEcAn.logger::logger.debug("Processing variable: ", 
                                     v)
        if (!(v %in% c(names(nc$var), names(nc$dim)))) {
          PEcAn.logger::logger.warn(paste(v, "missing in", 
                                          ncfile))
          next
        }
        newresult <- ncdf4::ncvar_get(nc, v, verbose = verbose)
        if ("pft" %in% sapply(nc$var[[v]]$dim, `[[`, 
                              "name")) {
          pft.string <- ncdf4::ncatt_get(nc, "PFT", verbose = verbose)
          pft.ind <- strsplit(pft.string$long_name, ",")[[1]] == 
            pft.name
          dim.check <- length(dim(newresult))
          if (any(pft.ind)) {
            if (dim.check == 1) {
              newresult <- newresult[pft.ind]
            }
            else {
              newresult <- newresult[, pft.ind]
            }
          }
          else {
            if (dim.check == 1) {
              newresult <- sum(newresult)
            }
            else {
              newresult <- apply(newresult, 1, sum)
            }
          }
        }
        result[[v]] <- abind::abind(result[[v]], newresult)
      }
      ncdf4::nc_close(nc)
    }
    if (print_summary) {
      result_means <- vapply(result, mean, numeric(1), 
                             na.rm = TRUE)
      result_medians <- vapply(result, stats::median, numeric(1), 
                               na.rm = TRUE)
      summary_matrix <- signif(cbind(Mean = result_means, 
                                     Median = result_medians), 3)
      rownames(summary_matrix) <- names(result)
      PEcAn.logger::logger.info("Result summary:\n", PEcAn.logger::print2string(summary_matrix), 
                                wrap = FALSE)
    }
  }
  if (!dataframe) 
    return(result)
  for (var in names(result)) {
    c <- dim(result[[var]])[2]
    r <- dim(result[[var]])[1]
    if (!is.na(c) & r > 1) {
      PEcAn.logger::logger.warn("Variable", var, "has", 
                                r, "dimensions,\n      it cannot be loaded and will be omitted.")
      result[[var]] <- NULL
    }
  }
  model <- as.data.frame(result)
  model[["posix"]] <- as.POSIXct(model[["posix"]], origin = run_origin, 
                                 tz = "UTC")
  model[["year"]] <- lubridate::year(model[["posix"]])
  return(model)
  
  # print("finish read.output")
}

##' @title analysis_sda_block (debug version)
##' @name  analysis_sda_block
##' @author 修改：谷洋 + ChatGPT
##' 
##' @description
##' 串行版 + 超啰嗦 debug 打印：
##'  - 打印每个 block 的关键结构信息（y.censored, H, aqq/bqq 等）
##'  - 检查明显的维度不一致
##'  - 如果 MCMC_block_function 出错，保存该 block 到 RDS 文件
##' 
##' @return 和原版一样：list(block.list.all, mu.f, Pf, mu.a, Pa, Y, R, analysis)
# analysis_sda_block <- function(settings, block.list.all, X, obs.mean, obs.cov,
#                                t, nt, MCMC.args, block.list.all.pre = NULL) {
#   
#   cores <- 1  # force serial for debug
#   
#   ## Step 1. Build block structure
#   block.results <- build.block.xy(
#     settings = settings,
#     block.list.all = block.list.all,
#     X = X,
#     obs.mean = obs.mean,
#     obs.cov = obs.cov,
#     t = t
#   )
#   
#   block.list.all <- block.results[[1]]
#   H  <- block.results[[2]]
#   Y  <- block.results[[3]]
#   R  <- block.results[[4]]
#   
#   ## Step 2. update q
#   block.list.all <- update_q(
#     block.list.all, t, nt,
#     aqq.Init = as.numeric(settings$state.data.assimilation$aqq.Init),
#     bqq.Init = as.numeric(settings$state.data.assimilation$bqq.Init),
#     MCMC_dat = NULL,
#     block.list.all.pre
#   )
#   
#   ## Step 3. add init
#   block.list.all[[t]] <- MCMC_Init(block.list.all[[t]], X)
#   
#   ## Step 4. assign MCMC args
#   block.list.all[[t]] <- purrr::map(block.list.all[[t]], function(l){
#     l$MCMC <- MCMC.args
#     return(l)
#   })
#   
#   PEcAn.logger::logger.info(
#     sprintf("Running MCMC (serial, debug) for %s blocks at t = %s",
#             length(block.list.all[[t]]), t)
#   )
#   
#   ## Step 5. Loop over blocks
#   for (b in seq_along(block.list.all[[t]])) {
#     block <- block.list.all[[t]][[b]]
#     
#     .print_block_info <- function(block, t, b) {
#       PEcAn.logger::logger.info(sprintf(
#         "Block %s at t = %s: len(y.censored) = %s YN = %s dim(H) = %s x %s constant$H = %s length(constant$H) = %s dim(aqq) = %s x %s dim(bqq) = %s x %s",
#         b, t, length(block$data$y.censored), block$constant$YN,
#         nrow(block$H), ncol(block$H),
#         paste(block$constant$H, collapse = " "),
#         length(block$constant$H),
#         nrow(block$aqq), ncol(block$aqq),
#         nrow(block$bqq), ncol(block$bqq)
#       ))
#     }
#     
#     PEcAn.logger::logger.info(sprintf("-> [DEBUG] Checking block %s / %s",
#                                       b, length(block.list.all[[t]])))
#     .print_block_info(block, t, b)
#     
#     ## Step 6. run MCMC for this block
#     tryres <- try( MCMC_block_function(block), silent = TRUE )
#     
#     if (inherits(tryres, "try-error")) {
#       PEcAn.logger::logger.severe(
#         paste0("MCMC_block_function FAILED for block ", b,
#                " at t = ", t, " with error: ", as.character(tryres))
#       )
#       stop(as.character(tryres))
#     }
#     
#     block.list.all[[t]][[b]] <- tryres
#   }
#   
#   # Step 7. assemble vectors
#   # V <- block.2.vector(block.list.all[[t]], X, H)
#   # 
#   # return(list(
#   #   block.list.all = block.list.all,
#   #   mu.f = V$mu.f,
#   #   Pf = V$Pf,
#   #   mu.a = V$mu.a,
#   #   Pa = V$Pa,
#   #   Y = Y,
#   #   R = R,
#   #   analysis = V$analysis
#   # ))
#   
#   ## Step 7. assemble vectors
#   V <- block.2.vector(block.list.all[[t]], X, H)
#   
#   ## Step 8. collect diagnostics
#   diag_table <- purrr::map_dfr(seq_along(block.list.all[[t]]), function(i) {
#     L <- block.list.all[[t]][[i]]
#     
#     data.frame(
#       block_id = i,
#       t_index = t,
#       site_ids = paste(L$site.ids, collapse = ","),
#       YN = L$constant$YN,
#       q_type = L$constant$q.type,
#       rhat_max = L$diag$summary$rhat_max,
#       rhat_median = L$diag$summary$rhat_median,
#       ess_min = L$diag$summary$ess_min,
#       ess_median = L$diag$summary$ess_median,
#       nchain = L$diag$summary$nchain,
#       niter = L$diag$summary$niter,
#       nburnin = L$diag$summary$nburnin,
#       nthin = L$diag$summary$nthin,
#       stringsAsFactors = FALSE
#     )
#   })
#   
#   return(list(
#     block.list.all = block.list.all,
#     mu.f = V$mu.f,
#     Pf = V$Pf,
#     mu.a = V$mu.a,
#     Pa = V$Pa,
#     Y = Y,
#     R = R,
#     analysis = V$analysis,
#     diag_table = diag_table
#   ))
# }



# ##' @title analysis_sda_block
# ##' @name  analysis_sda_block
# ##' @author Dongchen Zhang
# ##' 
# ##' @param settings  pecan standard multi-site settings list.  
# ##' @param block.list.all Lists of forecast and analysis outputs for each time point of each block. If t=1, we initialize those outputs of each block with NULL from the `sda.enkf.multisite` function.
# ##' @param X A matrix contains ensemble forecasts with the dimensions of `[ensemble number, site number * number of state variables]`. The columns are matched with the site.ids and state variable names of the inside the `FORECAST` object in the `sda.enkf.multisite` script. 
# ##' @param obs.mean Lists of date times named by time points, which contains lists of sites named by site ids, which contains observation means for each state variables of each site for each time point. 
# ##' @param obs.cov   Lists of date times named by time points, which contains lists of sites named by site ids, which contains observation covariances for all state variables of each site for each time point. 
# ##' @param t time point in format of YYYY-MM-DD.
# ##' @param nt total length of time steps, corresponding to the `nt` variable in the `sda.enkf.multisite` function.
# ##' @param MCMC.args arguments for the MCMC sampling, details can be found in the roxygen strucutre for control list in the `sda.enkf.multisite` function.
# ##' @param block.list.all.pre pre-existed block.list.all object for passing the aqq and bqq to the current SDA run, the default is NULL. Details can be found in the roxygen structure for `pre_enkf_params` of the `sda.enkf.multisite` function.
# ##' @details This function will add data and constants into each block that are needed for the MCMC sampling.
# ##'  
# ##' @description This function provides the block-based MCMC sampling approach.
# ##' 
# ##' @return It returns the `build.block.xy` object and the analysis results.
# ##' @importFrom dplyr %>%
analysis_sda_block <- function (settings, block.list.all, X, obs.mean, obs.cov,
                                t, nt, MCMC.args, block.list.all.pre = NULL) {
  # grab cores from settings.
  cores <- as.numeric(settings$state.data.assimilation$batch.settings$general.job$cores)
  # if we didn't assign number of CPUs in the settings.
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
    # if we only have one CPU.
    if (cores < 1) cores <- 1
  }
  
  # convert from vector values to block lists.
  if ("try-error" %in% class(try(
    block.results <- build.block.xy(
      settings      = settings,
      block.list.all = block.list.all,
      X             = X,
      obs.mean      = obs.mean,
      obs.cov       = obs.cov,
      t             = t
    )
  ))) {
    PEcAn.logger::logger.severe("Something wrong within the build.block.xy function.")
    return(0)
  }
  
  # grab block.list and H from the results.
  block.list.all <- block.results[[1]]
  H              <- block.results[[2]]
  Y              <- block.results[[3]]
  R              <- block.results[[4]]
  
  # update q.
  if ("try-error" %in% class(try(
    block.list.all <- update_q(
      block.list.all   = block.list.all,
      t                = t,
      nt               = nt,
      aqq.Init         = as.numeric(settings$state.data.assimilation$aqq.Init),
      bqq.Init         = as.numeric(settings$state.data.assimilation$bqq.Init),
      MCMC_dat         = NULL,
      block.list.all.pre = block.list.all.pre
    )
  ))) {
    PEcAn.logger::logger.severe("Something wrong within the update_q function.")
    return(0)
  }
  
  # add initial conditions for the MCMC sampling.
  if ("try-error" %in% class(try(
    block.list.all[[t]] <- MCMC_Init(block.list.all[[t]], X)
  ))) {
    PEcAn.logger::logger.severe("Something wrong within the MCMC_Init function.")
    return(0)
  }
  
  # update MCMC args.
  block.list.all[[t]] <- block.list.all[[t]] %>%
    purrr::map(function(l) {
      l$MCMC <- MCMC.args
      l
    })
  
  # parallel for loop over each block.
  PEcAn.logger::logger.info(
    paste0("Running MCMC ", "for ", length(block.list.all[[t]]), " blocks")
  )
  
  cl <- parallel::makeCluster(as.numeric(cores))
  doSNOW::registerDoSNOW(cl)
  l <- NULL
  
  if ("try-error" %in% class(try(
    block.list.all[[t]] <- foreach::foreach(
      l         = block.list.all[[t]],
      .packages = c("Kendall", "purrr", "nimble", "PEcAnAssimSequential"),
      # ★★ 关键修改：显式导出 MCMC_block_function 到每个 worker ★★
      .export   = c("MCMC_block_function")
    ) %dopar% {
      MCMC_block_function(l)
    }
  ))) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    PEcAn.logger::logger.severe("Something wrong within the MCMC_block_function function.")
    return(0)
  }
  
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  PEcAn.logger::logger.info("Completed!")
  
  # convert from block lists to vector values.
  if ("try-error" %in% class(try(
    V <- block.2.vector(block.list.all[[t]], X, H)
  ))) {
    PEcAn.logger::logger.severe("Something wrong within the block.2.vector function.")
    return(0)
  }
  
  # add diagnose table for each block. 
  diag_table <- purrr::map_dfr(seq_along(block.list.all[[t]]), function(i) {
    L <- block.list.all[[t]][[i]]
    
    data.frame(
      block_id    = i,
      t_index     = t,
      site_ids    = paste(L$site.ids, collapse = ","),
      YN          = L$constant$YN,
      q_type      = L$constant$q.type,
      rhat_max    = L$diag$summary$rhat_max,
      rhat_median = L$diag$summary$rhat_median,
      ess_min     = L$diag$summary$ess_min,
      ess_median  = L$diag$summary$ess_median,
      nchain      = L$diag$summary$nchain,
      niter       = L$diag$summary$niter,
      nburnin     = L$diag$summary$nburnin,
      nthin       = L$diag$summary$nthin,
      stringsAsFactors = FALSE
    )
  })
  
  # return values
  return(list(
    block.list.all = block.list.all,
    mu.f           = V$mu.f,
    Pf             = V$Pf,
    mu.a           = V$mu.a,
    Pa             = V$Pa,
    Y              = Y,
    R              = R,
    analysis       = V$analysis,
    diag_table     = diag_table
  ))
}

##' @title build.block.xy
##' @name  build.block.xy
##' @author Dongchen Zhang
##' 
##' @param settings  pecan standard multi-site settings list.  
##' @param block.list.all List contains nt empty sub-elements.
##' @param X A matrix contains ensemble forecasts.
##' @param obs.mean List of dataframe of observation means, named with observation datetime.
##' @param obs.cov   List of covariance matrices of state variables , named with observation datetime.
##' @param t time point.
##' @details This function will add data and constants into each block that are needed for the MCMC sampling.
##'  
##' @description This function split long vector and covariance matrix into blocks corresponding to the localization.
##' 
build.block.xy <- function(settings, block.list.all, X, obs.mean, obs.cov, t) {
  #set q.type from settings.
  if (settings$state.data.assimilation$q.type == "vector") {
    q.type <- 3
  } else if (settings$state.data.assimilation$q.type == "wishart") {
    q.type <- 4
  }
  #grab basic arguments based on X.
  site.ids <- unique(attributes(X)$Site)
  var.names <- unique(colnames(X))
  mu.f <- colMeans(X)
  Pf <- stats::cov(X)
  if (length(diag(Pf)[which(diag(Pf)==0)]) > 0) {
    diag(Pf)[which(diag(Pf)==0)] <- min(diag(Pf)[which(diag(Pf) != 0)])/5 #fixing det(Pf)==0
    PEcAn.logger::logger.warn("The zero variances in Pf is being replaced by one fifth of the minimum variance in those matrices respectively.")
  }
  #distance calculations and localization
  if (!is.null(settings$state.data.assimilation$Localization.FUN) && 
      ! as.numeric(settings$state.data.assimilation$scalef) == 0) {
    site.locs <- settings$run %>%
      purrr::map('site') %>%
      purrr::map_dfr(~c(.x[['lon']],.x[['lat']]) %>% as.numeric)%>%
      t %>%
      `colnames<-`(c("Lon","Lat")) %>%
      `rownames<-`(site.ids)
    #Finding the distance between the sites
    dis.matrix <- sp::spDists(site.locs, longlat = TRUE)
    Localization.FUN <- get(settings$state.data.assimilation$Localization.FUN)
    #turn that into a blocked matrix format
    blocked.dis <- block_matrix(dis.matrix %>% as.numeric(), rep(length(var.names), length(site.ids)))
    Pf <- Localization.FUN(Pf, blocked.dis, settings$state.data.assimilation$scalef %>% as.numeric())
  }
  #Handle observation
  #observation number per site
  #free run special case.
  if (is.null(obs.mean[[t]])) {
    obs_per_site <- rep(0, length(site.ids)) %>% purrr::set_names(site.ids)
  } else {
    obs_per_site <- purrr::map_int(obs.mean[[t]], length)
  }
  #if we do free run or the current obs.mean are all NULL.
  if (as.logical(settings$state.data.assimilation$free.run) | all(is.null(unlist(obs.mean[[t]])))) {
    H <- list(ind = seq_along(rep(var.names, length(site.ids))))
    Y <- rep(NA, length(H$ind))
    R <- diag(1, length(H$ind))
  } else if (!as.logical(settings$state.data.assimilation$free.run) && all(is.null(unlist(obs.mean[[t]])))) {
    PEcAn.logger::logger.error("Please set the settings$state.data.assimilation$free.run as TRUE if you don't have any observations!")
    return(0)
  } else {
    Obs.cons <- Construct.R(site.ids, var.names, obs.mean[[t]], obs.cov[[t]])
    Y <- Obs.cons$Y
    R <- Obs.cons$R
    if (length(Y) > 1) {
      if (length(diag(R)[which(diag(R)==0)]) > 0) {
        diag(R)[which(diag(R)==0)] <- min(diag(R)[which(diag(R) != 0)])/2
        PEcAn.logger::logger.warn("The zero variances in R is being replaced by half of the minimum variance in those matrices respectively.")
      }
    }
    #create matrix the describes the support for each observed state variable at time t
    min_max <- settings$state.data.assimilation$state.variables %>% 
      purrr::map(function(state.variable){
        c(as.numeric(state.variable$min_value),
          as.numeric(state.variable$max_value))
      }) %>% unlist() %>% as.vector() %>% 
      matrix(length(settings$state.data.assimilation$state.variables), 2, byrow = T) %>%
      `rownames<-`(var.names)
    #Create y.censored and y.ind
    #describing if the obs are within the defined range.
    y.ind <- y.censored <- c()
    for (i in seq_along(Y)) {
      if (Y[i] > min_max[names(Y[i]), 1]) {
        y.ind[i] = 1; y.censored[i] = Y[i]
      } else {y.ind[i] <- y.censored[i] <- 0}
    }
    #create H
    # if there is any site that has zero observation.
    if (any(obs_per_site == 0)) {
      #name matching between observation names and state variable names.
      f.2.y.ind <- obs.mean[[t]] %>%
        purrr::map(\(x)which(var.names %in% names(x))) %>%
        unlist %>%
        unique %>% 
        sort
      H <- list(ind = f.2.y.ind %>% purrr::map(function(start){
        seq(start, length(site.ids) * length(var.names), length(var.names))
      }) %>% unlist() %>% sort)
    } else {
      H <- construct_nimble_H(site.ids = site.ids,
                              var.names = var.names,
                              obs.t = obs.mean[[t]],
                              pft.path = settings[[1]]$run$inputs$pft.site$path,
                              by = "block_pft_var")
    }
  }
  #start the blocking process
  #should we consider interactions between sites?
  if(as.numeric(settings$state.data.assimilation$scalef) == 0){
    block.list <- vector("list", length(site.ids))
    #loop over sites
    for (i in seq_along(site.ids)) {
      #store which block contains which sites.
      block.list[[i]]$sites.per.block <- i
      block.list[[i]]$site.ids <- site.ids[i]
      block.list[[i]]$t <- t
      #fill in mu.f and Pf
      f.start <- (i - 1) * length(var.names) + 1
      f.end <- i * length(var.names)
      block.list[[i]]$data$muf <- mu.f[f.start:f.end]
      block.list[[i]]$data$pf <- Pf[f.start:f.end, f.start:f.end]
      #find indexs for Y.
      y.start <- sum(obs_per_site[1:i]) - obs_per_site[i] + 1
      y.end <- sum(obs_per_site[1:i])
      #fill in y and r
      #if there is no observation for this site.
      if (y.end < y.start) {
        #if every site has zero observation/free run.
        if (max(obs_per_site) == 0) {
          block.list[[i]]$data$y.censored <- rep(NA, length(var.names))
          block.list[[i]]$data$r <- diag(1, length(var.names))
          block.h <- matrix(1, 1, length(var.names))
        } else {
          block.list[[i]]$data$y.censored <- rep(NA, length(f.2.y.ind))
          block.list[[i]]$data$r <- diag(1, length(f.2.y.ind))
          block.h <- matrix(NA, 1, length(var.names))
          block.h[1, f.2.y.ind] <- 1
        }
      } else {
        block.list[[i]]$data$y.censored <- y.censored[y.start:y.end]
        block.list[[i]]$data$r <- solve(R[y.start:y.end, y.start:y.end])
        block.h <- Construct.H.multisite(site.ids[i], var.names, obs.mean[[t]])
      }
      #fill in constants.
      block.list[[i]]$H <- block.h
      block.list[[i]]$constant$H <- which(apply(block.h, 2, sum) == 1)
      block.list[[i]]$constant$N <- length(f.start:f.end)
      block.list[[i]]$constant$YN <- length(block.list[[i]]$data$y.censored)
      block.list[[i]]$constant$q.type <- q.type
    }
    names(block.list) <- site.ids
  } else {
    #find networks given TRUE/FALSE matrix representing sites' interactions.
    block.vec <- matrix_network(dis.matrix <= as.numeric(settings$state.data.assimilation$scalef))
    #check if the matrix_network function is working correctly.
    #check if the blocks are calculated correctly.
    if (block.vec %>% 
        purrr::map(function(l){length(l)}) %>%
        unlist %>%
        sum() != length(site.ids)) {
      PEcAn.logger::logger.severe("Block calculation failed, please check the matrix_network function!")
      return(0)
    }
    block.list <- vector("list", length(block.vec))
    #loop over sites
    for (i in seq_along(block.vec)) {#i is site index
      #store which block contains which sites.
      ids <- block.vec[[i]]
      block.list[[i]]$sites.per.block <- ids
      block.list[[i]]$site.ids <- site.ids[ids]
      block.list[[i]]$t <- t
      y.ind <- f.ind <- na.ind <- c()
      r.block <- y.block <- c()
      for (j in seq_along(ids)) {
        f.start <- (ids[j] - 1) * length(var.names) + 1
        f.end <- ids[j] * length(var.names)
        y.start <- sum(obs_per_site[1:ids[j]]) - obs_per_site[ids[j]] + 1
        y.end <- sum(obs_per_site[1:ids[j]])
        f.ind <- c(f.ind, f.start:f.end)
        #if the current site has greater or equal than 1 observation.
        if (y.end >= y.start) {
          # y.ind <- c(y.ind, y.start:y.end)
          y.block <- c(y.block, y.censored[y.start:y.end])
          r.block <- c(r.block, diag(R)[y.start:y.end])
        } else {
          #if the current site has zero observation.
          #if for free run.
          if (max(obs_per_site) == 0) {
            y.block <- c(y.block, rep(NA, length(var.names)))
            r.block <- c(r.block, rep(1, length(var.names)))
          } else {
            y.block <- c(y.block, rep(NA, max(obs_per_site)))
            r.block <- c(r.block, rep(1, max(obs_per_site)))
          }
        }
      }
      #if we have NA for y, we will build H differently.
      if (any(is.na(y.block))) {
        block.h <- matrix(0, 1, length(ids)*length(var.names))
        #if for free run.
        if (is.null(obs.mean[[t]])) {
          f.2.y.ind <- seq_along(var.names)
        } else {
          f.2.y.ind <- obs.mean[[t]] %>%
            purrr::map(\(x)which(var.names %in% names(x))) %>%
            unlist %>%
            unique
        }
        seq.ind <- f.2.y.ind %>% purrr::map(function(start){
          seq(start, dim(block.h)[2], length(var.names))
        }) %>% unlist()
        block.h[1, seq.ind] <- 1
      } else {
        block.h <- Construct.H.multisite(site.ids[ids], var.names, obs.mean[[t]])
      }
      #fill in  mu.f and Pf
      block.list[[i]]$data$muf <- mu.f[f.ind]
      block.list[[i]]$data$pf <- GrabFillMatrix(Pf, f.ind)
      #fill in y and R
      block.list[[i]]$data$y.censored <- y.block
      if (length(r.block)  == 1) {
        block.list[[i]]$data$r <- 1/r.block
      } else {
        block.list[[i]]$data$r <- solve(diag(r.block))
      }
      block.list[[i]]$H <- block.h
      block.list[[i]]$constant$H <- which(apply(block.h, 2, sum) == 1)
      block.list[[i]]$constant$N <- length(f.ind)
      block.list[[i]]$constant$YN <- length(y.block)
      block.list[[i]]$constant$q.type <- q.type
    }
  }
  #if it's Wishart Q, we need to replace any NA Y with corresponding muf, and r with Pf.
  #also, if length of observation is 1, the Wishart Q is not suitable for the MCMC.
  #we will then need to change the Q type to 3, which is the vector Q.
  #the wishart-MCMC is still under development, so I commented them out for now.
  if (q.type == 4) {
    for (i in seq_along(block.list)) {
      #check length.
      if (block.list[[i]]$constant$YN == 1) {
        block.list[[i]]$constant$q.type <- 3
        next
      }
      # #check NAs.
      # na.ind <- which(is.na(block.list[[i]]$data$y.censored))
      # if (length(na.ind) > 0) {
      #     block.list[[i]]$constant$YN <- block.list[[i]]$constant$YN - length(na.ind)
      #     block.list[[i]]$constant$H <- block.list[[i]]$constant$H[-na.ind]
      #     block.list[[i]]$data$y.censored <- block.list[[i]]$data$y.censored[-na.ind]
      #     block.list[[i]]$data$r <- diag(diag(block.list[[i]]$data$r)[-na.ind])
      # }
      # na.site.ind <- which(obs_per_site[block.list[[i]]$site.ids] == 0)
      # na.ind <- which(is.na(block.list[[i]]$data$y.censored))
      # if (length(na.site.ind) > 0) {
      #   site.inds <- block.list[[i]]$sites.per.block[na.site.ind]
      #   y.2.muf.ind <- f.2.y.ind %>% purrr::map(function(start){
      #     seq(start, length(mu.f), length(var.names))[site.inds]
      #   }) %>% unlist() %>% sort()
      #   block.list[[i]]$data$y.censored[na.ind] <- mu.f[y.2.muf.ind]
      #   block.list[[i]]$data$r[na.ind, na.ind] <- Pf[y.2.muf.ind, y.2.muf.ind]
      # }
    }
  }
  #return values.
  block.list.all[[t]] <- block.list
  return(list(block.list.all = block.list.all, H = H, Y = Y, R = R))
}

##' @title MCMC_Init
##' @name  MCMC_Init
##' @author Dongchen Zhang
##' 
##' @param block.list  lists of blocks generated by the `build.block.xy` function.
##' @param X A matrix contains ensemble forecasts.
##' @details This function helps create initial conditions for the MCMC sampling.
##' 
##' @return It returns the `block.list` object with initial conditions filled in.
MCMC_Init <- function (block.list, X) {
  var.names <- unique(names(X))
  #sample mu.f from X.
  sample.mu.f <- colMeans(X)
  for (i in seq_along(block.list)) {
    #number of observations.
    num.obs <- length(block.list[[i]]$data$y.censored)
    #loop over each site within each block
    for (j in seq_along(block.list[[i]]$sites.per.block)) {
      #initialize mu.f
      start <- (block.list[[i]]$sites.per.block[j] - 1) * length(var.names) + 1
      end <- (block.list[[i]]$sites.per.block[j]) * length(var.names)
      block.list[[i]]$Inits$X.mod <- c(block.list[[i]]$Inits$X.mod, sample.mu.f[start:end])
      #initialize X
      block.list[[i]]$Inits$X <- block.list[[i]]$data$y.censored
      #initialize Xs
      block.list[[i]]$Inits$Xs <- block.list[[i]]$Inits$X.mod[block.list[[i]]$constant$H]
    }
    #initialize q.
    #if we want the vector q.
    if (block.list[[i]]$constant$q.type == 3) {
      for (j in seq_along(block.list[[i]]$data$y.censored)) {
        temp.q <- stats::rgamma(1, shape = block.list[[i]]$data$aq[j], rate = block.list[[i]]$data$bq[j])
        if (temp.q < 0.001) {
          temp.q <- 0.001
        }
        block.list[[i]]$Inits$q <- c(block.list[[i]]$Inits$q, temp.q)
      }
    } else if (block.list[[i]]$constant$q.type == 4) {
      #if we want the wishart Q.
      if ("try-error" %in% class(try(block.list[[i]]$Inits$q <- 
                                     stats::rWishart(1, df = block.list[[i]]$data$bq, Sigma = block.list[[i]]$data$aq)[,,1], silent = T))) {
        block.list[[i]]$Inits$q <- 
          stats::rWishart(1, df = block.list[[i]]$data$bq, Sigma = stats::toeplitz((block.list[[i]]$constant$YN:1)/block.list[[i]]$constant$YN))[,,1]
      }
    }
  }
  #return values.
  return(block.list)
}

##' @title MCMC_block_function
##' @name  MCMC_block_function
##' @author Dongchen Zhang
##' 
##' @param block  each block within the `block.list` lists.
##' 
##' @return It returns the `block` object with analysis results filled in.
# MCMC_block_function <- function(block) {
#   
#   PEcAn.logger::logger.info(
#     sprintf("MCMC_block_function: t = %s, len(y.censored) = %s, YN = %s, q.type = %s",
#             block$t, length(block$data$y.censored), block$constant$YN, block$constant$q.type)
#   )
#   
#   nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE,
#                 checkNimbleFunction = FALSE, checkDuplicateNodeDefinitions = FALSE)
#   
#   ## ---- 1. Build nimble model ----
#   model_pred <- nimble::nimbleModel(
#     GEF.MultiSite.Nimble,
#     data = block$data,
#     inits = block$Inits,
#     constants = block$constant,
#     name = 'base'
#   )
#   
#   ## ---- 2. Configure MCMC ----
#   conf <- nimble::configureMCMC(model_pred, print = FALSE)
#   conf$setMonitors(c("X", "X.mod", "q"))
#   
#   samplerLists <- conf$getSamplers()
#   samplerNumberOffset <- length(samplerLists)
#   
#   ## Replace all RW samplers with ESS sampler if Wishart
#   if (block$constant$q.type == 4) {
#     samplerLists <- purrr::map(samplerLists, function(l) { l$setName("ess"); l })
#   }
#   
#   conf$setSamplers(samplerLists)
#   
#   ## ---- 3. Replace sampler for X.mod with ESS ----
#   xmod_targets <- purrr::map_chr(samplerLists, ~ .x$target)
#   xmod_ind <- grep("X.mod", xmod_targets)
#   
#   if (length(xmod_ind) == 1) {
#     conf$removeSampler(samplerLists[[xmod_ind]]$target)
#     conf$addSampler(
#       target = samplerLists[[xmod_ind]]$target,
#       type = "ess",
#       control = list(propCov = block$data$pf,
#                      adaptScaleOnly = TRUE,
#                      latents = "X",
#                      pfOptimizeNparticles = TRUE)
#     )
#   }
#   
#   PEcAn.logger::logger.info(
#     sprintf("MCMC_block_function: total samplers = %s; targets = [%s]",
#             length(samplerLists), paste(xmod_targets, collapse = ", "))
#   )
#   
#   ## ---- 4. Add toggle sampler for censored Y ----
#   for (i in 1:block$constant$YN) {
#     conf$addSampler(paste0("y.censored[", i, "]"), "toggle", control = list(type = "RW"))
#   }
#   
#   ## ---- 5. Compile ----
#   Rmcmc <- nimble::buildMCMC(conf)
#   Cmodel <- nimble::compileNimble(model_pred)
#   Cmcmc <- nimble::compileNimble(Rmcmc, project = model_pred, showCompilerOutput = FALSE)
#   
#   ## ---- 6. Disable toggle sampler if no NA ----
#   samplerFns <- Cmcmc$samplerFunctions
#   nsam <- length(samplerFns)
#   
#   PEcAn.logger::logger.info(sprintf(
#     "MCMC_block_function: samplerNumberOffset = %s, YN = %s, length(samplerFunctions) = %s",
#     samplerNumberOffset, block$constant$YN, nsam))
#   
#   if (nsam >= samplerNumberOffset + block$constant$YN) {
#     # safe to toggle
#     if (!any(is.na(block$data$y.censored))) {
#       for (i in 1:block$constant$YN) {
#         idx <- samplerNumberOffset + i
#         valueInCompiledNimbleFunction(samplerFns[[idx]], "toggle", 0)
#       }
#     }
#   } else {
#     PEcAn.logger::logger.warn(
#       sprintf("MCMC_block_function: samplerFunctions length (%s) < samplerNumberOffset + YN (%s); skip toggling.",
#               nsam, samplerNumberOffset + block$constant$YN))
#   }
#   
#   ## ---- 7. Run MCMC ----
#   dat <- runMCMC(Cmcmc,
#                  niter = block$MCMC$niter,
#                  nburnin = block$MCMC$nburnin,
#                  thin = block$MCMC$nthin,
#                  nchains = block$MCMC$nchain)
#   
#   PEcAn.logger::logger.info(
#     sprintf("MCMC_block_function: MCMC finished; dim(dat) = %s x %s",
#             nrow(dat), ncol(dat))
#   )
#   
#   ## ---- 8. Extract posterior for X and X.mod ----
#   iX      <- grep("^X\\[", colnames(dat))
#   iX_mod  <- grep("^X.mod\\[", colnames(dat))
#   
#   PEcAn.logger::logger.info(
#     sprintf("MCMC_block_function: iX length = %s, iX.mod length = %s",
#             length(iX), length(iX_mod))
#   )
#   
#   # Ensure non-empty
#   if (length(iX) == 0)      stop("ERROR: iX is empty — no X[] in MCMC output.")
#   if (length(iX_mod) == 0)  stop("ERROR: iX.mod is empty — no X.mod[] in MCMC output.")
#   
#   ## ---- 9. Compute mua, pa ----
#   if (length(iX) == 1) {
#     mua <- mean(dat[, iX])
#     pa  <- var(dat[, iX])
#   } else {
#     mua <- colMeans(dat[, iX])
#     pa  <- cov(dat[, iX])
#   }
#   
#   ## ---- 10. Compute mufa, pfa ----
#   mufa <- colMeans(dat[, iX_mod])
#   pfa  <- cov(dat[, iX_mod])
#   
#   ## ---- 11. Return ----
#   block$update <- list(
#     mua = mua,
#     pa  = pa,
#     mufa = mufa,
#     pfa = pfa,
#     aq = NA,
#     bq = NA
#   )
#   
#   return(block)
# }
#################MCMC_block_function with diagnose!!!!!
MCMC_block_function <- function(block) {
  
  PEcAn.logger::logger.info(
    sprintf("MCMC_block_function: t = %s, len(y.censored) = %s, YN = %s, q.type = %s",
            block$t, length(block$data$y.censored), block$constant$YN, block$constant$q.type)
  )
  
  nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE,
                checkNimbleFunction = FALSE, checkDuplicateNodeDefinitions = FALSE)
  
  ## ---- 1. Build nimble model ----
  model_pred <- nimble::nimbleModel(
    GEF.MultiSite.Nimble,
    data = block$data,
    inits = block$Inits,
    constants = block$constant,
    name = 'base'
  )
  
  ## ---- 2. Configure MCMC ----
  conf <- nimble::configureMCMC(model_pred, print = FALSE)
  conf$setMonitors(c("X", "X.mod", "q"))
  
  samplerLists <- conf$getSamplers()
  samplerNumberOffset <- length(samplerLists)
  
  ## Replace all RW samplers with ESS sampler if Wishart
  if (block$constant$q.type == 4) {
    samplerLists <- purrr::map(samplerLists, function(l) { l$setName("ess"); l })
  }
  
  conf$setSamplers(samplerLists)
  
  ## ---- 3. Replace sampler for X.mod with ESS ----
  xmod_targets <- purrr::map_chr(samplerLists, ~ .x$target)
  xmod_ind <- grep("X.mod", xmod_targets)
  
  if (length(xmod_ind) == 1) {
    conf$removeSampler(samplerLists[[xmod_ind]]$target)
    conf$addSampler(
      target = samplerLists[[xmod_ind]]$target,
      type = "ess",
      control = list(propCov = block$data$pf,
                     adaptScaleOnly = TRUE,
                     latents = "X",
                     pfOptimizeNparticles = TRUE)
    )
  }
  
  PEcAn.logger::logger.info(
    sprintf("MCMC_block_function: total samplers = %s; targets = [%s]",
            length(samplerLists), paste(xmod_targets, collapse = ", "))
  )
  
  ## ---- 4. Add toggle sampler for censored Y ----
  for (i in 1:block$constant$YN) {
    conf$addSampler(paste0("y.censored[", i, "]"), "toggle", control = list(type = "RW"))
  }
  
  ## ---- 5. Compile ----
  Rmcmc <- nimble::buildMCMC(conf)
  Cmodel <- nimble::compileNimble(model_pred)
  Cmcmc <- nimble::compileNimble(Rmcmc, project = model_pred, showCompilerOutput = FALSE)
  
  ## ---- 6. Disable toggle sampler if no NA ----
  samplerFns <- Cmcmc$samplerFunctions
  nsam <- length(samplerFns)
  
  PEcAn.logger::logger.info(sprintf(
    "MCMC_block_function: samplerNumberOffset = %s, YN = %s, length(samplerFunctions) = %s",
    samplerNumberOffset, block$constant$YN, nsam))
  
  if (nsam >= samplerNumberOffset + block$constant$YN) {
    if (!any(is.na(block$data$y.censored))) {
      for (i in 1:block$constant$YN) {
        idx <- samplerNumberOffset + i
        valueInCompiledNimbleFunction(samplerFns[[idx]], "toggle", 0)
      }
    }
  } else {
    PEcAn.logger::logger.warn(
      sprintf("MCMC_block_function: samplerFunctions length (%s) < samplerNumberOffset + YN (%s); skip toggling.",
              nsam, samplerNumberOffset + block$constant$YN))
  }
  
  ## ---- 7. Run MCMC (keep chains separately for diagnostics) ----
  dat <- runMCMC(
    Cmcmc,
    niter = block$MCMC$niter,
    nburnin = block$MCMC$nburnin,
    thin = block$MCMC$nthin,
    nchains = block$MCMC$nchain,
    samplesAsCodaMCMC = TRUE,
    summary = FALSE
  )
  
  PEcAn.logger::logger.info("MCMC finished.")
  
  ## ---- 8. Diagnostics ----
  rhat_mat <- tryCatch(
    coda::gelman.diag(dat, autoburnin = FALSE, multivariate = FALSE)$psrf,
    error = function(e) NULL
  )
  
  ess_vec <- tryCatch(
    coda::effectiveSize(dat),
    error = function(e) NULL
  )
  
  ## block-level 摘要
  diag_summary <- list(
    rhat_max = if (!is.null(rhat_mat)) max(rhat_mat[, 1], na.rm = TRUE) else NA_real_,
    rhat_median = if (!is.null(rhat_mat)) stats::median(rhat_mat[, 1], na.rm = TRUE) else NA_real_,
    ess_min = if (!is.null(ess_vec)) min(ess_vec, na.rm = TRUE) else NA_real_,
    ess_median = if (!is.null(ess_vec)) stats::median(ess_vec, na.rm = TRUE) else NA_real_,
    nchain = block$MCMC$nchain,
    niter = block$MCMC$niter,
    nburnin = block$MCMC$nburnin,
    nthin = block$MCMC$nthin
  )
  
  ## per-parameter diagnostics
  diag_by_param <- NULL
  if (!is.null(rhat_mat) || !is.null(ess_vec)) {
    
    par_names <- union(
      if (!is.null(rhat_mat)) rownames(rhat_mat) else character(0),
      if (!is.null(ess_vec)) names(ess_vec) else character(0)
    )
    
    diag_by_param <- data.frame(
      par_name = par_names,
      rhat = if (!is.null(rhat_mat)) rhat_mat[par_names, 1] else NA_real_,
      ess  = if (!is.null(ess_vec)) as.numeric(ess_vec[par_names]) else NA_real_,
      stringsAsFactors = FALSE
    )
    
    ## 参数映射回状态变量
    var_names <- names(block$data$muf)
    if (is.null(var_names)) {
      var_names <- paste0("state_", seq_along(block$data$muf))
    }
    
    obs_state_names <- var_names[block$constant$H]
    
    diag_by_param$state_var  <- NA_character_
    diag_by_param$param_type <- NA_character_
    
    ## X[k]
    ix <- grepl("^X\\[", diag_by_param$par_name)
    if (any(ix)) {
      k <- as.integer(sub("^X\\[(\\d+)\\]$", "\\1", diag_by_param$par_name[ix]))
      diag_by_param$state_var[ix]  <- var_names[k]
      diag_by_param$param_type[ix] <- "X"
    }
    
    ## X.mod[k]
    ixm <- grepl("^X\\.mod\\[", diag_by_param$par_name)
    if (any(ixm)) {
      k <- as.integer(sub("^X\\.mod\\[(\\d+)\\]$", "\\1", diag_by_param$par_name[ixm]))
      diag_by_param$state_var[ixm]  <- var_names[k]
      diag_by_param$param_type[ixm] <- "X.mod"
    }
    
    ## q[j]
    iq <- grepl("^q\\[", diag_by_param$par_name)
    if (any(iq)) {
      j <- as.integer(sub("^q\\[(\\d+)\\]$", "\\1", diag_by_param$par_name[iq]))
      diag_by_param$state_var[iq]  <- obs_state_names[j]
      diag_by_param$param_type[iq] <- "q"
    }
  }
  
  
  
  ## ---- 9. Merge chains into one matrix for posterior summaries ----
  dat_mat <- as.matrix(do.call(rbind, dat))
  
  PEcAn.logger::logger.info(
    sprintf("MCMC_block_function: combined sample dim(dat_mat) = %s x %s",
            nrow(dat_mat), ncol(dat_mat))
  )
  
  ## ---- 10. Extract posterior for X and X.mod ----
  iX      <- grep("^X\\[", colnames(dat_mat))
  iX_mod  <- grep("^X.mod\\[", colnames(dat_mat))
  
  PEcAn.logger::logger.info(
    sprintf("MCMC_block_function: iX length = %s, iX.mod length = %s",
            length(iX), length(iX_mod))
  )
  
  if (length(iX) == 0)     stop("ERROR: iX is empty — no X[] in MCMC output.")
  if (length(iX_mod) == 0) stop("ERROR: iX.mod is empty — no X.mod[] in MCMC output.")
  
  ## ---- 11. Compute mua, pa ----
  if (length(iX) == 1) {
    mua <- mean(dat_mat[, iX])
    pa  <- matrix(stats::var(dat_mat[, iX]), nrow = 1, ncol = 1)
  } else {
    mua <- colMeans(dat_mat[, iX, drop = FALSE])
    pa  <- stats::cov(dat_mat[, iX, drop = FALSE])
  }
  
  ## ---- 12. Compute mufa, pfa ----
  if (length(iX_mod) == 1) {
    mufa <- mean(dat_mat[, iX_mod])
    pfa  <- matrix(stats::var(dat_mat[, iX_mod]), nrow = 1, ncol = 1)
  } else {
    mufa <- colMeans(dat_mat[, iX_mod, drop = FALSE])
    pfa  <- stats::cov(dat_mat[, iX_mod, drop = FALSE])
  }
  
  ## ---- 13. Return ----
  block$update <- list(
    mua = mua,
    pa  = pa,
    mufa = mufa,
    pfa = pfa,
    aq = NA,
    bq = NA
  )
  
  ## 保存诊断对象
  block$diag <- list(
    summary = diag_summary,
    by_param = diag_by_param,
    rhat = rhat_mat,
    ess = ess_vec
  )
  
  return(block)
}

##' @title update_q
##' @name  update_q
##' @author Dongchen Zhang
##' 
##' @param block.list.all  each block within the `block.list` lists.
##' @param t time point.
##' @param nt total length of time steps.
##' @param aqq.Init the initial values of aqq, the default is NULL.
##' @param bqq.Init the initial values of bqq, the default is NULL.
##' @param MCMC_dat data frame of MCMC samples, the default it NULL.
##' @param block.list.all.pre pre-existed block.list.all object for passing the aqq and bqq to the current SDA run, the default is NULL.
##' 
##' @return It returns the `block.list.all` object with initialized/updated Q filled in.
update_q <- function (block.list.all, t, nt, aqq.Init = NULL, bqq.Init = NULL, MCMC_dat = NULL, block.list.all.pre = NULL) {
  block.list <- block.list.all[[t]]
  #if it's an update.
  if (is.null(MCMC_dat)) {
    #loop over blocks
    if (t == 1) {
      for (i in seq_along(block.list)) {
        nvar <- length(block.list[[i]]$data$muf)
        nobs <- length(block.list[[i]]$data$y.censored)
        if (block.list[[i]]$constant$q.type == 3) {
          #initialize aqq and bqq for nt
          if (!is.null(aqq.Init) && !is.null(bqq.Init)) {
            a_vec <- c(4, 2, 2, 3, 2, 4)
            b_vec <- c(40, 2, 1, 6, 1, 20)
            
            block.list[[i]]$aqq <- matrix(
              rep(a_vec, nt + 1),
              nrow = nvar,
              ncol = nt + 1,
              byrow = FALSE
            )
            
            block.list[[i]]$bqq <- matrix(
              rep(b_vec, nt + 1),
              nrow = nvar,
              ncol = nt + 1,
              byrow = FALSE
            )
          } else {
            block.list[[i]]$aqq <- array(1, dim = c(nvar, nt + 1))
            block.list[[i]]$bqq <- array(1, dim = c(nvar, nt + 1))
          }
          #update aq and bq based on aqq and bqq
          block.list[[i]]$data$aq <- block.list[[i]]$aqq[block.list[[i]]$constant$H, t]
          block.list[[i]]$data$bq <- block.list[[i]]$bqq[block.list[[i]]$constant$H, t]
        } else if (block.list[[i]]$constant$q.type == 4) {
          
          
          v_diag <- c(0.1, 1.0, 2.0, 0.5, 2.0, 0.2)
          
          if (nvar != length(v_diag)) {
            stop(paste0(
              "q.type == 4: nvar = ", nvar,
              " but v_diag length = ", length(v_diag),
              ". Check block state order."
            ))
          }
          
          block.list[[i]]$aqq <- array(NA_real_, dim = c(nvar, nvar, nt + 1))
          
          for (tt in 1:(nt + 1)) {
            block.list[[i]]$aqq[,,tt] <- diag(v_diag)
          }
          
          # df / strength
          block.list[[i]]$bqq <- rep(max(nobs, nvar + 2), nt + 1)
          
          block.list[[i]]$data$aq <- GrabFillMatrix(
            block.list[[i]]$aqq[,,t],
            block.list[[i]]$constant$H
          )
          block.list[[i]]$data$bq <- block.list[[i]]$bqq[t]
        }
      }
    } else if (t > 1) {
      if (!is.null(block.list.all.pre)) {
        block.list.pre <- block.list.all.pre[[t - 1]]
      } else {
        #if we want to update q from previous SDA runs.
        block.list.pre <- block.list.all[[t - 1]]
      }
      for (i in seq_along(block.list)) {
        nvar <- length(block.list[[i]]$data$muf)
        nobs <- length(block.list[[i]]$data$y.censored)
        if (block.list[[i]]$constant$q.type == 3) {
          #copy previous aqq and bqq to the current t
          block.list[[i]]$aqq <- block.list.pre[[i]]$aqq
          block.list[[i]]$bqq <- block.list.pre[[i]]$bqq
          #update aq and bq
          block.list[[i]]$data$aq <- block.list[[i]]$aqq[block.list[[i]]$constant$H, t]
          block.list[[i]]$data$bq <- block.list[[i]]$bqq[block.list[[i]]$constant$H, t]
        } else if (block.list[[i]]$constant$q.type == 4) {
          
          # ---------- state-specific marginal scales ----------
          v_diag <- c(0.1, 1.0, 2.0, 0.5, 2.0, 0.2)
          
          if (nvar != length(v_diag)) {
            stop(paste0(
              "q.type == 4: nvar = ", nvar,
              " but v_diag length = ", length(v_diag),
              ". Check block state order."
            ))
          }
          
          # ---------- weak process correlation structure ----------
          R <- diag(nvar)
          
          # assuming order:
          # 1 AbvGrndWood
          # 2 NEE
          # 3 Qle
          # 4 LAI
          # 5 SoilMoistFrac
          # 6 TotSoilCarb
          
          R[2,3] <- R[3,2] <- 0.30   # NEE - Qle
          R[3,5] <- R[5,3] <- 0.35   # Qle - SoilMoistFrac
          R[1,4] <- R[4,1] <- 0.20   # AbvGrndWood - LAI
          R[1,6] <- R[6,1] <- 0.10   # AbvGrndWood - TotSoilCarb
          R[4,5] <- R[5,4] <- 0.15   # LAI - SoilMoistFrac
          
          # ---------- build scale matrix V = D R D ----------
          sd_vec <- sqrt(v_diag)
          Dmat <- diag(sd_vec)
          Vmat <- Dmat %*% R %*% Dmat
          
          # numerical safety: force symmetry
          Vmat <- (Vmat + t(Vmat)) / 2
          
          # numerical safety: add tiny jitter if needed
          eig_min <- min(eigen(Vmat, symmetric = TRUE, only.values = TRUE)$values)
          if (eig_min <= 1e-8) {
            Vmat <- Vmat + diag(abs(eig_min) + 1e-6, nvar)
          }
          
          block.list[[i]]$aqq <- array(NA_real_, dim = c(nvar, nvar, nt + 1))
          for (tt in 1:(nt + 1)) {
            block.list[[i]]$aqq[,,tt] <- Vmat
          }
          
          # ---------- stronger df ----------
          df_val <- max(nobs, 2 * nvar + 5)
          block.list[[i]]$bqq <- rep(df_val, nt + 1)
          
          block.list[[i]]$data$aq <- GrabFillMatrix(
            block.list[[i]]$aqq[,,t],
            block.list[[i]]$constant$H
          )
          block.list[[i]]$data$bq <- block.list[[i]]$bqq[t]
        }
      }
    }
  } else {
    #TODO: Implement the feature that Q can be updated based on the pft types.
  }
  
  #return values.
  block.list.all[[t]] <- block.list
  return(block.list.all)
}

##' @title block.2.vector
##' @name  block.2.vector
##' @author Dongchen Zhang
##' 
##' @param block.list  lists of blocks generated by the `build.block.xy` function.
##' @param X A matrix contains ensemble forecasts.
##' @param H H index created by the `construct_nimble_H` function.
##' 
##' @return It returns a list of analysis results by MCMC sampling.
##' ####### Posterior changing from mufa and pfa to mua and pa
# block.2.vector <- function (block.list, X, H) {
#   site.ids <- attributes(X)$Site
#   mu.f <- mu.a <- c()
#   Pf <- Pa <- matrix(0, length(site.ids), length(site.ids))
#   analysis <- X
#   for (L in block.list) {
#     ind <- c()
#     for (id in L$site.ids) {
#       ind <- c(ind, which(site.ids == id))
#     }
#     #convert mu.f and pf
#     mu.a[ind] <- mu.f[ind] <- L$data$muf
#     Pa[ind, ind] <- Pf[ind, ind] <- L$data$pf
#     # MVN sample based on block.
#     ##### Change the posterior
#     sample <- as.data.frame(mvtnorm::rmvnorm(nrow(X),
#                                              L$update$mua,
#                                              L$update$pa,
#                                              method = "svd"))
#     analysis[,ind] <- sample
#     #convert mu.a and pa
#     ind <- intersect(ind, H$H.ind)
#     mu.a[ind] <- L$update$mua
#     Pa[ind, ind] <- L$update$pa
#   }
#   return(list(mu.f = mu.f,
#               Pf = Pf,
#               mu.a = mu.a,
#               Pa = Pa,
#               analysis = analysis))
# }

##' @title block.2.vector
##' @name  block.2.vector
##' @author Dongchen Zhang
##' 
##' @param block.list  lists of blocks generated by the `build.block.xy` function.
##' @param X A matrix contains ensemble forecasts.
##' @param H H index created by the `construct_nimble_H` function.
##' @param adjustment logical variable determine if we want to adjust the analysis ensembles based on likelihood.
##' 
##' @return It returns a list of analysis results by MCMC sampling.
# block.2.vector <- function (block.list, X, H, adjustment) {
#   # initialize site.ids, mu.a, mu.f, pa, and pf.
#   site.ids <- attributes(X)$Site
#   mu.f <- mu.a <- c()
#   Pf <- Pa <- matrix(0, length(site.ids), length(site.ids))
#   analysis <- X
#   # loop over blocks.
#   for (L in block.list) {
#     # grab index for the locations within the current block.
#     ind <- c()
#     for (id in L$site.ids) {
#       ind <- c(ind, which(site.ids == id))
#     }
#     # grab mu.a, mu.f, pa, and pf from the MCMC updates.
#     mu.a[ind] <- L$update$mufa
#     Pa[ind, ind] <- L$update$pfa
#     mu.f[ind] <- L$data$muf
#     Pf[ind, ind] <- L$data$pf
#     # adjustment.
#     if (as.logical(adjustment)) {
#       sample <- as.data.frame(adj.ens(Pf[ind, ind], X[,ind], mu.f[ind], mu.a[ind], Pa[ind, ind]))
#     } else {
#       sample <- as.data.frame(mvtnorm::rmvnorm(nrow(X), L$update$mufa, L$update$pfa, method = "svd"))
#     }
#     analysis[,ind] <- sample
#   }
#   return(list(mu.f = mu.f,
#               Pf = Pf,
#               mu.a = mu.a,
#               Pa = Pa,
#               analysis = analysis))
# }

########## Updated block.2.vector：：：：： mufa version
block.2.vector <- function(block.list, X, H = NULL) {
  site.ids_all <- attributes(X)$Site
  
  mu.f <- rep(NA_real_, ncol(X))
  mu.a <- rep(NA_real_, ncol(X))
  Pf <- matrix(0, ncol(X), ncol(X))
  Pa <- matrix(0, ncol(X), ncol(X))
  
  analysis <- X
  
  for (L in block.list) {
    ind <- c()
    for (id in L$site.ids) {
      ind <- c(ind, which(site.ids_all == id))
    }
    
    mu.f[ind] <- L$data$muf
    Pf[ind, ind] <- L$data$pf
    
    mu.a[ind] <- L$update$mufa  
    Pa[ind, ind] <- L$update$pfa 
    
    sample <- mvtnorm::rmvnorm(
      n = nrow(X),
      mean = L$update$mufa, 
      sigma = L$update$pfa, 
      method = "svd"
    )
    
    analysis[, ind] <- sample
  }
  
  return(list(
    mu.f = mu.f,
    Pf = Pf,
    mu.a = mu.a,
    Pa = Pa,
    analysis = analysis
  ))
}

sda.forecast.local <- function (settings, obs.mean, obs.cov, Q = NULL, pre_enkf_params = NULL, 
                                ensemble.samples = NULL, outdir = NULL, control = list(TimeseriesPlot = FALSE, 
                                                                                       OutlierDetection = FALSE, send_email = NULL, keepNC = TRUE, 
                                                                                       forceRun = TRUE, MCMC.args = NULL)) 
{
  if (future::supportsMulticore()) {
    future::plan(future::multicore)
  }
  else {
    future::plan(future::multisession)
  }
  cores <- as.numeric(settings$state.data.assimilation$batch.settings$general.job$cores)
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
    if (cores < 1) 
      cores <- 1
  }
  if (!is.null(outdir)) {
    PEcAn.logger::logger.info(paste0("Replacing model output directories with ", 
                                     outdir, "."))
    PEcAn.logger::logger.info("Please note that the workflow will only work locally.")
    PEcAn.logger::logger.info("Please swap the SDA function to `sda.enkf.multisite` function if you would like to run jobs through remote server.")
    settings$outdir <- outdir
    settings$rundir <- file.path(outdir, "run")
    settings$modeloutdir <- file.path(outdir, "out")
    settings$host$folder <- file.path(outdir, "out")
    settings$host$outdir <- file.path(outdir, "out")
    settings$host$rundir <- file.path(outdir, "run")
  }
  adjustment <- settings$state.data.assimilation$adjustment
  model <- settings$model$type
  defaults <- settings$pfts
  outdir <- settings$modeloutdir
  rundir <- settings$host$rundir
  nens <- as.numeric(settings$ensemble$size)
  var.names <- sapply(settings$state.data.assimilation$state.variable, 
                      "[[", "variable.name")
  names(var.names) <- NULL
  restart.list <- NULL
  if (!dir.exists(settings$outdir)) 
    dir.create(settings$outdir, showWarnings = FALSE)
  interval <- NULL
  state.interval <- cbind(as.numeric(lapply(settings$state.data.assimilation$state.variables, 
                                            "[[", "min_value")), as.numeric(lapply(settings$state.data.assimilation$state.variables, 
                                                                                   "[[", "max_value")))
  rownames(state.interval) <- var.names
  conf.settings <- settings
  site.ids <- conf.settings %>% purrr::map(~.x[["run"]]) %>% 
    purrr::map("site") %>% purrr::map("id") %>% base::unlist() %>% 
    base::as.character()
  # print(site.ids)
  # print("bk1")
  site.locs <- conf.settings %>% purrr::map(~.x[["run"]]) %>% 
    purrr::map("site") %>% purrr::map(function(s) {
      temp <- as.numeric(c(s$lon, s$lat))
      names(temp) <- c("Lon", "Lat")
      temp
    }) %>% dplyr::bind_rows() %>% as.data.frame() %>% `rownames<-`(site.ids)
  # print("bk1.5")
  start.cut <- lubridate::ymd_hms(settings$state.data.assimilation$start.date, 
                                  truncated = 3)
  Start.year <- (lubridate::year(settings$state.data.assimilation$start.date))
  End.year <- lubridate::year(settings$state.data.assimilation$end.date)
  assim.sda <- Start.year:End.year
  obs.mean <- obs.mean[sapply(lubridate::year(names(obs.mean)),
                              function(obs.year) obs.year %in% (assim.sda))]
  obs.cov <- obs.cov[sapply(lubridate::year(names(obs.cov)),
                            function(obs.year) obs.year %in% (assim.sda))]
  obs.times <- names(obs.mean)
  # print("bk2")
  obs.times.POSIX <- lubridate::ymd_hms(obs.times)
  for (i in seq_along(obs.times)) {
    if (is.na(obs.times.POSIX[i])) {
      if (is.na(lubridate::ymd(obs.times[i]))) {
        PEcAn.logger::logger.warn("Error: no dates associated with observations")
      }
      else {
        obs.times.POSIX[i] <- lubridate::ymd_hms(paste(obs.times[i], 
                                                       "23:59:59"))
      }
    }
  }
  # print("bk3")
  obs.times <- obs.times.POSIX
  read_restart_times <- c(lubridate::ymd_hms(start.cut, truncated = 3), 
                          obs.times)
  nt <- length(obs.times)
  if (nt == 0) 
    PEcAn.logger::logger.severe("There has to be at least one Obs.")
  do.call("library", list(paste0("PEcAn.", model)))
  my.write_restart <- paste0("write_restart.", model)
  ### Modify
  # my.read_restart <- paste0("read_restart.", model)
  my.read_restart <- read.restart.SIPNET
  
  my.split_inputs <- paste0("split_inputs.", model)
  if (is.null(adjustment)) 
    adjustment <- TRUE
  register.xml <- system.file(paste0("register.", model, ".xml"), 
                              package = paste0("PEcAn.", model))
  register <- XML::xmlToList(XML::xmlParse(register.xml))
  no_split <- !as.logical(register$exact.dates)
  if (!exists(my.split_inputs) & !no_split) {
    PEcAn.logger::logger.warn(my.split_inputs, "does not exist")
    PEcAn.logger::logger.severe("please make sure that the PEcAn interface is loaded for", 
                                model)
    PEcAn.logger::logger.warn(my.split_inputs, "If your model does not need the split function you can specify that in register.Model.xml in model's inst folder by adding <exact.dates>FALSE</exact.dates> tag.")
  }
  if (!file.exists(paste0(settings$outdir, "/Extracted_met/"))) {
    dir.create(paste0(settings$outdir, "/Extracted_met/"))
  }
  PEcAn.logger::logger.info("Splitting mets!")
  conf.settings <- conf.settings %>% `class<-`(c("list")) %>% 
    furrr::future_map(function(settings) {
      library(paste0("PEcAn.", settings$model$type), character.only = TRUE)
      inputs.split <- list()
      if (!no_split) {
        for (i in 1:length(settings$run$inputs$met$path)) {
          settings$run$inputs$met$path[[i]] <- do.call(my.split_inputs, 
                                                       args = list(settings = settings, start.time = lubridate::ymd_hms(settings$run$site$met.start, 
                                                                                                                        truncated = 3), stop.time = lubridate::ymd_hms(settings$run$site$met.end, 
                                                                                                                                                                       truncated = 3), inputs = settings$run$inputs$met$path[[i]], 
                                                                   outpath = paste0(paste0(settings$outdir, 
                                                                                           "/Extracted_met/"), settings$run$site$id), 
                                                                   overwrite = F))
          settings$run$start.date <- lubridate::ymd_hms(settings$state.data.assimilation$start.date, 
                                                        truncated = 3)
          settings$run$end.date <- lubridate::ymd_hms(settings$state.data.assimilation$end.date, 
                                                      truncated = 3)
        }
      }
      else {
        inputs.split <- inputs
      }
      settings
    }, .progress = F)
  conf.settings <- PEcAn.settings::as.MultiSettings(conf.settings)
  if (is.null(ensemble.samples)) {
    load(file.path(settings$outdir, "samples.Rdata"))
  }
  new.params <- sda_matchparam(settings, ensemble.samples, 
                               site.ids, nens)
  samp <- conf.settings$ensemble$samplingspace
  parents <- lapply(samp, "[[", "parent")
  order <- names(samp)[lapply(parents, function(tr) which(names(samp) %in% 
                                                            tr)) %>% unlist()]
  samp.ordered <- samp[c(order, names(samp)[!(names(samp) %in% 
                                                order)])]
  inputs <- vector("list", length(conf.settings))
  for (s in seq_along(conf.settings)) {
    if (is.null(inputs[[s]])) {
      inputs[[s]] <- list()
    }
    for (i in seq_along(samp.ordered)) {
      inputs[[s]][[names(samp.ordered)[i]]] <- PEcAn.uncertainty::input.ens.gen(settings = conf.settings[[s]], 
                                                                                input = names(samp.ordered)[i], method = samp.ordered[[i]]$method, 
                                                                                parent_ids = NULL)
    }
  }
  for (t in 1:nt) {
    sda.outputs <- FORECAST <- enkf.params <- ANALYSIS <- ens_weights <- list()
    obs.t <- as.character(lubridate::date(obs.times[t]))
    obs.year <- lubridate::year(obs.t)
    PEcAn.logger::logger.info(paste("Processing Year:", obs.year))
    if (t > 1) {
      PEcAn.logger::logger.info("Splitting mets!")
      inputs.split <- furrr::future_pmap(list(conf.settings %>% 
                                                `class<-`(c("list")), inputs, model), function(settings, 
                                                                                               inputs, model) {
                                                  library(paste0("PEcAn.", model), character.only = TRUE)
                                                  inputs.split <- inputs
                                                  if (!no_split) {
                                                    for (i in seq_len(nens)) {
                                                      inputs.split$met$samples[i] <- do.call(my.split_inputs, 
                                                                                             args = list(settings = settings, start.time = (lubridate::ymd_hms(obs.times[t - 
                                                                                                                                                                           1], truncated = 3) + lubridate::second(lubridate::hms("00:00:01"))), 
                                                                                                         stop.time = lubridate::ymd_hms(obs.times[t], 
                                                                                                                                        truncated = 3), inputs = inputs$met$samples[[i]]))
                                                    }
                                                  }
                                                  else {
                                                    inputs.split <- inputs
                                                  }
                                                  inputs.split
                                                })
      PEcAn.logger::logger.info("Collecting restart info!")
      restart.list <- furrr::future_pmap(list(out.configs, 
                                              conf.settings %>% `class<-`(c("list")), params.list, 
                                              inputs.split), function(configs, settings, new.params, 
                                                                      inputs) {
                                                new_state_site <- new.state[, which(attr(X, "Site") %in% 
                                                                                      settings$run$site$id)]
                                                if (is.vector(new_state_site)) {
                                                  new_state_site <- matrix(new_state_site)
                                                }
                                                list(runid = configs$runs$id, start.time = strptime(obs.times[t - 
                                                                                                                1], format = "%Y-%m-%d %H:%M:%S") + lubridate::second(lubridate::hms("00:00:01")), 
                                                     stop.time = strptime(obs.times[t], format = "%Y-%m-%d %H:%M:%S"), 
                                                     settings = settings, new.state = new_state_site, 
                                                     new.params = new.params, inputs = inputs, RENAME = TRUE, 
                                                     ensemble.id = settings$ensemble$ensemble.id)
                                              })
    }
    else {
      restart.list <- vector("list", length(conf.settings))
    }
    gc()
    PEcAn.logger::logger.info("Writting configs!")
    ####
    # print("Debugging out.configs")
    ####################This is supposed to be the bugging part
    out.configs <- furrr::future_pmap(list(conf.settings %>% 
                                             `class<-`(c("list")), restart.list, inputs), function(settings, 
                                                                                                   restart.arg, inputs) {
                                               library(paste0("PEcAn.", settings$model$type), character.only = TRUE)
                                               write.ensemble.configs(defaults = settings$pfts, 
                                                                      ensemble.samples = ensemble.samples, settings = settings, 
                                                                      model = settings$model$type, write.to.db = settings$database$bety$write, 
                                                                      restart = restart.arg, samples = inputs, rename = TRUE)
                                             }) %>% stats::setNames(site.ids)
    ####################End
    # print("finish de-bugging ")
    ensemble.ids <- site.ids %>% furrr::future_map(function(i) {
      run.list <- c()
      for (j in 1:nens) {
        run.list <- c(run.list, paste0("ENS-", sprintf("%05d", 
                                                       j), "-", i))
      }
      return(run.list)
    }, .progress = F) %>% unlist
    runs.tmp <- file.path(rundir, ensemble.ids)
    PEcAn.logger::logger.info("Running models!")
    job.files <- file.path(runs.tmp, "job.sh")
    temp <- job.files %>% furrr::future_map(function(f) {
      cmd <- paste0("cd ", dirname(f), ";./job.sh")
      system(cmd, intern = F, ignore.stdout = T, ignore.stderr = T)
    }, .progress = F)
    PEcAn.logger::logger.info("Reading forecast outputs!")
    
    #### Having the problems (Solved)
    reads <- build_X(out.configs = out.configs, settings = settings, 
                     new.params = new.params, nens = nens, read_restart_times = read_restart_times, 
                     outdir = outdir, t = t, var.names = var.names, my.read_restart = my.read_restart, 
                     restart_flag = FALSE)
    
    # print("breakpoint11")
    params.list <- reads %>% purrr::map(~.x %>% purrr::map("params"))
    X <- reads %>% furrr::future_map(function(r) {
      r %>% purrr::map_df(~.x[["X"]] %>% t %>% as.data.frame)
    })
    if (control$OutlierDetection) {
      X <- outlier.detector.boxplot(X)
      PEcAn.logger::logger.info("Outlier Detection.")
    }
    # print("breakpoint12")
    X <- seq_along(X) %>% furrr::future_map(function(i) {
      temp <- do.call(cbind, X[i])
      colnames(temp) <- paste0(var.names, ".", i)
      return(temp)
    }) %>% dplyr::bind_cols() %>% `colnames<-`(c(rep(var.names, 
                                                     length(X)))) %>% `attr<-`("Site", c(rep(site.ids, 
                                                                                             each = length(var.names))))
    FORECAST[[obs.t]] <- X
    gc()
    
    
    #### TO Forecasting or to reanalysis????
    if (!is.null(obs.mean[[t]][[1]]) || (as.logical(settings$state.data.assimilation$free.run) &
                                         control$forceRun)) {
      ##### Initialzation
      if (t == 1 || !exists("block.list.all")) {
        block.list.all <- obs.mean %>% purrr::map(function(l) {
          NULL
        })
      }
      if (is.null(control$MCMC.args)) {
        MCMC.args <- list(niter = 1e+05, nthin = 10,
                          nchain = 4, nburnin = 50000)
      }
      else {
        MCMC.args <- control$MCMC.args
      }
      settings$state.data.assimilation$batch.settings$analysis <- NULL
      ##### DO THE KALMAN Calculation
      # enkf.params[[obs.t]] <- analysis_sda_block(settings,
      #                                            block.list.all, X, obs.mean, obs.cov, t, nt,
      #                                            MCMC.args, pre_enkf_params)
      # enkf.params[[obs.t]] <- c(enkf.params[[obs.t]], RestartList = list(restart.list %>%
      #                                                                      stats::setNames(site.ids)))
      
      ##### Revised
      enkf.params[[obs.t]] <- analysis_sda_block(settings,
                                                 block.list.all, X, obs.mean, obs.cov, t, nt,
                                                 MCMC.args, pre_enkf_params)
      
      ## ---- write MCMC diagnostics to outdir ----
      diag_csv_path <- file.path(settings$outdir, paste0("mcmc_diag_t", t, ".csv"))
      utils::write.csv(enkf.params[[obs.t]]$diag_table, diag_csv_path, row.names = FALSE)
      PEcAn.logger::logger.info(paste0("Saved MCMC diagnostics to: ", diag_csv_path))
      
      enkf.params[[obs.t]] <- c(enkf.params[[obs.t]], RestartList = list(restart.list %>%
                                                                           stats::setNames(site.ids)))
      
      
      block.list.all <- enkf.params[[obs.t]]$block.list.all
      mu.f <- enkf.params[[obs.t]]$mu.f
      Pf <- enkf.params[[obs.t]]$Pf
      Pa <- enkf.params[[obs.t]]$Pa
      mu.a <- enkf.params[[obs.t]]$mu.a
    }
    #### Release some space
    if (!is.null(enkf.params[[obs.t]])) {
      # 有同化结果 → 用同化后的分析场
      analysis <- enkf.params[[obs.t]]$analysis
    } else {
      # 没有观测（或 free.run） → 退回纯 forecast
      analysis <- FORECAST[[obs.t]]
    }
    ## ===== OBS → SIPNET (before restart) =====
    NEE_cols <- which(colnames(analysis) == "NEE")
    if (length(NEE_cols) > 0) {
      analysis[, NEE_cols] <-
        nee_obs_to_model(analysis[, NEE_cols])
    }
    
    # enkf.params[[obs.t]]$analysis <- NULL
    ##### Make the analysis value within the range of the min and max
    for (i in 1:ncol(analysis)) {
      int.save <- state.interval[which(startsWith(colnames(analysis)[i], 
                                                  var.names)), ]
      analysis[analysis[, i] < int.save[1], i] <- int.save[1]
      analysis[analysis[, i] > int.save[2], i] <- int.save[2]
    }
    ##### Update the new state
    new.state <- as.data.frame(analysis)
    ANALYSIS[[obs.t]] <- analysis
    ### Calculate the ensemble weight
    # ens_weights[[obs.t]] <- PEcAnAssimSequential::sda_weights_site(FORECAST, 
    #                                                                ANALYSIS, 1, nens)
    ### Pure Forecasting, weight is averaged
    ens_weights[[obs.t]] <- matrix(1/nens, nrow = nens, ncol = ncol(ANALYSIS[[obs.t]]))
    ##### Save and output
    sda.outputs <- list(
      obs.mean = obs.mean[[t]],
      obs.cov = obs.cov[[t]],
      forecast = FORECAST[[obs.t]],
      analysis = ANALYSIS[[obs.t]], 
      enkf.params = enkf.params[[obs.t]],
      ens_weights = ens_weights[[obs.t]], 
      params.list = params.list, 
      restart.list = restart.list)
    save(sda.outputs, file = file.path(settings$outdir, paste0("sda.output", 
                                                               t, ".Rdata")))
    ##### Delete NC file to save some space if keepNC is FALSE
    if (!(control$keepNC) && t == 1) {
      PEcAn.logger::logger.info("Deleting NC files!")
      outs.tmp <- file.path(outdir, ensemble.ids)
      temp <- outs.tmp %>% furrr::future_map(function(f) {
        temp <- list.files(f, "*.nc", full.names = T)
        unlink(temp)
      }, .progress = F)
    }
    ##### send email to progress
    if (!is.null(control$send_email)) {
      sendmail <- Sys.which("sendmail")
      mailfile <- tempfile("mail")
      cat(paste0("From: ", control$send_email$from, "\n", 
                 "Subject: ", "SDA progress report", "\n", "To: ", 
                 control$send_email$to, "\n", "\n", paste("Time point:", 
                                                          obs.times[t], "has been completed!")), file = mailfile)
      system2(sendmail, c("-f", paste0("\"", control$send_email$from, 
                                       "\""), paste0("\"", control$send_email$to, "\""), 
                          "<", mailfile))
      unlink(mailfile)
    }
  }
  ##### After all time is looped, summarize all output to 2 lists
  sda.out.files <- file.path(settings$outdir, paste0("sda.output", 
                                                     1:nt, ".Rdata"))
  analysis.all <- forecast.all <- vector("list", nt)
  for (file in seq_along(sda.out.files)) {
    res_env <- new.env()
    load(sda.out.files[file], envir = res_env)
    analysis.all[[file]] <- res_env$sda.outputs$analysis
    forecast.all[[file]] <- res_env$sda.outputs$forecast
  }
  names(analysis.all) <- as.character(lubridate::date(obs.times))
  names(forecast.all) <- as.character(lubridate::date(obs.times))
  save(list = c("analysis.all", "forecast.all"), file = file.path(settings$outdir, 
                                                                  "sda.all.forecast.analysis.Rdata"))
  gc()
}

nee_model_to_obs <- function(x) {
  x * 1e8 / 1.157407
}

nee_obs_to_model <- function(x) {
  x * 1.157407 / 1e8
}