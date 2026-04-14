# loading libraries
library(dplyr)
library(xts)
library(PEcAn.all)
library(purrr)
library(furrr)
library(lubridate)
library(nimble)
library(ncdf4)
library(PEcAnAssimSequential)
library(dplyr)
library(sp)
library(raster)
library(zoo)
library(ggplot2)
library(mnormt)
library(sjmisc)
library(stringr)
library(doParallel)
library(doSNOW)
library(data.table)
library(Kendall)
library(lgarch)
library(parallel)
library(foreach)
library(terra)
setwd("/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/wishart_sda/")
source("/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/wishart_sda/functions_normalize.R")

## read settings xml file.
# settings_dir <- "/projectnb/dietzelab/guYANG/pecan/runners/test4/pecan_flux.xml"
# settings <- PEcAn.settings::read.settings(settings_dir)
load("/projectnb/dietzelab/guYANG/pecan/runners/test10/pecan_flux.RData")
# load observations.
load("/projectnb/dietzelab/guYANG/pecan/runners/test10/obs.mean.RData")
load("/projectnb/dietzelab/guYANG/pecan/runners/test10/obs.cov.RData")

#### Fix the multi output in one timestep bug
settings$model$jobtemplate <- "/projectnb/dietzelab/guYANG/pecan/runners/test7/sipnet_template.job"
settings$outdir <- "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q_3/"
settings$host$prerun <- "module load R/4.4.0"

load("/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/sda_idx.Rdata")
all_ids  <- vapply(settings, \(s) as.character(s$run$site$id), "")
settings <- settings[all_ids %in% keep_ids]
settings <- PEcAn.settings::as.MultiSettings(settings)
sub_obs <- function(L, keep) setNames(lapply(L, \(l) l[names(l) %in% keep]), names(L))
obs.mean <- sub_obs(obs.mean, keep_ids)
obs.cov  <- sub_obs(obs.cov,  keep_ids)

# replace zero observations and variances with small numbers.
for (i in 1:length(obs.mean)) {
  if(is.null(obs.mean[[i]][[1]])){
    next
  }
  for (j in 1:length(obs.mean[[i]])) {
    if (length(obs.mean[[i]][[j]])==0) {
      next
    }
    obs.mean[[i]][[j]][which(obs.mean[[i]][[j]]==0)] <- 0.01
    if(length(obs.cov[[i]][[j]]) > 1){
      diag(obs.cov[[i]][[j]])[which(diag(obs.cov[[i]][[j]]<=0.1))] <- 0.1
    }else{
      if(obs.cov[[i]][[j]] <= 0.1){
        obs.cov[[i]][[j]] <- 0.1
      }
    }
  }
}

# setup the batch job settings.
general.job <- list(cores = 28, folder.num = 80)
batch.settings = structure(list(
  general.job = general.job,
  qsub.cmd = "qsub -l h_rt=24:00:00 -l mem_per_core=4G -l buyin -pe omp @CORES@ -V -N @NAME@ -o @STDOUT@ -e @STDERR@ -S /bin/bash"
))
settings$state.data.assimilation$batch.settings <- batch.settings

# alter the ensemble size.
settings$ensemble$size <- 20

# update settings with the actual PFTs.
settings <- PEcAn.settings::prepare.settings(settings)


# load PFT parameter file.
load("/projectnb/dietzelab/dongchen/anchorSites/NA_runs/SDA_8k_site/samples.Rdata") ### Revised

settings$rundir <- "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q_3/run"

###### Change Q type
settings$state.data.assimilation$q.type <- "wishart"

###### free run or not?
# settings$state.data.assimilation$free.run <- "TRUE"
# obs.mean <- lapply(obs.mean, function(x) {
#   setNames(
#     vector("list", length(x)),  
#     names(x)
#   )
# })
# 
# obs.cov <- lapply(obs.cov, function(x) {
#   setNames(
#     vector("list", length(x)),
#     names(x)
#   )
# })

##########Problems update:
source("/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/wishart_sda/functions_normalize.R")
# execute the SDA.
sda.forecast.local(settings = settings,
                   obs.cov = obs.cov,
                   obs.mean = obs.mean, 
                   Q = NULL,
                   pre_enkf_params = NULL,
                   ensemble.samples = ensemble.samples,
                   outdir = "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/output_inter_q_3/",
                   control = list(TimeseriesPlot = FALSE,
                                  OutlierDetection=FALSE,
                                  send_email = NULL,
                                  debug.input.ens.gen = TRUE,
                                  keepNC = TRUE,
                                  forceRun = FALSE,
                                  MCMC.args = list(
                                    niter = 200000,   
                                    nburnin = 100000,
                                    nthin = 5,
                                    nchain = 4
                                  )))


# job_lines <- c(
#   "#!/bin/bash",
#   "module load R/4.4.0",
#   "Rscript /projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/wishart_sda/Block_sda_norm.R"
# )
# writeLines(job_lines, "/projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/logs/Block_sda_norm.sh")

# qsub -l h_rt=3:00:00 \
# -l buyin \
# -l mem_per_core=8G \
# -pe omp 28 \
# -V \
# -N Block_sda_norm \
# -o /projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/logs/Block_sda_norm.out \
# -e /projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/logs/Block_sda_norm.err \
# -M yanggu@bu.edu \
# -m abe \
# -S /bin/bash \
# /projectnb/dietzelab/guYANG/pecan/runners/wishart_sda/logs/Block_sda_norm.sh
