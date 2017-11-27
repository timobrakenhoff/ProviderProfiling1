## HPC script for extra simulations of MIR1
## 4 Scenarios for # of providers. 
## 02-09-2016

provs       <- 20
sims        <- 2000
seed        <- 200

#HPC
setwd("/home/julius_te/tbrakenhoff/Rcode/MIR1/AnalysisR")

#Load General functions
source("FunctionsGeneral/fun_recpred.R")
source("FunctionsGeneral/fun_spec_dec.R")
source("FunctionsGeneral/fun_Factorize.R")
source("FunctionsGeneral/fun_FACTOR.R")
source("FunctionsGeneral/fun_min_yint.R")

#LOAD FUNCTIONS for extra scenario of simulation (# of providers.)
source("fun_data_simulation_ex.R")
source("fun_LRana_svy_ex.R")
source("fun_PS_cov_ex.R")
source("fun_PS_wei_ex.R")
source("fun_PS_wei_trim_ex.R")
source("fun_zh_sampling_sim_ex2.R")

#Load ress for marginal true values
load(file="EX_SIM_REF_VALS_6000.Rdata")

prov.list.mar <- lapply(ress$mean.sd.list,function(x)x[c(1,1,2,2),])


### Set working directory
setwd("/home/julius_te/tbrakenhoff/Rcode/MIR1/OutputR")

#required packages
require("MASS")
require("nnet")
library(dummies)
require("survey")
require(microbenchmark)

### PARAMETERS ###
n.perzh     <- 3000 #Now specify amount of people per zh (n.max is dependent on sample size)
fit         <- T
ref         <- "A"
cors        <- 0
confs       <- 10
b.yi        <- 0.10 
trim        <- 0.02

sav.im <- paste0("HPC_MIR1_prov",provs,".Rdata")

### List of functions
funs         <- c(LR=LRana_svy_ex,
                  PS.cov=PS_cov_ex,
                  PS.wei=PS_wei_ex,
                  PS.wei.trim=PS_wei_trim_ex)

## RUN FUNCTION 
SIM_RES_EXTRA <- zh_sampling_ex(provs=provs,coefs.mar=prov.list.mar[[paste0("prov",provs)]],b.yi=b.yi,n.perzh=n.perzh,
                                funs=funs,sims=sims,seed=seed,trim=trim,ref="A",confs=confs,cors=cors)


#Save results
save.image(sav.im)
