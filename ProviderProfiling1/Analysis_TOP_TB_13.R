### TOP Analysis function for:

#Title:
#Investigating risk adjustment methods for health care provider profiling
#when observations are scarce or events rare

#Authors: 
#TB Brakenhoff, KGM Moons, J Kluin, RHH Groenwold


#Date created: 16-05-2015



### Set working directory

### Load  empirical data (ANALYSIS OF EMPIRICAL DATA NOT PROVIDED PUBLICLY)
#Obtained from the NVT (http://www.nvtnet.nl/)
#Cannot be made public due to privacy concerns.

### Wipe Workspace

rm(list = ls())

### SOURCE functions

#Load BASE functions
source("Functions\\Functions - General\\fun_recpred.R")
source("Functions\\Functions - General\\fun_spec_dec.R")
source("Functions\\Functions - General\\fun_Factorize.R")
source("Functions\\Functions - General\\fun_FACTOR.R")
source("Functions\\Functions - General\\fun_min_yint.R")

#Load Total simulation function
source("Functions\\Functions - Simulation\\fun_zh_sampling22_sim.R")

#Load Data simulations function (10 confs)
source("Functions\\Functions - Simulation\\\\fun_data_simulation8.R")

#Load Methods functions
source("Functions\\Functions - Simulation\\fun_LRana_raw.R")
source("Functions\\Functions - Simulation\\fun_LRana_firth.R")
source("Functions\\Functions - Simulation\\fun_LRana_glm.R")
source("Functions\\Functions - Simulation\\fun_LRana_svy.R")
source("Functions\\Functions - Simulation\\fun_PS_wei19.R")
source("Functions\\Functions - Simulation\\fun_PS_wei_trim19.R")
source("Functions\\Functions - Simulation\\fun_PS_cov16.R")
source("Functions\\Functions - Simulation\\fun_PS_mat22_alt.R")

#LOAD FUNCTIONS for extra scenario of simulation (# of providers.)
source("Functions\\Functions - Simulation\\HPC_MoreProviderScen08-2016\\fun_refvals_ex.R")
source("Functions\\Functions - Simulation\\HPC_MoreProviderScen08-2016\\fun_data_simulation_ex.R")
source("Functions\\Functions - Simulation\\HPC_MoreProviderScen08-2016\\fun_LRana_svy_ex.R")
source("Functions\\Functions - Simulation\\HPC_MoreProviderScen08-2016\\fun_PS_cov_ex.R")
source("Functions\\Functions - Simulation\\HPC_MoreProviderScen08-2016\\fun_PS_wei_ex.R")
source("Functions\\Functions - Simulation\\HPC_MoreProviderScen08-2016\\fun_PS_wei_trim_ex.R")
source("Functions\\Functions - Simulation\\HPC_MoreProviderScen08-2016\\fun_zh_sampling_sim_ex2.R")

#LOAD FUNCTIONS for Empirical Analysis with 16 Provider! (NOT PROVIDED IN REPOSITORY)

#Load plotting function
source("Functions\\Functions - Simulation\\fun_zh_plots25_sim.R")

#Load reference values function
source("Functions\\Functions - Simulation\\fun_refvals.R")

### Install Packages

#install.packages("nnet")
#install.packages("ipw")
#install.packages("data.table")
#install.packages("sampling")
#install.packages("microbenchmark")
#install.packages("MatchIt")
#install.packages("brglm")
#install.packages('sp')
#install.packages('rgeos')
#install.packages('logistf')
#install.packages('ggplot2')
#install.packages('mice')
#install.packages('rms')
#install.packages('concreg')
#install.packages("dummies")
#install.packages("MASS")
#install.packages("MatchIt")
#install.packages("optmatch")
#install.packages("clue")
#install.packages("fields")
#install.packages("plyr")
#install.packages("abind")
#install.packages("psych")
#install.packages("epitools")
#install.packages("prodlim")
#install.packages("Rmisc")
install.packages("extrafont")

### Require Packages

require(ModelGood)
require(ggplot2)
require(gmodels)
require(xtable)
require(survey)
require(plyr)   #Think if this package is actually necessary
require(dplyr)
require(nnet)
require(data.table)
require(sampling)
require(microbenchmark)
require(MatchIt)
require(optmatch)
require(brglm)
require(sp)
require(rgeos)
require(logistf)
library(dummies)
require(clue)
require(fields)
require(abind)
require(psych)
require(epitools)
require(prodlim)
require(Rmisc)
require(extrafont)
library(showtext)

#Make Times New Roman available
windowsFonts(Times=windowsFont("TT Times New Roman"))

#font.add("TT Times New Roman", regular = "times.ttf",
#         bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")


##############           ##############
#######        ANALYSIS         #######
##############           ##############

############################
## LIST OF SIM DATA CHARS ##
############################

#df of coefficients (updated 04-11-15)
b.coefs <- data.frame(b.y0 = 2.29,
                      b0A = 0.02,
                      b0C = 0.02,
                      b1A = 1,
                      b1C = 1,
                      byXA = -0.5,
                      byXB = -1,
                      byXC = 0.5,
                      byz  = 1)

#list of functions
funs         <- c(#LR.raw=LRana_raw,
                  LR=LRana_svy,
                  #LR.firth=LRana_firth,
                  PS.cov=PS_cov_new,
                  PS.wei=PS_wei_new,
                  PS.wei.trim=PS_wei_new_trim,
                  PS.mat=PS_mat)
                  

#df of scenarios
sim.scens <- data.frame(zhperc=c(0.05,0.1,0.2,0.5,rep(1,8)),
                        by0=c(rep(2.29,5),1,2,4,5,rep(2.29,3)),
                        b0s=c(rep(0.02,9),-2.34,-0.69,1.87))


#TRUE Reference values (different for by0s)
mat.con <- matrix(c(-0.5,0.5),ncol=nrow(sim.scens),nrow=2)

#MARGINAL for cor=0 FIXED 07-01-16
mat.mar <- matrix(c(rep(c(-0.497,0.494),5),
                    c(-0.491,0.489),
                    c(-0.496,0.495),
                    c(-0.497,0.500),
                    c(-0.496,0.500),
                    rep(c(-0.497,0.494),3)),
                  ncol=nrow(sim.scens),nrow=2)


#FIXED 07-01-16
mat.mar.01 <- matrix(c(rep(c(-0.492,0.490),5),
                       c(-0.483,0.481),
                       c(-0.491,0.488),
                       c(-0.498,0.498),
                       c(-0.503,0.500),
                       rep(c(-0.492,0.490),3)),
                     ncol=nrow(sim.scens),nrow=2)

#MATCHING  for cor=0 FIXED 07-01-16
mat.mat <- matrix(c(rep(c(-0.491,0.500),5),
                    c(-0.487,0.497),
                    c(-0.491,0.500),
                    c(-0.494,0.506),
                    c(-0.497,0.506),
                    rep(c(-0.491,0.500),3)),
                  ncol=nrow(sim.scens),nrow=2)


#List of true values (ORDER ACCORDING TO FUNS!)
zh.tot.coef.0 <- c(rep(list(mat.con),2),
                 rep(list(mat.mar),3))#,
                 #list(mat.mat))
                  
zh.tot.coef.01 <- c(rep(list(mat.con),2),
                   rep(list(mat.mar.01),3))#,


#parameters
n.max       <- 10000
n.zh        <- 3

############################################################
## LIST OF SIM DATA CHARS FOR EXTRA SCENARIO (#PROVIDERS) ##
############################################################

#required packages
require("MASS")
require("nnet")
library(dummies)
require("survey")
require(microbenchmark)

### PARAMETERS ###
n.perzh     <- 3000 #Now specify amount of people per zh (n.max is dependent on sample size)
seed        <- 100
fit         <- T
ref         <- "A"
cors        <- 0
confs       <- 10
provs       <- 10
b.yi        <- 0.10 
trim        <- 0.02
sims        <- 2

### List of functions
funs         <- c(LR=LRana_svy_ex,
                  PS.cov=PS_cov_ex,
                  PS.wei=PS_wei_ex,
                  PS.wei.trim=PS_wei_trim_ex)


### To determine marginal values for the marginal methods for different numbers of providers (n.perzh=10000, sims=6000)
ress <- refvals_ex(provs=c(5,10,15,20),n.perzh=10000,ref="A",cors=0,confs=10,b.yi=0.10,sims=6000)

  #Save results to file
  save(ress,file="EX_SIM_REF_VALS_6000.Rdata")
  load(file="EX_SIM_REF_VALS_6000.Rdata")
  
  #time function
  ll <-microbenchmark(r1=refvals_ex(provs=c(5,10,15,20),n.perzh=3000,ref="A",cors=0,confs=10,b.yi=0.10,sims=20),times=5L)

  ### From ress extra coefs.out data (LIST WITH EACH 4 rows corresponding to true values for each method)
  prov.list.mar <- lapply(ress$mean.sd.list,function(x)x[c(1,1,2,2),])

  
### RUN SIMULATION FUNCTION with marginal ref values from prov.list.mar 
  
#provs=10  
sim.ex2          <- zh_sampling_ex(provs=provs,coefs.mar=prov.list.mar$prov10,b.yi=b.yi,n.perzh=n.perzh,funs=funs,sims=sims,
                                   seed=seed,trim=trim,ref="A",confs=confs,cors=cors)
  
#For the HPC make 4 seperate files that act as input for each different prov  

#################################
## SIM DATA SAMPLING FUNCTION ###
#################################

#with corr=0
SIM.TOTS.0 <- zh_sampling_new(n.zh,n.max,b.coefs,funs,sim.scens,sims=2000,seed=13,trim=0.02,
                            l.plots=T,ref="B",sav=T,plot.t=1:4,zh.tot.coef=zh.tot.coef.0,confs=10,cors=0)


#SAVE IMAGE
save.image("SIMS_08_01_16_corr0.RData")


#With corr=0.1
SIM.TOTS.01 <- zh_sampling_new(n.zh,n.max,b.coefs,funs,sim.scens,sims=2000,seed=13,trim=0.02,
                               l.plots=T,ref="B",sav=T,plot.t=1:4,zh.tot.coef=zh.tot.coef.01,confs=10,cors=0.1)

#SAVE IMAGE
save.image("SIMS_08_01_16_corr01.RData")



# Benchmarking
ll <- microbenchmark(SIM.TOTS=zh_sampling_new(n.zh,n.max,b.coefs,funs,sim.scens,sims=10,seed=11,trim=0.02,
                                              l.plots=T,ref="B",sav=T,plot.t=1:4,zh.tot.coef=zh.tot.coef,bin.inc=0.01),times=1L)

#Check results
cbind(SIM.TOTS.0$sim.scens,SIM.TOTS.0$sim.chars)
cbind(SIM.TOTS.01$sim.scens,SIM.TOTS.01$sim.chars)

#X table results
xtable(cbind(SIM.TOTS.01$sim.scens,SIM.TOTS.01$sim.chars),digits = 2)
xtable(cbind(SIM.TOTS.0$sim.scens,SIM.TOTS.0$tot.bias),digits = 2)
xtable(cbind(SIM.TOTS.0$sim.scens,SIM.TOTS.0$MSE.mat),digits = 2)
xtable(cbind(SIM.TOTS.0$sim.scens,SIM.TOTS.0$cov.mat),digits = 2)
xtable(cbind(SIM.TOTS.0$sim.scens,SIM.TOTS.0$rat.se.sd),digits = 2)


#Different scenarios
xtable(cbind(SIM.TOTS$sim.scens,SIM.TOTS$sim.chars))

###############################
### SEPERATE PLOTS FUNCTION ###
###############################

#Results of simulation (.Rdata) can be found in: Final Results Files - R\SIMULATION

# INPUT PARS FOR SEPERATE PLOT FUNCTION corr=0
tot.bias=SIM.TOTS.0$tot.bias
MSE.mat=SIM.TOTS.0$MSE.mat
cov.mat=SIM.TOTS.0$cov.mat
rat.se.sd=SIM.TOTS.0$rat.se.sd
sim.scens=SIM.TOTS.0$sim.scens
sim.chars=SIM.TOTS.0$sim.chars
zh.afk=SIM.TOTS.0$zh.afk
n.max=10000
save=T
funs=SIM.TOTS.0$funs
cors=0

# INPUT PARS FOR SEPERATE PLOT FUNCTION corr=0.1
# tot.bias=SIM.TOTS.01$tot.bias
# MSE.mat=SIM.TOTS.01$MSE.mat
# cov.mat=SIM.TOTS.01$cov.mat
# rat.se.sd=SIM.TOTS.01$rat.se.sd
# sim.scens=SIM.TOTS.01$sim.scens
# sim.chars=SIM.TOTS.01$sim.chars
# zh.afk=SIM.TOTS.01$zh.afk
# n.max=10000
# save=T
# funs=SIM.TOTS.01$funs
# cors=0.1

#Run function
JOJO <- zh_plots(tot.bias,MSE.mat,cov.mat,rat.se.sd,sim.scens,sim.chars,n.max,zh.afk,save=T,funs,plot.types=1:4,cors=cors)

JOJO$Bias_b0s_0.pdf

#################################################################
## Simulation for raw effects provided in text Methods section ##
#################################################################

n.max=nmax
nmax=1e5
ref="B"
confs=10
cors=0.1

col.raw <- NULL

#Loop for raw effects with cors=0
for(g in 1:20){
  
  cat(g)
    
  dat.raw  <- sim.data(b.coefs,nmax,seed=NULL,fit=F,ref=ref,confs=10,cors=cors)
  
  ana.mar  <- coefficients(glm(overleden_zh~relevel(zh_afk,ref),data=dat.raw$df.reg,family=binomial))[2:3]
  
  col.raw <- rbind(col.raw,ana.mar)
}

col.raw01 <- NULL

cors=0.1

for(g in 1:20){
  
  cat(g)
  
  dat.raw  <- sim.data(b.coefs,nmax,seed=NULL,fit=F,ref=ref,confs=10,cors=cors)
  
  ana.mar  <- coefficients(glm(overleden_zh~relevel(zh_afk,ref),data=dat.raw$df.reg,family=binomial))[2:3]
  
  col.raw01 <- rbind(col.raw,ana.mar)
}


apply(col.raw,2,function(x) c(mean(x),sd(x)))
apply(col.raw01,2,function(x) c(mean(x),sd(x)))

##########################################################
## FUNCTION TO DETERMINE REFERENCE VALS MAR AND MAT     ##
##########################################################

#To check if values for cors=0.1 are the same
refs <- refvals(b.coefs,vals=c(2.29,1,2,4,5),nmax=1e6,ref="B",sims=100,cors=0.1)
refs

write.csv(refs,file="refs01.csv")

nn <- read.csv("refs00.csv")

colnames(nn) <- c("X",rep(c("A","C"),10))

nn[,1:11]
nn[,c(1,12:21)]