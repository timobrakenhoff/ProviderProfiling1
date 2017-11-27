#########################
## SAMPLING SIMULATION ##
#########################
## This function is now built for one scenario only ##
## Problem is that collection arrays are actually of different sizes for different # provs
## Now that problem is averted by just running scenario at a time
##


zh_sampling_ex <- function(provs,coefs.mar=NULL,b.yi,n.perzh,funs,sims=10,seed=123,trim=0.02,
                          ref="A",confs=10,cors=0){

#Set seed  
  if(is.null(seed)==F){set.seed(seed)}
  
#zh names
  zh.afk  <- LETTERS[LETTERS!=ref][1:(provs-1)]

#Initialize coefficient arrays in a list (Stil generalized for multiple scenarios)
  zh.coef <-  rep(list(array(NA,dim=c(provs-1,1,sims),dimnames=list(zh.afk,1:1,1:sims))),
                  length(funs))

#Initialize SE arrays in a list    
  zh.se   <-  rep(list(array(NA,dim=c(provs-1,1,sims),dimnames=list(zh.afk,1:1,1:sims))),
                  length(funs))

#Initialize CI arrays in a list    
  zh.ci   <-  rep(list(array(NA,dim=c((provs-1),2,sims),dimnames=list(zh.afk,c("lower","upper"),1:sims))),length(funs))
#Make list for chars of the samples
  nams.cols <- c(paste0("p.zh",LETTERS[1:provs]),paste0("py1.zh",LETTERS[1:provs]),"py1.tot")
  zh.chars  <-  array(NA,dim=c(1,provs*2+1,sims),dimnames=list(1:1,
                      nams.cols,1:sims))
#####
## START LOOP OVER SCENARIOS
#####

#specify coefficients
  #Coefficients for the multinom model with first provider as reference
  #Coefficients for the outcome model. Takes values from seq based on length provs 
  b.coefs.ex  <- c(0,rep(1,provs-1)) 
  coefs.out   <- c(1,seq(-1,1,length.out=provs-1)) 
  
  if(is.null(coefs.mar)){
  coefs.true <- matrix(coefs.out[-1],nrow = length(funs),ncol = length(coefs.out[-1]),byrow=T)  
  } else {coefs.true <- coefs.mar}
  
#Loop over sims  
  for (j in 1:sims){

#Percentage counter (progress bar)
  if(j%%(sims/10)==0){cat(j/sims*10)}


#####
## Sampling procedure per hospital (output is list)
#####

#Simulate data for specified scenario  
  samp.sim <- sim.data.ex(n.perzh=n.perzh,seed=NULL,fit=T,ref=ref,cors=cors,confs=confs,
                          provs=provs, b.yi=b.yi,b.coefs.ex=b.coefs.ex,coefs.out=coefs.out)

#Run analyses on samp.df
  zh.samp      <- lapply(funs,function(f) f(dat=samp.sim$df.reg,fit.ml=samp.sim$gen.ps,
                                            trim=trim,ref=ref,n.zh=provs-1))
#Collect coefficients   
  zh.samp.coef <- lapply(zh.samp, function(x) x$LR.coef)  

#Collect SEs
  zh.samp.SE   <- lapply(zh.samp, function(x) x$LR.se)

#Collect CIs
  zh.samp.CI   <- lapply(zh.samp, function(x) x$LR.ci)

#Store weight summaries for PS.wei methods IGNORE FOR NOW, SIZE OF MATRIX DEPENDS ON PROVS...
#  w.tab[i,,j]  <- c(zh.samp$PS.wei$sum.w,zh.samp$PS.wei.trim$sum.w)

#Store match chars
#  mat.char[i,,j] <- zh.samp$PS.mat$st.char

#Store chars of simulation
  zh.chars[1,,j] <- samp.sim$sim.chars

#Assign to elements of the list of coef and SE  
  for (k in 1:length(funs)){
    zh.coef[[k]][,1,j] <- zh.samp.coef[[k]]
    zh.se[[k]][,1,j]   <- zh.samp.SE[[k]]
    zh.ci[[k]][,,j]   <- zh.samp.CI[[k]]
  }   
  }



  
#####
## COLLECTION Results (optimized for only one row of sim.scens)
#####
  cat("\n","Collecting simulation results...","\n")

#Initialize collection matrices
  se.mat  <- coef.mat <- est.sd <- rat.se.sd <- est.ex <- tot.bias <- per.bias <- MSE.mat <- cov.mat <-
  matrix(NA,nrow=length(funs),ncol=provs-1,dimnames=list(names(funs),zh.afk))
  
#Loop over list elements to collect results  
  for (jj in 1:length(funs)){
  
    #Counter for filling of matrices
    #m.c <- seq(1,(length(funs)*2),2)[jj]
    
    #Mean SE, coef, and sd of coefficients
    se.mat[jj,]      <- apply(zh.se[[jj]],c(2,1),mean)
    coef.mat[jj,]    <- apply(zh.coef[[jj]],c(2,1),mean)
    est.sd[jj,]      <- apply(zh.coef[[jj]],c(2,1),sd)
    
    #Ratio of mean(SE) to SD of estimates
    rat.se.sd[jj,]   <- sqrt(apply(zh.se[[jj]]^2,c(2,1),mean))/est.sd[jj,]
    
    #Percentage of effect estimates extremer than 5 (methods are cbinded)
    est.ex[jj,]      <- apply(zh.coef[[jj]],c(2,1),function(x) sum(x>5 | x< -5)/length(x))
  
    #Bias per scenario and zh (one vector of coefs to compare with (NO DIFF BETWEEN REF FOR METHODS))
    tot.bias[jj,]    <- t(apply(zh.coef[[jj]],c(1,2),mean)-matrix(coefs.true[jj,]))
    
    #Percentage bias (WATCH OUT FOR SITUATION WHERE TRUE COEF IS 0...)
    per.bias[jj,]    <- t(100*tot.bias[jj,]/matrix(coefs.true[jj,]))
    
    #MSE matrix
    MSE.mat[jj,]     <- apply(zh.se[[jj]]^2,c(2,1),mean)+tot.bias[jj,]^2
    
    
    #Coverage matrix (cbind 2 vectors that calculate coverage for each coefficient over sims)
    cov.mat[jj,]     <- apply(zh.ci[[jj]][,1,]<= coefs.true[jj,] & coefs.true[jj,] <= zh.ci[[jj]][,2,],1,mean)
    
    }           

#Mean sim.chars over simulations
  sim.chars  <- apply(zh.chars,c(1,2),mean)  

#CIs of treatment effects  
  CI.up      <- coef.mat+1.96*(est.sd/sqrt(sims))
  CI.lo      <- coef.mat-1.96*(est.sd/sqrt(sims))
    
#Output list (FIX FOR MAT: add zh.PSmat.coef)
  out.list   <- list('est.sd'=est.sd,'est.ex'=est.ex,'zh.coef'=zh.coef,'zh.se'=zh.se,'zh.ci'=zh.ci,
                    'tot.bias'=tot.bias,'per.bias'=per.bias,'MSE.mat'=MSE.mat,'rat.se.sd'=rat.se.sd,
                    'se.mat'=se.mat,'coef.mat'=coef.mat,'cov.mat'=cov.mat,'provs'=provs,
                    'zh.afk'=zh.afk,'funs'=funs,"b.yi"=b.yi,'sim.chars'=sim.chars,
                    'CI.up'=CI.up,'CI.lo'=CI.lo,'coefs.true'=coefs.true)
                   
#####
## LINE GRAPHS
#####
  
#    if(l.plots==T){
#    
#      cat("\n","Making line graphs...","\n")
# 
#      zh.plots  <- zh_plots(tot.bias,MSE.mat,cov.mat,rat.se.sd,sim.scens,sim.chars,n.max,
#                            zh.afk,save=sav,funs=funs,plot.types=plot.t,cors=cors,colour=T)
#      
#      out.list  <- c(out.list,list(zh.plots=zh.plots))
#      }
  
###

  cat("\n","Simulation Finished.","\n")

  out.list
}