#########################
## SAMPLING SIMULATION ##
#########################

zh_sampling_new <- function(n.zh,n.max,b.coefs,funs,sim.scens,sims=10,seed=123,trim=0.02,
                            l.plots=F,ref="B",sav=T,plot.t=1:4,zh.tot.coef,confs=10,cors=0){
                            #,caliper=NULL,org.pts=5000){

#Set seed  
  if(is.null(seed)==F){set.seed(seed)}
  
#zh names
  zh.afk  <- LETTERS[LETTERS!=ref][1:(n.zh-1)]

#Initialize coefficient arrays in a list  
  zh.coef <-  rep(list(array(NA,dim=c(n.zh-1,nrow(sim.scens),sims),dimnames=list(zh.afk,1:nrow(sim.scens),1:sims))),
                  length(funs))

#Initialize SE arrays in a list    
  zh.se   <-  rep(list(array(NA,dim=c(n.zh-1,nrow(sim.scens),sims),dimnames=list(zh.afk,1:nrow(sim.scens),1:sims))),
                  length(funs))

#Initialize CI arrays in a list    
  zh.ci   <-  rep(list(array(NA,dim=c((n.zh-1)*2,nrow(sim.scens),sims),dimnames=list(paste0(rep(zh.afk,each=2),c(".L",".U")),
                                                                       1:nrow(sim.scens),1:sims))),length(funs))
#Make list for chars of the samples
  nams.cols <- c(paste0("p.zh",LETTERS[1:n.zh]),paste0("py1.zh",LETTERS[1:n.zh]),"py1.tot")
  zh.chars  <-  array(NA,dim=c(nrow(sim.scens),n.zh*2+1,sims),dimnames=list(1:nrow(sim.scens),
                      nams.cols,1:sims))
#Weight matrices
  w.tab     <- array(NA,dim=c(nrow(sim.scens),9*2,sims),
                   dimnames=list(1:nrow(sim.scens),
                   rep(c("Min","Q1","Med","Mean","Q3","Max","trim.A","trim.B","trim.C"),2),1:sims)) #Untrimmed and trimmed are rbinded
#Matching char matrix
  mat.char  <- array(NA,dim=c(nrow(sim.scens),2,sims),
                   dimnames=list(1:nrow(sim.scens),c("calp","mat.set"),1:sims)) #Untrimmed and trimmed are rbinded
  
#####
## START LOOP OVER SCENARIOS
#####

for (i in 1:nrow(sim.scens)){

#Print scenario number  
  cat("\n","Scenario",i,"of",nrow(sim.scens),"\n")  

#Loop over sims  
  for (j in 1:sims){

#Percentage counter (progress bar)
  if(j%%(sims/10)==0){cat(j/sims*10)}


#####
## Sampling procedure per hospital (output is list)
#####

#Update simulation parameters in b.coefs
  b.coefs$b.y0 <- sim.scens[i,2]
  b.coefs$b0A  <- sim.scens[i,3]
  b.coefs$b0C  <- sim.scens[i,3]

#Store simulated data
  samp.sim     <- sim.data(b.coefs=b.coefs,n.max=n.max,zh.perc=sim.scens[i,1],ref=ref,cors=cors,confs=confs)

#Run analyses on samp.df
  zh.samp      <- lapply(funs,function(f) f(dat=samp.sim$df.reg,fit.ml=samp.sim$gen.ps,
                                            trim=trim,ref=ref,n.zh=n.zh-1))
#Collect coefficients   
  zh.samp.coef <- lapply(zh.samp, function(x) x$LR.coef)  

#Collect SEs
  zh.samp.SE   <- lapply(zh.samp, function(x) x$LR.se)

#Collect CIs
  zh.samp.CI   <- lapply(zh.samp, function(x) x$LR.ci)

#Store weight summaries for PS.wei methods
  w.tab[i,,j]  <- c(zh.samp$PS.wei$sum.w,zh.samp$PS.wei.trim$sum.w)

#Store match chars
  mat.char[i,,j] <- zh.samp$PS.mat$st.char
  
#Store chars of simulation
  zh.chars[i,,j] <- samp.sim$sim.chars

#Assign to elements of the list of coef and SE  
  for (k in 1:length(funs)){
    zh.coef[[k]][,i,j] <- zh.samp.coef[[k]]
    zh.se[[k]][,i,j]   <- zh.samp.SE[[k]]
    zh.ci[[k]][,i,j]   <- t(zh.samp.CI[[k]])
  }   
  }
}


  
#####
## COLLECTION Results
#####
  cat("\n","Collecting simulation results...","\n")

#Initialize collection matrices
  se.mat  <- coef.mat <- est.sd <- rat.se.sd <- est.ex <- tot.bias <- per.bias <- MSE.mat <- cov.mat <-
  matrix(NA,ncol=length(funs)*2,nrow=nrow(sim.scens),dimnames=list(1:nrow(sim.scens),rep(zh.afk,length(funs))))
  
#Loop over list elements to collect results  
  for (jj in 1:length(funs)){
  
    #Counter for filling of matrices
    m.c <- seq(1,(length(funs)*2),2)[jj]
    
    #Mean SE, coef, and sd of coefficients
    se.mat[,m.c:(m.c+1)]      <- apply(zh.se[[jj]],c(2,1),mean)
    coef.mat[,m.c:(m.c+1)]    <- apply(zh.coef[[jj]],c(2,1),mean)
    est.sd[,m.c:(m.c+1)]      <- apply(zh.coef[[jj]],c(2,1),sd)
    
    #Ratio of mean(SE) to SD of estimates
    rat.se.sd[,m.c:(m.c+1)]   <- sqrt(apply(zh.se[[jj]]^2,c(2,1),mean))/est.sd[,m.c:(m.c+1)]
    
    #Percentage of effect estimates extremer than 5 (methods are cbinded)
    est.ex[,m.c:(m.c+1)]      <- apply(zh.coef[[jj]],c(2,1),function(x) sum(x>5 | x< -5)/length(x))
  
    #Bias per scenario and zh
    tot.bias[,m.c:(m.c+1)]    <- t(apply(zh.coef[[jj]],c(1,2),mean)-zh.tot.coef[[jj]])
    
    #Percentage bias
    per.bias[,m.c:(m.c+1)]    <- t(100*t(tot.bias[,m.c:(m.c+1)])/zh.tot.coef[[jj]])
    
    #MSE matrix
    MSE.mat[,m.c:(m.c+1)]     <- apply(zh.se[[jj]]^2,c(2,1),mean)+tot.bias[,m.c:(m.c+1)]^2
    
    #Coverage matrix (cbind 2 vectors that calculate coverage for each coefficient over sims)
    cov.mat[,m.c:(m.c+1)]     <- cbind(apply(zh.ci[[jj]][1,,] <= zh.tot.coef[[jj]][1] & 
                                             zh.tot.coef[[jj]][1] <= zh.ci[[jj]][2,,],1,mean), 
                                       apply(zh.ci[[jj]][3,,] <= zh.tot.coef[[jj]][2] & 
                                             zh.tot.coef[[jj]][2] <= zh.ci[[jj]][4,,],1,mean))
    }           

#Mean sim.chars over simulations
  sim.chars  <- apply(zh.chars,c(1,2),mean)  

#Mean weight summaries
  w.tab.sum  <- apply(w.tab,c(1,2),mean)
  
#Mat char summaries
  m.char.sum <- apply(mat.char,c(1,2),mean)

#CIs of treatment effects  
  CI.up      <- coef.mat+1.96*(est.sd/sqrt(sims))
  CI.lo      <- coef.mat-1.96*(est.sd/sqrt(sims))
    
#Output list (FIX FOR MAT: add zh.PSmat.coef)
  out.list   <- list('est.sd'=est.sd,'est.ex'=est.ex,'zh.coef'=zh.coef,'zh.se'=zh.se,'zh.ci'=zh.ci,
                    'tot.bias'=tot.bias,'per.bias'=per.bias,'MSE.mat'=MSE.mat,'rat.se.sd'=rat.se.sd,
                    'se.mat'=se.mat,'coef.mat'=coef.mat,'cov.mat'=cov.mat,'n.zh'=n.zh,
                    'zh.afk'=zh.afk,'funs'=funs,'sim.scens'=sim.scens,'sim.chars'=sim.chars,
                    'w.tab.sum'=w.tab.sum,'m.char.sum'=m.char.sum,'CI.up'=CI.up,'CI.lo'=CI.lo)
                   
#####
## LINE GRAPHS
#####
  
   if(l.plots==T){
   
     cat("\n","Making line graphs...","\n")

     zh.plots  <- zh_plots(tot.bias,MSE.mat,cov.mat,rat.se.sd,sim.scens,sim.chars,n.max,
                           zh.afk,save=sav,funs=funs,plot.types=plot.t,cors=cors,colour=T)
     
     out.list  <- c(out.list,list(zh.plots=zh.plots))
     }
  
###

  cat("\n","Simulation Finished.","\n")

  out.list
}