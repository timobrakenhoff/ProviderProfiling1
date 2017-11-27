#Reference values for marginal methods for different # of providers

refvals_ex <- function(provs,n.perzh,ref="A",cors=0,confs=10,b.yi=0.10,sims){
    
  #make empty list to collect marginal values
  mean.sd.list <- vector("list", length(provs))
  
  names(mean.sd.list) <- paste0("prov",provs)
  
  for(i in 1:length(provs)){
    
  #Counter
  cat("\n",paste0("prov",provs[i]),"\n")
    
  #Set confounding to 1 by putting al b.coefs for exposure model to 0
  b.coefs.ex  <- rep(0,provs[i])
  
  #Generate true coefficients
  coefs.out   <- c(1,seq(-1,1,length.out=provs[i]-1)) 
  
  #Make empty matrix for coefficients
  coef.col.mar <- matrix(NA,ncol=provs[i]-1,nrow=sims)
  
    for(j in 1:sims){
    
    #counter
    cat(j)
      
    #Generate data
    dat.mar     <- sim.data.ex(n.perzh=n.perzh,seed=NULL,fit=F,ref=ref,cors=cors,confs=confs,provs=provs[i],
                           b.yi=b.yi,b.coefs.ex=b.coefs.ex,coefs.out=coefs.out)
    
    #Perform marginal analysis and save to matrix
    coef.col.mar[j,]    <- coefficients(glm(overleden_zh~relevel(zh_afk,ref),data=dat.mar$df.reg,family=binomial))[-1]
    
    }
  
  #Save mean,sd,se into a list entry and gives names to col and rows
  mean.sd.list[[i]] <- rbind(coefs.out[-1],apply(coef.col.mar,2,function(x)c(mean(x),sd(x),sd(x)/sqrt(sims))))
  colnames(mean.sd.list[[i]]) <- LETTERS[2:provs[i]]
  rownames(mean.sd.list[[i]]) <- c("TRUE","Mean","SD","SE")
  
  }
  
#Gather output

list(mean.sd.list=mean.sd.list)
}
