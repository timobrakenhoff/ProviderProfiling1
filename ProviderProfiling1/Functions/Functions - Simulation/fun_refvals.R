#Function to determine reference values for marginal and matching treatment effects for different event rates

refvals <- function(b.coefs,vals=c(2.29,1,2,4,5),nmax=1e4,ref="B",sims=5,cors=0){

#Make b.coefs.mat for matching    
b.coefs.mat <- b.coefs  

#Set confounding to 0 for marginal
b.coefs[,4:5] <- 0

#Collection matrices
col.mar    <- matrix(NA,nrow=sims,ncol=2*length(vals))
col.mat    <- col.mar


#Loop over all vals
for(j in 1:length(vals)){
  
  cat(paste("\n",j,"\n"))
  
  #Set byo to one of vals
  b.coefs$b.y0     <- vals[j]
  b.coefs.mat$b.y0 <- vals[j]
  
  #Simulation for each vals
  for(l in 1:sims){
    
    cat(l)
    
    #Sim data
    dat.mar  <- sim.data(b.coefs,nmax,seed=NULL,fit=F,ref=ref,confs=10,cors=cors)
    dat.mat  <- sim.data(b.coefs.mat,nmax,seed=NULL,fit=F,ref=ref,confs=10,cors=cors) # FIT SET TO F SUCH THAT THERE IS VARIABILITY IN PS
    
    #Perform analysis
    ana.mar  <- coefficients(glm(overleden_zh~relevel(zh_afk,ref),data=dat.mar$df.reg,family=binomial))[2:3]
    ana.mat  <- PS_mat(dat=dat.mat$df.reg,remove=NA,seed=F,ref=ref,fit.ml=dat.mat$gen.ps,n.zh=2,plot.mat=F)
    
    #Collect coefficients
    col.mar[l,(j*2-1):(j*2)]  <- ana.mar
    col.mat[l,(j*2-1):(j*2)]  <- ana.mat$LR.coef
    }
}

mar.coef <- apply(col.mar,2,mean) 
mat.coef <- apply(col.mat,2,mean)

mar.sd   <- apply(col.mar,2,sd) 
mat.sd   <- apply(col.mat,2,sd)  

mar.se   <- mar.sd/sqrt(sims)
mat.se   <- mat.sd/sqrt(sims)

dat.mar  <- round(rbind(rep(vals,each=2),mar.coef,mar.sd,mar.se),3)
dat.mat  <- round(rbind(rep(vals,each=2),mat.coef,mat.sd,mat.se),3)

colnames(dat.mar) <- colnames(dat.mat) <- rep(c("A","C"),length(vals))

list('Marginal'=dat.mar,'Matching'=dat.mat)
}