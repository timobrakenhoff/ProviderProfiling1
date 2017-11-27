#################### SIMULATION FUNCTION ###################
## Edited to allow for differing amount of providers ###
## Start editing on 25-08-2016 ##
## Borrowed elements from fun_nvtdata_boot16.R ##

sim.data.ex <- function(n.perzh,seed=NULL,fit=T,ref="B",cors,confs,provs=5,
                        b.yi,b.coefs.ex,coefs.out){
  
  #Set seed
  if(is.null(seed)==F){set.seed(seed)}
  
  #Calculate maximum amount of people necessary
  n.max <- n.perzh*provs
  
  #Add sample size to be able to draw equal sized groups for providers
  n.ext <- n.max*2
  
  #Make varcov for mvrnorm
  Sigma       <- matrix(rep(cors,confs*confs),nrow=confs)
  diag(Sigma) <- 1
  
  #Sample variables into z
  z <- mvrnorm(n = n.ext, rep(0,confs), Sigma)
  
  #compute the probabilities of each outcome for the multinomial model (ref=1st zh)
  #Put these values into a matrix with providers in columns
  
  #Make probs matrix for multinom probabilities
  vProb        <-  matrix(NA,nrow=nrow(z),ncol=provs)
  
  #Calculate vector of probabilities (not really) for only coefficients (ASSUMES AN INTERCEPT of 0)
  #No intercept included in the equation (thus assumes an intercept of 0)
  for(j in 1:provs){
    vProb[,j]  <-  exp(z%*%rep(b.coefs.ex[j]/confs,confs)) 
  }
  
  #Make multinomial selection
  mChoices     <- t(apply(vProb, 1, rmultinom, n = 1, size = 1))
  
  #Select first x.sams of each zh and store as rownumbers in 1 vector
  sel.choice <- NULL
  for(h in 1:provs){
    sel.choice <-c(sel.choice,which(mChoices[,h]==1)[1:n.perzh])}
  
  
  #Check maximum rownumbers (how many total people you required)
  max.sel.choice <- max(sel.choice)
  
  #Take the selection of covariates (per column of zh) 
  sel.samp   <- z[sel.choice,] 
  
  #Draw categories of X from mChoices with rownumbers as selected in sel.choice (goes per column)
  X         <- as.factor(apply(mChoices[sel.choice,], 1, function(x) which(x==1)))
  
  #Change labels of X (Change to letters (there wont be any more categories than 20))
  levels(X) <- LETTERS[1:length(levels(X))]
  
  #Fit multinomial model for fitted values
  if(fit==T){
    ml.samp      <- multinom(X~sel.samp,trace=F)
    
    #Fitted values are the propensity scores (order columns according to levels of X)
    gen.ps       <- as.data.frame(fitted.values(ml.samp)[,levels(X)])
    
    #Give numbered names to the gen PS (NO more categories than 20 so LETTERS ARE FINE)
    #names(gen.ps) <- paste("zh",names(gen.ps),sep="_")
    
    } else {gen.ps <- vProb}

#Simulate y variable  
  
  #Store dummy matrix and change colnames
  dums    <- dummy(X)
  colnames(dums) <- LETTERS[1:ncol(dums)]

  #Select relevant dummies (no need for a reference...)
  dums1   <- dums[,!colnames(dums)==ref]

  #Combine dummy variable values with other covariates
  comb.set <- cbind(dums1,sel.samp)
  
  #Coefficients vector for all covariates (dummies + confounders)
  x.zhs.d <- c(coefs.out[-1],rep(1/confs,confs)) #Remove effect of reference prov from coefs.out
  
  #Create linear predictor (No intercept/ or assumed to be 0)
  lin.pred.cov.mat <- comb.set %*% x.zhs.d
  
  #Run uniroot to find overall intercept (ALLOW INTERVAL TO BE BROADENED)
  y.int            <- uniroot(f=min.yint,interval=c(-5,2),lin.pred=lin.pred.cov.mat,
                              y.i=b.yi,extendInt = "yes")$root
  
  # Calculate predicted probability to have an event using y.i.bnt
  y.pr             <- as.vector(1/(1+exp(-(lin.pred.cov.mat+y.int)))) 

  #Use probabilities to draw from binomial distribution
  y                <- rbinom(length(y.pr),1,y.pr)
  
  #True Incidence of outcome per X (IN PERCENTAGE) theoretical prob in dataset (not observed)
  ind.zh           <- round(unlist(by(y.pr,X,mean,simplify=F))*100,digits=2)
  
  #Observed incidence of outcome
  obs.ind          <- tapply(y,X,mean)
  
  #Combine into final dataset
  df.reg <- data.frame(zh_afk=X,sel.samp,overleden_zh=y)

  #Check if coefficients are reproduced when running log reg  
  #out.mod <- glm(overleden_zh~zh_afk+.,data=df.reg,family="binomial")
  
  #Collect simulation chars
  sim.chars <- round(c(colMeans(gen.ps),ind.zh,mean(y.pr)),3)
  
  list(df.reg=df.reg,sim.chars=sim.chars,gen.ps=gen.ps)#,out.mod=out.mod)
}


