#################### SIMULATION FUNCTION ###################
sim.data <- function(b.coefs,n.max,zh.perc=NULL,seed=NULL,fit=T,ref="B",cors,confs){
  
  #Set seed
  if(is.null(seed)==F){set.seed(seed)}
  
  #Set n if type is observations (then the amount of observations is simply diminished)
  n  <- ifelse(is.null(zh.perc),n.max,n.max*zh.perc)
  
  #Make varcov for mvrnorm
  Sigma       <- matrix(rep(cors,confs*confs),nrow=confs)
  diag(Sigma) <- 1
  
  #Sample variables into z
  z <- mvrnorm(n = n, rep(0,confs), Sigma)
  
  #compute the probabilities of each outcome for the multinomial model (ref="C")
  pA <- exp(b.coefs$b0A+rowSums((b.coefs$b1A/confs)*z))/
        (1+exp(b.coefs$b0A+rowSums((b.coefs$b1A/confs)*z))+
           exp(b.coefs$b0C+rowSums((b.coefs$b1C/confs)*z)))
  pC <- exp(b.coefs$b0C+rowSums((b.coefs$b1C/confs)*z))/
        (1+exp(b.coefs$b0A+rowSums((b.coefs$b1A/confs)*z))+
           exp(b.coefs$b0C+rowSums((b.coefs$b1C/confs)*z)))
  pB <- 1-pA-pC

  #Sample out of three categories using probs calculated above
  X <- as.factor(apply(cbind(pA,pB,pC),1,function(x) sample(c("A","B","C"),1,replace=T,prob=x)))

  #Fit multinomial model for fitted values
  if(fit==T){
    ml.samp      <- multinom(relevel(X,ref=ref)~z,trace=F)
    
    #Fitted values are the propensity scores (order columns according to alphabet)
    gen.ps       <- fitted.values(ml.samp)[,LETTERS[1:3]]
    } else {gen.ps <- cbind(pA,pB,pC)}

#Simulate y variable  
  #set reference zh
  zh.ref  <- paste0("dummy(X)",ref)
  
  #Store dummy matrix and change colnames
  dums    <- dummy(X)
  colnames(dums) <- LETTERS[1:ncol(dums)]

  #Select relevant dummies
  dums1   <- dums[,!colnames(dums)==ref]

  #Put effects in the right order
  x.zhs   <- c(b.coefs$byXA,b.coefs$byXB,b.coefs$byXC)
  x.zhs.d <- as.vector(c(x.zhs[which(LETTERS %in% ref)],x.zhs[-which(LETTERS %in% ref)]))

#Make linear predictor for X (include !ref dummies and add an intercept)
  linpred <- cbind(b.coefs$b.y0,dums1) %*% x.zhs.d
  
  #Make linear predictor
  y.reg   <- linpred + rowSums((b.coefs$byz/confs)*z)
  
  #Make y
  y.pr    <- as.vector(1/(1+exp(-y.reg)))
  y       <- rbinom(n,1,y.pr)

  #Incidence per zh
  ind.zh <- unlist(by(y.pr,X,mean,simplify=F))
  
  #Combine into data frame
  df.reg <- data.frame(zh_afk=X,z,overleden_zh=y)
  
  #Collect simulation chars
  sim.chars <- round(c(colMeans(gen.ps),ind.zh,mean(y.pr)),3)
  
  list(df.reg=df.reg,sim.chars=sim.chars,gen.ps=gen.ps)
}

### TEST OF OUTCOME MODEL ###
#ll <- glm(overleden_zh~relevel(zh_afk,ref)+.-zh_afk,data=df.reg,family="binomial")
###
