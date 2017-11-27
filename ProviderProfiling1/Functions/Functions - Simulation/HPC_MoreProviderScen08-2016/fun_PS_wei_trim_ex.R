##############################
# Propensity Score weighting #
##############################
##INPUT:
#DATA     frame with last col  = outcome, first col = indicator, All predictors (centered)
#REMOVE   Variable #s that need to be removed before analysis
#SEED     Select to set Seed 
#PLOTS    Select true if you want PS overlap plots 
#PREC     Percentage of total sample size (necessary in function for plots; (REMOVED)
#         For data is specified outside function) (REMOVED)

##OUTPUT:
#GLM            GLM object with binomial and logit    
#BRIER SCORE    Calculated Brier Score

PS_wei_trim_ex <- function(dat,remove=NA,seed=F,ref="A",fit.ml=NULL,n.zh=2,trim=0.02,...) {
  
  #set seed  
  if (seed==T){set.seed(123)}  
  
  #Remove variables indicated in function
  if (is.na(remove)==F){dat <- dat[,-remove]}
  
  #Define indicator variable
  ind   <- colnames(dat)[1]
  
  #Define outcome variable
  out   <- tail(colnames(dat),1)  

  #Variable with zh names
  nam.zh <- colnames(fit.ml)
  
  #Append column with unstabalized weights (1/ PS of received treatment) (See McCaffrey 2013)
  weig <- 0
  for(i in 1:ncol(fit.ml)){
    weig <- weig +1/fit.ml[,nam.zh[i]]*(dat[,1]==nam.zh[i])  
  }
  
  #Append weights to dataset
  dat2 <- dat
  dat2$w <- weig
  
  #Identity of trimmed weights
  id.w <- table(dat2[dat2$w>quantile(dat2$w,1-trim),][,1])
  
  #Trim weights (Mean is allowed to shift)
  dat2$w[dat2$w>quantile(dat2$w,1-trim)] <-quantile(dat2$w,1-trim)

  #Summary of weights distribution (including id.w)
  sum.w       <- c(summary(dat2$w),round(id.w/sum(id.w),3))
  
  #Make formula for regression
  form3       <- as.formula(paste0(out,"~1+ relevel(",ind,",\"",ref,"\")"))

  #Make weight design
  design.c    <- svydesign(ids=~1,weights=~w,data=dat2)
  
  #Run IPTW svyGLM
  GLM         <- svyglm(form3,design=design.c,family=quasibinomial())

  #Brier score
  BS       <- sum((GLM$fitted.values-dat2[,out])^2)/nrow(dat2)
  
  #Extract coefficients
  LR.coef <- GLM$coefficients[1:n.zh+1]
  
  #Extract SEs
  LR.se   <- summary(GLM)$coefficients[1:n.zh+1,2]
  
  #Extract CIs
  LR.ci   <- confint(GLM,parm=1:n.zh+1)
  
  list('GLM'=GLM,'BS'=BS,'sum.w'=sum.w,'fit.ml'=fit.ml,"LR.coef"=LR.coef,"LR.se"=LR.se,"LR.ci"=LR.ci)
}
