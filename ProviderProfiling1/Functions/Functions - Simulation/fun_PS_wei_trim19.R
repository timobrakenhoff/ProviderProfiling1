##############################
# Propensity Score weighting #
##############################
##INPUT:
#DATA     frame with last col  = outcome, first col = indicator, All predictors (centered)
#REMOVE   Variable #s that need to be removed before analysis
#SEED     Select to set Seed 
#PLOTS    Select true if you want PS overlap plots 
#PREC     Percentage of total sample size (necessary in function for plots;
#         For data is specified outside function)

##OUTPUT:
#GLM            GLM object with binomial and logit    
#BRIER SCORE    Calculated Brier Score

PS_wei_new_trim <- function(dat,remove=NA,seed=F,trim=0.02,ref="A",fit.ml=NULL,n.zh=2,...) {
  
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
  dat$w <-  1/fit.ml[,nam.zh[1]]*(dat[,1]==nam.zh[1]) + 
            1/fit.ml[,nam.zh[2]]*(dat[,1]==nam.zh[2]) + 
            1/fit.ml[,nam.zh[3]]*(dat[,1]==nam.zh[3])
  
  #Identity of trimmed weights
  id.w <- table(dat[dat$w>quantile(dat$w,1-trim),][,1])
  
  #Trim weights (Mean is allowed to shift)
  dat$w[dat$w>quantile(dat$w,1-trim)] <-quantile(dat$w,1-trim)
  
  #Summary of weights distribution (including id.w)
  sum.w       <- c(summary(dat$w),round(id.w/sum(id.w),3))
  
  #Make formula for regression
  form3    <- as.formula(paste0(out,"~1+ relevel(",ind,",\"",ref,"\")"))
  
  #Make weight design
  design.c <- svydesign(ids=~1,weights=~w,data=dat)
  
  #Run IPTW svyGLM
  GLM      <- svyglm(form3,design=design.c,family=quasibinomial())
  
  #Brier score
  BS       <- sum((GLM$fitted.values-dat[,out])^2)/nrow(dat)
  
  #Extract coefficients
  LR.coef  <- GLM$coefficients[1:n.zh+1]
  
  #Extract SEs
  LR.se    <- summary(GLM)$coefficients[1:n.zh+1,2]
  
  #Extract CIs
  LR.ci    <- confint(GLM,parm=1:n.zh+1)
  
  list('GLM'=GLM,'BS'=BS,'sum.w'=sum.w,'fit.ml'=fit.ml,"LR.coef"=LR.coef,"LR.se"=LR.se,"LR.ci"=LR.ci)
}
