###################################################
# Logistic regression Function with fixed effects #
###################################################
##INPUT:
#dat     frame with last col  = outcome, first col = indicator, All predictors (centered)
#REMOVE   Variable #s that need to be removed before analysis
#SEED     Select to set seed 

##OUTPUT:
#GLM            GLM object with binomial and logit    
#BRIER SCORE    Calculated Brier Score

LRana_svy <- function(dat,remove=NA,seed=F,ref="B",n.zh=2,...){
  
  #set seed  
  if (seed==T){set.seed(123)}  
  
  #Remove variables indicated in function
  if (is.na(remove)==F){dat <- dat[,-remove]} else {dat <- dat}  
  
  #Define outcome variable name
  out   <- tail(colnames(dat),1)
  
  #Preds for model with ind
  preds <- colnames(dat)[-ncol(dat)]
  
  #Make formula for GLM
  form  <- as.formula(paste0(out,"~1+ relevel(",preds[1],",\"",
                            ref,"\") + ",paste(preds[-1],collapse="+")))
  
  #Make design for svyglm
  des   <- svydesign(ids=~1,weights=rep(1,nrow(dat)),data=dat)
  
  #Fit risk adjustment model  
  ml    <- svyglm(formula=form, design=des,data=dat,family=binomial())
  
  #Extract BS  
  BS    <- sum((ml$fitted.values-dat[,ncol(dat)])^2)/nrow(dat)
    
  #Extract coefficients
  LR.coef <- ml$coefficients[1:n.zh+1]
  
  #Extract SEs
  LR.se   <- summary(ml)$coefficients[1:n.zh+1,2]
  
  #Extract CIs
  LR.ci   <- confint(ml,parm=1:n.zh+1)
  
  #List GLM model and Brier Score
  list('GLM'=ml,"BS"=BS,"LR.coef"=LR.coef,"LR.se"=LR.se,"LR.ci"=LR.ci)
} 


