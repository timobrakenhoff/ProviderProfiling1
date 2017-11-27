###############################################################
# Propensity Score covariate adjustment -- Spreeuwenberg 2010 #
###############################################################
##INPUT:
#DATA     frame with last col  = outcome, first col = indicator, All predictors (centered)
#REMOVE   Variable #s that need to be removed before analysis
#SEED     Select to set Seed 
#PLOTS    Select true if you want PS overlap plots(REMOVED)
#PREC     Percentage of total sample size (necessary in function for plots (REMOVED);
#         For data is specified outside function)

##OUTPUT:
#GLM            GLM object with binomial and logit    
#BRIER SCORE    Calculated Brier Score

PS_cov_new <- function(dat,remove=NA,seed=F,scal=F,ref="A",fit.ml=NULL,n.zh=2,...) {

#set seed  
if (seed==T){set.seed(123)}    

#Remove variables indicated in function
if (is.na(remove)==F){dat <- dat[,-remove]}

#Define indicator variable
ind   <- colnames(dat)[1]

#Define outcome variable
out   <- tail(colnames(dat),1)  

#If scale, scale fitted values
if(scal==T){fit.ml <- scale(fit.ml)}

#Append fitted values to dataframe
dat <- cbind(dat,fit.ml)

###Effect Estimation After Correction
  #PS variables
  preds.ps <- levels(dat[,1])[-length(levels(dat[,1]))]
  
  #Formula
  form3    <- as.formula(paste0(out,"~ 1 + relevel(",ind,",\"",ref,"\") + ",
                                paste(preds.ps,collapse="+")))   
  
  #Post correction logistic model using svyglm()
  des      <- svydesign(ids=~1,weights=rep(1,nrow(dat)),data=dat)
  GLM      <- svyglm(formula=form3, design=des,data=dat,family=binomial())

###
# Brier score
BS       <- sum((GLM$fitted.values-dat[,out])^2)/nrow(dat)

#Extract coefficients
LR.coef <- GLM$coefficients[1:n.zh+1]

#Extract SEs
LR.se   <- summary(GLM)$coefficients[1:n.zh+1,2]

#Extract CIs
LR.ci   <- confint(GLM,parm=1:n.zh+1)

list('GLM'=GLM,"BS"=BS,'fit.ml'=fit.ml,"LR.coef"=LR.coef,"LR.se"=LR.se,"LR.ci"=LR.ci)
}