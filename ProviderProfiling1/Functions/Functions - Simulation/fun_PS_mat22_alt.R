#PS MATCHING
# Extra arguments: bin.inc -- bin increments for matching PSs on x and y
#                  plot.mat -- Whether matching plot should be produced and saved.
#

PS_mat <- function(dat,remove=NA,seed=F,ref="B",fit.ml=NULL,n.zh=2,plot.mat=F,...){

#set seed  
if (seed==T){set.seed(123)}  

#Remove variables indicated in function
if (is.na(remove)==F){dat <- dat[,-remove]}

#Define indicator variable
ind   <- colnames(dat)[1]

#Define outcome variable
out   <- tail(colnames(dat),1)  

#Make fit.ml into dataframe to ensure rownumbers
fit.ml <- as.data.frame(fit.ml)

## Determine caliper width

#take log of fit.ml
log.fit  <- log(fit.ml)

#element of pool.sd
sd.tes   <- nrow(log.fit)-1

#pooled sd of PS A and B
pool.sd  <- sqrt((sd.tes*sd(log.fit[,1])+ sd.tes*sd(log.fit[,2]))/(nrow(log.fit)*3-1))

#caliper rounded to 3 decs
calp     <- round(0.2*pool.sd,3)

#Seperate PS into 3 groups according to values in ind
x1 <- fit.ml[dat$zh_afk==sort(unique(dat[,ind]))[1],1:2] 
x2 <- fit.ml[dat$zh_afk==sort(unique(dat[,ind]))[2],1:2] 
x3 <- fit.ml[dat$zh_afk==sort(unique(dat[,ind]))[3],1:2] 

#Round of PSs and make df
rx1 <- rounder(x1,calp)
rx2 <- rounder(x2,calp)
rx3 <- rounder(x3,calp)

#Make columns into factors
f.rx1 <- lapply(rx1,function(x) factor(x,levels=seq(0,1,calp)))
f.rx2 <- lapply(rx2,function(x) factor(x,levels=seq(0,1,calp)))
f.rx3 <- lapply(rx3,function(x) factor(x,levels=seq(0,1,calp)))

#Make contigency table of counts and make into vectors
t.rx1 <- as.vector(table(f.rx1))
t.rx2 <- as.vector(table(f.rx2))
t.rx3 <- as.vector(table(f.rx3))

#Bind vectors together into matrix
a.rx      <- cbind(t.rx1,t.rx2,t.rx3)

#Find minimum value per group for each bin (row)
min.a.rx  <- apply(a.rx,1,min)

#Find complement per group (amount of non matches per bin)
C1 	<- t.rx1 - min.a.rx		
C2 	<- t.rx2 - min.a.rx		
C3 	<- t.rx3 - min.a.rx		

#Rbind match number and non match number for each bin for each group
counts.1 <- as.vector(rbind(min.a.rx,C1))
counts.2 <- as.vector(rbind(min.a.rx,C2))
counts.3 <- as.vector(rbind(min.a.rx,C3))

#Make vector with indexes to choose
S1 <- rep(rep(c(1,0),length(min.a.rx)),times=counts.1)
S2 <- rep(rep(c(1,0),length(min.a.rx)),times=counts.2)
S3 <- rep(rep(c(1,0),length(min.a.rx)),times=counts.3)

#Properly arranged variables with rownames
arr.rx1 <- rx1[with(rx1, order(rx1[,2], rx1[,1])), ]
arr.rx2 <- rx2[with(rx2, order(rx2[,2], rx2[,1])), ]
arr.rx3 <- rx3[with(rx3, order(rx3[,2], rx3[,1])), ]

#Selected rows of arr.rx1
sel.rx1 <- arr.rx1[as.logical(S1),]
sel.rx2 <- arr.rx2[as.logical(S2),]
sel.rx3 <- arr.rx3[as.logical(S3),]

#Indices of rows that should be selected
ind.rx1 <- as.numeric(rownames(sel.rx1))
ind.rx2 <- as.numeric(rownames(sel.rx2))
ind.rx3 <- as.numeric(rownames(sel.rx3))

#Plot if plot.mat==T
if(plot.mat==T){

  #Max limits for plot
  xy.lim <- round_any(max(fit.ml[,1:2]), 0.01, f = ceiling)
  
  #Save plot if plot.mat==T
  png(filename="match_plot.png")
    
  #Make plot with axis limits = to max PS score in both groups
  plot(x1[,1:2],col="red",pch=16,xlim=c(0,xy.lim),ylim=c(0,xy.lim))
  points(x2[,1:2],col="blue",pch=16)
  points(x3[,1:2],col="green",pch=16)
  
  segments(fit.ml[ind.rx1,1],fit.ml[ind.rx1,2],fit.ml[ind.rx2,1],fit.ml[ind.rx2,2],col="red")
  segments(fit.ml[ind.rx1,1],fit.ml[ind.rx1,2],fit.ml[ind.rx3,1],fit.ml[ind.rx3,2],col="green")
  segments(fit.ml[ind.rx2,1],fit.ml[ind.rx2,2],fit.ml[ind.rx3,1],fit.ml[ind.rx3,2],col="blue")
  abline(v=(seq(calp/2,xy.lim,calp)), col="lightgray", lty="dotted")
  abline(h=(seq(calp/2,xy.lim,calp)), col="lightgray", lty="dotted")
  
  dev.off()
}

#### Analysis ####

#Final Dataset with matched individuals. Take indexes from original dat
fin.df   <- dat[c(ind.rx1,ind.rx2,ind.rx3),]

#Store calp and matched sets
st.char  <- c(calp,nrow(fin.df)/3)

#Make formula for regression
form3    <- as.formula(paste0(out,"~1+ relevel(",ind,",\"",ref,"\")"))

des      <- svydesign(ids=~1,weights=rep(1,nrow(fin.df)),data=fin.df)
GLM      <- svyglm(formula=form3, design=des,data=fin.df,family=binomial())

#Brier score
BS       <- sum((GLM$fitted.values-fin.df[,out])^2)/nrow(fin.df)

#Extract coefficients
LR.coef  <- GLM$coefficients[1:n.zh+1]

#Extract SEs
LR.se    <- summary(GLM)$coefficients[1:n.zh+1,2]

#Extract CIs
LR.ci    <- confint(GLM,parm=1:n.zh+1)

list('GLM'=GLM,"BS"=BS,'fit.ml'=fit.ml,"LR.coef"=LR.coef,"LR.se"=LR.se,"LR.ci"=LR.ci,
     "st.char"=st.char,"fin.df"=fin.df)
}