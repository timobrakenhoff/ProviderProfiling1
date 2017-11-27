#Function necessary for uniroot in fun_nvtdata_boot#.R

min.yint <- function(lin.pred,y.i,x){
  mean(as.vector(1/(1+exp(-(lin.pred+x)))))-y.i
}