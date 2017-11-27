#Function to center and recode covariates
#Continuous variables are centered
#binary variables are recoded to +1 and -1

recpred <- function(data, vars, type=c("cent","cod"),scal=F){

if (type=="cod"){ 
for (i in vars){  
  if(length(table(data[,i]))==2){
    data[data[,i]==0,i] <- -1
    } else {data[,i] <- as.numeric(scale(data[,i],scale=scal))}}
} else {if(type=="cent"){
  for (i in vars){  
    data[,i] <- scale(data[,i],scale=scal)}}
}

data
}
