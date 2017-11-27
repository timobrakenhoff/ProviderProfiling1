#Function for rbind of tables of binary variables
#Date created: 25-03-2015


tabsum <- function(data,dich.vars){
  
  tabs <- NULL

for (i in dich.vars){
  ff   <- table(data[,i])
  if(ff[[1]]==nrow(data)){ff<-c(ff,0)}
  tabs <- rbind(tabs,ff)
}

rownames(tabs) <- colnames(data)[dich.vars]

tabs}