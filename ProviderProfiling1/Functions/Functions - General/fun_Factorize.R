# FUNCTION that factorizes multiple numeric variables in a data frame
# with certain levels and labels
# Notee that NAs have already been removed

factorize <- function(dat,colIndex,levs=c(1,2),labs=c("Ja","Nee")){
  
  sums <- NULL
  for (i in 1:length(colIndex)){
    dat[,colIndex[i]] <- factor(dat[,colIndex[i]],levels=levs,labels=labs)
    sums              <- rbind(sums,summary(dat[,colIndex[i]]))
  }
  list(DATA=dat, SUMS=sums)
}

