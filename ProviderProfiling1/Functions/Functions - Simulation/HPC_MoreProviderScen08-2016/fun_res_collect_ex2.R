#################### Extra Simulation Collection ###################
## Collect results of extra simulation

#packages
install.packages("tidyr")
#install.packages("rounder")
require(tidyr)
require(ggplot2)
require(Rmisc)
#require(rounder)

#Setwd
setwd("F:\\Google Drive\\PhD UMC\\MIR\\MIR Project 1 Risk Adjustment\\\\Final Images - R\\SIMULATION\\Extra_scen_16")
setwd("C:\\Users\\Timo\\Google Drive\\PhD UMC\\MIR\\MIR Project 1 Risk Adjustment\\Final Images - R\\SIMULATION\\Extra_scen_16")

#Matrices to include for plotting
nam.mats      <- c("tot.bias","cov.mat","rat.se.sd","MSE.mat")
plot.nam.mats <- c("Bias","Coverage","SE/SD", "MSE")

#List all files
fils <- list.files(pattern="HPC")

sim_plot_ex <- function(nam.mats,plot.nam.mats,fils,sav=F){

#Give manual names for the methods for plotting  
nam.meths <- c("LR", "PS C", "PS W", "PS WT")
nam.list  <- list(expression("LR"),
                expression("PS"[C]),
                expression("PS"[W]),
                expression("PS"[WT]))

#Default GGPLOT COLOURS
gg_color_hue <- function(n) {
hues = seq(15, 375, length=n+1)
hcl(h=hues, l=65, c=100)[1:n]
}

#Specify colour codes for fills per method
cols.lines <- gg_color_hue(length(nam.meths)+1)[1:length(nam.meths)]

#Extract name for plot naming    
fils.sim <- gsub("HPC_MIR1_","",fils)
fils.sim.c <- gsub(".Rdata","",fils.sim)

#Make list to save all plots in
p.list    <- vector("list", length(nam.mats)/2*length(fils))

#Empty list for collection of all sim results
sim.res.list <- vector("list", length(fils))

#Limits matrix
lim.mat <- matrix(NA,nrow=length(nam.mats),ncol=2)

#Load all data to find limits for the graphs
for(k in 1:length(fils)){

  #Load files
  load(fils[k])
  
  #Save file in list
  sim.res.list[[k]] <- SIM_RES_EXTRA
}

#Determine limits for each outcome measure
for(i in seq_along(nam.mats)){
  
  #Bind all results from outcome measure together
  all.tog <- do.call(cbind,lapply(sim.res.list,function(x)x[[nam.mats[i]]]))  
  
  #Determine limits by taking min and max
  lim.mat[i,] <- c(rounder(min(all.tog),0.02,fun="floor"),
                   rounder(max(all.tog),0.02,fun="ceiling"))
}


#Loop over the different output files
for(j in 1:length(fils)){

  #Extract outcome tables
  out_mat.b <- as.data.frame(sim.res.list[[j]][[nam.mats[1]]])
  out_mat.c <- as.data.frame(sim.res.list[[j]][[nam.mats[2]]])
  out_mat.s <- as.data.frame(sim.res.list[[j]][[nam.mats[3]]])
  out_mat.m <- as.data.frame(sim.res.list[[j]][[nam.mats[4]]])

  
  #Add method names as column in df
  out_mat.b$Method <- out_mat.c$Method <- out_mat.s$Method <- 
  out_mat.m$Method <- rownames(out_mat.b)
  
  #Add outcome measure type
  out_mat.b$mes <- plot.nam.mats[1]
  out_mat.c$mes <- plot.nam.mats[2]
  out_mat.s$mes <- plot.nam.mats[3]
  out_mat.m$mes <- plot.nam.mats[4]
  
  #Reshape df into long
  out_mat.b.long   <- gather_(out_mat.b, "Provider", "out",
                            colnames(out_mat.b)[-c(length(colnames(out_mat.b))-1,
                                                   length(colnames(out_mat.b)))])
  out_mat.c.long   <- gather_(out_mat.c, "Provider", "out",
                              colnames(out_mat.c)[-c(length(colnames(out_mat.c))-1,
                                                     length(colnames(out_mat.c)))])
  out_mat.s.long   <- gather_(out_mat.s, "Provider", "out",
                              colnames(out_mat.s)[-c(length(colnames(out_mat.s))-1,
                                                     length(colnames(out_mat.s)))])
  out_mat.m.long   <- gather_(out_mat.m, "Provider", "out",
                              colnames(out_mat.m)[-c(length(colnames(out_mat.m))-1,
                                                     length(colnames(out_mat.m)))])
  
  
  #COmbine bias + coverage and MSE + sesd
  out.mat.bc <- rbind(out_mat.b.long,out_mat.c.long)
  out.mat.ms <- rbind(out_mat.s.long,out_mat.m.long)
  
  #Add abline to facets
  ab.bc <- data.frame(mes = "Coverage", Z = 0.95)
  ab.ms <- data.frame(mes = "SE/SD", Z = 1)
  
  
  ## ADD DATAPOINTS FOR FIXED SCALES
  out.mat.bc.bias <- out.mat.bc[which(out.mat.bc$mes=="Bias")[1:2],]
  out.mat.bc.cov  <- out.mat.bc[which(out.mat.bc$mes=="Coverage")[1:2],]
  out.mat.ms.mse  <- out.mat.ms[which(out.mat.ms$mes=="MSE")[1:2],]
  out.mat.ms.sesd <- out.mat.ms[which(out.mat.ms$mes=="SE/SD")[1:2],]
  
  
  out.mat.bc.bias$Provider <- c("X","X")
  out.mat.bc.bias$out <- lim.mat[1,]
  out.mat.bc.cov$Provider <- c("X","X")
  out.mat.bc.cov$out  <- lim.mat[2,]
  
  out.mat.ms.mse$Provider <- c("X","X")
  out.mat.ms.mse$out  <- lim.mat[4,]
  out.mat.ms.sesd$Provider <- c("X","X")
  out.mat.ms.sesd$out <- c(lim.mat[3,1]-0.005,lim.mat[3,2]+0.005)
  
  out.eve.bc <- rbind(out.mat.bc,out.mat.bc.bias,out.mat.bc.cov)
  out.eve.ms <- rbind(out.mat.ms,out.mat.ms.mse,out.mat.ms.sesd)
  
  
 #PLot for bias and coverage with 2 plot types in one facet grid
  plot.bc <- ggplot() + geom_bar(data=subset(out.eve.bc,mes=="Bias"),aes_string(x="Provider", y="out", fill="Method"),
                      stat="identity", position=position_dodge())+
              geom_point(data=subset(out.eve.bc,mes=="Coverage"),
                        aes_string(x="Provider", y="out", colour="Method",shape="Method"),
                        position=position_dodge(width=0.15),size=4)+
              scale_fill_manual(values=cols.lines,labels=nam.list)+          
              guides(fill=FALSE)+          
              scale_colour_manual(values=cols.lines,labels=nam.list)+
              scale_shape_manual(values=c(15:17,12),labels=nam.list)+
              facet_grid(mes~.,scales="free_y")+
              geom_hline(data = ab.bc, aes(yintercept = Z),linetype=2)+
              scale_x_discrete(limits=unique(out.mat.bc$Provider))+
              ylab("")+
              theme_bw() +
              theme(text=element_text(size=24,family="Times"))
    
  
  #Plot for mse and sesd
  plot.ms <- ggplot() + geom_bar(data=subset(out.eve.ms,mes=="MSE"),aes_string(x="Provider", y="out", fill="Method"),
                                 stat="identity", position=position_dodge())+
                        geom_point(data=subset(out.eve.ms,mes=="SE/SD"),
                                    aes_string(x="Provider", y="out", colour="Method",shape="Method"),
                                    position=position_dodge(width=0.15),size=4)+
                        scale_fill_manual(values=cols.lines,labels=nam.list)+          
                        guides(fill=FALSE)+   
                        scale_colour_manual(values=cols.lines,labels=nam.list)+
                        scale_shape_manual(values=c(15:17,12),labels=nam.list)+
                        facet_grid(mes~.,scales="free_y")+
                        geom_hline(data = ab.ms, aes(yintercept = Z),linetype=2)+
                        scale_x_discrete(limits=unique(out.mat.ms$Provider))+
                        ylab("")+
                        theme_bw() +
                        theme(text=element_text(size=24,family="Times"))
                      
  
  #Make name of plot
  nam.plot.bc <- paste0(fils.sim.c[j],"_BC",".pdf")
  nam.plot.ms <- paste0(fils.sim.c[j],"_MS",".pdf")
  
  #Save plots
  if(sav==T){
    ggsave(plot.bc,filename=nam.plot.bc,width=14.3,height=8)
    ggsave(plot.ms,filename=nam.plot.ms,width=14.3,height=8)
  }
  
  
  #Save plots to out.list and give them a name
  p.list[[(j*2)-1]]      <- plot.bc
  p.list[[(j*2)]]        <- plot.ms
  names(p.list)[(j*2-1):(j*2)] <- c(nam.plot.bc,nam.plot.ms)
}

p.list
}

sim_plot_ex(nam.mats,plot.nam.mats,fils,sav=T)

