## FUNCTION TO MAKE GRAPHS USING OUTPUT FROM func_zh_sampling.R
## INPUT is object from previous function
# Fixed function to combine similar plots in one total figure. 



zh_plots <- function(tot.bias,MSE.mat,cov.mat,rat.se.sd,sim.scens,sim.chars,
                     n.max,zh.afk,save=T,funs,plot.types=1:4,cors,colour=T){

  #Reference scenario
  bas       <- c(1,sim.scens$by0[1],sim.scens$b0s[1])
  
  #All plot types
  plot.all  <- list("Bias"=tot.bias, "MSE"=MSE.mat, "Coverage"=cov.mat, "SE/SD"=rat.se.sd)
  
  #Which plots to include
  plot.list <- plot.all[plot.types]
  
  #To catch plots
  p.list    <- vector("list", length(plot.list)*ncol(sim.scens)/2)
  
  #Determine min and max for ylim on plot
  lims.bias <- c(rounder(min(tot.bias),0.02,fun="floor"),rounder(max(tot.bias),0.02,fun="ceiling"))
  lims.cov  <- c(rounder(min(cov.mat),0.01,fun="floor"),rounder(max(cov.mat),0.01,fun="ceiling"))

  lims.mse  <- c(0,rounder(max(MSE.mat),0.05,fun="ceiling"))
  lims.sesd <- c(rounder(min(rat.se.sd),0.1,fun="floor"),rounder(max(rat.se.sd),0.1,fun="ceiling"))
  
  #Determine tick marks using factors (NOT REALLY USED AT THE MOMENT)
  # fac.bias <- FACTOR(diff(lims.bias)*100)$pos
  # len.bias <- fac.bias[which.min(abs(fac.bias - 8))]+1
  # 
  # fac.cov  <- FACTOR(diff(lims.cov*100))$pos
  # len.cov  <- fac.cov[which.min(abs(fac.cov - 8))]+1
  # 
  # fac.mse  <- FACTOR(diff(lims.mse)*100)$pos
  # len.mse  <- fac.mse[which.min(abs(fac.mse - 8))]+1
  # 
  # fac.sesd <- FACTOR(diff(lims.sesd)*10)$pos
  # len.sesd <- fac.sesd[which.min(abs(fac.sesd - 8))]+1
  # 
  # #Determine tickmarks
  # ticks.bias <- seq(lims.bias[1],lims.bias[2],length.out=len.bias)
  # ticks.cov  <- seq(lims.cov[1],lims.cov[2],length.out=len.cov)
  # ticks.mse  <- seq(lims.mse[1],lims.mse[2],length.out=len.mse)
  # ticks.sesd <- seq(lims.sesd[1],lims.sesd[2],length.out=len.sesd)
  
  #X lab simulation matrix
  sim.lab   <- cbind(sim.scens[,1]*n.max,sim.chars[,ncol(sim.chars)]*n.max,sim.chars[,"p.zhB"]*n.max)
  colnames(sim.lab) <- c("Total sample size","Total amount of events","Provider volumes")
  
  #Determine ticks for x axis
  ticks.x <- list(zhperc=seq(0,10000,2000),
                  by0=seq(rounder(min(sim.lab[,2]),500,"floor"),rounder(max(sim.lab[,2]),500,"ceiling"),500),
                  b0s=seq(rounder(min(sim.lab[,3]),2000,"floor"),rounder(max(sim.lab[,3]),2000,"ceiling"),2000))
  
  
  # METHOD NAMES (MANUAL FOR PROPER LABELLING IN GGPLOT)
  #nam.meths <- gsub("[.]", " ", names(funs))
  nam.meths <- c("LR", "PS C", "PS W", "PS WT", "PS M") # MANUAL NAME GIVING...
  nam.list  <- list(expression("LR"),
                    expression("PS"[C]),
                    expression("PS"[W]),
                    expression("PS"[WT]),
                    expression("PS"[M]))
  
  #Axis titles for plots
  nam.tits <- c("total sample size","total amount of events","Provider volumes")
  
  #Functions for facet grid naming (NO LONGER WORKS IN NEW GGPLOT, FIXED BELOW WITH RENAMING ZH VARIABLE)
  # make_label <- function(value) {
  #   x <- as.character(value[1])
  #   y <- as.character(value[2])
  #   
  #   c(expression(bquote(B[.(x)])),expression(bquote(B[.(y)])))
  #   }
  # 
  # plot_labeller <- function(value) {
  #   do.call(expression, lapply(levels(value), make_label))
  #   }
  
  #Default GGPLOT COLOURS
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  cols.lines <- gg_color_hue(length(nam.meths))
  
################
## ALL PLOTS  ##
################

#LOOP OVER SCENARIO TYPES
for(l in seq_along(colnames(sim.scens))){

  #Location of reference in scenario matrix (sim.scens)
  ref.bas   <- which(colSums(t(sim.scens) == bas) == ncol(sim.scens))
  
  #Search for other scenarios excluding ref.bas
  ref.ot    <- match(unique(sim.scens[,l])[-which(unique(sim.scens[,l])==bas[l])],sim.scens[,l])
  
  #Bind all relevant scenarios
  colss     <- sort(c(ref.bas,ref.ot))
  
  #Make skeleton of df needed for plot using expand grid
  sim.grid  <- expand.grid(list(obs=sim.lab[,l][colss],zh=zh.afk,Method=nam.meths))
  
  #Fix zh labels for ggplotting
  sim.grid$zh <- paste0("B[",sim.grid$zh,"]")
  
  #Abline for ggplot
  #g.ab      <- geom_vline(xintercept = as.numeric(sim.lab[ref.bas,l]),linetype = "longdash")
  x.lab     <- colnames(sim.lab)[l]

  #####
  ### MAKE MANUAL GRAPHS FOR THE COMBINATIONS OF PERFORMANCE MEASURES
  #####
  
  #Make dataframes for each list element and cbind to grid
  plot.ress <- lapply(plot.list,function(x) cbind(sim.grid,as.vector(x[colss,])))
  
  #Now rbind all elements together
  plot.tots <- cbind(do.call("rbind",plot.ress),rep(names(plot.list),each=nrow(sim.grid)))
  
  #Fix colnames of matrix
  colnames(plot.tots) <- c("obs","zh","Method","out","mes")
  
  #####
  #Data for plotting bias + coverage AND mse + se/sd
  #####
  
  plot.b.c <- plot.tots %>% filter(.,mes=="Bias" | mes=="Coverage")
  plot.m.s <- plot.tots %>% filter(.,mes=="MSE" | mes=="SE/SD")
  
  #Workaround to create stable limits for plots by adding datapoints to ggplot 
  #outside of x margins. Add new datapoints with max and min for cov and bias
  plot.b.c.bias <- plot.b.c[which(plot.b.c$mes=="Bias")[1:2],]
  plot.b.c.cov <- plot.b.c[which(plot.b.c$mes=="Coverage")[1:2],]
  
  plot.m.s.mse <- plot.m.s[which(plot.m.s$mes=="MSE")[1:2],]
  plot.m.s.sesd <- plot.m.s[which(plot.m.s$mes=="SE/SD")[1:2],]
  
  
  plot.b.c.bias$obs <- c(-1000,-1000)
  plot.b.c.bias$out <- lims.bias
  plot.b.c.cov$obs <- c(-1000,-1000)
  plot.b.c.cov$out <- lims.cov
  
  plot.m.s.mse$obs <- c(-1000,-1000)
  plot.m.s.mse$out <- lims.mse
  plot.m.s.sesd$obs <- c(-1000,-1000)
  plot.m.s.sesd$out <- lims.sesd
  
  
  #Total plotting dataset
  plot.eve.b.c <- rbind(plot.b.c,plot.b.c.bias,plot.b.c.cov)
  plot.eve.m.s <- rbind(plot.m.s,plot.m.s.mse,plot.m.s.sesd)
  
  #Make condition for b0s where there is a double X axis
  if(colnames(sim.scens)[l]=="b0s"){
    
    lab.c  <- c(paste("     B=",ticks.x[[l]][2]),ticks.x[[l]][-c(1,2)])
    lab.ab <- c(paste("A=C=",(n.max-ticks.x[[l]][2])/2),(n.max-ticks.x[[l]][-c(1,2)])/2)
    
    #Labels for x axis
    labs.p <- paste(lab.c,lab.ab,sep="\n")
    
    #Plotting of bias + coverage
    sim.plot.bc <- ggplot(data=plot.eve.b.c, aes(x=obs, y=out,colour=Method)) +
                    geom_line(size=1.3) +
                    geom_point(size=4,aes(shape=Method)) + 
                    theme_bw()  +
                    facet_grid(mes~zh,scales="free_y", labeller=label_parsed) +
                    xlab(x.lab) + 
                    ylab("") +
                    scale_shape_manual(values = c(15:17,12,4),labels = nam.list) +
                    scale_colour_manual(values = cols.lines,labels = nam.list) +
                    scale_x_continuous(limits=ticks.x[[l]][c(1,length(ticks.x[[l]]))],
                                       breaks=ticks.x[[l]][-1],labels=labs.p)+
                    theme(text=element_text(size=24,family="Times"))  
        
    #plotting of MSE and SESD
    sim.plot.ms <- ggplot(data=plot.eve.m.s, aes(x=obs, y=out,colour=Method)) +
                    geom_line(size=1.3) +
                    geom_point(size=4,aes(shape=Method)) + 
                    theme_bw()  +
                    facet_grid(mes~zh,scales="free_y", labeller=label_parsed) +
                    xlab(x.lab) + 
                    ylab("") +
                    scale_shape_manual(values = c(15:17,12,4),labels = nam.list) +
                    scale_colour_manual(values = cols.lines,labels = nam.list) +
                    scale_x_continuous(limits=ticks.x[[l]][c(1,length(ticks.x[[l]]))],
                                       breaks=ticks.x[[l]][-1],labels=labs.p)+
                    theme(text=element_text(size=24,family="Times"))
    
  } else {
  
    #Plotting of bias + coverage
    sim.plot.bc <- ggplot(data=plot.eve.b.c, aes(x=obs, y=out,colour=Method)) +
                    geom_line(size=1.3) +
                    geom_point(size=4,aes(shape=Method)) + 
                    theme_bw()  +
                    facet_grid(mes~zh,scales="free_y", labeller=label_parsed) +
                    xlab(x.lab) + 
                    ylab("") +
                    scale_shape_manual(values = c(15:17,12,4),labels = nam.list) +
                    scale_colour_manual(values = cols.lines,labels = nam.list) +
                    scale_x_continuous(limits=ticks.x[[l]][c(1,length(ticks.x[[l]]))],
                                       breaks=ticks.x[[l]])+
                    theme(text=element_text(size=24,family="Times"))  
              
    #plotting of MSE and SESD
    sim.plot.ms <- ggplot(data=plot.eve.m.s, aes(x=obs, y=out,colour=Method)) +
                    geom_line(size=1.3) +
                    geom_point(size=4,aes(shape=Method)) + 
                    theme_bw()  +
                    facet_grid(mes~zh,scales="free_y", labeller=label_parsed) +
                    xlab(x.lab) + 
                    ylab("") +
                    scale_shape_manual(values = c(15:17,12,4),labels = nam.list) +
                    scale_colour_manual(values = cols.lines,labels = nam.list) +
                    scale_x_continuous(limits=ticks.x[[l]][c(1,length(ticks.x[[l]]))],
                                       breaks=ticks.x[[l]])+
                    theme(text=element_text(size=24,family="Times"))
  
    }
  
    
  
  #Name of saved pdf
  nam.plots  <- c(paste0("BC_",colnames(sim.scens)[l],"_",cors,".pdf"),
                  paste0("MS_",colnames(sim.scens)[l],"_",cors,".pdf"))
  
  
  #Save plot
  if(save==T){
    ggsave(sim.plot.bc,filename=nam.plots[1],width=14.3,height=8)
    ggsave(sim.plot.ms,filename=nam.plots[2],width=14.3,height=8)
  }
  
  #Save plots to out.list and give them a name
  p.list[[(l*2)-1]]      <- sim.plot.bc
  p.list[[(l*2)]]      <- sim.plot.ms
  names(p.list)[(l*2-1):(l*2)] <- nam.plots
  
}

p.list
}

