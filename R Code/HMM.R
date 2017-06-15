######################################################
#
# File containing the function definition pertaining to
# the markov chain analysis. Contains data loading,
# Manipulation, plotting and CTMC fitting
#
#######################################################



CreateMarkovChainList <- function(TimeSeries){

  # create the list of states along with the time spent to fit a Markov chain from the time series
  
  datatest <- as.data.table(TimeSeries)
  setnames(datatest,"V1")
  datatest[, count := .I]
  datatest[, mincount := min(count)-1,.(rleid(V1))]
  datatest <- unique(datatest[,.(V1,mincount)])
  
  listtest <- list(datatest$V1,datatest$mincount)
  
  return(listtest)
}

CreateTimeSeries <- function(mcdf){

  # recreate the time series from the list obtained by simulating the markov chain

  mcdf2<-as.data.table(mcdf)
  mcdf2[,time := floor(as.numeric(as.character(time)))]
  mcdf2[,states := as.numeric(as.character(states))]
  mcdf2$time2 <- mcdf2$time[-1] - mcdf2$time[-length(mcdf2$time)] 
  mcdf3 <- rep(mcdf2$states,mcdf2$time2)
  return(mcdf3)
}

GetSummarydf <-function(listtest){

  # generate a data frame containing a number of useful summary statistics of the markov chain/time series

  df.test <- as.data.frame(listtest)
  setnames(df.test,c("state","time"))
  df.test <- as.data.table(cbind(df.test[-nrow(df.test),],df.test[-1,]))
  setnames(df.test,c("Statebegin","time1","Stateend","time2"))
  df.test[,TimeSpent := time2-time1]
  df.test[,MeanTimeSpent := mean(TimeSpent),.(Statebegin)]
  df.test[,TimeVisited := .N,(Statebegin)]
  
  stat <- list()
  
  stat$NumberOfTransitions <- nrow(df.test)
  
  df.test <- unique(df.test[,.(TimeVisited,MeanTimeSpent,Statebegin)])
  stat$summary <- df.test
  return(stat)
}

TestMarkov <- function(listtest){

 # Generate a data frame containing a number of useful statistics to test whether or not the time series is actually Markovian.

  df.test <- as.data.frame(listtest)
  setnames(df.test,c("state","time"))
  df.test <- as.data.table(cbind(df.test[-c(nrow(df.test),nrow(df.test)-1),],df.test[-c(1,nrow(df.test)),],df.test[-c(1,2),]))
  setnames(df.test,c("State2","time1","State1","time2","Stateactual","time3"))
  test <- df.test[,.(State2,State1,Stateactual)]

  test2 <- test[,.N,.(State2,State1,Stateactual)]
  return(test2)
  
}

EstimateMarkovChain <- function(listtest){

  #fits one COMPLETE markov chain, including transition probabilities using the list. 
  #Coefficients are computed using the MLE.
  
  #setting up df.
  df.test <- as.data.frame(listtest)
  setnames(df.test,c("state","time"))
  TransitionMatrix <- matrix(0,nrow=3,ncol=3)
  
  df.test <- as.data.table(cbind(df.test[-nrow(df.test),],df.test[-1,]))
  setnames(df.test,c("Statebegin","time1","Stateend","time2"))
  df.test[,TimeSpent := time2-time1]
  
  
  #Compute Exponential Proba
  
  df.test.exp <- copy(df.test)
  df.test.exp[,TimeSpent := sum(TimeSpent),.(Statebegin)]
  df.test.exp[,Timevisited := .N,.(Statebegin)]
  df.test.exp <- unique(df.test.exp[,.(TimeSpent,Statebegin,Timevisited)])
  df.test.exp[,lambda:=Timevisited/TimeSpent]
  
  #Compute transitional proba
  
  df.test[, numberoftrans := .N,.(Statebegin,Stateend)]
  df.test[, totalnumberoftrans := .N,.(Statebegin)]
  df.test[,ProbaTrans := numberoftrans/totalnumberoftrans]
  df.test <- unique(df.test[,.(ProbaTrans,Statebegin,Stateend)])
  
  #filling up matrix
  
  dimnames(TransitionMatrix) <-list( c("-1","0","1"),c("-1","0","1"))
  for(i in c("-1","0","1")){
  TransitionMatrix[i,i] <- ifelse(length(df.test.exp[Statebegin==i]$lambda)==0,NA,df.test.exp[Statebegin==i]$lambda)
    for(j in c("-1","0","1")){
      if(i!=j & nrow(df.test[Statebegin==i & Stateend ==j])>0){
        TransitionMatrix[i,j] <- ifelse(length(df.test[Statebegin==i & Stateend ==j]$ProbaTrans)==0,0,df.test[Statebegin==i & Stateend ==j]$ProbaTrans)
      }
    }
  }
  return(TransitionMatrix)
}

FormatTransitionMatrix <- function(matrix){

  # function to format the transition matrix obtained by the MLI estimation to the one that rctmc can use
  # i.e. one the actual correct transition matrix, with negative coeff on the diagonal and the qik/qii otherwise

  for(i in 1:3){
    for(j in 1:3){
      if(i!=j){
        matrix[i,j] <- matrix[i,i]*matrix[i,j]
      }
    }
  }
  diag(matrix) <- -1*diag(matrix)
  return(matrix)
}

BootstrapSimulation <- function(matrix,n){

  # function to simulate $n$ number of time the markov chain mctest.
  # returns a list with useful statistics about the bootstrap

  bootsumm <- list()
  bootsumm$NumberOfTransitions <- c()
  bootsumm$summary <- data.frame()
  mctest <- new("ctmc", states = c("-1","0","1"), byrow = T, generator = matrix)
  for(i in 1:n){
    mcdf <-rctmc(200,mctest,T=1190,out.type = "df")
    summ <- GetSummarydf(list(as.character(mcdf$states),as.numeric(as.character(mcdf$time))))
    bootsumm$NumberOfTransitions <- c(bootsumm$NumberOfTransitions,summ$NumberOfTransitions)
    bootsumm$summary <- rbind(bootsumm$summary,summ$summary)
  }
  bootsumm$MeanTimeVisited <-  as.data.table(summ$summary)[,mean(TimeVisited),.(Statebegin)]
  bootsumm$MeanTimeSpent  <- as.data.table(summ$summary)[,mean(MeanTimeSpent),.(Statebegin)]
  bootsumm$MeanTransitions <- mean(summ$NumberOfTransitions)
  return(bootsumm)
}

EstimateSeveralMarkovChains <- function(regions = 1:90,patient=1,mean=T,data=HCP.bin){

  #fits several COMPLETE markov chains, including transition probabilities. 
  #Coefficients are computed using the MLE.
  #if mean is set to TRUE, coefficients are average across patients

  coeff <- array(0,dim=c(3,3,length(regions),length(patient)))
  for(j in 1:length(patient)){
    cat(paste("doing patient",j,"\n"))
    for(i in 1:length(regions)){
      listtest <- CreateMarkovChainList(data[regions[i],,patient[j]])
      EasyMatrix <- EstimateMarkovChain(listtest)
      coeff[,,i,j] <- as.array(EasyMatrix)
    }
  }
  if(mean)  coeff <- apply(coeff,c(1,2,3),mean)
  return(coeff)
}


DotPlotCluster<- function(coef){

  # plots the dot plot colored by cluster of the coefficients of the markov chain stored in coef.
  # usually used for one single patient or the mean of the coefficients.
  
  # I - setting up the dataframe
  coef.names <- c("Lambda_-1","p_0-1","p_1-1","p_-10","Lambda_0","p_10","p_-11","p_01","Lambda_1")
  plot.df <- data.table()
  
  for(i in 1:90){
    plot.df.temp <- cbind(rep(Nodes[index==i]$name,9),rep(Nodes[index==i]$cluster,9))
    coef.vector <- c(coef[,,i])
    plot.df.temp <- cbind(plot.df.temp,coef.vector,coef.names)
    plot.df <- rbind(plot.df,plot.df.temp)
  }
  
  setnames(plot.df,c("Region","Cluster","Coefc","CoefName"))
  plot.df[,CoefValue := as.numeric(Coefc)]
  plot.df <- plot.df[CoefName %in% c("Lambda_-1","p_01","Lambda_0","Lambda_1")]
  
  # II - plotting
  plot1 <- ggplot(data=plot.df) + geom_dotplot(aes(x=1,y=CoefValue,fill=as.factor(Cluster)),stackgroups=T,
                    binpositions="all",binaxis = "y",stackdir="center",dotsize = 1) + 
   scale_fill_brewer(palette="Set1", name = "Cluster") + facet_wrap(~CoefName,scales="free") + xlab('') + 
    ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
  ggsave(filename = paste("./Plots/dotPlot_MarkovCoefCluster",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 9,
         height = 8,
         units = "in")
   print(plot1)
   return(plot.df)
}



DotPlotClusterCut<- function(coef,coef.name,ncut){

  # plots the dot plot and the line plots colored by cluster, faceted by sections (with $ncut$ cuts), of the $coef.name$ coefficient of the markov chain stored in coef.
  # usually used for one single patient or the mean of the coefficient.
  # the coef.name can be one of ("Lambda_-1","p_01","Lambda_0","Lambda_1")
  
  # I - setting up the dataframe
  coef.names <- c("Lambda_-1","p_0-1","p_1-1","p_-10","Lambda_0","p_10","p_-11","p_01","Lambda_1")
  plot.df <- data.table()
  
  for(i in 1:90){
    for(j in 1:ncut){
      plot.df.temp <- cbind(rep(Nodes[index==i]$name,9),rep(Nodes[index==i]$cluster,9),rep(j,9))
      coef.vector <- c(coef[,,i,j])
      plot.df.temp <- cbind(plot.df.temp,coef.vector,coef.names)
      plot.df <- rbind(plot.df,plot.df.temp)
    }
  }
  
  setnames(plot.df,c("Region","Cluster","Section","Coefc","CoefName"))
  plot.df[,CoefValue := as.numeric(Coefc)]
  plot.df <- plot.df[CoefName == coef.name]
  
  # II - plotting
  
   plot1 <- ggplot(data=plot.df) + geom_dotplot(aes(x=1,y=CoefValue,fill=as.factor(Cluster)),stackgroups=T,
                                                binpositions="all",binaxis = "y",stackdir="center",dotsize = 1) +
     scale_fill_brewer(palette="Set1", name = "Cluster") + facet_grid(~Section,scales="free") + xlab('') +
     ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
   
  
  plot2 <- ggplot(data=plot.df) + geom_line(aes(x=Section,y=CoefValue,color=as.factor(Cluster),group=Region),size=1) +
    scale_color_brewer(palette="Set1", name = "Cluster") + xlab('') +
    ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
  ggsave(filename = paste("./Plots/dotPlot_MarkovCoefClusterCut",coef.name,format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 12,
         height = 8,
         units = "in")
  
  ggsave(filename = paste("./Plots/linePlot_MarkovCoefClusterCut",coef.name,format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot2,
         width = 12,
         height = 8,
         units = "in")
  print(plot1)
  print(plot2)
  return(plot.df)
}


ComputeBootstrapTestTopBottom <- function(i,j){

  # function to compute the boostrap tests between the top8 and bottom 8 regions
  # returns a $p$ value matrix

  top.regions <- c("OLF.L","OLF.R","AMYG.L","AMYG.R","PAL.L","PAL.R","THA.L","THA.R")
  bottom.regions <- c("PreCG.L","PreCG.R", "SOG.R","IOG.R","PoCG.L", "PoCG.R", "PCL.L","PCL.R")
  top.index<-Nodes[name %in% top.regions]$index
  bottom.index<-Nodes[name %in% bottom.regions]$index
  
  p.matrix <- matrix(0,nrow=length(top.index),ncol=length(bottom.index))
  
  
  for(k in 1:length(top.index)){
    for(l in 1:length(bottom.index)){
      p.matrix[k,l] <- BootstrapTest(coeff[i,j,top.index[k],],coeff[i,j,bottom.index[l],])
    }
  }
  return(p.matrix)
}


ComputetTestSection <- function(i,j,coeff){

  #function to compute the bootrap test p value and generate the statistics to test whether or not there is a
  # significant differences in coef value between section i and section j
  
  p.vector <- as.data.table(1:90)
  setnames(p.vector,"region")
  p.vector[,cluster := Nodes[index==region]$cluster,.(region)]
  p.vector[,`Lambda_-1` := BootstrapTest(coeff[1,1,region,,i],coeff[1,1,region,,j]),.(region)]
  p.vector[,`Lambda_1` := BootstrapTest(coeff[3,3,region,,i],coeff[3,3,region,,j]),.(region)]
  p.vector[,`Lambda_0` := BootstrapTest(coeff[2,2,region,,i],coeff[2,2,region,,j]),.(region)]
  
  p.vector[,`:=`(Max_0=max(Lambda_0),Max_1=max(Lambda_1),`Max_-1`=max(`Lambda_-1` )),.(cluster)]
  p.vector[,`:=`(Min_0=min(Lambda_0),Min_1=min(Lambda_1),`Min_-1`=min(`Lambda_-1` )),.(cluster)]
  p.vector[`Lambda_-1`< 0.05,`NumberSign_-1` := as.numeric(.N),.(cluster)]
  p.vector[`Lambda_0`< 0.05,`NumberSign_0` :=as.numeric(.N),.(cluster)]
  p.vector[`Lambda_1`< 0.05,`NumberSign_1` := as.numeric(.N),.(cluster)]
  p.vector[,numbersign:=(`Lambda_-1`< 0.05)+(`Lambda_0`< 0.05)+(`Lambda_1`< 0.05)]
  p.vector[,`NumberSign_-1` := `NumberSign_-1`/.N *100,.(cluster)]
  p.vector[,`NumberSign_0` := `NumberSign_0`/.N *100,.(cluster)]
  p.vector[,`NumberSign_1` := `NumberSign_1`/.N *100,.(cluster)]
  cat("TotalNumberOfSign ",sum(p.vector$numbersign)/(nrow(p.vector)*3)*100,"%")
  p.vector <- unique(na.omit(p.vector[,.(cluster,Max_0,Max_1,`Max_-1`,Min_0,Min_1,`Min_-1`,
                                 `NumberSign_1`,`NumberSign_-1`,`NumberSign_0`)]))
  
 return(p.vector)
}


ComputeTestCluster <- function(i,j,clusteri){

  #this function compute the bootstrap test to check wheter or not the coefficient are the same within clusters.
  # remove region 36 if selected cluster is 2 as an outlier
  # select the coefficient to test with i and j, and the cluster with clusteri

  index <- Nodes[cluster==clusteri & index != 36]$index #cluster 2
  
  p.matrix <- matrix(0,nrow=length(index),ncol=length(index))
  
  
  for(k in 1:length(index)){
    for(l in 1:length(index)){
      if(k!=l)  p.matrix[k,l] <- BootstrapTest(coeff[i,j,index[k],],coeff[i,j,index[l],])
    }
  }
  return(p.matrix)
}

OverallSummary <- function(){

  # returns some summary statistics about all the regions' transitions and time spent

  summ.df <-data.table()
  summ.vec <- c()
  for(j in 1:90){
    cat(paste("doing patient",j,"\n"))
    for(i in 1:90){
      summtemp <- GetSummarydf(CreateMarkovChainList(HCP.bin[i,,j])) 
      summ.df <- rbind(summ.df,cbind(summtemp$summary,rep(i,3),rep(j,3)))
      summ.vec <- c(summ.vec,summtemp$NumberOfTransitions)
    }
  } 
  summary <- list()
  
  setnames(summ.df,c("TimeVisited","MeanTimeSpent","Statebegin","RegionIndex","PatientIndex"))
  summ.df <- merge(summ.df,Nodes[,.(name,index,cluster)],by.x="RegionIndex",by.y="index",all.x=T)
  summary$NumberOfTransitions <- summ.vec
  summary$df <- summ.df
  return(summary)
}

FitExp <- function(region,patient,state){

  # tests whether the exponential distribution is a good fit for the time spent in a given state for region and patient
  # plots the QQ-plots and the density plot as diagnostic

  listtest <- CreateMarkovChainList(HCP.bin[region,,patient])
  df.test <- as.data.frame(listtest)
  setnames(df.test,c("state","time"))
  df.test <- as.data.table(cbind(df.test[-nrow(df.test),],df.test[-1,]))
  setnames(df.test,c("Statebegin","time1","Stateend","time2"))
  df.test[,TimeSpent := time2-time1]
  df.test <- df.test[Statebegin == state,.(TimeSpent)]
  setkey(df.test,"TimeSpent")
  
  n <- nrow(df.test)
  fit1 <-fitdist(df.test$TimeSpent,distr="exp")
  qq <- qexp((1:n)/(n+1),rate=fit1$estimate[["rate"]])
  vec <- 1:(max(qq)+2)
  vec2 <- 1+1:300*(max(df.test$TimeSpent)-1)/300

  density2 <- dexp(vec2,rate=fit1$estimate[["rate"]])
  #qqplot
  plot1 <- ggplot() + geom_point(aes(y=df.test$TimeSpent,x=qq)) + geom_line(aes(x=vec,y=vec),color="red") +
    xlab("Theoretical Quantiles") + ylab("Actual Quantiles") + theme_classic()
  ggsave(filename = paste("./Plots/qqplot_exp_",state,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 9,
         height = 8,
         units = "in")
  print(plot1)
  
  #density
  plot2 <- ggplot() + geom_density(aes(x=df.test$TimeSpent)) + geom_line(aes(x=vec2,y=density2),color="red") + 
    xlab("TimeSpent") + ylab("Density") + theme_classic()
  ggsave(filename = paste("./Plots/density_exp_",state,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot2,
         width = 9,
         height = 8,
         units = "in")
  
  print(plot2)
  cat("aic for exp", fit1$aic,"\n")
  cat("bic for exp", fit1$bic,"\n")
  
}

LogNormalfit <- function(region,patient,state){

  # tests whether the lognormal distribution is a good fit for the time spent in a given state for region and patient
  # plots the QQ-plots and the density plot as diagnostic

  listtest <- CreateMarkovChainList(HCP.bin[region,,patient])
  df.test <- as.data.frame(listtest)
  setnames(df.test,c("state","time"))
  df.test <- as.data.table(cbind(df.test[-nrow(df.test),],df.test[-1,]))
  setnames(df.test,c("Statebegin","time1","Stateend","time2"))
  df.test[,TimeSpent := time2-time1]
  df.test <- df.test[Statebegin == state,.(TimeSpent)]
  setkey(df.test,"TimeSpent")
  
  n <- nrow(df.test)
  
  #lognormal
  fit1 <- fitdist(df.test$TimeSpent,dist="lnorm")
  qq <- qlnorm((1:n)/(n+1),meanlog=fit1$estimate[["meanlog"]],sdlog=fit1$estimate[["sdlog"]])
  vec <- 1:(max(qq)+2)
  vec2 <- 1+1:300*(max(df.test$TimeSpent)-1)/300
  density2 <- dlnorm(vec2,meanlog=fit1$estimate[["meanlog"]],sdlog=fit1$estimate[["sdlog"]])
  
  #qqplot
  plot1 <- ggplot() + geom_point(aes(y=df.test$TimeSpent,x=qq)) + geom_line(aes(x=vec,y=vec),color="red") +
    xlab("Theoretical Quantiles") + ylab("Actual Quantiles") + theme_classic()
  ggsave(filename = paste("./Plots/otherqqplot_lognormal_",state,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 9,
         height = 8,
         units = "in")
  print(plot1)
  
  #density
  plot2 <- ggplot() + geom_density(aes(x=df.test$TimeSpent)) + geom_line(aes(x=vec2,y=density2),color="red") + 
    xlab("TimeSpent") + ylab("Density") + theme_classic()
  ggsave(filename = paste("./Plots/otherdensity_lognormal_",state,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot2,
         width = 9,
         height = 8,
         units = "in")
  
  print(plot2)
  
  cat("aic for lnorm", fit1$aic,"\n")
  cat("bic for lnorm", fit1$bic,"\n")
}

Gammafit <- function(region,patient,state){

  # tests whether the Gamma distribution is a good fit for the time spent in a given state for region and patient
  # plots the QQ-plots and the density plot as diagnostic 

  listtest <- CreateMarkovChainList(HCP.bin[region,,patient])
  df.test <- as.data.frame(listtest)
  setnames(df.test,c("state","time"))
  df.test <- as.data.table(cbind(df.test[-nrow(df.test),],df.test[-1,]))
  setnames(df.test,c("Statebegin","time1","Stateend","time2"))
  df.test[,TimeSpent := time2-time1]
  df.test <- df.test[Statebegin == state,.(TimeSpent)]
  setkey(df.test,"TimeSpent")
  
  n <- nrow(df.test)
  
  #lognormal
  fit1 <- fitdist(df.test$TimeSpent,dist="gamma")
  qq <- qgamma((1:n)/(n+1),shape=fit1$estimate[["shape"]],rate=fit1$estimate[["rate"]])
  vec <- 1:(max(qq)+2)
  vec2 <- 1+1:300*(max(df.test$TimeSpent)-1)/300
  density2 <- dgamma(vec2,shape=fit1$estimate[["shape"]],rate=fit1$estimate[["rate"]])
  
  #qqplot
  plot1 <- ggplot() + geom_point(aes(y=df.test$TimeSpent,x=qq)) + geom_line(aes(x=vec,y=vec),color="red") +
    xlab("Theoretical Quantiles") + ylab("Actual Quantiles") + theme_classic()
  ggsave(filename = paste("./Plots/otherqqplot_lognormal_",state,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 9,
         height = 8,
         units = "in")
  print(plot1)
  
  #density
  plot2 <- ggplot() + geom_density(aes(x=df.test$TimeSpent)) + geom_line(aes(x=vec2,y=density2),color="red") + 
    xlab("TimeSpent") + ylab("Density") + theme_classic()
  ggsave(filename = paste("./Plots/otherdensity_lognormal_",state,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot2,
         width = 9,
         height = 8,
         units = "in")
  
  print(plot2)
  
  cat("aic for gamma", fit1$aic,"\n")
  cat("bic for gamma", fit1$bic,"\n")
}

CompareRatioTimeSpent <- function(summ){
 
 # plots the ratio of time spent in 0 compared to time spent in 1 and -1.
 # largely unused. takes the output of OverallSummary() as input.

  summ.df <- summ$df
  summ.df <- dcast(summ.df,formula = RegionIndex + PatientIndex + cluster ~ Statebegin, value.var="MeanTimeSpent") 
  summ.df[,c("ratio-1","ratio1") := .(`0`/`-1`,`0`/`1`)]
  
  plot1 <- ggplot(data=summ.df[PatientIndex==1]) + geom_dotplot(aes(x=1,y=`ratio1`,fill=as.factor(cluster)),stackgroups=T,
                                               binpositions="all",binaxis = "y",stackdir="center",dotsize = 1) + 
    scale_fill_brewer(palette="Set1", name = "Cluster") + xlab('') + 
    ylab("Ratio Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
  ggsave(filename = paste("./Plots/dotPlot_RatioTimeSpent",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 9,
         height = 8,
         units = "in")
  print(plot1)
  return(plot.df)
} 

WeibullFit<- function(region,patient,state){

  # tests whether the Weibull distribution is a good fit for the time spent in a given state for region and patient
  # plots the QQ-plots and the density plot as diagnostic

  listtest <- CreateMarkovChainList(HCP.bin[region,,patient])
  #setting up df.
  df.test <- as.data.frame(listtest)
  setnames(df.test,c("state","time"))
  df.test <- as.data.table(cbind(df.test[-nrow(df.test),],df.test[-1,]))
  setnames(df.test,c("Statebegin","time1","Stateend","time2"))
  df.test[,TimeSpent := time2-time1]
  df.test <- df.test[Statebegin == state,.(TimeSpent)]
  setkey(df.test,"TimeSpent")
  
  fit1 <- fitdist(df.test$TimeSpent,distr="weibull")
  n <- nrow(df.test)
  #weibull
  qq <- qweibull((1:n)/(n+1),shape=fit1$estimate[["shape"]],scale=fit1$estimate[["scale"]])
  vec <- 1:(max(qq)+2)
  vec2 <- 1+1:300*(max(df.test$TimeSpent)-1)/300
  density2 <- dweibull(vec2,shape=fit1$estimate[["shape"]],scale=fit1$estimate[["scale"]])
  
  #qqplot
  plot1 <- ggplot() + geom_point(aes(y=df.test$TimeSpent,x=qq)) + geom_line(aes(x=vec,y=vec),color="red") +
    xlab("Theoretical Quantiles") + ylab("Actual Quantiles") + theme_classic()
  ggsave(filename = paste("./Plots/otherqqplot_weibull_",state,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 9,
         height = 8,
         units = "in")
  print(plot1)
  
  #density
  plot2 <- ggplot() + geom_density(aes(x=df.test$TimeSpent)) + geom_line(aes(x=vec2,y=density2),color="red") + 
    xlab("TimeSpent") + ylab("Density") + theme_classic()
  ggsave(filename = paste("./Plots/otherdensity_weibull_",state,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot2,
         width = 9,
         height = 8,
         units = "in")
  
  print(plot2)
  cat("aic for weibull", fit1$aic,"\n")
  cat("bic for weibull", fit1$bic,"\n")
  
}

CompareFit <- function(region,patient,state){

  # Wrapper to test all 3 fits at the same time

  LogNormalfit(region,patient,state)
  WeibullFit(region,patient,state)
  FitExp(region,patient,state)
}

MarkovCut <- function(timeseries,numberofcut,plot=F){

  # generate the list of states and time spent for the stationarity analysis for one time series.
  # i.e. cuts the series into $numberofcut$ sections.
  # The function then fits a markov chain for each section. outputs an array with all the coefficients
  # if plot is set to TRUE, then it also plots the time series along with vertical cut at the sections breakpoints 

  dimn <- list(c("-1","0","1"),c("-1","0","1"),NULL)
  #generate the list containt the n cut time series
  len <- floor(length(timeseries)/numberofcut)
  ts <- list()
  for(i in 0:(numberofcut-2)){
    ts[[i+1]] <- timeseries[(1+i*len):((i+1)*len)]
  }
  ts[[numberofcut]] <- timeseries[((numberofcut-1)*len):length(timeseries)]
  
  #generates the list containing the markov chain coeff
  markovcut <- array(0,dim=c(3,3,numberofcut))
  for(i in 1:numberofcut){
    markovcut[,,i] <- as.array(EstimateMarkovChain(CreateMarkovChainList(ts[[i]])))
  }
  if(plot){
    plot1 <- PlotTimeSeries(timeseries,save=F)
    df <- as.data.table(cbind(rep(max(timeseries*1.1),numberofcut-1),rep(min(timeseries*1.1),numberofcut-1),
                              length(timeseries)/numberofcut*1:(numberofcut-1)))
    plot1 <- plot1+geom_segment(data=df,aes(x=V3,xend=V3,y=V1,yend=V2),color="blue",size=2)
    ggsave(filename = paste("./Plots/PlotTimeSeriesCut",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
           plot = plot1,
           width = 12,
           height = 8,
           units = "in") 
    print(plot1)
  }
  dimnames(markovcut) <- dimn
  return(markovcut)
}

MarkovCutAll <- function(regions = 1:90,patient=1,numberofcut,mean=T){

  #fits several markov chains, cutting the time series in $numberofcut$ sections each time.
  #if mean is set to TRUE, returns the mean of the coefficient across patients, otherwise returns the full array.

  coeff <- array(0,dim=c(3,3,length(regions),length(patient),numberofcut))
  coefflist <- list()
  for(j in 1:length(patient)){
    cat(paste("doing patient",j," "))
    for(i in 1:length(regions)){
      cat(".")
      mc <- MarkovCut(HCP.bin[regions[i],,patient[j]],numberofcut,F)
      coeff[,,i,j,] <- mc
    }
    cat("\n")
  }
  if(mean)  coeff <- apply(coeff,c(1,2,3,5),function(x) mean(x,na.rm=T))
  return(coeff)
}



DotPlotClusterMCI<- function(MarkovCoeffMCI){

  # plots the dot plot colored by cluster of the coefficients for the MCI dataset of the markov chain stored in MarkovCoeffMCI.
  # usually used for one single patient or the mean of the coefficients.

  # I - setting up the dataframe
  coef.names <- c("Lambda_-1","p_0-1","p_1-1","p_-10","Lambda_0","p_10","p_-11","p_01","Lambda_1")
  plot.df <- data.table()
  Nodes2 <- Nodes[!(index %in% c(75,76))]
  setkey(Nodes2,index)
  
  coef <- array(0,dim=c(3,3,88,3))
  coef[,,,1] <- apply(MarkovCoeffMCI$sCTL,c(1,2,3),mean)
  coef[,,,2] <- apply(MarkovCoeffMCI$dCTL,c(1,2,3),mean)
  coef[,,,3] <- apply(MarkovCoeffMCI$MCI,c(1,2,3),mean)
  
  for(j in 1:3){
    for(i in 1:88){
      ii <- c(1:74,77:90)[i]
      plot.df.temp <- cbind(rep(j,9),rep(Nodes2[index==ii]$name,9),rep(Nodes2[index==ii]$cluster,9))
      coef.vector <- c(coef[,,i,j])
      plot.df.temp <- cbind(plot.df.temp,coef.vector,coef.names)
      plot.df <- rbind(plot.df,plot.df.temp)
    }
  }
  
  setnames(plot.df,c("Group","Region","Cluster","Coefc","CoefName"))
  plot.df[Group=="1",Group:="sCTL"]
  plot.df[Group=="2",Group:="dCTL"]
  plot.df[Group=="3",Group:="MCI"]
  plot.df[,Group:=factor(Group,levels = c("sCTL","dCTL","MCI"))]
  plot.df[,CoefValue := as.numeric(Coefc)]
  plot.df <- plot.df[CoefName %in% c("Lambda_-1","p_01","Lambda_0","Lambda_1")]
  
  # II - plotting
  for(cc in c("Lambda_-1","p_01","Lambda_0","Lambda_1")){
  plot1 <- ggplot(data=plot.df[CoefName==cc]) + geom_dotplot(aes(x=1,y=CoefValue,fill=as.factor(Cluster)),stackgroups=T,
                                               binpositions="all",binaxis = "y",stackdir="center",dotsize = 1) + 
    scale_fill_brewer(palette="Set1", name = "Cluster") + facet_wrap(~Group) + xlab('') + 
    ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
  ggsave(filename = paste("./Plots/dotPlot_MarkovCoefCluster",cc,format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 12,
         height = 8,
         units = "in")
  print(plot1)
  }
  return(plot.df)
}

ComputetTestMCI <- function(MarkovCoeffMCI){

  #function to compute the bootstrap test p value and generate the statistics for the MCI dataset
  # test whether or not the coefficients of regions vary from patient pool to patient pool.

  setkey(Nodes,index)
  Nodes2<-Nodes[!(index%in%c(75,76))]
  p.vector <- as.data.table(1:88)
  setnames(p.vector,"region")
  p.vector[,cluster := Nodes[index==region]$cluster,.(region)]
  p.vector[,`Lambda_-1` := mean(BootstrapTest(MarkovCoeffMCI$MCI[1,1,region,],MarkovCoeffMCI$dCTL[1,1,region,])),.(region)]
  p.vector[,`Lambda_1` :=  mean(BootstrapTest(MarkovCoeffMCI$MCI[3,3,region,],MarkovCoeffMCI$dCTL[3,3,region,])),.(region)]
  p.vector[,`Lambda_0` :=  mean(BootstrapTest(MarkovCoeffMCI$MCI[2,2,region,],MarkovCoeffMCI$dCTL[2,2,region,])),.(region)]
  p.vector <- merge(p.vector,Nodes[,.(index,name)],by.x="region",by.y="index")


  return(p.vector)
}

MarkovMovingLambda <- function(region = 43,patient=1:15,length=200,mean=T){

  # computes the moving lambdas, cutting each time the time series into section of $length$ unit of time
  # if mean is set to TRUE, returns the mean of the coefficient across patients, otherwise, returns the full array
  # this is used for Movie watching data only ! change the dataset inside CreateMarkovChainList to adapt

  coeff <- array(0,dim=c(3,3,length(patient),1290-length))
  dimnames(coeff) <-list( c("-1","0","1"),c("-1","0","1"),NULL,NULL)
  for(j in 1:length(patient)){
    cat(paste("doing patient",j," "))
    for(i in (length+1):1290){
      cat(".")
      mc <- EstimateMarkovChain(CreateMarkovChainList(MovieData$fullbinarized[region,(i-length):i,patient[j]]))
      coeff[,,j,(i-length)] <- mc
    }
    cat("\n")
  }
  if(mean)  coeff <- apply(coeff,c(1,2,4),function(x) mean(x,na.rm=T))
  return(drop(coeff))
}

MovingCoeffPlot <- function(MovingCoeff,length,name){

  # creates the Moving lambda plots. Input MovingCoeff is obtained from the MarkovMovingLambda() function
  # This functions is meant to be used with the Moviewatching dataset.

  dimss <- dim(MovingCoeff)
  index <- 1:(dimss[3])
  for(i in c("1","-1","0")){
    plot.df <- as.data.table(cbind(index, MovingCoeff[i,i,]))
    setnames(plot.df,c("Time","Value"))
    plot.df[,group:=1]
    plot.df[,Type := "Transition"]
    plot.df[(Time+length) %in% MovieData$TimePointsMovie,Type := "Movie"]
    plot.df[(Time+length) %in% MovieData$TimePointsRest,Type := "Rest"]
    plot1 <- ggplot(data=plot.df) + geom_line(aes(x=Time,y=Value,color=Type,group=group),size=1) +
      scale_color_manual(values = c(brewer.pal(3,"Set1")[1:2],"black")) + theme_minimal() +
      ggtitle(paste("Moving Lambda",i))
    ggsave(filename = paste("./Plots/PlotMovingCoeff_",name,"_",i,".png", sep = ""),
           plot = plot1,
           width = 12,
           height = 8,
           units = "in") 
    print(plot1)
  }
}

####EXECUTION EXAMPLES #####

## Generate the Markov Coeff lists for the MCI dataset
#MarkovCoeffMCI <- list()
#MarkovCoeffMCI$sCTL <- EstimateSeveralMarkovChains(regions=1:88,patient=1:77,mean=F,data=MCIData$sCTLbin)
#MarkovCoeffMCI$dCTL <- EstimateSeveralMarkovChains(regions=1:88,patient=1:66,mean=F,data=MCIData$dCTLbin)
#MarkovCoeffMCI$MCI <- EstimateSeveralMarkovChains(regions=1:88,patient=1:69,mean=F,data=MCIData$MCIbin)
#save(MarkovCoeffMCI, file = "./Data/MarkovCoeffMCI.rda")

## plots the moving lambdas for the movie watching for patient 1, region 22 and length of 50
#MovingCoeff <- MarkovMovingLambda(22,1,50,F)
#MovingCoeffPlot(MovingCoeff,50,22)

## Fits and simulate a Markov chain to region 58 of patient 27
#listtest <- CreateMarkovChainList(HCP.bin[58,,27])
#EasyMatrix <- EstimateMarkovChain(listtest)
#MatrixEstimate <- FormatTransitionMatrix(EasyMatrix)
#Summ <- BootstrapSimulation(MatrixEstimate,50)

## Plots the dot plots for the stationarity analysis for the mean patient
#mcall <- MarkovCutAll(1:90,1:90,4,T)
#DotPlotClusterCut(mcall,"Lambda_0",4)
#DotPlotClusterCut(mcall,"Lambda_1",4)
#DotPlotClusterCut(mcall,"Lambda_-1",4)
#DotPlotClusterCut(mcall,"p_01",4)
