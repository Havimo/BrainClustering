# 
# 
# This .R file is the base file for the Thesis' R code
# It contains data acquisition and transformation as well as simple plotting function.
# Source this file first in order to use the other R files pertaining to this project
# 
# The init() function will automatically source all of the other .R files containing the functions
# each of these files have some simple working example commented out to get a quick understanding on 
# what they do
#
#
#

#### Initilazing ####


init <- function(workingdir=""){

  #packages
  library(R.matlab)
  library(randomForest)
  library(data.table)
  library(ggplot2)
  library(corrplot)
  library(ppcor)
  library(RColorBrewer)
  library(Matrix)
  library(qgraph)
  library(cluster)
  library(NbClust)
  library(factoextra)
  library(network)
  library(ggdendro)
  library(markovchain)
  library(fitdistrplus)
  
  #working directory
  if(workingdir=="" & dir.exists("D:/Google Drive/")){
  workingdir <- "D:/Google Drive/Master/Thesis"
  } else if (workingdir==""){
  workingdir <- "D:/GoogleDrive/Master/Thesis"
  }
  setwd(workingdir)

  #sources the other .R files containing the functions
  source("./R Code/HMM.R")
  source("./R Code/ClusterAnalysis.R")
  source("./R Code/HMM_connection_6_states.R")
  source("./R Code/MCIClassification.R")
  source("./R Code/ClusterConnections.R")


  #set some global variables used throughout the project and the functions.
  #the .rda files are loaded for reasons of speed but can be recreated through the functions 
  index <<- unique(1:89%/%2)*2+1
  Nodes <<- as.data.table(read.csv("./Data/Nodes_MetaData.csv",stringsAsFactors = F))
  load("./Data/ListofAllPartCorrBin.rda",envir=.GlobalEnv)
  load("./Data/ListofAllPartCorr.rda",envir=.GlobalEnv)
  load("./Data/FinalBoostrapClustering.rda",envir=.GlobalEnv)
  load("./Data/MarkovChainLambdas4CutCoeffConnections.rda",envir=.GlobalEnv)
  
  #merge the final cluster distribution with the Nodes meta data. Allows for easier plotting and stats.
  Nodes <<- merge(Nodes,data.table(FinalBoostrapClustering$cluster,row.names(as.data.frame(FinalBoostrapClustering$cluster)))
                 ,by.x="name",by.y="V2")
  setnames(Nodes,"V1","cluster")
}


#### General Utiliy Functions ####

is.even <- function(num){
  bool <- (num/2 == floor(num/2))
  return(bool)
}

randn <- function(lower,upper,n=1){
  return(floor(runif(n)*(upper-lower + 1))+1)
}

Extract_Symmpair_Coef <- function(CorrMatrix){
  coef <- diag(CorrMatrix[index,(index+1)])
  names(coef) <- Nodes$name[index]
}

#### Data Loading ####

LoadData <- function(){

  #Main dataset loading

  files=list.files("./Data/Mat/")
  N <<- 90 #Number of Region
  Tm <<- 1190 #Number of time points
  n <<- 90 #Number of patients
  HCP=array(NA, dim=c(N,Tm,n),dimnames = list(Region=Nodes$name,Time = 1:1190,Patient=1:90))
  
  for(k in 1:n) {
      HCP[,,k] <- readMat(paste("./Data/Mat/", files[k], sep=""),sparseMatrixClass="Matrix")$TCS
  }
  cat(paste("data was loaded successfuly with ", sum(is.na(HCP)), "NA values", "\n"))
  return(HCP)
}


LoadAdditionalData_MovieWatching <- function(){

  # Load movie watching data

  path="./Data/RS MW/"
  Names=list.files(paste(path,"Time courses/",sep=""))
  n=length(Names)
  
  data=array(0,dim=c(88,1290,n))
  for (k in 1:n){
    if(k !=15){
      data[,,k]<-as.matrix(read.table(paste(path,"Time courses/",Names[k], sep=""), sep=","))[-c(75,76),]
    }else{ #problem with 15 having one less time point.
      data[,1:1289,k]<-as.matrix(read.table(paste(path,"Time courses/",Names[k], sep=""), sep=","))[-c(75,76),]
      
    }
  }
  TimeInfo=readMat(paste(path, "Matlab data/scansIdx_rest_movies.mat",sep=""))
  
  nInov=9
  #TimeInfo$scansIdx.movies[[1]]
  dataMW=array(NA,dim=c(88,368,n)) #88regions, 368 timepoints, n patients
  MWTPs=NULL
  WinSizesMW=NULL
  InovMW=0
  
  for(Inov_ind in 1:nInov){
    MWTPs=c(MWTPs, as.numeric(unlist(TimeInfo$scansIdx.movies[[Inov_ind]]))) 
    WinSizesMW=c(WinSizesMW,length(unlist(TimeInfo$scansIdx.movies[[Inov_ind]])))
    InovMW=c(InovMW, sum(WinSizesMW))
  }
  InovMW=rev(rev(InovMW)[-1])
  dataMW=data[,MWTPs,]
  
  nInov=9
  dataRest=array(NA,dim=c(88,720,n)) #88regions, 720 timepoints, n patients
  RestTPs=NULL
  WinSizesRest=NULL
  InovRest=0
  
  for(Inov_ind in 1:nInov){
    RestTPs=c(RestTPs, as.numeric(unlist(TimeInfo$scansIdx.rest[[Inov_ind]]))) 
    WinSizesRest=c(WinSizesRest,length(unlist(TimeInfo$scansIdx.movies[[Inov_ind]])))
    InovRest=c(InovRest, sum(WinSizesRest))
  }
  dataRest=data[,RestTPs,]
  
  dataf<-list()
  dataf$full <- data
  dataf$fullbinarized <- BinarizeDataset(dataf$full,0.84,0.16)
  dataf$rest <- dataRest
  dataf$movie <- dataMW
  dataf$WinSizesRest <-WinSizesRest
  dataf$TimePointsMovie <- MWTPs
  dataf$TimePointsRest <- RestTPs
  dataf$InovRest <- InovRest
  dataf$WinSizesMW <- WinSizesMW
  dataf$InovMW <- InovMW

  return(dataf)
}

LoadAdditionalData_MCI <- function(){

  # Load time courses of MCI data

  matrices=readMat("./Data/MCI time courses/TCS.mat",sparseMatrixClass="Matrix")[[1]]
  sCTL=array(NA, dim=c(88,170,77))#88regions, 170 time points,77 patient
  dCTL=array(NA, dim=c(88,170,66))
  MCI=array(NA,  dim=c(88,170,69))
  # remove 2 regions
  RegiosToDiscard=c(75,76)
  for(k in   1:77) sCTL[,,k]=(matrices[[k]])[[1]][-RegiosToDiscard,]
  for(k in   1:66) dCTL[,,k]=(matrices[[(k+77)]])[[1]][-RegiosToDiscard,]
  for(k in   1:69)  MCI [,,k]=(matrices[[(k+143)]])[[1]][-RegiosToDiscard,]
  MCIf <- list()
  MCIf$sCTL <- sCTL
  MCIf$dCTL <- dCTL
  MCIf$MCI <- MCI
  MCIf$sCTLbin <- BinarizeDataset(sCTL,0.84,0.16)
  MCIf$dCTLbin <- BinarizeDataset(dCTL,0.84,0.16)
  MCIf$MCIbin <- BinarizeDataset(MCI,0.84,0.16)
  ConcordanceMCI <<- list()
  ConcordanceMCI$sCTL <<- ComputeConcordance(MCIf$sCTLbin)
  ConcordanceMCI$dCTL <<- ComputeConcordance(MCIf$dCTLbin)
  ConcordanceMCI$MCI <<- ComputeConcordance(MCIf$MCIbin)
  return(MCIf)
}

CreateDistanceMatrix <- function(){

  # Returns the physical distance matrix between the nodes

  Distance <- matrix(NA,nrow=90,ncol=90)
  for(i in 1:90){
    for(j in 1:90){
      Distance[i,j] <- sqrt( (Nodes$x[i]-Nodes$x[j])^2 + (Nodes$y[i]-Nodes$y[j])^2 + (Nodes$z[i]-Nodes$z[j])^2 )
    }
  }
  return(Distance)
}

#### Transform Data ####

#Transform data by normalzing it and applying the binary threshold transform

NormalizeDataset <- function(HCP){
  
  normalize <- function(X){  #create the normalization function
    X = X - mean(X)
    X = X/sd(X)
    return(X)
  }
  
#HCP.normalized <- apply(HCP,c(1,3),normalize) switches dimension around.
  
 for(i in 1:N){
    for(j in 1:n){
    #call normalization function for each line of HCP array, i.e
    #the normalization is done for each individual and is region specific
      HCP[i,,j]=normalize(HCP[i,,j])
    }
 }
  return(HCP)
}


BinarizeDataset <- function(HCP, upperb = NULL, lowerb = NULL){

  # Binarize the Dataset

  N <- dim(HCP)[1]
  n <- dim(HCP)[3]
  for(i in 1:N){
    for(j in 1:n){
      
      qt <- as.numeric(quantile(HCP[i,,j],c(lowerb,upperb)))
      l<- qt[1]
      u<- qt[2]
      HCP[i,,j]=ifelse(HCP[i,,j]<=l,-1,0)+ifelse(HCP[i,,j]>=u,1,0)
    }
  }
  return(HCP)
}

ThresholdDataset <- function(HCP, upperb = 0.8, lowerb = 0.2){

  # Threshold the dataset

  for(i in 1:N){
    for(j in 1:n){
      
      qt <- as.numeric(quantile(HCP[i,,j],c(lowerb,upperb)))
      l<- qt[1]
      u<- qt[2]
      HCP[i,,j]=ifelse(HCP[i,,j]<=l,HCP[i,,j],0)+ifelse(HCP[i,,j]>=u,HCP[i,,j],0)
    }
  }
  return(HCP)
}

ComputeConcordance<-function(HCPbin){

  # Computing accordance and discodance from the binarized dataset HCPbin

  dims <- dim(HCPbin)
  N <- dims[1]
  n <- dims[3]
  acc <- array(NA, dim = c(N,N,n))
  dis <- array(NA, dim = c(N,N,n))
  
  ComputeConcordanceForOnePatient <- function(HCPbin,patient,plot=T){
 
    X <- HCPbin[,,patient]
    Xu <- ifelse(X>0,1,0)
    Xl <- ifelse(X<0,-1,0)
    
    Scav=1/Tm*(Xu%*%t(Xu))
    Scdv=1/Tm*(Xl%*%t(Xl))
    Sa=Scav+Scdv
    E=diag(Sa)^(-0.5)           # Calculate the energy
    Acc=diag(E)%*%Sa%*%diag(E)  # Calculate the Accordance matrix
    
    Sd=1/Tm*(Xu%*%t(Xl)+Xl%*%t(Xu))
    Dis=diag(E)%*%Sd%*%diag(E)  # Calculate the Discordance matrix
    
    con=NULL
    
    con$acc=Acc
    con$dis=Dis
    
    if(plot){
      corrplot(con$acc)
      corrplot(con$dis)
    }
    
    return(con)
  }
  
  for(i in 1:n){
    conc <-  ComputeConcordanceForOnePatient(HCPbin,i,F)
    acc[,,i] <-  conc$acc
    dis[,,i] <- conc$dis
  }
  conc$acc <- acc
  conc$dis <- dis
  return(conc)
}

ComputePartialCorrelation <- function(HCP,patient){

  # Function to compute the partial correlation from the residuals of two linear regressions
  # this version takes ALL other variables into accout for the computation

  data <- HCP[,,patient]
  data <- as.data.table(t(data))
  
  partial.cor <- function(i,j){
    lm.region.1 <- lm(data=data[,-j,with=FALSE], formula = paste(names(data)[i]," ~ .",sep=""))
    residuals.region.1 <- residuals(lm.region.1)
    
    lm.region.2 <- lm(data=data[,-i,with=FALSE], formula =paste(names(data)[j]," ~ .",sep=""))
    residuals.region.2 <- residuals(lm.region.2)
    
    pcor = cor(residuals.region.1, residuals.region.2, method = 'pearson')
    
    return(pcor)
  }
  
  pcor <- matrix(0,ncol = 90,nrow = 90)  
  #loop
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      cat(paste("Computing Patient ",patient," cor for pair",i,j,"\n"))
      pcor[i,j] <- partial.cor(i,j)
    }
  }
  
  diag(pcor) <- NA
  pcor <- as.matrix(forceSymmetric(pcor))
  return(pcor)
}

ComputePartialCorrelation.Cluster <- function(HCP,patient){

  # Function to compute the partial correlation from the residuals of two linear regressions
  # this version takes only the variables into  the same cluster into accout for the computation

  data <- HCP[,,patient]
  data <- as.data.table(t(data))
  ClustIndex <- list()
  for(i in 1:7){
    ClustIndex[[i]]<- Nodes[cluster==i]$index
  }
  
  partial.cor <- function(i,j,dataclust){
    lm.region.1 <- lm(data=dataclust[,-j,with=FALSE], formula = paste(names(dataclust)[i]," ~ .",sep=""))
    residuals.region.1 <- residuals(lm.region.1)
    
    lm.region.2 <- lm(data=dataclust[,-i,with=FALSE], formula =paste(names(dataclust)[j]," ~ .",sep=""))
    residuals.region.2 <- residuals(lm.region.2)
    
    pcor = cor(residuals.region.1, residuals.region.2, method = 'pearson')
    
    return(pcor)
  }
  
  pcor <- matrix(0,ncol = 90,nrow = 90)  
  #loop
  for (i in 1:n){
    clusternum <- Nodes[index ==i]$cluster
    for (j in ClustIndex[[clusternum]]){
      if(i!=j){
        #get "new indices" of the regions inside its cluster df.
        ii <- which(ClustIndex[[clusternum]] == i)
        jj <- which(ClustIndex[[clusternum]] == j)
        cat(paste("Computing Patient ",patient," cor for cluster", clusternum, "and pair",i,j,"\n"))
        pcor[i,j] <- partial.cor(ii,jj,dataclust = data[, ClustIndex[[clusternum]],with=FALSE])
      }
    }
  }
  
  diag(pcor) <- NA
  pcor <- as.matrix(forceSymmetric(pcor))
  return(pcor)
}


MergePcorMatrices <- function(pcorcluster,pcorall){
  
  # merges the Partial Correlation Matrices, one only computed with the cluster, the other with everything.
  # used in the double correlation circle plot

  returnm<-ifelse(pcorcluster != 0,pcorcluster,pcorall)  
  return(returnm)
}


MeanPearson <- function(HCP){

  # Returns the mean Pearson correlation of the HCP dataset

  meanarray <- array(NA,dim=c(90,90,90))
  for(i in 1:90){
    meanarray[,,i] <- cor(t(HCP[,,i]))
  }
  return(apply(meanarray,c(1,2),mean))
}

BootstrapTest <- function(sample1,sample2,n=10000){

  #function to compute the bootstrap hypothesis testing 
  # test whether the mean of sample 1 is equal to the mean of sample 2
  # runs for $n$ iterations, return the first approximated p value

  sample1new <- sample1 - mean(sample1) + mean(sample2) #make sample 1 have mean of sample2
  sample2new <- sample2 - mean(sample2) + mean(sample1) #make sample 2 have mean of sample1
  bstrap1 <- bstrap2 <- c()
  
  for (i in 1:n){ #bootstrap
    bstrap1 <- c(bstrap1, mean(sample(sample1new,length(sample1),replace=T)))
    bstrap2 <- c(bstrap2, mean(sample(sample2new,length(sample2),replace=T)))
  }
  
  diff <- abs(mean(sample1)-mean(sample2))
  p.value.1 <- sum(bstrap1 >= mean(sample1))/n#proportion of more extreme cases
  p.value.2 <- sum(bstrap2 <= mean(sample2))/n
  return(p.value.1)
  
}

#### Plotting Functions ####

PlotDistanceCorrelation <- function(Corr,name,ylabel="Correlation"){

# Plots the Correlation Coefficient stored in Corr with respect to distance.

  dist <- CreateDistanceMatrix()
  dataplot <- as.data.table(cbind(c(dist),c(Corr), c(0:8099 %/% 90 + 1), c(0:8099 %% 90 +1)))
  setnames(dataplot,c('Distance','Correlation',"RA","RB"))
  dataplot[, Correlation := as.numeric(Correlation)]
  dataplot[, Distance := as.numeric(Distance)]
  dataplot[Correlation >= 0.8*max(Correlation) & Correlation <1 & RA<RB, Label := paste(RA,RB,sep=":")]
  plot1 <- ggplot(data=dataplot) + geom_point(aes(x=Distance,y=Correlation)) + 
     ylab(ylabel) # + geom_text(aes(x=Distance,y=Correlation,label=Label),angle = 45,nudge_x = 4, nudge_y = 0.05)
   
  ggsave(filename = paste("./Plots/DistanceCorrPlot_",name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 8,
         height = 8,
         units = "in")
  print(plot1)
}


PlotTimeSeries <- function(HCP,region=NULL,patient=NULL,save=T,name=""){

  #Plots the time series of HCP for the selected patient and region

  if(!is.null(region)){ 
    tp <- length(HCP[region,,patient])
    Dataplot <- as.data.table(cbind(HCP[region,,patient],1:tp))
  } else {
    tp <- length(HCP)
    Dataplot <- as.data.table(cbind(HCP,1:tp))
  }
  setnames(Dataplot, c("BOLD","Time"))
  plot1 <- ggplot(data = Dataplot) + geom_line(aes(x=Time, y= BOLD),size=1)
  if(save){
    ggsave(filename = paste("./Plots/PlotTimeSeries",name,".png", sep = ""),
           plot = plot1,
           width = 12,
           height = 8,
           units = "in") 
  }
  print(plot1)
  return(plot1)
}


PlotTimeSeriesMovieWatching <- function(HCP,region=1,patient=1,save=T,name=""){

  #Plots the time series of the movie watching dataset (inputed in HCP) for the selected patient and region.
  #Colors the periods
  
  if(!is.null(region)){ 
    tp <- length(HCP[region,,patient])
    Dataplot <- as.data.table(cbind(HCP[region,,patient],1:tp))
  } else {
    tp <- length(HCP)
    Dataplot <- as.data.table(cbind(HCP,1:tp))
  }
  setnames(Dataplot, c("BOLD","Time"))
  Dataplot[,Type := "Transition"]
  Dataplot[Time %in% MovieData$TimePointsMovie,Type := "Movie"]
  Dataplot[Time %in% MovieData$TimePointsRest,Type := "Rest"]
  Dataplot[,group:=1]
  plot1 <- ggplot(data = Dataplot) + geom_line(aes(x=Time,y=BOLD,color=Type,group=group),size=1) +
    scale_color_manual(values = c(brewer.pal(3,"Set1")[1:2],"black")) + theme_minimal()
  if(save){
    ggsave(filename = paste("./Plots/PlotMovieTimeSeries",name,".png", sep = ""),
           plot = plot1,
           width = 12,
           height = 8,
           units = "in") 
  }
  print(plot1)
  return(plot1)
}

PlotBothTimeSeries <- function(HCP,HCPbin,region=NULL,patient=NULL,name1="Normalized",name2="Thresholded",save=T){

  # Plots two time series at the same time. Input is relaxed. If patient and region are selected, the whole dataset can be provided in HCP and HCPbin,
  # and it will plot for the same region and patient in both cases
  # if not, then only provide the corresponding vectors in both.

  if(!is.null(region)){ 
  Dataplot1 <- as.data.table(cbind(HCP[region,,patient],1:1190))
  Dataplot2 <- as.data.table(cbind(HCPbin[region,,patient],1:1190))
  } else {
  Dataplot1 <- as.data.table(cbind(HCP,1:1190))
  Dataplot2 <- as.data.table(cbind(HCPbin,1:1190))
    
  }
  setnames(Dataplot1, c("BOLD","Time"))
  setnames(Dataplot2, c("BOLD","Time"))
  Dataplot1[, signal := name1]
  Dataplot2[, signal := name2]
  Dataplot <- rbind(Dataplot1,Dataplot2)
  plot1 <- ggplot(data = Dataplot) + geom_line(aes(x=Time, y= BOLD, color = signal), size = 1) + 
    scale_colour_manual(values = c('black','#d7191c'))
  if(save){
   ggsave(filename = paste("./Plots/PlotBothTimeSeries_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 12,
         height = 8,
         units = "in") 
  }
  print(plot1)
}

PlotMeanConnectivityMeasures <- function(){

  # convenient wrapper to plot mean accordance, mean discordance and mean pearson for the binarized dataset

  MeanAccordance <- apply(Concordance$acc,c(1,2),mean)
  PlotCorr(MeanAccordance,name="MeanAccordance")
  MeanDiscordance <- apply(Concordance$dis,c(1,2),mean)
  PlotCorr(abs(MeanDiscordance),name="MeanDiscordance")
  
  PlotCorr(MeanPearson(HCP.bin),name="MeanPearson")
}



PlotConcordance <- function(Conc,patient=1,name){

  # Plot the concordance for the selected patient

  PlotCorr(Conc$acc[,,patient],name=paste("Accordance",patient,name,sep=""))
  PlotCorr(abs(Conc$dis[,,patient]),name=paste("Discordance",patient,name,sep=""))

}

PlotPartialCorr <- function(HCP,patient,name=""){

  #Plots the partial correlation, using the generalized inverse method

  data <- pcor(t(HCP[,,patient]),method='pearson')$estimate
  data <- ifelse(data >= 0.99,NA,data)
  PlotCorr(data,paste('partial_',name,sep=""))
}


PlotCorr <- function(CorrMatrix,name="",save=T){

  # general purpose plotting function to plot any connectivity metric, stored in the matrix CorrMatrix

  dt <- as.data.table(CorrMatrix)
  m <- dim(CorrMatrix)[1]
  setnames(dt,as.character(paste(1:m,sep="")))
  dt[,RegionA:= as.character(.I)]
  dt <- melt(dt, id.vars='RegionA',value.name = "Corr",variable.name = 'RegionB')
  dt[,RegionA := factor(RegionA, levels = paste(c(1:m),sep=""))]
  dt[,RegionB := factor(RegionB, levels = paste(c(m:1),sep=""))]
  
  plot1 <- ggplot(data=dt) + geom_tile(aes(x=RegionA,y=RegionB,fill=Corr),color='grey') + 
    scale_fill_gradient2(low="#E41A1C", mid = 'white', high = "#377EB8") + labs(x='Region A',y='Region B') + 
    theme_classic() + 
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
  if(save){
    ggsave(filename = paste("./Plots/CorrPlot_",name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
           plot = plot1,
           width = 8,
           height = 8,
           units = "in")
  }
  print(plot1)
}




#### EXECUTION EXAMPLES ####

##Set up Data to use every functions
 init(workingdir = "") #input here the working directory
 HCP <- LoadData()
 MovieData <- LoadAdditionalData_MovieWatching()
 MCIData <- LoadAdditionalData_MCI()
 load("./Data/MarkovCoeffMCI.rda")
 HCP.normal <- NormalizeDataset(HCP)
 HCP.bin <- BinarizeDataset(HCP.normal,0.84,0.16)
 Concordance <- ComputeConcordance(HCP.bin)

##Compute the partial correlation using the residuals for patient 1
#p.cor <- ComputePartialCorrelation(HCP.bin,1)

## Some Easy Plots
#PlotBothTimeSeries(HCP.normal,HCP.bin,1,1)
#PlotBothTimeSeries(HCP.normal,HCP.Thresh,34,65)
#PlotConcordance(Concordance,53)
#PlotCorr(cor(t(HCP.Thresh[,,53])),name='BaseThresh53')
#PlotCorr(cor(t(HCP.normal[,,53])),name='normal53')
#DistanceCorrelationPlot(HCP)
#PlotMeanConnectivityMeasures()

