######################################################
#
# File containing the function definition pertaining to
# the markov chain connection analysis using 6 states ! 
# Contains data loading, Manipulation, plotting and CTMC fitting
#
#
# These functions were used to study the 9 states case, but were
# not used in the clustering approach. All MC coef are computed
#
#######################################################


####Function Def####

CreateConnectionTimeSeries <- function(){

  # create the array containing the binarized time series
  # unused in the project, but discussed as a possible next step

  pichu <- array(0,dim=c(90,90,1190,90))
   for(p in 1:90){
     for(j in 1:90){
       for(i in 1:90){
         cat(i,j,p,"\n")
     pichu[i,j,,p] <- 1*((HCP.bin[i,,p]==1 & HCP.bin[j,,p]==1) | (HCP.bin[i,,p]==-1 & HCP.bin[j,,p]==-1)) - 
                     1*((HCP.bin[i,,p]==1 & HCP.bin[j,,p]==-1) | (HCP.bin[i,,p]==-1 & HCP.bin[j,,p]==1))
       }
     }
   }
}

MapState <- function(a,b){

  # maps the value of two time series to the proper state of the markov chain
  # this version contains all 9 states! ud and du are considered different 


  dt<-as.data.table(cbind(a,b))
  setnames(dt,c("a","b"))
  dt[a==1 & b==1,d:="uu"]
  dt[a==1 & b==0, d:="uz"]
  dt[a==1 & b==-1, d:="ud"]
  dt[a==0 & b==1, d:="zu"]
  dt[a==0 & b==0, d:="zz"]
  dt[a==0 & b==-1, d:="zd"]
  dt[a==-1 & b==1, d:="du"]
  dt[a==-1 & b==0, d:="dz"]
  dt[a==-1 & b==-1, d:="dd"]
  return(dt[,.(d)])
}

CreateMarkovChainList_Conn <- function(TimeSeries1,TimeSeries2){

  # create the list of states along with the time spent to fit a Markov chain from the two time series
   
  datatest <- MapState(TimeSeries1,TimeSeries2)
  setnames(datatest,"V1")
  datatest[, count := .I]
  datatest[, mincount := min(count)-1,.(rleid(V1))]
  datatest <- unique(datatest[,.(V1,mincount)])
  
  listtest <- list(datatest$V1,datatest$mincount)
  
  return(listtest)
}

EstimateMarkovChain_Conn <- function(listtest){

  #fits one COMPLETE markov chain, including transition probabilities using the list. 
  #Coefficients are computed using the MLE.
  
  #setting up df.
  df.test <- as.data.frame(listtest)
  setnames(df.test,c("state","time"))
  statename <- c("uu","uz","ud","zu","zz","zd","du","dz","dd")
  nstate <- 9
  TransitionMatrix <- matrix(0,nrow=nstate,ncol=nstate)
  
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

  dimnames(TransitionMatrix) <-list(statename,statename)
  for(i in statename){
    TransitionMatrix[i,i] <- ifelse(length(df.test.exp[Statebegin==i]$lambda)==0,NA,df.test.exp[Statebegin==i]$lambda)
    for(j in statename){
      if(i!=j & nrow(df.test[Statebegin==i & Stateend ==j])>0){
        TransitionMatrix[i,j] <- ifelse(length(df.test[Statebegin==i & Stateend ==j]$ProbaTrans)==0,0,df.test[Statebegin==i & Stateend ==j]$ProbaTrans)
      }
    }
  }
  return(TransitionMatrix)
}

EstimateSeveralMarkovChains_Conn <- function(region1 = 1,region2=1:90,patient1=1,mean=F){
  
  #fits all coefficients of several markov chains for a chosen number of patients and regions
  # note : this function is very slow to compute the coefficients for the whole and complete dataset
  #Coefficients are computed using the MLE.
  # if mean is set to TRUE, returns the mean coefficients across patient. Otherwise returns the full array.
  
  coeff <- array(0,dim=c(length(region1),length(region2),9,9,length(patient1)))
  setkey(Nodes,index)
  names <- Nodes$name
  names1 <- names[region1]
  names2 <- names[region2]
  statename <- c("uu","uz","ud","zu","zz","zd","du","dz","dd")
  dimnames(coeff)<-list(names1,names2,statename,statename,NULL)
  
  for(j in 1:length(patient1)){
    for(i in 1:length(region1)){
      for(ii in 1:length(region2)){
        listtest <- CreateMarkovChainList_Conn(HCP.bin[region1[i],,patient1[j]],HCP.bin[region2[ii],,patient1[j]])
        EasyMatrix <- EstimateMarkovChain_Conn(listtest)
        coeff[i,ii,,,j] <- as.array(EasyMatrix)
      }
    }
  }
  if(mean)  coeff <- apply(coeff,c(1,2,3,4),function(x) mean(x,na.rm=T))
  return(coeff)
}

DotPlotClusterLambda_Conn<- function(region1 = 1,region2=1:90,patient1=1){
  
  #Plots the dotplots of the lambda coefficient by fixing one region of the pair.
  # the coefficient are normalized by the coefficient of the markov chain of region1.
  # this function compute the coefficient again for the selected region & patient combination.

  # I - setting up the dataframe
  coef <- EstimateSeveralMarkovChains_Conn(region1,region2,patient1,F)
  matrix1 <- EstimateMarkovChain(CreateMarkovChainList(HCP.bin[region1,,patient1]))
  plot.df <- data.table()
  statename <- c("uu","uz","ud","zu","zz","zd","du","dz","dd")
  for(i in 1:90){
    matrix2 <- EstimateMarkovChain(CreateMarkovChainList(HCP.bin[i,,patient1]))
    plot.df.temp <- cbind(rep(Nodes[index==i]$name,9),rep(Nodes[index==i]$cluster,9))
    coef.vector <- diag(coef[1,i,,,1])
    normalizing.vector <- c(max(matrix1["1","1"],matrix2["1","1"]),
                            max(matrix1["1","1"],matrix2["0","0"]),
                            max(matrix1["1","1"],matrix2["-1","-1"]),
                            max(matrix1["0","0"],matrix2["1","1"]),
                            max(matrix1["0","0"],matrix2["0","0"]),
                            max(matrix1["0","0"],matrix2["-1","-1"]),
                            max(matrix1["-1","-1"],matrix2["1","1"]),
                            max(matrix1["-1","-1"],matrix2["0","0"]),
                            max(matrix1["-1","-1"],matrix2["-1","-1"]))
    plot.df.temp <- cbind(plot.df.temp,coef.vector,normalizing.vector,statename)
    plot.df <- rbind(plot.df,plot.df.temp)
  }
  
  setnames(plot.df,c("Region","Cluster","Coefc","Normal","CoefName"))
  plot.df[,CoefValue := as.numeric(Coefc)]
  plot.df[,NormalValue := as.numeric(Normal)]
  #plot.df <- plot.df[CoefName %in% c("Lambda_-1","p_01","Lambda_0","Lambda_1")]
  
  # II - plotting
  plot1 <- ggplot(data=plot.df) + geom_dotplot(aes(x=1,y=CoefValue/NormalValue,fill=as.factor(Cluster)),stackgroups=T,
                                               binpositions="all",binaxis = "y",stackdir="center",dotsize = 1) + 
    scale_fill_brewer(palette="Set1", name = "Cluster") + facet_wrap(~CoefName,scales="free") + xlab('') + 
    ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
  ggsave(filename = paste("./Plots/dotPlot_MarkovCoefCluster_Conn",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 9,
         height = 8,
         units = "in")
  print(plot1)
  return(plot.df)
}

DotPlotClusterLambdaCut_Conn<- function(region1 = 1,region2=1:90,patient1=1,ncut=4,mean=F){
  
  #Plots the dotplots and the line plots of the lambda coefficient over sections (time) by fixing one region of the pair, faceted by clusters.
  # the coefficient are NOT normalized
  # this function compute the coefficient again for the selected region & patient combination.

  # I - setting up the dataframe
  patient.name=patient1
  if(length(patient1)>1) patient.name="All"
  coef <- EstimateSeveralMarkovChainsCut_Conn(region1,region2,patient1,ncut,mean)
  #matrix1 <- MarkovCut(HCP.bin[region1,,patient1],ncut,F)
  plot.df <- data.table()
  statename <- c("uu","uz","ud","zu","zz","zd","du","dz","dd")
  for(i in 1:90){
    #matrix2 <- MarkovCut(HCP.bin[i,,patient1],ncut,F)
    for(j in 1:ncut){
    plot.df.temp <- cbind(rep(Nodes[index==i]$name,9),rep(Nodes[index==i]$cluster,9),rep(j,9))
    if(mean) {
      coef.vector <- diag(coef[1,i,,,j])
    } else {
      coef.vector <- diag(coef[1,i,,,1,j])
    }
    
    normalizing.vector <- rep(1,9)
    plot.df.temp <- cbind(plot.df.temp,coef.vector,normalizing.vector,statename)
    plot.df <- rbind(plot.df,plot.df.temp)
    }
  }
  
  setnames(plot.df,c("Region","Cluster","Section","Coefc","Normal","CoefName"))
  plot.df[,CoefValue := as.numeric(Coefc)]
  plot.df[,NormalValue := as.numeric(Normal)]
  
  # II - plotting
  for(coef.name in statename){
    plot.df.temp <- plot.df[CoefName == coef.name]
    plot1 <- ggplot(data=plot.df.temp) + geom_dotplot(aes(x=1,y=CoefValue/NormalValue,fill=as.factor(Cluster)),stackgroups=T,
                                               binpositions="all",binaxis = "y",stackdir="center",dotsize = 1) + 
    scale_fill_brewer(palette="Set1", name = "Cluster") + facet_grid(~as.factor(Section),scales="free") + xlab('') + 
    ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
    ggsave(filename = paste("./Plots/dotPlot_MarkovCoefClusterCut_Conn_",region1,"_",patient.name,"_",coef.name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
           plot = plot1,
           width = 9,
           height = 8,
           units = "in")
    print(plot1)
    
    plot2 <- ggplot(data=plot.df.temp) + geom_line(aes(x=Section,y=CoefValue/NormalValue,color=as.factor(Cluster),group=Region)) + 
      scale_color_brewer(palette="Set1", name = "Cluster") + xlab('Section') + 
      ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    ggsave(filename = paste("./Plots/linePlot_MarkovCoefClusterCut_Conn_",region1,"_",patient.name,"_",coef.name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
           plot = plot2,
           width = 9,
           height = 8,
           units = "in")
    print(plot2)
  }
  return(plot.df)
}



FitExp_Conn <- function(region1,region2,patient,state){

    #Test if the time spent is actually exponential distributed


  test <- CreateMarkovChainList_Conn(HCP.bin[region1,,patient],HCP.bin[region2,,patient])
  matrixt <- EstimateMarkovChain_Conn(test)
  df.test <- as.data.frame(test)
  setnames(df.test,c("state","time"))
  df.test <- as.data.table(cbind(df.test[-nrow(df.test),],df.test[-1,]))
  setnames(df.test,c("Statebegin","time1","Stateend","time2"))
  df.test[,TimeSpent := time2-time1]
  df.test <- df.test[Statebegin == state,.(TimeSpent)]
  setkey(df.test,"TimeSpent")
  
  n <- nrow(df.test)
  if(n<=1){
    cat("not enough state in the time series !")
    return(FALSE)
  }
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
CompareCoeffWithSingleRegion <- function(ts1,ts2){

   # test function to compare the lambda coeff of the pair to the one of each of its componenents.


  test <- CreateMarkovChainList_Conn(ts1,ts2)
  matrixt <- EstimateMarkovChain_Conn(test)
  matrix1 <- EstimateMarkovChain(CreateMarkovChainList(ts1))
  matrix2 <- EstimateMarkovChain(CreateMarkovChainList(ts2))
  cat("uu/lambda_1^1" ,matrixt["uu","uu"]/matrix1["1","1"],
      " uu/lambda_1^2 ",matrixt["uu","uu"]/matrix2["1","1"],"\n")
  
  cat("dd/lambda_-1^1" ,matrixt["dd","dd"]/matrix1["-1","-1"],
      " dd/lambda_-1^2 ",matrixt["dd","dd"]/matrix2["-1","-1"],"\n")
  
  cat("zz/lambda_0^1" ,matrixt["zz","zz"]/matrix1["0","0"],
      " zz/lambda_0^2 ",matrixt["zz","zz"]/matrix2["0","0"],"\n")

}

PlotCoefficientTest <- function(region1,coef.name="zu"){


    # test plotting function. Largely unused.

  setkey(Nodes,"index")
  coef <- EstimateSeveralMarkovChains_Conn(region1,1:90,1,F)
  matrix1 <-EstimateMarkovChain(CreateMarkovChainList(HCP.bin[region1,,1]))
  matrix2 <-EstimateSeveralMarkovChains(regions=1:90,patient=1,mean=F)
  matrix2 <- matrix2[3,3,,1]
  #create plot data frame
  plot.df <- as.data.table(cbind(coef[1,,"uu","uu",1],coef[1,,"uz","uz",1],coef[1,,"ud","ud",1],
                   coef[1,,"du","du",1],coef[1,,"zu","zu",1],Nodes$cluster))
  setnames(plot.df,c("uu","uz","ud","du","zu","cluster"))
  plot.df[,uu:=uu/matrix1["1","1"]]
  plot.df[,zu:=zu/matrix2]
  # 
  # plot.df <- as.data.table(cbind(coef[1,,"zz","zz",1],coef[1,,"uz","uz",1],coef[1,,"ud","ud",1],
  #                                coef[1,,"du","du",1],coef[1,,"zu","zu",1],Nodes$cluster))
  # setnames(plot.df,c("zz","uz","ud","du","zu","cluster"))
  # plot.df[,zz:=zz/matrix1["0","0"]]
  
  plot1 <- ggplot(data=plot.df) + geom_point(aes(x=uu,y=get(coef.name),color=as.factor(cluster))) +
    geom_smooth(aes(x=uu,y=get(coef.name))) + scale_color_brewer(palette="Set1")
    # geom_point(aes(x=uu,y=ud),color=green) + geom_line(aes(x=uu,y=uz),color=red)
  print(plot1)

}

MarkovCut_Conn <- function(ts1,ts2,numberofcut){

  #Fits one markov chain for the stationarity analysis. The two time series are cut $numberofcut$ times
  # and a complete CTMC is fitted for each of the sections

  #generate the list containt the n cut time series
  len <- floor(length(ts1)/numberofcut)
  ts1_c <- list()
  ts2_c <- list()
  for(i in 0:(numberofcut-2)){
    ts1_c[[i+1]] <- ts1[(1+i*len):((i+1)*len)]
    ts2_c[[i+1]] <- ts2[(1+i*len):((i+1)*len)]
  }
  ts1_c[[numberofcut]] <- ts1[((numberofcut-1)*len):length(ts1)]
  ts2_c[[numberofcut]] <- ts2[((numberofcut-1)*len):length(ts2)]
  #generates the list containing the markov chain coeff
  markovcut <- array(0,dim=c(9,9,numberofcut))
  for(i in 1:numberofcut){
    markovcut[,,i] <- as.array(EstimateMarkovChain_Conn(CreateMarkovChainList_Conn(ts1_c[[i]],ts2_c[[i]])))
  }
  return(markovcut)
}

EstimateSeveralMarkovChainsCut_Conn <- function(region1 = 1,region2=1:90,patient1=1,numberofcut,mean=F,coeff){


  #Fits several markov chains for the stationarity analysis. The time series are cut $numberofcut$ times
  # and a CTMC are fitted for each of the sections
  # if mean is set to TRUE, returns the mean coefficients across patient. Otherwise returns the full array.
  # also saves the data as .rda for easy retrieving

    for(j in 1:length(patient1)){
      cat("computing patient",j,"at",format(Sys.time(),"%H:%M"))
      for(i in 1:length(region1)){
        for(ii in 1:length(region2)){
          cat(".")
          EasyMatrix <- MarkovCut_Conn(HCP.bin[region1[i],,patient1[j]],HCP.bin[region2[ii],,patient1[j]],numberofcut)
          coeff[i,ii,,,j,] <- as.array(EasyMatrix)
        }
      }
      save(coeff, file = 'MarkovChain4CutCoeffConnections.rda')
      cat("\n")
    }
    if(mean)  coeff <- apply(coeff,c(1,2,3,4,6),function(x) mean(x,na.rm=T))
    return(coeff)
}