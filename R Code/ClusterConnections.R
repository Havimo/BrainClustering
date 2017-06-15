######################################################
#
# File containing the function definition pertaining to
# the two dimensional cluster analysis. Contains data loading,
# Manipulation, plotting and cluster fitting.
#
#######################################################


ClusterOnePatient <- function(array,indexmap,nc=9,standalone=F){

	# cluster the pairs for one patient.
	# if standalone is selected returns the complete data frame
	# if not, only returns the vector of clusters

  lambdavector.matrix<-matrix(0,nrow=nrow(indexmap), ncol=6*4)
  cat("creating distance matrix at ",format(Sys.time(),"%H:%M"))
  for(i in 1:nrow(indexmap)){
    ii <- as.numeric(indexmap[index==i,.(lhs)])
    jj <- as.numeric(indexmap[index==i,.(rhs)])
    lambdavector.matrix[i,] <- c(array[ii,jj,,1],array[ii,jj,,2],array[ii,jj,,3],array[ii,jj,,4])
  }
  cat("\n")
  dimnames(lambdavector.matrix) <- list(indexmap$name,NULL)
  cat("computing distance matrix at",format(Sys.time(),"%H:%M"))
  distance.matrix <- dist(lambdavector.matrix)
  cat("\n")
  tree <- hclust(d= distance.matrix,method="complete")
  cluster <- cutree(tree,k=nc)
  if(standalone){
    clusterf <- merge(indexmap,data.table(cluster,row.names(as.data.frame(cluster))),by.x="name",by.y="V2")
    clusterf[,cluster.rhs := Nodes[name==rhsname,.(cluster)],.(rhs,index)]
    clusterf[,cluster.lhs := Nodes[name==lhsname,.(cluster)],.(lhs,index)]
  } else {
    clusterf <- cluster
  }
  return(clusterf)
}

CreateRegionIndexMap <- function(){

  #this function creates the linear mapping to the pair index from the lexico graphic order of the pairs 
  #note that due to algorithm running time and size of the dataset, we only consider 4005 pairs.
  
   lhs <- rep(1:90,each = 90)
   rhs <- rep(1:90,90)
   test <- as.data.table(cbind(lhs,rhs))
   setnames(test,c("lhs","rhs"))
   test <- test[lhs!=rhs]
   test[,lhsname := Nodes[index==lhs]$name,.(lhs)]
   test[,rhsname := Nodes[index==rhs]$name,.(rhs)]
   test[,name := paste("(",lhsname,",",rhsname,")",sep="")]
   test[lhs<rhs,orderedpair := paste(lhs,rhs,sep="-")]
   test[rhs<lhs,orderedpair := paste(rhs,lhs,sep="-")]
   setkey(test,lhs)
   test[, index2 := rowid(orderedpair)]
   test <- test[index2==1]
   setkey(test,lhs)
   test[,index := .I]
   setkey(test,index)
   return(test)
}

CreateClusterMatrix_Conn <- function(cluster,name="",save=T,plot=T){

	# create the cluster matrix using a vector product instead of a for loop for reasons of speed

  cat("Creating Cluster Matrix at",format(Sys.time(),"%H:%M"))
  clustermatrix <- matrix(0,nrow=length(cluster),ncol=length(cluster))
  for(i in 1:(length(unique(cluster)))){
    cluster_temp <- ifelse(cluster==i,1,0)
    temp_matrix <- cluster_temp %*% t(cluster_temp)
    clustermatrix = clustermatrix + temp_matrix
  }
  if(plot) PlotCluster(clustermatrix,name=name,save=save)
  cat("\n")
  return(clustermatrix)
}

CreateBootsrapConnectivity_Conn <- function(array,indexmap,nc=9){
	
	# computes the bootstrap connectivity metric for the two dimensional case
	
  npatient <- 90
  ClusterMatrix <- matrix(0,nrow=nrow(indexmap),ncol=nrow(indexmap))
  
  for(i in 1:npatient){
    cat(paste("clustering for patient",i,"\n"))
    time1<- Sys.time()
    cluster <- ClusterOnePatient(array[,,,i,],indexmap,nc,F)
    temp <- CreateClusterMatrix_Conn(cluster,save=F,plot=F)
    temp[temp>0] <- 1
    temp[is.na(temp)] <- 0
    ClusterMatrix <- ClusterMatrix + temp
    cat("patient done in ")
    print(Sys.time()-time1)
    cat("\n")
  }
  
  dimnames(ClusterMatrix) <-list(indexmap$name,indexmap$name)
  #PlotCorr(ClusterMatrix/(npatient),name=paste("Cluster7tree",method,sep=""))
  return(ClusterMatrix/(npatient))
  
}

FinalBoostrapClustering_Conn <- function(Bootstrap.Conn,indexmap,Nodes,name=""){
	
	# computes the final cluster distribution for the two dimensional case

  tree <- hclust(d=as.dist(1-Bootstrap.Conn),method="ward.D2")
  cluster <- cutree(tree,k=9)
  clusterf <- merge(indexmap,data.table(cluster,row.names(as.data.frame(cluster))),by.x="name",by.y="V2")
  clusterf[,cluster.rhs := Nodes[name==rhsname,.(cluster)],.(rhs,index)]
  clusterf[,cluster.lhs := Nodes[name==lhsname,.(cluster)],.(lhs,index)]
  return(clusterf)
  
}


DotPlotClusterLambdaCutCluster_Conn<- function(cluster,coeff,patient1=1,mean = F){

	# plots the markov chain coefficients across the 4 sectoins
	# only line plot is not commented as this the most interesting one 
	# but box plots and dot plots are also implemented
  # the inputs are as follows
  #   -  cluster is a dataframe containing the cluster distribution and the nodes metadata (ouput of FinalBoostrapClustering_Conn())
  #   -  coeff is the dataframe containing all the coefficients (loaded from load("MarkovChainLambdas4CutCoeffConnections.rda") )
  # if mean is set to TRUE, then plots the mean coefficients. Otherwise plots for the selected patient.  
  
  
  # I - setting up the dataframe
  statename <- c("uu","uz","ud","zz","zd","dd")
  plot.df <- data.table()
  if(mean){
    coeff <- apply(coeff,c(1,2,3,4,6),function(x) mean(x,na.rm=T))
    for(j in 1:nrow(cluster)){
      ii <- as.numeric(indexmap[index==j,.(lhs)])
      jj <- as.numeric(indexmap[index==j,.(rhs)])
      plot.df.t <- cbind(rep(cluster[index==j]$name,6*4),rep(cluster[index==j]$cluster,6*4),rep(1:4,each=6),rep(statename,4))
      plot.df.t <- as.data.table(cbind(plot.df.t,c(coeff[ii,jj,,patient1,1],coeff[ii,jj,,patient1,2],
                                                   coeff[ii,jj,,patient1,3],coeff[ii,jj,,patient1,4])))
      plot.df <- rbind(plot.df,plot.df.t)
    }

  } else {
    for(j in 1:nrow(cluster)){
      ii <- as.numeric(indexmap[index==j,.(lhs)])
      jj <- as.numeric(indexmap[index==j,.(rhs)])
      plot.df.t <- cbind(rep(cluster[index==j]$name,6*4),rep(cluster[index==j]$cluster,6*4),rep(1:4,each=6),rep(statename,4))
      plot.df.t <- as.data.table(cbind(plot.df.t,c(coeff[ii,jj,,patient1,1],coeff[ii,jj,,patient1,2],
                                                   coeff[ii,jj,,patient1,3],coeff[ii,jj,,patient1,4])))
      plot.df <- rbind(plot.df,plot.df.t)
    }
  }
  
  setnames(plot.df,c("Pair","Cluster","Section","CoefName","Coefc"))
  plot.df[,CoefValue := as.numeric(Coefc)]
  
  # II - plotting
  for(coef.name in statename){
    plot.df.temp <- plot.df[CoefName == coef.name]
    # plot1 <- ggplot(data=plot.df.temp) + geom_dotplot(aes(x=1,y=CoefValue,fill=as.factor(Cluster)),stackgroups=T,
    #                                                   binpositions="all",binwidth=100,binaxis = "y",stackdir="center",dotsize = 1) + 
    #   scale_fill_brewer(palette="Set1", name = "Cluster") + facet_grid(~as.factor(Section),scales="free") + xlab('') + 
    #   ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    # 
    # ggsave(filename = paste("./Plots/dotPlot_MarkovCoefClusterCut_Conn_",coef.name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
    #        plot = plot1,
    #        width = 9,
    #        height = 8,
    #        units = "in")
    # print(plot1)
    
    plot2 <- ggplot(data=plot.df.temp) + geom_line(aes(x=Section,y=CoefValue,color=as.factor(Cluster),group=Pair)) + 
      scale_color_brewer(palette="Set1", name = "Cluster") + xlab('Section') + facet_wrap(~Cluster)
    ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    ggsave(filename = paste("./Plots/linePlot_MarkovCoefClusterCut_Conn_",coef.name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
           plot = plot2,
           width = 9,
           height = 8,
           units = "in")
    print(plot2)
    
    #plot3 <- ggplot(data=plot.df.temp) + geom_boxplot(aes(x=Section,y=CoefValue,fill=as.factor(Cluster))) + 
    #  scale_fill_brewer(palette="Set1", name = "Cluster") + facet_grid(~Cluster) + xlab('') + 
    #  ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    #
    #ggsave(filename = paste("./Plots/boxPlot_MarkovCoefClusterCut_Conn_",coef.name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
    #       plot = plot3,
    #       width = 9,
    #       height = 8,
    #       units = "in")
    #print(plot3)
  }
  return(plot.df)
}

LinePlotClusterLambdaCutClusterRegion_Conn<- function(cluster,coeff,patient=1,region=1,mean=F){

	# plots the box plot for the cut into 4 markov chain coefficients 
  # the inputs are as follows
  #   -  cluster is a dataframe containing the cluster distribution and the nodes metadata (ouput of FinalBoostrapClustering_Conn())
  #   -  coeff is the dataframe containing all the coefficients (loaded from load("MarkovChainLambdas4CutCoeffConnections.rda") )
  # if mean is set to TRUE, then plots the mean coefficients. Otherwise plots for the selected patient.  
  
  # I - setting up the dataframe
  statename <- c("uu","uz","ud","zz","zd","dd")
  plot.df <- data.table()
  if(mean==F){
    for(j in 1:nrow(cluster)){
      ii <- as.numeric(indexmap[index==j,.(lhs)])
      jj <- as.numeric(indexmap[index==j,.(rhs)])
      plot.df.t <- cbind(rep(ii,6*4),rep(jj,6*4),rep(cluster[index==j]$name,6*4),rep(cluster[index==j]$cluster,6*4),rep(1:4,each=6),rep(statename,4))
      plot.df.t <- as.data.table(cbind(plot.df.t,c(coeff[ii,jj,,patient,1],coeff[ii,jj,,patient,2],
                                                   coeff[ii,jj,,patient,3],coeff[ii,jj,,patient,4])))
      plot.df <- rbind(plot.df,plot.df.t)
    }
  } else {
    coeff <- apply(coeff,c(1,2,3,4,6),function(x) mean(x,na.rm=T))
    for(j in 1:nrow(cluster)){
      ii <- as.numeric(indexmap[index==j,.(lhs)])
      jj <- as.numeric(indexmap[index==j,.(rhs)])
      plot.df.t <- cbind(rep(ii,6*4),rep(jj,6*4),rep(cluster[index==j]$name,6*4),rep(cluster[index==j]$cluster,6*4),rep(1:4,each=6),rep(statename,4))
      plot.df.t <- as.data.table(cbind(plot.df.t,c(coeff[ii,jj,,1],coeff[ii,jj,,2],
                                                   coeff[ii,jj,,3],coeff[ii,jj,,4])))
      plot.df <- rbind(plot.df,plot.df.t)
    }
  }
  
  setnames(plot.df,c("lhs","rhs","Pair","Cluster","Section","CoefName","Coefc"))
  plot.df[,CoefValue := as.numeric(Coefc)]
  plot.df <-plot.df[rhs==region | lhs ==region]
  # II - plotting
  for(coef.name in statename){
    plot.df.temp <- plot.df[CoefName == coef.name]
    plot2 <- ggplot(data=plot.df.temp) + geom_line(aes(x=Section,y=CoefValue,color=as.factor(Cluster),group=Pair)) +
      scale_color_brewer(palette="Set1", name = "Cluster") + xlab('Section') +
      ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    ggsave(filename = paste("./Plots/linePlot_MarkovCoefClusterCutRegion_Conn_",coef.name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
           plot = plot2,
           width = 9,
           height = 8,
           units = "in")
    print(plot2)
  }
  return(plot.df)
}

BoxPlotClusterLambdaCutCluster_Conn<- function(cluster,coeff,patient=1:90){

	# plots the box plot for the cut into 4 markov chain coefficients 
  # the inputs are as follows
  #   -  cluster is a dataframe containing the cluster distribution and the nodes metadata (ouput of FinalBoostrapClustering_Conn())
  #   -  coeff is the dataframe containing all the coefficients (loaded from load("MarkovChainLambdas4CutCoeffConnections.rda") )
  
  # I - setting up the dataframe
  statename <- c("uu","uz","ud","zz","zd","dd")
  plot.df <- data.table()
  for(p in patient){
    cat(p)
    for(j in 1:nrow(cluster)){
      cat(".")
      ii <- as.numeric(indexmap[index==j,.(lhs)])
      jj <- as.numeric(indexmap[index==j,.(rhs)])
      plot.df.t <- cbind(rep(p,6*4),rep(ii,6*4),rep(jj,6*4),rep(cluster[index==j]$name,6*4),rep(cluster[index==j]$cluster,6*4),rep(1:4,each=6),rep(statename,4))
      plot.df.t <- as.data.table(cbind(plot.df.t,c(coeff[ii,jj,,p,1],coeff[ii,jj,,p,2],
                                                   coeff[ii,jj,,p,3],coeff[ii,jj,,p,4])))
      plot.df <- rbind(plot.df,plot.df.t)
    }
    cat("\n")
  }
  
  setnames(plot.df,c("patient","lhs","rhs","Pair","Cluster","Section","CoefName","Coefc"))
  plot.df[,CoefValue := as.numeric(Coefc)]
  # II - plotting
  for(coef.name in statename){
    plot.df.temp <- plot.df[CoefName == coef.name]
    plot3 <- ggplot(data=plot.df.temp) + geom_boxplot(aes(x=Section,y=CoefValue,fill=as.factor(Cluster))) + 
      scale_fill_brewer(palette="Set1", name = "Cluster") + facet_grid(~Cluster) + xlab('') + 
      ylab("Coef. Value") +  theme_classic() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    ggsave(filename = paste("./Plots/boxPlot_MarkovCoefClusterCut_Conn_",coef.name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
           plot = plot3,
           width = 9,
           height = 8,
           units = "in")
    print(plot3)
  }
  return(plot.df)
}

####Execution####

##Preliminary Loading
#indexmap <- CreateRegionIndexMap()
#load("MarkovChainLambdas4CutCoeffConnections.rda")

##Cluster Patient 1
#test <- ClusterOnePatient(coeff[,,,1,],indexmap,standalone = T)

##Create the Booststrap Connectivity for the Connections for 9 clusters
#Bootstrap.Conn <- CreateBootsrapConnectivity_Conn(coeff,indexmap,9)
#ClusterConn <- FinalBoostrapClustering_Conn(Bootstrap.Conn,indexmap)
#DotPlotClusterLambdaCutCluster_Conn(ClusterConn,coeff,mean=T)
