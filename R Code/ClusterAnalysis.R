######################################################
#
# File containing the function definition pertaining to
# the one dimensional cluster analysis. Contains data loading,
# Manipulation, plotting and cluster fitting
#
#######################################################

####Function Def####

#Prior data handling and loading

AverageData <- function(HCP,dimension = "region"){

	# this function averages the dataset across patients or regions
  # not really used furter on

  dataset <- HCP
  if(dimension == "region"){
    dataset <- apply(HCP,c(1,2),mean)
  } else if(dimension == "patient") {
    dataset <- apply(HCP,c(2,3),mean)
  } 
if(dim(dataset)[1]==1190) dataset <- t(dataset)
return(dataset)
}

ShoulderCriterion <- function(data, nc=15, seed=1234){

  #creates and plots the shoulder criterion graph used to graphically determine the number of clusters.
  # data should be a matrix or 2d array, hence either the mean patient or one patient

  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)
  }
  plot1 <- ggplot(data=as.data.table(cbind(1:nc,wss)),aes(x=V1, y= wss)) + geom_point() + geom_line() + xlab("Number of Clusters") +
       geom_segment(aes(x=8,xend=8,y=1.8e7,yend=4.5e7),color="red") + ylab("Within groups sum of squares") + theme_classic()
  ggsave(filename = paste("./Plots/shouldercrit.png", sep = ""),
         plot = plot1,
         width = 8,
         height = 8,
         units = "in") 
  print(plot1)
}


CreateClusterMatrix <- function(cluster,name="",save=T,plot=T){

	# creates the cluster matrix to plot later on using PlotCluster

  cluster <- as.factor(cluster)
  n <- length(cluster)
  clustermatrix <- matrix(NA,nrow=n,ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      if(cluster[i]==cluster[j]) clustermatrix[i,j] <- cluster[i]
    }
  }
  if(plot) PlotCluster(clustermatrix,name=name,save=save)
  return(clustermatrix)
}

PlotCluster <- function(CorrMatrix,name="",pValueMatrix=NULL,save=T){

	# Plots the cluster matrix from the CreateClusterMatrix function taken as input.

  dt <- as.data.table(CorrMatrix)
  n <- dim(CorrMatrix)[1]
  setnames(dt,as.character(paste(1:n,sep="")))
  dt[,RegionA:= as.character(.I)]
  dt <- melt(dt, id.vars='RegionA',value.name = "Cluster",variable.name = 'RegionB')
  dt[,RegionA := factor(RegionA, levels = paste(c(1:n),sep=""))]
  dt[,RegionB := factor(RegionB, levels = paste(c(n:1),sep=""))]
  dt[,Cluster := as.factor(Cluster)]
  
  plot1 <- ggplot(data=dt) + geom_tile(aes(x=RegionA,y=RegionB,fill=Cluster),color='grey') + 
    scale_fill_brewer(palette = 'Set1') + labs(x='Region A',y='Region B') + 
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

PlotClusteredTimeSeries <- function(){

# Unused plots

  dataset.dt <- as.data.table(dataset)  
  setnames(dataset.dt,paste('Region',1:90,sep=""))
  Dataplot <- as.data.table(cbind(dataset.dt,1:1190))
  setnames(Dataplot, "V2","Time")
  plot1 <- ggplot(data = Dataplot) + geom_line(aes(x=Time, y= Region1, color = cluster[1]))  +
  geom_line(aes(x=Time, y= Region50, color = cluster[50])) +
    geom_line(aes(x=Time, y= Region33, color = cluster[33])) +
    geom_line(aes(x=Time, y= Region72, color = cluster[72])) +
    geom_line(aes(x=Time, y= Region49, color = cluster[49])) +
    geom_line(aes(x=Time, y= Region10, color = cluster[10]))
  ggsave(filename = paste("./Plots/ClusterTimeSeries.png", sep = ""),
         plot = plot1,
         width = 8,
         height = 8,
         units = "in") 
  print(plot1)
}
  
PlotCorrCluster <- function(CorrMatrix,name="",clusterMatrix=NULL,save=T){

	# plots cluster matrix and correlation at the same time on the same plot.
	# unused in the report

    dt <- as.data.table(CorrMatrix)
    setnames(dt,as.character(paste(1:90,sep="")))
    dt[,RegionA:= as.character(.I)]
    dt <- melt(dt, id.vars='RegionA',value.name = "Corr",variable.name = 'RegionB')
    dt[,RegionA := factor(RegionA, levels = paste(c(1:90),sep=""))]
    dt[,RegionB := factor(RegionB, levels = paste(c(90:1),sep=""))]
    
    dt2 <- as.data.table(clusterMatrix)
    setnames(dt2,as.character(paste(1:90,sep="")))
    dt2[,RegionA:= as.character(.I)]
    dt2 <- melt(dt2, id.vars='RegionA',value.name = "Cluster",variable.name = 'RegionB')
    dt2[,RegionA := factor(RegionA, levels = paste(c(1:90),sep=""))]
    dt2[,RegionB := factor(RegionB, levels = paste(c(90:1),sep=""))]
    dt2[Cluster=='NA',Cluster:=""]
    dt2[,Cluster := factor(Cluster)]
    
    dt <- merge(dt,dt2,by=c('RegionA','RegionB'))
    
    plot1 <- ggplot(data=dt) + geom_tile(aes(x=RegionA,y=RegionB,fill=Corr),color='grey') + 
      geom_point(aes(x=RegionA,y=RegionB,color=Cluster),alpha=0.8) +
      scale_fill_gradient2(low="#E41A1C", mid = 'white', high = "#377EB8") + labs(x='Region A',y='Region B') + 
      #scale_color_manual(values = c('#E41A1C','#4DAF4A','#984EA3','#FF7F00')) + #Brewer pal Set1 
      scale_color_brewer(palette = 'Spectral') + 
      theme_classic() + 
      theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())

    print(plot1)
    ggsave(filename = paste("./Plots/ClusterCorrPlot_",name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
           plot = plot1,
           width = 8,
           height = 8,
           units = "in")

    return(plot1)
}

CirclePlot <- function(Part.Corr,kmeans,name){

	# generates the simple correlation circle plot. 

  max <- max(ifelse(is.na(Part.Corr),0,Part.Corr))
  min <- min(ifelse(is.na(Part.Corr),0,Part.Corr))
  
  Part.Corr.tresholded <- ifelse(Part.Corr <= 0.5*max &  Part.Corr >= 0.5*min,NA,Part.Corr)
  
  qgraph(Part.Corr.tresholded, groups = factor(kmeans$cluster,levels=1:8),layout='circle',color = brewer.pal(8,"Set1"),
         filetype='png',filename = paste("./Plots/CirclePlot_",name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),sep = ""),
         width = 12,height = 12,normalize=F)

}

CreateBootsrapConnectivity.Tree <- function(HCP,method=c("pearson","accordance","discordance"),k=7){

	# runs the hierachical dendrogram analysis using k clusters. 
	# choice of dissimilarity of pearson, accordance of discordance
	# returns the bootstrap connectivity as a matrix
	
  dims <- dim(HCP)
  n<- dims[1]
  N <- dims[3]
  ClusterMatrix <- matrix(0,nrow=n,ncol=n)
  
  for(i in 1:N){
    cat(paste("clustering for patient",i,"\n"))
    if(method=="pearson"){
      cor <- cor(t(HCP[,,i]))
      dimnames(cor) <- list(Nodes$name,Nodes$name)
      tree <- hclust(d= as.dist(1-cor),method="complete")
    } else if(method=="accordance"){
      cor <- Concordance$acc[,,i]
      dimnames(cor) <- list(Nodes$name,Nodes$name)
      tree <- hclust(d= as.dist(1-cor),method="complete")
    } else if(method=="discordance"){
      cor <- Concordance$dis[,,i]
      dimnames(cor) <- list(Nodes$name,Nodes$name)
      tree <- hclust(d= as.dist(abs(cor)),method="complete")
    }
    cluster <- cutree(tree,k=k)
    temp <- CreateClusterMatrix(cluster,name=paste("Cluster",k,"tree",method,sep=""),save=F,plot=F)
    temp[temp>0] <- 1
    temp[is.na(temp)] <- 0
    ClusterMatrix <- ClusterMatrix + temp
  }
  
  dimnames(ClusterMatrix) <-list(Nodes$name,Nodes$name)
  PlotCorr(ClusterMatrix/90,name=paste("Cluster",k,"tree",method,sep=""))
  return(ClusterMatrix/90)
  
}

CreateBootsrapConnectivity.Tree.MovieWatching <- function(HCP,method=c("pearson","accordance","discordance"),k=7){

	# runs the hierachical dendrogram analysis using k clusters for the movie watching data !
	# choice of dissimilarity of pearson, accordance of discordance
    # returns the bootstrap connectivity as a matrix
	
  dims <- dim(HCP)
  n<- dims[1]
  N <- dims[3]
  setkey(Nodes,index)
  namer <- Nodes[!(index %in% c(75,76))]$name
  ClusterMatrix <- matrix(0,nrow=n,ncol=n)
  
  for(i in 1:N){
    cat(paste("clustering for patient",i,"\n"))
    if(method=="pearson"){
      cor <- cor(t(HCP[,,i]))
      dimnames(cor) <- list(namer,namer)
      tree <- hclust(d= as.dist(1-cor),method="complete")
    } else if(method=="accordance"){
      cor <- ConcordanceMW$acc[,,i]
      dimnames(cor) <- list(namer,namer)
      tree <- hclust(d= as.dist(1-cor),method="complete")
    } else if(method=="discordance"){
      cor <- ConcordanceMW$dis[,,i]
      dimnames(cor) <- list(Nnamer,namer)
      tree <- hclust(d= as.dist(abs(cor)),method="complete")
    }
    cluster <- cutree(tree,k=k)
    temp <- CreateClusterMatrix(cluster,name=paste("Cluster",k,"tree",method,sep=""),save=F,plot=F)
    temp[temp>0] <- 1
    temp[is.na(temp)] <- 0
    ClusterMatrix <- ClusterMatrix + temp
  }
  
  dimnames(ClusterMatrix) <-list(namer,namer)
  PlotCorr(ClusterMatrix/90,name=paste("Cluster",k,"tree",method,sep=""))
  return(ClusterMatrix/90)
  
}


CreateBootsrapConnectivity.kmeans <- function(HCP){

	# runs the k-means using 7 clusters for each patient individually
  # returns the bootstrap connectivity as a matrix
  
  ClusterMatrix <- matrix(0,nrow=90,ncol=90)
  
  for(i in 1:90){
    cat(paste("clustering for patient",i,"\n"))
    kmean <- kmeans(HCP[,,i],center=7)
    temp <- CreateClusterMatrix(kmean$cluster,name="Cluster7tree",save=F,plot=F)
    temp[temp>0] <- 1
    temp[is.na(temp)] <- 0
    ClusterMatrix <- ClusterMatrix + temp
  }
  PlotCorr(ClusterMatrix/90,name=paste("Cluster7tree",method,sep=""))
  return(ClusterMatrix/90)
  
}

SplittingAnalysis.kmeans<- function(data,k=8,name="",diss=c("euclidian","correlation")){

  # test function to see if it made sense to only cluster on one half of the symmetric pairs.
  # using k means
  #first do the clustering on the complete dataset
  kmeans <- kmeans(x=data,centers=k)
  switch(diss,
         euclidian = {
           sil <- silhouette(kmeans$cluster, dist(data, method = "euclidian"))
         },
         correlation = {
           sil <- silhouette(kmeans$cluster, dmatrix = (1-cor(t(data))))
         })
  png(paste("./Plots/Silhouette/Full_kmean_silhouette_",name,format(Sys.time(),"%d-%m-%Y-%H-%M"),".png",sep=""),width=8,height=8,units="in",res=720)
  plot(sil,col=brewer.pal(8,"Spectral"))
  dev.off()
  
  #second do the clustering on the halved dataset and then retroactively add in the rest.
  index <- unique(1:89%/%2)*2+1
  dataset <-data[index,]
  kmeans2 <- kmeans(x=dataset,centers=8)
  cluster2 <- rep(kmeans2$cluster,each=2)
  switch(diss,
         euclidian = {
           sil2 <- silhouette(cluster2, dist(data, method = "euclidian"))
         },
         correlation = {
           sil2 <- silhouette(cluster2, dmatrix = (1-cor(t(data)))) #choose distance here of course.
         })
  png(paste("./Plots/Silhouette/split_kmeans_silhouette_",name,format(Sys.time(),"%d-%m-%Y-%H-%M"),".png",sep=""),width=8,height=8,units="in",res=720)
  plot(sil2,col=brewer.pal(8,"Spectral"))
  dev.off()
}

SplittingAnalysis.pam<- function(data,k=8,name="",diss=c("euclidian","correlation")){


  # test function to see if it made sense to only cluster on one half of the symmetric pairs.
  # using k medoids
  
  #first do the clustering on the complete dataset
  kmeans <- pam(x=data,k=k)
  switch(diss,
         euclidian = {
           sil <- silhouette(kmeans$cluster, dist(data, method = "euclidian"))
         },
         correlation = {
           sil <- silhouette(kmeans$cluster, dmatrix = (1-cor(t(data))))
         })
  png(paste("./Plots/Silhouette/Full_pam_silhouette_",name,format(Sys.time(),"%d-%m-%Y-%H-%M"),".png",sep=""),width=8,height=8,units="in",res=720)
  plot(sil,col=brewer.pal(8,"Spectral"))
  dev.off()
  
  #second do the clustering on the halved dataset and then retroactively add in the rest.
  index <- unique(1:89%/%2)*2+1
  dataset <-data[index,]
  kmeans2 <- pam(x=dataset,k=8)
  cluster2 <- rep(kmeans2$cluster,each=2)
  switch(diss,
         euclidian = {
           sil2 <- silhouette(cluster2, dist(data, method = "euclidian"))
         },
         correlation = {
           sil2 <- silhouette(cluster2, dmatrix = (1-cor(t(data)))) #choose distance here of course.
         })
  png(paste("./Plots/Silhouette/split_pam_silhouette_",name,format(Sys.time(),"%d-%m-%Y-%H-%M"),".png",sep=""),width=8,height=8,units="in",res=720)
  plot(sil2,col=brewer.pal(8,"Spectral"))
  dev.off()
}

FinalBoostrapClustering <- function(Bootstrap.Conn,method="complete",k=7,plot=T,save=T,name=""){

  # computes the final cluster distribution for the one dimensional case
  # can plot the silhouette but is commented out since the plot is not very useful
  # method used in the report was "complete" but others can be used 

  #the functions plots the final dendrogram and the silhouettes and returns the 
  #clustermatrix and the cluster vector

  tree <- hclust(d= as.dist(1-Bootstrap.Conn),method=method)
  cluster <- cutree(tree,k=k)
  clustermatrix <- CreateClusterMatrix(cluster,name="Cluster7tree",save,plot)
  
  png(paste("./Plots/",name,".png",sep=""))
  sil <- silhouette(cluster, dmatrix = 1-Bootstrap.Conn)
  plot(sil,col=brewer.pal(k,"Set1"),main="7 clusters",do.n.k=F,do.clus.stat=F)
  dev.off()

  PlotTree(tree,cluster,name)
 

  return(list("cluster"=cluster,"ClusterMatrix"=clustermatrix))
  
}

DoubleCorrelationCirclePlot <- function(Nodes,pcormatrix){
  
  #creates and plots the double correlation circle in the appendix B
  #inputs are straighforward. Note that the Nodes dataframe must contain a column
  #with the cluster attribute !

  
  #create the coordinates for the vertices
  k <- length(unique(Nodes$cluster))
  xc <- sin(2*pi/k*Nodes$cluster)#+runif(length(cluster))*0.7
  yc <- cos(2*pi/k*Nodes$cluster)#+runif(length(cluster))*0.7
  
  
  #create vertex data frame
  vertex.data.frame <- as.data.table(cbind(xc,yc,Nodes))
  
  vertex.data.frame[,Tot := .N,.(cluster)]
  vertex.data.frame[,Nb := rowid(cluster)]
  
  vertex.data.frame[,xv := xc+sin(2*pi/Tot*Nb)*0.1*log(Tot)]
  vertex.data.frame[,yv := yc+cos(2*pi/Tot*Nb)*0.1*log(Tot)]
  vertex.data.frame[,xtext := xc+sin(2*pi/Tot*Nb)*(0.1*log(Tot)+0.05)]
  vertex.data.frame[,ytext := yc+cos(2*pi/Tot*Nb)*(0.1*log(Tot)+0.05)]
  
  #create the edges data frame
  formvector<-function(){
    #this function creates the vector to create the dataframe of "edges" in the pcor network plot
    n <- length(Nodes$cluster)
    temp <- cbind(rep(1,n),1:n)
    for(i in 2:90){
      temp <- rbind(temp,cbind(rep(i,n),1:n))
    }
    return(temp)
  }
  
  NewVector <- formvector()
  #get pcor corresponding to the relation
  edge.data.frame <- as.data.table(cbind(NewVector,c(pcormatrix),1:nrow(NewVector)))
  setnames(edge.data.frame,c("nodes1","nodes2","pcor","SortingIndex"))
  #get begin coor
  edge.data.frame <- merge(edge.data.frame,vertex.data.frame,by.x="nodes1",by.y="index")
  setnames(edge.data.frame,c("xv","yv","xc","yc"),c("xbegin","ybegin","xcbegin","ycbegin"))
  #get end coor
  edge.data.frame <- merge(edge.data.frame,vertex.data.frame,by.x="nodes2",by.y="index")
  setnames(edge.data.frame,c("xv","yv","xc","yc","cluster.x","name.x","cluster.y","name.y"),
           c("xend","yend","xcend","ycend","clusterbegin","namebegin","clusterend","name.end"))
  setkey(edge.data.frame,"SortingIndex")
  
  #set inclustser and oucluster edge data frame
  incluster.edge.data.frame <- edge.data.frame[clusterbegin==clusterend]
  outcluster.edge.data.frame <- edge.data.frame[clusterbegin!=clusterend]
  
  #keep only high pcor for incluster
  incluster.edge.data.frame[is.na(pcor),pcor:=0]
  max <- max(incluster.edge.data.frame$pcor)
  min <- min(incluster.edge.data.frame$pcor)
  maxout <- max(outcluster.edge.data.frame$pcor)
  minout <- min(outcluster.edge.data.frame$pcor)
  
  incluster.edge.data.frame.positive <- incluster.edge.data.frame[pcor >= 0.3*max]
  incluster.edge.data.frame.negative <- incluster.edge.data.frame[pcor <= 0.3 * min]
  
  #keep only max pcor for outcluster
  outcluster.edge.data.frame <- outcluster.edge.data.frame[,c("pcormax","pcormin") := .(max(pcor),min(pcor)),.(xcend,ycend,xcbegin,ycbegin,clusterbegin,clusterend)]
  
  outcluster.edge.data.frame <- outcluster.edge.data.frame[pcormax >= 0.75*maxout | pcormin <= 0.75*minout]
  
  #plot
  plot1 <- ggplot() + 
    geom_segment(data=incluster.edge.data.frame.positive,aes(x=xbegin,y=ybegin,xend=xend,yend=yend,size=abs(pcor)),color="darkgreen")+
    geom_segment(data=incluster.edge.data.frame.negative,aes(x=xbegin,y=ybegin,xend=xend,yend=yend,size=abs(pcor)),color="darkred")+
    geom_point(data=vertex.data.frame,aes(xv,yv,color=as.factor(cluster)),size=5) + 
    geom_point(data=vertex.data.frame,aes(xc/2,yc/2,color=as.factor(cluster)),size=3) +
    geom_segment(data=outcluster.edge.data.frame,aes(x=xcbegin/2,y=ycbegin/2,xend=xcend/2,yend=ycend/2,size=abs(pcormax)),color="grey")+
    geom_text(data=vertex.data.frame,aes(xtext,ytext,label=name))  +
    scale_size(range = c(0, 2)) + scale_color_brewer(palette="Set1") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  ggsave(filename = paste("./Plots/CirclePlot_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = plot1,
         width = 11,
         height = 8,
         units = "in")
  
  print(plot1)
}

PlotTree <- function(tree,cluster,name=""){

	# plots the dendrogram using ggdendro package as a base

  dend_data <- dendro_data(tree, type = "rectangle")
  dend_data$labels <- merge(dend_data$labels, data.table(cluster,row.names(as.data.frame(cluster))),by.x="label",by.y="V2")
  if(row.names(as.data.frame(cluster))[1]=="1"){
    dend_data$labels <- as.data.table(merge(dend_data$labels, Nodes,by.x="label",by.y="index"))
    dend_data$labels[,label := name]
    dend_data$labels[,x := x.x]
    dend_data$labels[,y := y.x]
  }
  p <- ggplot(dend_data$segments) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    geom_text(data = dend_data$labels, aes(x, y, label = label),
              nudge_y = -0.1, angle = 90, size = 3)+
    geom_point(data = dend_data$labels, aes(x, y, color = as.factor(cluster)),size=2) +
    ylim(-0.25, 1) + guides(color=guide_legend(title="Cluster")) + scale_color_brewer(palette =  "Set1") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +
    ylab("Height (Pearson Correlation)")
  ggsave(filename = paste("./Plots/Dendrogram_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png", sep = ""),
         plot = p,
         width = 11,
         height = 8,
         units = "in")
  print(p)
}

silhouette_MovieWatching <- function(){

	# plots the silhouette for the movie watching data
	
  setkey(Nodes,index)
  cluster <- Nodes[!(index %in% c(75,76))]$cluster
  distance <- apply(ConcordanceMW$acc,c(1,2),mean)
  sil <- silhouette(cluster, dmatrix = as.matrix(1-distance))
  png(paste("./Plots/","Sil_MovieWatchingAccordance",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png",sep=""),width=8,height=8,units="in",res=720)
  plot(sil,col=brewer.pal(7,"Set1"),main="7 clusters",do.n.k=F,do.clus.stat=F)
  dev.off()
}

silhouette_MCI <- function(ConcordanceMW,name="",patient=1,mean=T){

	# plots the silhouette for the MCI data
	
  setkey(Nodes,index)
  cluster <- Nodes[!(index %in% c(75,76))]$cluster
  if(mean){
    distance <- apply(ConcordanceMW$acc,c(1,2),mean)
    distance2 <- abs(apply(ConcordanceMW$dis,c(1,2),mean))
    name<-paste(name,"mean",sep="")
  }else{
    distance <- ConcordanceMW$acc[,,patient]
    distance2 <- abs(ConcordanceMW$dis[,,patient])
    name<-paste(name,patient,sep="")
  }
  sil <- silhouette(cluster, dmatrix = as.matrix(1-distance))
  png(paste("./Plots/","Sil_MCIAccordance_",name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png",sep=""),width=8,height=8,units="in",res=720)
  plot(sil,col=brewer.pal(7,"Set1"),main="7 clusters",do.n.k=F,do.clus.stat=F)
  dev.off()
  sil <- silhouette(cluster, dmatrix = as.matrix(distance2))
  png(paste("./Plots/","Sil_MCIdiscordance_",name,"_",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png",sep=""),width=8,height=8,units="in",res=720)
  plot(sil,col=brewer.pal(7,"Set1"),main="7 clusters",do.n.k=F,do.clus.stat=F)
  dev.off()
}

silhouette_MCI_test <- function(ConcordanceMCI,patient){

	# This function returns some statistics to play around the silhouette width with the MCI data
	
  setkey(Nodes,index)
  cluster <- Nodes[!(index %in% c(75,76))]$cluster
  distance <- ConcordanceMCI$acc[,,patient]
  distance2 <- abs(ConcordanceMCI$dis[,,patient])
  sil <- as.data.table(silhouette(cluster, dmatrix = as.matrix(1-distance))[1:88,])
  sil2 <- as.data.table(silhouette(cluster, dmatrix = as.matrix(distance2))[1:88,])
  result1 <- sum(sil$sil_width>0)/88
  result2 <- sum(sil2$sil_width>0)/88
  result <- c(result1,mean(sil$sil_width),result2,mean(sil2$sil_width))
  return(result)
}

silhouette_All <- function(){

	# plots the silhouette for the main dataset
		
  setkey(Nodes,index)
  cluster <- Nodes$cluster
  distance <- apply(Concordance$acc,c(1,2),mean)
  sil <- silhouette(cluster, dmatrix = as.matrix(1-distance))
  png(paste("./Plots/","Sil_MainDataAccordance",format(Sys.time(),"%d-%m-%Y-%H-%M"),".png",sep=""),width=8,height=8,units="in",res=720)
  plot(sil,col=brewer.pal(7,"Set1"),main="7 clusters",do.n.k=F,do.clus.stat=F)
  dev.off()
}

#### EXECUTION EXAMPLES ####


## Script to get the final clusters and save them as a file
## 7 clusters, binarized dataset and accordance based dendrogram are used as input parameters
#BoostrapConnectivity <- CreateBootsrapConnectivity.Tree(HCP.bin,method="accordance",k=7)
#FinalBoostrapCluster <- FinalBoostrapClustering(BoostrapConnectivity,method="complete",name="accordanceboostrap")
#FinalBoostrapCluster$connectivity <- BoostrapConnectivity
#save(FinalBoostrapCluster, file="FinalBoostrapClustering.rda")

## Plots the double correlation circle, with only one node per symmetrical pair
#Nodes2 <- Nodes[index,]
#Nodes2$name <- gsub(".L","",Nodes2$name)
#DoubleCorrelationCirclePlot(Nodes2,apply(ListofAllPart.Corrbin,c(1,2),mean)[index,index])
  

