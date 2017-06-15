######################################################
#
# File containing the function definition pertaining to
# the MCI application section. Contains data loading,
# Manipulation, and random forest predictive model.
#
# Note : need to load both LoadAdditionalData_MCI() and "./Data/MarkovCoeffMCI.rda"
# to the global environment
#
#######################################################

CreateDataFrame <- function(MarkovCoeffMCI){

  #function to create a workable, comprehensive data table of the Markov Coeff MCI dataset to plot it or build predictive model
  # contains the lambda coefficients of the CTMC, the accordance and concordance coefficients, and the silhouette witdh.

  final.df <- data.table()
  for(i in 1:dim(MarkovCoeffMCI$sCTL)[4]){ #patient loop
    patient.df <- c()
    for(j in 1:dim(MarkovCoeffMCI$sCTL)[3]){ #region loop
      patient.df <- c(patient.df,diag(MarkovCoeffMCI$sCTL[,,j,i]))
    }
    patient.df <- c(patient.df,mean(ConcordanceMCI$sCTL$acc[,,i]),mean(ConcordanceMCI$sCTL$dis[,,i]),
                    silhouette_MCI_test(ConcordanceMCI$sCTL,i)) 
    final.df<- rbind(final.df,t(patient.df))
  }
  final.df[,group:="sCTL"]
  
  for(i in 1:dim(MarkovCoeffMCI$dCTL)[4]){ #patient loop
    patient.df <- c()
    for(j in 1:dim(MarkovCoeffMCI$dCTL)[3]){ #region loop
      patient.df <- c(patient.df,diag(MarkovCoeffMCI$dCTL[,,j,i]))
    }
    patient.df <- c(patient.df,mean(ConcordanceMCI$sCTL$acc[,,i]),mean(ConcordanceMCI$sCTL$dis[,,i]),
                    silhouette_MCI_test(ConcordanceMCI$sCTL,i)) 
    final.df<- rbind(final.df,t(patient.df),fill=T)
  }
  final.df[is.na(group),group:="dCTL"]
  
  for(i in 1:dim(MarkovCoeffMCI$MCI)[4]){ #patient loop
    patient.df <- c()
    for(j in 1:dim(MarkovCoeffMCI$MCI)[3]){ #region loop
      patient.df <- c(patient.df,diag(MarkovCoeffMCI$MCI[,,j,i]))
    }
    patient.df <- c(patient.df,mean(ConcordanceMCI$sCTL$acc[,,i]),mean(ConcordanceMCI$sCTL$dis[,,i]),
                    silhouette_MCI_test(ConcordanceMCI$sCTL,i)) 
    final.df<- rbind(final.df,t(patient.df),fill=T)
  }
  final.df[is.na(group),group:="MCI"]

  Nodes2 <- Nodes[!(index %in% c(75,76))]
  setkey(Nodes2,index)
  setnames(final.df,c(paste(rep(1:88,each=3),rep(c("L_-1","L_0","L_1"),88),sep=""),"meanAcc","meanDis",
                      "silAcc","meanSilAcc","silDis","meanSilDis","group"))
  final.df <- final.df[group!="dCTL"]
  final.df[,group:=as.factor(group)]
}

ComputeNumberOfSignificantlyDifferentSisters <- function(data){

  #function to compute the number of symetric pairs with statically different coefficient for each of the patient pools

  setkey(Nodes,index)
  p.vector.f <- data.table()
  index2 <- unique(1:87%/%2)*2+1
  for(i in index2){
    ii <- ifelse(i<=74,i,i+2)
    p.vector <- c(ii,Nodes[index==ii]$name)
    p.vector <- c(p.vector,
               mean(BootstrapTestTwoSided(data[1,1,i,],data[1,1,i+1,])),
               mean(BootstrapTestTwoSided(data[2,2,i,],data[2,2,i+1,])),
               mean(BootstrapTestTwoSided(data[3,3,i,],data[3,3,i+1,])))
    p.vector.f <- rbind(p.vector.f,t(p.vector))
  }
  setnames(p.vector.f,c("index","name","-1","0","1"))
  return(p.vector.f)
}


#### EXECUTION EXAMPLES ####

##Data loading if not done before
# MCIData <- LoadAdditionalData_MCI()
# load("./Data/MarkovCoeffMCI.rda")

##run a random forest on the most significant regions
# Complete.df <- CreateDataFrame(MarkovCoeffMCI)
# Complete.df <- Complete.df[,.(`34L_-1`,`42L_-1`,`33L_-1`,`70L_-1`,`17L_-1`,`78L_-1`,`61L_-1`,
#                               `56L_1`,`42L_1`,`66L_1`,`81L_1`,`20L_1`,
#                               `17L_0`,`66L_0`,`18L_0`,`19L_0`,`70L_0`,`34L_0`,group)]
# rf <- randomForest(x= Complete.df[,!("group"),with=F],y= Complete.df$group)

##run a PCA on the dataset
# pca <- prcomp(Complete.df[,!"group",with=F])
# plot1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
#               groups = Complete.df$group, ellipse = TRUE, 
#               circle = TRUE)  + scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', 
#               legend.position = 'top')
# print(plot1)


## compute the different sisters for all 3 patient pools
# sCTLdiff<-ComputeNumberOfSignificantlyDifferentSisters(MarkovCoeffMCI$sCTL)
# dCTLdiff <- ComputeNumberOfSignificantlyDifferentSisters(MarkovCoeffMCI$dCTL)
# MCIdiff <- ComputeNumberOfSignificantlyDifferentSisters(MarkovCoeffMCI$MCI)
