# Author: Matthew Deitz

# Editor: Xiuxia Du

# Started in June 2015





rm(list=ls())



PPMTolerance=10





if (!is.element("xlsx", installed.packages()[,1])) {
    install.packages("xlsx")
}

library(xlsx)




# ==============================================================
# set working directory here
# ==============================================================


setwd("C:/Users/matt/Desktop/MetaboliteDB/")




# ==============================================
# !!! set working directory using setwd() here
# ==============================================

data=read.csv("HMDBMetaboliteList.csv")

#Reads in the data from the worksheet
mzvalues <- sort(data$Monisotopic_weight)

mzvalues <- mzvalues[ with(mzvalues,order(Query_mz)),]
#Orders the mzvalues by the Query value


results<-list()
endFlag=T

while (endFlag) {
  index <- findSmallestDist_2(mzvalues)
  
  if (((mzvalues[index+1]-mzvalues[index])/mzvalues[index]*1000000)<PPMTolerance){
    
    indexValues<-extendCluster(mzvalues,index,index+1,PPMTolerance)
    
    smallIndex=indexValues[1]
    largeIndex=indexValues[2]
    cluster = mzvalues[smallIndex:largeIndex]
    
    results[[length(results)+1]]<-cluster
    
    removeValues=c(smallIndex:largeIndex)
    mzvalues<-mzvalues[-removeValues]
    
    
  }
  else{
    endFlag=F
    results<-c(results,mzvalues)
  }
}





findSmallestDist_2 <- function(x) {
    x1 <- x[-length(x)]
    x2 <- x[-1]
    all_distance <- x2 - x1
    II <- which(all_distance == min(all_distance))
    return(II)
}


extendCluster <- function(mzvalues,smallIndex,largeIndex,PPMTolerance){
  returnSmallIndex <- smallIndex
  returnLargeIndex <- largeIndex

  if (returnSmallIndex-1==0){
    if (((mzvalues[returnSmallIndex]+mzvalues[returnLargeIndex+1])/2
         -mzvalues[returnSmallIndex])/mzvalues[returnSmallIndex]
        *1000000<PPMTolerance){
      return(extendCluster(mzvalues,returnSmallIndex,returnLargeIndex+1,PPMTolerance))
    }
    else{
      return(c(returnSmallIndex,returnLargeIndex))
    }
  }
  
  if (returnLargeIndex+1==length(mzvalues)){
    if (((mzvalues[returnSmallIndex-1]+mzvalues[returnLargeIndex])/2
         -mzvalues[returnSmallIndex-1])/mzvalues[returnSmallIndex-1]
        *1000000<PPMTolerance){
      return(extendCluster(mzvalues,returnSmallIndex-1,returnLargeIndex,PPMTolerance))
    }
    else{
      return(c(returnSmallIndex,returnLargeIndex))
    }
  }
  
  lowerDist=mzvalues[returnSmallIndex]-mzvalues[returnSmallIndex-1]
  upperDist=mzvalues[returnLargeIndex+1]-mzvalues[returnLargeIndex]
  
  if (lowerDist<upperDist){
    if (((mzvalues[returnSmallIndex-1]+mzvalues[returnLargeIndex])/2
         -mzvalues[returnSmallIndex-1])/mzvalues[returnSmallIndex-1]
        *1000000<PPMTolerance){
      return(extendCluster(mzvalues,returnSmallIndex-1,returnLargeIndex,PPMTolerance))
    }
    else{
      return(c(returnSmallIndex,returnLargeIndex))
    }
  }
  else{
    if (((mzvalues[returnSmallIndex]+mzvalues[returnLargeIndex+1])/2
         -mzvalues[returnSmallIndex])/mzvalues[returnSmallIndex]
        *1000000<PPMTolerance){
      return(extendCluster(mzvalues,returnSmallIndex,returnLargeIndex+1,PPMTolerance))
    }
    else{
      return(c(returnSmallIndex,returnLargeIndex))
    }
  }
}

xx<-lapply(results,unlist)
max<-max(sapply(xx,length))


write.csv(do.call(rbind,lapply(xx,function(z)c(z,rep(NA,max-length(z))))),"QTClusterResults.csv")
