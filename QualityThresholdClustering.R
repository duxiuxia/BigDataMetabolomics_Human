# Author: Matthew Deitz

# Editor: Xiuxia Du

# Started in June 2015






PPMTolerance=10



if ( !is.element("xlsx", installed.packages()[,1]) ) {
    install.packages("xlsx")
}

library(xlsx)




# ==============================================================
# set working directory here
# ==============================================================






data1=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV1")
data2=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV2")
data3=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV3")
data4=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV4")
data5=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV5")
data6=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV6")
data7=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV7")
#Reads in the data from the worksheet





mzvalues <- rbind(data1,data2,data3,data4,data5,data6,data7)
#combines all the rows into one dataset



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




findSmallestDist<-function(mzvalues){
  i=2
  smallestDist=mzvalues[2]-mzvalues[1]
  smallestIndex=1
  while (i<(length(mzvalues)-1)){
    if (smallestDist>(mzvalues[i+1]-mzvalues[i])){
      smallestDist = mzvalues[i+1]-mzvalues[i]
      smallestIndex=i
    }
    i=i+1
  }
  return(smallestIndex)
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
