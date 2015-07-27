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



data=read.csv("NISTDatabaseResults.csv")
data$X=NULL
#Reads in the data from the worksheet
mzvalues <- sort(data$ExactMass)
#Orders the mzvalues by the Query value


results<-data.frame(Name=character(),InChIKey=character(),Formula=character(),
                    MW=character(),ExactMass=character(),CASNO=character())
endFlag=T

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

while (endFlag) {
    index <- findSmallestDist_2(mzvalues)
    if(length(index)>1){
        index <- index[1]
    }
    
    if (((mzvalues[index+1]-mzvalues[index])/mzvalues[index]*1000000)<PPMTolerance){
        
        indexValues<-extendCluster(mzvalues,index,index+1,PPMTolerance)
        
        smallIndex=indexValues[1]
        largeIndex=indexValues[2]
        
        
        
        removeValues=c(smallIndex:largeIndex)
        cName=""
        cInChIKey=""
        cFormula=""
        cMW=""
        cExactMass=""
        cCASNO=""
        
        for(cValue in removeValues){
            cName<-paste(cName,data[which(data$ExactMass==mzvalues[cValue]),]$Name,sep=" ")
            cInChIKey<-paste(cInChIKey,data[which(data$ExactMass==mzvalues[cValue]),]$InChIKey,sep=" ")
            cFormula<-paste(cFormula,data[which(data$ExactMass==mzvalues[cValue]),]$Formula,sep=" ")
            cMW<-paste(cMW,data[which(data$ExactMass==mzvalues[cValue]),]$MW,sep=" ")
            cExactMass<-(mzvalues[smallIndex]+mzvalues[largeIndex])/2
            cCASNO<-paste(cCASNO,data[which(data$ExactMass==mzvalues[cValue]),]$CASNO,sep=" ")
        }
        mzvalues<-mzvalues[-removeValues]
        
        results<-rbind(results,data.frame(Name=cName,InChIKey=cInChIKey,Formula=cFormula,MW=cMW,ExactMass=cExactMass,CASNO=cCASNO))
    }
    else{
        endFlag=F

        results<-rbind(results,data[which(data$ExactMass==as.character(mzvalues)),])
    }
}

write.csv(results,"NISTQTClusterResults.csv")
