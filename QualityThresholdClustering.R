# Author: Matthew Deitz

# Editor: Xiuxia Du

# Started in June 2015





rm(list=ls())



PPMTolerance=3

if (!is.element("rJava", installed.packages()[,1])) {
    install.packages("rJava")
}

library(rJava)

if (!is.element("xlsxjars", installed.packages()[,1])) {
    install.packages("xlsxjars")
}

library(xlsxjars)

if (!is.element("xlsx", installed.packages()[,1])) {
    install.packages("xlsx")
}

library(xlsx)

if (!is.element("kulife", installed.packages()[,1])) {
    install.packages("kulife")
}

library(kulife)


# ==============================================================
# set working directory here
# ==============================================================


setwd("/users/matthewdeitz/Desktop/MetaboliteDB/")




# ==============================================
# !!! set working directory using setwd() here
# ==============================================



data=read.csv("NISTDatabaseResults.csv")
data$X=NULL
#Reads in the data from the worksheet
mzvalues <- sort(data$ExactMass)
data <- data[order(data$ExactMass),]
#Orders the mzvalues by the Query value


results<-data.frame(CenterMass = character(), Name=character(),InChIKey=character(),Formula=character(),
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
        cCenterMass = (mzvalues[smallIndex]+mzvalues[largeIndex])/2
        for(cValue in removeValues){
            cName<-paste(cName,data[cValue,]$Name,sep=", ")
            cInChIKey<-paste(cInChIKey,data[cValue,]$InChIKey,sep=", ")
            cFormula<-paste(cFormula,data[cValue,]$Formula,sep=", ")
            cMW<-paste(cMW,data[cValue,]$MW,sep=", ")
            cExactMass<-paste(cExactMass,data[cValue,]$ExactMass,sep=", ")
            cCASNO<-paste(cCASNO,data[cValue,]$CASNO,sep=", ")
        }
        data <- data[-removeValues,]
        mzvalues<-mzvalues[-removeValues]
        
        results<-rbind(results,data.frame(CenterMass = cCenterMass, Name=cName,InChIKey=cInChIKey,Formula=cFormula,MW=cMW,ExactMass=cExactMass,CASNO=cCASNO))
    }
    else{
        endFlag=F
        for(cValue in mzvalues){
            temdf <- data
            templist <- data$ExactMass
            tempdf$CenterMass <- termplist
            results<-rbind(results,tempdf)
        }
    }
}
write.xml(results, "NistClusters.xml")
