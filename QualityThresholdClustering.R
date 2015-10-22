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


#Define which Database you want to read
data=read.csv("KEGGDatabaseResults.csv")
data$X=NULL
#Reads in the data from the worksheet
mzvalues <- sort(data$Exact_Mass)
data <- data[order(data$Exact_Mass),]
#Orders the mzvalues by the Query value


results<-data.frame(CenterMass = character(), Name=character(),Entry=character(),Formula=character(),
                    Exact_Mass=character())
endFlag=T


#A function that finds the compounds with the smallest mass difference
findSmallestDist_2 <- function(x) {
    x1 <- x[-length(x)]
    x2 <- x[-1]
    all_distance <- x2 - x1
    II <- which(all_distance == min(all_distance))
    return(II)
}

# a recursive function that will take the cluster indexes and try to extend them as 
#much as possible without going over the user defined PPMTolerance
extendCluster <- function(mzvalues,smallIndex,largeIndex,PPMTolerance){
  returnSmallIndex <- smallIndex
  returnLargeIndex <- largeIndex

  #if the cluster is the first element in the array
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
  
  #if it is the last element in the array
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
  
  #if not first or last then get the distance in both directions
  lowerDist=mzvalues[returnSmallIndex]-mzvalues[returnSmallIndex-1]
  upperDist=mzvalues[returnLargeIndex+1]-mzvalues[returnLargeIndex]
  
  #whichever distance is smaller try adding it to the cluster
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
    #get the index for the smallest distance compounds
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
        cEntry=""
        cFormula=""
        cExact_Mass=""
        cCenterMass = (mzvalues[smallIndex]+mzvalues[largeIndex])/2
        # for all the vales in the cluster
        for(cValue in removeValues){
            cName<-paste(cName,data[cValue,]$Name,sep=", ")
            cEntry<-paste(cEntry,data[cValue,]$Entry,sep=", ")
            cFormula<-paste(cFormula,data[cValue,]$Formula,sep=", ")
            cExact_Mass<-paste(cExact_Mass,data[cValue,]$Exact_Mass,sep=", ")
            #cMonisotopic_weight<-paste(cMonisotopic_weight,data[cValue,]$Monisotopic_weight,sep=", ")
        }
        #remove the values from the list so they can't be added to any other clusters
        data <- data[-removeValues,]
        mzvalues<-mzvalues[-removeValues]
        
        #add the results
        results<-rbind(results,data.frame(CenterMass = cCenterMass, Name=cName,Entry=cEntry,Formula=cFormula,Exact_Mass=cExact_Mass))
    }
    else{
        endFlag=F
        #add all the results that aren't clustered as single clusters
        tempdf <- data.frame(lapply(data, as.character), stringsAsFactors = FALSE)
        templist <- as.character(data$Exact_Mass)
        tempdf$CenterMass <- templist
        results<-rbind(results,tempdf)
    }
}
write.xml(results, file="KEGGClusterResults.xml")

