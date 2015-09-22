# Author: Matthew Deitz

# Editor: Xiuxia Du

# Started in June 2015





rm(list=ls())



PPMTolerance=3

if (!is.element("rJava", installed.packages()[,1])) {
    install.packages("rJava", lib="/projects/dulab_research/R_Libs", repos="http://cran.r-project.org")
}

library(rJava, lib.loc="/projects/dulab_research/R_Libs")

if (!is.element("xlsxjars", installed.packages()[,1])) {
    install.packages("xlsxjars", lib="/projects/dulab_research/R_Libs", repos="http://cran.r-project.org")
}

library(xlsxjars, lib.loc="/projects/dulab_research/R_Libs")

if (!is.element("xlsx", installed.packages()[,1])) {
    install.packages("xlsx", lib="/projects/dulab_research/R_Libs", repos="http://cran.r-project.org")
}

library(xlsx, lib.loc="/projects/dulab_research/R_Libs")

if (!is.element("kulife", installed.packages()[,1])) {
    install.packages("kulife", lib="/projects/dulab_research/R_Libs", repos="http://cran.r-project.org")
}

library(kulife, lib.loc="/projects/dulab_research/R_Libs")


# ==============================================================
# set working directory here
# ==============================================================


#setwd("/projects/dulab_research/Datasets/PubChem/ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/CompletedFiles/")
setwd("/Users/matthewdeitz/Desktop")



# ==============================================
# !!! set working directory using setwd() here
# ==============================================
data<-data.frame(CenterMass = character(), IUPAC_Name=character(),PubChem_ID=character(),Formula=character(),
                    Exact_Mass=character(), InChiKey = character())
temp = list.files(pattern="*.csv")
for (filename in temp)
{
    tempdata=read.csv(filename)
    tempdata$X=NULL
    tempdata<-tempdata[-which(is.na(tempdata$Exact_Mass)),]
    data<-rbind(data,tempdata)
}


#Reads in the data from the worksheet
mzvalues <- sort(data$Exact_Mass)
data <- data[order(data$Exact_Mass),]
#Orders the mzvalues by the Query value


results<-data.frame(CenterMass = character(), IUPAC_Name=character(),PubChem_ID=character(),Formula=character(),
                    Exact_Mass=character(), InChIKey = character(), Molecular_Weight = character())
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
        cIUPAC_Name=""
        cPubChem_ID=""
        cInChIKey=""
        cFormula=""
        cExact_Mass=""
        cMolecular_Weight=""
        cCenterMass = (mzvalues[smallIndex]+mzvalues[largeIndex])/2
        for(cValue in removeValues){
            cIUPAC_Name<-paste(cIUPAC_Name,data[cValue,]$IUPAC_Name,sep=", ")
            cPubChem_ID<-paste(cPubChem_ID,data[cValue,]$PubChem_ID,sep=", ")
            cInChIKey<-paste(cInChIKey,data[cValue,]$InChIKey,sep=", ")
            cFormula<-paste(cFormula,data[cValue,]$Formula,sep=", ")
            cExact_Mass<-paste(cExact_Mass,data[cValue,]$Exact_Mass,sep=", ")
            cMolecular_Weight<-paste(cMolecular_Weight,data[cValue,]$Molecular_Weight,sep=", ")
            #cMonisotopic_weight<-paste(cMonisotopic_weight,data[cValue,]$Monisotopic_weight,sep=", ")
        }
        data <- data[-removeValues,]
        mzvalues<-mzvalues[-removeValues]
        
        results<-rbind(results,data.frame(CenterMass = cCenterMass, IUPAC_Name=cIUPAC_Name,PubChem_ID=cPubChem_ID,Formula=cFormula,Exact_Mass=cExact_Mass, InChIKey=cInChIKey, Molecular_Weight = cMolecular_Weight))
    }
    else{
        endFlag=F
        tempdf <- data.frame(lapply(data, as.character), stringsAsFactors = FALSE)
        templist <- as.character(data$Exact_Mass)
        tempdf$CenterMass <- templist
        results<-rbind(results,tempdf)
    }
}
write.xml(results, file="PubChemClusterResults.xml")

