


######################################
#Set User Directory Here

setwd("C:/Users/matt/Desktop/NISTMetabolites")
######################################


results<-data.frame(Name=character(),InChIKey=character(),Formula=character(),MW=character(),ExactMass=character(),CASNO=character())
#create the results data frame and open the database file .MSP
fileLines<-readLines(file("mainlib.MSP",open="r"))
i=1
while(i<length(fileLines)){
    #Follows the following format to get the name and then check to see if it has the 
    #following exceptions Related CAS#, Salt, or Known Impurity
    
  if(length(grep("Name:",fileLines[i]))>0){
    cpdName=substr(gsub(" ","",fileLines[i]),6,nchar(gsub(" ","",fileLines[i])))
    if(length(grep("Related_CAS#",fileLines[i+1]))>0){
      i=i+1
    }
    if(length(grep("Salt",fileLines[i+1]))>0){
      i=i+1
    }
    if(length(grep("Known_impurity",fileLines[i+1]))>0){
      i=i+1
    }
    
    #gets the rest of the information which is the same format for every compound
    cpdInChIKey=substr(gsub(" ","",fileLines[i+1]),10,nchar(gsub(" ","",fileLines[i+1])))
    cpdFormula=substr(gsub(" ","",fileLines[i+2]),9,nchar(gsub(" ","",fileLines[i+2])))
    cpdMW=substr(gsub(" ","",fileLines[i+3]),4,nchar(gsub(" ","",fileLines[i+3])))
    
    
    if(length(grep("ExactMass:",fileLines[i+4]))>0){
      cpdExactMass=substr(gsub(" ","",fileLines[i+4]),11,nchar(gsub(" ","",fileLines[i+4])))
    }
    else{
      cpdExactMass="NA"
      i=i-1
    }
    
    
    
    cpdCASNO=substr(gsub(" ","",fileLines[i+5]),7,nchar(gsub(" ","",fileLines[i+5])))
    #add the results to the dataframe results
    tmp<-data.frame(Name=cpdName,InChIKey=cpdInChIKey,Formula=cpdFormula,MW=cpdMW,ExactMass=cpdExactMass,CASNO=cpdCASNO)
    results<-rbind(results,tmp)
    
    #The following code is to speed up the code through the while loop
    if(length(grep("Comment:",fileLines[i+7]))>0){
      numPeaks=substr(gsub(" ","",fileLines[i+8]),10,nchar(gsub(" ","",fileLines[i+8])))
      i=i+8+as.integer(numPeaks)
    }
    else{
      numPeaks=substr(gsub(" ","",fileLines[i+7]),10,nchar(gsub(" ","",fileLines[i+7])))
      i=i+7+as.integer(numPeaks)
    }
  }
  i=i+1
}

write.csv(results,file="NISTDatabaseResults.csv")
