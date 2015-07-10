
setwd("C:/Users/matt/Desktop/NISTMetabolites")

results<-data.frame(Name=character(),InChIKey=character(),Formula=character(),MW=character(),ExactMass=character(),CASNO=character())

fileLines<-readLines(file("mainlib.MSP",open="r"))
i=1
while (i < length(fileLines)){
  if(length(grep("Name:",fileLines[i]))>0){
    cpdName=substr(gsub(" ","",fileLines[i]),6,nchar(gsub(" ","",fileLines[i])))
    if(length(grep("Related_CAS#:",fileLines[i+1]))>0){
      i=i+1
    }
    cpdInChIKey=substr(gsub(" ","",fileLines[i+1]),10,nchar(gsub(" ","",fileLines[i+1])))
    cpdFormula=substr(gsub(" ","",fileLines[i+2]),9,nchar(gsub(" ","",fileLines[i+2])))
    cpdMW=substr(gsub(" ","",fileLines[i+3]),4,nchar(gsub(" ","",fileLines[i+3])))
    cpdExactMass=substr(gsub(" ","",fileLines[i+4]),11,nchar(gsub(" ","",fileLines[i+4])))
    cpdCASNO=substr(gsub(" ","",fileLines[i+5]),7,nchar(gsub(" ","",fileLines[i+5])))
    
    tmp<-data.frame(Name=cpdName,InChIKey=cpdInChIKey,Formula=cpdFormula,MW=cpdMW,ExactMass=cpdExactMass,CASNO=cpdCASNO)
    results<-rbind(results,tmp)
    
    numPeaks = substr(gsub(" ","",fileLines[i+8]),11,nchar(gsub(" ","",fileLines[i+8])))
    
    i=i+8+strtoi(numPeaks)
  }
  i=i+1
}

write.csv(results,file="NISTDatabaseResults.csv")
