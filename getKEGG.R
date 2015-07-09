# getKEGG.R

# Written in June 2015.
setwd("C:/Users/matt/Desktop/KEGGMetabolites")

results<-data.frame(Entry=character(),Name=character(),Formula=character(),Exact_Mass=character())

masses <- seq(from=50, to=1500, by=0.01)

for (i in 1:length(masses)) {
    current_mass <- as.character(masses[i])
    
    query_1 <- paste("http://rest.kegg.jp", "find", "compound", 
                   current_mass, "exact_mass", sep=.Platform$file.sep)
    
    tmp<-try(read.table(query_1),silent=T)
   
    if (inherits(tmp,'try-error')){
      #If an empty file do nothing
      
      
      }
    else{
      dataIn <- read.table(query_1)
    
      compounds <- as.character(dataIn$V1)
    
      for (c in 1:length(compounds)) {
        cpdID <- compounds[c]
        filename<-paste(substring(cpdID,5,nchar(cpdID)),"txt",sep=".")
        
        if (!file.exists(filename)){
          query_2 <- paste("http://rest.kegg.jp", "get", cpdID, sep=.Platform$file.sep)
        
          download.file(query_2, destfile=filename, method="auto",quiet=T)
        
          # then import the test.txt file, parse it, and store the Entry, Name, formula, exact_mass, 
          # just as what you have done for parsing HMDB file
          fileinfo<-paste(readLines(filename))
        
          for(line in fileinfo){
            if (length(grep("ENTRY",line))>0){
              cpdEntry<-substr(gsub(" ","",line),6,nchar(gsub(" ","",line))-8)
            }
            if(length(grep("NAME",line))>0){
              cpdName<-substr(gsub(" ","",line),5,nchar(gsub(" ","",line)))
            }
            if(length(grep("FORMULA",line))>0){
              cpdFormula<-substr(gsub(" ","",line),8,nchar(gsub(" ","",line)))
            }
            if(length(grep("EXACT_MASS",line))>0){
              cpdMass<-substr(gsub(" ","",line),11,nchar(gsub(" ","",line)))
            }
          }
        newcpd<-data.frame(Entry=cpdEntry,Name=cpdName,Formula=cpdFormula,Exact_Mass=cpdMass)
        results<-rbind(unique(results),newcpd)
        }
      }
    }
}


write.csv(unique(results),file="KEGGDatabaseResults.csv")
