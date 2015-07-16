# getKEGG.R

# Author: Matthew Deitz

# Editor: Xiuxia Du

# Started in June 2015.




rm(list=ls())






# ================================================
# set working directory in console using setwd()
# ================================================
setwd("C:/Users/matt/Desktop/KEGGMetabolites")



results <- data.frame(Entry=character(), Name=character(), Formula=character(), Exact_Mass=character())


masses <- seq(from=50, to=1500, by=0.01) 
# /find/compound/174.05/exact_mass   for 174.045 =< exact_mass =< 174.055, so use by=0.01 



fileName <- "warning_message.txt"
if (!file.exists(fileName)) {
    file.create(fileName)
}

error_message=""


for (i in 1:length(masses)) {
    current_mass <- as.character(masses[i])
    
    query_1 <- paste("http://rest.kegg.jp", "find", "compound", 
                     current_mass, "exact_mass", sep=.Platform$file.sep)
    
    tmp <- try(read.table(query_1), silent=T)
   
    if (inherits(tmp, 'try-error')) {
      #If an empty file do nothing
      
        warning_message <- paste("No compounds in KEGG for mass = ", current_mass, ".")
        error_message <- paste(error_message,warning_message,sep="\n")
    }
    else{
        dataIn <- read.table(query_1)
        
        compounds <- as.character(dataIn$V1)
        
        for (c in 1:length(compounds)) {
            cpdID <- compounds[c]
            
            filename <- paste(substr(x=cpdID, start=5, stop=nchar(cpdID)), "txt", sep=".")
            
            if (!file.exists(filename)) {
                query_2 <- paste("http://rest.kegg.jp", "get", cpdID, sep=.Platform$file.sep)
                
                # retrieve info for the specific compound and store the info in a .txt file
                download.file(query_2, destfile=filename, method="auto",quiet=T)
            }
            # parse the .txt file
            fileinfo <- paste(readLines(filename)) # Is paste necessary here?
                
            for (line in fileinfo) {
                  if (length(grep("ENTRY", line)) > 0) {
                      cpdEntry <- substr(x=gsub(pattern=" ", replacement="", x=line), 
                                         start=6, 
                                         stop=nchar(gsub(pattern=" ", replacement="", x=line))-8)
                  }
                    
                  if (length(grep("NAME",line)) > 0) {
                      cpdName<-substr(x=gsub(pattern=" ", replacement="", x=line),
                                      start=5,
                                      stop=nchar(gsub(pattern=" ", replacement="", x=line)))
                  }
                    
                  if (length(grep("FORMULA",line)) > 0) {
                      cpdFormula<-substr(x=gsub(pattern=" ", replacement="", x=line),
                                         start=8,
                                         stop=nchar(gsub(pattern=" ", replacement="", x=line)))
                  }
                    
                  if (length(grep("EXACT_MASS",line)) > 0) {
                      cpdMass<-substr(x=gsub(pattern=" ", replacement="", x=line),
                                      start=11,
                                      stop=nchar(gsub(pattern=" ", replacement="", x=line)))
                  }
                    
              }
                
              newcpd <- data.frame(Entry=cpdEntry, Name=cpdName, Formula=cpdFormula, Exact_Mass=cpdMass)
                
              results <- rbind(unique(results), newcpd)
        }
    }
}


cat(error_message, file="ErrorFile.txt", sep="\n", append=T)
write.csv(unique(results), file="KEGGDatabaseResults.csv")
