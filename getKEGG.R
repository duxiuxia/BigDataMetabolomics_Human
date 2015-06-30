# getKEGG.R

# Written in June 2015.


masses <- seq(from=50, to=1500, by=0.01)

for (i in 1:length(mass_range)) {
    current_mass <- as.character(masses[i])
    
    query_1 <- paste("http://rest.kegg.jp", "find", "compound", 
                   current_mass, "exact_mass", sep=.Platform$file.sep)
    
    dataIn <- read.table(query_1)
    
    compounds <- as.character(dataIn$V1)
    
    for (c in 1:length(compounds)) {
        cpdID <- compounds[c]
        
        query_2 <- paste("http://rest.kegg.jp", "get", cpdID, sep=.Platform$file.sep)
        
        download.file(query_2, destfile="test.txt", method="auto")
        
        # then import the test.txt file, parse it, and store the Entry, Name, formula, exact_mass, 
        # just as what you have done for parsing HMDB file
        
    }
}



