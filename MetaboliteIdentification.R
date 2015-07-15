setwd("C:/Users/matt/Desktop/HMDBMetabolites")
HMDBdb<-read.table("")
setwd("C:/Users/matt/Desktop/KEGGMetabolites")
KEGGdb<-read.table("")
setwd("C:/Users/matt/Desktop/NISTMetabolites")
NISTdb<-read.table("")
setwd("C:/Users/matt/Desktop/PubChemMetaboblites")
PubChemdb<-read.table("")

setwd("C:/Users/matt/Desktop")
unknowns <- read.table("")

PPMTolerance=10

for(metabolite in unknowns){
  for(mzvalue in HMDBdb$V1){
    if(<PPMTolerance){
      results #Identified
    }
  }
  for(mzvalue in KEGGdb$V1){
    if(<PPMTolerance){
      results #Identified
    }
  }
  for(mzvalue in NISTdb$V1){
    if(<PPMTolerance){
      results #Identified
    }
  }
  for(mzvalue in PubChemdb$V1){
    if(<PPMTolerance){
      results #Identified
    }
  }
}