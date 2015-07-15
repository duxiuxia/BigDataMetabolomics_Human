setwd("C:/Users/matt/Desktop/MetaboliteDB")
HMDBdb<-read.table("")
KEGGdb<-read.table("")
NISTdb<-read.table("")
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