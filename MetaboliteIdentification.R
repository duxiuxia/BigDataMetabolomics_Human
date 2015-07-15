setwd("C:/Users/matt/Desktop/MetaboliteDB")
HMDBdb<-read.table("")
KEGGdb<-read.table("")
NISTdb<-read.table("")
PubChemdb<-read.table("")

setwd("C:/Users/matt/Desktop")
unknowns <- read.table("")

PPMTolerance=10
results=data.frame(MetaboliteMZ=character(),DataBase=character(),FoundID=character())

HMDBList=sort(HMDBdb$V1)
KEGGList=sort(KEGGdb$V1)
NISTList=sort(NISTdb$V1)
PubChemList=sort(PubChemdb$V1)


for(metabolite in sort(unknowns)){
  
  for(mzvalue in HMDBList){
    if(((mzvalue-metabolite)/metabolite*1000000)<-PPMTolerance){
      HMDBList=HMDBList[HMDBList != mzvalue]
      #reupdate list in HMDB for future values since they are larger
    }    
    if((((mzvalue-metabolite)/metabolite*1000000)<PPMTolerance)&(((mzvalue-metabolite)/metabolite*1000000)>-PPMTolerance)){
      results =rbind(results,data.frame(MetaboliteMZ=metabolite,DataBase="HMDB",FoundID=paste(HMDBdb[which(HMDBdb$V1==mzvalue)],sep=" ")))
      #update results for a found metabolite
    }
    if (((mzvalue-metabolite)/metabolite*1000000)>PPMTolerance){
      break #get out of for loop because all other values are larger
      #go to next database
    }
  }
  
  
  for(mzvalue in KEGGList){
    if(((mzvalue-metabolite)/metabolite*1000000)<-PPMTolerance){
      KEGGList=KEGGList[KEGGList != mzvalue]
      #reupdate list in HMDB for future values since they are larger
    }    
    if((((mzvalue-metabolite)/metabolite*1000000)<PPMTolerance)&(((mzvalue-metabolite)/metabolite*1000000)>-PPMTolerance)){
      results =rbind(results,data.frame(MetaboliteMZ=metabolite,DataBase="KEGG",FoundID=paste(KEGGdb[which(KEGGdb$V1==mzvalue)],sep=" ")))
      #update results for a found metabolite
    }
    if (((mzvalue-metabolite)/metabolite*1000000)>PPMTolerance){
      break #get out of for loop because all other values are larger
      #go to next database
    }
  }
  
  
  for(mzvalue in NISTList){
    if(((mzvalue-metabolite)/metabolite*1000000)<-PPMTolerance){
      NISTList=NISTList[NISTList != mzvalue]
      #reupdate list in HMDB for future values since they are larger
    }    
    if((((mzvalue-metabolite)/metabolite*1000000)<PPMTolerance)&(((mzvalue-metabolite)/metabolite*1000000)>-PPMTolerance)){
      results =rbind(results,data.frame(MetaboliteMZ=metabolite,DataBase="NIST",FoundID=paste(NISTdb[which(NISTdb$V1==mzvalue)],sep=" ")))
      #update results for a found metabolite
    }
    if (((mzvalue-metabolite)/metabolite*1000000)>PPMTolerance){
      break #get out of for loop because all other values are larger
      #go to next database
    }
  }
  
  
  for(mzvalue in PubChemList){
    if(((mzvalue-metabolite)/metabolite*1000000)<-PPMTolerance){
      PubChemList=PubChemList[PubChemList != mzvalue]
      #reupdate list in HMDB for future values since they are larger
    }    
    if((((mzvalue-metabolite)/metabolite*1000000)<PPMTolerance)&(((mzvalue-metabolite)/metabolite*1000000)>-PPMTolerance)){
      results =rbind(results,data.frame(MetaboliteMZ=metabolite,DataBase="PubChem",FoundID=paste(PubChemdb[which(PubChemdb$V1==mzvalue)],sep=" ")))
      #update results for a found metabolite
    }
    if (((mzvalue-metabolite)/metabolite*1000000)>PPMTolerance){
      break #get out of for loop because all other values are larger
      #go to next database
    }
  }
}
