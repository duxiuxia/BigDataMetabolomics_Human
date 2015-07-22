
#Needs xlsx package
library(xlsx)

###########################################
#   Set working directory for Databases
###########################################
setwd("C:/Users/matt/Desktop/MetaboliteDB/")
HMDBdb<-read.csv("HMDBMetaboliteList.csv")
KEGGdb<-read.csv("KEGGDatabaseResults.csv")
NISTdb<-read.csv("NISTDatabaseResults.csv")


#####PubChemdb<-read.csv("")
#----- Temporarily removed



#####################################
#   Set Directory for unknowns to be run
#####################################

setwd("C:/Users/matt/Desktop/XCMSResultsNutritionStudy/")
unknowns <- read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName = "PV1")



PPMTolerance=10


results=data.frame(MetaboliteMZ=character(),DataBase=character(),FoundID=character())


HMDBList=sort(HMDBdb$Monisotopic_weight)
KEGGList=sort(KEGGdb$Exact_Mass)
NISTList=sort(NISTdb$ExactMass)

#PubChemList=sort(PubChemdb$V1)


#gets and sorts the list of mz values


for(metabolite in sort(unknowns$query_mz)){
    Flag = 0
    for(mzvalue in HMDBList){      
        if(((mzvalue-metabolite)/metabolite*1000000)<-PPMTolerance){
            HMDBList=HMDBList[HMDBList != mzvalue]
            #reupdate list in HMDB for future values since they are larger
        }
        if((((mzvalue-metabolite)/metabolite*1000000)<PPMTolerance)&(((mzvalue-metabolite)/metabolite*1000000)>-PPMTolerance)){
            results =rbind(results,data.frame(MetaboliteMZ=metabolite,DataBase="HMDB",FoundID=paste(HMDBdb[which(HMDBdb$V1==mzvalue)],sep=" ")))
            Flag = 1
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
            Flag = 1
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
            Flag = 1
            #update results for a found metabolite
        }
        if (((mzvalue-metabolite)/metabolite*1000000)>PPMTolerance){
            break #get out of for loop because all other values are larger
            #go to next database
        }
    }
  
  
#    for(mzvalue in PubChemList){
#        if(((mzvalue-metabolite)/metabolite*1000000)<-PPMTolerance){
#            PubChemList=PubChemList[PubChemList != mzvalue]
#            #reupdate list in HMDB for future values since they are larger
#        }    
#        if((((mzvalue-metabolite)/metabolite*1000000)<PPMTolerance)&(((mzvalue-metabolite)/metabolite*1000000)>-PPMTolerance)){
#            results =rbind(results,data.frame(MetaboliteMZ=metabolite,DataBase="PubChem",FoundID=paste(PubChemdb[which(PubChemdb$V1==mzvalue)],sep=" ")))
#            Flag = 1
#            #update results for a found metabolite
#        }
#        if (((mzvalue-metabolite)/metabolite*1000000)>PPMTolerance){
#            break #get out of for loop because all other values are larger
#            #go to next database
#        }
#    }
    
    #Determine if no match was found then if not record it as unknown
    if(Flag != 1){
        results = rbind(results,data.frame(MetaboliteMZ=metabolite,DataBase="NA",FoundID="NA"))
    }
}

write.csv("IdentifiedMetabolites.csv",results)
