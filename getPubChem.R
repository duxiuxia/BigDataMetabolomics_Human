# getPubChem.R


# Written by Xiuxia Du and Matthew Deitz.


# Started in July 2015.





rm(list=ls())




# ===========================================================
# set working directory in CMD Line by passing dir as variable
# ===========================================================
args <- commandArgs(trailingOnly = TRUE)
workingDir <- as.character(args[1])




setwd(workingDir)

MainDirectory <- getwd()





if (!is.element("XML", installed.packages()[,1])) {
    install.packages("XML")
}

library(XML)



if (!is.element("R.utils", installed.packages()[,1])) {
    install.packages("R.utils")
}

library(R.utils)


ZipFiles <- list.files(path=getwd(), pattern="*.gz")
#get all of the gzipped files


dir.create(file.path(MainDirectory, "Temp", fsep=.Platform$file.sep), showWarnings=T)
pathToTemp <- paste(MainDirectory, "Temp", sep=.Platform$file.sep)
#Create a temp dir to store unzipped files and then remove
#record the path to that temp dir

all_compounds <- data.frame(PubChem_ID=numeric(), 
                            IUPAC_Name=character(), 
                            InChIKey=character(), 
                            Formula=character(),
                            Exact_Mass=numeric(),
                            Molecular_Weight=numeric()
)


for(zipfile in ZipFiles){
#for every zipped file in the directory
    
    
    file.copy(from=paste(MainDirectory, zipfile, sep=.Platform$file.sep), 
              to=paste(pathToTemp, zipfile, sep=.Platform$file.sep))
    #copy the .gz file to Temp directory
    
    setwd(pathToTemp)
    gunzip(zipfile)
    #go to that temp dir and unzip file
    fileName <- substr(zipfile,1,nchar(zipfile)-3)
    #unzip the files and get the name from that file without the .gz
    
    
    
    inFile <- fileName
    
    dataIn <- xmlParse(inFile)
    
    xmlData <- xmlToList(dataIn)
    #get the file and parse it to a list
    
    
    for (i in 1:length(xmlData)) {
        current_compound <- xmlData[[i]]
        
        if(is.recursive(current_compound)){
            current_compound_id <- as.numeric(current_compound$"PC-Compound_id"$"PC-CompoundType"$"PC-CompoundType_id"$"PC-CompoundType_id_cid")
        }else{
            next
        }
        
    

  
            current_compound_properties <- current_compound[["PC-Compound_props"]]
            
            #First if statement checks for index out of bounds error
            #nested if statement checks for corresponding value ex. IUPAC Name
            if (class(try(grepl(pattern="IUPAC Name", x=current_compound_properties[[11]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label"),silent=TRUE))!="try-error") {
                if(grepl(pattern="IUPAC Name", x=current_compound_properties[[11]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label")){
                    current_IUPAC_name <- current_compound_properties[[11]]$"PC-InfoData_value"$"PC-InfoData_value_sval"
                }else{
                    current_IUPAC_name <- NA
                }
                
            }else{
                current_IUPAC_name <- NA
            }
            
            
            
            
            if (class(try(grepl(pattern="InChIKey", x=current_compound_properties[[13]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label"),silent=TRUE))!="try-error"){
                if(grepl(pattern="InChIKey", x=current_compound_properties[[13]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label")){
                    current_InChIKey <- current_compound_properties[[13]]$"PC-InfoData_value"$"PC-InfoData_value_sval"
                }else{
                    current_InChIKey <- NA
                }
            }else{
                current_InChIKey <- NA
            }
            
            
            
            
            if ((class(try(grepl(pattern="Mass", x=current_compound_properties[[15]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label" ),silent=TRUE))!="try-error") && 
                    (class(try(grepl(pattern="Exact", x=current_compound_properties[[15]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_name" ),silent=TRUE))!="try-error")) {
                
                if(grepl(pattern="Mass", x=current_compound_properties[[15]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label" )&&
                       grepl(pattern="Exact", x=current_compound_properties[[15]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_name" )){
                    
                    current_exact_mass <- as.numeric(current_compound_properties[[15]]$"PC-InfoData_value"$"PC-InfoData_value_fval")
                }else{
                    current_exact_mass <- NA
                }
            }else{
                current_exact_mass <- NA
            }
            
            
            
            
            
            if (class(try(grepl(pattern="Molecular Formula", x=current_compound_properties[[16]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label" ),silent=TRUE))!="try-error"){
                if(grepl(pattern="Molecular Formula", x=current_compound_properties[[16]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label")){
                    current_formula <- current_compound_properties[[16]]$"PC-InfoData_value"$"PC-InfoData_value_sval"
                }else{
                    current_formula <- NA
                }
            }else{
                current_formula <- NA
            }
            
            
            
            
            
            if (class(try(grepl(pattern="Molecular Weight", x=current_compound_properties[[17]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label" ) ,silent=TRUE))!="try-error"){
                if(grepl(pattern="Molecular Weight", x=current_compound_properties[[17]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label" )){
                    current_molecular_weight <- as.numeric(current_compound_properties[[17]]$"PC-InfoData_value"$"PC-InfoData_value_fval")
                }else{
                    current_molecular_weight <- NA
                }
            }else{
                current_molecular_weight <- NA
            }
            
            
            new_compound <- data.frame(PubChem_ID=current_compound_id, 
                                       IUPAC_Name=current_IUPAC_name, 
                                       InChIKey=current_InChIKey, 
                                       Formula=current_formula,
                                       Exact_Mass=current_exact_mass,
                                       Molecular_Weight=current_molecular_weight
            )
            
            all_compounds <- rbind(all_compounds, new_compound)
        }
        
        out_file_name <- paste(substr(x=fileName, start=1, stop=nchar(fileName)-4), ".csv", sep="")
        
        write.csv(x=all_compounds, file=out_file_name)

        
        file.remove(fileName)
        #remove the unzipped file for storage space
        
        setwd(MainDirectory)
        
    }
