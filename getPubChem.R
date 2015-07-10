# getPubChem.R


# Written by Xiuxia Du.


# Started in July 2015.





rm(list=ls())




# ===========================================================
# set working directory in RStudio console using setwd() 
# ===========================================================





if (!is.element("XML", installed.packages()[,1])) {
    install.packages("XML")
}

library(XML)




all_compounds <- data.frame(IUPAC_Name=character(), 
                           InChIKey=character(), 
                           Formula=character(),
                           Exact_Mass=numeric(),
                           Molecular_Weight=numeric()
                           )


inFile <- "Compound_033850001_033875000.xml"

dataIn <- xmlParse(inFile)

xmlData <- xmlToList(dataIn)


for (i in 1:length(xmlData)) {
    current_compound <- xmlData[[i]]
    
    current_compound_properties <- current_compound[["PC-Compound_props"]]
    

    if (grepl(pattern="IUPAC Name", x=current_compound_properties[[11]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label")) {
        current_IUPAC_name <- current_compound_properties[[11]]$"PC-InfoData_value"$"PC-InfoData_value_sval"
    }else{
        current_IUPAC_name <- NA
    }
    
    
    if (grepl(pattern="InChIKey", x=current_compound_properties[[13]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label")) {
        current_InChIKey <- current_compound_properties[[13]]$"PC-InfoData_value"$"PC-InfoData_value_sval"
    }else{
        current_InChIKey <- NA
    }
        
    if ( (grepl(pattern="Mass", x=current_compound_properties[[15]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label" )) && 
             (grepl(pattern="Exact", x=current_compound_properties[[15]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_name" )) 
    ) {
        current_exact_mass <- as.numeric(current_compound_properties[[15]]$"PC-InfoData_value"$"PC-InfoData_value_fval")
    }else{
        current_exact_mass <- NA
    }
    
    
    if ( grepl(pattern="Molecular Formula", x=current_compound_properties[[16]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label" ) ) {
        current_formula <- current_compound_properties[[16]]$"PC-InfoData_value"$"PC-InfoData_value_sval"
    }else{
        current_formula <- NA
    }
    
    
    if ( grepl(pattern="Molecular Weight", x=current_compound_properties[[17]]$"PC-InfoData_urn"$"PC-Urn"$"PC-Urn_label" ) ) {
        current_molecular_weight <- as.numeric(current_compound_properties[[17]]$"PC-InfoData_value"$"PC-InfoData_value_fval")
    }else{
        current_molecular_weight <- NA
    }
    
    
    new_compound <- data.frame(IUPAC_Name=current_IUPAC_name, 
                               InChIKey=current_InChIKey, 
                               Formula=current_formula,
                               Exact_Mass=current_exact_mass,
                               Molecular_Weight=current_molecular_weight
                               )
    
    all_compounds <- rbind(all_compounds, new_compound)
}



index <- regexpr(pattern=".xml", text=inFile)
out_file_name <- paste(substr(x=inFile, start=1, stop=index-1), ".csv", sep="")

write.csv(x=all_compounds, file=out_file_name)