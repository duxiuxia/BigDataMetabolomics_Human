# Author: Matthew Deitz

# Editor: Xiuxia Du

# Started in June 2015




if (!is.element("XML", installed.packages()[,1])) {
    install.packages("XML")
}


library(XML)



metaboliteList<-data.frame(Accession=character(0),
                           Name=character(0),
                           Chemical_Formula=character(0),
                           Avg_Molecular_Weight=character(0),
                           Monisotopic_weight=character(0))




for (infile in dir(getwd(),pattern="*.xml")){
  data<-xmlParse(infile)
  
  xml_data<-xmlToList(data)
  
  metaboliteList<-rbind(metaboliteList,
                        data.frame(Accession=xml_data[["accession"]],
                                   Name=xml_data[["name"]],
                                   Chemical_Formula=xml_data[["chemical_formula"]],
                                   Avg_Molecular_Weight=xml_data[["average_molecular_weight"]],
                                   Monisotopic_weight=xml_data[["monisotopic_moleculate_weight"]]))
}


write.csv(metaboliteList,"HMDBMetaboliteList.csv")