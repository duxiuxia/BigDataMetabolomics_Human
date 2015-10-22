# Author: Matthew Deitz

# Editor: Xiuxia Du

# Started in June 2015


#Run this program in the same directory as the HMDB XML database file

if (!is.element("XML", installed.packages()[,1])) {
    install.packages("XML")
}


library(XML)


#data frame of the results
metaboliteList<-data.frame(Accession=character(0),
                           Name=character(0),
                           Chemical_Formula=character(0),
                           Avg_Molecular_Weight=character(0),
                           Monisotopic_weight=character(0))



#for the xml file in your current working directory
for (infile in dir(getwd(),pattern="*.xml")){
  data<-xmlParse(infile)
  #parse the file to a list
  xml_data<-xmlToList(data)
  
  metaboliteList<-rbind(metaboliteList,
                        data.frame(Accession=toString(xml_data[["accession"]]),
                                   Name=toString(xml_data[["name"]]),
                                   Chemical_Formula=toString(xml_data[["chemical_formula"]]),
                                   Avg_Molecular_Weight=toString(xml_data[["average_molecular_weight"]]),
                                   Monisotopic_weight=toString(xml_data[["monisotopic_moleculate_weight"]])))
}


write.csv(metaboliteList,"HMDBMetaboliteList.csv")
