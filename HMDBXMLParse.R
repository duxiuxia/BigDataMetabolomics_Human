library(XML)
setwd("C:/Users/matt/Desktop/HMDBMetabolites")

metaboliteList<-data.frame(Accession=character(0),Name=character(0),Chemical_Formula=character(0),Avg_Molecular_Weight=character(0),Monisotopic_weight=character(0))
for (infile in dir(getwd(),pattern="*.xml")){
  data<-xmlParse(infile)
  xml_data<-xmlToList(data)
  metaboliteList<-rbind(metaboliteList,data.frame(Accession=xml_data[[4]],Name=xml_data[[6]],Chemical_Formula=xml_data[[9]],Avg_Molecular_Weight=xml_data[[10]],Monisotopic_weight=xml_data[[11]]))
}
write.csv(metaboliteList,"HMDBMetaboliteList.csv")