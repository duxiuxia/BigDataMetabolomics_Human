if (!is.element("XML", installed.packages()[,1])) {
    install.packages("XML")
}


library(XML)

infile <- "Compound_033850001_033875000.xml"

dataIn<-xmlParse(infile)

xml_data<-xmlToList(dataIn)
