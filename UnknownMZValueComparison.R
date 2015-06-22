PPMTolerance=10

library(xlsx)
setwd("C:/Users/matt/Desktop/XCMSResultsNutritionStudy")
data1=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV1")
data2=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV2")
data3=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV3")
data4=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV4")
data5=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV5")
data6=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV6")
data7=read.xlsx("UnknownMetabolitesMZValues.xlsx",sheetName="NV7")
#Reads in the data from the worksheet

mzvalues = rbind(data1,data2,data3,data4,data5,data6,data7)
#combines all the rows into one dataset

mzvalues = mzvalues[ with(mzvalues,order(Query_mz)),]
#Orders the mzvalues by the Query value

m=matrix(0,ncol=1,nrow=length(mzvalues))

mzvalues = cbind(mzvalues,m)
#adds 0 to the data frame for each value

for (i in 1:length(mzvalues[,1])){
  value1 = mzvalues[i,1]
  tempFlag=T
  j=i+1
  #Goes from the begining of the array to the end
  #tempFlag for determing PPM tolerance in the while loop
  #j is the comparison value to i
  while(tempFlag & (j<length(mzvalues[,1]))){
    #while true and < continue looking for similar values
    value2 = mzvalues[j,1]
    #the comparison value
    if ((value2-value1)/value2*1000000>PPMTolerance){
      tempFlag=F
      #if the value is greater than the ppm tolerance set the flag to false
    }
    else{
      #else if within PPM tolernce
      mzvalues[i,2]=mzvalues[i,2]+1
      mzvalues[j,2]=mzvalues[j,2]+1
      j=j+1
      #increase the value for determining similar values
      #by one for each so we don't have to go backwards
      #go to the next j value
    }
  }
}
write.xlsx(mzvalues,"C:/Users/matt/Desktop/NegativeResults.xlsx")
