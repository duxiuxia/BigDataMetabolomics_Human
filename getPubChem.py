import xml.etree.cElementTree as ET
import os

def process_buffer(buf):
    tnode = ET.fromstring(buf)
    # for child in tnode:
    #     print child.tag, child.attrib
    return tnode



#   Set the working directory with all the pubchem .xml files as path


path = '/Users/matthewdeitz/Desktop/PubChem'
for filename in os.listdir(path):
    if not filename.endswith('.xml'): continue
    #gets all the files and only checks for those is .msl ending


    fullname = os.path.join(path, filename)
# full file path

    Results = list()
    inputbuffer = ''
    with open(fullname,'rb') as inputfile:
        append = False
        for line in inputfile:
            if '<PC-Compound>' in line:#parsing by line for faster speed and read in
                inputbuffer = line
                append = True
            elif '</PC-Compound>' in line: #end of compund pass the lines into process buffer to string
                inputbuffer += line
                append = False
                current_compound=process_buffer(inputbuffer)
                inputbuffer = None
                del inputbuffer #probably redundant...

                PC_Compound_props = current_compound.find("PC-Compound_props")#values we are looking for in the compound

                PC_Compound = list(current_compound)

                current_compound_id = PC_Compound[0][0][0][0].text

                for PC_Info_Data in list(PC_Compound_props):
                    if PC_Info_Data[0][0][0].text == 'IUPAC Name':
                        current_IUPAC_name = PC_Info_Data[1][0].text
                    elif PC_Info_Data[0][0][0].text == 'InChIKey':
                        current_InChIKey = PC_Info_Data[1][0].text
                    elif PC_Info_Data[0][0][0].text == 'Mass':
                        current_exact_mass = PC_Info_Data[1][0].text
                    elif PC_Info_Data[0][0][0].text == 'Molecular Formula':
                        current_formula = PC_Info_Data[1][0].text
                    elif PC_Info_Data[0][0][0].text == 'Molecular Weight':
                        current_molecular_weight = PC_Info_Data[1][0].text
                        
#append the values to Results and then reset them for next itteration
                Results.append(current_compound_id+', '+current_IUPAC_name+', '+current_InChIKey+', '+current_exact_mass+', '+current_formula+', '+current_molecular_weight)
                current_compound_id = ""
                current_IUPAC_name = ""
                current_InChIKey = ""
                current_exact_mass = ""
                current_formula = ""
                current_molecular_weight = ""
            elif append:
                inputbuffer += line
    inputfile.close()
    #close the input file


# make corrresponding output file and write the results
    base, ext = os.path.splitext(filename)
    outputfilename = base + '.csv'
    outputfile = os.path.join(path, outputfilename)
    with open(outputfile,'wb') as f:
        for s in Results:
            f.write(s)
            f.write("\n")
    f.close()
