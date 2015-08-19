import xml.etree.cElementTree as ET

def process_buffer(buf):
    tnode = ET.fromstring(buf)
    # for child in tnode:
    #     print child.tag, child.attrib
    return tnode


inputbuffer = ''
with open('Compound_015775001_015800000.xml','rb') as inputfile:
    append = False
    for line in inputfile:
        if '<PC-Compound>' in line:
            inputbuffer = line
            append = True
        elif '</PC-Compound>' in line:
            inputbuffer += line
            append = False
            current_compound=process_buffer(inputbuffer)
            inputbuffer = None
            del inputbuffer #probably redundant...

            PC_Compound = list(current_compound)

            current_compound_id = PC_Compound[0][0][0][0].text
            print current_compound_id

            current_IUPAC_name = PC_Compound[5][10][1][0].text
            print current_IUPAC_name

            current_InChIKey = PC_Compound[5][12][1][0].text
            print current_InChIKey

            current_exact_mass = PC_Compound[5][14][1][0].text
            print current_exact_mass

            current_formula = PC_Compound[5][15][1][0].text
            print current_formula

            current_molecular_weight = PC_Compound[5][16][1][0].text
            print current_molecular_weight

        elif append:
            inputbuffer += line

