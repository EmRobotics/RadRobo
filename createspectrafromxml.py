import os
import xmlwrapper
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET 

def xmlread(file):
    # A file containing:
    # <XMLname attrib1="Some value">
    #   <Element>Some text</Element>
    #   <DifferentElement attrib2="2">Some more text</Element>
    #   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
    # </XMLname>
    #
    # Can be accessed by:
    # XMLname.Attributes.attrib1 = "Some value";
    # XMLname.Element.Text = "Some text";
    # XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
    # XMLname.DifferentElement{1}.Text = "Some more text";
    # XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
    # XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
    # XMLname.DifferentElement{2}.Text = "Even more text";
    #
    # Adapted from a MATLAB code originally
    # Written by W. Falkena, ASTI, TUDelft, 21-08-2010
    # Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
    # Added CDATA support by I. Smirnov, 20-3-2012
    # Modified by X. Mo, University of Wisconsin, 12-5-2012
    #
    # Python code written by Emeline Hanna 2022


    if os.path.exists(file):

        if file.endswith('.n42'):

            with open(file) as f:
                xmlstring = f.read()

            colon_bool = False 
            for element in range(0, len(xmlstring)):
                if (':' == xmlstring[element]):
                    colon_bool = True

            
            if colon_bool:
                ind = xmlstring.index('xmlns') + len('xmlns') 

                xmlstring = xmlstring[:ind] + ':H3D' + xmlstring[ind:]

            #tree = ET.parse(data)
            tree = ET.ElementTree(ET.fromstring(xmlstring))

            # Getting the parent tag of the xml document
            root = tree.getroot()

            return xmlwrapper.XmlParse(root)

        else:
            print('Not compatible file type.')
    else:
        print('File {} not found.'.format(file))


def createspectra(file,counts):
    counts = counts.split()

    counts = [int(el) for el in counts]
    zeros = []
    
    #0 compression
    for i in range(0,len(counts)):
        count = counts[i]
    
        if count == 0:
            zeros.append([i, counts[i+1]])     
    
    for i in range(len(zeros)-1,-1,-1):
    
        ind,num = zeros[i]
    
        counts = counts[0:ind] + [0]*num + counts[ind+2:]

    x = range(0,len(counts))

    slash_index = file.rindex('/')+1
    rename_file = file.replace('.n42','')

    #create spectra
    plt.step(x,counts)
    plt.xlabel('Channel number')
    plt.ylabel('Counts')
    plt.title(rename_file[slash_index:])
    #plt.show()
    plt.savefig(str(rename_file)+'.png')
    plt.clf()

#
parent_path = "/Users/emeline/Documents/NERS 491"

# Define path where xml files can be found
#path = "/Users/emeline/Documents/NERS 491/IsotopeID_n42/"
data = "IsotopeID_n42"

data_path = os.path.join(parent_path, data)
figure_path = os.path.join(parent_path, "IsotopeID_spectra")

if os.path.exists(figure_path):
    for root, dirs, files in os.walk(figure_path):
        for file in files:
            os.remove(os.path.join(root, file))  
else:
    os.mkdir(figure_path)     


# Create list of the names of all xml files in this path
files = [x for x in os.listdir(data_path) if x.endswith('.n42')]

#
for file in files:

    RadInstrumentData = xmlread(os.path.join(data_path,file))
    
    if not (RadInstrumentData.RadMeasurement.Spectrum.ChannelData.text == None):
    #i don't think this is actually doing anything

        counts = RadInstrumentData.RadMeasurement.Spectrum.ChannelData.text
        #print(counts)
        createspectra(os.path.join(figure_path,file),counts)
