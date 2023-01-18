import xml.etree.ElementTree as ET
import os

# Object that abstracts the dict['key'] operation using the '.' operator
class Attributes:
    def __init__(self, dict):
        self.dict = dict
    
    def __getattr__(self, key):
        return self.dict[key]

# Object that abstracts etree.findall, etree.attrib, and etree.text 
# Attributes in one place using the '.' operator
class XmlParse:
    def __init__(self, root):
        self.root = root
        self.xmlns = root.tag[:root.tag.find('}')+1]
        
    def __getattr__(self, key):
        # Text and Attrib
        if (key == 'text'):
            return self.root.text
        elif (key == 'Attributes'):
            return Attributes(self.root.attrib)
        
        #Other Children
        to_search_for = self.xmlns + key
        new_etrees = self.root.findall(to_search_for)


        if len(new_etrees) == 0:
            raise Exception('{} has no attribute {}'.format(get_tag(self.root),key))
        elif len(new_etrees) == 1:
            return XmlParse(new_etrees[0]) 
        else:
            arr = [0]
            for etree in new_etrees:
                arr.append(XmlParse(etree))
            return arr

def get_tag(child):
    return child.tag[child.tag.find('}')+1:]      


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

            #return xmlwrapper.XmlParse(root)
            return XmlParse(root)

        else:
            print('Not compatible file type.')
    else:
        print('File {} not found.'.format(file))
    