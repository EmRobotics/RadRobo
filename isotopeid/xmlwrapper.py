import xml.etree.ElementTree as ET

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



    