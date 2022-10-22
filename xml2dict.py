import os
import xml.etree.ElementTree as ET
import json

list_tree = []

#returns a dict
def recursive_dict(node):
    arr = []
    for child in node:        
        
        if len(child) > 0 and not any(get_tag(child) in i for i in arr):
            arr.append((get_tag(child), dict(recursive_dict(child))))
        elif len(child) > 0:
            print('help')
        else:
            arr.append((get_tag(child), child.text))
    
    return arr
    


def recurse_tree_print(root):
    dict = {get_tag(root): {}}
    for child in root:
        dict[get_tag(root)] = recurse_tree_print(child)
        #tag = get_tag(child)
        #for i in range(level):
        #    print("\t", end="")
        
        #dict[tag]
        #print(tag, child.attrib,len(child))

        #ist_tree.append([tag, level, len(child)])

        #recurse_tree(child,level+1)

    return dict  




def recurse_tree(root):
    #parent = get_tag(root)
    dict = {}
    
    for child in root: 
        tag = get_tag(child)

        if len(child) > 0:
            #parent = tag
            dict[tag] = ''
            recurse_tree(child)
        else:
            #inner dict instead of vec
            vec = [get_tag(child), child.text]
            dict.update({get_tag(root): vec})

def recurse_tree(root):
    dict = {}
    
    for child in root: 
        tag = get_tag(child)

        if len(child) > 0:
            dict[tag] = ''
            recurse_tree(child)
        else:
            inner_dict = {get_tag(child): child.text}
            dict.update({get_tag(root): inner_dict})

            return(dict)



           
            

#
    #    tag_list.append(tag)
    #    recurse_tree(child)
    #    tag = get_tag(child)
    #    tag_list.append(tag)
    #
    #print(tag_list)
#
    #dict[get_tag(root)] = tag_list
#
    #for child in root: 
    #    recurse_tree(child)
#
 #   return dict
    
#def print_dictionary(dictionary, level = 0):
#    if (len(dictionary.keys()) > 1):
#        for key,value in dictionary.items():
#            for i in range(level):
#                print("\t", end="")
#    	        print('key: ', key, ', value: ', value)   
#
#        print_dictionary(value, level+1)
#  


    #list_tree.find
    # find if two+ tags in the list at the same level are the same

    # put all of those attributes in a single dictionary

def get_tag(child):
    return child.tag[child.tag.find('}')+1:]

def xmlread(file):
    
    # Passing the path of the xml file to enable parsing process
    tree = ET.parse(file)
    
    # getting the parent tag of the xml document
    root = tree.getroot()
    #print(get_tag(root))
    d = {get_tag(root): dict(recursive_dict(root))}
    print(json.dumps(d, indent = 4))



# to do
# add .attributes and .text
# fix overwritting if same key


    #print(d['RadInstrumentData'].keys())
    #print(d.RadInstrumentData.RadInstrumentInformation.RadInstrumentManufacturerName)


    #print(root.__dict__)
    #dict = {}
    #dict = recurse_tree(dict,root)
    #dict = {get_tag(root), dict} 
    #print_dictionary(list_tree)

    #create_dictionary(root,list_tree)



#def xml2struct(file):
    #Convert xml file into a python structure

    # A file containing:
    # <XMLname attrib1="Some value">
    #   <Element>Some text</Element>
    #   <DifferentElement attrib2="2">Some more text</Element>
    #   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
    # </XMLname>
    #
    # Will produce:
    # s.XMLname.Attributes.attrib1 = "Some value";
    # s.XMLname.Element.Text = "Some text";
    # s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
    # s.XMLname.DifferentElement{1}.Text = "Some more text";
    # s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
    # s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
    # s.XMLname.DifferentElement{2}.Text = "Even more text";
    #
    # Please note that the following characters are substituted
    # '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
    #
    # Written by W. Falkena, ASTI, TUDelft, 21-08-2010
    # Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
    # Added CDATA support by I. Smirnov, 20-3-2012
    #
    # Modified by X. Mo, University of Wisconsin, 12-5-2012
    #
    # Turned into python from MATLAB code by Emeline Hanna 2022
file = '/Users/emeline/Documents/NERS 491/IsotopeID_n42/am241_200.n42'

if os.path.exists(file):
    
    #this needs to be better
    if file.endswith('.n42'):
        xmlread(file)
    else:
        print('Not compatible file type.')
else:
    print('File not found.')
  
