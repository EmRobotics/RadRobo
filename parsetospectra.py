# Imports
import os
import xml2dict
import matplotlib.pyplot as plt

#
parent_path = "/Users/emeline/Documents/NERS 491"

# Define path where xml files can be found
#path = "/Users/emeline/Documents/NERS 491/IsotopeID_n42/"
data = "IsotopeID_n42"

data_path = os.path.join(parent_path, data)
figure_path = os.path.join(parent_path, "IsotopeID_spectra")
os.mkdir(figure_path)

# Create list of the names of all xml files in this path
files = [x for x in os.listdir(data_path) if x.endswith('.n42')]

#
for file in files:

    #
    struct = xml2dict(file)

    #RadInstrumentData/RadMeasurement/Spectrum/ChannelData/Text
    #ans.RadInstrumentData.RadMeasurement.Spectrum.ChannelData.Text 
    struct_str = 

    string_list = struct_str.split() #include space in ()?

    counts = [int(el) for el in string_list]

    zeros = []

    #0 compression
    for i in range(len(counts)):
        count = counts[i]
      
        if count == 0:
            zeros.add([i, counts(i+1)])

    for i in range(len(zeros)):

        [ind, num] = zeros(i)
        
        counts.pop(ind+1)
        
        zeros = [0] * (num-1)
        counts.insert(ind, zeros)


    #create spectra

    plt.hist(counts)

    plt.xlabel('x - axis')
    plt.ylabel('y - axis')
    plt.title('My first graph!')
    plt.show()

    plt.savefig(file.replace('.n42','') + '.jpg')


    #save as file with corresponding filename




#print(files[0])