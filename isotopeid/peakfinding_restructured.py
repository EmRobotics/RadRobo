import os
import xmlwrapper
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks,savgol_filter

# Create class for describing spectra
class Spectra:
    def __init__(self, file, counts,compressed):

        # File name
        self.iso_name = file[file.rindex('/')+1:file.rindex('_')]

        # Uncompress zeros
        if compressed:
            self.counts = expand_zeros(counts)
        
        # Extract numbers from string
        else:
            self.counts = [float(el) for el in counts.split(" ")[:-1]]
    
        # Define arrays of important parameters
        self.wdt_arr = np.array([max(i,2) for i in np.arange(0,len(self.counts),1)*0.005])
        self.height_vec = find_height(self.counts)
        self.prom = 0.01*max(self.counts)

# Expand the zeros in compressed spectra
def expand_zeros(counts):
    counts = counts.split()

    counts = [int(el) for el in counts]
    zeros = []
    
    # Zeros compression
    for i in range(0,len(counts)):
        count = counts[i]
    
        if count == 0:
            zeros.append([i, counts[i+1]])     
    
    for i in range(len(zeros)-1,-1,-1):
    
        ind,num = zeros[i]
        counts = counts[0:ind] + [0]*num + counts[ind+2:]

    return counts   

# Determine a height threshold
def find_height(counts):

    # Find standard deviation of specta
    std_dev = []
    for i in range(1,len(counts)-1):
        std_dev.append(np.sqrt(counts[i-1] + 4*counts[i] + counts[i+1]))
    std_dev.append(0)
    std_dev.insert(0,0)  

    # Filter counts using Savitzky–Golay filter
    counts_filtered = savgol_filter(counts, 50, 2)

    # Calculate height threshold
    height = np.array([int(max(std_dev[i]+0.98*counts_filtered[i],1)) for i in range(len(std_dev))])

    return height 

# Determine isotopes given the name
def determine_name(iso_name):

    iso_name = iso_name.lower()

    # Redefine isotope name to follow isotope-number naming convention
    if "heu" in iso_name:
        iso_name = iso_name.replace("heu","u235u-238")
    elif "du" in iso_name:
        iso_name = iso_name.replace("du","u235u-238")
    elif "wgpu" in iso_name:
        iso_name = iso_name.replace("wgpu","pu239")

    # Split up name into list
    old_nm = iso_name
    for x in range(len(old_nm)):
        if iso_name[x].isalpha() and iso_name[x+1].isdigit():
            iso_name = iso_name[:x+1] + '-' + iso_name[x+1:]

        elif iso_name[x].isdigit() and iso_name[x+1].isalpha():
            iso_name = iso_name[:x+1] + ',' + iso_name[x+1:]

    # Redefine isotope name due to isotope database
    if "tc-99" in iso_name:
        iso_name += "m"

    return iso_name.split(',')

# Find true isotope peaks
def isotope_peaks(iso_name,isotope_dict):

    # Determine true isotope
    iso_list = determine_name(iso_name)

    peaks = ""
    prob = ""

    # Loop through isotope dictionary
    for x in range(len(iso_list)):
        for key,val in isotope_dict.items():
            if key[3:] == iso_list[x]:

                # Add peaks and probability to list
                if x > 0:
                    peaks += "," + val['peaks']
                    prob += "," + val['probability']
                else:
                    peaks +=  val['peaks']
                    prob +=  val['probability']
   
    # Return peaks, probability
    return [float(i) for i in peaks.split(',')],[float(i) for i in prob.split(',')]

# Score isotope identification
def score_isotope(guess_name,true_name):
    
    score = 0

    # Spectra with isotopes
    if not true_name[0]=="":
        if len(guess_name) > len(true_name):
            return 1/len(guess_name)

        # Compare true isotope and guess isotopes
        for i in range(len(true_name)):
            for j in range(len(guess_name)):
                if true_name[i] == guess_name[j]:
                    score += 1        

        return score/len(true_name)
    
    # Background specta
    else:
        if len(guess_name) > 0:
            
            # Expected background isotopes
            if len(guess_name) == 1 and guess_name[0] == 'ra-226':
                return 1
            else:
                return 0
        else:
            return 1

# Define xrays from file
def xray_list(isotopes):

    xrays = np.array([])
    for key,val in isotopes.items():
        if val['xray']:
            for i in val['xray'].split(','):
                xrays = np.append(xrays,float(i))

    # Return xray peak list
    return xrays

# Guess isotope based on peaks found
def determine_isotope(all_peaks,all_prominences,isotopes,iso_name):

    # Initialize variables
    score = np.zeros(len(isotopes))
    name = np.array([])
    index = 0
    sim_peaks = [93, 186]

    # Loop through all possible isotopes
    for key,val in isotopes.items():

        # Initialize score boost for each isotope
        score_boost_og = 1

        # Check for x-ray peaks
        if val['xray']:
            xrays = np.array([])

            for i in val['xray'].split(','):
                xrays = np.append(xrays,float(i))

            peaks = np.array([])
            prominences = np.array([])

            for x in range(len(all_peaks)):

                # Only include gamma peaks
                if not any(abs(xray - all_peaks[x]) < 2 for xray in xrays):
                    peaks = np.append(peaks,all_peaks[x])
                    prominences = np.append(prominences,all_prominences[x])
                
                # If x-ray peaks are found, increase score
                else:
                    score_boost_og *= 2 

        else:
            peaks = all_peaks
            prominences = all_prominences 

        # Define isotope database parameters
        name = np.append(name,key)
        iso_peaks = np.array([float(i) for i in val['peaks'].split(',')])
        prob_peaks = np.array([float(i) for i in val['probability'].split(',')])

        # Define sorted arrays for isotope database
        sorted_index_array = np.argsort(prob_peaks)
        sorted_iso_peaks = np.flip(iso_peaks[sorted_index_array],0)
        sorted_prob_peaks = np.flip(prob_peaks[sorted_index_array],0)

        # Define sorted arrays for unkown isotope
        sorted_ind_array = np.argsort(prominences)
        sorted_peaks = np.flip(peaks[sorted_ind_array],0)
        sorted_prom_peaks = np.flip(prominences[sorted_ind_array],0)

        # Loop through peaks and compare
        for i in range(len(sorted_peaks)):
            for j in range(len(sorted_iso_peaks)):

                peak = sorted_peaks[i]
                iso_peak = sorted_iso_peaks[j]

                # If a peak is within error of isotope peak
                if abs(peak - iso_peak) < 2.5:

                    # Factor based on distance between peaks
                    dist_factor = 1-(abs(peak - iso_peak)/2)*0.1

                    # Parameters based on peak iteration
                    iso_peak_prob = sorted_prob_peaks[j]
                    iso_peak_prob_frac = iso_peak_prob/max(prob_peaks)
                    prom_peak_frac = sorted_prom_peaks[i]/max(sorted_prom_peaks)

                    # Ensure that score boost from x-ray peaks is not too large
                    score_boost = min(5,score_boost_og) 

                    # Increase score based on prominances and emission probability 
                    if abs(prom_peak_frac - iso_peak_prob_frac) < .25:
                        score_boost *= 10
                    elif abs(prom_peak_frac - iso_peak_prob_frac) < .35:
                        score_boost *= 8
                    elif abs(prom_peak_frac - iso_peak_prob_frac) < .38:
                        score_boost *= 5
                    elif abs(prom_peak_frac - iso_peak_prob_frac) < .46:
                        score_boost *= 2

                    # Increase score based on closeness of i and j
                    # i and j related to prominances and emission probability 
                    if abs(i - j) < 3:
                        score_boost *= 10

                    # Increase score if higher prominances
                    if prom_peak_frac > .65:
                        score_boost *= 5
                    elif prom_peak_frac > .6:
                        score_boost += 30

                    # Additionally boost if middle probability
                    if iso_peak_prob_frac > .4 and iso_peak_prob_frac < 1:
                        score_boost *= 6

                    # Decrease score if low energy due to noise
                    if peak < 90:
                        score_boost *= 0.1
                    elif peak < 100:
                        score_boost *= 0.6 

                    # Increase score based on prominance
                    if i == 1 and prom_peak_frac > 0.5:
                        score_boost *= 10
                    elif i == 1 and prom_peak_frac > 0.3:
                        score_boost *= 2

                    # Decrease score if low prominance
                    if sorted_prom_peaks[i] < 25:
                        score_boost *= 0.001
                    elif sorted_prom_peaks[i] < 40:
                        score_boost *= 0.05

                    # Look through peaks that are similar in multiple isotopes
                    if abs(prom_peak_frac - iso_peak_prob_frac) < 0.1 and any(abs(peak-pk)<2 for pk in sim_peaks) and len(sorted_prom_peaks)>1:
                        if sorted_prom_peaks[0]/sorted_prom_peaks[1] < 5:
                            score_boost *= 0.06

                    # Based on similar peak
                    if abs(peak - 186) < 3 and iso_peak_prob_frac > 0.9 and prom_peak_frac < 0.25:
                         score_boost *= 10

                    # Calculate score
                    scr = score_boost*dist_factor*iso_peak_prob_frac
                    score[index] += scr 

        index+=1

    # Threshold for possible isotopes
    if np.amax(score) > 40:
        ind = np.where(score > 0.4*np.amax(score))
    else:
        ind = []

    # Convert index to isotope names
    guess_name = []
    for i in range(len(name[ind])):
        guess_name.append(name[ind][i][3:])

    # Determine actual isotope
    true_name = determine_name(iso_name)

    # Print guess isotopes
    if len(guess_name) > 0:
        print("The isotope identified: ",guess_name," true isotope: ",true_name," score: ",str(np.round(score[ind])))

    # Score guess vs true isotopes
    scored_iso = score_isotope(guess_name,true_name)

    return scored_iso, peaks, guess_name

# Plot spectra
def plot_spectra(file,counts,peaks,true_peaks=None,height_vec=None,x_lim=None,guess_nm=None):

    # File name
    slash_index = file.rindex('/')+1
    rename_file = file.replace('.n42','')

    x = range(0,len(counts))

    # Create spectra
    plt.step(x,counts)

    # Height threshold
    if height_vec:
        plt.plot(x,height_vec, label = "Threshold")

    # Peaks from isotope database
    if true_peaks:
        plt.vlines(x=true_peaks, ymin = -0.01*max(counts), ymax = max(counts)*1.2, linestyles = 'dashed', colors = 'purple',alpha=0.4, label=rename_file[slash_index:-2].capitalize() + " true peaks")

    # Scatter plot
    plt.scatter([x[int(i)] for i in peaks], [counts[int(i)] for i in peaks],color='deepskyblue',label="Main peaks")
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    plt.xlim(left=0)

    # Create x-axis limit
    if x_lim:
        plt.xlim(right=x_lim)

    # Title based on isotope identification
    name = ""
    for i in range(len(guess_nm)):
        if i == 0:
            name = "Isotope identified: " + guess_nm[0].capitalize()
        else:
            name += ", " + guess_nm[i].capitalize()

    # Define title
    if len(name)>0:
        plt.title(name) 
    else:   
        plt.title("Isotope identified: None") 
    
    plt.legend()
    
    # Extra plot options
    #plt.title(rename_file[slash_index:])
    #plt.rc('font', size=22) 
    #plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
    #plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
    #plt.rc('ytick', labelsize=20)    # fontsize of the tick labels

    # Save figure
    plt.savefig(str(rename_file)+str((np.random.random()))+'.png')
    plt.clf()

# Create class and plot spectra
def isotopeID(file,counts,isotopes,compressed):
    
    # Create spectra
    spect = Spectra(file,counts,compressed)
   
    # Find peaks in spectra
    [peaks,peaks_dict] = find_peaks(spect.counts,prominence=spect.prom,width=spect.wdt_arr,rel_height=0.6,height=spect.height_vec,distance=4)

    # Isotope identification
    score,final_peaks,guess = determine_isotope(peaks[peaks>42],peaks_dict['prominences'][peaks>42],isotopes,spect.iso_name)

    # Only use more prominent peaks in plot
    plt_peaks = peaks[peaks_dict['prominences']>0.45*np.mean(peaks_dict['prominences'])]

    # Create a plot of spectra
    plot_spectra(file,spect.counts,plt_peaks,guess_nm=guess)

    return score

##### ------------------------------------------------------------ #####

parent_path = "/Users/emeline/Documents/NERS 491"

# Define path where xml files can be found
data = "IsotopeID_n42"
#data = "120SecondBkg"
#data = "efficiency data/actual"

# Define path where data is and path where figures are created
data_path = os.path.join(parent_path, data)
figure_path = os.path.join(parent_path, "IsotopeID_spectra")

# Remove figure path if it exists
if os.path.exists(figure_path):
    for root, dirs, files in os.walk(figure_path):
        for file in files:
            os.remove(os.path.join(root, file))  

# Create figure path
else:
    os.mkdir(figure_path)     

isotopes = {}

# Read in isotope database
with open("/Users/emeline/Dev/RadRobo/isotopeid/isotopes_xray.txt") as f:
    for lines in f.readlines():
        line = lines.split(';')

        # X-ray peaks
        xray = None
        if len(line[3]) > 0:
            xray = line[3].replace('\n','')

        # Populate dictionary
        isotopes[line[0].lower()] = {'peaks':line[1],'probability':line[2].replace('\n',''),'xray':xray}

f.close()

# Create list of the names of all xml files in this path
files = [x for x in os.listdir(data_path) if x.endswith('.n42')]

score = 0
length = 0

# Loop through all files in data folder
for file in files:

    # Make sure data in right format otherwise code will fail
    #with open(os.path.join(data_path,file)) as reader:
    #    data = reader.readlines()
    #    if data[1].count("xmlns") > 1:
    #        data[1] = "<RadInstrumentData xmlns='http://physics.nist.gov/N42/2011/N42'>\n"

    # Write back to file
    #with open(os.path.join(data_path,file), 'w') as writer:
    #    writer.writelines(data)

    # Can set a condition here if only want to read in certain files
    if "" in file:
        
        # Read in file
        RadInstrumentData = xmlwrapper.xmlread(os.path.join(data_path,file))
 
        # Read in counts and perform isotope identification
        # Different loops are used due to slighly differentn setups in .n42 files
        if data == "IsotopeID_n42" or "longSpectrum" in file:
            counts = RadInstrumentData.RadMeasurement.Spectrum.ChannelData.text
            compressed = True
            score += isotopeID(os.path.join(figure_path,file),counts,isotopes,compressed)
            length += 1
        elif "actual" in data:
            
            for i in range(1,len(RadInstrumentData.RadMeasurement.Spectrum)):
                counts = RadInstrumentData.RadMeasurement.Spectrum[i].ChannelData.text
                compressed = True
                score += isotopeID(os.path.join(figure_path,file),counts,isotopes,compressed)
                length += 1
        else:
            counts = RadInstrumentData.Measurement.Spectrum.ChannelData.text
            compressed = False
            score += isotopeID(os.path.join(figure_path,file),counts,isotopes,compressed)
            length += 1

# Determine and print final score
score = score/length
print('You scored ' + str(score*100) + ' %!' + " Out of " + str(length) + " files")
