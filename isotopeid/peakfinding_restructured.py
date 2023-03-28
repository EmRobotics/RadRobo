import os
import xmlwrapper
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks,savgol_filter
      
class Spectra:
    def __init__(self, file, counts,compressed):
        self.iso_name = file[file.rindex('/')+1:file.rindex('_')]

        if compressed:
            self.counts = expand_zeros(counts)
        else:
            #self.counts = counts = [float(el) for el in counts.split(" ")[:-1]]
            self.counts = [float(el) for el in counts.split(" ")[:-1]]
            #self.counts =  savgol_filter([float(el) for el in counts.split(" ")[:-1]],4,2)
        self.wdt_arr = np.array([max(i,2) for i in np.arange(0,len(self.counts),1)*0.005])
        self.height_vec = find_height(self.counts)
        self.prom = 0.01*max(self.counts)

# Expand the zeros in compressed spectra
def expand_zeros(counts):
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

    return counts   

# Determine a height threshold
def find_height(counts):
    der = []
    std_dev = []
    for i in range(1,len(counts)-1):
        der.append(counts[i-1] - 2*counts[i] + counts[i+1])
        std_dev.append(np.sqrt(counts[i-1] + 4*counts[i] + counts[i+1]))
    std_dev.append(0)
    std_dev.insert(0,0)  

    counts_filtered = savgol_filter(counts, 50, 2)
    height = np.array([int(max(std_dev[i]+0.98*counts_filtered[i],1)) for i in range(len(std_dev))])

    return height 

# Determine isotopes given the name
def determine_name(iso_name):

    iso_name = iso_name.lower()

    if "heu" in iso_name:
        iso_name = iso_name.replace("heu","u235u-238")
    elif "du" in iso_name:
        #iso_name = iso_name.replace("du","u238")
        iso_name = iso_name.replace("du","u235u-238")
    elif "wgpu" in iso_name:
        iso_name = iso_name.replace("wgpu","pu239")

    old_nm = iso_name
    for x in range(len(old_nm)):
        if iso_name[x].isalpha() and iso_name[x+1].isdigit():
            iso_name = iso_name[:x+1] + '-' + iso_name[x+1:]

        elif iso_name[x].isdigit() and iso_name[x+1].isalpha():
            iso_name = iso_name[:x+1] + ',' + iso_name[x+1:]

    if "tc-99" in iso_name:
        iso_name += "m"

    return iso_name.split(',')

# Find true isotope peaks
def isotope_peaks(iso_name,isotope_dict):

    iso_list = determine_name(iso_name)
    peaks = ""
    prob = ""

    for x in range(len(iso_list)):
        for key,val in isotope_dict.items():
            if key[3:] == iso_list[x]:

                if x > 0:
                    peaks += "," + val['peaks']
                    prob += "," + val['probability']
                else:
                    peaks +=  val['peaks']
                    prob +=  val['probability']
   
    return [float(i) for i in peaks.split(',')],[float(i) for i in prob.split(',')]

# Score isotope id
def score_isotope(guess_name,true_name):
    score = 0

    #if len(guess_name) > 0:
    #    return 0
    #else:
    #    return 1

    
    #use this for normal spectra
    if not true_name[0]=="":
        if len(guess_name) > len(true_name):
            return 1/len(guess_name)

        for i in range(len(true_name)):
            for j in range(len(guess_name)):
                if true_name[i] == guess_name[j]:
                    score += 1        

        return score/len(true_name) #fraction of 1
    else:
        if len(guess_name) > 0:
            #need to check this
            if len(guess_name) == 1 and guess_name[0] == 'ra-226':
                return 1
            else:
                return 0
        else:
            return 1

def xray_list(isotopes):

    xrays = np.array([])
    for key,val in isotopes.items():
        if val['xray']:
            for i in val['xray'].split(','):
                xrays = np.append(xrays,float(i))

    return xrays

# Guess isotope based on peaks found
def determine_isotope(all_peaks,all_prominences,isotopes,iso_name):

    score = np.zeros(len(isotopes))
    name = np.array([])
    index = 0
    sim_peaks = [93, 186]

    for key,val in isotopes.items():

        score_boost_og = 1
        if val['xray']:
            xrays = np.array([])

            for i in val['xray'].split(','):
                xrays = np.append(xrays,float(i))

            peaks = np.array([])
            prominences = np.array([])

            for x in range(len(all_peaks)):
                if not any(abs(xray - all_peaks[x]) < 2 for xray in xrays): #or abs(all_peaks[x] -94)<0.05:
                    peaks = np.append(peaks,all_peaks[x])
                    prominences = np.append(prominences,all_prominences[x])
                else:
                    score_boost_og *= 2 #2?
        else:
            peaks = all_peaks
            prominences = all_prominences 


        #if "u" in key:
        #print(key,score_boost_og)

        name = np.append(name,key)

        iso_peaks = np.array([float(i) for i in val['peaks'].split(',')])
        prob_peaks = np.array([float(i) for i in val['probability'].split(',')])

        sorted_index_array = np.argsort(prob_peaks)
        sorted_iso_peaks = np.flip(iso_peaks[sorted_index_array],0)
        sorted_prob_peaks = np.flip(prob_peaks[sorted_index_array],0)


        #print(key, sorted_iso_peaks, sorted_prob_peaks)

        sorted_ind_array = np.argsort(prominences)
        sorted_peaks = np.flip(peaks[sorted_ind_array],0)
        sorted_prom_peaks = np.flip(prominences[sorted_ind_array],0)

        #print(sorted_peaks, sorted_prom_peaks)
        #if any(sorted_prom_peaks > 40):
        #print(sorted_prom_peaks)

        for i in range(len(sorted_peaks)):
            for j in range(len(sorted_iso_peaks)):

                peak = sorted_peaks[i]
                iso_peak = sorted_iso_peaks[j]

                #if abs(peak - iso_peak) < 2.5:
                if abs(peak - iso_peak) < 2.5:

                    dist_factor = 1-(abs(peak - iso_peak)/2)*0.1

                    iso_peak_prob = sorted_prob_peaks[j]
                    iso_peak_prob_frac = iso_peak_prob/max(prob_peaks)

                    prom_peak_frac = sorted_prom_peaks[i]/max(sorted_prom_peaks)

                    score_boost = min(5,score_boost_og) #1 #score_boost_og
                    if abs(prom_peak_frac - iso_peak_prob_frac) < .25:
                        score_boost *= 10
                    elif abs(prom_peak_frac - iso_peak_prob_frac) < .35:
                        score_boost *= 8   #5 is working
                    elif abs(prom_peak_frac - iso_peak_prob_frac) < .38:
                        score_boost *= 5 
                    elif abs(prom_peak_frac - iso_peak_prob_frac) < .46:
                        score_boost *= 2      

                    if abs(i - j) < 3:
                        score_boost *= 10

                    if prom_peak_frac > .65: #.6 is working
                        score_boost *= 5
                    elif prom_peak_frac > .6 :#and score_boost < 20: 
                        score_boost += 30

                    if iso_peak_prob_frac > .4 and iso_peak_prob_frac < 1:
                        score_boost *= 6

                    if peak < 90:
                        score_boost *= 0.1
                    elif peak < 100:
                        score_boost *= 0.6 

                    if i == 1 and prom_peak_frac > 0.5:
                        score_boost *= 10
                    elif i == 1 and prom_peak_frac > 0.3:
                        score_boost *= 2

                    # FOR BACKGROUND
                    #if sorted_prom_peaks[i] < 40:
                    # 40 99% accuracy for background, 94% accuracy for isotopes
                    # 25 100% accuracy for isotopes
                    # 27 99% accuracy for isotopes

                    # 100% accuracy for isotopes, 99.34% accuracy for background
                    if sorted_prom_peaks[i] < 25:
                        score_boost *= 0.001
                    elif sorted_prom_peaks[i] < 40:
                        score_boost *= 0.05

                    #if sorted_prom_peaks[i] < 45:
                    #    score_boost *= 0.0001
                    #maybe do something where if no peaks are greater than 20 than just... idk


                    #if i == 0 and abs(peak - 93) < 2:
                    #    print(sorted_prom_peaks[0]/sorted_prom_peaks[1] )

                    if abs(prom_peak_frac - iso_peak_prob_frac) < 0.1 and any(abs(peak-pk)<2 for pk in sim_peaks) and len(sorted_prom_peaks)>1:
                        if sorted_prom_peaks[0]/sorted_prom_peaks[1] < 5:
                            score_boost *= 0.06 #0.1 was working
                        #print(i,iso_name,"isotope: ", key, sorted_prom_peaks[0]/sorted_prom_peaks[1])

                    if abs(peak - 186) < 3 and iso_peak_prob_frac > 0.9 and prom_peak_frac < 0.25:
                         score_boost *= 10



                    #scr = score_boost*sorted_prom_peaks[i]/max(sorted_prom_peaks)*dist_factor*iso_peak_prob_frac #*sorted_prom_peaks[i]
                    #score_boost*dist_factor*iso_peak_prob_frac*100/(i+1)/(j+1)

                    scr = score_boost*dist_factor*iso_peak_prob_frac
                    score[index] += scr 
                     
                    
                    #print("peak:",peak,"iso-peak:",iso_peak,"iso:",key)
                    #print("i:",i,"j:",j,"dist:",dist_factor,"prob frac:",iso_peak_prob_frac)
                    #print("prob peak:",iso_peak_prob,"prom peak:",sorted_prom_peaks[i],"prom frac:",sorted_prom_peaks[i]/max(sorted_prom_peaks))
                    #print("score:",scr,"\n")

        index+=1

    # u-238
    #for peak in peaks:
    #    if abs(peak - 1001) < 3:
    #        score[np.where(name == "92-u-238")] += 1000

    #ind = np.where(score > 0.25*np.amax(score))
    if np.amax(score) > 40:
        ind = np.where(score > 0.4*np.amax(score))
    else:
        ind = []

    # newest is 0.4 is working
    #0.6 is working

    guess_name = []
    for i in range(len(name[ind])):
        guess_name.append(name[ind][i][3:])

    true_name = determine_name(iso_name)

    print("The isotope identified: ",guess_name," true isotope: ",true_name," score: ",str(np.round(score[ind])))

    scored_iso = score_isotope(guess_name,true_name)

    return scored_iso, peaks,guess_name

# Plot spectra
def plot_spectra(file,counts,peaks,true_peaks=None,height_vec=None,x_lim=None,guess_nm=None):

    slash_index = file.rindex('/')+1
    rename_file = file.replace('.n42','')

    x = range(0,len(counts))

    #create spectra
    plt.step(x,counts)

    if height_vec:
        plt.plot(x,height_vec, label = "Threshold")

    if true_peaks:
        plt.vlines(x=true_peaks, ymin = -0.01*max(counts), ymax = max(counts)*1.2, linestyles = 'dashed', colors = 'purple',alpha=0.4, label=rename_file[slash_index:-2].capitalize() + " true peaks")

    plt.scatter([x[int(i)] for i in peaks], [counts[int(i)] for i in peaks],color='deepskyblue',label="Main peaks")
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    plt.xlim(left=0)

    if x_lim:
        plt.xlim(right=x_lim)

    name = ""
    for i in range(len(guess_nm)):
        if i == 0:
            name = "Isotope identified: " + guess_nm[0].capitalize()
        else:
            name += ", " + guess_nm[i].capitalize()

    if len(name)>0:
        plt.title(name) # + str(peaks))
    else:   
        plt.title("Isotope identified: None") 
        #plt.title(str(peaks))
    #plt.title(rename_file[slash_index:].capitalize() + " Energy Spectra")
    plt.legend()
    #plt.title(rename_file[slash_index:])

    plt.savefig(str(rename_file)+str((np.random.random()))+'.png')
    plt.clf()

# Create class and plot spectra
def isotopeID(file,counts,isotopes,compressed):
    
    spect = Spectra(file,counts,compressed)
   
    [peaks,peaks_dict] = find_peaks(spect.counts,prominence=spect.prom,width=spect.wdt_arr,rel_height=0.6,height=spect.height_vec,distance=4)

    #print(peaks[peaks>42],peaks_dict['prominences'][peaks>42])
    #true_peaks,branching_ratios = isotope_peaks(spect.iso_name,isotopes)

    #score,final_peaks,guess = determine_isotope(peaks,peaks_dict['prominences'],isotopes,spect.iso_name)
    score,final_peaks,guess = determine_isotope(peaks[peaks>42],peaks_dict['prominences'][peaks>42],isotopes,spect.iso_name)

    #print(peaks_dict['prominences'])
    #plot_spectra(file,spect.counts,final_peaks,true_peaks=true_peaks,guess_nm=guess)
    #plot_spectra(file,spect.counts,final_peaks,guess_nm=guess) #true_peaks=true_peaks,guess_nm=guess) #height=counts.height_vec,x_lim=
    # 
    plt_peaks = peaks[peaks_dict['prominences']>0.5*np.mean(peaks_dict['prominences'])]

    plot_spectra(file,spect.counts,plt_peaks,guess_nm=guess)

    return score

##### ----------------------------- CODE----------------------------- #####

parent_path = "/Users/emeline/Documents/NERS 491"

# Define path where xml files can be found
#path = "/Users/emeline/Documents/NERS 491/IsotopeID_n42/"
data = "IsotopeID_n42"

data = "120SecondBkg"

#data = "efficiency data/actual"

data_path = os.path.join(parent_path, data)
figure_path = os.path.join(parent_path, "IsotopeID_spectra")

if os.path.exists(figure_path):
    for root, dirs, files in os.walk(figure_path):
        for file in files:
            os.remove(os.path.join(root, file))  
else:
    os.mkdir(figure_path)     

isotopes = {}

#with open("/Users/emeline/Dev/RadRobo/isotopeid/isotopes_limited.txt") as f:
with open("/Users/emeline/Dev/RadRobo/isotopeid/isotopes_xray.txt") as f:
    for lines in f.readlines():
        line = lines.split(';')

        xray = None
        if len(line[3]) > 0:
            xray = line[3].replace('\n','')

        isotopes[line[0].lower()] = {'peaks':line[1],'probability':line[2].replace('\n',''),'xray':xray}

f.close()

# Create list of the names of all xml files in this path
files = [x for x in os.listdir(data_path) if x.endswith('.n42')]

score = 0
length = 0

for file in files:
    #print(file)
    if "" in file:
        RadInstrumentData = xmlwrapper.xmlread(os.path.join(data_path,file))
 
        if data == "IsotopeID_n42" or "longSpectrum" in file:
            counts = RadInstrumentData.RadMeasurement.Spectrum.ChannelData.text
            compressed = True
            score += isotopeID(os.path.join(figure_path,file),counts,isotopes,compressed)
            length += 1
        elif "actual" in data: #efficiency
            
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



        #try: 
        #    if (RadInstrumentData.RadMeasurement.Spectrum.ChannelData.Attributes.compressionCode):
        #        print(RadInstrumentData.RadMeasurement.Spectrum.ChannelData.Attributes)
        #        compressed = True
        #        counts = RadInstrumentData.RadMeasurement.Spectrum.ChannelData.text
        #except KeyError:
        #    counts = RadInstrumentData.Measurement.Spectrum.ChannelData.text

        #score += isotopeID(os.path.join(figure_path,file),counts,isotopes,compressed)
        #length += 1

score = score/length

print('You scored ' + str(score*100) + ' %!' + " Out of " + str(length) + " files")
