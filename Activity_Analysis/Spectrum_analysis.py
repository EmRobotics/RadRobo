from typing import List, Tuple
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import pdb
import os
import glob
from scipy.signal import find_peaks,savgol_filter

"""
Combined script for spectrum analysis
Find_peaks method adopted and integrated directly from isotope_algorithm.py
This file serves as the main file for maintaining information regarding the .n42
measurement file pre and post analysis. The class dependency is as follows:

N42
Maintains overhead information of detector measurement such as livetime, filename, raw counts etc

Spectrum
Maintains spectrum information such as identified peaks, guess of isotope, spectrum by channel,
plot method, background estimation etc
Has Peak and N42 class objects

Peak
Maintains individual peak information such as guess of the source (if multiple sources), energy,
background count, total count, estimated count etc

This script reads in a .n42 file and performs spectrum analysis on the data. It finds peaks in the
spectrum and identifies the most likely isotope responsible for the observed peaks. The script
also provides various methods for plotting and analyzing the spectrum data, such as calculating
the full width at half maximum (FWHM) for each peak, estimating background counts, and
subtracting background to obtain clean peak counts. The script relies on an external file
(isotopes_xray.txt) for information about the isotopes and their corresponding x-ray peaks.
"""

#File path to .n42 file
file_path = ""

# INCLUDE PATH TO isotopes_xray.txt HERE
isotopes_xray_path = ""

isotopes = {}

with open(isotopes_xray_path) as f:
    for lines in f.readlines():
        line = lines.split(';')

        xray = None
        if len(line[3]) > 0:
            xray = line[3].replace('\n','')

        isotopes[line[0].lower()] = {'peaks':line[1],'probability':line[2].replace('\n',''),'xray':xray}

f.close()

class N42:
    def __init__(self, filepath: str):
        self.filepath = filepath
        with open(filepath, "r") as openFile:
            self.spec_text = openFile.read()
            
        self.spec_id = ['RadMeasurement-1-Spectrum-1-PUR', 'RadMeasurement-1-Spectrum-1']
        self.namespace = {'' : 'http://physics.nist.gov/N42/2011/N42'}
        self.parse_header()
        

    def parse_header(self):
        # Parse header information and store as instance variables  
        xml_text = ET.fromstring(self.spec_text)
        for child in xml_text.findall('RadMeasurement', self.namespace):
            #Q: Is "Foreground" text always in tag "MeasurementClassCode"? If yes, does it occur only once?          
            if (child.find('MeasurementClassCode', self.namespace) == None or child.find('MeasurementClassCode', self.namespace).text != "Foreground"):
                continue     
            self.real_time = float(child.find('RealTimeDuration', self.namespace).text.split('T')[1].split('S')[0])
            
            for spec in child.findall('Spectrum', self.namespace):
                if spec.attrib['id'] in self.spec_id:
                    self.live_time = float(spec.find('LiveTimeDuration', self.namespace).text.split('T')[1].split('S')[0])
                    self.parse_data(spec.find('ChannelData', self.namespace).text)
        

    def parse_data(self, text, bins_per_kev = 1):
        # Parse data section and store as instance variables
        try:
            text_spectrum = text.split()
            
            spectrum = []
            propagate_zeros = False
            for i in text_spectrum:
                if float(i) == 0:
                    propagate_zeros = True
                    continue
                elif propagate_zeros:
                    propagate_zeros = False
                    for o in range(int(i)):
                        spectrum.append(0)
                else:
                    spectrum.append(float(i))
                    
        except AttributeError:
            spectrum = 3000*[0]
    
        # We only really care about 1 keV bins for this code so make it so
        kev_per_bin_spectrum = []
        if bins_per_kev == 1:
            first_time = False
        else:
            first_time = True
        bin_sum = 0
        for indx, bin in enumerate(spectrum):
            if indx % bins_per_kev == 0:
                if not first_time:
                    bin_sum += bin
                    kev_per_bin_spectrum.append(bin_sum)
                    bin_sum = 0
                else:
                    first_time = False
            else:
                bin_sum += bin
                    
        kev_per_bin_spectrum.append(bin_sum)
            
        self.count_spectrum = np.array(kev_per_bin_spectrum)

def parse_isotope_xrays(xrays_file_path):
    with open(xrays_file_path, 'r') as f:
        isotopes_xrays = {}
        for line in f.readlines():
            tokens = line.split(';')
            
            isotope_id = tokens[0].split("-")
            source_id = f"{isotope_id[1]}-{isotope_id[2]}"
            
            photopeak_energies = list(map(float, tokens[1].split(',')))
            branching_ratios = list(map(float, tokens[2].split(',')))
            
            isotopes_xrays[source_id] = {
                'photopeak_energies': photopeak_energies,
                'branching_ratios': branching_ratios,
                'energy_and_ratio': list(zip(photopeak_energies, branching_ratios))
            }
        
    #print(isotopes_xrays)
    return isotopes_xrays


isotope_xray = parse_isotope_xrays(isotopes_xray_path)


class Spectrum:
    def __init__(self, n42: N42):
        self.n42 = n42
        self.counts = n42.count_spectrum
        
        
        self.wdt_arr = np.array([max(i,2) for i in np.arange(0,len(self.counts),1)*0.005])
        self.height_vec = find_height(self.counts)
        self.prom = 0.01*max(self.counts)
        
        
        self.peaks = []
        self.guess = []

    def find_peaks(self):
        # Find peaks in the spectrum and store as instances of Peak
        
        [peaks, peaks_dict] = find_peaks(self.counts,prominence=self.prom,width=self.wdt_arr,rel_height=0.6,height=self.height_vec,distance=4)
        peaks, guess, real_peaks = determine_isotope(peaks,peaks_dict['prominences'],isotopes)
        
        #Uncomment if manually setting isotope guess
        #guess = "Cs-137"
        
        #If using expected peak energies instead of detected
        real_peaks = (isotope_xray[guess]['photopeak_energies'])
        peaks = (isotope_xray[guess])
        
        print(f"All detected peak(s) {peaks}, guess Isotope is {guess}, using line(s) {real_peaks}")
        #init peaks
        for i, peak in enumerate(real_peaks):
            print(peak)
            self.peaks.append(Peak(int(peak), guess, self.counts))
        return peaks, guess
        

    def calculate_fwhm(self):
        fwhms = []
        for peak in self.peaks:
            fwhms.append(peak.findFWHM())
        return fwhms
            

    def net_count(self):
        # Calculate net counts for each peak
        peak_counts = []
        for peak in self.peaks:
            peak_counts.append(peak.total_peak_count)
        return peak_counts

    def background_count(self):
        bgs = []
        for peak in self.peaks:
            bgs.append(peak.total_bg)
        return bgs

    def clean_peak_count(self):
        clean_peak_count = []
        for peak in self.peaks:
            clean_peak_count.append(peak.clean_peak_count)
        return clean_peak_count
    
    def plot_spectrum(self):
        plt.figure()
        plt.plot(self.counts, marker = ".")
        plt.title(self.n42.filepath)
        plt.xlabel('Energy (keV)')
        plt.ylabel('Count')
        
        for i in self.peaks:
            plt.plot(i.peak, self.counts[i.peak], marker="x", color = 'r', label = 'peak')
            
            i.print_peak_info()
            
            #Background
            plt.fill_between(range(i.bg_left_idx, i.total_peak_left_idx + 1), (self.counts[i.bg_left_idx:i.total_peak_left_idx + 1]), color = 'c', alpha = 0.6)
            plt.fill_between(range(i.total_peak_right_idx, i.bg_right_idx + 1), (self.counts[i.total_peak_right_idx:i.bg_right_idx + 1]), color = 'c', alpha = 0.6, label = 'background')
           
            #total peak 
            plt.fill_between(range(i.total_peak_left_idx, i.total_peak_right_idx+ 1), (self.counts[i.total_peak_left_idx:i.total_peak_right_idx+ 1]), color = 'b', label = 'net count')
            
            #clean peak
            plt.plot(i.clean_peak_x , i.clean_peak_y, color = 'r', linestyle = '--', linewidth = 0.6, label = 'background subtracted counts') 
        
        plt.legend()
        plt.show()
                      
        


class Peak(Spectrum):
    def __init__(self, energy, name, counts,  net_counts=None, fwhm=None, background=None):
        self.counts = counts
        self.peak = energy
        self.name = name
        
        self.findFWHM()
        self.estimate_peak_counts()
        
        
        
    def findFWHM(self):
               
            half_max = self.counts[(self.peak)]/2
            left_idx = self.peak
            right_idx = self.peak
            while left_idx > 0 and self.counts[left_idx] > half_max:
                left_idx -= 1
            while right_idx < len(self.counts) - 1 and self.counts[right_idx] > half_max:
                right_idx += 1
            self.fwhm_left_idx = left_idx
            self.fwhm_right_idx = right_idx
            self.fwhm = right_idx - left_idx

    def estimate_peak_counts(self):
        # Estimate background for the peak
        self.total_peak_left_idx  =  max(round(self.peak - self.fwhm * 1.5), 0)
        self.total_peak_right_idx =  min(round(self.peak + self.fwhm * 1.5), len(self.counts))
        
         
        self.total_peak_count = 0 if (self.total_peak_left_idx < 0 or self.total_peak_right_idx > len(self.counts)) else sum(self.counts[self.total_peak_left_idx:self.total_peak_right_idx+1])
       
        #background 
        self.bg_left_idx = max(round(self.peak - self.fwhm * 2.5), 0) #self.total_peak_left_idx)
        self.bg_right_idx = min(round(self.peak + self.fwhm * 2.5), len(self.counts)) #self.total_peak_right_idx))
        
        
        self.left_bg_count = sum(self.counts[self.bg_left_idx:self.total_peak_left_idx+1])
        self.right_bg_count = 0 if (self.bg_left_idx < 0 or self.bg_right_idx > len(self.counts)) else sum(self.counts[self.total_peak_right_idx:self.bg_right_idx+1])
    
      
        x_min = self.total_peak_left_idx
        x_max = self.total_peak_right_idx
        y_min = self.counts[x_min]
        y_max = self.counts[x_max]
        x_val = np.arange(x_min, x_max + 1)
        x_mixing = (x_val - x_min) / (x_max - x_min)
        y_val = (1 - x_mixing) * y_min + x_mixing * y_max
        
        self.total_bg_count = sum(y_val)
        self.total_bg = y_val
        self.total_clean_peak_count = sum((self.counts[x_min:x_max + 1] - y_val))
        self.clean_peak_y = self.counts[x_min:x_max + 1] - y_val
        self.clean_peak_x = x_val
        
    def print_peak_info(self):
        print(f"Peak is at energy {self.peak}, potential isotope {self.name}")
        print(f"Peak net count is {self.total_peak_count}, from idx {self.total_peak_left_idx} to {self.total_peak_right_idx}")
        print(f"The FWHM for this peak is {self.fwhm}")
        print(f"BG count left to peak is {self.left_bg_count}, from idx {self.bg_left_idx} to {self.total_peak_left_idx}")
        print(f"BG count right to peak {self.right_bg_count}, from idx {self.total_peak_right_idx} to {self.bg_right_idx}")
        print(f"BG count estimate is {self.total_bg_count}")
        print(f"Total clean peak count is estimate {self.total_clean_peak_count}")
        
        


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


#Guess isotope based on peaks found
def determine_isotope(all_peaks,all_prominences,isotopes):

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


        name = np.append(name,key)
        

        iso_peaks = np.array([float(i) for i in val['peaks'].split(',')])
        prob_peaks = np.array([float(i) for i in val['probability'].split(',')])

        sorted_index_array = np.argsort(prob_peaks)
        sorted_iso_peaks = np.flip(iso_peaks[sorted_index_array],0)
        sorted_prob_peaks = np.flip(prob_peaks[sorted_index_array],0)

        sorted_ind_array = np.argsort(prominences)
        sorted_peaks = np.flip(peaks[sorted_ind_array],0)
        sorted_prom_peaks = np.flip(prominences[sorted_ind_array],0)

       
        for i in range(len(sorted_peaks)):
            for j in range(len(sorted_iso_peaks)):

                peak = sorted_peaks[i]
                iso_peak = sorted_iso_peaks[j]

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

                    if sorted_prom_peaks[i] < 25:
                        score_boost *= 0.001
                    elif sorted_prom_peaks[i] < 40:
                        score_boost *= 0.05

                    if abs(prom_peak_frac - iso_peak_prob_frac) < 0.1 and any(abs(peak-pk)<2 for pk in sim_peaks) and len(sorted_prom_peaks)>1:
                        if sorted_prom_peaks[0]/sorted_prom_peaks[1] < 5:
                            score_boost *= 0.06 #0.1 was working
                        
                    if abs(peak - 186) < 3 and iso_peak_prob_frac > 0.9 and prom_peak_frac < 0.25:
                         score_boost *= 10

                    scr = score_boost*dist_factor*iso_peak_prob_frac
                    score[index] += scr 
                     

        index+=1


    if np.amax(score) > 40:
        ind = np.where(score > 0.4*np.amax(score))
    else:
        ind = []
        guess_name = None

    real_lines = []

    for i in range(len(name[ind])):
        if i == 0:
            guess_name = str(name[ind][i][3:]).capitalize()
          
            real_lines = isotopes[name[ind][i]]['peaks']
            
        else:
            guess_name += ", " + str(name[ind][i][3:]).capitalize()
            real_lines.append(isotopes[name[ind][i]]['peaks'])

    real_lines = real_lines.split(',')
    
    real_peaks = []
    
    for peak in all_peaks:
        if any(abs(float(line) - peak) < 2 for line in real_lines): 
           real_peaks.append(peak)
           
    print(f"The isotope(s) identified: {guess_name} with expected peak(s) {real_peaks}")

    return peaks,guess_name, real_peaks




    


    
    
    
   