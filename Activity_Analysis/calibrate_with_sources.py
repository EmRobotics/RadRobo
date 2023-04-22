
import json
import datetime
import numpy as np
import Spectrum_analysis as sp
import csv
import glob
import os


"""
Calculation of detector efficiency given folder with .n42 files

Requires: isotope_xray.txt path, folder with .n42 files that contains the source name
and detector angle, .json file with calibration source information

Modifies: writes to efficiency.csv file with energy, value, degree format

Output: calculation of efficiency of detector
"""


#isotope xray file from Isotope ID
isotopes_xray_path = ""
#json file to mantain source information
source_info_path = ""
path = ''


isotopes = {}


def parse_json(json_file_path):
    with open(json_file_path, 'r') as f:
        isotopes = json.load(f)
    return isotopes

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
        
    print(isotopes_xrays)
    return isotopes_xrays



output_file_path = "/Users/undarmaa/Desktop/Documents/RADical Robotics repo/RadRobo/eff.CVS"


def write_efficiencies_to_csv(efficiencies, output_file_path):
    with open(output_file_path, 'a', newline='') as csvfile:
        fieldnames = ['energy', 'efficiency', 'degree']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for efficiency_data in efficiencies:
            writer.writerow(efficiency_data)


def calculate_efficiency_for_spectrum(spectrum, degree):
    live_time = spectrum.n42.live_time
    
        
    efficiencies = []
    
    print(spectrum.peaks)
    
    for peak in spectrum.peaks:
        
        source_id = peak.name 
        source_count = peak.total_clean_peak_count
        peak_energy = peak.peak
        
        
        source_start_date = isotopes[source_id]["date"]
        source_start_activity = isotopes[source_id]["activity_start"]
        source_half_life = isotopes[source_id]["half_life"]
        
        start_datetime = datetime.datetime.strptime(source_start_date, "%B %d %Y %I:%M%p")
        elapsed_time = datetime.datetime.now() - start_datetime
        elapsed_time_in_half_lives = elapsed_time.total_seconds() / (source_half_life * 60 * 60 * 24 * 365.25)
        current_activity = source_start_activity * (0.5 ** elapsed_time_in_half_lives)
        
        
        threshold = 2.0
        
        for i, energy in enumerate(isotope_xray[source_id]["photopeak_energies"]):
            if abs(energy - peak_energy) < threshold:
                expected_energy = energy
                index = i
                break
        else:
            index = None

        branching_ratio = isotope_xray[source_id]["branching_ratios"][index]
        
        #number of count emitted by source, muCi converted to Bq
        N_emitted = current_activity * 3.7 * 10**4 * (branching_ratio /100) * live_time 
        
        #Omega
        solid_angle = 22 * 22 / 1000**2
        
        #number of counts detected by the detector under full energy peak
        #source count was found by finding the sum of the count from -1.5fwhm to 1.5fwhm 
        #from the peak energy and backgground subtracted 
        # by using min max interpolation on both sides of the peak 
        N_detected = source_count 
        
        # intrinic peak eff = abs_eff * 4pi/Omega
        intrinsic_eff = N_detected/N_emitted *  (4 * 3.142) / solid_angle
        
        print(intrinsic_eff)
        efficiencies.append({'energy': expected_energy,
                             'efficiency': intrinsic_eff,
                             #'uncertainty': uncertainty,
                             'degree': degree
                             })
    
    return efficiencies 
        
        
#Path to files for calibration measurements
main_path = ''

n42_files = []
for dirpath, dirnames, filenames in os.walk(main_path):
    if dirpath != main_path:
        continue
    for filename in filenames:
        if filename.endswith('.n42'):
            full_path = os.path.join(dirpath, filename)
            n42_files.append(full_path)



isotope_xray = parse_isotope_xrays(isotopes_xray_path)
isotopes = parse_json(source_info_path)


n42_files = glob.glob(path + '/**/*.n42', recursive=False)

eff = []

for file in n42_files:

    #requires that .n42 file has detector face angle in filename
    if file.endswith(".n42"):
        if "_0" in file:
            degree = 0
        elif "_90" in file:
            degree = 90
        elif "_180" in file:
            degree = 180
        elif "_270" in file:
            degree = 270

    print(file)

    n42_file = sp.N42(file)
    spectrum = sp.Spectrum(n42_file)
    
    spectrum.find_peaks()

    eff = calculate_efficiency_for_spectrum(spectrum, degree)
    print(eff)
    write_efficiencies_to_csv(eff, output_file_path)


