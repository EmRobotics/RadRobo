
import json
import datetime
import numpy as np
import Spectrum_analysis as sp
import csv
import glob
import os
import pandas as pd



"""
Driver for activity analysis

Requires: estimated distance, estimated direction of source, isotope_xray.txt, .n42 file path,
efficiency.csv path

Dependencies: Spectrum_analysis.py

Output: estimated current activity of source
"""

# read in csv file and create dataframe
df = pd.read_csv('efficiency.csv', header=None, names=['Energy', 'Value', 'Degree'])

# create lookup table
lookup_table = {}
for degree in [0, 90, 180, 270]:
    degree_df = df[df['Degree'] == degree]
    lookup_table[degree] = dict(zip(degree_df['Energy'], degree_df['Value']))



# Function to interpolate efficiency
def interpolated_efficiency(angle, energy):
    # Find the two closest angles in the lookup table
    closest_angles = sorted(lookup_table.keys(), key=lambda x: abs(x - angle))[:2]

    # Calculate the weights for each angle based on their distance to the input angle
    angle_weights = [1 - abs(closest_angle - angle) / 90 for closest_angle in closest_angles]

    interpolated_efficiencies = []

    for closest_angle, angle_weight in zip(closest_angles, angle_weights):
        energy_efficiency_dict = lookup_table[closest_angle]

        # Find the two closest energies in the lookup table for the current angle
        closest_energies = sorted(energy_efficiency_dict.keys(), key=lambda x: abs(x - energy))[:2]

        # Calculate the weights for each energy based on their distance to the input energy
        energy_weights = [1 - abs(closest_energy - energy) / abs(closest_energies[1] - closest_energies[0])
                          for closest_energy in closest_energies]

        # Interpolate the efficiency for the current angle using the closest energy values
        efficiency = np.dot([energy_efficiency_dict[e] for e in closest_energies], energy_weights)

        # Weight the interpolated efficiency by the angle weight
        weighted_efficiency = efficiency * angle_weight
        interpolated_efficiencies.append(weighted_efficiency)

    # Combine the weighted efficiencies to get the final interpolated efficiency
    return sum(interpolated_efficiencies)



if __name__ == "__main__":
    # Instantiate N42File object and pass to Spectrum object
    file_path = ""
    n42_file = sp.N42(file_path)
    spectrum = sp.Spectrum(n42_file)
    
    # Analyze spectrum
    spectrum.find_peaks()

    
    distance = 0  #REPLACE WITH DISTANCE FROM SOURCE TO DETECTOR OR CALL FUNCTION (cm)
    angle =  0 #REPLACE WITH SOURCE DIRECTION OR CALL FUNCTION 
    
    
    
    for peak in spectrum.peaks:
        energy = peak.energy
        counts = peak.counts
        name = peak.name
        
            
        eff = interpolated_efficiency(angle, energy)
        branching_ratio = 0 #obtain from isotope_xray.txt
        
        solid_angle = 22 * 22 / (distance * 10) **2
        
        N_detected = peak.total_clean_peak_count
        N_detected_uncertainty = np.sqrt(N_detected)  # Assuming Poisson statistics for the counts

    
        N_emitted = (N_detected / eff * (4* 3.142/solid_angle)) 
        N_emitted_uncertainty = N_emitted * np.sqrt((N_detected_uncertainty / N_detected) ** 2 + (eff_uncertainty / eff) ** 2)
        
        
        print(f"N_detected {N_detected} and N emitted: {N_emitted}")
        
        current_activity = N_emitted / (3.7 * 10**4 * (branching_ratio /100) * spectrum.n42.live_time )
        
        current_activity_uncertainty = current_activity * np.sqrt((N_emitted_uncertainty / N_emitted) ** 2 )
       
        print(f"Estimated activity of source is {current_activity}")
    


