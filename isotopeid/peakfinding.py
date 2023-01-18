import os
import xmlwrapper
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks,savgol_filter
      

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

def isotope_peaks(iso_name,isotope_dict):

    iso_name = iso_name.lower()

    if "heu" in iso_name:
        iso_name = iso_name.replace("heu","u235")
    elif "du" in iso_name:
        iso_name = iso_name.replace("du","u238")
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

    iso_list = iso_name.split(',')
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

    #31-Ga-67;93.311,184.577,208.951,300.219,393.529,494.169,703.110,794.386,887.693;39.2,21.2,2.40,16.80,4.68,0.0691,0.0106,0.0540,0.149

def score_peaks(peaks,true_peaks):
    print('good job!')

def plot_spectra(file,x,counts,peaks,true_peaks):

    slash_index = file.rindex('/')+1
    rename_file = file.replace('.n42','')

    #create spectra
    plt.step(x,counts)
    #plt.step(x,counts_filtered)
    #print(peaks)

    #plt.hlines(y=0.1*max(counts), xmin = 0, xmax=max(x), linestyles = 'dashed', colors = 'green')
    plt.vlines(x=true_peaks, ymin = -0.2*max(counts), ymax = max(counts)*1.2, linestyles = 'dashed', colors = 'purple',alpha=0.2)

    for i in peaks:
        plt.scatter(x[i], counts[i])

    #plt.xlabel('Channel number')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    #plt.xlim(right=250)

    #if "am241" in rename_file[slash_index:]:
    #    plt.yscale('log')

    plt.title(rename_file[slash_index:])
    plt.savefig(str(rename_file)+'.png')
    plt.clf()


def createspectra(file,counts,isotopes):
    counts = expand_zeros(counts)

    x = range(0,len(counts))

    rename_file = file.replace('.n42','')
    tru_peaks,branching_ratios = isotope_peaks(rename_file[file.rindex('/')+1:rename_file.rindex('_')],isotopes)

    #print(branching_ratios)

    true_peaks = []

    for i in range(len(tru_peaks)):
        if branching_ratios[i] > 0: #max(branching_ratios)/100:
            true_peaks.append(tru_peaks[i])

    #print(true_peaks)
    #wdt_arr = np.arange(0,len(counts),1)*0.01

    counts_filtered  = savgol_filter(counts, 40, 5)
    #counts_sm = filtfilt(counts)

    #thresh = np.zeros(len(counts))

    #for i in range(0,len(counts)):
        #thresh[i] = 2*counts_sm[i]
        #sum(counts[i-round(wdt_arr[i]):i+round(wdt_arr[i])])/wdt_arr[i]
        #print(thresh[i])

    #[peaks,peak_dict] = find_peaks(counts, width = wdt_arr, distance = 4,rel_height=0.75, height = 2*counts_sm) #, prominence = 4, threshold = 40)
    #[peaks,peak_dict] = find_peaks(counts, width=2, height=0.1*max(counts))
    [peaks,peak_dict] = find_peaks(counts,width=2,prominence=100,threshold=20)
    #[peaks,peak_dict] = find_peaks(counts,width=2,prominence=30,threshold=5)

    #determine_isotope(points, isotopes,rename_file[slash_index:])

    # these can all be arrays!! (except distance)
    #height -> Required height of peaks
    #distance -> Required minimal horizontal distance (>= 1) in samples between neighbouring peaks
    #prominence -> Required prominence of peaks.
    #width -> Required width of peaks in samples.
    #threshold -> Required threshold of peaks, the vertical distance to its neighbouring samples

    #poisson limit is 2.35/sqrt(N)*energy
    #(2.35/(.01/667))^2

    plot_spectra(file,x,counts,peaks,true_peaks)



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

isotopes = {}

with open("/Users/emeline/Dev/RadRobo/isotopeid/isotopes.txt") as f:
    for lines in f.readlines():
        line = lines.split(';')

        isotopes[line[0].lower()] = {'peaks':line[1],'probability':line[2].replace('\n','')}

# Create list of the names of all xml files in this path
files = [x for x in os.listdir(data_path) if x.endswith('.n42')]



for file in files:

    #if "du" in file:

    RadInstrumentData = xmlwrapper.xmlread(os.path.join(data_path,file))
    counts = RadInstrumentData.RadMeasurement.Spectrum.ChannelData.text
    createspectra(os.path.join(figure_path,file),counts,isotopes)


    #second derivative?
    # in this gross count method, the spectrum is split into discrete energy ‘bins’ and the counts in each 
    # are then summed. If the summed counts in a specific bin are found to be over a predetermined threshold 
    # value relative to the neighbouring bins, then the existence of a peak is consequently flagged.


    #—--------------------------------------------
    # WHAT ARE THE PEAKS AT 77 AND 88 FOR RA-226?????
    # https://pdf.sciencedirectassets.com/272707/1-s2.0-S0020708X00X03316/1-s2.0-0020708X83902752/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEOj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQCplm9GHN0n5gjJAWPXO8hHxHh%2FosanfRRZDLE5mFgxYgIgYQGHOGNWMzR5rI6DXPDp%2BYLbD%2Ffc3aDRKbK12Qfq7CAqzAQIURAFGgwwNTkwMDM1NDY4NjUiDAs7HRpuAfydkWLaKCqpBAZTaGvnPNS3WwFPCt2dFfIDL7MpDadZpwK%2B8k7gw1ItbvC8MeeQHGJddaZeGlk4C6c62x525WX1hRIE2cbmsmNLA6MGpemKytynWp48ycRsXAdWimwb86WkZ0d2heojSrSJFLh7WU16SIO3pHXV8WUbkMLUeZ0hCJIQq5lQ38joR1MrSJJOyMiNLJLlwz4EnDV9tfxQ3C5W5OHpH051%2FspHr7t202Vsm67lpY9BeVRsP%2Fa5M88VwEm%2B33vkgHuLBJGjo6iUqWXPmwftWMuoVnIE21NeFCDQcx1ZG%2BHgcHPoZtiNXrPuOIjXkFaE8KZDgrP3ZBDBmuHbkdmNVogi4D0ktOeF6qZ5oGj0P60evk1y6bOzdMI8TKG7ooGr8eSUquUTs650LlmIFp6LhRYvMsOE%2BY7MBw%2FR0gdnGH1DgX2HpsqnxT7e3l5eKh3atfKbpbudAcPnxPjbNBUzQ1Sf3z1i%2F%2B6%2BdrN2Y8IXcCoTkUy6uAvv%2Btd0KN8vLWI4ALlJdE0I5RsuDwvZTIsRGi6puzS236ZqXL9AXx7LVjb1pl2dgpuw444g8McLAw4R%2Fe4pV9XPcgL%2F9ej1WVL8D9%2Fihrw3ZImvyRADW2vvdpX8LEy1qJfNFPc9Vue8ew9E4Hsuzj9s3oXqx6h86%2BxPmzqjn%2BXFdq3E%2FaW9S5ASfGy2JSe%2F4qKbYOlMcUgJ7a0okR%2B1YtfWQmsHwu%2F45DPe319AReeE%2Botl2oF9vLMwrOCcngY6qQHszJeuSQBMl1I6GrtV3Y2MwNz%2BJbrdmVtaAbC%2B67zSDbb7vi7gXTUuFgfFF8a7hmN3c%2Fjx9NqNkwyerH6Bn9Wu3Ygb%2FR4u7DmLRlBQqObH6XAd%2Fn4pIZxzashT%2F8LuvarpQIuwtgEZ3aFcErks9ehP3V4MUl7sbYDZDbMHTvfLYBB22oy0A%2FTTv6fNWJwroeB%2FfBo%2FJZsIm2MZx2HEbdEOLq7hwTMj1%2FOu&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230118T003549Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYTIYNG23C%2F20230118%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=f63a7cc89c4dae307f1e185c9604a5c0ce862d24c34f65210b309c0966bba7cf&hash=246e82663ea6e9aa73da731f12522e5d3d39019c9c27582a6002504a89cbf087&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0020708X83902752&tid=spdf-9dfbe0e8-84a1-4eb1-a3d0-ef659e4b9233&sid=1a7d585752f89741f17a0ad0d0eebebbf064gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=131052590c5c51570156&rr=78b340b8a85d294c&cc=us
    #
    # FIX TH-232
    #—--------------------------------------------
    #
    #Subset of peaks, check to see if peaks are feasible to detect
    #Need to redo library a bit
    #–
    #
    #0.2% enriched DU
    #Need to check U235 lines (1001, 700)
    #
    #HEU 238 and 235
    #Height of 235 peaks to determine if DU or HEU
    #
    #Peak ratio of 186 keV to 10001 keV
    #(see LOS ALAMOS)
    #https://www.lanl.gov/org/ddste/aldgs/sst-training/_assets/docs/PANDA/The%20Measurement%20of%20Uranimum%20Enrichment%20Ch.%207%20p.%20195-220.pdf
    #
    #–
    #
    #WGPU mostly 239 peaks, 
    #See useful gammas in below link
    #
    #https://www.lanl.gov/org/ddste/aldgs/sst-training/_assets/docs/PANDA/Plutonium%20Isotopic%20Composition%20by%20Gamma-Ray%20Spectroscopy%20Ch.%208%20p.%20221-272.pdf
    #
    #–
    #
    #Width should just be 1% of E or if < 2kev, set to 2
    #
    #PROMINENCE AND THRESHOLD for scipy find peaks
    #Need to do more parameter tuing
    #
    #Maybe try smooth? But if filter less than FWHM then peak will be lowered
    #
    #Peak centroid could be 0.05% different than “true” peak
    #
    #–
    #If peak matches, add to “score” of an isotope
    #
    #Maybe make it so if matches really well, >1
    #If ambiguous, add <1
    #etc.
    #


    # Need to do:
    # only plot certain branching ratios
    # re-look at all of the peaks and branching ratios of isotopes in dataset (16 isotopes)
    # tune parameters based on counts
    # write something to score each peak based on true peaks
    # after doing that contiune to tune parameters
    # then work on identification...