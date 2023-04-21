#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from numpy import rot90

from BIF import BIFReader
import pdb

bif_filename = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20221102-154955_ImagingTest(-1,0)\\20221102-160044_SS\\RadImage.bif"
bif_filename2 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20221102-153230_ImagingTest(0,1)\\20221102-154242_SS\\RadImage.bif"
bif_filename3 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20221102-154453_ImagingTest(1,0)\\20221102-154528_SS\\RadImage.bif"
image = BIFReader(bif_filename)
image2 = BIFReader(bif_filename2)
image3 = BIFReader(bif_filename3)
# Lets collapse the polar bins into one single bin per azimuthal direction
image_t = image.compton.T
image2_t = image2.compton.T
image3_t = image3.compton.T

new_image = []
for polar_bin in image_t:
    new_image.append(np.sum(polar_bin))
new_image = np.array(new_image)
normalized_line_image = new_image / np.sum(new_image)
spectrum = plt.plot(normalized_line_image)
xvalues = spectrum[0].get_xdata()
yvalues = spectrum[0].get_ydata()
probability_at_degree = []
for i in range (0,180):

    idx = np.where(xvalues==xvalues[i])
    ylocation = yvalues[idx]
    probability_at_degree.append(ylocation)

arr_1 = np.array(probability_at_degree)
list_1 = arr_1.tolist()
rad_data_test = np.savetxt("rad_data_test_1.txt", arr_1)

new_image2 = []
for polar_bin in image2_t:
    new_image2.append(np.sum(polar_bin))
new_image2 = np.array(new_image2)
normalized_line_image2 = new_image2 / np.sum(new_image2)
spectrum2 = plt.plot(normalized_line_image2)
xvalues2 = spectrum2[0].get_xdata()
yvalues2 = spectrum2[0].get_ydata()
probability_at_degree2 = []
for i in range (0,180):
    #probability_at_degree = i
    idx2 = np.where(xvalues2==xvalues2[i])
    ylocation2 = yvalues2[idx2]
    probability_at_degree2.append(ylocation2)
    #print(ylocation)
arr_2 = np.array(probability_at_degree2)
list_1= arr_1.tolist()
np.savetxt("rad_data2.txt", arr_2)
#print(probability_at_degree[0])

#firstvalue = arr_1[0]

#print(firstvalue)

matrixmap = np.zeros((20,20),dtype=float)
#o degrees
for i in range (0,10):
    matrixmap[i][0] = arr_1[0]
#90 degrees
for i in range(0,20):
    matrixmap[10][i] = arr_1[45]
#180 degrees
for i in range(11,20):
    matrixmap[i][0] = arr_1[90]

#0 degrees
matrixmap[0][0] = arr_1[1]
matrixmap[0][1] = arr_1[1]


#5 degrees
matrixmap[0][3] = arr_1[3]
matrixmap[0][2] = arr_1[3]
matrixmap[1][1] = arr_1[3]
matrixmap[1][0] = arr_1[3]
#10 degrees
matrixmap[0][5] = arr_1[5]
matrixmap[0][4] = arr_1[5]
matrixmap[1][3] = arr_1[5]
matrixmap[1][2] = arr_1[5]
matrixmap[2][1] = arr_1[5]
matrixmap[2][0] = arr_1[5]
#15 degrees
matrixmap[0][7] = arr_1[8]
matrixmap[0][6] = arr_1[8]
matrixmap[1][5] = arr_1[8]
matrixmap[1][4] = arr_1[8]
matrixmap[2][3] = arr_1[8]
matrixmap[2][2] = arr_1[8]
matrixmap[3][1] = arr_1[8]
matrixmap[3][0] = arr_1[8]

#20 degrees
matrixmap[0][9] = arr_1[10]
matrixmap[0][8] = arr_1[10]
matrixmap[1][7] = arr_1[10]
matrixmap[1][6] = arr_1[10]
matrixmap[2][5] = arr_1[10]
matrixmap[2][4] = arr_1[10]
matrixmap[3][3] = arr_1[10]
matrixmap[3][2] = arr_1[10]
matrixmap[4][1] = arr_1[10]
matrixmap[4][0] = arr_1[10]

#25 degrees
matrixmap[0][11] = arr_1[13]
matrixmap[0][10] = arr_1[13]
matrixmap[1][9] = arr_1[13]
matrixmap[1][8] = arr_1[13]
matrixmap[2][7] = arr_1[13]
matrixmap[2][6] = arr_1[13]
matrixmap[3][5] = arr_1[13]
matrixmap[3][4] = arr_1[13]
matrixmap[4][3] = arr_1[13]
matrixmap[4][2] = arr_1[13]
matrixmap[5][1] = arr_1[13]
matrixmap[5][0] = arr_1[13]

#30 degrees
matrixmap[0][13] = arr_1[15]
matrixmap[0][12] = arr_1[15]
matrixmap[1][11] = arr_1[15]
matrixmap[1][10] = arr_1[15]
matrixmap[2][9] = arr_1[15]
matrixmap[2][8] = arr_1[15]
matrixmap[3][7] = arr_1[15]
matrixmap[3][6] = arr_1[15]
matrixmap[4][5] = arr_1[15]
matrixmap[4][4] = arr_1[15]
matrixmap[5][3] = arr_1[15]
matrixmap[5][2] = arr_1[15]
matrixmap[6][1] = arr_1[15]
matrixmap[6][0] = arr_1[15]

#35 degrees
matrixmap[0][15] = arr_1[18]
matrixmap[0][14] = arr_1[18]
matrixmap[1][13] = arr_1[18]
matrixmap[1][12] = arr_1[18]
matrixmap[2][11] = arr_1[18]
matrixmap[2][10] = arr_1[18]
matrixmap[3][9] = arr_1[18]
matrixmap[3][8] = arr_1[18]
matrixmap[4][7] = arr_1[18]
matrixmap[4][6] = arr_1[18]
matrixmap[5][5] = arr_1[18]
matrixmap[5][4] = arr_1[18]
matrixmap[6][3] = arr_1[18]
matrixmap[6][2] = arr_1[18]
matrixmap[7][1] = arr_1[18]
matrixmap[7][0] = arr_1[18]


#40 degrees
matrixmap[0][17] = arr_1[20]
matrixmap[0][16] = arr_1[20]
matrixmap[1][15] = arr_1[20]
matrixmap[1][14] = arr_1[20]
matrixmap[2][13] = arr_1[20]
matrixmap[2][12] = arr_1[20]
matrixmap[3][11] = arr_1[20]
matrixmap[3][10] = arr_1[20]
matrixmap[4][9] = arr_1[20]
matrixmap[4][8] = arr_1[20]
matrixmap[5][7] = arr_1[20]
matrixmap[5][6] = arr_1[20]
matrixmap[6][5] = arr_1[20]
matrixmap[6][4] = arr_1[20]
matrixmap[7][3] = arr_1[20]
matrixmap[7][2] = arr_1[20]
matrixmap[8][1] = arr_1[20]
matrixmap[8][0] = arr_1[20]

#45 degrees
matrixmap[0][19] = arr_1[23]
matrixmap[0][18] = arr_1[23]
matrixmap[1][17] = arr_1[23]
matrixmap[1][16] = arr_1[23]
matrixmap[2][15] = arr_1[23]
matrixmap[2][14] = arr_1[23]
matrixmap[3][13] = arr_1[23]
matrixmap[3][12] = arr_1[23]
matrixmap[4][11] = arr_1[23]
matrixmap[4][10] = arr_1[23]
matrixmap[5][9] = arr_1[23]
matrixmap[5][8] = arr_1[23]
matrixmap[6][7] = arr_1[23]
matrixmap[6][6] = arr_1[23]
matrixmap[7][5] = arr_1[23]
matrixmap[7][4] = arr_1[23]
matrixmap[8][3] = arr_1[23]
matrixmap[8][2] = arr_1[23]
matrixmap[9][1] = arr_1[23]
matrixmap[9][0] = arr_1[23]

#50 degrees
matrixmap[1][19] = arr_1[25]
matrixmap[1][18] = arr_1[25]
matrixmap[2][17] = arr_1[25]
matrixmap[2][16] = arr_1[25]
matrixmap[3][15] = arr_1[25]
matrixmap[3][14] = arr_1[25]
matrixmap[4][13] = arr_1[25]
matrixmap[4][12] = arr_1[25]
matrixmap[5][11] = arr_1[25]
matrixmap[5][10] = arr_1[25]
matrixmap[6][9] = arr_1[25]
matrixmap[6][8] = arr_1[25]
matrixmap[7][7] = arr_1[25]
matrixmap[7][6] = arr_1[25]
matrixmap[8][5] = arr_1[25]
matrixmap[8][4] = arr_1[25]
matrixmap[9][3] = arr_1[25]
matrixmap[9][2] = arr_1[25]

#55 degrees
matrixmap[2][19] = arr_1[28]
matrixmap[2][18] = arr_1[28]
matrixmap[3][17] = arr_1[28]
matrixmap[3][16] = arr_1[28]
matrixmap[4][15] = arr_1[28]
matrixmap[4][14] = arr_1[28]
matrixmap[5][13] = arr_1[28]
matrixmap[5][12] = arr_1[28]
matrixmap[6][11] = arr_1[28]
matrixmap[6][10] = arr_1[28]
matrixmap[7][9] = arr_1[28]
matrixmap[7][8] = arr_1[28]
matrixmap[8][7] = arr_1[28]
matrixmap[8][6] = arr_1[28]
matrixmap[9][5] = arr_1[28]
matrixmap[9][4] = arr_1[28]

#60 degrees
matrixmap[3][19] = arr_1[30]
matrixmap[3][18] = arr_1[30]
matrixmap[4][17] = arr_1[30]
matrixmap[4][16] = arr_1[30]
matrixmap[5][15] = arr_1[30]
matrixmap[5][14] = arr_1[30]
matrixmap[6][13] = arr_1[30]
matrixmap[6][12] = arr_1[30]
matrixmap[7][11] = arr_1[30]
matrixmap[7][10] = arr_1[30]
matrixmap[8][9] = arr_1[30]
matrixmap[8][8] = arr_1[30]
matrixmap[9][7] = arr_1[30]
matrixmap[9][6] = arr_1[30]

#65 degrees
matrixmap[4][19] = arr_1[33]
matrixmap[4][18] = arr_1[33]
matrixmap[5][17] = arr_1[33]
matrixmap[5][16] = arr_1[33]
matrixmap[6][15] = arr_1[33]
matrixmap[6][14] = arr_1[33]
matrixmap[7][13] = arr_1[33]
matrixmap[7][12] = arr_1[33]
matrixmap[8][11] = arr_1[33]
matrixmap[8][10] = arr_1[33]
matrixmap[9][9] = arr_1[33]
matrixmap[9][8] = arr_1[33]

#70 degrees
matrixmap[5][19] = arr_1[35]
matrixmap[5][18] = arr_1[35]
matrixmap[6][17] = arr_1[35]
matrixmap[6][16] = arr_1[35]
matrixmap[7][15] = arr_1[35]
matrixmap[7][14] = arr_1[35]
matrixmap[8][13] = arr_1[35]
matrixmap[8][12] = arr_1[35]
matrixmap[9][11] = arr_1[35]
matrixmap[9][10] = arr_1[35]

#75 degrees
matrixmap[6][19] = arr_1[38]
matrixmap[6][18] = arr_1[38]
matrixmap[7][17] = arr_1[38]
matrixmap[7][16] = arr_1[38]
matrixmap[8][15] = arr_1[38]
matrixmap[8][14] = arr_1[38]
matrixmap[9][13] = arr_1[38]
matrixmap[9][12] = arr_1[38]

#80 degrees
matrixmap[7][19] = arr_1[40]
matrixmap[7][18] = arr_1[40]
matrixmap[8][17] = arr_1[40]
matrixmap[8][16] = arr_1[40]
matrixmap[9][15] = arr_1[40]
matrixmap[9][14] = arr_1[40]

#85 degrees
matrixmap[8][19] = arr_1[43]
matrixmap[8][18] = arr_1[43]
matrixmap[9][17] = arr_1[43]
matrixmap[9][16] = arr_1[43]

#90 degrees
matrixmap[9][19] = arr_1[45]
matrixmap[9][18] = arr_1[45]

#95 degrees
matrixmap[11][19] = arr_1[48]
matrixmap[11][18] = arr_1[48]
matrixmap[10][17] = arr_1[48]
matrixmap[10][16] = arr_1[48]

#100 degrees
matrixmap[12][19] = arr_1[50]
matrixmap[12][18] = arr_1[50]
matrixmap[11][17] = arr_1[50]
matrixmap[11][16] = arr_1[50]
matrixmap[10][15] = arr_1[50]
matrixmap[10][14] = arr_1[50]

#105 degrees
matrixmap[13][19] = arr_1[53]
matrixmap[13][18] = arr_1[53]
matrixmap[12][17] = arr_1[53]
matrixmap[12][16] = arr_1[53]
matrixmap[11][15] = arr_1[53]
matrixmap[11][14] = arr_1[53]
matrixmap[10][13] = arr_1[53]
matrixmap[10][12] = arr_1[53]

#110 degrees
matrixmap[14][19] = arr_1[55]
matrixmap[14][18] = arr_1[55]
matrixmap[13][17] = arr_1[55]
matrixmap[13][16] = arr_1[55]
matrixmap[12][15] = arr_1[55]
matrixmap[12][14] = arr_1[55]
matrixmap[11][13] = arr_1[55]
matrixmap[11][12] = arr_1[55]
matrixmap[10][11] = arr_1[55]
matrixmap[10][10] = arr_1[55]
#115 degrees
matrixmap[15][19] = arr_1[58]
matrixmap[15][18] = arr_1[58]
matrixmap[14][17] = arr_1[58]
matrixmap[14][16] = arr_1[58]
matrixmap[13][15] = arr_1[58]
matrixmap[13][14] = arr_1[58]
matrixmap[12][13] = arr_1[58]
matrixmap[12][12] = arr_1[58]
matrixmap[11][11] = arr_1[58]
matrixmap[11][10] = arr_1[58]
matrixmap[10][9] = arr_1[58]
matrixmap[10][8] = arr_1[58]

#120 degrees
matrixmap[16][19] = arr_1[60]
matrixmap[16][18] = arr_1[60]
matrixmap[15][17] = arr_1[60]
matrixmap[15][16] = arr_1[60]
matrixmap[14][15] = arr_1[60]
matrixmap[14][14] = arr_1[60]
matrixmap[13][13] = arr_1[60]
matrixmap[13][12] = arr_1[60]
matrixmap[12][11] = arr_1[60]
matrixmap[12][10] = arr_1[60]
matrixmap[11][9] = arr_1[60]
matrixmap[11][8] = arr_1[60]
matrixmap[10][7] = arr_1[60]
matrixmap[10][6] = arr_1[60]

#125 degrees
matrixmap[17][19] = arr_1[63]
matrixmap[17][18] = arr_1[63]
matrixmap[16][17] = arr_1[63]
matrixmap[16][16] = arr_1[63]
matrixmap[15][15] = arr_1[63]
matrixmap[15][14] = arr_1[63]
matrixmap[14][13] = arr_1[63]
matrixmap[14][12] = arr_1[63]
matrixmap[13][11] = arr_1[63]
matrixmap[13][10] = arr_1[63]
matrixmap[12][9] = arr_1[63]
matrixmap[12][8] = arr_1[63]
matrixmap[11][7] = arr_1[63]
matrixmap[11][6] = arr_1[63]
matrixmap[10][5] = arr_1[63]
matrixmap[10][4] = arr_1[63]
#130 degrees
matrixmap[18][19] = arr_1[65]
matrixmap[18][18] = arr_1[65]
matrixmap[17][17] = arr_1[65]
matrixmap[17][16] = arr_1[65]
matrixmap[16][15] = arr_1[65]
matrixmap[16][14] = arr_1[65]
matrixmap[15][13] = arr_1[65]
matrixmap[15][12] = arr_1[65]
matrixmap[14][11] = arr_1[65]
matrixmap[14][10] = arr_1[65]
matrixmap[13][9] = arr_1[65]
matrixmap[13][8] = arr_1[65]
matrixmap[12][7] = arr_1[65]
matrixmap[12][6] = arr_1[65]
matrixmap[11][5] = arr_1[65]
matrixmap[11][4] = arr_1[65]
matrixmap[10][3] = arr_1[65]
matrixmap[10][2] = arr_1[65]
#135 degrees
matrixmap[19][19] = arr_1[68]
matrixmap[19][18] = arr_1[68]
matrixmap[18][17] = arr_1[68]
matrixmap[18][16] = arr_1[68]
matrixmap[17][15] = arr_1[68]
matrixmap[17][14] = arr_1[68]
matrixmap[16][13] = arr_1[68]
matrixmap[16][12] = arr_1[68]
matrixmap[15][11] = arr_1[68]
matrixmap[15][10] = arr_1[68]
matrixmap[14][9] = arr_1[68]
matrixmap[14][8] = arr_1[68]
matrixmap[13][7] = arr_1[68]
matrixmap[13][6] = arr_1[68]
matrixmap[12][5] = arr_1[68]
matrixmap[12][4] = arr_1[68]
matrixmap[11][3] = arr_1[68]
matrixmap[11][2] = arr_1[68]
matrixmap[10][1] = arr_1[68]
matrixmap[10][0] = arr_1[68]

#140 degrees
matrixmap[19][17] = arr_1[70]
matrixmap[19][16] = arr_1[70]
matrixmap[18][15] = arr_1[70]
matrixmap[18][14] = arr_1[70]
matrixmap[17][13] = arr_1[70]
matrixmap[17][12] = arr_1[70]
matrixmap[16][11] = arr_1[70]
matrixmap[16][10] = arr_1[70]
matrixmap[15][9] = arr_1[70]
matrixmap[15][8] = arr_1[70]
matrixmap[14][7] = arr_1[70]
matrixmap[14][6] = arr_1[70]
matrixmap[13][5] = arr_1[70]
matrixmap[13][4] = arr_1[70]
matrixmap[12][3] = arr_1[70]
matrixmap[12][2] = arr_1[70]
matrixmap[11][1] = arr_1[70]
matrixmap[11][0] = arr_1[70]

#145 degrees
matrixmap[19][15] = arr_1[73]
matrixmap[19][14] = arr_1[73]
matrixmap[18][13] = arr_1[73]
matrixmap[18][12] = arr_1[73]
matrixmap[17][11] = arr_1[73]
matrixmap[17][10] = arr_1[73]
matrixmap[16][9] = arr_1[73]
matrixmap[16][8] = arr_1[73]
matrixmap[15][7] = arr_1[73]
matrixmap[15][6] = arr_1[73]
matrixmap[14][5] = arr_1[73]
matrixmap[14][4] = arr_1[73]
matrixmap[13][3] = arr_1[73]
matrixmap[13][2] = arr_1[73]
matrixmap[12][1] = arr_1[73]
matrixmap[12][0] = arr_1[73]
#150 degrees
matrixmap[19][13] = arr_1[75]
matrixmap[19][12] = arr_1[75]
matrixmap[18][11] = arr_1[75]
matrixmap[18][10] = arr_1[75]
matrixmap[17][9] = arr_1[75]
matrixmap[17][8] = arr_1[75]
matrixmap[16][7] = arr_1[75]
matrixmap[16][6] = arr_1[75]
matrixmap[15][5] = arr_1[75]
matrixmap[15][4] = arr_1[75]
matrixmap[14][3] = arr_1[75]
matrixmap[14][2] = arr_1[75]
matrixmap[13][1] = arr_1[75]
matrixmap[13][0] = arr_1[75]

#155 degrees
matrixmap[19][11] = arr_1[78]
matrixmap[19][10] = arr_1[78]
matrixmap[18][9] = arr_1[78]
matrixmap[18][8] = arr_1[78]
matrixmap[17][7] = arr_1[78]
matrixmap[17][6] = arr_1[78]
matrixmap[16][5] = arr_1[78]
matrixmap[16][4] = arr_1[78]
matrixmap[15][3] = arr_1[78]
matrixmap[15][2] = arr_1[78]
matrixmap[14][1] = arr_1[78]
matrixmap[14][0] = arr_1[78]

#160 degrees
matrixmap[19][9] = arr_1[80]
matrixmap[19][8] = arr_1[80]
matrixmap[18][7] = arr_1[80]
matrixmap[18][6] = arr_1[80]
matrixmap[17][5] = arr_1[80]
matrixmap[17][4] = arr_1[80]
matrixmap[16][3] = arr_1[80]
matrixmap[16][2] = arr_1[80]
matrixmap[15][1] = arr_1[80]
matrixmap[15][0] = arr_1[80]

#165 degrees
matrixmap[19][7] = arr_1[83]
matrixmap[19][6] = arr_1[83]
matrixmap[18][5] = arr_1[83]
matrixmap[18][4] = arr_1[83]
matrixmap[17][3] = arr_1[83]
matrixmap[17][2] = arr_1[83]
matrixmap[16][1] = arr_1[83]
matrixmap[16][0] = arr_1[83]

#170 degrees
matrixmap[19][5] = arr_1[85]
matrixmap[19][4] = arr_1[85]
matrixmap[18][3] = arr_1[85]
matrixmap[18][2] = arr_1[85]
matrixmap[17][1] = arr_1[85]
matrixmap[17][0] = arr_1[85]

#175 degrees
matrixmap[19][3] = arr_1[88]
matrixmap[19][2] = arr_1[88]
matrixmap[18][1] = arr_1[88]
matrixmap[18][0] = arr_1[88]

#180 degrees
matrixmap[19][0:19] = arr_1[90]
#matrixmap[19][0] = arr_1[90]


#print(matrixmap)
plt.matshow(matrixmap, cmap=plt.cm.jet)
plt.show()

newmatrix = rot90(matrixmap, k=1, axes=(0,1))
plt.matshow(newmatrix, cmap=plt.cm.jet)

addedmatrix = matrixmap + newmatrix
plt.matshow(addedmatrix, cmap=plt.cm.jet)

newmatrix2 = rot90(newmatrix, k=1, axes=(0,1))

#combinedmatrix= np.zeros((40,60),dtype=float)
#combinedmatrix[10:10+matrixmap.shape[0], 0:0+matrixmap.shape[1]] +=matrixmap
#combinedmatrix[20:20+newmatrix.shape[0], 20:20+newmatrix.shape[1]] +=newmatrix
#combinedmatrix[10:10+newmatrix2.shape[0], 59:59+newmatrix2.shape[1]] +=newmatrix2
#plt.matshow(combinedmatrix, cmap=plt.cm.jet)

#newmatrix2 = rot90(newmatrix, k=1, axes=(0,1))
plt.matshow(newmatrix2, cmap=plt.cm.jet)

newmatrix3 = rot90(newmatrix2, k=1, axes=(0,1))
plt.matshow(newmatrix3, cmap=plt.cm.jet)


#matrixmap_01 = matrixmap.transpose()

#plt.matshow(matrixmap_01,cmap = plt.cm.jet)

#matrixmap_10 = matrixmap_01.transpose()
#plt.matshow(matrixmap_10,cmap=plt.cm.jet)
#print(matrixmap)
x,y,z = [0, 0 ,0]
yaw = 0

bin_width = 0.05 # m

x_bins = 100
y_bins = 100

map_image = np.zeros((y_bins, x_bins))

map_image[20][20] = 100

plt.figure()
plt.imshow(map_image)
plt.show()

# pdb.set_trace()
# tets= 1