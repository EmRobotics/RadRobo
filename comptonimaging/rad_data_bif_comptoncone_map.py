#!/usr/bin/python3
import textwrap

import numpy as np
import matplotlib.pyplot as plt

from BIF import BIFReader
import pdb

bif1 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-175334_SS\\RadImage.bif"
bif2 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-175916_SS\\RadImage.bif"
bif3 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-180602_SS\\RadImage.bif"
bif4 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-181641_SS\\RadImage.bif"
bif5 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-182351_SS\\RadImage.bif"
def rot_mat(deg):
    theta = deg / 180 * np.pi
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, -s], [s, c]])



def load_rad_image(bif_filename):
    image = BIFReader(bif_filename)

    # Lets collapse the polar bins into one single bin per azimuthal direction
    image_t = image.compton.T
    new_image = []
    for polar_bin in image_t:
        new_image.append(np.sum(polar_bin))

    new_image = np.array(new_image)

    normalized_line_image = new_image / np.sum(new_image)

    spectrum = plt.plot(normalized_line_image)
    xvalues = spectrum[0].get_xdata()
    yvalues = spectrum[0].get_ydata()
    probability_at_degree = []
    for i in range(0, 180):
        idx = np.where(xvalues == xvalues[i])
        ylocation = yvalues[idx]
        probability_at_degree.append(ylocation)

    arr_1 = np.array(probability_at_degree)
    list_1 = arr_1.tolist()
    # np.savetxt("rad_data_test.txt", arr_1)

    return normalized_line_image, list_1


rad_data_1 = np.loadtxt('rad_data_april_7_1.txt')
rad_data_2 = np.loadtxt('rad_data_april_7_2.txt')
rad_data_3 = np.loadtxt('rad_data_april_7_3.txt')
rad_data_4 = np.loadtxt('rad_data_april_7_4.txt')
rad_data_5 = np.loadtxt('rad_data_april_7_5.txt')

# rad_data_3 = np.loadtxt('rad_data_april_5_3.txt')
# rad_data_4 = np.loadtxt('rad_data_april_5_4.txt')
# rad_data_5 = np.loadtxt('rad_data_april_5_5.txt')
def covert_to_index(x_min, x_max, y_min, y_max, x_size, y_size, pos):
    x_ind = x_size - ((x_max - pos[0]) / (x_max - x_min)) * x_size
    y_ind = y_size - ((y_max - pos[1]) / (y_max - y_min)) * y_size

    return int(y_ind), int(x_ind)

def project_image(x_size, y_size, start_point, yaw, rad_data):
    projected_image = np.zeros((y_size, x_size))

    for y_ind in range(y_size):
        for x_ind in range(x_size):

            deg = np.arctan2(x_ind - start_point[1], y_ind - start_point[0]) * 180 / np.pi

            deg += yaw

            if deg >= 360:
                deg -= 360

            int_deg = int(deg / 2)

            try:
                projected_image[y_ind, x_ind] = rad_data[int_deg]
            except IndexError:
                pdb.set_trace()
                test = 1

    return projected_image


# we have four measurements.  (x (m), y (m)) we want to convert to integer positions based on a total grid size (arbitrary)
# We want to get (y, x) indices in span position

# total grid size

x_size = 496
y_size = 353

y_min = -11.8472
y_max = y_min+0.05*y_size
x_min = -12.0194
x_max = x_min+0.05*x_size

bif1 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-175334_SS\\RadImage.bif"
bif2 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-175916_SS\\RadImage.bif"
bif3 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-180602_SS\\RadImage.bif"
bif4 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-181641_SS\\RadImage.bif"
bif5 = "C:\\Users\\katie\\PycharmProjects\\NERS491Project\\Compton Data\\ComptonImaging\\20230407-182351_SS\\RadImage.bif"

raw_img1 = load_rad_image(bif1)
yaw1 = 0
img_pos1 = [0.66307, -0.0843]


raw_img2 = load_rad_image(bif2)
yaw2 = 85.5662
img_pos2 = [1.86653, -1.35223]

raw_img3 = load_rad_image(bif3)
yaw3 = 182.834
img_pos3 = [3.34589, .059664]

raw_img4 = load_rad_image(bif4)
yaw4 = 276.272
img_pos4 = [1.89838, .985561]

raw_img5 = load_rad_image(bif5)
yaw5 = 181.435
img_pos5 = [.64764, .111667]

start_point_1 = covert_to_index(x_min, x_max, y_min, y_max, x_size, y_size, img_pos1)
start_point_2 = covert_to_index(x_min, x_max, y_min, y_max, x_size, y_size, img_pos2)
start_point_3 = covert_to_index(x_min, x_max, y_min, y_max, x_size, y_size, img_pos3)
start_point_4 = covert_to_index(x_min, x_max, y_min, y_max, x_size, y_size, img_pos4)
start_point_5 = covert_to_index(x_min, x_max, y_min, y_max, x_size, y_size, img_pos5)


img1 = project_image(x_size, y_size, start_point_1, 0, rad_data_1)
img2 = project_image(x_size, y_size, start_point_2, 85.5662, rad_data_2)
img3 = project_image(x_size, y_size, start_point_2, 182.834, rad_data_3)
img4 = project_image(x_size, y_size, start_point_2, 276.272, rad_data_4)
img5 = project_image(x_size, y_size, start_point_2, 181.435, rad_data_5)

cmap = plt.jet()

plt.imshow(img1)
plt.gca().invert_yaxis()
plt.title('Source Localization Map 1')
plt.xlabel('X dimension')
plt.ylabel('Y dimension')
plt.show()

#img2_max = np.max(img2)
# img2[img2 < (0.9*img2_max)] = 0
#normalized_image_matrix[normalized_image_matrix < (0.75*normalized_image_matrix_max)] = 0
plt.imshow(img2)
plt.gca().invert_yaxis()
plt.title('Source Localization Map 2')
plt.xlabel('X dimension')
plt.ylabel('Y dimension')
plt.show()

plt.imshow(img3)
plt.gca().invert_yaxis()
plt.title('Source Localization Map 3')
plt.xlabel('X dimension')
plt.ylabel('Y dimension')
plt.show()

plt.imshow(img4)
plt.gca().invert_yaxis()
plt.title('Source Localization Map 4')
plt.xlabel('X dimension')
plt.ylabel('Y dimension')
plt.show()

# img5[0:353, 0:200]= 0
# img5[0:170, 200:496] = 0
# img5_max = np.max(img5)
# img5[img5 < (0.9*img5_max)] = 0
plt.imshow(img5)
plt.gca().invert_yaxis()
plt.title('Source Localization Map 5')
plt.xlabel('X dimension')
plt.ylabel('Y dimension')
plt.show()

image_matrix = img1 + img2
plt.imshow(img1+img2)
#x = [0,100,200,300,400,496]
#xlabels = [0,5,10,15,20,24.8]
#y = [0,50,100,150,200,250,300,353]
#ylabels = [0,2.5,5,7.5,10,12.5,15,17.65]
# plt.xticks(x, xlabels)
# plt.yticks(y, ylabels)
plt.gca().invert_yaxis()
plt.title('Source Localization Map')
plt.xlabel('X dimension (m)')
plt.ylabel('Y dimension (m)')
plt.savefig('sourcelocalization_pixelated.png')

plt.show()

print('start point 1',start_point_1)
print('start point 2', start_point_2)
print('start point 3',start_point_3)
print('start point 4', start_point_4)
print('start point 5', start_point_5)




cmap = plt.gray()

print('The image matrix size is:', image_matrix.shape)
image_matrix_min = np.min(image_matrix)
print('Minimum value of Image Matrix:',image_matrix_min)
image_matrix_max = np.max(image_matrix)
print('Maximum value of Image Matrix:',image_matrix_max)
image_matrix_std = np.std(image_matrix)


print('Standard Deviation of Image Matrix:',image_matrix_std)

def normalize(x):
   return (x - image_matrix_min) / (image_matrix_max - image_matrix_min)


normalized_image_matrix = normalize(image_matrix)

normalized_image_matrix_min = np.min(normalized_image_matrix)
print('Minimum value of Image Matrix:',normalized_image_matrix_min)
normalized_image_matrix_max = np.max(normalized_image_matrix)
print('Maximum value of Image Matrix:',normalized_image_matrix_max)
normalized_image_matrix[normalized_image_matrix < (0.75*normalized_image_matrix_max)] = 0
print(normalized_image_matrix)
from scipy import ndimage
source_location = ndimage.center_of_mass(normalized_image_matrix)
plt.imshow(normalized_image_matrix)
plt.gca().invert_yaxis()
# x = [0,100,200,300,400,496]
# xlabels = [0,5,10,15,20,24.8]
# y = [0,50,100,150,200,250,300,353]
#ylabels = [0,2.5,5,7.5,10,12.5,15,17.65]
#plt.xticks(x, xlabels)
# plt.yticks(y, ylabels)
plt.title('Source Localization Map')
plt.xlabel('X Dimension (m)')
plt.ylabel('Y Dimension (m)')
plt.colorbar(label="Source Probability", orientation="horizontal")
plt.savefig('NormalizedImage_location_correctdims2.png')
plt.show()
plt.show()

print('The source is located at (row,column) position:', source_location)
source_locationx = source_location[1]*0.05
source_locationy = source_location[0]*0.05
print(source_locationx)
print(source_locationy)

#this tests to see where the source is relative to the detector
x_ind = source_location[1]
y_ind = source_location[0]
deg_test = np.arctan2(x_ind - start_point_5[1], y_ind - start_point_5[0]) * 180 / np.pi
deg_test_1 = deg_test + yaw5
print(deg_test_1)

