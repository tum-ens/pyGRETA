import numpy as np
import time

print("Started...")
t_0 = time.time()

xmin = 3
xmax = 8
ymin = 1
ymax = 6

GWA_array = np.random.rand(3999, 8052)*10
x_gwa = np.random.rand(3999, 8052)*10
y_gwa = np.random.rand(3999, 8052)*10

gwa_rows, gwa_cols = GWA_array.shape
gwa_value = np.zeros([gwa_rows, gwa_cols])

total_num_cells = 0
value_num_cells = 0

for k in range(gwa_rows):
    for l in range(gwa_cols):
        if xmin <= x_gwa[k,l] < xmax and ymin <= y_gwa[k,l] < ymax:
            gwa_value[k,l] = GWA_array[k,l]
            total_num_cells = total_num_cells + 1
            if gwa_value[k,l]:
                value_num_cells = value_num_cells+1
            if total_num_cells == 1:
                K_first = k
                L_first = l
            K_last = k
            L_last = l

t_1 = time.time()
print("1) Loop: " + str(t_1-t_0))


# 2) Selection
selection_index = (xmin <= x_gwa) & (x_gwa < xmax) & (ymin <= y_gwa) & (y_gwa < ymax)  # Determine the pixels that are inbetween the range
GWA_array[np.invert(selection_index)] = 0  # Set pixel not covered by the shapfile to zero
value_num_cells = np.count_nonzero(GWA_array)  # Number of non-zero pixels
# print (value_num_cells)

coordinates_nonzero_pixels_x, coordinates_nonzero_pixels_y = selection_index.nonzero()  # Determine the coordinates of non-zero pixels in order to determine the outer rows and columns. Tuple(x,y)
K_first = coordinates_nonzero_pixels_x[0]  # First x-coordinate of non-zero pixels
L_first = coordinates_nonzero_pixels_y[0]  # First y-coordinate of non-zero pixels
K_last = coordinates_nonzero_pixels_x[-1]  # Last x-coordinate of non-zero pixels
L_last = coordinates_nonzero_pixels_y[-1]  # Last y-coordinate of non-zero pixels

t_2 = time.time()
print("1) Indexing (Numpy): " + str(t_2-t_1))









i_offset = 0
j_offset = 0

K_first = 0
K_last = 199
L_last = 249
L_first = 0

value_num_cells = 5
Num_pix = 10

reMerra = np.zeros([200, 250, 8760])
merraData = np.random.rand(8760)*10
gwa_cut = GWA_array[K_first:K_last+1, L_first:L_last+1]


for h in range(8760):
    V = merraData[h]
    E = V**3
    E = E*value_num_cells/Num_pix
    merra_cut = gwa_cut/np.sum(gwa_cut)*E
    merra_cut = np.cbrt(merra_cut*Num_pix)
    reMerra[i_offset:i_offset+K_last+1-K_first,j_offset:j_offset+L_last+1-L_first,h] = merra_cut

t_3 = time.time()
print("2) Loop: " + str(t_3-t_2))

gwa_cut = GWA_array[K_first:K_last + 1, L_first:L_last + 1]
# Todo: Calculate the ratio with zero pixel of gwa_cut instead of GWA_array
gwa_cut_energy = np.power(gwa_cut, 3)  # Convert from wind speed to wind energy
merra_cut_energy_weighting = gwa_cut / np.sum(gwa_cut_energy)  # Redistribution
merra_cut_speed_weighting = np.cbrt(merra_cut_energy_weighting)
# E = (merraData ** 3) * value_num_cells / Num_pix  # Contains the percentage of energy for each pixel of
merra_cut_speed_redistributed = np.repeat(merra_cut_speed_weighting[..., None], 8760, axis=2) * merraData  # Expand the array along a new axis and multiply it with vector E (Broadcasting)
# merra_cut_redistributed_speed = np.cbrt(merra_cut_energy)  # Convert back from energy to wind speed
reMerra[i_offset:i_offset + K_last + 1 - K_first, j_offset:j_offset + L_last + 1 - L_first, :8760] = merra_cut_speed_redistributed


t_4 = time.time()
print("2) Matrix multiplication (Numpy) and math: " + str(t_4-t_3))

print('Finished!')
