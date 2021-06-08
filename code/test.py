import numpy as np
# xmin = 3
# xmax = 8
# ymin = 1
# ymax = 6
#
# x_gwa = np.random.rand(10,10)*10
# y_gwa = np.random.rand(10,10)*10
#
# GWA_array = np.ones((10, 10))*5
#
# gwa_value = np.zeros(GWA_array.shape)
#
# index = (xmin <= x_gwa) & (x_gwa < xmax) & (ymin <= y_gwa) & (y_gwa < ymax)
# print(index)
#
# a, b = index.nonzero()
#
# nonzero = index.nonzero()
# nonzero_x = nonzero[0]
# nonzero_y = nonzero[1]
#
# K_first = nonzero_x[0]
# L_first = nonzero_y[0]
# K_last = nonzero_x[-1]
# L_last = nonzero_y[-1]
#
# index_xy = np.unravel_index(index, GWA_array.shape)
# print(index_xy)
#
# gwa_value[index] = GWA_array[index]
# total_num_cells = len(index)
# print(gwa_value)
# print(total_num_cells)


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
GWA_array = np.random.rand(3999, 8052)*10

#gwa_cut = np.zeros([200, 250, 8760])
gwa_cut = GWA_array[K_first:K_last+1, L_first:L_last+1]

E = (merraData ** 3) * value_num_cells / Num_pix # What? E[h]
merra_cut_energy = gwa_cut / np.sum(gwa_cut)    #
merra_cut_energy_redistributed = np.repeat(merra_cut_energy[..., None], 8760, axis=2) * E   # Expand the array along a new axis and multiply it to vector E (Broadcasting)
merra_cut_redistributed = np.cbrt(merra_cut_energy_redistributed * Num_pix)     # Convert back from energy to wind speed
reMerra[i_offset:i_offset + K_last + 1 - K_first, j_offset:j_offset + L_last + 1 - L_first, :] = merra_cut_redistributed

print('Finished!')
