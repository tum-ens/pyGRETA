import os
# os.environ['MKL_NUM_THREADS'] = '8'
import numpy as np



def calc_CF_wind(w50m_h, rasterData, turbine, tech):
    # This function calculates the capacity factor for a given wind speed at 50m
    # Calculate the wind speed a the desired height

    w_new_h = w50m_h * rasterData["A_cf"]

    # retreive turbine parameters
    turbine = turbine[tech]

    # Calculate the capacity factor

    a = turbine["w_in"] ** 3 / (turbine["w_in"] ** 3 - turbine["w_r"] ** 3)
    b = 1 / (turbine["w_r"] ** 3 - turbine["w_in"] ** 3)

    CF = np.zeros(w_new_h.shape)
    # Case 1 : above the cut-in speed and below the rated speed
    idx1 = np.logical_and(turbine["w_in"] < w_new_h, w_new_h < turbine["w_r"])
    CF[idx1] = a + b * w_new_h[idx1] ** 3
    # Case 2 : above the rated wind speed and below the cut_off speed
    idx2 = np.logical_and(turbine["w_r"] <= w_new_h, w_new_h <= turbine["w_off"])
    CF[idx2] = 1
    # Other cases (below cut-in or above cut-off
    CF[np.logical_not(np.logical_or(idx1, idx2))] = 0

    return CF

