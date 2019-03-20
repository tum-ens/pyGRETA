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
	
	
def calc_FLH_wind(hours, args):
    # Decomposing the tuple args
    reg = args[0]
    nRegions = args[1]
    region_name = args[2]
    Ind = args[3]
    Crd = args[4]
    res = args[5]
    merraData = args[6]
    rasterData = args[7]
    m = args[8]
    n = args[9]
    W50M = args[10]
    turbine = args[11]
    tech = args[12]
	
    TS = np.zeros((8760, 1))
    FLH = np.zeros((m[1, reg], n[1, reg]))
	
    for hour in hours:
        # Show progress of the simulation
        print(str(reg) + '/' + str(nRegions - 1) + ' ' + region_name + ' ' + str(hour + 1))
        
        # Load MERRA data, increase its resolution, and fit it to the extent
        w50m_h = resizem(W50M[Ind[0, reg, 2] - 1:Ind[0, reg, 0], Ind[0, reg, 3] - 1:Ind[0, reg, 1], hour],
                         m[1, reg], n[1, reg])
        # Calculate hourly capacity factor
        CF = WindFunc.calc_CF_wind(w50m_h, rasterData, turbine, tech)
        # Aggregates CF to obtain the yearly FLH
        
        CF = CF * rasterData["A_region"]
        CF[np.isnan(CF)] = 0
        FLH = FLH + CF
        # Time series for the mean
        TS[hour] = np.mean(CF[rasterData["A_region"] == 1])
    return FLH, TS

