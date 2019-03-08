import datetime
import numpy as np
import os


# User preferences
# np.__config__.show()
region = 'Germany'
year = '2015'
technology = 'PV'  # 'PV','Wind','CSP'
windtechnology = ['Offshore', 'Onshore']  # 'Offshore', 'Onshore'

quantiles = np.array([100, 97, 95, 90, 75, 67, 50, 30])
correction = 1
savetiff = 1  # Save geotiff files of mask and weight rasters
# Correction factors
a = 2
b = 0.1

# Buffer amount around urban areas
buffer_pixel_amount = 0

# MERRA_Centroid_Extent = [74.5, 45, 19, -20.625]  # EUMENA
# MERRA_Centroid_Extent = [74.5, 36.25, 33.5, -16.25]  # Europe
# MERRA_Centroid_Extent = [49, -103.75, 28, -129.375]  # California
# MERRA_Centroid_Extent = np.array([56.25, 15.3125, 47.25, 2.8125])  # Germany

res = np.array([[1 / 2, 5 / 8], [1 / 240, 1 / 240]])
landuse = {"type": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])}
cost = {}
description = {}
if technology == 'PV':
    # Landuse reclassification
    landuse["Ross_coeff"] = np.array([208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208, 208,
                                      208])
    landuse["albedo"] = np.array([0, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 0, 25, 25, 25, 0, 25])
    landuse["suit_s"] = np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1])
    landuse["avail_s"] = np.array([0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 2, 2, 2, 0, 2])
    landuse["cost_s"] = np.array([100, 50, 50, 50, 50, 50, 20, 20, 20, 20, 10, 80, 10, 15, 10, 100, 80])
    # Technology-related parameters Yingli PV module
    pv = {"T_r": 25,
          "loss_coeff": 0.37,
          "tracking": 0}
    # For the GCR calculation
    GCR = {"shadefree_period": 6,
           "day_north": 79,
           "day_south": 263}
    # Mask
    mask = {"slope_s": 20,
            "suit_s": landuse["suit_s"]}
    # weight
    weight = {"GCR": GCR, "avail_s": landuse["avail_s"],
              "f_pd_pv": 0.000160,
              "f_performance_pv": 0.75,
              "f_pd_csp": 0.00016,
              "f_performance_csp": 0.9 * 0.75}
    # Cost
    cost = {"lu_s": landuse["cost_s"],
            "unit_pv": 1 * 10 ** 7,
            "unit_csp": 1 * 10 ** 7}
    description = {"pv": pv}

if technology == 'Wind':
    # Landuse reclassification
    landuse["hellmann"] = np.array([10, 25, 25, 25, 25, 25, 20, 20, 25, 25, 15, 15, 20, 40, 20, 15, 15])
    landuse["height"] = np.array([213, 366, 366, 366, 366, 366, 320, 320, 366, 366, 274, 274, 320, 457, 320, 274, 274])

    # Technology-related parameters
    turbine = {}
    weight = {}
    if 'Onshore' in windtechnology:
        # Onshore specific parameters
        landuse["Onshore"] = {"suit_w": np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1]),
                              "avail_w": np.array([0, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 0, 10, 0, 10, 0, 10]),
                              "cost_w": np.array(
                                  [100, 20, 20, 20, 20, 20, 20, 20, 20, 20, 10, 80, 10, 100, 10, 100, 10])
                              }

        turbine["Onshore"] = {"w_in": 4,
                              "w_r": 15,
                              "w_off": 25,
                              "P_r": 3,
                              "hub_height": 80}
        # Weight
        weight["Onshore"] = {"avail_w": landuse["Onshore"]["avail_w"],
                             "f_pd_w": 0.000008,
                             "f_performance_w": 0.87}

    if 'Offshore' in windtechnology:
        # Offshore specific parameters
        landuse["Offshore"] = {"suit_w": np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                                "avail_w": np.array([10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                                "cost_w": np.array(
                                    [100, 20, 20, 20, 20, 20, 20, 20, 20, 20, 10, 80, 10, 100, 10, 100, 10])
                                }

        turbine["Offshore"] = {"w_in": 3,
                               "w_r": 16.5,
                               "w_off": 34,
                               "P_r": 7.58,
                               "hub_height": 135}
        # Weight
        weight["Offshore"] = {"avail_w": landuse["Offshore"]["avail_w"],
                              "f_pd_w": 0.000020,
                              "f_performance_w": 0.87}
    # Mask
    mask = {"slope_w": 20,
            "depth": -40
            }
    # costs

    cost = {"depth_thres": 40,
            "unit_w": 1.37 * 10 ** 7,
            "offshore": 25}


    def depth(x):
        return np.exp((float(x) - cost["depth_thres"]) / 5)


    cost["depth"] = depth

    # Output description
    description = {"turbine": turbine}
# No distinction between technologies
cost["protected"] = 100
cost["slope_thres"] = 3


def slope(x):
    return np.exp((x - cost["slope_thres"]) / 10)


cost["slope"] = slope
description["region"] = region
description["landuse"] = landuse

# Protected Areas

protected_areas = {"pa_type": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
                   "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
                   "pa_availability": np.array([1, 0, 0, 0, 0, 0, 0.25, 1, 1, 1, 1]),
                   "IUCN_Category": np.array(['Not Protected',  # 0
                                              'Ia',  # 1
                                              'Ib',  # 2
                                              'II',  # 3
                                              'III',  # 4
                                              'IV',  # 5
                                              'V',  # 6
                                              'VI',  # 7
                                              'Not Applicable',  # 8
                                              'Not Assigned',  # 9
                                              'Not Reported'  # 10
                                              ])}

# Initialization - Define Paths
fs = os.path.sep
root = os.path.dirname(os.path.abspath(__file__)) + fs
paths = {}

# Shapefiles
PathTemp = root + "INPUTS" + fs + region + fs + "Shapefile" + fs + region
paths["SHP"] = PathTemp + "_NUTS0_wo_Balkans_with_EEZ.shp"

# MERRA2
PathTemp = root + "INPUTS" + fs + region + fs + "MERRA2 " + year + fs

paths["U50M"] = PathTemp + "u50m_" + year + ".mat"
paths["V50M"] = PathTemp + "v50m_" + year + ".mat"
paths["GHI"] = PathTemp + "swgdn_" + year + ".mat"
paths["TOA"] = PathTemp + "swtdn_" + year + ".mat"
paths["T2M"] = PathTemp + "t2m_" + year + ".mat"  # Temperature

# Global maps
PathTemp = root + "INPUTS" + fs + "Global maps" + fs

paths["LU_global"] = PathTemp + "Landuse" + fs + "LCType.tif"
paths["Topo_tiles"] = PathTemp + "Topography" + fs
paths["Pop_tiles"] = PathTemp + "Population" + fs
paths["Bathym_global"] = PathTemp + fs + "Bathymetry" + "ETOPO1_Ice_c_geotiff.tif"
paths["Protected"] = PathTemp + "Protected Areas" + fs + "WDPA_Jan2018-shapefile-polygons.shp"

# Local mapsTOP
PathTemp = root + "INPUTS" + fs + region + fs + region

paths["LAND"] = PathTemp + "_Land.tif"  # Land pixels
paths["EEZ"] = PathTemp + "_EEZ.tif"  # Sea pixels
paths["LU"] = PathTemp + "_Landuse.tif"  # Land use types
paths["TOPO"] = PathTemp + "_Topography.tif"  # Topography
paths["PA"] = PathTemp + "_Protected_areas.tif"  # Protected areas
paths["SLOPE"] = PathTemp + "_Slope.tif"  # Slope
paths["BATH"] = PathTemp + "_Bathymetry.tif"  # Bathymetry
paths["POP"] = PathTemp + "_Population.tif"  # Population
paths["PopBuff"] = PathTemp + "_Population_Buffered.tif"  # Buffered population

del PathTemp
# Define Output Folder Function

# Ouput Folders
timestamp = str(datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
paths["OUT"] = root + "OUTPUTS" + fs + region + fs + timestamp + fs
if technology == "Wind":
    paths["OUT"] = root + "OUTPUT" + fs + region + fs + str(turbine["hub_height"]) + "m_" + str(correction) + "corr_" + timestamp
else:
    paths["OUT"] = root + "OUTPUT" + fs + region + fs + str(pv["tracking"]) + "axis_" + timestamp)
paths["mask"] = paths["OUT"] + region + "_" + technology + "_mask_" + year + ".mat"
paths["FLH_mask"] = paths["OUT"] + region + "_" + technology + "_FLH_mask_" + year + ".mat"
paths["area"] = paths["OUT"] + region + "_" + technology + "_area_" + year + ".mat"
paths["weight"] = paths["OUT"] + region + "_" + technology + "_weight_" + year + ".mat"
paths["FLH_weight"] = paths["OUT"] + region + "_" + technology + "_FLH_weight_" + year + ".mat"

description["paths"] = paths
del root
