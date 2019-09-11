import datetime
import numpy as np
import os

###########################
#### User preferences #####
###########################

param = {}
param["region"] = 'Germany'
param["year"] = '2015'
param["technology"] = ['WindOn']  # ['PV', 'CSP', 'WindOn', 'WindOff']
param["quantiles"] = np.array([100, 97, 95, 90, 75, 67, 50, 30])
param["savetiff"] = 1  # Save geotiff files of mask and weight rasters
param["nproc"] = 6
param["CPU_limit"] = True
param["report_sampling"] = 100


# Regression Coefficient
param["solver"] = 'gurobi'
param["hub_heights"] = np.array([60, 80, 100])
param["no_solution"] = '**'

# MERRA_Centroid_Extent = [74.5, 45, 19, -20.625]  # EUMENA
# MERRA_Centroid_Extent = [74.5, 36.25, 33.5, -16.25]  # Europe
# MERRA_Centroid_Extent = [49, -103.75, 28, -129.375]  # California
# MERRA_Centroid_Extent = np.array([56.25, 15.3125, 47.25, 2.8125])  # Germany

param["res_weather"] = np.array([ 1/2,   5/8])
param["res_desired"] = np.array([1/240, 1/240])

# Landuse reclassification
# A_lu matrix element values range from 0 to 16:
# 0   -- Water
# 1   -- Evergreen needle leaf forest
# 2   -- Evergreen broad leaf forest
# 3   -- Deciduous needle leaf forest
# 4   -- deciduous broad leaf forest
# 5   -- Mixed forests
# 6   -- Closed shrublands
# 7   -- Open shrublands
# 8   -- Woody savannas
# 9   -- Savannas
# 10  -- Grasslands
# 11  -- Permanent wetland
# 12  -- Croplands
# 13  -- URBAN AND BUILT-UP
# 14  -- Croplands / natural vegetation mosaic
# 15  -- Snow and ice
# 16  -- Barren or sparsely vegetated

landuse = {"type": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]),
           "type_urban": 13,
           "Ross_coeff": np.array(
               [0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208]),
           "albedo": np.array([0.00, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.00, 0.20, 0.20, 0.20, 0.00, 0.20]),
           "hellmann": np.array([0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.20, 0.20, 0.25, 0.25, 0.15, 0.15, 0.20, 0.40, 0.20, 0.15, 0.15]),
           "height": np.array([213, 366, 366, 366, 366, 366, 320, 320, 366, 366, 274, 274, 320, 457, 320, 274, 274])
           }

# Protected Areas
protected_areas = {"type": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
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
                                              ])
                   }
param["landuse"] = landuse
param["protected_areas"] = protected_areas
del landuse, protected_areas

# Parameters related to PV
pv = {}
pv["resource"] = {"clearness_correction": 0.85
                  }
pv["technical"] = {"T_r": 25,
                   "loss_coeff": 0.37,
                   "tracking": 0
                   }
pv["mask"] = {"slope": 20,
              "lu_suitability": np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1]),
              "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
              }
GCR = {"shadefree_period": 6,
       "day_north": 79,
       "day_south": 263
       }
pv["weight"] = {"GCR": GCR,
                "lu_availability": np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.02, 0.00, 0.02, 0.02, 0.02, 0.00, 0.02]),
                "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
                "power_density": 0.000160,
                "f_performance": 0.75
                }
del GCR

# Parameters related to CSP
csp = {}
csp["mask"] = {"slope": 20,
               "lu_suitability": np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1]),
               "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
               }
csp["weight"] = {"lu_availability": np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.02, 0.00, 0.02, 0.02, 0.02, 0.00, 0.02]),
                 "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
                 "power_density": 0.000160,
                 "f_performance": 0.9 * 0.75
                 }

# Parameters related to onshore wind
windon = {}
windon["resource"] = {"res_correction": 1,
                      "topo_correction": 1,
                      "topo_weight": 'capacity' # 'none' or 'size' or 'capacity'
                      }
windon["technical"] = {"w_in": 4,
                       "w_r": 15,
                       "w_off": 25,
                       "P_r": 3,
                       "hub_height": 100
                       }
windon["mask"] = {"slope": 20,
                  "lu_suitability": np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1]),
                  "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
                  "buffer_pixel_amount": 1
                  }
windon["weight"] = {"lu_availability": np.array([0.00, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.10, 0.10, 0.10, 0.10, 0.00, 0.10, 0.00, 0.10, 0.00, 0.10]),
                    "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
                    "power_density": 0.000008,
                    "f_performance": 0.87
                    }

# Parameters related to offshore wind
windoff = {}
windoff["resource"] = {"res_correction": 1,
                       }
windoff["technical"] = {"w_in": 3,
                        "w_r": 16.5,
                        "w_off": 34,
                        "P_r": 7.58,
                        "hub_height": 135
                        }
windoff["mask"] = {"depth": -40,
                   "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
                   }
windoff["weight"] = {"lu_availability": np.array([0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]),
                     "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
                     "power_density": 0.000020,
                     "f_performance": 0.87
                     }

param["PV"] = pv
param["CSP"] = csp
param["WindOn"] = windon
param["WindOff"] = windoff
del pv, csp, windon, windoff

###########################
##### Define Paths ########
###########################

fs = os.path.sep
root = os.path.dirname(os.path.abspath(__file__)) + fs + ".." + fs
region = param["region"]
year = str(param["year"])

paths = {}

# Shapefiles
PathTemp = root + "01 Raw inputs" + fs + "Maps" + fs + "Shapefiles" + fs
paths["SHP"] = PathTemp + "Germany_with_EEZ.shp"
paths["Countries"] = PathTemp + region + "_with_EEZ.shp" # for eventual correction with the Global Wind Atlas

# MERRA2
paths["MERRA_IN"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + region + " " + year + fs
PathTemp = root + "02 Intermediate files" + fs + "Files " + region + fs + "Renewable energy" + fs + "MERRA2 " + year + fs
paths["U50M"] = PathTemp + "u50m_" + year + ".mat"
paths["V50M"] = PathTemp + "v50m_" + year + ".mat"
paths["W50M"] = PathTemp + "w50m_" + year + ".mat"
paths["GHI"] = PathTemp + "swgdn_" + year + ".mat"
paths["TOA"] = PathTemp + "swtdn_" + year + ".mat"
paths["CLEARNESS"] = PathTemp + "clearness_" + year + ".mat"
paths["T2M"] = PathTemp + "t2m_" + year + ".mat"

# IRENA
paths["inst-cap"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "IRENA " + year + fs + "inst_cap_" + year + ".csv"

# Global maps
PathTemp = root + "01 Raw inputs" + fs + "Maps" + fs + "Global maps" + fs
paths["LU_global"] = PathTemp + "Landuse" + fs + "LCType.tif"
paths["Topo_tiles"] = PathTemp + "Topography" + fs
paths["Pop_tiles"] = PathTemp + "Population" + fs
paths["Bathym_global"] = PathTemp + "Bathymetry" + fs + "ETOPO1_Ice_c_geotiff.tif"
paths["Protected"] = PathTemp + "Protected Areas" + fs + "WDPA_Nov2018-shapefile-polygons.shp"
paths["GWA"] = PathTemp + "Global Wind Atlas" + fs + fs + "windSpeed.csv"

# Local maps
PathTemp = root + "02 Intermediate files" + fs + "Files " + region + fs + "Maps" + fs + region
paths["LAND"] = PathTemp + "_Land.tif"  # Land pixels
paths["EEZ"] = PathTemp + "_EEZ.tif"  # Sea pixels
paths["LU"] = PathTemp + "_Landuse.tif"  # Land use types
paths["TOPO"] = PathTemp + "_Topography.tif"  # Topography
paths["PA"] = PathTemp + "_Protected_areas.tif"  # Protected areas
paths["SLOPE"] = PathTemp + "_Slope.tif"  # Slope
paths["BATH"] = PathTemp + "_Bathymetry.tif"  # Bathymetry
paths["POP"] = PathTemp + "_Population.tif"  # Population
paths["BUFFER"] = PathTemp + "_Population_Buffered.tif"  # Buffered population
paths["CORR_GWA"] = PathTemp + "_GWA_Correction.mat"  # Correction factors based on the GWA
paths["CORR"] = PathTemp + "_Wind_Correction.tif"  # Correction factors for wind speeds

# Ouput Folders
timestamp = str(datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
timestamp = timestamp
# timestamp = "test"
paths["OUT"] = root + "OUTPUTS" + fs + region + fs + timestamp + fs
if not os.path.isdir(paths["OUT"]):
    os.mkdir(paths["OUT"])

# Regression input
paths["IRENA"] = root + "INPUTS" + fs + region + fs + "EMHIRES_IRENA" + fs + "IRENA_FLH.txt"
paths["Reg_RM"] = root + "INPUTS" + fs + region + fs + "EMHIRES_IRENA" + fs + "README.txt"


# Regression folders
paths["regression"] = root + "OUTPUTS" + fs + region + fs + "Regression" + fs
paths["regression_in"] = paths["regression"] + "INPUTS" + fs
paths["regression_out"] = paths["regression"] + "OUTPUTS" + fs

for tech in param["technology"]:
    paths[tech] = {}

    if tech in ['WindOn', 'WindOff']:
        hubheight = str(param[tech]["technical"]["hub_height"])
        paths[tech]["FLH"] = paths["OUT"] + region + '_' + tech + '_' + hubheight + '_FLH_' + year + '.mat'
        paths[tech]["mask"] = paths["OUT"] + region + "_" + tech + '_' + hubheight + "_mask_" + year + ".mat"
        paths[tech]["FLH_mask"] = paths["OUT"] + region + "_" + tech + '_' + hubheight + "_FLH_mask_" + year + ".mat"
        paths[tech]["area"] = paths["OUT"] + region + "_" + tech + '_' + hubheight + "_area_" + year + ".mat"
        paths[tech]["weight"] = paths["OUT"] + region + "_" + tech + '_' + hubheight + "_weight_" + year + ".mat"
        paths[tech]["FLH_weight"] = paths["OUT"] + region + "_" + tech + '_' + hubheight + "_FLH_weight_" + year + ".mat"
        paths[tech]["Locations"] = paths["OUT"] + region + "_" + tech + '_' + hubheight + '_Locations.shp'
        paths[tech]["TS"] = paths["OUT"] + region + '_' + tech + '_' + hubheight + '_TS_' + year + '.csv'
        paths[tech]["TS_height"] = paths["regression_in"] + region + '_' + tech + '_'
        paths[tech]["Regression_summary"] = paths["regression_out"] + region + '_' + tech + '_reg_coefficients_' + timestamp + '.csv'
        paths[tech]["Region_Stats"] = paths["OUT"] + region + '_' + tech + '_' + hubheight + '_Region_stats_' + year + '.csv'
        paths[tech]["Sorted_FLH"] = paths["OUT"] + region + '_' + tech + '_' + hubheight + '_sorted_FLH_sampled_' + year + '.mat'
        paths[tech]["EMHIRES"] = root + "INPUTS" + fs + region + fs + "EMHIRES_IRENA" + fs + "EMHIRES_" + tech + '.txt'
    else:
        paths[tech]["FLH"] = paths["OUT"] + region + '_' + tech + '_FLH_' + year + '.mat'
        paths[tech]["mask"] = paths["OUT"] + region + "_" + tech + "_mask_" + year + ".mat"
        paths[tech]["FLH_mask"] = paths["OUT"] + region + "_" + tech + "_FLH_mask_" + year + ".mat"
        paths[tech]["area"] = paths["OUT"] + region + "_" + tech + "_area_" + year + ".mat"
        paths[tech]["weight"] = paths["OUT"] + region + "_" + tech + "_weight_" + year + ".mat"
        paths[tech]["FLH_weight"] = paths["OUT"] + region + "_" + tech + "_FLH_weight_" + year + ".mat"
        paths[tech]["Locations"] = paths["OUT"] + region + "_" + tech + '_Locations.shp'
        paths[tech]["TS"] = paths["OUT"] + region + '_' + tech + '_TS_' + year + '.csv'
        paths[tech]["Region_Stats"] = paths["OUT"] + region + '_' + tech + '_Region_stats_' + year + '.csv'
        paths[tech]["Sorted_FLH"] = paths["OUT"] + region + '_' + tech + '_sorted_FLH_sampled_' + year + '.mat'

del root, PathTemp, fs
