import datetime
import numpy as np
import os
from pathlib import Path

param = {}
paths = {}
fs = os.path.sep
git_RT_folder = os.path.dirname(os.path.abspath(__file__))
# For personal Computer:
# root = str(Path(git_RT_folder).parent.parent) + fs + "Database_KS" + fs
# For Server Computer:
root = str(Path(git_RT_folder).parent.parent) + "Database_KS" + fs

###########################
#### User preferences #####
###########################

# Custom timestamp
timestamp = str(datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
timestamp = '4NEMO'

# Important settings
param["region"] = 'Europe_wo_balkans'
param["year"] = 2015
param["technology"] = ['WindOff']  # ['PV', 'CSP', 'WindOn', 'WindOff']

# Computation
param["nproc"] = 18
param["CPU_limit"] = True

# Data resolution
param["res_weather"] = np.array([1 / 2, 5 / 8])
param["res_desired"] = np.array([1 / 240, 1 / 240])

# Weather
param["MERRA_coverage"] = 'World'

# Mask / Weight
param["savetiff"] = 1  # Save geotiff files of mask and weight rasters

# Reporting
param["report_sampling"] = 100

# Time series
param["quantiles"] = np.array([100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0])

# Regression
regression = {
    "solver": 'gurobi',
    "hub_heights": [80],
    "orientations": []}
param["regression"] = regression

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
               [0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208,
                0.0208, 0.0208, 0.0208, 0.0208]),
           "albedo": np.array(
               [0.00, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.00, 0.20, 0.20, 0.20, 0.00, 0.20]),
           "hellmann": np.array(
               [0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.20, 0.20, 0.25, 0.25, 0.15, 0.15, 0.20, 0.40, 0.20, 0.15, 0.15]),
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
pv["resource"] = {"clearness_correction": 1
                  }
pv["technical"] = {"T_r": 25,
                   "loss_coeff": 0.37,
                   "tracking": 0,
                   "orientation": -90  # | 0: South | 90: West | 180: North | -90: East |
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
                "lu_availability": np.array(
                    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.02, 0.00, 0.02, 0.02, 0.02, 0.00,
                     0.02]),
                "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
                "power_density": 0.000160,
                "f_performance": 0.75
                }
del GCR
if pv["technical"]["tracking"] != 0 and pv["technical"]["orientation"] not in [0, 180]:
    print('WARNING: ' + str(pv["technical"]["tracking"]) + ' axis tracking, overwrites orientation input: '
          + str(pv["technical"]["orientation"]))
    pv["technical"]["orientation"] = 'track_' + str(pv["technical"]["tracking"])
    
# Parameters related to CSP
csp = {}
csp["technical"] = {"T_avg_HTF": 350,  # °C
                    "loss_coeff": 1.06,  # Heat Loss coefficient W/m².K Independent on Wind Speed
                    "loss_coeff_wind": 1.19,  # Multiplied with (Wind speed)^(0.6)
                    "Flow_coeff": 0.95,  # heat transfer to the HTF factor (Flow or heat removal factor)
                    "AbRe_ratio": 0.00079,  # Receiver Area / Concentrator aperture (90mm diameter/8m aperture)
                    "Wind_cutoff": 22  # (m/s)  Maximum Wind Speed for effective tracking (Source:NREL)
                    }
csp["mask"] = {"slope": 20,
               "lu_suitability": np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1]),
               "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
               }
csp["weight"] = {"lu_availability": np.array(
    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.02, 0.00, 0.02, 0.02, 0.02, 0.00, 0.02]),
                 "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
                 "power_density": 0.000160,
                 "f_performance": 0.9 * 0.75
                 }

# Parameters related to onshore wind
windon = {}
windon["resource"] = {"res_correction": 1,
                      "topo_correction": 1,
                      "topo_weight": 'capacity'  # 'none' or 'size' or 'capacity'
                      }
windon["technical"] = {"w_in": 4,
                       "w_r": 15,
                       "w_off": 25,
                       "P_r": 3,
                       "hub_height": 140
                       }
windon["mask"] = {"slope": 20,
                  "lu_suitability": np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1]),
                  "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
                  "buffer_pixel_amount": 1
                  }
windon["weight"] = {"lu_availability": np.array(
    [0.00, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.10, 0.10, 0.10, 0.10, 0.00, 0.10, 0.00, 0.10, 0.00, 0.10]),
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
                        "hub_height": 120
                        }
windoff["mask"] = {"depth": -40,
                   "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
                   }
windoff["weight"] = {"lu_availability": np.array(
    [0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]),
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

region = param["region"]
MERRA_coverage = param["MERRA_coverage"]
year = str(param["year"])

# Regionalization
paths["region"] = root + "03 Intermediate files" + fs + "Files " + region + fs

# MERRA2
paths["MERRA_IN"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + MERRA_coverage + " " + year + fs
paths["weather_dat"] = paths["region"] + "Renewable energy" + fs + "MERRA2 " + year + fs
paths["U50M"] = paths["weather_dat"] + "u50m_" + year + ".mat"
paths["V50M"] = paths["weather_dat"] + "v50m_" + year + ".mat"
paths["W50M"] = paths["weather_dat"] + "w50m_" + year + ".mat"
paths["GHI"] = paths["weather_dat"] + "swgdn_" + year + ".mat"
paths["GHI_net"] = paths["weather_dat"] + "swgnt_" + year + ".mat"
paths["TOA"] = paths["weather_dat"] + "swtdn_" + year + ".mat"
paths["TOA_net"] = paths["weather_dat"] + "swtnt_" + year + ".mat"
paths["CLEARNESS"] = paths["weather_dat"] + "clearness_" + year + ".mat"
paths["CLEARNESS_net"] = paths["weather_dat"] + "clearness_net_" + year + ".mat"
paths["T2M"] = paths["weather_dat"] + "t2m_" + year + ".mat"

# IRENA
paths[
    "inst-cap"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "IRENA " + year + fs + "inst_cap_" + year + ".csv"
paths[
    "IRENA_FLH"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "IRENA " + year + fs + "IRENA_FLH_" + year + ".csv"

# Regression input
paths["Reg_RM"] = git_RT_folder + fs + "Regression_coef" + fs + "README.txt"
paths["inst-cap_example"] = git_RT_folder + fs + "Regression_coef" + fs + "IRENA_FLH_example.csv"

# Global maps
PathTemp = root + "01 Raw inputs" + fs + "Maps" + fs
paths["LU_global"] = PathTemp + "Landuse" + fs + "LCType.tif"
paths["Topo_tiles"] = PathTemp + "Topography" + fs
paths["Pop_tiles"] = PathTemp + "Population" + fs
paths["Bathym_global"] = PathTemp + "Bathymetry" + fs + "ETOPO1_Ice_c_geotiff.tif"
paths["Protected"] = PathTemp + "Protected Areas" + fs + "WDPA_Nov2018-shapefile-polygons.shp"
paths["GWA"] = PathTemp + "Global Wind Atlas" + fs + fs + "windSpeed.csv"

# Shapefiles
PathTemp = root + "02 Shapefiles for regions" + fs + "User-defined" + fs
paths["Countries"] = PathTemp + "Europe_NUTS0_wo_Balkans_with_EEZ.shp"
paths["SHP"] = paths["Countries"]

# Local maps
PathTemp = paths["region"] + "Maps" + fs + param["region"]
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

# Correction factors for wind speeds
turbine_height_on = str(param["WindOn"]["technical"]["hub_height"])
turbine_height_off = str(param["WindOff"]["technical"]["hub_height"])
paths["CORR"] = PathTemp + "_Wind_Correction_" + turbine_height_on + '_' + turbine_height_off + '.tif'

# Ouput Folders
paths["OUT"] = root + "03 Intermediate files" + fs + "Files " + region + fs + "Renewable energy" + fs + timestamp + fs

# Regression folders
paths["regression_in"] = paths["OUT"]
paths["regression_out"] = paths["OUT"] + "Regression_Outputs" + fs

for tech in param["technology"]:
    paths[tech] = {}
    if tech == 'WindOn':
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES " + year + fs + \
                                 "TS.CF.COUNTRY.30yr.date.txt"
    elif tech == 'WindOff':
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES " + year + fs + \
                                 "TS.CF.OFFSHORE.30yr.date.txt"
    elif tech == 'PV':
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES " + year + fs + \
                                 "EMHIRESPV_TSh_CF_Country_19862015.txt"

    if tech in ['WindOn', 'WindOff']:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["OUT"] + region + '_' + tech + '_' + hubheight
    elif tech in ['PV']:
        if 'orientation' in param["PV"]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = '0'
        PathTemp = paths["OUT"] + region + '_' + tech + '_' + orientation
    else:
        PathTemp = paths["OUT"] + region + '_' + tech

    paths[tech]["area"] = paths["OUT"] + region + "_area_" + year + ".mat"
    paths[tech]["FLH"] = PathTemp + '_FLH_' + year + '.mat'
    paths[tech]["mask"] = PathTemp + "_mask_" + year + ".mat"
    paths[tech]["FLH_mask"] = PathTemp + "_FLH_mask_" + year + ".mat"
    paths[tech]["weight"] = PathTemp + "_weight_" + year + ".mat"
    paths[tech]["FLH_weight"] = PathTemp + "_FLH_weight_" + year + ".mat"
    paths[tech]["Locations"] = PathTemp + '_Locations.shp'
    paths[tech]["TS"] = PathTemp + '_TS_' + year + '.csv'
    paths[tech]["Region_Stats"] = PathTemp + '_Region_stats_' + year + '.csv'
    paths[tech]["Sorted_FLH"] = PathTemp + '_sorted_FLH_sampled_' + year + '.mat'
    paths[tech]["TS_param"] = paths["regression_in"] + region + '_' + tech
    paths[tech]["Regression_summary"] = paths["regression_out"] + region + '_' + tech + '_reg_coefficients_'
    paths[tech]["Regression_TS"] = paths["regression_out"] + region + '_' + tech + '_reg_TimeSeries_'


        
del root, PathTemp, fs
