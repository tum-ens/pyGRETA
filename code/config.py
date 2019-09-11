import datetime
import os
from pathlib import Path
import numpy as np

def config():
    
    param, paths, root, current_folder, timestamp, fs = general_settings()
    param = scope_parameters(param)
    paths = shapefiles_input_paths(paths, fs, root)
    param = computation_parameters(param)
    param = resolution_parameters(param)
    return paths, param
    
    
    
    
    
def general_settings():
    ''' '''
    param = {}
    paths = {}
    fs = os.path.sep
    current_folder = os.path.dirname(os.path.abspath(__file__))
    # For personal Computer:
    # root = str(Path(current_folder).parent.parent.parent) + fs + "Database_KS" + fs
    # For Server Computer:
    root = str(Path(current_folder).parent.parent.parent) + "Database_KS" + fs
    
    # Custom timestamp
    timestamp = str(datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
    timestamp = 'test'
    
    return param, paths, root, current_folder, timestamp, fs
    
###########################
#### User preferences #####
###########################
    
def scope_parameters(param):
    '''
    '''
    param["region"] = 'Europe_wo_Balkans'  #'Europe_wo_Balkans'  # Name of the spatial scope, define path to shapefile below!
    param["year"] = 2015
    param["technology"] = ['PV']  # ['PV', 'CSP', 'WindOn', 'WindOff']
    return param
    
def shapefiles_input_paths(paths, fs, root):
    '''
    '''
    
    PathTemp = root + "02 Shapefiles for regions" + fs + "User-defined" + fs
    paths["spatial_scope"] = PathTemp + "Europe_NUTS0_wo_Balkans_with_EEZ.shp"
    paths["subregions"] = PathTemp + "Europe_NUTS0_wo_Balkans_with_EEZ.shp"
    return paths
    
def computation_parameters(param):
    '''
    '''
    param["nproc"] = 36
    param["CPU_limit"] = True
    return param
    
def resolution_parameters(paths, param, root, current_folder, timestamp, fs):
    '''
    '''
    
    param["res_weather"] = np.array([1 / 2, 5 / 8])
    param["res_desired"] = np.array([1 / 240, 1 / 240])
    return param
    
def weather_data_parameters(param):
    '''
    '''
    
    param["MERRA_coverage"] = 'World'
    param["MERRA_correction"] = 0.35
    return param
    
def file_saving_options(param):
    '''
    '''
    # Mask / Weight
    param["savetiff"] = 1  # Save geotiff files of mask and weight rasters
    
    # Reporting
    param["report_sampling"] = 100
    return param
    
def time_series_parameters(param):
    '''
    '''
    # Quantiles for time series
    param["quantiles"] = np.array([100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0])
    
    # Regression
    regression = {
        "solver": 'gurobi',
        "hub_heights": [80],
        "orientations": []}
    param["regression"] = regression
    return param
    
def landuse_parameters(param):
    '''
    '''
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
    param["landuse"] = landuse
    return param
    
def protected_areas_parameters(param):
    '''
    '''
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
    
    param["protected_areas"] = protected_areas
    return param
    
def pv_parameters(param):
    '''
    '''
    pv = {}
    pv["resource"] = {"clearness_correction": 1
                    }
    pv["technical"] = {"T_r": 25,
                    "loss_coeff": 0.37,
                    "tracking": 0,
                    "orientation": 180  # | 0: South | 90: West | 180: North | -90: East |
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
    param["PV"] = pv
    return param
        
def csp_parameters(param):
    '''
    '''
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
    param["CSP"] = csp
    return param
    
def onshore_wind_parameters(param):
    windon = {}
    windon["resource"] = {"res_correction": 1,
                        "topo_correction": 1,
                        "topo_weight": 'capacity'  # 'none' or 'size' or 'capacity'
                        }
    windon["technical"] = {"w_in": 4,
                        "w_r": 13,
                        "w_off": 25,
                        "P_r": 3,
                        "hub_height": 80
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
    param["WindOn"] = windon
    return param
    
def offshore_wind_paramters(param):
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
    param["WindOff"] = windoff
    return param
    
###########################
##### Define Paths ########
###########################
    
    region = param["region"]
    MERRA_coverage = param["MERRA_coverage"]
    year = str(param["year"])
    
    
    
    # MERRA2
    paths["MERRA_IN"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + MERRA_coverage + " " + year + fs
    paths["weather_dat"] = paths["region"] + "Renewable energy" + fs + "MERRA2 " + year + fs
    paths["W50M"] = paths["weather_dat"] + "w50m_" + year + ".mat"
    paths["CLEARNESS"] = paths["weather_dat"] + "clearness_" + year + ".mat"
    paths["T2M"] = paths["weather_dat"] + "t2m_" + year + ".mat"
    
    # IRENA input
    paths["IRENA"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "IRENA" + fs + "IRENA_RE_electricity_statistics_allcountries_alltech_" + year + ".csv"
    paths["IRENA_dict"] = root + "00 Assumptions" + fs + "dict_IRENA_countries.csv"
    
    # IRENA output
    paths["IRENA_out"] = paths["region"] + "Renewable energy" + fs + "IRENA_summary_" + year + ".csv"
    
    # Regression input
    paths["Reg_RM"] = current_folder + fs + "Regression_coef" + fs + "README.txt"
    paths["inst-cap_example"] = current_folder + fs + "Regression_coef" + fs + "IRENA_FLH_example.csv"
    
    # Global maps
    PathTemp = root + "01 Raw inputs" + fs + "Maps" + fs
    paths["LU_global"] = PathTemp + "Landuse" + fs + "LCType.tif"
    paths["Topo_tiles"] = PathTemp + "Topography" + fs
    paths["Pop_tiles"] = PathTemp + "Population" + fs
    paths["Bathym_global"] = PathTemp + "Bathymetry" + fs + "ETOPO1_Ice_c_geotiff.tif"
    paths["Protected"] = PathTemp + "Protected Areas" + fs + "WDPA_Nov2018-shapefile-polygons.shp"
    paths["GWA"] = PathTemp + "Global Wind Atlas" + fs + fs + "windSpeed.csv"
    paths["Countries"] = PathTemp + "Countries" + fs + "gadm36_0.shp"
    paths["EEZ_global"] = PathTemp + "EEZ" + fs + "eez_v10.shp"
    
    # Local maps
    PathTemp = paths["region"] + "Maps" + fs + param["region"]
    paths["LAND"] = PathTemp + "_Land.tif"  # Land pixels
    paths["EEZ"] = PathTemp + "_EEZ.tif"  # Sea pixels
    paths["SUB"] = PathTemp + "_Subregions.tif"  # Subregions pixels
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
    paths["CORR_ON"] = PathTemp + "_WindOn_Correction_" + turbine_height_on + '.tif'
    paths["CORR_OFF"] = PathTemp + "_WindOff_Correction_" + turbine_height_off + '.tif'
    
    
    
    # Regression folders
    paths["regression_in"] = paths["OUT"]
    paths["regression_out"] = paths["OUT"] + "Regression_Outputs" + fs
    
    for tech in param["technology"]:
        paths[tech] = {}
        paths = emhires_input_paths(paths, tech, root, fs, year)
        paths = technology_specific_output_paths(paths, param, tech, region, year)
    return paths, param
    
def output_folders(paths, root, fs):
    '''
    '''
    # Main output folder
    paths["region"] = root + "03 Intermediate files" + fs + "Files " + region + fs
    if not os.path.isdir(paths["region"]):
        os.makedirs(paths["region"] + "Renewable energy" + fs)
        os.makedirs(paths["region"] + "Maps" + fs)
        
    # Ouput Folders
    paths["OUT"] = root + "03 Intermediate files" + fs + "Files " + region + fs + "Renewable energy" + fs + timestamp + fs
    
def create_folders(paths):
    fs = os.path.sep
    
    if not os.path.isdir(paths["weather_dat"]):
        os.makedirs(paths["weather_dat"])
    if not os.path.isdir(paths["OUT"]):
        os.makedirs(paths["OUT"])

def emhires_input_paths(paths, tech, root, fs, year):
    '''
    '''
    if tech == 'WindOn':
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES " + year + fs + \
                                "TS.CF.COUNTRY.30yr.date.txt"
    elif tech == 'WindOff':
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES " + year + fs + \
                                "TS.CF.OFFSHORE.30yr.date.txt"
    elif tech == 'PV':
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES " + year + fs + \
                                "EMHIRESPV_TSh_CF_Country_19862015.txt"
    return paths
    
def technology_specific_output_paths(paths, param, tech, region, year):
    '''
    '''
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
    return paths