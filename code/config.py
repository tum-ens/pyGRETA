import os
from pathlib import Path
import numpy as np

def config():
    """
    This function is the main configuration function that calls all the other modules in the code.

    :return: The dictionary param containing all the user preferences, and the dictionary path containing all the paths to inputs and outputs.
    :rtype: tuple of dict
    """
    paths, param = general_settings()
    paths, param = scope_paths_and_parameters(paths, param)
    
    param = computation_parameters(param)
    param = resolution_parameters(param)
    param = weather_data_parameters(param)
    param = file_saving_options(param)
    param = time_series_parameters(param)
    param = landuse_parameters(param)
    param = protected_areas_parameters(param)
    param = pv_parameters(param)
    param = csp_parameters(param)
    param = onshore_wind_parameters(param)
    param = offshore_wind_paramters(param)
    
    paths = weather_input_folder(paths, param)
    paths = global_maps_input_paths(paths)
    paths = output_folders(paths, param)
    paths = weather_output_paths(paths, param)
    paths = local_maps_paths(paths, param)
    paths = irena_paths(paths, param)
    paths = regression_input_paths(paths)
    
    for tech in param["technology"]:
        paths[tech] = {}
        paths = emhires_input_paths(paths, param, tech)
        paths = potential_output_paths(paths, param, tech)
        paths = regional_analysis_output_paths(paths, param, tech)
    return paths, param
    
    
def general_settings():
    """
    This function creates and initializes the dictionaries param and paths. It also creates global variables for the root folder (root),
    the current folder (current_folder), and the system-dependent file separator (fs).

    :param paths: Dictionary including the path to the land use raster for the scope, and to the output path BUFFER.
    :type paths: dict
    :param param: Dictionary including the user-defined buffer (buffer_pixel_amount), the urban type within the land use map (type_urban), and the georeference dictionary.
    :type param: dict
    :return: The empty dictionary paths, and the dictionary param including some general information.
    :rtype: tuple of dict
    """
    # These variables will be initialized here, then read in other modules without modifying them.
    global fs
    global root
    global current_folder
    
    param = {}
    param["author"] = 'Kais Siala' # the name of the person running the script
    param["comment"] = 'Testing JSON'
    
    paths = {}
    fs = os.path.sep
    current_folder = os.path.dirname(os.path.abspath(__file__))
    # For personal Computer:
    # root = str(Path(current_folder).parent.parent.parent) + fs + "Database_KS" + fs
    # For Server Computer:
    root = str(Path(current_folder).parent.parent.parent) + "Database_KS" + fs
    
    return paths, param
    
###########################
#### User preferences #####
###########################
    
def scope_paths_and_parameters(paths, param):
    """
    This function defines the path of the geographic scope of the output (spatial_scope) and of the subregions of interest (subregions).
    It also associates two name tags for them, respectively (region_name) and (subregions_name), which define the names of output folders.
    Both paths should point to shapefiles of polygons or multipolygons.
    
    For (spatial_scope), only the bounding box around all the features matters.
    Example: In case of Europe, whether a shapefile of Europe as one multipolygon, or as a set of multiple features (countries, states, etc.) is used, does not make a difference.
    Potential maps (theoretical and technical) will be later generated for the whole scope of the bounding box.
    
    For (subregions), the shapes of the individual features matter, but not their scope.
    For each individual feature that lies within the scope, you can later generate a summary report and time series.
    The shapefile of (subregions) does not have to have the same bounding box as (spatial_scope).
    In case it is larger, features that lie completely outside the scope will be ignored, whereas those that lie partly inside it will be cropped using the bounding box
    of (spatial_scope). In case it is smaller, all features are used with no modification.
    
    (year) defines the year of the weather data, and (technology) the list of technologies that you are interested in.
    Currently, four technologies are defined: onshore wind ``WindOn``, offshore wind ``WindOff``, photovoltaics ``PV``, concentrated solar power ``CSP``.

    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict
    :return: The updated dictionaries paths and param.
    :rtype: tuple of dict
    """
    # Paths to the shapefiles
    PathTemp = root + "02 Shapefiles for regions" + fs + "User-defined" + fs
    paths["spatial_scope"] = PathTemp + "gadm36_SGP_0.shp"
    paths["subregions"] = PathTemp + "gadm36_SGP_1.shp"
    
    # Name tags for the scope and the subregions
    param["region_name"] = 'Singapore'  # Name tag of the spatial scope
    param["subregions_name"] = 'Singapore_districts' # Name tag of the subregions
    
    # Year
    param["year"] = 2015
    
    # Technologies
    param["technology"] = ['WindOn', 'PV', 'WindOff', 'CSP']  # ['PV', 'CSP', 'WindOn', 'WindOff']
    return paths, param
  
    
def computation_parameters(param):
    """
    This function defines parameters related to the processing.
    Some modules in ``Master.py`` allow parallel processing. The key (nproc), which takes an integer as a value, limits the number of parallel processes.
    (CPU_limit) is a boolean parameter ?????

    :param param: Dictionary including the user preferences.
    :type param: dict
    :return: The updated dictionary param.
    :rtype: dict
    """
    param["nproc"] = 36
    param["CPU_limit"] = True
    return param
    
def resolution_parameters(param):
    """
    This function defines the resolution of weather data (low resolution), and the desired resolution of output rasters (high resolution).
    Both are numpy array with two numbers. The first number is the resolution in the vertical dimension (in degrees of latitude),
    the second is for the horizontal dimension (in degrees of longitude).

    :param param: Dictionary including the user preferences.
    :type param: dict
    :return: The updated dictionary param.
    :rtype: dict
    """
    
    param["res_weather"] = np.array([1 / 2, 5 / 8])
    param["res_desired"] = np.array([1 / 240, 1 / 240])
    return param
    
def weather_data_parameters(param):
    """
    This function defines the coverage of the weather data (MERRA_coverage), and how outliers should be corrected (MERRA_correction).
    If you have downloaded the MERRA-2 data for the world, enter the name tag 'World'. The code will later search for the data in the corresponding folder.
    It is possible to download the MERRA-2 just for the geographic scope of the analysis. In that case, enter another name tag (we recommend using the same one as the spatial scope).
    
    MERRA-2 contains some outliers, especially in the wind data. (MERRA_correction) sets the threshold of the relative distance between the yearly mean of the data point
    to the yearly mean of its neighbors. 

    :param param: Dictionary including the user preferences.
    :type param: dict
    :return: The updated dictionary param.
    :rtype: dict
    """
    
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

def weather_input_folder(paths, param):
    '''
    '''
    global root
    global fs
    
    MERRA_coverage = param["MERRA_coverage"]
    year = str(param["year"])
    paths["MERRA_IN"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + MERRA_coverage + " " + year + fs
    
    return paths
    
def global_maps_input_paths(paths):
    '''
    '''
    global root
    global fs
    
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
    
    return paths

def output_folders(paths, param):
    '''
    '''
    global root
    global fs
    
    region = param["region_name"]
    subregions = param["subregions_name"]
    
    # Main output folder
    paths["region"] = root + "03 Intermediate files" + fs + "Files " + region + fs
    
    # Output folder for weather data
    paths["weather_data"] = paths["region"] + "Weather data" + fs
    if not os.path.isdir(paths["weather_data"]):
        os.makedirs(paths["weather_data"])
        
    # Output folder for local maps of the scope
    paths["local_maps"] = paths["region"] + "Maps" + fs
    if not os.path.isdir(paths["local_maps"]):
        os.makedirs(paths["local_maps"])
        
    # Output folder for the potential of the scope
    paths["potential"] = paths["region"] + "Renewable energy" + fs + "Potential" + fs
    if not os.path.isdir(paths["potential"]):
        os.makedirs(paths["potential"])
        
    # Output folder for the regional analysis
    paths["regional_analysis"] = paths["region"] + "Renewable energy" + fs + "Regional analysis" + fs + subregions + fs
    if not os.path.isdir(paths["regional_analysis"]):
        os.makedirs(paths["regional_analysis"])
        
    # Regression output
    paths["regression_out"] = paths["regional_analysis"] + "Regression outputs" + fs
    if not os.path.isdir(paths["regression_out"]):
        os.makedirs(paths["regression_out"])

    return paths
    

    
def weather_output_paths(paths, param):
    '''
    '''
    year = str(param["year"])
    paths["W50M"] = paths["weather_data"] + "w50m_" + year + ".mat"
    paths["CLEARNESS"] = paths["weather_data"] + "clearness_" + year + ".mat"
    paths["T2M"] = paths["weather_data"] + "t2m_" + year + ".mat"
    
    return paths
    
def local_maps_paths(paths, param):
    '''
    '''
    # Local maps
    PathTemp = paths["local_maps"] + param["region_name"]
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
    paths["AREA"] = PathTemp + "_Area.mat" # Area per pixel in m²
    
    # Correction factors for wind speeds
    turbine_height_on = str(param["WindOn"]["technical"]["hub_height"])
    turbine_height_off = str(param["WindOff"]["technical"]["hub_height"])
    paths["CORR_ON"] = PathTemp + "_WindOn_Correction_" + turbine_height_on + '.tif'
    paths["CORR_OFF"] = PathTemp + "_WindOff_Correction_" + turbine_height_off + '.tif'
    
    return paths
    
def irena_paths(paths, param):
    '''
    '''
    global root
    global fs
    
    year = str(param["year"])
    # IRENA input
    paths["IRENA"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "IRENA" + fs + "IRENA_RE_electricity_statistics_allcountries_alltech_" + year + ".csv"
    paths["IRENA_dict"] = root + "00 Assumptions" + fs + "dict_IRENA_countries.csv"
    
    # IRENA output
    paths["IRENA_out"] = paths["region"] + "Renewable energy" + fs + "IRENA_summary_" + year + ".csv"
    
    return paths
    
def regression_input_paths(paths):
    '''
    '''
    global current_folder
    global fs
    
    paths["Reg_RM"] = current_folder + fs + "Regression_coef" + fs + "README.txt"
    paths["inst-cap_example"] = current_folder + fs + "Regression_coef" + fs + "IRENA_FLH_example.csv"
    paths["regression_in"] = paths["regional_analysis"]
    
    return paths
    
def emhires_input_paths(paths, param, tech):
    '''
    '''
    global root
    global fs
    
    year = str(param["year"])
    
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
    
def potential_output_paths(paths, param, tech):
    '''
    '''
    
    region = param["region_name"]
    year = str(param["year"])
    
    if tech in ['WindOn', 'WindOff']:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["potential"] + region + '_' + tech + '_' + hubheight
    elif tech in ['PV']:
        if 'orientation' in param["PV"]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = '0'
        PathTemp = paths["potential"] + region + '_' + tech + '_' + orientation
    else:
        PathTemp = paths["potential"] + region + '_' + tech
    
    paths[tech]["FLH"] = PathTemp + '_FLH_' + year + '.mat'
    paths[tech]["mask"] = PathTemp + "_mask_" + year + ".mat"
    paths[tech]["FLH_mask"] = PathTemp + "_FLH_mask_" + year + ".mat"
    paths[tech]["weight"] = PathTemp + "_weight_" + year + ".mat"
    paths[tech]["FLH_weight"] = PathTemp + "_FLH_weight_" + year + ".mat"

    return paths
    
def regional_analysis_output_paths(paths, param, tech):
    '''
    '''
    
    subregions = param["subregions_name"]
    year = str(param["year"])
    
    if tech in ['WindOn', 'WindOff']:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["regional_analysis"] + subregions + '_' + tech + '_' + hubheight
    elif tech in ['PV']:
        if 'orientation' in param["PV"]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = '0'
        PathTemp = paths["regional_analysis"] + subregions + '_' + tech + '_' + orientation
    else:
        PathTemp = paths["regional_analysis"] + subregions + '_' + tech
    
    
    paths[tech]["FLH"] = PathTemp + '_FLH_' + year + '.mat'
    paths[tech]["mask"] = PathTemp + "_mask_" + year + ".mat"
    paths[tech]["FLH_mask"] = PathTemp + "_FLH_mask_" + year + ".mat"
    paths[tech]["weight"] = PathTemp + "_weight_" + year + ".mat"
    paths[tech]["FLH_weight"] = PathTemp + "_FLH_weight_" + year + ".mat"
    paths[tech]["Locations"] = PathTemp + '_Locations.shp'
    paths[tech]["TS"] = PathTemp + '_TS_' + year + '.csv'
    paths[tech]["Region_Stats"] = PathTemp + '_Region_stats_' + year + '.csv'
    paths[tech]["Sorted_FLH"] = PathTemp + '_sorted_FLH_sampled_' + year + '.mat'
    paths[tech]["TS_param"] = paths["regression_in"] + subregions + '_' + tech
    paths[tech]["Regression_summary"] = paths["regression_out"] + subregions + '_' + tech + '_reg_coefficients_'
    paths[tech]["Regression_TS"] = paths["regression_out"] + subregions + '_' + tech + '_reg_TimeSeries_'
    return paths