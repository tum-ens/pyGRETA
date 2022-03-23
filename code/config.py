from pathlib import Path
from warnings import warn
import pandas as pd
import numpy as np
import os
import sys

def configuration(config_file):
    """
    This function is the main configuration function that calls all the other modules in the code.

    :return (paths, param): The dictionary paths containing all the paths to inputs and outputs, and the dictionary param containing all the user preferences.
    :rtype: tuple(dict, dict)
    """
    paths, param = general_settings()
    paths, param = scope_paths_and_parameters(paths, param, config_file)

    param = computation_parameters(param)
    param = resolution_parameters(param)
    param = weather_data_parameters(param)
    param = file_saving_options(param)
    param = time_series_parameters(param)
    param = landuse_parameters(param)
    param = protected_areas_parameters(param)
    param = osm_areas(param)
    param = buffers(param)
    param = openfieldpv_parameters(param)
    param = rooftoppv_parameters(param)
    param = csp_parameters(param)
    param = onshore_wind_parameters(param)
    param = offshore_wind_paramters(param)
    param = biomass_parameters(param)

    paths = weather_input_folder(paths, param)
    paths, param = global_maps_input_paths(paths, param)
    paths = output_folders(paths, param)
    paths = weather_output_paths(paths, param)
    paths = local_maps_paths(paths, param)
    paths = irena_paths(paths, param)

    for tech in param["technology"]:
        paths[tech] = {}
        paths = regression_paths(paths, param, tech)
        paths = emhires_input_paths(paths, tech)
        paths = potential_output_paths(paths, param, tech)
        paths = regional_analysis_output_paths(paths, param, tech)
        paths = discrete_output_paths(paths, param, tech)
    return paths, param


def general_settings():
    """
    This function creates and initializes the dictionaries param and paths. It also creates global variables for the root folder ``root``,
    and the system-dependent file separator ``fs``.

    :return (paths, param): The empty dictionary paths, and the dictionary param including some general information.
    :rtype: tuple(dict, dict)
    """
    # These variables will be initialized here, then read in other modules without modifying them.
    global fs
    global root
    fs = os.path.sep

    paths = {}

    param = {}
    param["author"] = "Thushara Addanki"  # the name of the person running the script
    param["comment"] = "Potential Analysis"
    param["path_database_windows"] = "..\..\pyGRETA\Database_KS"    # specify the relative (from the runme.py file) or the absolute path to the database folder
    param["path_database_linux"] = os.path.expanduser("~") + fs + "pygreta" + fs + "0_database"  # specify path to database on linux

    if sys.platform.startswith("win"):
        root = param["path_database_windows"]       # use windows location of database folder
        # print('Working on windows')
    elif sys.platform.startswith("linux"):
        root = param["path_database_linux"]         # use linux location of database folder
        # print('Working on linux')

    if not root.endswith(fs):
        root += fs      # add final slash('\' or '/') if not already there

    return paths, param


###########################
#### User preferences #####
###########################


def scope_paths_and_parameters(paths, param, config_file):
    """
    This function defines the path of the geographic scope of the output *spatial_scope* and of the subregions of interest *subregions*.
    Both paths should point to shapefiles of polygons or multipolygons.
    It also associates two name tags for them, respectively *region_name* and *subregions_name*, which define the names of output folders.
    
    * For *spatial_scope*, only the bounding box around all the features matters.
      Example: In case of Europe, whether a shapefile of Europe as one multipolygon, or as a set of multiple features (countries, states, etc.) is used, does not make a difference.
      Potential maps (theoretical and technical) will be later generated for the whole scope of the bounding box.
    
    * For *subregions*, the shapes of the individual features matter, but not their scope.
      For each individual feature that lies within the scope, you can later generate a summary report and time series.
      The shapefile of *subregions* does not have to have the same bounding box as *spatial_scope*.
      In case it is larger, features that lie completely outside the scope will be ignored, whereas those that lie partly inside it will be cropped using the bounding box
      of *spatial_scope*. In case it is smaller, all features are used with no modification.
    
    * *year* defines the year of the input data. 
    
    * *technology* defines the list of technologies that you are interested in.
      Currently, four technologies are defined: onshore wind ``'WindOn'``, offshore wind ``'WindOff'``, photovoltaics ``'PV'``, concentrated solar power ``'CSP'``.

    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return (paths, param): The updated dictionaries paths and param.
    :rtype: tuple(dict, dict)
    """
    # Paths to the shapefiles
    PathTemp = root + "02 Shapefiles for regions" + fs + "User-defined" + fs

    input_df = pd.read_csv('../configs' + fs + config_file, delimiter=':', comment='#', header=None, index_col=0,
                           skip_blank_lines=True, )  # Import parameters from config_files in folder 'configs'
    input_dict = input_df[1].to_dict()  # Convert dataframe to dict with values from the first column

    param["region_name"] = input_dict["region_name"].replace(" ", "")  # Name tag of the spatial scope
    param["subregions_name"] = input_dict["subregions_name"] .replace(" ", "") # Name tag of the subregions
    param["country_code"] = input_dict["country_code"].replace(" ", "")
    paths["subregions"] = PathTemp + "gadm40_" + param["country_code"] + "_shp" + fs + input_dict["regions"].replace(
        " ", "")
    param["year"] = int(input_dict["year"].replace(" ", ""))  # Convert string 'xxxx' to int
    param["technology"] = input_dict["technology"].replace(" ", "").split(',')  # Creat array by comma separated string
    param["gid"] = input_dict["gid"].replace(" ", "")  # Define spatial level according to GID

    return paths, param

def computation_parameters(param):
    """
    This function defines parameters related to the processing:
    
    * *nproc* is an integer that limits the number of parallel processes (some modules in ``potential.py`` and ``time_series.py`` allow parallel processing).
    
    * *CPU_limit* is a boolean parameter that sets the level of priority for all processes in the multiprocessesing.
      Leave ``True`` if you plan on using the computer while FLH and TS are being computed, ``False`` for fastest computation time.

    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    param["nproc"] = 10
    param["CPU_limit"] = True
    return param


def resolution_parameters(param):
    """
    This function defines the resolution of weather data (low resolution), and the desired resolution of output rasters (high resolution).
    Both are numpy arrays with two numbers. The first number is the resolution in the vertical dimension (in degrees of latitude),
    the second is for the horizontal dimension (in degrees of longitude).

    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    param["res_weather"] = np.array([1 / 2, 5 / 8])
    param["res_desired"] = np.array([1 / 400, 1 / 400])
    return param


def weather_data_parameters(param):
    """
    This function defines the coverage of the weather data *MERRA_coverage*, and how outliers should be corrected using *MERRA_correction*:
    
    * *MERRA_coverage*: If you have downloaded the MERRA-2 data for the world, enter the name tag ``'World'``. The code will later search for the data in the corresponding folder.
      It is possible to download the MERRA-2 just for the geographic scope of the analysis. In that case, enter another name tag (we recommend using the same one as the spatial scope).
    
    * *MERRA_correction*: MERRA-2 contains some outliers, especially in the wind data. *MERRA_correction* sets the threshold of the relative distance between the yearly mean of the data point
      to the yearly mean of its neighbors. 

    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    param["MERRA_coverage"] = "World"
    param["MERRA_correction"] = True
    param["MERRA_correction_factor"] = {"W50M": 0.35, "W50M_offshore": 0.35,
                                        "CLEARNESS": 0.35, "T2M": 0.35}  # Wind Speed  # Clearness index  # Temperature at 2 m
    return param


def file_saving_options(param):
    """
    This function sets some options for saving files.
    
    * *savetiff* is a boolean that determines whether tif rasters for the potentials are saved (``True``), or whether only mat files are saved (``False``).
      The latter are saved in any case.
    
    *  *report_sampling* is an integer that sets the sample size for the sorted FLH values per region (relevant for :mod:`potential.reporting`).
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    # imput maps
    param["savetiff_inputmaps"] = False

    # Mask / Weight
    param["savetiff_potentials"] = True  # Save geotiff files of mask and weight rasters

    # Reporting
    param["report_sampling"] = 100
    return param


def time_series_parameters(param):
    """
    This function determines the time series that will be created.
    
    * *quantiles* is a list of floats between 100 and 0. Within each subregion, the FLH values will be sorted,
      and points with FLH values at a certain quantile will be later selected. The time series will be created for these points.
      The value 100 corresponds to the maximum, 50 to the median, and 0 to the minimum.
      
    * *regression* is a dictionary of options for :mod:`regression.regression_coefficients`:
      
      - *solver* is the name of the solver for the regression.
      - *WindOn* is a dictionary containing a list of hub heights that will be considered in the regression, with a name tag for the list.
      - *WindOff* is a dictionary containing a list of hub heights that will be considered in the regression, with a name tag for the list.
      - *PV* is a dictionary containing a list of orientations that will be considered in the regression, with a name tag for the list.
      - *CSP* is a dictionary containing a list of settings that will be considered in the regression, with a name tag for the list.
    
      If all the available settings should be used, you can leave an empty list.
      
    * *modes* is a dictionary that groups the quantiles and assigns names for each subgroup. You can define the groups as you wish.
      If you want to use all the quantiles in one group without splitting them in subgroups, you can write::
      
       param["modes"] = {"all": param["quantiles"]}
      
    * *combo* is a dictionary of options for :mod:`time_series.generate_stratified_timeseries`:
    
      - *WindOn* is a dictionary containing the different combinations of hub heights for which stratified time series should be generated, with a name tag for each list.
      - *WindOff* is a dictionary containing the different combinations of hub heights for which stratified time series should be generated, with a name tag for each list.
      - *PV* is a dictionary containing the different combinations of orientations for which stratified time series should be generated, with a name tag for each list.
      - *CSP* is a dictionary containing the different combinations of settings for which stratified time series should be generated, with a name tag for each list.
      
      If all the available settings should be used, you can leave an empty list.
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    # Quantiles for time series
    param["quantiles"] = [90,70,50,30]

    # User defined locations
    param["useloc"] = {
        "Augsburg": (48.370, 10.898), "Forchheim": (49.720, 11.059),
        "Kempten": (47.725, 10.315), "München": (48.144, 11.571),
        "Neuburg": (48.736, 11.181), "Nürnberg": (49.453, 11.077),
        "Rosenheim": (47.855, 12.125), "Straubing": (48.882, 12.570),
        "Waldmünchen": (49.378, 12.706)
    }  # {"point name": (latitude, longitude),...}

    # Regression
    param["regression"] = {
        "solver": "gurobi",  # string
        "WindOn": {"all": []},  # dictionary of hub height combinations
        "WindOff": {"80m": []},  # dictionary of hub height combinations
        "OpenFieldPV": {"all": []},  # list of orientation combinations
        "RoofTopPV": {"all": []},  # list of orientation combinations
        "CSP": {"all": []},
    }

    # Stratified time series
    param["modes"] = {"high": [90, 70], "mid": [60, 40], "low": [], "all": param["quantiles"]}
    param["combo"] = {
        # dictionary of hub height and orientation combinations
        "WindOn": {"2015": [60, 80, 100], "2030": [80, 100, 120], "2050": [100, 120, 140]},
        "WindOff": {"80m": [80], "100m": [100], "120m": [120]},
        "OpenFieldPV": {"Solar": [0, 180, -90, 90]},
        "RoofTopPV": {"Solar": [0, 180, -90, 90]},
        "CSP": {"all": []},
    }

    return param


def landuse_parameters(param):
    """
    This function sets the land use parameters in the dictionary *landuse* inside param:
    
      * *type* is a numpy array of integers that associates a number to each land use type.
      * *type_urban* is the number associated to urban areas (useful for :mod:`input_maps.generate_buffered_population`).
      * *Ross_coeff* is a numpy array of Ross coefficients associated to each land use type (relevant for :mod:`physical_models.loss`).
      * *albedo* is a numpy array of albedo coefficients between 0 and 1 associated to each land use type (relevant for reflected irradiation, see :mod:`physical_models.calc_CF_solar`).
      * *hellmann* is a numpy array of Hellmann coefficients associated to each land use type (relevant for :mod:`correction_functions.generate_wind_correction`).

    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict

# Landuse reclassification
        # 0	    No data
        # 10	Cropland, rain-fed
        # 11	Herbaceous cover
        # 12	Tree or shrub cover
        # 20	Cropland, irrigated or post-flooding
        # 30	Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)
        # 40	Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)
        # 50	Tree cover, broadleaved, evergreen, closed to open (>15%)
        # 60	Tree cover, broadleaved, deciduous, closed to open (>15%)
        # 61	Tree cover, broadleaved, deciduous, closed (>40%)
        # 62	Tree cover, broadleaved, deciduous, open (15-40%)
        # 70	Tree cover, needleleaved, evergreen, closed to open (>15%)
        # 71	Tree cover, needleleaved, evergreen, closed (>40%)
        # 72	Tree cover, needleleaved, evergreen, open (15-40%)
        # 80	Tree cover, needleleaved, deciduous, closed to open (>15%)
        # 81	Tree cover, needleleaved, deciduous, closed (>40%)
        # 82	Tree cover, needleleaved, deciduous, open (15-40%)
        # 90	Tree cover, mixed leaf type (broadleaved and needleleaved)
        # 100	Mosaic tree and shrub (>50%) / herbaceous cover (<50%)
        # 110	Mosaic herbaceous cover (>50%) / tree and shrub (<50%)
        # 120	Shrubland
        # 121	Shrubland evergreen
        # 122	Shrubland deciduous
        # 130	Grassland
        # 140	Lichens and mosses
        # 150	Sparse vegetation (tree, shrub, herbaceous cover) (<15%)
        # 151	Sparse tree (<15%)
        # 152	Sparse shrub (<15%)
        # 153	Sparse herbaceous cover (<15%)
        # 160	Tree cover, flooded, fresh or brakish water
        # 170	Tree cover, flooded, saline water
        # 180	Shrub or herbaceous cover, flooded, fresh/saline/brakish water
        # 190	Urban areas
        # 200	Bare areas
        # 201	Consolidated bare areas
        # 202	Unconsolidated bare areas
        # 210	Water bodies
        # 220	Permanent snow and ice


    """
    landuse = {
        "type": np.array([0, 10, 11, 12, 20, 30, 40, 50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100, 110, 120,
                        121, 122, 130, 140, 150, 151, 152, 153, 160, 170, 180, 190, 200, 201, 202, 210, 220]),
        "Ross_coeff": np.array([0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208,
                                0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208,
                                0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208,
                                0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208]),
        "albedo": np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                            0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0]),
        "hellmann": np.array([0.25, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                              0.25, 0.25, 0.25, 0.2, 0.2, 0.2, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.25, 0.25, 0.25,
                              0.4, 0.15, 0.15, 0.15, 0.1, 0.15])
    }
    param["landuse"] = landuse
    return param


def protected_areas_parameters(param):
    """
    This function sets the parameters for protected areas in the dictionary *protected_areas* inside param:
    
      * *type* is a numpy array of integers that associates a number to each protection type.
      * *IUCN_Category* is an array of strings with names associated to each protection type (for your information).
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    protected_areas = {
        "type": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        "IUCN_Category": np.array(
            [
                "Not Protected",  # 0
                "Ia",  # 1
                "Ib",  # 2
                "II",  # 3
                "III",  # 4
                "IV",  # 5
                "V",  # 6
                "VI",  # 7
                "Not Applicable",  # 8
                "Not Assigned",  # 9
                "Not Reported",  # 10
            ]
        ),
    }

    param["protected_areas"] = protected_areas
    return param

def osm_areas(param):

    #"fclass" ILIKE 'commercial' OR 'industrial' OR 'quarry' OR 'military' OR 'park' OR 'recreation_ground'
    osm_areas = {
        "type": np.array([1, 2, 3, 4, 5, 6]),
        "Category": np.array(
            [
                "commercial",  # 1
                "industrial",  # 2
                "quarry",  # 3
                "military",  # 4
                "park",  # 5
                "recreation_ground" #6
            ]
        ),
    }

    param["osm_areas"] = osm_areas
    return param

def buffers(param):

    buffer = {
        "snow": 4, #Landuse
        "water": 1, #Landuse
        "wetland": 1, #Landuse

        "protected_areas_pv": 1, #ProtectedAreas
        "protected_areas_windon": 2, #ProtectedAreas
        "protected_areas_windoff": 4, #ProtectedAreas

        "airport_windon": 16, #Airports
        "boarder": 2, #gadm

        "commercial_windon" : 2, #osm
        "industrial_windon" : 1, #osm
        "mining" : 1, #osm
        "military_windon" : 2, #osm
        "park_pv" : 1, #osm
        "park_windon" : 3, #osm
        "recreation_windon" : 1, #osm

        "settlement_pv": 1, #WSF
        "settlement_windon" : 4, #WSF

        "hydrolakes" : 1,
        "hydrorivers_pv" : 1,
    }

    param["buffer"] = buffer
    return param

def openfieldpv_parameters(param):
    """
    This function sets the parameters for photovoltaics in the dictionary *pv* inside param:
    
    * *resource* is a dictionary including the parameters related to the resource potential:
    
      * *clearness_correction* is a factor that will be multiplied with the clearness index matrix to correct it. If no correction is required, leave it equal to 1.
    
    * *technical* is a dictionary including the parameters related to the module:
    
      * *T_r* is the rated temperature in °C.
      * *loss_coeff* is the loss coefficient (relevant for :mod:`physical_models.loss`).
      * *tracking* is either 0 for no tracking, 1 for one-axis tracking, or 2 for two-axes tracking.
      * *orientation* is the azimuth orientation of the module in degrees relative to the equator.
      * The tilt angle from the horizontal is chosen optimally in the code, see :mod:`physical_models.angles`.
    
    * *mask* is a dictionary including the parameters related to the masking:
    
      * *slope* is the threshold slope in percent. Areas with a larger slope are excluded.
      * *lu_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of land use types.
      * *pa_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of protected area categories.
      
    * *weight* is a dictionary including the parameters related to the weighting:
    
      * *GCR* is a dictionary of design settings for the ground cover ratio:
      
        * *shadefree_period* is the number of shadefree hours in the design day.
        * *day_north* is the design day for the northern hemisphere.
        * *day_south* is the design day for the southern hemisphere.
        
      * *lu_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of land use types.
      * *pa_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of protected area categories.
      * *power_density* is the power density of PV projects, assuming a GCR = 1, in MW/m².
      * *f_performance* is a number smaller than 1, taking into account all the other losses from the module until the AC substation.
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    :raise Tracking Warning: If *tracking* is not set to 0 and *orientation* is given as a value other than 0 or 180 (South or North), the orientation is ignored.

    """
    openfieldpv = {}
    openfieldpv["resource"] = {"clearness_correction": 1}
    openfieldpv["technical"] = {
        "T_r": 25,  # °C
        "loss_coeff": 0.37,
        "tracking": 0,  # 0 for no tracking, 1 for one-axis tracking, 2 for two-axes tracking
        "orientation": 0,  # | 0: Towards equator | 90: West | 180: Away from equator | -90: East |
    }
    openfieldpv["mask"] = {
        "slope": 10,
        "lu_suitability": np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,1,1,1,0,0]),
        "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
    GCR = {"shadefree_period": 6, "day_north": 79, "day_south": 263}
    openfieldpv["weight"] = {
        "GCR": GCR,
        "power_density": 0.000160,
        "f_performance": 0.75,
    }
    del GCR
    if openfieldpv["technical"]["tracking"] != 0 and openfieldpv["technical"]["orientation"] not in [0, 180]:
        warn("WARNING: " + str(openfieldpv["technical"]["tracking"]) + " axis tracking, overwrites orientation input: " + str(
            openfieldpv["technical"]["orientation"]))
        openfieldpv["technical"]["orientation"] = "track_" + str(openfieldpv["technical"]["tracking"])
    param["OpenFieldPV"] = openfieldpv
    return param


def rooftoppv_parameters(param):
    """
    This function sets the parameters for photovoltaics in the dictionary *pv* inside param:

    * *resource* is a dictionary including the parameters related to the resource potential:

      * *clearness_correction* is a factor that will be multiplied with the clearness index matrix to correct it. If no correction is required, leave it equal to 1.

    * *technical* is a dictionary including the parameters related to the module:

      * *T_r* is the rated temperature in °C.
      * *loss_coeff* is the loss coefficient (relevant for :mod:`physical_models.loss`).
      * *tracking* is either 0 for no tracking, 1 for one-axis tracking, or 2 for two-axes tracking.
      * *orientation* is the azimuth orientation of the module in degrees relative to the equator.
      * The tilt angle from the horizontal is chosen optimally in the code, see :mod:`physical_models.angles`.

    * *mask* is a dictionary including the parameters related to the masking:

      * *slope* is the threshold slope in percent. Areas with a larger slope are excluded.
      * *lu_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of land use types.
      * *pa_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of protected area categories.

    * *weight* is a dictionary including the parameters related to the weighting:

      * *GCR* is a dictionary of design settings for the ground cover ratio:

        * *shadefree_period* is the number of shadefree hours in the design day.
        * *day_north* is the design day for the northern hemisphere.
        * *day_south* is the design day for the southern hemisphere.

      * *lu_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of land use types.
      * *pa_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of protected area categories.
      * *power_density* is the power density of PV projects, assuming a GCR = 1, in MW/m².
      * *f_performance* is a number smaller than 1, taking into account all the other losses from the module until the AC substation.

    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    :raise Tracking Warning: If *tracking* is not set to 0 and *orientation* is given as a value other than 0 or 180 (South or North), the orientation is ignored.

    """
    rooftoppv = {}
    rooftoppv["resource"] = {"clearness_correction": 1}
    rooftoppv["technical"] = {
        "T_r": 25,  # °C
        "loss_coeff": 0.37,
        "tracking": 0,  # 0 for no tracking, 1 for one-axis tracking, 2 for two-axes tracking
        "orientation": 0,  # | 0: Towards equator | 90: West | 180: Away from equator | -90: East |
    }
    rooftoppv["mask"] = {
        "slope": 10,
        "lu_suitability": np.array(
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1,
             0, 0]),
        "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
    GCR = {"shadefree_period": 6, "day_north": 79, "day_south": 263}
    rooftoppv["weight"] = {
        "GCR": GCR,
        "suitable_roofs": 0.4,  # Percentage of south faced roofs. 40% should work for Europe!
        "power_density": 0.000160,
        "f_performance": 0.75,
    }
    del GCR
    if rooftoppv["technical"]["tracking"] != 0 and rooftoppv["technical"]["orientation"] not in [0, 180]:
        warn("WARNING: " + str(
            rooftoppv["technical"]["tracking"]) + " axis tracking, overwrites orientation input: " + str(
            rooftoppv["technical"]["orientation"]))
        rooftoppv["technical"]["orientation"] = "track_" + str(rooftoppv["technical"]["tracking"])
    param["RoofTopPV"] = rooftoppv
    return param


def csp_parameters(param):
    """
    This function sets the parameters for concentrated solar power in the dictionary *csp* inside param:
    
    * *resource* is a dictionary including the parameters related to the resource potential:
    
      * *clearness_correction* is a factor that will be multiplied with the clearness index matrix to correct it. If no correction is required, leave it equal to 1.
    
    * *technical* is a dictionary including the parameters related to the module:
    
      * *T_avg_HTF* is the average temperature in °C of the heat transfer fluid between the inlet and outlet of the solar field.
      * *loss_coeff* is the the heat loss coefficient in W/(m²K), which does not depend on wind speed (relevant for :mod:`physical_models.calc_CF_solar`).
      * *loss_coeff_wind* is the the heat loss coefficient in W/(m²K(m/s)^0.6), which depends on wind speed (relevant for :mod:`physical_models.calc_CF_solar`).
      * *Flow_coeff* is a factor smaller than 1 for the heat transfer to the HTF (Flow or heat removal factor).
      * *AbRe_ratio* is the ratio between the receiver area and the concentrator aperture.
      * *Wind_cutoff* is the maximum wind speed for effective tracking in m/s.
    
    * *mask* is a dictionary including the parameters related to the masking:
    
      * *slope* is the threshold slope in percent. Areas with a larger slope are excluded.
      * *lu_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of land use types.
      * *pa_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of protected area categories.
      
    * *weight* is a dictionary including the parameters related to the weighting:
        
      * *lu_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of land use types.
      * *pa_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of protected area categories.
      * *power_density* is the power density of CSP projects in MW/m².
      * *f_performance* is a number smaller than 1, taking into account all the other losses from the CSP module until the AC substation.
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    csp = {}
    csp["resource"] = {"clearness_correction": 1}
    csp["technical"] = {
        "T_avg_HTF": 350,  # °C
        "loss_coeff": 1.06,  # Heat Loss coefficient W/(m².K) independent of Wind Speed
        "loss_coeff_wind": 1.19,  # Multiplied with (Wind speed)^(0.6)
        "Flow_coeff": 0.95,  # heat transfer to the HTF factor (Flow or heat removal factor)
        "AbRe_ratio": 0.00079,  # Receiver Area / Concentrator aperture (90mm diameter/8m aperture)
        "Wind_cutoff": 22,  # (m/s)  Maximum Wind Speed for effective tracking (Source:NREL)
    }
    csp["mask"] = {
        "slope": 20,
        "lu_suitability": np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,1,1,1,0,0]),
        "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
    csp["weight"] = {
        "power_density": 0.000160,
        "f_performance": 0.9 * 0.75,
    }
    param["CSP"] = csp
    return param


def onshore_wind_parameters(param):
    """
    This function sets the parameters for onshore wind in the dictionary *windon* inside param:
    
    * *resource* is a dictionary including the parameters related to the resource potential:
    
      * *res_correction* is either 1 (perform a redistribution of wind speed when increasing the resolution) or 0 (repeat the same value from the low resolution data). It is relevant for :mod:`correction_functions.generate_wind_correction`.
      * *topo_correction* is either 1 (perform a correction of wind speed based on the altitude and the Global Wind Atlas) or 0 (no correction based on altitude).
      * *topo_weight* is only relevant if *topo_correction* = 1. It defines how to weight the correction factors of each country. There are three options: ``'none'`` (all countries have the same weight), ``'size'`` (larger countries have a higher weight), or ``'capacity'`` (countries with a higher installed capacity according to IRENA have a higher weight).
    
    * *technical* is a dictionary including the parameters related to the wind turbine:
    
      * *w_in* is the cut-in speed in m/s.
      * *w_r* is the rated wind speed in m/s.
      * *w_off* is the cut-off wind speed in m/s.
      * *P_r* is the rated power output in MW.
      * *hub_height* is the hub height in m.
    
    * *mask* is a dictionary including the parameters related to the masking:
    
      * *slope* is the threshold slope in percent. Areas with a larger slope are excluded.
      * *lu_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of land use types.
      * *pa_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of protected area categories.
      * *buffer_pixel_amount* is an integer that defines the number of pixels making a buffer of exclusion around urban areas.
      
    * *weight* is a dictionary including the parameters related to the weighting:
        
      * *lu_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of land use types.
      * *pa_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of protected area categories.
      * *power_density* is the power density of onshore wind projects in MW/m².
      * *f_performance* is a number smaller than 1, taking into account all the other losses from the turbine generator until the AC substation.
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    windon = {}
    # windon["technical"] = {"w_in": 2, "w_r": 14, "w_off": 25, "P_r": 2.3, "hub_height": 99}
    windon["technical"] = {"w_in": 4, "w_r": 13, "w_off": 25, "P_r": 3, "hub_height": 80}
    # windon["technical"] = {"w_in": 3, "w_r": 10.7, "w_off": 22, "P_r": 3.5, "hub_height": 100}
    windon["mask"] = {
        "slope": 17,
        "lu_suitability": np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,0,0]),
        "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
    windon["weight"] = {
        "power_density": 0.000008,
        "f_performance": 0.87,
    }
    param["WindOn"] = windon
    return param


def offshore_wind_paramters(param):
    """
    This function sets the parameters for offshore wind in the dictionary *windoff* inside param:
    
    * *resource* is a dictionary including the parameters related to the resource potential:
    
      * *res_correction* is either 1 (perform a redistribution of wind speed when increasing the resolution) or 0 (repeat the same value from the low resolution data).
        It is relevant for :mod:`correction_functions.generate_wind_correction`.
    
    * *technical* is a dictionary including the parameters related to the wind turbine:
    
      * *w_in* is the cut-in speed in m/s.
      * *w_r* is the rated wind speed in m/s.
      * *w_off* is the cut-off wind speed in m/s.
      * *P_r* is the rated power output in MW.
      * *hub_height* is the hub height in m.
    
    * *mask* is a dictionary including the parameters related to the masking:
    
      * *depth* is the threshold depth in meter (negative number). Areas that are deeper are excluded.
      * *pa_suitability* is a numpy array of values 0 (unsuitable) or 1 (suitable). It has the same size as the array of protected area categories.
      
    * *weight* is a dictionary including the parameters related to the weighting:

      * *power_density* is the power density of offshore wind projects in MW/m².
      * *f_performance* is a number smaller than 1, taking into account all the other losses from the turbine generator until the AC substation.
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    windoff = {}
    # windoff["technical"] = {"w_in": 3, "w_r": 16.5, "w_off": 34, "P_r": 7.58, "hub_height": 100}
    windoff["technical"] = {"w_in": 3, "w_r": 13, "w_off": 25, "P_r": 7, "hub_height": 150}
    windoff["mask"] = {
        "depth": -50,
        "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    }
    windoff["weight"] = {
        "power_density": 0.000003,
        "f_performance": 0.87,
    }
    param["WindOff"] = windoff
    return param


def biomass_parameters(param):
    """
    This function sets the parameters for biomass in the dictionary *Biomass* inside param:
   
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    biomass = {}
    biomass["agriculture"] = {
        "crops": np.array(
            [
                "Rice, paddy",
                "Maize",
                "Sugar cane",
                "Oil palm fruit",
                "Cassava",
                "Coconuts",
                "Coffee, green",
                "Groundnuts, with shell",
                "Sugar beet",
                "Wheat",
                "Barley",
                "Rye",
                "Rapeseed",
                "Oats"
            ]
        ),
        "residue": {
            "Rice, paddy": ["Straw", "Husk"],
            "Maize": ["Stalk", "Cob"],
            "Sugar cane": ["Bagasse", "Top&Leave"],
            "Oil palm fruit": ["Shell", "Fiber", "EFB", "POME", "Frond"],
            "Cassava": ["Stalk"],
            "Coconuts": ["Husk", "Shell", "Frond"],
            "Coffee, green": ["Husk"],
            "Groundnuts, with shell": ["Shell", "Straw"],
            "Sugar beet": ["Leaves"],
            "Wheat": ["Straw"],
            "Barley": ["Straw"],
            "Rye": ["Straw"],
            "Rapeseed": ["Straw"],
            "Oats": ["Straw"]
        },
        "rpr": {
            "Rice, paddy": {"Straw": 1.00, "Husk": 0.27},
            "Maize": {"Stalk": 1.00, "Cob": 0.25},
            "Sugar cane": {"Bagasse": 0.25, "Top&Leave": 0.30},
            "Oil palm fruit": {"Shell": 0.07, "Fiber": 0.13, "EFB": 0.23, "POME": 0.67, "Frond": 0.55},
            "Cassava": {"Stalk": 0.09},
            "Coconuts": {"Husk": 0.36, "Shell": 0.16, "Frond": 0.23},
            "Coffee, green": {"Husk": 2.10},
            "Groundnuts, with shell": {"Shell": 0.32, "Straw": 2.30},
            "Sugar beet": {"Leaves": 0.3},
            "Wheat": {"Straw": 1.2},
            "Barley": {"Straw": 1.1},
            "Rye": {"Straw": 1.3},
            "Rapeseed": {"Straw": 1.7},
            "Oats": {"Straw": 1.2}
        },
        "af": {
            "Rice, paddy": {"Straw": 0.5, "Husk": 0.47},
            "Maize": {"Stalk": 0.33, "Cob": 0.67},
            "Sugar cane": {"Bagasse": 0.21, "Top&Leave": 0.99},
            "Oil palm fruit": {"Shell": 0.04, "Fiber": 0.13, "EFB": 0.58, "POME": 0.65, "Frond": 0.05},
            "Cassava": {"Stalk": 0.41},
            "Coconuts": {"Husk": 0.60, "Shell": 0.38, "Frond": 0.81},
            "Coffee, green": {"Husk": 0.33},
            "Groundnuts, with shell": {"Shell": 1.00, "Straw": 0.33},
            "Sugar beet": {"Leaves": 0.4},
            "Wheat": {"Straw": 0.4},
            "Barley": {"Straw": 0.4},
            "Rye": {"Straw": 0.4},
            "Rapeseed": {"Straw": 0.5},
            "Oats": {"Straw": 0.4}
        },
        "lhv": {
            "Rice, paddy": {"Straw": 3.89, "Husk": 3.57},
            "Maize": {"Stalk": 3.97, "Cob": 4.62},
            "Sugar cane": {"Bagasse": 1.79, "Top&Leave": 1.89},
            "Oil palm fruit": {"Shell": 4.72, "Fiber": 3.08, "EFB": 1.69, "POME": 0.17, "Frond": 2.21},
            "Cassava": {"Stalk": 4.72},
            "Coconuts": {"Husk": 4.09, "Shell": 4.58, "Frond": 4.04},
            "Coffee, green": {"Husk": 3.44},
            "Groundnuts, with shell": {"Shell": 3.12, "Straw": 4.88},
            "Sugar beet": {"Leaves": 3.7},
            "Wheat": {"Straw": 4},
            "Barley": {"Straw": 4},
            "Rye": {"Straw": 4},
            "Rapeseed": {"Straw": 4},
            "Oats": {"Straw": 4}
        },  # MWh/ton
        "emission factor" : 1585  # kg/ton of crop residues
    }

    biomass["forest"] = {
        "woods" : np.array(["Wood fuel, coniferous","Wood fuel, non-coniferous"]),
        "density" : {"Wood fuel, coniferous": 0.75,"Wood fuel, non-coniferous":  0.85},  # tons/m3
        "rpr": 0.67,
        "af": 0.4,
        "lhv": 4.31,  # MWh/ton
        "emission factor": 1500  # kg/ton of wood
    }

    biomass["livestock"] = {
        "animal": np.array(["Cattle", "Sheep", "Goats", "Chickens", "Pigs"]),
        "rpr": np.array([0.84, 0.12, 0.13, 0.007, 0.11]),  # ton of manure per animal
        "af": np.array([0.02, 0.02, 0.02, 0.47, 0.47]),
        "lhv": np.array([1.01, 1.31, 1.31, 2.43, 2.93]),  # MWh/ton of manure
        "emission factor": 859  # kg/ton of manure
    }

    param["Biomass"] = biomass
    return param


###########################
##### Define Paths ########
###########################


def weather_input_folder(paths, param):
    """
    This function defines the path *MERRA_IN* where the MERRA-2 data is saved. It depends on the coverage of the data and the year.
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    global root
    global fs

    MERRA_coverage = param["MERRA_coverage"]
    year = str(param["year"])
    paths["MERRA_IN"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + MERRA_coverage + " " + year + fs
    paths["GWA_global"] = root + "04 Global Wind Atlas Wind Speed" + fs + param["country_code"] + "_wind-speed_50m.tif"

    return paths


def global_maps_input_paths(paths, param):
    """
    This function defines the paths where the global maps are saved:
    
      * *LU_global* for the land use raster
      * *Topo_tiles* for the topography tiles (rasters)
      * *Pop_global* for the global population raster
      * *Bathym_global* for the bathymetry raster
      * *Protected* for the shapefile of protected areas
      * *GWA* for the country data retrieved from the Global Wind Atlas (missing the country code, which will be filled in a for-loop in :mod:correction_functions.calc_gwa_correction)
      * *Countries* for the shapefiles of countries
      * *EEZ_global* for the shapefile of exclusive economic zones of countries
    
    :param paths: Dictionary including the paths.
    :type paths: dict

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    global root
    global fs

    # Global maps
    PathTemp = root + "01 Raw inputs" + fs + "Maps" + fs
    
    paths["Countries"] = PathTemp + "Countries" + fs + "gadm36_0.shp"
    
    paths["EEZ_global"] = PathTemp + "EEZ" + fs + "eez_v10.shp"
    paths["Internal_waters_global"] = PathTemp + "Internal Waters" + fs + "eez_internal_waters_v3.shp"
    paths["Topo_global"] = PathTemp + "Topography" + fs + "Topography_250m.tif"  # [1/400 , 1/400] resampled in QGIS
    paths["Bathym_global"] = PathTemp + "Bathymetry" + fs + "ETOPO1_Ice_c_geotiff_250m.tif" # [1/400 , 1/400] resampled in QGIS
    paths["LU_global"] = PathTemp + "Landuse" + fs + "CCI-250m.tif" #[1/400 , 1/400] resampled in QGIS
    paths["Protected"] = PathTemp + "Protected Areas" + fs + "WDPA_Nov2018-shapefile-polygons.shp"
    paths["Airports"] = PathTemp + "Openflights" + fs + "airports.csv"
    paths["OSM_Roads"] = PathTemp + "OSM" + fs + param["country_code"] + "-roads.shp"
    paths["OSM_Railways"] = PathTemp + "OSM" + fs + param["country_code"] + "-railways.shp"
    paths["OSM_Landuse"] = PathTemp + "OSM" + fs + param["country_code"] + "-landuse.shp"
    paths["WSF_global"] = PathTemp + "WSF" + fs + "WSF2015_Full.tif" #[1/400 , 1/400] resampled in QGIS
    paths["HydroLakes"] = PathTemp + "HydroLakes" + fs + "HydroLAKES_polys_v10.shp"
    paths["HydroRivers"] = PathTemp + "HydroRivers" + fs + "HydroRIVERS_v10.shp"

    paths["LS_global"] = PathTemp + "Livestock" + fs + "Glb_"
    param["res_livestock"] = np.array([1 / 120, 1 / 120])

    # paths["Pop_global"] = PathTemp + "Population" + fs + "gpw_v4_population_count_rev10_2015_30_sec.tif"
    # param["res_population"] = np.array([1 / 120, 1 / 120])

    paths["Biomass_Crops"] = PathTemp + "FAOSTAT" + fs + "FAOSTAT_data_9-28-2021.csv"
    paths["Biomass_Forestry"] = PathTemp + "FAOSTAT" + fs + "FAOSTAT_Forestry_data_6-2-2021.csv"

    return paths, param


def output_folders(paths, param):
    """
    This function defines the paths to multiple output folders:
    
      * *region* is the main output folder.
      * *weather_data* is the output folder for the weather data of the spatial scope.
      * *local_maps* is the output folder for the local maps of the spatial scope.
      * *potential* is the output folder for the ressource and technical potential maps.
      * *regional_analysis* is the output folder for the time series and the report of the subregions.
      * *regression_in* is the folder where the regression parameters (FLH, fitting time series) are saved.
      * *regression_out* is the output folder for the regression results.
      
    All the folders are created at the beginning of the calculation, if they do not already exist,
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    global root
    global fs

    region = param["region_name"]
    subregions = param["subregions_name"]

    # Main output folder
    # paths["region"] = root + "03 Intermediate files" + fs + "Files " + region + fs        # old structure of database
    paths["region"] = os.path.join(root, "..", "1_results", "Files " + region, "")

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

    # Output folder for discrete analysis
    paths["discrete_analysis"] = paths["region"] + "Renewable energy" + fs + "Discrete analysis" + fs
    if not os.path.isdir(paths["discrete_analysis"]):
        os.makedirs(paths["discrete_analysis"])

    # Regression parameters
    paths["regression_in"] = paths["regional_analysis"] + "Regression outputs" + fs + "Parameters" + fs
    if not os.path.isdir(paths["regression_in"]):
        os.makedirs(paths["regression_in"])

    # Regression output
    paths["regression_out"] = paths["regional_analysis"] + "Regression outputs" + fs

    return paths


def weather_output_paths(paths, param):
    """
    This function defines the paths to weather filesfor a specific *year*:
    
      * *W50M* is the file for the wind speed at 50m in m/s.
      * *CLEARNESS* is the file for the clearness index, e.g. the ratio between total ground horizontal radiation and total top-of-the-atmosphere horizontal radiation.
      * *T2M* is the file for the temperature at 2m in Kelvin.
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    year = str(param["year"])
    paths["W50M"] = paths["weather_data"] + "w50m_" + year + ".mat"
    paths["CLEARNESS"] = paths["weather_data"] + "clearness_" + year + ".mat"
    paths["SWGDN"] = paths["weather_data"] + "SWGDN_" + year + ".mat"
    paths["SWTDN"] = paths["weather_data"] + "SWTDN_" + year + ".mat"
    paths["SWGDN_real"] = paths["weather_data"] + "SWGDN_real_" + year + ".mat"
    paths["T2M"] = paths["weather_data"] + "t2m_" + year + ".mat"
    paths["W50M_offshore"] = paths["weather_data"] + "w50m_offshore" + year + ".mat"

    paths["MERRA_XMIN"] = paths["weather_data"] + "w50m_xmin" + year + ".mat"
    paths["MERRA_XMAX"] = paths["weather_data"] + "w50m_xmax" + year + ".mat"
    paths["MERRA_YMIN"] = paths["weather_data"] + "w50m_ymin" + year + ".mat"
    paths["MERRA_YMAX"] = paths["weather_data"] + "w50m_ymax" + year + ".mat"

    paths["GWA_X"] = paths["weather_data"] + "gwa_x" + year + ".mat"
    paths["GWA_Y"] = paths["weather_data"] + "gwa_y" + year + ".mat"

    paths["IND_GWA_MERRA_X"] = paths["weather_data"] + "ind_gwa_merra_x" + ".mat"
    paths["IND_GWA_MERRA_Y"] = paths["weather_data"] + "ind_gwa_merra_y" + ".mat"

    return paths


def local_maps_paths(paths, param):
    """
    This function defines the paths where the local maps will be saved:
    
      * *LAND* for the raster of land areas within the scope
      * *EEZ* for the raster of sea areas within the scope
      * *SUB* for the raster of areas covered by subregions (both land and sea) within the scope
      * *LU* for the land use raster within the scope
      * *BATH* for the bathymetry raster within the scope
      * *TOPO* for the topography raster within the scope
      * *SLOPE* for the slope raster within the scope
      * *PA* for the raster of protected areas within the scope
      * *POP* for the population raster within the scope
      * *BUFFER* for the raster of population buffer areas within the scope
      * *CORR_GWA* for correction factors based on the Global Wind Atlas (mat file)
      * *CORR_ON* for the onshore wind correction factors (raster)
      * *CORR_OFF* for the offshore wind correction factors (raster)
      * *AREA* for the area per pixel in m² (mat file)
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    # Local maps
    PathTemp = paths["local_maps"] + param["region_name"]
    paths["LAND"] = PathTemp + "_Land.tif"  # Land pixels
    paths["EEZ"] = PathTemp + "_EEZ.tif"  # Sea pixels
    paths["INT_WATER"] = PathTemp + "_INT_WATER.tif"  # Internal Waters pixels
    paths["AREA"] = PathTemp + "_Area.mat"  # Area per pixel in m²
    paths["AREA_offshore"] = PathTemp + "_Area_offshore.mat"  # Area per pixel in m²
    paths["TOPO"] = PathTemp + "_Topography.mat"  # Topography
    paths["TOPO_low"] = PathTemp + "_Topography_low.mat"
    paths["SLOPE"] = PathTemp + "_Slope.mat"  # Slope
    paths["BATH"] = PathTemp + "_Bathymetry.mat"  # Bathymetry

    paths["LU"] = PathTemp + "_Landuse.mat"  # Land use types
    paths["WATER_BUFFER"] = PathTemp + "_Water_Buffered.mat"  # Buffered Water
    paths["WETLAND_BUFFER"] = PathTemp + "_Wetland_Buffered.mat"  # Buffered Wetlands
    paths["SNOW_BUFFER"] = PathTemp + "_Snow_Buffered.mat"  # Buffered Snow

    paths["PA"] = PathTemp + "_Protected_areas.tif"  # Protected areas
    paths["PV_PA_BUFFER"] = PathTemp + "_PV_Protected_areas_Buffered.mat"  # Buffered Protected areas for PV
    paths["WINDON_PA_BUFFER"] = PathTemp + "_WindOn_Protected_areas_Buffered.mat"  # Buffered Protected areas for Wind Onshore
    paths["PA_offshore"] = PathTemp + "_Protected_areas_offshore.tif"  # Protected areas
    paths["WINDOFF_PA_BUFFER"] = PathTemp + "_WindOff_Protected_areas_Buffered.mat"  # Buffered Protected areas for Wind Offshore

    paths["AIRPORTS"] = PathTemp + "_Airports.mat"  # Buffered Airports
    paths["BOARDERS"] = PathTemp + "_Boarders.mat"  # Buffered Boarders

    paths["ROADS"] = PathTemp + "_Roads.tif"
    paths["RAILS"] = PathTemp + "_Rails.tif"
    paths["OSM_AREAS"] = PathTemp + "_OSM_areas.tif"
    paths["OSM_COM_BUFFER"] = PathTemp + "_OSM_Commercial_Buffered.mat"
    paths["OSM_IND_BUFFER"] = PathTemp + "_OSM_Industrial_Buffered.mat"
    paths["OSM_MINE_BUFFER"] = PathTemp + "_OSM_Mining_Buffered.mat"
    paths["OSM_MIL_BUFFER"] = PathTemp + "_OSM_Military_Buffered.mat"
    paths["PV_OSM_PARK_BUFFER"] = PathTemp + "_PV_OSM_Parks_Buffered.mat"
    paths["WINDON_OSM_PARK_BUFFER"] = PathTemp + "_WindOn_OSM_Parks_Buffered.mat"
    paths["OSM_REC_BUFFER"] = PathTemp + "_OSM_Recreation_Buffered.mat"

    paths["WSF"] = PathTemp + "_Settlements.mat"
    paths["PV_WSF_BUFFER"] = PathTemp + "_PV_Settlements_Buffered.mat"
    paths["WINDON_WSF_BUFFER"] = PathTemp + "_WindOn_Settlements_Buffered.mat"

    paths["HYDROLAKES"] = PathTemp + "_HydroLakes.tif"
    paths["HYDROLAKES_BUFFER"] = PathTemp + "_HydroLakes_Buffered.mat"
    paths["HYDRORIVERS"] = PathTemp + "_HydroRivers.tif"
    paths["HYDRORIVERS_BUFFER"] = PathTemp + "_HydroRivers_Buffered.mat"

    paths["LS"] = PathTemp + "_Livestock_"  # Livestock density per animal type
    # paths["POP"] = PathTemp + "_Population.tif"  # Population

    # Correction factors for wind speeds
    turbine_height_on = str(param["WindOn"]["technical"]["hub_height"])
    turbine_height_off = str(param["WindOff"]["technical"]["hub_height"])
    paths["CORR_ON"] = PathTemp + "_WindOn_Correction_" + turbine_height_on + ".tif"
    paths["CORR_OFF"] = PathTemp + "_WindOff_Correction_" + turbine_height_off + ".tif"

    return paths


def irena_paths(paths, param):
    """
    This function defines the paths for the IRENA inputs and outputs:
        
      * *IRENA* is a csv file containing statistics for all countries and technologies for a specific *year*, created using a query tool of IRENA.
      * *IRENA_dict* is a csv file to convert the code names of countries from the IRENA database to the database of the shapefile of countries.
      * *IRENA_summary* is a csv file with a summary of renewable energy statistics for the countries within the scope.
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    global root
    global fs

    year = str(param["year"])

    # IRENA input
    paths["IRENA"] = (
            root + "01 Raw inputs" + fs + "Renewable energy" + fs + "IRENA" + fs + "IRENA_RE_electricity_statistics_allcountries_alltech_" + year + ".csv"
    )

    current_folder = os.path.dirname(os.path.abspath(__file__))
    PathTemp = str(Path(current_folder).parent)
    if PathTemp[-1] != fs:
        PathTemp = PathTemp + fs
    PathTemp = PathTemp + "assumptions"

    paths["IRENA_dict"] = PathTemp + fs + "dict_countries.csv"

    # IRENA output
    paths["IRENA_summary"] = paths["region"] + "Renewable energy" + fs + "IRENA_summary_" + year + ".csv"

    return paths


def regression_paths(paths, param, tech):
    """
    This function defines the paths for the regression parameters:
    
      * *FLH_regression* is a csv file containing FLH statistics for the subregions and the four technologies for a specific *year*, based on the previously created *IRENA_summary*.
      * *TS_regression* is a csv file containing time series to be match for each subregion and technology, based on EMHIRES time series if available.
    
    :param paths: Dictionary including the paths.
    :type paths: dict

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    year = str(param["year"])

    paths["FLH_regression"] = paths["regression_in"] + "FLH_regression_" + year + ".csv"
    paths[tech]["TS_regression"] = paths["regression_in"] + "TimeSeries_regression_" + tech + "_" + year + ".csv"

    return paths


def emhires_input_paths(paths, tech):
    """
    This function defines the path to the EMHIRES input file for each technology (only ``'WindOn'``,
    ``'WindOff'``, and ``'PV'`` are supported by EMHIRES).
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict
    :param tech: Name of the technology.
    :type tech: string

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    global root
    global fs

    if tech == "WindOn":
        paths[tech][
            "EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES" + fs + "TS.CF.COUNTRY.30yr.date.txt"
    elif tech == "WindOff":
        paths[tech][
            "EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES" + fs + "TS.CF.OFFSHORE.30yr.date.txt"
    elif tech == "PV":
        paths[tech][
            "EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES" + fs + "EMHIRESPV_TSh_CF_Country_19862015.txt"

    return paths


def potential_output_paths(paths, param, tech):
    """
    This function defines the paths of the files that will be saved in the folder for the potential outputs:
    
      * *FLH* is the file with the full-load hours for all pixels within the scope (mat file).
      * *mask* is the file with the suitable pixels within the scope (mat file).
      * *FLH_mask* is the file with the full-load hours for the suitable pixels within the scope (mat file).
      * *weight* is the power density for all the pixels in the scope (mat file).
      * *FLH_weight* is the potential energy output for all the pixels in the scope (mat file).
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict
    :param tech: Name of the technology.
    :type tech: string

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    region = param["region_name"]
    year = str(param["year"])

    PathTemp = ''
    if tech in ["WindOn", "WindOff"]:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["potential"] + region + "_" + tech + "_" + hubheight
    elif tech in ["OpenFieldPV", "RoofTopPV"]:
        if "orientation" in param[tech]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = "0"
        PathTemp = paths["potential"] + region + "_" + tech + "_" + orientation
    elif tech in ["CSP"]:
        orientation = "0"
        PathTemp = paths["potential"] + region + "_" + tech + "_" + orientation
    elif tech in ["Biomass"]:
        PathTemp = paths["potential"] + region + "_" + tech

    # File name for Biomass
    if tech in ["Biomass"]:
        paths[tech]["BIOMASS_CROPS"] = PathTemp + "_Biomass_Crops.mat"
        paths[tech]["BIOMASS_WOOD"] = PathTemp + "_Biomass_Wood.mat"
        paths[tech]["BIOMASS_MANURE"] = PathTemp + "_Biomass_Manure.mat"
        paths[tech]["BIOMASS_ENERGY"] = PathTemp + "_Biomass_Energy.mat"
        paths[tech]["BIOMASS_CO2"] = PathTemp + "_Biomass_CO2.mat"
    # File name for all other tech
    else:
        paths[tech]["FLH"] = PathTemp + "_FLH_" + year + ".mat"
        paths[tech]["mask"] = PathTemp + "_mask_" + year + ".mat"
        paths[tech]["FLH_mask"] = PathTemp + "_FLH_mask_" + year + ".mat"
        paths[tech]["weight"] = PathTemp + "_weight_" + year + ".mat"
        paths[tech]["FLH_weight"] = PathTemp + "_FLH_weight_" + year + ".mat"

    return paths


def regional_analysis_output_paths(paths, param, tech):
    """
    This function defines the paths of the files that will be saved in the folder for the regional analysis outputs:
    
      * *Locations* is the shapefile of points that correspond to the selected quantiles in each subregion, for which the time series will be generated.
      * *TS* is the csv file with the time series for all subregions and quantiles.
      * *Region_Stats* is the csv file with the summary report for all subregions.
      * *Sorted_FLH* is the mat file with the sorted samples of FLH for each subregion.
      * *Regression_coefficients* is the path format for a csv files containing the regression coefficients found by the solver
      * *Regression_TS* is the path format for a csv files with the regression resulting timeseries for the tech and settings
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict
    :param tech: Name of the technology.
    :type tech: string

    :return paths: The updated dictionary paths.
    :rtype: dict
    """
    subregions = param["subregions_name"]
    year = str(param["year"])

    PathTemp = ''
    if tech in ["WindOn", "WindOff"]:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["regional_analysis"] + subregions + "_" + tech + "_" + hubheight
    elif tech in ["OpenFieldPV", "RoofTopPV"]:
        if "orientation" in param[tech]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = "0"
        PathTemp = paths["regional_analysis"] + subregions + "_" + tech + "_" + orientation
    elif tech in ["CSP"]:
        orientation = "0"
        PathTemp = paths["regional_analysis"] + subregions + "_" + tech + "_" + orientation
    elif tech in ["Biomass"]:
        PathTemp = paths["regional_analysis"] + subregions + "_" + tech

    paths[tech]["Region_Stats"] = PathTemp + "_Region_stats_" + year + ".csv"

    if tech != "Biomass":
        paths[tech]["Locations"] = PathTemp + "_Locations_" + year + ".shp"
        paths[tech]["TS"] = PathTemp + "_TS_" + year + ".csv"
        paths[tech]["Sorted_FLH"] = PathTemp + "_sorted_FLH_sampled_" + year + ".mat"
        paths[tech]["Regression_coefficients"] = paths["regression_out"] + subregions + "_" + tech + "_reg_coefficients_"
        paths[tech]["Regression_TS"] = paths["regression_out"] + subregions + "_" + tech + "_reg_TimeSeries_"

    return paths


def discrete_output_paths(paths, param, tech):
    region = param["region_name"]
    year = str(param["year"])

    PathTemp = ''
    if tech in ["WindOn", "WindOff"]:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["discrete_analysis"] + region + "_" + tech + "_" + hubheight
    elif tech in ["OpenFieldPV", "RoofTopPV"]:
        if "orientation" in param[tech]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = "0"
        PathTemp = paths["discrete_analysis"] + region + "_" + tech + "_" + orientation
    elif tech in ["CSP"]:
        orientation = "0"
        PathTemp = paths["discrete_analysis"] + region + "_" + tech + "_" + orientation

    if tech != "Biomass":
        paths[tech]["TS_discrete"] = PathTemp + "_TS_" + year + ".csv"

    return paths
