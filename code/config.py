import os
from pathlib import Path
from warnings import warn

import numpy as np


def configuration():
    """
    This function is the main configuration function that calls all the other modules in the code.

    :return (paths, param): The dictionary paths containing all the paths to inputs and outputs, and the dictionary param containing all the user preferences.
    :rtype: tuple(dict, dict)
    """
    paths, param = general_settings()
    paths, param = scope_paths_and_parameters(paths, param)

    param = computation_parameters(param)
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
    paths, param = global_maps_input_paths(paths, param)
    paths = output_folders(paths, param)
    paths = weather_output_paths(paths, param)
    paths = local_maps_paths(paths, param)
    paths = irena_paths(paths, param)

    for tech in param["technology"]:
        paths[tech] = {}
        paths = regression_paths(paths, param, tech)
        paths = emhires_input_paths(paths, param, tech)
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

    param = {}
    param["author"] = "User"  # the name of the person running the script
    param["comment"] = "Tutorial"

    paths = {}
    fs = os.path.sep
    current_folder = os.path.dirname(os.path.abspath(__file__))
    root = str(Path(current_folder).parent.parent)
    # For use at TUM ENS
    if root[-1] != fs:
        root = root + fs + "Database" + fs
    else:
        root = root + "Database" + fs

    return paths, param


###########################
#### User preferences #####
###########################


def scope_paths_and_parameters(paths, param):
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
    
    * *res_desired* is a numpy array with two numbers. The first number is the resolution in the vertical dimension (in degrees of latitude),
      the second is for the horizontal dimension (in degrees of longitude).

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
    PathTemp = root + "02 Shapefiles for regions" + fs

    paths["spatial_scope"] = PathTemp + "User-defined" + fs + "gadm36_AUT_0.shp"
    paths["subregions"] = PathTemp + "Clustering outputs" + fs + "Austria" + fs + "Wind_FLH - Solar_FLH" + fs + "05 final_output" + fs + "final_result.shp"

    # Name tags for the scope and the subregions
    param["region_name"] = "Austria"  # Name tag of the spatial scope
    param["subregions_name"] = "Austria"  # Name tag of the subregions

    # Desired resolution
    param["res_desired"] = np.array([1 / 240, 1 / 240])

    # Year
    param["year"] = 2015

    # Technologies
    param["technology"] = ["WindOn", "PV"]  # ["PV", "CSP", "WindOn", "WindOff"]

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
    param["nproc"] = 6
    param["CPU_limit"] = True
    return param


def weather_data_parameters(param):
    """
    This function defines the coverage of the weather data *MERRA_coverage*, and how outliers should be corrected using *MERRA_correction*:
    
    * *MERRA_coverage*: If you have downloaded the MERRA-2 data for the world, enter the name tag ``'World'``. The code will later search for the data in the corresponding folder.
      It is possible to download the MERRA-2 just for the geographic scope of the analysis. In that case, enter another name tag (we recommend using the same one as the spatial scope).
    
    * *MERRA_correction*: MERRA-2 contains some outliers, especially in the wind data. *MERRA_correction* decides whether these outliers are dealt with.

    * *MERRA_correction_factor*: if *MERRA_correction* is ``'True'``, this sets the threshold of the relative distance between the yearly mean of the data point
      to the yearly mean of its neighbors.

    * *res_weather*: defines the resolution of weather data using a numpy array with two numbers. The first number is the resolution in the vertical dimension (in degrees of latitude),
    the second is for the horizontal dimension (in degrees of longitude).

    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    param["MERRA_coverage"] = "World"
    param["MERRA_correction"] = True
    param["MERRA_correction_factor"] = {"W50M": 0.35, "CLEARNESS": 0.35, "T2M": 0.35}  # Wind Speed  # Clearness index  # Temperature at 2 m
    param["res_weather"] = np.array([1 / 2, 5 / 8])
    return param


def file_saving_options(param):
    """
    This function sets some options for saving files.
    
    * *savetiff* is a boolean that determines whether tif rasters for the potentials are saved (``True``), or whether only mat files are saved (``False``).
      The latter are saved in any case.
    
    * *report_sampling* is an integer that sets the sample size for the sorted FLH values per region (relevant for :mod:`potential.report_potentials`).
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    # Mask / Weight
    param["savetiff"] = True  # Save geotiff files of mask and weight rasters

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
    param["quantiles"] = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0]

    # User defined locations
    param["useloc"] = {"Point1": (0, -80), "Point2": (1, 1)}  # {"point name": (latitude, longitude),...}

    # Regression
    param["regression"] = {
        "solver": "gurobi",  # string
        "WindOn": {"all": []},  # dictionary of hub height combinations
        "WindOff": {"80m": []},  # dictionary of hub height combinations
        "PV": {"all": []},  # list of orientation combinations
        "CSP": {"all": []},
    }

    # Stratified time series
    param["modes"] = {"all": param["quantiles"]}
    param["combo"] = {
        # dictionary of hub height and orientation combinations
        "WindOn": {"2015": [80]},
        "WindOff": {"80m": [80], "100m": [100], "120m": [120]},
        "PV": {"Solar": [0]},
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
      * *height* is a numpy array of gradient heights in meter associated to each land use type (relevant for :mod:`correction_functions.generate_wind_correction`).
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict

    Land use reclassification::

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
        # 13  -- Urban and built-up
        # 14  -- Croplands / natural vegetation mosaic
        # 15  -- Snow and ice
        # 16  -- Barren or sparsely vegetated
    """
    landuse = {
        "type": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]),
        "type_urban": 13,
        "Ross_coeff": np.array(
            [0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208, 0.0208]
        ),
        "albedo": np.array([0.00, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.00, 0.20, 0.20, 0.20, 0.00, 0.20]),
        "hellmann": np.array([0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.20, 0.20, 0.25, 0.25, 0.15, 0.15, 0.20, 0.40, 0.20, 0.15, 0.15]),
        "height": np.array([213, 366, 366, 366, 366, 366, 320, 320, 366, 366, 274, 274, 320, 457, 320, 274, 274]),
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


def pv_parameters(param):
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
    pv = {}
    pv["resource"] = {"clearness_correction": 1}
    pv["technical"] = {
        "T_r": 25,  # °C
        "loss_coeff": 0.37,
        "tracking": 0,  # 0 for no tracking, 1 for one-axis tracking, 2 for two-axes tracking
        "orientation": 0,  # | 0: Towards equator | 90: West | 180: Away from equator | -90: East |
    }
    pv["mask"] = {
        "slope": 20,
        "lu_suitability": np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1]),
        "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
    }
    GCR = {"shadefree_period": 6, "day_north": 79, "day_south": 263}
    pv["weight"] = {
        "GCR": GCR,
        "lu_availability": np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.02, 0.00, 0.02, 0.02, 0.02, 0.00, 0.02]),
        "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
        "power_density": 0.000160,
        "f_performance": 0.75,
    }
    del GCR
    if pv["technical"]["tracking"] != 0 and pv["technical"]["orientation"] not in [0, 180]:
        warn("WARNING: " + str(pv["technical"]["tracking"]) + " axis tracking, overwrites orientation input: " + str(pv["technical"]["orientation"]))
        pv["technical"]["orientation"] = "track_" + str(pv["technical"]["tracking"])
    param["PV"] = pv
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
        "lu_suitability": np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1]),
        "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
    }
    csp["weight"] = {
        "lu_availability": np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.02, 0.00, 0.02, 0.02, 0.02, 0.00, 0.02]),
        "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
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
    windon["resource"] = {"res_correction": 1, "topo_correction": 1, "topo_weight": "capacity"}  # 'none' or 'size' or 'capacity'
    windon["technical"] = {"w_in": 4, "w_r": 13, "w_off": 25, "P_r": 3, "hub_height": 80}
    windon["mask"] = {
        "slope": 20,
        "lu_suitability": np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1]),
        "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),
        "buffer_pixel_amount": 1,
    }
    windon["weight"] = {
        "lu_availability": np.array([0.00, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.10, 0.10, 0.10, 0.10, 0.00, 0.10, 0.00, 0.10, 0.00, 0.10]),
        "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
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
        
      * *lu_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of land use types.
      * *pa_availability* is a numpy array of values between 0 (completely not available) and 1 (completely available). It has the same size as the array of protected area categories.
      * *power_density* is the power density of offshore wind projects in MW/m².
      * *f_performance* is a number smaller than 1, taking into account all the other losses from the turbine generator until the AC substation.
    
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return param: The updated dictionary param.
    :rtype: dict
    """
    windoff = {}
    windoff["resource"] = {"res_correction": 1}
    windoff["technical"] = {"w_in": 3, "w_r": 16.5, "w_off": 34, "P_r": 7.58, "hub_height": 100}
    windoff["mask"] = {"depth": -40, "pa_suitability": np.array([1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1])}
    windoff["weight"] = {
        "lu_availability": np.array([0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]),
        "pa_availability": np.array([1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 1.00, 1.00, 1.00, 1.00]),
        "power_density": 0.000020,
        "f_performance": 0.87,
    }
    param["WindOff"] = windoff
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

    It also defines the resolution of the input rasters (*res_landuse*, *res_topography*, *res_population*, *res_bathymetry*) using numpy arrays with two numbers. The first number is the resolution in the vertical dimension (in degrees of latitude),
    the second is for the horizontal dimension (in degrees of longitude).

    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return (paths, param): The updated dictionaries paths and param.
    :rtype: tuple(dict, dict)
    """
    global root
    global fs

    # Global maps
    PathTemp = root + "01 Raw inputs" + fs + "Maps" + fs
    paths["LU_global"] = PathTemp + "Landuse" + fs + "LCType.tif"
    param["res_landuse"] = np.array([1 / 240, 1 / 240])

    paths["Topo_tiles"] = PathTemp + "Topography" + fs
    param["res_topography"] = np.array([1 / 240, 1 / 240])

    paths["Pop_global"] = PathTemp + "Population" + fs + "gpw_v4_population_count_rev10_2015_30_sec.tif"
    param["res_population"] = np.array([1 / 120, 1 / 120])

    paths["Bathym_global"] = PathTemp + "Bathymetry" + fs + "ETOPO1_Ice_c_geotiff.tif"
    param["res_bathymetry"] = np.array([1 / 60, 1 / 60])

    paths["Protected"] = PathTemp + "Protected Areas" + fs + "WDPA_Nov2018-shapefile-polygons.shp"
    paths["GWA"] = PathTemp + "Global Wind Atlas" + fs + fs + "windSpeed.csv"
    paths["Countries"] = PathTemp + "Countries" + fs + "gadm36_0.shp"
    paths["EEZ_global"] = PathTemp + "EEZ" + fs + "eez_v10.shp"

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
    paths["T2M"] = paths["weather_data"] + "t2m_" + year + ".mat"

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
    paths["SUB"] = PathTemp + "_Subregions.tif"  # Subregions pixels
    paths["LU"] = PathTemp + "_Landuse.tif"  # Land use types
    paths["TOPO"] = PathTemp + "_Topography.tif"  # Topography
    paths["PA"] = PathTemp + "_Protected_areas.tif"  # Protected areas
    paths["SLOPE"] = PathTemp + "_Slope.tif"  # Slope
    paths["BATH"] = PathTemp + "_Bathymetry.tif"  # Bathymetry
    paths["POP"] = PathTemp + "_Population.tif"  # Population
    paths["BUFFER"] = PathTemp + "_Population_Buffered.tif"  # Buffered population
    paths["CORR_GWA"] = PathTemp + "_GWA_Correction.mat"  # Correction factors based on the GWA
    paths["AREA"] = PathTemp + "_Area.mat"  # Area per pixel in m²

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


def emhires_input_paths(paths, param, tech):
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
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES" + fs + "TS.CF.COUNTRY.30yr.date.txt"
    elif tech == "WindOff":
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES" + fs + "TS.CF.OFFSHORE.30yr.date.txt"
    elif tech == "PV":
        paths[tech]["EMHIRES"] = root + "01 Raw inputs" + fs + "Renewable energy" + fs + "EMHIRES" + fs + "EMHIRESPV_TSh_CF_Country_19862015.txt"

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

    if tech in ["WindOn", "WindOff"]:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["potential"] + region + "_" + tech + "_" + hubheight
    elif tech in ["PV"]:
        if "orientation" in param["PV"]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = "0"
        PathTemp = paths["potential"] + region + "_" + tech + "_" + orientation
    elif tech in ["CSP"]:
        orientation = "0"
        PathTemp = paths["potential"] + region + "_" + tech + "_" + orientation

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

    if tech in ["WindOn", "WindOff"]:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["regional_analysis"] + subregions + "_" + tech + "_" + hubheight
    elif tech in ["PV"]:
        if "orientation" in param["PV"]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = "0"
        PathTemp = paths["regional_analysis"] + subregions + "_" + tech + "_" + orientation
    elif tech in ["CSP"]:
        orientation = "0"
        PathTemp = paths["regional_analysis"] + subregions + "_" + tech + "_" + orientation

    paths[tech]["Locations"] = PathTemp + "_Locations.shp"
    paths[tech]["TS"] = PathTemp + "_TS_" + year + ".csv"
    paths[tech]["Region_Stats"] = PathTemp + "_Region_stats_" + year + ".csv"
    paths[tech]["Sorted_FLH"] = PathTemp + "_sorted_FLH_sampled_" + year + ".mat"

    paths[tech]["Regression_coefficients"] = paths["regression_out"] + subregions + "_" + tech + "_reg_coefficients_"
    paths[tech]["Regression_TS"] = paths["regression_out"] + subregions + "_" + tech + "_reg_TimeSeries_"

    return paths


def discrete_output_paths(paths, param, tech):

    region = param["region_name"]
    year = str(param["year"])

    if tech in ["WindOn", "WindOff"]:
        hubheight = str(param[tech]["technical"]["hub_height"])
        PathTemp = paths["discrete_analysis"] + region + "_" + tech + "_" + hubheight
    elif tech in ["PV"]:
        if "orientation" in param["PV"]["technical"].keys():
            orientation = str(param[tech]["technical"]["orientation"])
        else:
            orientation = "0"
        PathTemp = paths["discrete_analysis"] + region + "_" + tech + "_" + orientation
    elif tech in ["CSP"]:
        orientation = "0"
        PathTemp = paths["discrete_analysis"] + region + "_" + tech + "_" + orientation

    paths[tech]["TS_discrete"] = PathTemp + "_TS_" + year + ".csv"

    return paths
