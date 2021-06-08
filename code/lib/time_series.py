from .spatial_functions import *
from .physical_models import calc_CF_solar, calc_CF_wind
from .potential import get_merra_raster_data, redistribution_array


def find_representative_locations(paths, param, tech):
    """
    This function reads the masked FLH raster and finds the coordinates and indices of the pixels for the user-defined quantiles for each region.
    It creates a shapefile containing the position of those points for each region, and two MAT files with their
    coordinates and indices.

    :param paths: Dictionary of dictionaries containing path values for FLH MAT files, region statistics, and output paths.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the user-defined quantiles, FLH resolution, and spatial scope.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: The shapefile with the locations and the two MAT files for the coordinates and the indices are saved
             directly in the given paths, along with their corresponding metadata in JSON files.
    :rtype: None
    """
    timecheck("Start")
    FLH_mask = hdf5storage.read("FLH_mask", paths[tech]["FLH_mask"])
    quantiles = param["quantiles"]
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    # Select only indices in the report
    filter = pd.read_csv(paths[tech]["Region_Stats"], sep=";", decimal=",", index_col=0).index
    regions_shp = param["regions_sub"].loc[filter]
    nRegions = len(regions_shp)
    Crd_regions = param["Crd_subregions"]
    Ind = ind_merra(Crd_regions, Crd_all, res_desired)

    reg_ind = np.zeros((nRegions, len(quantiles), 2))
    k = 0
    list_names = []
    list_quantiles = []
    for reg in filter:
        # A_region
        A_region = calc_region(regions_shp.loc[reg], Crd_regions[reg, :], res_desired, GeoRef)

        FLH_reg = A_region * FLH_mask[Ind[reg, 2] - 1 : Ind[reg, 0], Ind[reg, 3] - 1 : Ind[reg, 1]]
        FLH_reg[FLH_reg == 0] = np.nan
        X = FLH_reg.flatten(order="F")
        I_old = np.argsort(X)

        # Escape loop if intersection only yields NaN
        if np.isnan(X).all():
            # do nothing
            continue

        q_rank = 0
        for q in quantiles:
            list_names.append(regions_shp["NAME_SHORT"].loc[reg])
            list_quantiles.append("q" + str(q))
            if q == 100:
                I = I_old[(len(X) - 1) - sum(np.isnan(X).astype(int))]
            elif q == 0:
                I = I_old[0]
            else:
                I = I_old[int(np.round(q / 100 * (len(X) - 1 - sum(np.isnan(X).astype(int)))))]
            # Convert the indices to row-column indices
            I, J = ind2sub(FLH_reg.shape, I)
            reg_ind[k, q_rank, :] = np.array([I + Ind[reg, 2], J + Ind[reg, 3]]).astype(int)
            q_rank = q_rank + 1
        k = k + 1

    reg_ind = np.reshape(reg_ind, (-1, 2), "C").astype(int)
    reg_ind = (reg_ind[:, 0], reg_ind[:, 1])
    param[tech]["Ind_points"] = reg_ind
    param[tech]["Crd_points"] = crd_exact_points(reg_ind, Crd_all, res_desired)
    param[tech]["Crd_points"] = (param[tech]["Crd_points"][0], param[tech]["Crd_points"][1], list_names, list_quantiles)

    # Format point locations
    points = [(param[tech]["Crd_points"][1][i], param[tech]["Crd_points"][0][i]) for i in range(0, len(param[tech]["Crd_points"][0]))]

    # Create shapefile
    schema = {"geometry": "Point", "properties": {"NAME_SHORT": "str", "quantile": "str"}}
    with fiona.open(paths[tech]["Locations"], "w", "ESRI Shapefile", schema) as c:
        c.writerecords(
            [
                {"geometry": mapping(Point(points[i])), "properties": {"NAME_SHORT": list_names[i], "quantile": list_quantiles[i]}}
                for i in range(0, len(points))
            ]
        )
    hdf5storage.writes(
        {"Ind_points": param[tech]["Ind_points"]}, paths[tech]["Locations"][:-4] + "_Ind.mat", store_python_metadata=True, matlab_compatible=True
    )
    hdf5storage.writes(
        {"Crd_points": param[tech]["Crd_points"]}, paths[tech]["Locations"][:-4] + "_Crd.mat", store_python_metadata=True, matlab_compatible=True
    )
    create_json(
        paths[tech]["Locations"],
        param,
        ["author", "comment", tech, "region_name", "subregions_name", "quantiles", "Crd_all"],
        paths,
        ["spatial_scope", "subregions"],
    )
    print("files saved: " + paths[tech]["Locations"])
    timecheck("End")


def generate_time_series_for_representative_locations(paths, param, tech):
    """
    This function generates yearly capacity factor time-series for the technology of choice at quantile locations
    generated in find_locations_quantiles.
    The timeseries are saved in CSV files.

    :param paths: Dictionary of dictionaries containing paths to coordinate and indices of the quantile locations.
    :type paths: dict
    :param param: Dictionary of dictionaries containing processing parameters.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: The CSV file with the time series for all subregions and quantiles is saved directly in the given path,
             along with the corresponding metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    nproc = param["nproc"]
    CPU_limit = np.full((1, nproc), param["CPU_limit"])
    param[tech]["Crd_points"] = hdf5storage.read("Crd_points", paths[tech]["Locations"][:-4] + "_Crd.mat")
    param[tech]["Ind_points"] = hdf5storage.read("Ind_points", paths[tech]["Locations"][:-4] + "_Ind.mat")
    list_names = param[tech]["Crd_points"][2]
    list_quantiles = param[tech]["Crd_points"][3]
    if tech in ["PV", "CSP"]:
        res_weather = param["res_weather"]
        Crd_all = param["Crd_all"]
        Ind = ind_merra(Crd_all, Crd_all, res_weather)[0]

    list_hours = np.arange(0, 8760)
    param["status_bar_limit"] = list_hours[-1]

    # Obtain weather and correction matrices
    param["Ind_nz"] = param[tech]["Ind_points"]
    merraData, rasterData = get_merra_raster_data(paths, param, tech)

    if tech in ["PV", "CSP"]:

        day_filter = np.nonzero(merraData["CLEARNESS"][Ind[2] - 1 : Ind[0], Ind[3] - 1 : Ind[1], :].sum(axis=(0, 1)))
        list_hours = np.arange(0, 8760)
        if nproc == 1:
            param["status_bar_limit"] = list_hours[-1]
            results = calc_TS_solar(list_hours[day_filter], [param, tech, rasterData, merraData])
        else:
            list_hours = np.array_split(list_hours[day_filter], nproc)
            param["status_bar_limit"] = list_hours[0][-1]
            results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
                calc_TS_solar, product(list_hours, [[param, tech, rasterData, merraData]])
            )

    elif tech in ["WindOn", "WindOff"]:

        list_hours = np.array_split(np.arange(0, 8760), nproc)
        param["status_bar_limit"] = list_hours[0][-1]
        results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
            calc_TS_wind, product(list_hours, [[param, tech, rasterData, merraData]])
        )
    print("\n")

    # Collecting results
    TS = np.zeros((len(param[tech]["Ind_points"][0]), 8760))
    if nproc > 1:
        for p in range(len(results)):
            TS = TS + results[p]
    else:
        TS = results

    # Restructuring results
    tuples = list(zip(list_names, list_quantiles))
    column_names = pd.MultiIndex.from_tuples(tuples, names=["NAME_SHORT", "Quantile"])
    results = pd.DataFrame(TS.transpose(), columns=column_names)
    results.to_csv(paths[tech]["TS"], sep=";", decimal=",")
    create_json(
        paths[tech]["TS"],
        param,
        ["author", "comment", tech, "quantiles", "region_name", "subregions_name", "year", "Crd_all"],
        paths,
        [tech, "spatial_scope", "subregions"],
    )
    print("files saved: " + paths[tech]["TS"])
    timecheck("End")


def generate_time_series_for_specific_locations(paths, param, tech):


    """
    This function generates yearly capacity factor time-series for the technology of choice at user defined locations.
    The timeseries are saved in CSV files.

    :param paths: Dictionary of dictionaries containing paths output desired locations.
    :type paths: dict
    :param param: Dictionary of dictionaries containing processing parameters, and user-defined locations.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: The CSV file with the time series for all subregions and quantiles is saved directly in the given path,
             along with the corresponding metadata in a JSON file.
    :rtype: None
    :raise Point locations not found: Is raised when the dictionary containing the points names and locations is empty.
    :raise Points outside spatial scope: Some points are not located inside of the spatial scope, therefore no input maps are available for the calculations
    """
    timecheck("Start")

    nproc = param["nproc"]
    CPU_limit = np.full((1, nproc), param["CPU_limit"])
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    m_high = param["m_high"]
    
    ind_gwa_merra_x = hdf5storage.read("IND_GWA_MERRA_X", paths["IND_GWA_MERRA_X"])
    ind_gwa_merra_y = hdf5storage.read("IND_GWA_MERRA_Y", paths["IND_GWA_MERRA_Y"])
    ind_gwa_merra_x = np.flipud(ind_gwa_merra_x)
    ind_gwa_merra_y = np.flipud(ind_gwa_merra_y)

    # Read user defined locations dictionary
    if not param["useloc"]:
        warn(
            "Point locations not found: Please fill in the name and locations of the points in config.py prior to executing this function",
            UserWarning,
        )
        timecheck("End")
        return
    points_df = pd.DataFrame.from_dict(param["useloc"], orient="index", columns=["lat", "lon"])

    # Filter points outside spatial scope
    lat_max, lon_max, lat_min, lon_min = param["spatial_scope"][0]
    # Points outside the scope bound
    out_scope_df = points_df.loc[
        (lat_min > points_df["lat"]) | (lat_max < points_df["lat"]) | (lon_min > points_df["lon"]) | (lon_max < points_df["lon"])
    ].copy()
    if not out_scope_df.empty:
        out_points = list(out_scope_df.index)
        print("WARNING: The following points are located outside of the spatial scope " + str(param["spatial_scope"][0]) + ": \n" + str(out_scope_df))
        warn("Points located outside spatial scope", UserWarning)
    # Points inside the scope bounds
    points_df = points_df.loc[
        (lat_min <= points_df["lat"]) & (lat_max >= points_df["lat"]) & (lon_min <= points_df["lon"]) & (lon_max >= points_df["lon"])
    ].copy()
    if not points_df.empty:
        # Prepare input for calc_TS functions
        crd = (points_df["lat"].to_numpy(), points_df["lon"].to_numpy())
        ind = ind_exact_points(crd, Crd_all, res_desired)
        list_names = ["UD"] * len(crd[0])
        list_points = list(points_df.index)

        param[tech]["Crd_points"] = (crd[0], crd[1], list_names, list_points)
        param[tech]["Ind_points"] = ind_exact_points(crd, Crd_all, res_desired)
        param[tech]["Ind_merra_points"] = ind
        
        print (ind[0])
        print (ind[1])
     
        for p in range(len(crd[0])):
            param[tech]["Ind_merra_points"][0][p] = (m_high-ind[0][p]-1)/200
            param[tech]["Ind_merra_points"][1][p] = (ind[1][p]+1)/250

        print (param[tech]["Ind_merra_points"][0])
        print (param[tech]["Ind_merra_points"][1])

        # Obtain weather and correction matrices
        param["Ind_nz"] = param[tech]["Ind_points"]
        merraData, rasterData = get_merra_raster_data(paths, param, tech)
        
        for p in range(len(crd[0])):
            print (rasterData["A_cf"][ind[0][p],ind[1][p]])

        if tech in ["PV", "CSP"]:
            # Set up day_filter
            res_weather = param["res_weather"]
            Ind = ind_merra(Crd_all, Crd_all, res_weather)[0]
            day_filter = np.nonzero(merraData["CLEARNESS"][Ind[2] - 1 : Ind[0], Ind[3] - 1 : Ind[1], :].sum(axis=(0, 1)))

            list_hours = np.arange(0, 8760)
            if nproc == 1:
                param["status_bar_limit"] = list_hours[-1]
                results = calc_TS_solar(list_hours[day_filter], [param, tech, rasterData, merraData])
            else:
                list_hours = np.array_split(list_hours[day_filter], nproc)
                param["status_bar_limit"] = list_hours[0][-1]
                results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
                    calc_TS_solar, product(list_hours, [[param, tech, rasterData, merraData]])
                )

        elif tech in ["WindOn", "WindOff"]:
            TS = calc_TS_wind(param, paths, tech, rasterData, merraData["W50M"])

            #list_hours = np.array_split(np.arange(0, 8760), nproc)
            #param["status_bar_limit"] = list_hours[0][-1]
            #results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
                #calc_TS_wind, product(list_hours, [[param, tech, rasterData, merraData]])
            #)
        print("\n")

        # Collecting results
        #TS = np.zeros((len(param[tech]["Ind_points"][0]), 8760))
        #if nproc > 1:
            #for p in range(len(results)):
                #TS = TS + results[p]
        #else:
            #TS = results

        # Restructuring results
        results = pd.DataFrame(TS.transpose(), columns=list_points).rename_axis("Points", axis="columns")
        results.to_csv(paths[tech]["TS_discrete"], sep=";", decimal=",")
        create_json(
            paths[tech]["TS_discrete"],
            param,
            ["author", "comment", tech, "useloc", "region_name", "subregions_name", "year", "Crd_all"],
            paths,
            [tech, "spatial_scope", "subregions"],
        )
        print("files saved: " + paths[tech]["TS_discrete"])
    timecheck("End")


def generate_time_series_for_full_scope(paths, param, tech):
    
    # Obtain weather and correction matrices
    merraData, rasterData = get_merra_raster_data(paths, param, tech)
    
    if tech in ["WindOn", "WindOff"]:
        TS = calc_TS_wind(param, paths, tech, rasterData, merraData["W50M"])
        #TS = np.flipud(TS)
        #hdf5storage.writes({"FLH": FLH}, paths[tech]["FLH"], store_python_metadata=True, matlab_compatible=True)

        # Restructuring results
        results = pd.DataFrame(TS.transpose()).rename_axis("Points", axis="columns")
        results.to_csv(paths[tech]["TS_discrete"], sep=";", decimal=",")
        create_json(
            paths[tech]["TS_discrete"],
            param,
            ["author", "comment", tech, "useloc", "region_name", "subregions_name", "year", "Crd_all"],
            paths,
            [tech, "spatial_scope", "subregions"],
        )
        print("files saved: " + paths[tech]["TS_discrete"])
 
 
def calc_TS_solar(hours, args):
    """
    This function computes the hourly PV and CSP capacity factor for the desired quantiles.

    :param hours: Hour ranks of the year (from 0 to 8759).
    :type hours: numpy array
    :param args: List of arguments:
    
        * *param* (dict): Dictionary including multiple parameters such as the status bar limit, the name of the region, and
          others for calculating the hourly capacity factors.
        * *tech* (str): Name of the technology.
        * *rasterData* (dict): Dictionary of numpy arrays containing land use types, Ross coefficients, albedo coefficients,
          and wind speed correction for every point in *reg_ind*.
        * *merraData* (dict): Dictionary of numpy arrays containing the weather data for every point in *reg_ind*.

    :type args: list

    :return TS: Array of time series for the desired quantiles for each subregion.
    :rtype: numpy array
    """
    # Decomposing the list args
    param = args[0]
    tech = args[1]
    rasterData = args[2]
    merraData = args[3]
    reg_ind = param[tech]["Ind_points"]

    TS = np.zeros((len(reg_ind[0]), 8760))
    status = 0
    for hour in hours:
        if hour <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            display_progress(tech + " " + param["subregions_name"] + " ", (len(hours), status))

        if tech == "PV":
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech)[0]
        elif tech == "CSP":
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech)[1]

        # Aggregates CF to obtain the time series
        CF[np.isnan(CF)] = 0
        TS[:, hour] = CF
    return TS


def calc_TS_wind(param, paths, tech, rasterData, merraData):
    """
    This function computes the hourly onshore and offshore wind capacity factor for desired quantiles.

    :param hours: Hour ranks of the year (from 0 to 8759).
    :type hours: numpy array
    :param args: List of arguments:
    
        * *param* (dict): Dictionary including multiple parameters such as the status bar limit, the name of the region, and
          others for calculating the hourly capacity factors.
        * *tech* (str): Name of the technology.
        * *rasterData* (dict): Dictionary of numpy arrays containing the wind speed correction for every point in *reg_ind*.
        * *merraData* (dict): Dictionary of numpy arrays containing the weather data for every point in *reg_ind*.

    :type args: list

    :return TS: Array of time series for the desired quantiles for each subregion.
    :rtype: numpy array
    """
    # Decomposing the tuple args
    m_high = param["m_high"]
    n_high = param["n_high"]
    m_low = param["m_low"]
    n_low = param["n_low"]
    #reg_ind = param["Ind_nz"]
    res_weather = param["res_weather"]
    res_desired = param["res_desired"]
    
    turbine = param[tech]["technical"]
    Ind_merra_points = param[tech]["Ind_merra_points"]
    Ind_points = param[tech]["Ind_points"]

    TS = np.zeros((len(Ind_points[0]),8760))
    
    #No of pixels from gwa for each pixel from merra
    Num_pix = (res_weather[0] * res_weather[1]) / (res_desired[0] * res_desired[1])

    b_xmin = hdf5storage.read("MERRA_XMIN", paths["MERRA_XMIN"])
    b_xmax = hdf5storage.read("MERRA_XMAX", paths["MERRA_XMAX"])
    b_ymin = hdf5storage.read("MERRA_YMIN", paths["MERRA_YMIN"])
    b_ymax = hdf5storage.read("MERRA_YMAX", paths["MERRA_YMAX"])
    
    merraData = merraData[::-1,:,:]


    for p in range(len(Ind_points[0])):
        #for j in range(n_low):
        print(str(Ind_merra_points[0][p])+"_"+str(Ind_merra_points[1][p]))
        i = Ind_merra_points[0][p]
        j = Ind_merra_points[1][p]
        
        reMerra = redistribution_array(param, paths, merraData[i,j,:],i,j,b_xmin[i,j],b_xmax[i,j],b_ymin[i,j],b_ymax[i,j],Num_pix)
        reRaster = np.flipud(rasterData["A_cf"])
        if np.sum(reMerra):
            for hour in range(8760):
                # Calculate hourly capacity factor
                CF = calc_CF_wind(hour, turbine, reMerra, reRaster[i*200:((i+1)*200),j*250:((j+1)*250)])
                #CF[np.isnan(CF)] = 0
                TS[p, hour] = CF[(m_high-Ind_points[0][p])%200,Ind_points[1][p]%250]
    return TS


def combinations_for_time_series(paths, param, tech):
    """
    This function reads the list of generated regression coefficients for different hub heights and orientations,
    compares it to the user-defined modes and combos and returns a list of lists containing all the available
    combinations. The function will return a warning if the user input and the available time series are not congruent.

    :param paths: Dictionary of dictionaries containing the paths to the regression output folder.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the year, the user defined combos, and subregions name.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return combinations: List of combinations of settings to be used in stratified time series.
    :return inputfiles: List of regression outputs to be used in generating the stratified time series.
    :rtype: tuple (list, list)
    :raise No coefficients: If regression coefficients are not available, a warning is raised.
    :raise Missing coefficients: If regression coefficients are missing based on user-defined combos and mode, a warning is raised.
    """
    subregions = param["subregions_name"]
    year = str(param["year"])

    bef_setting = paths["regression_out"] + subregions + "_" + tech + "_reg_coefficients_"
    aft_setting = "_" + year + ".csv"

    # Reads the files present in input folder
    inputfiles = glob(paths["regression_out"] + subregions + "_" + tech + "_reg_coefficients*" + year + ".csv")

    # Case 1: no files existing
    if len(inputfiles) == 0:
        warn("Run the regression first, before creating stratified time series!", UserWarning)
        return

    # Get existing settings
    settings_existing = []
    for filename in inputfiles:
        list_settings = filename.replace(bef_setting, "").replace(aft_setting, "").split("_")
        settings_existing = settings_existing + [sorted([int(x) for x in list_settings], reverse=True)]
    settings_sorted = set(sorted(map(tuple, settings_existing)))
    print("\nFor technology " + tech + ", regression coefficients for the following combinations have been detected: ", settings_existing)

    # Get required settings
    combinations = param["combo"][tech].values()
    combinations_sorted = []
    for combi in combinations:
        if combi == []:
            combi = sorted(set([item for sublist in settings_sorted for item in sublist]), reverse=True)
        combinations_sorted = combinations_sorted + [tuple(sorted(combi, reverse=True))]
    combinations = set(sorted(combinations_sorted))

    # Case 2: some files are missing
    if not combinations.issubset(settings_sorted):
        print("\nFor technology " + tech + ", regression coefficients for the following combinations are required: ", combinations)
        warn(
            "Not all regression coefficients are available! Generate the missing regression coefficients first, then create stratified time series.",
            UserWarning,
        )
        return

    # Create inputfiles list matching combinations
    combifiles = []
    for combi in combinations:
        combifiles = combifiles + [bef_setting + "_".join(map(str, list(combi))) + aft_setting]

    return list(combinations), combifiles


def generate_time_series_for_regions(paths, param, tech):
    """
    This function reads the coefficients obtained from the regression function as well as the generated time series for
    the combinations of hub heights / orientations and quantiles, to combine them according to user-defined
    *modes* (quantile combination) and *combos* (hub heights / orientation combinations) and saves the results (time series) 
    in a CSV file.

    :param paths: Dictionary of dictionaries containing the paths to the regression coefficients and the time series.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the list of subregions, the modes, and the combos.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: The stratified time series for each region, mode, and combo are saved directly in the given path, along with the metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    modes = param["modes"]
    subregions = param["subregions_name"]
    year = str(param["year"])

    try:
        combinations, inputfiles = combinations_for_time_series(paths, param, tech)
    except UserWarning:
        timecheck("End")
        return

    # Display the combinations of settings to be used
    if tech in ["WindOn", "WindOff"]:
        print("Combinations of hub heights to be used for the stratified time series: ", combinations)
    elif tech in ["PV"]:
        print("Orientations to be used for the stratified time series: ", combinations)

    for tag, combo in param["combo"][tech].items():
        if combo == []:
            combo = sorted(set([item for sublist in combinations for item in sublist]))
        ind = combinations.index(tuple(sorted(combo, reverse=True)))
        coef = pd.read_csv(inputfiles[ind], sep=";", decimal=",", index_col=[0])

        # Extract names of regions
        regions = sorted(list(set([col.split("_")[0] for col in coef.columns])))

        # Load the TS files
        TS_files = {}
        for setting in combo:
            setting_path = paths["regional_analysis"] + subregions + "_" + tech + "_" + str(setting) + "_TS_" + year + ".csv"
            TS_files[setting] = pd.read_csv(setting_path, sep=";", decimal=",", header=[0, 1], index_col=[0])
        quantiles_existing = list(map(int, [s.strip("q") for s in list(TS_files[list(TS_files.keys())[0]].columns.levels[1])]))

        # Loop over modes and regions
        TS_df = pd.DataFrame(index=range(8760), dtype="float16")
        for mode_tag, quantiles in modes.items():
            # Check if quantiles are available
            if not set(quantiles).issubset(set(quantiles_existing)):
                warn("\nSet quantiles " + str(quantiles) + " do not match available quantiles from input files: " + str(quantiles_existing))
                timecheck("End")
                return

            for reg in regions:
                col_name = reg + "_" + tech + "_" + tag + "_" + mode_tag
                TS_df[col_name] = np.zeros((8760, 1))
                filter_reg = [col for col in coef if col.startswith(reg)]
                for setting in combo:
                    sum_quantiles = coef.loc[quantiles, filter_reg].sum().sum()
                    for quantile in quantiles:
                        if sum_quantiles:
                            TS_df[col_name] = (
                                TS_df[col_name]
                                + TS_files[setting][reg, "q" + str(quantile)] * coef.loc[quantile, reg + "_" + str(setting)] / sum_quantiles
                            )
                        else:
                            TS_df[col_name] = TS_df[col_name] + TS_files[setting][reg, "q" + str(quantile)] / len(quantiles) / len(combo)

        st = ""
        for setting in combo:
            st = st + str(setting) + "_"
        param["st"] = st
        TS_df.to_csv(paths[tech]["Regression_TS"] + st + year + ".csv", sep=";", decimal=",")
        create_json(
            paths[tech]["Regression_TS"] + st + year + ".csv",
            param,
            ["author", "comment", tech, "quantiles", "modes", "combo", "region_name", "subregions_name", "year"],
            paths,
            ["spatial_scope", "subregions"],
        )
        print("File Saved: " + paths[tech]["Regression_TS"] + st + year + ".csv")
    timecheck("End")
