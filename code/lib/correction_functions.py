from .spatial_functions import calc_region, array2raster
from .util import *


def clean_weather_data(paths, param):
    """
    This function detects data outliers in the weather input .mat files. An outlier is a data point, for which
    the absolute value of the difference between the yearly average value and the mean of the direct neighbors
    (Moore neighborhood) is higher than a user-defined threshold *MERRA_correction_factor*. It replaces the hourly values
    with the hourly values of the mean of the neighbors, and overwrites the original .mat file.

    :param paths: Dictionary including the path to the weather .mat files.
    :type paths: dict
    :param param: Dictionary including the threshold value *MERRA_correction_factor*.
    :type param: dict

    :return: The file weather .mat files are overwritten after the correction.
    :rtype: None
    """
    timecheck("Start")
    for p in ["W50M", "CLEARNESS", "T2M"]:

        # Read Weather Data
        weather = hdf5storage.read(p, paths[p])
        mean = np.mean(weather, 2)

        # Set convolution mask
        kernel = np.ones((3, 3))
        kernel[1, 1] = 0

        # Compute average Convolution
        neighbors = generic_filter(mean, np.nanmean, footprint=kernel, mode="constant", cval=np.NaN)
        ratio = mean / neighbors

        # Extract over threshold Points
        points = np.where(abs(ratio - np.mean(ratio)) > param["MERRA_correction_factor"][p])

        # Correct points hourly
        for t in range(weather.shape[2]):
            weather[points[0], points[1], t] = weather[points[0], points[1], t] / ratio[points[0], points[1]]

        # Save corrected Wind
        hdf5storage.writes({p: weather}, paths[p], store_python_metadata=True, matlab_compatible=True)
    timecheck("End")


def generate_wind_correction(paths, param):
    """
    This function creates a matrix of correction factors for onshore and/or offshore wind.
    There are different types of correction:

    * Gradient correction: Adjusts for the hub height of the wind turbines, based on the Hellmann coefficients of each land use type.
      This correction applies always.
    * Resolution correction: Performs a redistribution of wind speed when increasing the resolution based on land use types, while ensuring that
      the average of each MERRA-2 cell at 50m is still the same. This correction is optional, and is activated if *res_correction* is 1.
      If not activated, the same value from the low resolution is repeated.
    * Topographic/Orographic correction: Takes into account the elevation of the terrain, because MERRA-2 usually underestimates
      the wind speed in mountains. This correction is optional, uses data from the Global Wind Atlas for all countries in the scope,
      and is activated only for onshore wind if *topo_correction* is 1

    :param paths: Dictionary of dictionaries containing the paths to the land, land use, and topography rasters, and to the output files CORR_ON and CORR_OFF.
    :type paths: dict
    :param param: Dictionary of dictionaries containing user-preferences regarding the wind correction, landuse, hub height, weather and desired resolutions.
    :type param: dict

    :return: The rasters for wind correction CORR_ON and/or CORR_OFF are saved directly in the user-defined paths, along with their metadata in JSON files.
    :rtype: None
    """
    timecheck("Start")
    res_correction_on = param["WindOn"]["resource"]["res_correction"]
    res_correction_off = param["WindOff"]["resource"]["res_correction"]
    topo_correction = param["WindOn"]["resource"]["topo_correction"]
    GeoRef = param["GeoRef"]
    landuse = param["landuse"]
    with rasterio.open(paths["LU"]) as src:
        A_lu = np.flipud(src.read(1)).astype(int)
    A_hellmann = changem(A_lu, landuse["hellmann"], landuse["type"]).astype(float)

    # Onshore resolution correction
    if "WindOn" in param["technology"]:
        turbine_height_on = param["WindOn"]["technical"]["hub_height"]

        if res_correction_on:
            m_low = param["m_low"]
            n_low = param["n_low"]
            m_high = param["m_high"]
            n_high = param["n_high"]
            res_weather = param["res_weather"]
            res_desired = param["res_desired"]
            A_gradient_height = changem(A_lu.astype(float), landuse["height"], landuse["type"])
            Sigma = sumnorm_MERRA2((50 / A_gradient_height) ** A_hellmann, m_low, n_low, res_weather, res_desired)
            A_cf_on = ((turbine_height_on / 50) * turbine_height_on / A_gradient_height) ** A_hellmann / resizem(Sigma, m_high, n_high)
            del A_gradient_height, Sigma
        else:
            A_cf_on = (turbine_height_on / 50) ** A_hellmann

        # Topographic correction (only onshore)
        with rasterio.open(paths["LAND"]) as src:
            A_land = np.flipud(src.read(1)).astype(int)
        A_cf_on = A_cf_on * A_land
        del A_land
        if topo_correction:
            if not os.path.isfile(paths["CORR_GWA"]):
                calc_gwa_correction(paths, param)
            gwa_correction = hdf5storage.read("correction_" + param["WindOn"]["resource"]["topo_weight"], paths["CORR_GWA"])
            A_cf_on = A_cf_on * gwa_correction
        array2raster(paths["CORR_ON"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_cf_on)
        create_json(paths["CORR_ON"], param, ["region_name", "year", "WindOn", "landuse", "res_weather", "res_desired"], paths, ["LAND", "CORR_GWA"])
        print("\nfiles saved: " + paths["CORR_ON"])

    # Offshore resolution correction
    if "WindOff" in param["technology"]:
        turbine_height_off = param["WindOff"]["technical"]["hub_height"]

        if res_correction_off:
            m_low = param["m_low"]
            n_low = param["n_low"]
            m_high = param["m_high"]
            n_high = param["n_high"]
            res_weather = param["res_weather"]
            res_desired = param["res_desired"]
            A_gradient_height = changem(A_lu.astype(float), landuse["height"], landuse["type"])
            Sigma = sumnorm_MERRA2((50 / A_gradient_height) ** A_hellmann, m_low, n_low, res_weather, res_desired)
            A_cf_off = ((turbine_height_off / 50) * turbine_height_off / A_gradient_height) ** A_hellmann / resizem(Sigma, m_high, n_high)
            del A_gradient_height, Sigma
        else:
            A_cf_off = (turbine_height_off / 50) ** A_hellmann
        del A_hellmann
        with rasterio.open(paths["EEZ"]) as src:
            A_eez = np.flipud(src.read(1)).astype(int)
        A_cf_off = A_cf_off * A_eez

        array2raster(paths["CORR_OFF"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_cf_off)
        create_json(paths["CORR_ON"], param, ["region_name", "year", "WindOff", "landuse", "res_weather", "res_desired"], paths, ["CORR_GWA"])
        print("\nfiles saved: " + paths["CORR_OFF"])
    timecheck("End")


def calc_gwa_correction(paths, param):
    """
    This function creates a correction matrix for onshore wind based on the topography and the frequency distribution of wind speed in each country in the Global
    Wind Atlas.
    We first read the MERRA-2 data for wind and increase its resolution without editing it. We also read the topographic data for the whole scope.
    For each country, we filter the two datasets based on valid pixels and obtain *w50m_reg* and *topo_reg*. We correct the wind data based on the following formula:
    
        :math:`w50m_{corrected} = w50m_{reg} * min(exp(ai * topo_{reg} + bi), 3.5)`
        
    where *ai* and *bi* are two parameters that have to be determined, so that the error (difference to the sorted frequencies of wind speeds from the GWA) is minimalized
    for the whole scope. Instead of using a nonlinear optimization, we simply iterate over a range of discrete possibilities for *ai* and *bi*, save the errors, then pick
    the combinations that minimize the error.
    It is possible to weight the error of each country based on its area or on its installed onshore wind capacity, or to give the same weight to all the countries.
    Finally, the three possible correction matrices are saved.

    :param paths: Dictionary that contains the paths to wind speed data, topography, IRENA summary, Global Wind Atlass folder, and to the output location.
    :type paths: dict
    :param param: Dictionary that contains assumptions about the desired resolution, size of the output, shapefile of countries, and georeference dictionary.
    :type param: dict

    :return: The correction matrices are saved in the same MAT file, directly in the given path, along with the metadata in the corresponding JSON file.
    :rtype: None
    """
    m_high = param["m_high"]
    n_high = param["n_high"]
    res_desired = param["res_desired"]
    nCountries = param["nRegions_land"]
    countries_shp = param["regions_land"]
    Crd_countries = param["Crd_regions"][0:nCountries, :]
    GeoRef = param["GeoRef"]

    # Obtain wind speed at 50m
    W50M = hdf5storage.read("W50M", paths["W50M"])
    W50M = np.mean(W50M, 2)
    W50M = resizem(W50M, m_high, n_high)

    # Obtain topography
    with rasterio.open(paths["TOPO"]) as src:
        w = src.read(1)
    TOPO = np.flipud(w)

    # Clean IRENA data and filter them for desired scope
    if not os.path.isfile(paths["IRENA_summary"]):
        clean_IRENA_summary(paths, param)

    # Get the installed capacities
    inst_cap = pd.read_csv(paths["IRENA_summary"], sep=";", decimal=",", index_col=0, usecols=[0, 1, 2])
    inst_cap = inst_cap.loc[inst_cap["Technology"] == "Onshore wind energy"]

    w_size = np.zeros((nCountries, 1))
    w_cap = np.zeros((nCountries, 1))
    # Try different combinations of (a, b)
    combi_list = list(product(np.arange(0.00046, 0.00066, 0.00002), np.arange(-0.3, 0, 0.025)))
    errors = np.zeros((len(combi_list), nCountries))
    status = 0
    for reg in range(0, nCountries):
        # Show status bar
        status = status + 1
        display_progress("Finding wind correction factors", (nCountries, status))

        A_region = calc_region(countries_shp.iloc[reg], Crd_countries[reg, :], res_desired, GeoRef)
        reg_name = countries_shp.iloc[reg]["GID_0"]
        Ind_reg = np.nonzero(A_region)
        w_size[reg] = len(Ind_reg[0])
        try:
            w_cap[reg] = inst_cap.loc[reg_name, "inst-cap (MW)"]
        except KeyError:
            w_cap[reg] = 0

        # Load MERRA data, increase its resolution, and fit it to the extent
        w50m_reg = W50M[Ind_reg]
        topo_reg = TOPO[Ind_reg]

        # Get the sampled frequencies from the GWA
        try:
            try:
                w50m_gwa = pd.read_csv(paths["GWA"][:-14] + reg_name + paths["GWA"][-14:], usecols=["gwa_ws"]).to_numpy()[:, 0]
            except:
                w50m_gwa = pd.read_csv(paths["GWA"][:-14] + reg_name + paths["GWA"][-14:], usecols=["val"]).to_numpy()[:, 0]
        except:
            w50m_gwa = pd.read_csv(paths["GWA"][:-14] + reg_name + paths["GWA"][-14:], usecols=[0]).to_numpy()[:, 0]

        i = 0
        for combi in combi_list:
            ai, bi = combi
            w50m_corrected = w50m_reg * np.minimum(np.exp(ai * topo_reg + bi), 3.5)
            w50m_sorted = np.sort(w50m_corrected)
            w50m_sampled = np.flipud(w50m_sorted[:: (len(w50m_sorted) // len(w50m_gwa) + 1)])
            if len(w50m_sampled) != len(w50m_gwa):
                len_diff = len(w50m_gwa) - len(w50m_sampled)
                w50m_sampled = np.append(w50m_sampled, w50m_sorted[:len_diff])
            try:
                w50m_diff = w50m_sampled - w50m_gwa
                errors[i, reg] = np.sqrt((w50m_diff ** 2).sum())
            except ValueError:
                errors[i, reg] = 0
            i = i + 1

    w_size = np.tile(w_size / w_size.sum(), (1, len(combi_list))).transpose()
    w_cap = np.tile(w_cap / w_cap.sum(), (1, len(combi_list))).transpose()

    ae, be = combi_list[np.argmin(np.sum(errors / nCountries, 1))]
    correction_none = np.zeros(TOPO.shape)
    correction_none = np.minimum(np.exp(ae * TOPO + be), 3.5)

    a_size, b_size = combi_list[np.argmin(np.sum(errors * w_size, 1))]
    correction_size = np.zeros(TOPO.shape)
    correction_size = np.minimum(np.exp(a_size * TOPO + b_size), 3.5)

    a_cap, b_cap = combi_list[np.argmin(np.sum(errors * w_cap, 1))]
    correction_capacity = np.zeros(TOPO.shape)
    correction_capacity = np.minimum(np.exp(a_cap * TOPO + b_cap), 3.5)

    hdf5storage.writes(
        {"correction_none": correction_none, "correction_size": correction_size, "correction_capacity": correction_capacity},
        paths["CORR_GWA"],
        store_python_metadata=True,
        matlab_compatible=True,
    )
    create_json(
        paths["CORR_GWA"],
        param,
        ["author", "comment", "region_name", "subregions_name", "year", "Crd_all", "res_desired", "GeoRef"],
        paths,
        ["W50M", "TOPO", "IRENA_summary"],
    )
    return


def clean_IRENA_summary(paths, param):
    """
    This function reads the IRENA database, format the output for selected regions and computes the FLH based on the
    installed capacity and yearly energy production. The results are saved in CSV file.

    :param param: Dictionary of dictionaries containing list of subregions, and year.
    :type param: dict
    :param paths: Dictionary of dictionaries containing the paths to the IRENA country name dictionary, and IRENA database.
    :type paths: dict

    :return: The CSV file containing the summary of IRENA data for the countries within the scope is saved directly in the desired path, along with the corresponding metadata in a JSON file.
    :rtype: None
    """
    year = str(param["year"])
    filter_countries = param["regions_land"]["GID_0"].to_list()
    IRENA_dict = pd.read_csv(paths["IRENA_dict"], sep=";", index_col=0)
    IRENA_dict = IRENA_dict["Countries shapefile"].to_dict()
    IRENA = pd.read_csv(paths["IRENA"], skiprows=7, sep=";", index_col=False, usecols=[0, 1, 2, 3])
    for i in IRENA.index:
        if pd.isnull(IRENA.loc[i, "Country/area"]):
            IRENA.loc[i, "Country/area"] = IRENA.loc[i - 1, "Country/area"]
        if pd.isnull(IRENA.loc[i, "Technology"]):
            IRENA.loc[i, "Technology"] = IRENA.loc[i - 1, "Technology"]

    for c in IRENA["Country/area"].unique():
        IRENA.loc[IRENA["Country/area"] == c, "Country/area"] = IRENA_dict[c]

    IRENA = IRENA.set_index(["Country/area", "Technology"])

    IRENA = IRENA.fillna(0).sort_index()

    for (c, t) in IRENA.index.unique():
        sub_df = IRENA.loc[(c, t), :]
        inst_cap = sub_df.loc[sub_df["Indicator"] == "Electricity capacity (MW)", year][0]
        if isinstance(inst_cap, str):
            inst_cap = int(inst_cap.replace(" ", ""))
            IRENA.loc[(IRENA.index.isin([(c, t)])) & (IRENA["Indicator"] == "Electricity capacity (MW)"), year] = inst_cap
        gen_prod = sub_df.loc[sub_df["Indicator"] == "Electricity generation (GWh)", year][0]
        if isinstance(gen_prod, str):
            gen_prod = 1000 * int(gen_prod.replace(" ", ""))
            IRENA.loc[(IRENA.index.isin([(c, t)])) & (IRENA["Indicator"] == "Electricity generation (GWh)"), year] = gen_prod
        if inst_cap == 0:
            FLH = 0
        else:
            FLH = gen_prod / inst_cap
        IRENA = IRENA.append(pd.DataFrame([["FLH (h)", FLH]], index=[(c, t)], columns=["Indicator", year])).sort_index()

    # Filter countries
    IRENA = IRENA.reset_index()
    IRENA = IRENA.set_index(["Country/area"]).sort_index()
    IRENA = IRENA.loc[IRENA.index.isin(filter_countries)]
    # Reshape
    IRENA = IRENA.reset_index()
    IRENA = IRENA.set_index(["Country/area", "Technology"])
    IRENA = IRENA.pivot(columns="Indicator")[year].rename(
        columns={"Electricity capacity (MW)": "inst-cap (MW)", "Electricity generation (GWh)": "prod (MWh)"}
    )
    IRENA = IRENA.astype(float)
    IRENA.to_csv(paths["IRENA_summary"], sep=";", decimal=",", index=True)
    create_json(paths["IRENA_summary"], param, ["author", "comment", "region_name", "year"], paths, ["Countries", "IRENA", "IRENA_dict"])
    print("files saved: " + paths["IRENA_summary"])
