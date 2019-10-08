from util import *


def define_spatial_scope(scope_shp):
    """
    This function reads the spatial scope shapefile and returns its bounding box.

    :param scope_shp: Spatial scope shapefile
    :type scope_shp: Geopandas dataframe

    :return box: list of the bounding box coordinates
    :rtype box: list
    """
    scope_shp = scope_shp.to_crs({"init": "epsg:4326"})
    r = scope_shp.total_bounds
    box = r[::-1][np.newaxis]
    return box


def crd_merra(Crd_regions, res_weather):
    """
    This function calculates coordinates of the bounding box covering MERRA-2 data.

    :param Crd_regions: Coordinates of the bounding boxes of the regions.
    :type Crd_regions: numpy array
    :param res_weather: Weather data resolution.
    :type res_weather: list

    :return Crd: Coordinates of the bounding box covering MERRA-2 data for each region.
    :rtype: numpy array
    """
    Crd = np.array(
        [
            np.ceil((Crd_regions[:, 0] + res_weather[0] / 2) / res_weather[0]) * res_weather[0] - res_weather[0] / 2,
            np.ceil(Crd_regions[:, 1] / res_weather[1]) * res_weather[1],
            np.floor((Crd_regions[:, 2] + res_weather[0] / 2) / res_weather[0]) * res_weather[0] - res_weather[0] / 2,
            np.floor(Crd_regions[:, 3] / res_weather[1]) * res_weather[1],
        ]
    )
    Crd = Crd.T
    return Crd


def crd_exact_points(Ind_points, Crd_all, res):
    """
    This function converts indices of points in high resolution rasters into longitude and latitude coordinates.

    :param Ind_points: Tuple of arrays of indices in the vertical and horizontal axes.
    :type Ind_points: tuple of arrays
    :param Crd_all: Array of coordinates of the bounding box of the spatial scope.
    :type Crd_all: numpy array
    :param res: Data resolution in the vertical and horizontal dimensions.
    :type res: list
    
    :return Crd_points: Coordinates of the points in the vertical and horizontal dimensions.
    :rtype: list of arrays
    """
    Crd_points = [Ind_points[0] * res[0] + Crd_all[2], Ind_points[1] * res[1] + Crd_all[3]]
    return Crd_points


def subset(A, param):
    """
    This function retrieves a subset of the global MERRA-2 coverage based on weather resolution and the
    bounding box coordinates of the spatial scope.

    :param A: Weather data on a global scale.
    :type A: numpy array
    :param param: Dictionary of parameters containing MERRA-2 coverage and the name of the region.
    :type param: dict

    :return subset: The subset of the weather data contained in the bounding box of *spatial_scope*.
    :rtype: numpy array
    """
    if param["MERRA_coverage"] == "World" and param["region_name"] != "World":
        crd = param["Crd_all"]
        res = param["res_weather"]
        southlim = int(math.floor((crd[2] + res[0] / 10 + 90 + res[0] / 2) / res[0]))
        northlim = int(math.ceil((crd[0] - res[0] / 10 + 90 + res[0] / 2) / res[0]))
        westlim = int(math.floor((crd[3] + res[1] / 10 + 180) / res[1]))
        eastlim = int(math.ceil((crd[1] - res[1] / 10 + 180) / res[1]))
        subset = A[:, southlim:northlim, westlim:eastlim]
    else:
        subset = A
    return subset


def ind_merra(Crd, Crd_all, res):
    """
    This function converts longitude and latitude coordinates into indices within the spatial scope of MERRA-2 data.

    :param Crd: Coordinates to be converted into indices.
    :type Crd: numpy array
    :param Crd_all: Coordinates of the bounding box of the spatial scope.
    :type Crd_all: numpy array
    :param res: Resolution of the data, for which the indices are produced.
    :type res: list
    
    :return Ind: Indices within the spatial scope of MERRA-2 data.
    :rtype: 
    """
    if len(Crd.shape) == 1:
        Crd = Crd[np.newaxis]
    Ind = np.array(
        [
            (Crd[:, 0] - Crd_all[2]) / res[0],
            (Crd[:, 1] - Crd_all[3]) / res[1],
            (Crd[:, 2] - Crd_all[2]) / res[0] + 1,
            (Crd[:, 3] - Crd_all[3]) / res[1] + 1,
        ]
    )
    Ind = np.transpose(Ind).astype(int)
    return Ind


def ind_global(Crd, res_desired):
    """
    This function converts longitude and latitude coordinates into indices on a global data scope, where the origin is at (-90, -180).

    :param Crd: Coordinates to be converted into indices.
    :type Crd: numpy array
    :param res_desired: Desired resolution in the vertical and horizontal dimensions.
    :type res_desired: list
    
    :return Ind: Indices on a global data scope.
    :rtype: numpy array
    """
    if len(Crd.shape) == 1:
        Crd = Crd[np.newaxis]
    Ind = np.array(
        [
            np.round((90 - Crd[:, 0]) / res_desired[0]) + 1,
            np.round((180 + Crd[:, 1]) / res_desired[1]),
            np.round((90 - Crd[:, 2]) / res_desired[0]),
            np.round((180 + Crd[:, 3]) / res_desired[1]) + 1,
        ]
    )
    Ind = np.transpose(Ind.astype(int))
    return Ind


def calc_geotiff(Crd_all, res_desired):
    """
    This function returns dictionary containing the Georefferencing parameters for geotiff creation,
    based on the desired extent and resolution.

    :param Crd_all: Coordinates of the bounding box of the spatial scope.
    :type Crd_all: numpy array
    :param res_desired: Desired data resolution in the vertical and horizontal dimensions.
    :type res_desired: list

    :return GeoRef: Georeference dictionary containing *RasterOrigin*, *RasterOrigin_alt*, *pixelWidth*, and *pixelHeight*.
    :rtype: dict
    """
    GeoRef = {
        "RasterOrigin": [Crd_all[3], Crd_all[0]],
        "RasterOrigin_alt": [Crd_all[3], Crd_all[2]],
        "pixelWidth": res_desired[1],
        "pixelHeight": -res_desired[0],
    }
    return GeoRef


def calc_region(region, Crd_reg, res_desired, GeoRef):
    """
    This function reads the region geometry, and returns a masking raster equal to 1 for pixels within and 0 outside of
    the region.

    :param region: Region geometry
    :type region: Geopandas series
    :param Crd_reg: Coordinates of the region
    :type Crd_reg: list
    :param res_desired: Desired high resolution of the output raster
    :type res_desired: list
    :param GeoRef: Georeference dictionary containing *RasterOrigin*, *RasterOrigin_alt*, *pixelWidth*, and *pixelHeight*.
    :type GeoRef: dict

    :return A_region: Masking raster of the region.
    :rtype A_region: numpy array
    """
    latlim = Crd_reg[2] - Crd_reg[0]
    lonlim = Crd_reg[3] - Crd_reg[1]
    M = int(math.fabs(latlim) / res_desired[0])
    N = int(math.fabs(lonlim) / res_desired[1])
    A_region = np.ones((M, N))
    origin = [Crd_reg[3], Crd_reg[2]]

    if region.geometry.geom_type == "MultiPolygon":
        features = [feature for feature in region.geometry]
    else:
        features = [region.geometry]
    west = origin[0]
    south = origin[1]
    profile = {
        "driver": "GTiff",
        "height": M,
        "width": N,
        "count": 1,
        "dtype": rasterio.float64,
        "crs": "EPSG:4326",
        "transform": rasterio.transform.from_origin(west, south, GeoRef["pixelWidth"], GeoRef["pixelHeight"]),
    }

    with MemoryFile() as memfile:
        with memfile.open(**profile) as f:
            f.write(A_region, 1)
            out_image, out_transform = mask.mask(f, features, crop=False, nodata=0, all_touched=False, filled=True)
        A_region = out_image[0]

    return A_region


def clean_IRENA_summary(param, paths):
    """
    This function reads the IRENA database, format the output for selected regions and computes the FLH based on the
    installed capacity and yearly energy production. The results are saved in .csv file.

    :param param: Dictionary of dictionaries containing list of subregions, and year.
    :type param: dict
    :param paths: Dictionary of dictionaries containing the paths to the IRENA country name dictionary, and IRENA database.
    :type paths: dict

    :return: None
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
            sub_df.loc[sub_df["Indicator"] == "Electricity capacity (MW)", year] = inst_cap
        gen_prod = sub_df.loc[sub_df["Indicator"] == "Electricity generation (GWh)", year][0]
        if isinstance(gen_prod, str):
            gen_prod = 1000 * int(gen_prod.replace(" ", ""))
            sub_df.loc[sub_df["Indicator"] == "Electricity generation (GWh)", year] = gen_prod
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
    print("files saved: " + paths["IRENA_summary"])


def clean_FLH_regression(param, paths):
    """
    This function creates a .csv file containing the model FLH used for regression. If the region is present in the
    IRENA database then the FLH are extracted directly from there. In case it is not present, a place holder for the
    regions is written in the csv file and it is the user's responsibility to fill in an appropriate value.
    The function will warn the user, and print all regions that are left blank.

    :param param: Dictionary of dictionaries containing the list of regions.
    :type param: dict
    :param paths: Dictionary of dictionaries containing the paths to IRENA_summary, IRENA_dict.
    :type paths: dict

    :return: list of string of the missing regions.
    :rtype: list of str
    """
    # Read IRENA summary
    if not os.path.isfile(paths["IRENA_summary"]):
        clean_IRENA_summary(param, paths)
    IRENA_summary = pd.read_csv(paths["IRENA_summary"], sep=";", decimal=",", index_col=[0, 1])

    # Load IRENA dictionary
    IRENA_dict = pd.read_csv(paths["IRENA_dict"], sep=";")
    IRENA_dict.dropna(inplace=True)
    IRENA_dict.set_index(["NAME_SHORT"], inplace=True)
    IRENA_dict = IRENA_dict["Countries shapefile"].to_dict()

    # Setup FLH_regression dataframe
    list_regions = param["regions_sub"]["NAME_SHORT"].values.tolist()
    FLH_regression = pd.DataFrame(columns=["WindOn", "WindOff", "PV", "CSP"], index=list_regions, dtype=float)

    # Fill in FLH_regression dataframe
    missing = []
    for country in list_regions:
        # Country/region found in IRENA_dict and IRENA_summary
        if country in IRENA_dict.keys():
            FLH_regression.loc[country] = [
                IRENA_summary.loc[(IRENA_dict[country], "Onshore wind energy"), "FLH (h)"],
                IRENA_summary.loc[(IRENA_dict[country], "Offshore wind energy"), "FLH (h)"],
                IRENA_summary.loc[(IRENA_dict[country], "Solar photovoltaic"), "FLH (h)"],
                IRENA_summary.loc[(IRENA_dict[country], "Concentrated solar power"), "FLH (h)"],
            ]
        # Missing country/region, require user input
        else:
            missing.append(country)

    # Save FLH_regression
    FLH_regression.to_csv(paths["FLH_regression"], sep=";", decimal=",", index=True)
    print("files saved: " + paths["FLH_regression"])

    # Return Missing countries/regions
    warn(
        "The following countries/regions are not present in the IRENA Database: "
        + ",".join(missing)
        + ".\nTheir corresponding FLH have been left blank.",
        UserWarning,
    )
    return missing

  
def clean_TS_regression(param, paths, tech):
    """
    This function creates a .csv file containing the model Time-series used for regression. If the region is present in
    the EMHIRES text files then the TS is extracted directly from it. If the region is not present in the EMHIRES text
    files the highest FLH generated TS is used instead and is scaled to match IRENA FLH.

    :param param: Dictionary containing technologies under study, *FLH_regression* dataframe, list of subregions contained in shapefile, and year.
    :type param: dict
    :param paths: Dictionary containing paths to EMHIRES text files
    :type paths: dict

    :return: None
    """

    # load IRENA FLH data
    irena = pd.read_csv(paths['FLH_regression'], sep=';', decimal=',', index_col=0)
    technologies = param["technology"]
    # Find intersection between desired regions and irena regions
    list_regions = param["regions_sub"]["NAME_SHORT"].values.tolist()
    list_regions = sorted(list(set(list_regions).intersection(set(irena.index))))

    # Create TS_regression dataframe
    TS_regression = pd.DataFrame(index=range(1, 8761), columns=list_regions)

    # Load EMHIRES data for desired year
    if tech in ["PV", "CSP"]:
        date_index = pd.date_range(start="1/1/1986", end="1/1/2016", freq="H", closed="left")
        EMHIRES = pd.read_csv(paths[tech]["EMHIRES"], " ")
        EMHIRES = EMHIRES.set_index(date_index)
        EMHIRES = EMHIRES.loc["1/1/" + str(param["year"]) : "1/1/" + str(param["year"] + 1)]
    else:
        EMHIRES = pd.read_csv(paths[tech]["EMHIRES"], "\t")
        EMHIRES = EMHIRES[EMHIRES["Year"] == param["year"]].reset_index()
        EMHIRES = EMHIRES.drop(["index", "Time step", "Date", "Year", "Month", "Day", "Hour"], axis=1)

    # Find intersection between EMHIRES and list_regions
    intersect_regions = sorted(list((set(list_regions).intersection(set(EMHIRES.columns)))))

    # Load setting combinations
    settings = combinations_for_regression(paths, param, tech)
    if not settings[0]:
        settings = settings[1]
    else:
        settings = settings[0]

    for region in list_regions:
        IRENA_FLH = irena.loc[region, tech]
        # Region is present in both EMHIRES and IRENA
        if region in intersect_regions:
            # Scale EMHIRES TS to IRENA
            TS_regression[region] = (EMHIRES[region] * (IRENA_FLH / sum(EMHIRES[region]))).values
        # Region is not present in EMHIRES, Use scaled generated TS instead
        else:
            # Load Generated TS and Scale it with IRENA FLH
            GenTS = read_generated_TS(paths, param, tech, settings, region)
            # Find highest FLH valued TS
            settings_sorted = np.array(
                pd.DataFrame((np.nansum(GenTS[key]) for key in GenTS.keys()), index=settings, columns=["FLH_all_quant"])
                .sort_values(by="FLH_all_quant", ascending=0)
                .index
            )
            GenTS["TS_Max"] = GenTS[str(settings_sorted[0])]["q" + str(np.max(param["quantiles"]))]
            # Scale max TS to IRENA
            TS_regression[region] = (GenTS["TS_Max"] * (IRENA_FLH / GenTS["TS_Max"].sum())).values


    # Save TS_regression as .csv
    TS_regression.to_csv(paths[tech]["TS_regression"], sep=";", decimal=",", index=True)
    print("files saved: " + paths[tech]["TS_regression"])


def calc_gwa_correction(param, paths):
    """
    Missing description

    :param param:
    :param paths:

    :return:
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
        clean_IRENA_summary(param, paths)

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
        sys.stdout.write("\rFinding wind correction factors " + "[%-50s] %d%%" % ("=" * ((status * 50) // nCountries), (status * 100) // nCountries))
        sys.stdout.flush()

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
        w50m_gwa = pd.read_csv(paths["GWA"][:-14] + reg_name + paths["GWA"][-14:], usecols=["gwa_ws"]).to_numpy()[:, 0]

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
    return


def calc_gcr(Crd_all, m_high, n_high, res_desired, GCR):
    """
    This function creates a GCR weighting matrix for the desired geographic extent.
    The sizing of the PV system is conducted on a user-defined day for a shade-free exposure
    to the sun during a given number of hours.

    :param Crd_all: desired geographic extent of the whole region (north, east, south, west)
    :type Crd_all: list
    :param m_high: number of rows
    :type m_high: int
    :param n_high: number of columns
    :type n_high: int
    :param res_desired: map's high resolution
    :type res_desired: list
    :param GCR: includes the user-defined day and the duration of the shade-free period

    :return: GCR raster
    :rtype: numpy array
    """
    # Vector of latitudes between (south) and (north), with resolution (res_should) degrees
    lat = np.arange((Crd_all[2] + res_desired[0] / 2), Crd_all[0], res_desired[0])[np.newaxis]
    lon = np.arange((Crd_all[3] + res_desired[1] / 2), Crd_all[1], res_desired[1])[np.newaxis]

    # Repeating for all longitudes/latitudes
    lat = repmat(lat.transpose(), 1, n_high)
    lon = repmat(lon, m_high, 1)

    # Solar time where shade-free exposure starts
    omegast = 12 - GCR["shadefree_period"] / 2

    # Calculation
    omega = 15 * (omegast - 12)  # Hour angle
    phi = abs(lat)  # Latitude angle

    beta = np.maximum(phi, 15)  # Tilt angle = latitude, but at least 15 degrees
    # Optimal tilt angle (loosely based on Breyer 2010)
    beta = np.minimum(np.abs(phi), 55)  # The tilt angle is preferably equal to the latitude
    range_lat = np.logical_and(np.abs(phi) >= 35, np.abs(phi) < 65)
    beta[range_lat] = (beta[range_lat] - 35) / 65 * 55 + 35  # Tilt angle does not increase very quickly
    range_lat = np.logical_and(lat >= 35, lat < 65)
    range_lon = np.logical_and(lon >= -20, lon < 30)
    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat, range_lon)] - 35) / 65 * 45 + 35  # Europe
    range_lat = np.logical_and(lat >= 20, lat < 65)
    range_lon = np.logical_and(lon >= 75, lon < 140)
    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat, range_lon)] - 20) / 65 * 60 + 20  # Asia/China

    if Crd_all[2] > 0:
        day = GCR["day_north"]
        # Declination angle
        delta = repmat(arcsind(0.3978) * sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), m_high, 1)

    if Crd_all[0] < 0:
        day = GCR["day_south"]
        # Declination angle
        delta = repmat(arcsind(0.3978) * sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), m_high, 1)

    if (Crd_all[2] * Crd_all[0]) < 0:
        lat_pos = int(np.sum(lat >= 0, axis=0)[0])
        day = GCR["day_north"]
        # Declination angle
        delta_pos = repmat(arcsind(0.3978) * sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), lat_pos, 1)

        lat_neg = int(np.sum(lat < 0, axis=0)[0])
        day = GCR["day_south"]
        # Declination angle
        delta_neg = repmat(arcsind(0.3978) * sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), lat_neg, 1)
        delta = np.append(delta_neg, delta_pos, axis=0)

    # Elevation angle
    alpha = arcsind(sind(delta) * sind(phi) + cosd(delta) * cosd(phi) * cosd(omega))

    # Azimuth angle
    azi = arccosd((sind(delta) * cosd(phi) - cosd(delta) * sind(phi) * cosd(omega)) / cosd(alpha))

    # The GCR
    A_GCR = 1 / (cosd(beta) + np.abs(cosd(azi)) * sind(beta) / tand(alpha))

    # Fix too large and too small values of GCR
    A_GCR[A_GCR < 0.2] = 0.2
    A_GCR[A_GCR > 0.9] = 0.9

    return A_GCR


def sampled_sorting(Raster, sampling):
    """
    This function returns a list with a defined length of sorted values sampled from a numpy array.

    :param Raster: Input raster to be sorted
    :type Raster: numpy array
    :param sampling: Number of values to be sampled from the raster, defines length of outputted list
    :type sampling: int

    :return: List of sorted values sampled from Raster.
    :rtype: List
    """
    # Flatten the raster and sort raster from highest to lowest
    Sorted_FLH = np.sort(Raster.flatten(order="F"))
    Sorted_FLH = np.flipud(Sorted_FLH)

    # Loop over list with sampling increment

    s = Sorted_FLH[0]  # Highest value
    for n in np.arange(sampling, len(Sorted_FLH), sampling):
        s = np.append(s, Sorted_FLH[n])
    s = np.append(s, Sorted_FLH[-1])  # Lowest value

    return s


def read_generated_TS(paths, param, tech, settings, subregion):
    """
    This function returns a dictionary containing the available Time series generated by the script based on
    the desired technology and settings.

    :param paths: Dictionary including output folder for regional_analysis
    :type paths: dict
    :param param: Dictionary including list of subregions and year
    :type param: dict
    :param tech: Technology
    :type tech: str
    :param settings: list of lists containing setting combinations
    :type settigns: list of list
    :param subregion: name of the subregion
    :type subregion: str

    :return GenTS: Dictionary of timeseries indexed by setting and quantile
    :rtype GenTS: dict
    """
    subregions = param["subregions_name"]
    year = str(param["year"])

    bef_setting = paths["regional_analysis"] + subregions + "_" + tech + "_"
    if tech == "CSP":
        bef_setting = paths["regional_analysis"] + subregions + "_" + tech
    aft_setting = "_TS_" + year + ".csv"

    # Setup the data dictionary for generated TS for each quantile
    GenTS = {}

    for setting in settings:
        TS_Temp = pd.read_csv(bef_setting + str(setting) + aft_setting, sep=";", decimal=",", dtype=str)

        filter_reg = [col for col in TS_Temp if col.startswith(subregion)]
        # Remove undesired regions
        TS_Temp = TS_Temp[filter_reg]

        # Exit function if subregion is not present in TS files
        if TS_Temp.empty:
            return None

        TS_Temp.columns = TS_Temp.iloc[0]
        TS_Temp = TS_Temp.drop(0)
        # Replace ',' with '.' for float conversion
        for q in range(0, len(TS_Temp.columns)):
            TS_Temp.iloc[:, q] = TS_Temp.iloc[:, q].apply(lambda x: x.replace(",", ".")).astype(float)
            TS_Temp.columns.name = ""
            TS_Temp.reset_index(inplace=True, drop=True)
        GenTS[str(setting)] = TS_Temp

    return GenTS


def regmodel_load_data(paths, param, tech, settings, subregion):
    """
    This function returns a dictionary used to initialize a pyomo abstract model for the regression analysis
    of each region.

    :param paths: dictionary of dictionaries containing the paths to the Timeseries csv files
    :param param: dictionry of dictionaries contating IRENA's region list, FLHs and EMHIRES model timeseries.
    :param tech: name of the technology under study
    :type tech: str
    :param settings: list of all the settings (hub heights/orientations) to be used in the regression
    :type settings: list
    :param subregion: code name of region
    :type subregion: str

    :return: Dictionary containing regression parameters
    :rtype: dict
    """
    subregions = param["subregions_name"]
    year = str(param["year"])
    FLH = param["FLH_regression"]
    time = range(1, 8761)

    # Setup dataframe for IRENA
    model_FLH = FLH.loc[subregion, tech]
    if np.isnan(model_FLH):
        return None

    # Read data from output folder
    GenTS = read_generated_TS(paths, param, tech, settings, subregion)
    if GenTS is None:
        return None

    # reorder hubheights to go from max TS to min TS:
    settings_sorted = np.array(
        pd.DataFrame((np.nansum(GenTS[key]) for key in GenTS.keys()), index=settings, columns=["FLH_all_quant"])
        .sort_values(by="FLH_all_quant", ascending=0)
        .index
    )
    GenTS["TS_Max"] = np.nansum(GenTS[str(settings_sorted[0])]["q" + str(np.max(param["quantiles"]))])
    GenTS["TS_Min"] = np.nansum(GenTS[str(settings_sorted[-1])]["q" + str(np.min(param["quantiles"]))])

    # Check if solution exists
    solution_check = (GenTS["TS_Max"] > model_FLH, GenTS["TS_Min"] < model_FLH)

    # Prepare Timeseries dictionary indexing by height and quantile

    if solution_check == (False, True):
        Timeseries = GenTS[str(settings_sorted[0])]["q" + str(np.max(param["quantiles"]))]

    elif solution_check == (True, False):
        Timeseries = GenTS[str(settings_sorted[-1])]["q" + str(np.min(param["quantiles"]))]

    elif solution_check == (True, True):
        Timeseries = {}

        for s in settings_sorted:
            for q in param["quantiles"]:
                for t in time:
                    Timeseries[(s, q, t)] = np.array(GenTS[str(s)]["q" + str(q)])[t - 1]

    # Setup dataframe for TS Models
    TS_reg = param["TS_regression"]
    ts = np.array(TS_reg[subregion].values)

    TS = {}
    for t in time:
        TS[(t,)] = ts[t - 1]
    # Create data_input dictionary
    data = {
        None: {
            "s": {None: settings_sorted},
            "q": {None: param["quantiles"]},
            "FLH": {None: model_FLH},
            "shape": TS,
            "t": {None: np.array(time)},
            "TS": Timeseries,
            "IRENA_best_worst": solution_check,
            "GenTS": GenTS,
        }
    }
    return data


def combinations_for_regression(paths, param, tech):
    """
    This function reads the list of generated time-series for different hub-heights and orientations, compares it to the
    user defined combinations and returns a list of lists containing all the available combinations. The function will
    returns a warning if the user input and the available time-series are not congruent.

    :param paths: Dictionary of dictionaries containing the paths to the regional analysis output folder.
    :type paths: dict
    :param param: Dictionary of dictionaries containing a list of sub-regions' name, user defined combinations.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: List of combinations
    :rtype: list of lists
    :raise missing data: If no time series are available for this technology
    :raise missing combination: If a hub-height or orientation is missing based on user defined combinations
    """
    subregions = param["subregions_name"]
    year = str(param["year"])

    # Reads the files present in input folder
    inputfiles = glob(paths["regional_analysis"] + subregions + "_" + tech + "*_TS_" + year + ".csv")

    # Case 1: no files existing
    if len(inputfiles) == 0:
        warn("Generate time series first, before doing the regression!", UserWarning)
        return

    # Get existing settings
    settings_existing = []
    for filename in inputfiles:
        bef_setting = paths["regional_analysis"] + subregions + "_" + tech + "_"
        aft_setting = "_TS_" + year + ".csv"
        settings_existing = settings_existing + [int(filename.replace(bef_setting, "").replace(aft_setting, ""))]
    settings_existing = set(settings_existing)
    print("\nFor technology " + tech + ", time series for the following settings have been detected: ", settings_existing)

    # Get required settings
    combinations = param["regression"][tech].values()
    combinations_sorted = []
    for combi in combinations:
        combinations_sorted = combinations_sorted + [sorted(combi)]
    combinations = combinations_sorted
    settings_required = set([item for sublist in combinations for item in sublist])

    # Case 2: some files are missing
    if not settings_required.issubset(settings_existing):
        print("\nFor technology " + tech + ", time series for the following settings are required: ", settings_required)
        warn("Not all time series are available! Generate the missing time series first, then do the regression.", UserWarning)
        return
    # Case 3: all the required files exist (and nothing more)
    if settings_existing.issubset(settings_required):
        if [] in combinations:
            full_combi = sorted(list(settings_existing))
            if not (full_combi in combinations):
                combinations = combinations + [full_combi]
            combinations.remove([])
    # Case 4: more files exist than those required
    else:
        if [] in combinations:
            full_combi = sorted(list(settings_existing))
            if not (full_combi in combinations):
                combinations = combinations + [full_combi]
            combinations.remove([])
    return combinations


def combinations_for_stratified_timeseries(paths, param, tech):
    """
    This function reads the list of generated regression coefficients for different hub-heights and orientations,
    compares it to the user defined modes and combos and returns a list of lists containing all the available
    combinations. The function will returns a warning if the user input and the available time series are not congruent.

    :param paths: Dictionary of dictionaries containing the paths to the regression output folder.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the year, the user defined combos, and list of sub-regions' name.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: List of combinations
    :rtype: list of lists
    :raise No coefficients: If regression coefficients are not available.
    :raise Missing coefficients: If regression coefficients are missing based on user defined combos and mode.
    """
    subregions = param["subregions_name"]
    year = str(param["year"])

    # Reads the files present in input folder
    inputfiles = glob(paths["regression_out"] + subregions + "_" + tech + "_reg_coefficients*" + year + ".csv")

    # Case 1: no files existing
    if len(inputfiles) == 0:
        warn("Run the regression first, before creating stratified time series!", UserWarning)
        return

    # Get existing settings
    settings_existing = []
    for filename in inputfiles:
        bef_setting = paths["regression_out"] + subregions + "_" + tech + "_reg_coefficients_"
        aft_setting = "_" + year + ".csv"
        list_settings = filename.replace(bef_setting, "").replace(aft_setting, "").split("_")
        settings_existing = settings_existing + [sorted([int(x) for x in list_settings])]
    settings_sorted = sorted(settings_existing)
    print("\nFor technology " + tech + ", regression coefficients for the following combinations have been detected: ", settings_existing)

    # Get required settings
    combinations = param["combo"][tech].values()
    combinations_sorted = []
    for combi in combinations:
        if combi == []:
            combi = list(set([item for sublist in combinations for item in sublist]))
        combinations_sorted = combinations_sorted + [sorted(combi)]
    combinations = sorted(combinations_sorted)

    # Case 2: some files are missing
    if combinations != settings_sorted:
        print("\nFor technology " + tech + ", regression coefficients for the following combinations are required: ", combinations)
        warn(
            "Not all regression coefficients are available! Generate the missing regression coefficients first, then create stratified time series.",
            UserWarning,
        )
        return
    return settings_existing, inputfiles


def get_merra_raster_Data(paths, param, tech):
    """
    Returns tuple of two dictionaries containing weather and correction rasters for specified technology.

    :param paths: dictionary of dictionaries containing the paths to the input weather and raster data.
    :type paths: dict
    :param param: dictionary of dictionaries containing landuse, Ross coefficients, Albedo, and hellman coefficient.
    :type param: dict
    :param tech: Technology under study
    :type tech: str

    :return: tuple of dictionaries for the weather and correction data
    :rtype: tuple (dict, dict)
    """
    landuse = param["landuse"]
    merraData = {}
    rasterData = {}
    # Wind Speed Data
    merraData["W50M"] = hdf5storage.read("W50M", paths["W50M"])
    if tech in ["PV", "CSP"]:

        # Other weather Data
        # Clearness index - stored variable CLEARNESS
        merraData["CLEARNESS"] = hdf5storage.read("CLEARNESS", paths["CLEARNESS"])
        # Temperature 2m above the ground - stored variable T2M
        merraData["T2M"] = hdf5storage.read("T2M", paths["T2M"])

        # Calculate A matrices correction
        # A_lu
        with rasterio.open(paths["LU"]) as src:
            w = src.read(1)
        rasterData["A_lu"] = np.flipud(w)
        # A_Ross (Temperature coefficients for heating losses)
        rasterData["A_Ross"] = changem(rasterData["A_lu"], param["landuse"]["Ross_coeff"], param["landuse"]["type"]).astype("float16")
        # A_albedo (Reflectivity coefficients)
        rasterData["A_albedo"] = changem(rasterData["A_lu"], param["landuse"]["albedo"], param["landuse"]["type"]).astype("float16")
        # A_WS_Coef wind Speed at 2m above the ground
        A_hellmann = changem(rasterData["A_lu"], landuse["hellmann"], landuse["type"])
        rasterData["A_WindSpeed_Corr"] = ((2 / 50) ** A_hellmann).astype("float16")
        del A_hellmann

    elif tech in ["WindOn", "WindOff"]:
        reg_ind = param["Ind_nz"]
        # A_cf
        if tech == "WindOn":
            paths_corr = paths["CORR_ON"]
        else:
            paths_corr = paths["CORR_OFF"]
        with rasterio.open(paths_corr) as src:
            w = src.read(1)
        rasterData["A_cf"] = np.flipud(w).astype("float16")
        rasterData["A_cf"] = rasterData["A_cf"][reg_ind]
        del w
    return merraData, rasterData
