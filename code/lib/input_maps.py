from .correction_functions import clean_weather_data
from .spatial_functions import *


def generate_maps_for_scope(paths, param):
    """
    This function calls the individual functions that generate the maps for the geographic scope.
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict
    
    :return: The maps are saved directly in the desired paths.
    :rtype: None
    """
    
    generate_weather_files(paths, param)  # MERRA Weather data
    generate_landsea(paths, param)  # Land and Sea
    generate_subregions(paths, param)  # Subregions
    generate_area(paths, param)  # Area Gradient
    generate_landuse(paths, param)  # Landuse
    generate_bathymetry(paths, param)  # Bathymetry
    generate_topography(paths, param)  # Topography
    generate_slope(paths, param)  # Slope
    generate_population(paths, param)  # Population
    generate_protected_areas(paths, param)  # Protected areas
    generate_buffered_population(paths, param)  # Buffered Population


def generate_weather_files(paths, param):
    """
    This function reads the daily NetCDF data (from MERRA-2) for SWGDN, SWTDN, T2M, U50m, and V50m,
    and saves them in matrices with yearly time series with low spatial resolution. Depending on the *MERRA_correction*
    parameter this function will also call clean_weather_data() to remove data outliers.
    This function has to be run only once.

    :param paths: Dictionary including the paths to the MERRA-2 input files *MERRA_IN*, and to the desired output locations for *T2M*, *W50M* and *CLEARNESS*.
    :type paths: dict
    :param param: Dictionary including the year, the spatial scope, and the MERRA_correction parameter.
    :type param: dict

    :return: The files T2M.mat, W50M.mat, and CLEARNESS.mat are saved directly in the defined paths, along with their metadata in JSON files.
    :rtype: None
    """
    timecheck("Start")
    start = datetime.date(param["year"], 1, 1)
    end = datetime.date(param["year"], 12, 31)

    SWGDN = np.array([])
    SWTDN = np.array([])
    T2M = np.array([])
    U50M = np.array([])
    V50M = np.array([])
    status = 0
    delta = (end - start).days + 1
    for date in pd.date_range(start, end):
        # Show status bar
        status = status + 1
        sys.stdout.write("\r")
        sys.stdout.write("Reading NetCDF files " + "[%-50s] %d%%" % ("=" * ((status * 50) // delta), (status * 100) // delta))
        sys.stdout.flush()

        tomorrow = date + pd.Timedelta("1 day")
        if date.day == 29 and date.month == 2:
            continue

        # Name and path of the NetCDF file to be read
        name = paths["MERRA_IN"] + "MERRA2_400.tavg1_2d_rad_Nx." + date.strftime("%Y%m%d") + ".SUB.nc"
        name2 = paths["MERRA_IN"] + "MERRA2_400.tavg1_2d_slv_Nx." + date.strftime("%Y%m%d") + ".SUB.nc"

        # Read NetCDF file, extract hourly tables
        with h5netcdf.File(name, "r") as f:
            # [time, lat 361, lon 576]
            swgdn = np.transpose(subset(f["SWGDN"], param), [1, 2, 0])
            if SWGDN.size == 0:
                SWGDN = swgdn
            else:
                SWGDN = np.concatenate((SWGDN, swgdn), axis=2)

            swtdn = np.transpose(subset(f["SWTDN"], param), [1, 2, 0])
            if SWTDN.size == 0:
                SWTDN = swtdn
            else:
                SWTDN = np.concatenate((SWTDN, swtdn), axis=2)

        with h5netcdf.File(name2, "r") as f:
            t2m = np.transpose(subset(f["T2M"], param), [1, 2, 0])
            if T2M.size == 0:
                T2M = t2m
            else:
                T2M = np.concatenate((T2M, t2m), axis=2)

            u50m = np.transpose(subset(f["U50M"], param), [1, 2, 0])
            if U50M.size == 0:
                U50M = u50m
            else:
                U50M = np.concatenate((U50M, u50m), axis=2)

            v50m = np.transpose(subset(f["V50M"], param), [1, 2, 0])
            if V50M.size == 0:
                V50M = v50m
            else:
                V50M = np.concatenate((V50M, v50m), axis=2)
        if date.year != tomorrow.year:
            # Create the overall wind speed
            W50M = abs(U50M + (1j * V50M))
            # Calculate the clearness index
            CLEARNESS = np.zeros(SWGDN.shape)
            CLEARNESS = np.divide(SWGDN, SWTDN, where=SWTDN != 0)

            sys.stdout.write("\n")
            timecheck("Writing Files: T2M, W50M, CLEARNESS")
            hdf5storage.writes({"T2M": T2M}, paths["T2M"], store_python_metadata=True, matlab_compatible=True)
            hdf5storage.writes({"W50M": W50M}, paths["W50M"], store_python_metadata=True, matlab_compatible=True)
            hdf5storage.writes({"CLEARNESS": CLEARNESS}, paths["CLEARNESS"], store_python_metadata=True, matlab_compatible=True)

            if param["MERRA_correction"]:
                clean_weather_data(paths, param)

            create_json(
                paths["W50M"],
                param,
                ["MERRA_coverage", "region_name", "Crd_all", "res_weather", "MERRA_correction", "MERRA_correction_factor"],
                paths,
                ["MERRA_IN", "W50M"],
            )
            create_json(
                paths["T2M"],
                param,
                ["MERRA_coverage", "region_name", "Crd_all", "res_weather", "MERRA_correction", "MERRA_correction_factor"],
                paths,
                ["MERRA_IN", "T2M"],
            )
            create_json(
                paths["CLEARNESS"],
                param,
                ["MERRA_coverage", "region_name", "Crd_all", "res_weather", "MERRA_correction", "MERRA_correction_factor"],
                paths,
                ["MERRA_IN", "CLEARNESS"],
            )
    timecheck("End")


def generate_landsea(paths, param):
    """
    This function reads the shapefiles of the countries (land areas) and of the exclusive economic zones (sea areas)
    within the scope, and creates two rasters out of them.

    :param paths: Dictionary including the paths *LAND* and *EEZ*.
    :type paths: dict
    :param param: Dictionary including the geodataframes of the shapefiles, the number of features, the coordinates of the bounding box of the spatial scope, and the number of rows and columns.
    :type param: dict

    :return: The tif files for *LAND* and *EEZ* are saved in their respective paths, along with their metadata in JSON files.
    :rtype: None
    """
    m_high = param["m_high"]
    n_high = param["n_high"]
    Crd_all = param["Crd_all"]
    res_desired = param["res_desired"]
    GeoRef = param["GeoRef"]
    nRegions_land = param["nRegions_land"]
    nRegions_sea = param["nRegions_sea"]

    timecheck("Start")
    timecheck("Start Land")
    # Extract land areas
    countries_shp = param["regions_land"]
    Crd_regions_land = param["Crd_regions"][:nRegions_land]
    Ind = ind_merra(Crd_regions_land, Crd_all, res_desired)
    A_land = np.zeros((m_high, n_high))
    status = 0
    for reg in range(0, param["nRegions_land"]):
        # Show status bar
        status = status + 1
        sys.stdout.write("\r")
        sys.stdout.write("Creating A_land " + "[%-50s] %d%%" % ("=" * ((status * 50) // nRegions_land), (status * 100) // nRegions_land))
        sys.stdout.flush()

        # Calculate A_region
        try:
            A_region = calc_region(countries_shp.iloc[reg], Crd_regions_land[reg, :], res_desired, GeoRef)

            # Include A_region in A_land
            A_land[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] = (
                A_land[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] + A_region
            )
        except:
            pass

    # Saving file
    array2raster(paths["LAND"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_land)
    print("\nfiles saved: " + paths["LAND"])
    create_json(
        paths["LAND"], param, ["region_name", "m_high", "n_high", "Crd_all", "res_desired", "GeoRef", "nRegions_land"], paths, ["Countries", "LAND"]
    )
    timecheck("Finish Land")

    timecheck("Start Sea")
    # Extract sea areas
    eez_shp = param["regions_sea"]
    Crd_regions_sea = param["Crd_regions"][-nRegions_sea:]
    Ind = ind_merra(Crd_regions_sea, Crd_all, res_desired)
    A_sea = np.zeros((m_high, n_high))
    status = 0
    for reg in range(0, param["nRegions_sea"]):
        # Show status bar
        status = status + 1
        sys.stdout.write("\r")
        sys.stdout.write(
            "Creating A_sea " + "[%-50s] %d%%" % ("=" * ((status * 50) // param["nRegions_sea"]), (status * 100) // param["nRegions_sea"])
        )
        sys.stdout.flush()

        # Calculate A_region
        A_region = calc_region(eez_shp.iloc[reg], Crd_regions_sea[reg, :], res_desired, GeoRef)

        # Include A_region in A_sea
        A_sea[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] = (
            A_sea[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] + A_region
        )

    # Fixing pixels on the borders to avoid duplicates
    A_sea[A_sea > 0] = 1
    A_sea[A_land > 0] = 0
    # Saving file
    array2raster(paths["EEZ"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_sea)
    print("\nfiles saved: " + paths["EEZ"])
    create_json(
        paths["EEZ"], param, ["region_name", "m_high", "n_high", "Crd_all", "res_desired", "GeoRef", "nRegions_sea"], paths, ["EEZ_global", "EEZ"]
    )
    timecheck("Finish Sea")


def generate_subregions(paths, param):
    """
    This function reads the shapefile of the subregions within the scope, and creates a raster out of it.

    :param paths: Dictionary including the paths *SUB*, *LAND*, *EEZ*.
    :type paths: dict
    :param param: Dictionary including the geodataframe of the shapefile, the number of features, the coordinates of the bounding box of the spatial scope, and the number of rows and columns.
    :type param: dict

    :return: The tif file for *SUB* is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    m_high = param["m_high"]
    n_high = param["n_high"]
    Crd_all = param["Crd_all"]
    res_desired = param["res_desired"]
    GeoRef = param["GeoRef"]
    nRegions_sub = param["nRegions_sub"]

    timecheck("Start Subregions")
    # Read shapefile of regions
    regions_shp = param["regions_sub"]
    Crd_regions_sub = param["Crd_subregions"]
    Ind = ind_merra(Crd_regions_sub, Crd_all, res_desired)
    A_sub = np.zeros((m_high, n_high))
    status = 0
    for reg in range(0, nRegions_sub):
        # Show status bar
        status = status + 1
        sys.stdout.write("\r")
        sys.stdout.write("Creating A_subregions " + "[%-50s] %d%%" % ("=" * ((status * 50) // nRegions_sub), (status * 100) // nRegions_sub))
        sys.stdout.flush()

        # Calculate A_region
        A_region = calc_region(regions_shp.iloc[reg], Crd_regions_sub[reg, :], res_desired, GeoRef)

        # Include A_region in A_sub
        A_sub[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] = (
            A_sub[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] + A_region
        )

    # Fixing pixels on the borders
    with rasterio.open(paths["EEZ"]) as src:
        A_sea = np.flipud(src.read(1)).astype(int)
    with rasterio.open(paths["LAND"]) as src:
        A_land = np.flipud(src.read(1)).astype(int)
    A_sub = A_sub * (A_land + A_sea)

    # Saving file
    array2raster(paths["SUB"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_sub)
    print("\nfiles saved: " + paths["SUB"])
    create_json(
        paths["SUB"], param, ["subregions_name", "m_high", "n_high", "Crd_all", "res_desired", "GeoRef", "nRegions_sea"], paths, ["subregions", "SUB"]
    )
    timecheck("Finish Subregions")

    timecheck("End")


def generate_landuse(paths, param):
    """
    This function reads the global map of land use, and creates a raster out of it for the desired scope.
    There are 17 discrete possible values from 0 to 16, corresponding to different land use classes.
    See :mod:`config.py` for more information on the land use map.

    :param paths: Dictionary including the paths to the global land use raster *LU_global* and to the output path *LU*.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict

    :return: The tif file for *LU* is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, param["res_landuse"])[0]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["LU_global"]) as src:
        w = src.read(1, window=windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]), slice(Ind[3] - 1, Ind[1])))
        w = np.flipud(w)
    w = adjust_resolution(w, param["res_landuse"], param["res_desired"], "category")
    array2raster(paths["LU"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], w)
    print("files saved: " + paths["LU"])
    create_json(paths["LU"], param, ["region_name", "Crd_all", "res_landuse", "res_desired", "GeoRef"], paths, ["LU_global", "LU"])
    timecheck("End")


def generate_bathymetry(paths, param):
    """
    This function reads the global map of bathymetry, resizes it, and creates a raster out of it for the desired scope.
    The values are in meter (negative in the sea).

    :param paths: Dictionary including the paths to the global bathymetry raster *Bathym_global* and to the output path *BATH*.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict

    :return: The tif file for *BATH* is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, param["res_desired"])[0]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["Bathym_global"]) as src:
        A_BATH = src.read(1)
    #A_BATH = resizem(A_BATH, 180 * 240, 360 * 240)
    A_BATH = adjust_resolution(A_BATH, param["res_bathymetry"], param["res_desired"], "mean")
    A_BATH = np.flipud(A_BATH[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]])
    
    array2raster(paths["BATH"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_BATH)
    create_json(paths["BATH"], param, ["region_name", "Crd_all", "res_bathymetry", "res_desired", "GeoRef"], paths, ["Bathym_global", "BATH"])
    print("files saved: " + paths["BATH"])
    timecheck("End")


def generate_topography(paths, param):
    """
    This function reads the tiles that make the global map of topography, picks those that lie completely or partially in the scope,
    and creates a raster out of them for the desired scope. The values are in meter.

    :param paths: Dictionary including the paths to the tiles of the global topography raster *Topo_tiles* and to the output path *TOPO*.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict

    :return: The tif file for *TOPO* is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, param["res_topography"])[0]
    GeoRef = param["GeoRef"]
    Topo = np.zeros((int(180 / param["res_topography"][0]), int(360 / param["res_topography"][1])))
    tile_extents = np.zeros((24, 4), dtype=int)
    i = 1
    j = 1
    for letter in char_range("A", "X"):
        north = (i - 1) * 45 / param["res_topography"][0] + 1
        east = j * 60 / param["res_topography"][1]
        south = i * 45 / param["res_topography"][0]
        west = (j - 1) * 60 / param["res_topography"][1] + 1
        tile_extents[ord(letter) - ord("A"), :] = [north, east, south, west]
        j = j + 1
        if j == 7:
            i = i + 1
            j = 1
    n_min = (Ind[0] // (45 * 240)) * 45 / param["res_topography"][0] + 1
    e_max = (Ind[1] // (60 * 240) + 1) * 60 / param["res_topography"][1]
    s_max = (Ind[2] // (45 * 240) + 1) * 45 / param["res_topography"][0]
    w_min = (Ind[3] // (60 * 240)) * 60 / param["res_topography"][1] + 1

    need = np.logical_and(
        (np.logical_and((tile_extents[:, 0] >= n_min), (tile_extents[:, 1] <= e_max))),
        np.logical_and((tile_extents[:, 2] <= s_max), (tile_extents[:, 3] >= w_min)),
    )

    status = 0
    for letter in char_range("A", "X"):
        index = ord(letter) - ord("A")
        if need[index]:
            # Show status bar
            status = status + 1
            sys.stdout.write("\r")
            sys.stdout.write(
                "Generating topography map from tiles " + "[%-50s] %d%%" % ("=" * ((status * 50) // sum(need)), (status * 100) // sum(need))
            )
            sys.stdout.flush()

            with rasterio.open(paths["Topo_tiles"] + "15-" + letter + ".tif") as src:
                tile = src.read()
            Topo[tile_extents[index, 0] - 1 : tile_extents[index, 2], tile_extents[index, 3] - 1 : tile_extents[index, 1]] = tile[0, 0:-1, 0:-1]

    A_TOPO = np.flipud(Topo[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]])
    A_TOPO = adjust_resolution(A_TOPO, param["res_topography"], param["res_desired"], "mean")
    array2raster(paths["TOPO"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_TOPO)
    print("\nfiles saved: " + paths["TOPO"])
    create_json(paths["TOPO"], param, ["region_name", "Crd_all", "res_topography", "res_desired", "GeoRef"], paths, ["Topo_tiles", "TOPO"])
    timecheck("End")


def generate_slope(paths, param):
    """
    This function reads the topography raster for the scope, and creates a raster of slope out of it. The slope is calculated in
    percentage, although this can be changed easily at the end of the code.

    :param paths: Dictionary including the paths to the topography map of the scope *TOPO* and to the output path *SLOPE*.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict

    :return: The tif file for SLOPE is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, res_desired)[0]
    GeoRef = param["GeoRef"]
    Lat1 = np.arange(-90, 90, res_desired[0])
    Lat2 = np.arange(-90 + res_desired[0], 90 + res_desired[0], res_desired[0])
    latMid = (Lat1 + Lat2) / 2
    deltaLat = abs(Lat1 - Lat2)

    Lat1 = np.arange(-90, 90, res_desired[0])
    Lat2 = np.arange(-90 + res_desired[0], 90 + res_desired[0], res_desired[0])
    latMid_2 = (Lat1 + Lat2) / 2

    Lon1 = np.arange(-180, 180, res_desired[1])
    Lon2 = np.arange(-180 + res_desired[1], 180 + res_desired[1], res_desired[1])
    deltaLon = abs(Lon1 - Lon2)

    m_per_deg_lat = 111132.954 - 559.822 * cos(np.deg2rad(2 * latMid)) + 1.175 * cos(np.deg2rad(4 * latMid))
    m_per_deg_lon = (np.pi / 180) * 6367449 * cos(np.deg2rad(latMid_2))

    x_cell = repmat(deltaLon, int(180 / res_desired[1]), 1) * repmat(m_per_deg_lon, int(360 / res_desired[1]), 1).T
    x_cell = x_cell[Ind[0] - 1 : Ind[2], Ind[3] -1: Ind[1]]
    x_cell = np.flipud(x_cell)

    y_cell = repmat((deltaLat * m_per_deg_lat), int(360 / res_desired[0]), 1).T
    y_cell = y_cell[Ind[0] - 1 : Ind[2], Ind[3] -1: Ind[1]]
    y_cell = np.flipud(y_cell)

    with rasterio.open(paths["TOPO"]) as src:
        A_TOPO = src.read(1)

    kernel = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]) / 8
    dzdx = convolve(A_TOPO, kernel) / x_cell
    kernel = np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]]) / 8
    dzdy = convolve(A_TOPO, kernel) / y_cell

    slope_deg = arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / np.pi
    slope_pc = tan(np.deg2rad(slope_deg)) * 100

    A_SLP = np.flipud(slope_pc)
    array2raster(paths["SLOPE"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_SLP)
    print("files saved: " + paths["SLOPE"])
    create_json(paths["SLOPE"], param, ["region_name", "Crd_all", "res_topography", "res_desired", "GeoRef"], paths, ["TOPO", "SLOPE"])
    timecheck("End")


def generate_population(paths, param):
    """
    This function reads the global map of population density, resizes it, and creates a raster out of it for the desired scope.
    The values are in population per pixel.

    :param paths: Dictionary including the paths to the global population raster *Pop_global* and to the output path *POP*.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict

    :return: The tif file for *POP* is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, res_desired)[0]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["Pop_global"]) as src:
        A_POP_part = src.read(1)  # map is only between latitudes -60 and 85
    A_POP = np.zeros((21600, 43200))
    A_POP[600:18000, :] = A_POP_part
    A_POP = adjust_resolution(A_POP, param["res_population"], param["res_desired"], "sum")
    #A_POP = resizem(A_POP, 180 * 240, 360 * 240) / 4  # density is divided by 4
    A_POP = np.flipud(A_POP[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]])
    array2raster(paths["POP"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_POP)
    print("\nfiles saved: " + paths["POP"])
    create_json(paths["POP"], param, ["region_name", "Crd_all", "res_population", "res_desired", "GeoRef"], paths, ["Pop_global", "POP"])
    timecheck("End")


def generate_protected_areas(paths, param):
    """
    This function reads the shapefile of the globally protected areas, adds an attribute whose values are based on the dictionary 
    of conversion (protected_areas) to identify the protection category, then converts the shapefile into a raster for the scope.
    The values are integers from 0 to 10.

    :param paths: Dictionary including the paths to the shapefile of the globally protected areas, to the landuse raster of the scope, and to the output path PA.
    :type paths: dict
    :param param: Dictionary including the dictionary of conversion of protection categories (protected_areas).
    :type param: dict
    :return: The tif file for PA is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """

    timecheck("Start")
    protected_areas = param["protected_areas"]
    # set up protected areas dictionary
    protection_type = dict(zip(protected_areas["IUCN_Category"], protected_areas["type"]))

    # First we will open our raster image, to understand how we will want to rasterize our vector
    raster_ds = gdal.Open(paths["LU"], gdal.GA_ReadOnly)

    # Fetch number of rows and columns
    ncol = raster_ds.RasterXSize
    nrow = raster_ds.RasterYSize

    # Fetch projection and extent
    proj = raster_ds.GetProjectionRef()
    ext = raster_ds.GetGeoTransform()

    raster_ds = None
    shp_path = paths["Protected"]
    # Open the dataset from the file
    dataset = ogr.Open(shp_path, 1)
    layer = dataset.GetLayerByIndex(0)

    # Add a new field
    if not field_exists("Raster", shp_path):
        new_field = ogr.FieldDefn("Raster", ogr.OFTInteger)
        layer.CreateField(new_field)

        for feat in layer:
            pt = feat.GetField("IUCN_CAT")
            feat.SetField("Raster", protection_type[pt])
            layer.SetFeature(feat)
            feat = None

    # Create a second (modified) layer
    outdriver = ogr.GetDriverByName("MEMORY")
    source = outdriver.CreateDataSource("memData")

    # Create the raster dataset
    memory_driver = gdal.GetDriverByName("GTiff")
    out_raster_ds = memory_driver.Create(paths["PA"], ncol, nrow, 1, gdal.GDT_Byte)

    # Set the ROI image's projection and extent to our input raster's projection and extent
    out_raster_ds.SetProjection(proj)
    out_raster_ds.SetGeoTransform(ext)

    # Fill our output band with the 0 blank, no class label, value
    b = out_raster_ds.GetRasterBand(1)
    b.Fill(0)

    # Rasterize the shapefile layer to our new dataset
    gdal.RasterizeLayer(
        out_raster_ds,  # output to our new dataset
        [1],  # output to our new dataset's first band
        layer,  # rasterize this layer
        None,
        None,  # don't worry about transformations since we're in same projection
        [0],  # burn value 0
        [
            "ALL_TOUCHED=FALSE",  # rasterize all pixels touched by polygons
            "ATTRIBUTE=Raster",
        ],  # put raster values according to the 'Raster' field values
    )
    create_json(paths["PA"], param, ["region_name", "protected_areas", "Crd_all", "res_desired", "GeoRef"], paths, ["Protected", "PA"])

    # Close dataset
    out_raster_ds = None
    print("files saved: " + paths["PA"])
    timecheck("End")


def generate_buffered_population(paths, param):
    """
    This function reads the land use raster, identifies urban areas, and excludes pixels around them based on a
    user-defined buffer *buffer_pixel_amount*. It creates a masking raster of boolean values (0 or 1) for the scope.
    Zero means the pixel is excluded, one means it is suitable.
    The function is useful in case there is a policy to exclude renewable energy projects next to urban settlements.

    :param paths: Dictionary including the path to the land use raster for the scope, and to the output path BUFFER.
    :type paths: dict
    :param param: Dictionary including the user-defined buffer (buffer_pixel_amount), the urban type within the land use map (type_urban), and the georeference dictionary.
    :type param: dict

    :return: The tif file for BUFFER is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    buffer_pixel_amount = param["WindOn"]["mask"]["buffer_pixel_amount"]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
    A_lu = np.flipud(A_lu).astype(int)
    A_lu = A_lu == param["landuse"]["type_urban"]  # Land use type for Urban and built-up
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_lu_buffered = generic_filter(A_lu, np.nanmean, footprint=kernel, mode="constant", cval=np.NaN)
    A_notPopulated = (~A_lu_buffered).astype(int)

    array2raster(paths["BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notPopulated)
    print("files saved: " + paths["BUFFER"])
    create_json(paths["BUFFER"], param, ["region_name", "landuse", "WindOn", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "BUFFER"])
    timecheck("End")

    
def generate_area(paths, param):
    """
    This function retreives the coordinates of the spatial scope and computes the pixel area gradient of the corresponding
    raster.

    :param paths: Dictionary of dictionaries containing the path to the output file.
    :type paths: dict
    :param param: Dictionary of dictionaries containing spatial scope coordinates and desired resolution.
    :type param: dict

    :return: The mat file for AREA is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    Crd_all = param["Crd_all"]
    n_high = param["n_high"]
    res_desired = param["res_desired"]

    # Calculate available area
    # WSG84 ellipsoid constants
    a = 6378137  # major axis
    b = 6356752.3142  # minor axis
    e = np.sqrt(1 - (b / a) ** 2)

    # Lower pixel latitudes
    lat_vec = np.arange(Crd_all[2], Crd_all[0], res_desired[0])
    lat_vec = lat_vec[np.newaxis]

    # Lower slice areas
    # Areas between the equator and the lower pixel latitudes circling the globe
    f_lower = np.deg2rad(lat_vec)
    zm_lower = 1 - (e * sin(f_lower))
    zp_lower = 1 + (e * sin(f_lower))

    lowerSliceAreas = np.pi * b ** 2 * ((2 * np.arctanh(e * sin(f_lower))) / (2 * e) + (sin(f_lower) / (zp_lower * zm_lower)))

    # Upper slice areas
    # Areas between the equator and the upper pixel latitudes circling the globe
    f_upper = np.deg2rad(lat_vec + res_desired[0])

    zm_upper = 1 - (e * sin(f_upper))
    zp_upper = 1 + (e * sin(f_upper))

    upperSliceAreas = np.pi * b ** 2 * ((2 * np.arctanh((e * sin(f_upper)))) / (2 * e) + (sin(f_upper) / (zp_upper * zm_upper)))

    # Pixel areas
    # Finding the latitudinal pixel-sized globe slice areas then dividing them by the longitudinal pixel size
    area_vec = ((upperSliceAreas - lowerSliceAreas) * res_desired[1] / 360).T
    A_area = np.tile(area_vec, (1, n_high))

    # Save to HDF File
    hdf5storage.writes({"A_area": A_area}, paths["AREA"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths["AREA"])
    create_json(paths["AREA"], param, ["Crd_all", "res_desired", "n_high"], paths, [])

    timecheck("End")
