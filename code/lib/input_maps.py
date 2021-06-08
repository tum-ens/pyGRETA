from .correction_functions import clean_weather_data
from .spatial_functions import *
#import cython


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
    #only for windoff #generate_bathymetry(paths, param)  # Bathymetry
    generate_topography(paths, param)  # Topography
    generate_slope(paths, param)  # Slope
    #generate_population(paths, param)  # Population #not used anywhere?
    generate_protected_areas(paths, param)  # Protected areas
    generate_buffered_population(paths, param)  # Buffered Population #Not related to population? Why?
    generate_array_coordinates(paths, param)
    generate_livestock(paths,param)
    generate_biomass_production(paths, param)


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
        name = paths["MERRA_IN"] + "MERRA2_400.tavg1_2d_rad_Nx." + date.strftime("%Y%m%d") + ".nc4.nc4"
        name2 = paths["MERRA_IN"] + "MERRA2_400.tavg1_2d_slv_Nx." + date.strftime("%Y%m%d") + ".nc4.nc4"

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
    #res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, param["res_landuse"])[0]
    GeoRef = param["GeoRef"]
    lu_a = param["WindOn"]["weight"]["lu_availability"]
    with rasterio.open(paths["LU_global"]) as src:
        w = src.read(1, window=windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]), slice(Ind[3] - 1, Ind[1])))
        w = np.flipud(w)
    w = adjust_resolution(w, param["res_landuse"], param["res_desired"], "category")
    #if "WindOn" in param["technology"]:
    w = recalc_lu_resolution(w, param["res_landuse"], param["res_desired"], lu_a)
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
    #res_desired = param["res_desired"]
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
    #res_desired = param["res_desired"]
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
    #if "WindOn" in param["technology"]:
    A_TOPO = recalc_topo_resolution(A_TOPO, param["res_topography"], param["res_desired"])
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

    #x_cell = repmat(deltaLon, int(180 / res_desired[1]), 1) * repmat(m_per_deg_lon, int(360 / res_desired[1]), 1).T
    #x_cell = x_cell[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]]
    x_cell = repmat(deltaLon[Ind[3] - 1 : Ind[1]], Ind[2]-Ind[0]+1, 1) * repmat(m_per_deg_lon[Ind[0] - 1 : Ind[2]], Ind[1]-Ind[3]+1, 1).T
    x_cell = np.flipud(x_cell)

    #y_cell = repmat((deltaLat * m_per_deg_lat), int(360 / res_desired[0]), 1).T
    #y_cell = y_cell[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]]
    y_cell = repmat((deltaLat[Ind[0] - 1 : Ind[2]] * m_per_deg_lat[Ind[0] - 1 : Ind[2]]), Ind[1]-Ind[3]+1, 1).T
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
    Ind = ind_global(Crd_all, param["res_desired"])[0]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["Pop_global"]) as src:
        A_POP_part = src.read(1)  # map is only between latitudes -60 and 85
    A_POP = np.zeros((21600, 43200))
    A_POP[600:18000, :] = A_POP_part
    #A_POP = adjust_resolution(A_POP, param["res_population"], param["res_desired"], "sum")
    #A_POP = resizem(A_POP, 180 * 240, 360 * 240) / 4  # density is divided by 4
    A_POP = np.flipud(A_POP[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]])
    #if "WindOn" in param["technology"]:
    A_POP = recalc_topo_resolution(A_POP, param["res_landuse"], param["res_desired"])
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


def generate_array_coordinates(paths, param):
    
    timecheck("Start")
    
    Crd_all = param["Crd_all"]
    ymax, xmax, ymin, xmin = Crd_all
    res_weather = param["res_weather"]

    W50M = hdf5storage.read("W50M", paths["W50M"])
    w50m_shape = W50M.shape
    
    #bounding box coordinates of each pixel in merra
    b_xmin = np.zeros([w50m_shape[0],w50m_shape[1]])
    b_ymin = np.zeros([w50m_shape[0],w50m_shape[1]])
    b_xmax = np.zeros([w50m_shape[0],w50m_shape[1]])
    b_ymax = np.zeros([w50m_shape[0],w50m_shape[1]])
    for i in range(w50m_shape[0]):
        for j in range(w50m_shape[1]):
            b_xmin[i,j] = xmin + j * res_weather[1]
            b_xmax[i,j] = xmin + (j+1) * res_weather[1]
            b_ymin[i,j] = ymax - (i+1) * res_weather[0]
            b_ymax[i,j] = ymax - i * res_weather[0]
            
    hdf5storage.writes({"MERRA_XMIN": b_xmin}, paths["MERRA_XMIN"], store_python_metadata=True, matlab_compatible=True)
    hdf5storage.writes({"MERRA_XMAX": b_xmax}, paths["MERRA_XMAX"], store_python_metadata=True, matlab_compatible=True)
    hdf5storage.writes({"MERRA_YMIN": b_ymin}, paths["MERRA_YMIN"], store_python_metadata=True, matlab_compatible=True)
    hdf5storage.writes({"MERRA_YMAX": b_ymax}, paths["MERRA_YMAX"], store_python_metadata=True, matlab_compatible=True)
            
    GWA_speed = rasterio.open(paths["GWA_global"])
    GWA_array = GWA_speed.read(1)
    gwa_rows, gwa_cols = GWA_array.shape
    
    #coordinates for center of each pixel in GWA
    x_gwa = np.zeros([gwa_rows,gwa_cols])
    y_gwa = np.zeros([gwa_rows,gwa_cols])
    for k in range(gwa_rows):
        for l in range(gwa_cols):
            x_gwa[k,l] = GWA_speed.xy(k,l,offset='center')[0]
            y_gwa[k,l] = GWA_speed.xy(k,l,offset='center')[1]
    
    hdf5storage.writes({"GWA_X": x_gwa}, paths["GWA_X"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths["GWA_X"])
    hdf5storage.writes({"GWA_Y": y_gwa}, paths["GWA_Y"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths["GWA_Y"])
    timecheck("End")


def generate_livestock(paths, param):
    """
    This function reads the global maps of each livestock density, resizes it, and creates a raster out of it for the desired scope.
    The values are in number of animals per sq.km.

    :param paths: Dictionary including the paths to the global livestock rasters *LS_global* and to the output path *LS*.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict

    :return: The tif files for *LS* is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, param["res_livestock"])[0]
    GeoRef = param["GeoRef"]
    
    A_area = hdf5storage.read("A_area", paths["AREA"])
    
    for animal in param["Biomass"]["livestock"]["animal"]:
        with rasterio.open(paths["LS_global"]+animal+"_2006.tif") as src:
            A_LS = src.read(1, window=windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]), slice(Ind[3] - 1, Ind[1])))
        A_LS = np.flipud(A_LS)
        A_LS = recalc_livestock_resolution(A_LS, param["res_livestock"], param["res_desired"])
        #print (np.size(A_LS))
        A_LS[A_LS<0]=float(0)
        A_LS = np.multiply(A_LS, A_area) / (10 ** 6)
        array2raster(paths["LS"]+animal+".tif", GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_LS)
        print("\nfiles saved: " + paths["LS"]+animal+".tif")
        create_json(paths["LS"]+animal+".tif", param, ["region_name", "Crd_all", "res_livestock", "res_desired", "GeoRef"], paths, ["LS_global", "LS"])
    
    timecheck("End")
    
    
def generate_biomass_production(paths, param):
    
    timecheck("Start")
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    nRegions = param["nRegions_sub"]
    regions_shp = param["regions_sub"]
    countries_shp = param["regions_land"]
    
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
    A_lu = np.flipud(A_lu).astype(int)
    
    with rasterio.open(paths["PA"]) as src:
        A_protect = src.read(1)
    A_protect = np.flipud(A_protect).astype(int)
    A_include = np.zeros(A_lu.shape) 
    val_include = [0,5,6,7,8,9,10]
    for i in val_include:
        A_i = A_protect==i
        A_include = A_include + A_i
    
    A_lu_forest = np.zeros(A_lu.shape)
    val_forest = [1,2,3,4,5]
    for j in val_forest:
        A_j = A_lu == j
        A_lu_forest = A_lu_forest + A_j
        
    A_Bioenergy = np.zeros(A_lu.shape)
    A_Bioco2 = np.zeros(A_lu.shape)
    A_country_area = np.zeros(A_lu.shape)
    
    # for reg in range(0, param["nRegions_land"]):
        # r = calc_region(countries_shp.loc[reg], Crd_all, res_desired, GeoRef)  
        # A_country_area = A_country_area + r
      
    for reg in range(0, nRegions):
        A_region_extended = calc_region(regions_shp.loc[reg], Crd_all, res_desired, GeoRef)
        A_country_area = A_country_area + A_region_extended    
        
        A_lu_crop = A_lu == param["landuse"]["type_croplands"]
        A_lu_veg = A_lu == param["landuse"]["type_vegetation"]
        
        A_lu_crop = np.multiply(A_lu_crop,A_region_extended)
        n_lu_crop = np.sum(A_lu_crop)
        #print (n_lu_crop)
        A_lu_veg = np.multiply(A_lu_veg,A_region_extended)
        n_lu_veg = np.sum(A_lu_veg)
        #print (n_lu_veg)

        if n_lu_crop or n_lu_veg:        
            region_code = regions_shp.loc[reg]["GID_2"][:-2]
            if len(region_code)==11:
                region_code = region_code[0:3]+region_code[4:7]+region_code[8:]
            elif len(region_code) == 10:
                if region_code[7] == ".":
                    region_code = region_code[0:3]+region_code[4:7]+"0"+region_code[8:]
                else:
                    region_code = region_code[0:3]+"0"+region_code[4:6]+region_code[7:]
            elif len(region_code) == 9:
                if region_code[7] == ".":
                    region_code = region_code[0:3]+region_code[4:7]+"00"+region_code[8:]
                elif region_code[6] == ".":
                    region_code = region_code[0:3]+"0"+region_code[4:6]+"0"+region_code[7:]
                else:
                    region_code = region_code[0:3]+"00"+region_code[4]+region_code[6:]
            elif len(region_code) == 8:
                if region_code[6] == ".":
                    region_code = region_code[0:3]+"0"+region_code[4:6]+"00"+region_code[7:]
                else:
                    region_code = region_code[0:3]+"00"+region_code[4]+"0"+region_code[6:]
            else:
                region_code = region_code[0:3]+"00"+region_code[4]+"00"+region_code[6:]
            print (region_code)
            bio_energy = 0
            bio_co2 = 0
            for crop in param["Biomass"]["agriculture"]["crops"]:
                production_crop = pd.read_csv(paths["Biomass"]+crop+".csv", index_col = ["agromap_area_2"],usecols=["agromap_area_2","Production"])
                if region_code in production_crop.index:
                    production_crop = production_crop[production_crop.index==region_code]
                    for residue in param["Biomass"]["agriculture"]["residue"][crop]:
                        if isinstance(production_crop["Production"].iloc[0], str):
                            bio_energy = bio_energy + float(production_crop["Production"].iloc[0].replace(',',''))* param["Biomass"]["agro_rpr"][crop][residue]*param["Biomass"]["agro_af"][crop][residue]*param["Biomass"]["agro_lhv"][crop][residue]
                            bio_co2 = bio_co2 + float(production_crop["Production"].iloc[0].replace(',',''))* param["Biomass"]["agro_rpr"][crop][residue]*param["Biomass"]["agro_af"][crop][residue]*param["Biomass"]["agro_emission factor"]
                        else:
                            bio_energy = bio_energy + float(production_crop["Production"].iloc[0])* param["Biomass"]["agro_rpr"][crop][residue]*param["Biomass"]["agro_af"][crop][residue]*param["Biomass"]["agro_lhv"][crop][residue]
                            bio_co2 = bio_co2 + float(production_crop["Production"].iloc[0])* param["Biomass"]["agro_rpr"][crop][residue]*param["Biomass"]["agro_af"][crop][residue]*param["Biomass"]["agro_emission factor"]

            energy_lu_crop = A_lu_crop * 2 * bio_energy / (2*n_lu_crop+n_lu_veg)
            energy_lu_veg = A_lu_veg * bio_energy / (2*n_lu_crop+n_lu_veg)
            
            co2_lu_crop = A_lu_crop * 2 * bio_co2 / (2*n_lu_crop+n_lu_veg)
            co2_lu_veg = A_lu_veg * bio_co2 / (2*n_lu_crop+n_lu_veg)
        
            energy_reg = energy_lu_crop + energy_lu_veg
            A_Bioenergy = A_Bioenergy + energy_reg

            co2_reg = co2_lu_crop + co2_lu_veg
            A_Bioco2 = A_Bioco2 + co2_reg
     
    A_lu_forest = np.multiply(A_lu_forest, A_country_area)
    A_lu_forest = np.multiply(A_lu_forest, A_include)
    n_lu_forest = np.sum(A_lu_forest)
    forest = param["Biomass"]["forest"]
    
    if n_lu_forest:
        energy_lu_forest = A_lu_forest * (forest["Wood fuel, coniferous"]*forest["density, coniferous"] + forest["Wood fuel, non-coniferous"]*forest["density, non-coniferous"])*forest["rpr"]*forest["af"]*forest["lhv"] / n_lu_forest
        co2_lu_forest = A_lu_forest * (forest["Wood fuel, coniferous"]*forest["density, coniferous"] + forest["Wood fuel, non-coniferous"]*forest["density, non-coniferous"])*forest["rpr"]*forest["af"]*forest["emission factor"] / n_lu_forest
        
        A_Bioenergy = A_Bioenergy + energy_lu_forest
        A_Bioco2 = A_Bioco2 + co2_lu_forest
    
    n_animal = 0
    for animal in param["Biomass"]["livestock"]["animal"]:
        with rasterio.open(paths["LS"]+animal+".tif") as src:
            A_LS_animal = src.read(1)
        A_LS_animal = np.flipud(A_LS_animal)
        A_LS_animal = np.multiply(A_LS_animal,A_protect)
        A_LS_animal = np.multiply(A_LS_animal,A_country_area)
        energy_ls_animal = A_LS_animal * param["Biomass"]["livestock"]["rpr"][n_animal] * param["Biomass"]["livestock"]["af"][n_animal] * param["Biomass"]["livestock"]["lhv"][n_animal]
        co2_ls_animal = A_LS_animal * param["Biomass"]["livestock"]["rpr"][n_animal] * param["Biomass"]["livestock"]["af"][n_animal] * param["Biomass"]["livestock"]["emission factor"]
        
        A_Bioenergy = A_Bioenergy + energy_ls_animal
        A_Bioco2 = A_Bioco2 + co2_ls_animal
        
        n_animal = n_animal+1
        
    array2raster(paths["BIOMASS_ENERGY"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_Bioenergy)
    print("files saved: " + paths["BIOMASS_ENERGY"])
    create_json(paths["BIOMASS_ENERGY"], param, ["region_name", "landuse", "Biomass", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "BIOMASS_ENERGY"])
    
    array2raster(paths["BIOMASS_CO2"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_Bioco2)
    print("files saved: " + paths["BIOMASS_CO2"])
    create_json(paths["BIOMASS_CO2"], param, ["region_name", "landuse", "Biomass", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "BIOMASS_CO2"])
    
    timecheck("End") 