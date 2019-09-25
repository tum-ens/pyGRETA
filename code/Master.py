from model_functions import *


def initialization():
    """
    This function reads the user-defined parameters and paths from :mod:`config.py`, then adds additional parameters related
    to the shapefiles. 
    First, it saves the spatial scope of the problem.
    Then, it distinguishes between countries, exclusive economic zones and subregions. For each one of them, 
    it saves the geodataframes, the number of features, and the coordinates of the bounding boxes of each feature.
    Finally, it saves the number of rows and columns in the low and righ resolution, and a georeference dictionary
    used for saving tif files.

    :return: The updated dictionaries param and paths.
    :rtype: tuple (dict, dict)
    """
    timecheck('Start')

    # import param and paths
    from config import configuration
    paths, param = configuration()

    # Read shapefile of scope
    scope_shp = gpd.read_file(paths["spatial_scope"])
    param["spatial_scope"] = define_spatial_scope(scope_shp)

    res_weather = param["res_weather"]
    res_desired = param["res_desired"]
    Crd_all = crd_merra(param["spatial_scope"], res_weather)[0]
    param["Crd_all"] = Crd_all
    ymax, xmax, ymin, xmin = Crd_all
    bounds_box = Polygon([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)])

    timecheck('Read shapefile of countries')
    # Extract land areas
    countries_shp = gpd.read_file(paths["Countries"], bbox=scope_shp)
    countries_shp = countries_shp.to_crs({'init': 'epsg:4326'})

    # Crop all polygons and take the part inside the bounding box
    countries_shp['geometry'] = countries_shp['geometry'].intersection(bounds_box)
    countries_shp = countries_shp[countries_shp.geometry.area > 0]
    param["regions_land"] = countries_shp
    param["nRegions_land"] = len(param["regions_land"])
    Crd_regions_land = np.zeros((param["nRegions_land"], 4))

    for reg in range(0, param["nRegions_land"]):
        # Box coordinates for MERRA2 data
        r = countries_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd_regions_land[reg, :] = crd_merra(box, res_weather)

    timecheck('Read shapefile of EEZ')
    # Extract sea areas
    eez_shp = gpd.read_file(paths["EEZ_global"], bbox=scope_shp)
    eez_shp = eez_shp.to_crs({'init': 'epsg:4326'})

    # Crop all polygons and take the part inside the bounding box
    eez_shp['geometry'] = eez_shp['geometry'].intersection(bounds_box)
    eez_shp = eez_shp[eez_shp.geometry.area > 0]
    param["regions_sea"] = eez_shp
    param["nRegions_sea"] = len(param["regions_sea"])
    Crd_regions_sea = np.zeros((param["nRegions_sea"], 4))

    for reg in range(0, param["nRegions_sea"]):
        # Box coordinates for MERRA2 data
        r = eez_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd_regions_sea[reg, :] = crd_merra(box, res_weather)

    timecheck('Read shapefile of subregions')
    # Read shapefile of regions
    regions_shp = gpd.read_file(paths["subregions"], bbox=scope_shp)
    regions_shp = regions_shp.to_crs({'init': 'epsg:4326'})

    # Crop all polygons and take the part inside the bounding box
    regions_shp['geometry'] = regions_shp['geometry'].intersection(bounds_box)
    regions_shp = regions_shp[regions_shp.geometry.area > 0]
    regions_shp.sort_values(by=['NAME_SHORT'], inplace=True)
    regions_shp.reset_index(inplace=True)
    param["regions_sub"] = regions_shp
    param["nRegions_sub"] = len(param["regions_sub"])
    Crd_regions_sub = np.zeros((param["nRegions_sub"], 4))

    for reg in range(0, param["nRegions_sub"]):
        # Box coordinates for MERRA2 data
        r = regions_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd_regions_sub[reg, :] = crd_merra(box, res_weather)

    # Saving parameters
    param["Crd_subregions"] = Crd_regions_sub
    param["Crd_regions"] = np.concatenate((Crd_regions_land, Crd_regions_sea), axis=0)

    # Indices and matrix dimensions
    Ind_all_low = ind_merra(Crd_all, Crd_all, res_weather)
    Ind_all_high = ind_merra(Crd_all, Crd_all, res_desired)

    param["m_high"] = int((Ind_all_high[:, 0] - Ind_all_high[:, 2] + 1)[0])  # number of rows
    param["n_high"] = int((Ind_all_high[:, 1] - Ind_all_high[:, 3] + 1)[0])  # number of columns
    param["m_low"] = int((Ind_all_low[:, 0] - Ind_all_low[:, 2] + 1)[0])  # number of rows
    param["n_low"] = int((Ind_all_low[:, 1] - Ind_all_low[:, 3] + 1)[0])  # number of columns
    param["GeoRef"] = calc_geotiff(Crd_all, res_desired)
    timecheck('End')

    # Display initial information
    print('\nRegion: ' + param["region_name"] + ' - Year: ' + str(param["year"]))
    print('Folder Path: ' + paths["region"] + '\n')

    return paths, param


def generate_weather_files(paths, param):
    """
    This function reads the daily NetCDF data (from MERRA-2) for SWGDN, SWTDN, T2M, U50m, and V50m,
    and saves them in matrices with yearly time series with low spatial resolution.
    This function has to be run only once.

    :param paths: Dictionary including the paths to the MERRA-2 input files *MERRA_IN*, and to the desired output locations for *T2M*, *W50M* and *CLEARNESS*.
    :type paths: dict
    :param param: Dictionary including the year and the spatial scope.
    :type param: dict
    :return: The files T2M.mat, W50M.mat, and CLEARNESS.mat are saved directly in the defined paths, along with their metadata in JSON files.
    :rtype: None
    """

    timecheck('Start')
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
        sys.stdout.write('\r')
        sys.stdout.write('Reading NetCDF files ' + '[%-50s] %d%%' % (
            '=' * ((status * 50) // delta), (status * 100) // delta))
        sys.stdout.flush()

        tomorrow = date + pd.Timedelta('1 day')
        if date.day == 29 and date.month == 2:
            continue

        # Name and path of the NetCDF file to be read
        name = paths["MERRA_IN"] + 'MERRA2_400.tavg1_2d_rad_Nx.' + date.strftime('%Y%m%d') + '.SUB.nc'
        name2 = paths["MERRA_IN"] + 'MERRA2_400.tavg1_2d_slv_Nx.' + date.strftime('%Y%m%d') + '.SUB.nc'

        # Read NetCDF file, extract hourly tables
        with h5netcdf.File(name, 'r') as f:
            # [time, lat 361, lon 576]
            swgdn = np.transpose(subset(f['SWGDN'], param), [1, 2, 0])
            if SWGDN.size == 0:
                SWGDN = swgdn
            else:
                SWGDN = np.concatenate((SWGDN, swgdn), axis=2)

            swtdn = np.transpose(subset(f['SWTDN'], param), [1, 2, 0])
            if SWTDN.size == 0:
                SWTDN = swtdn
            else:
                SWTDN = np.concatenate((SWTDN, swtdn), axis=2)

        with h5netcdf.File(name2, 'r') as f:
            t2m = np.transpose(subset(f['T2M'], param), [1, 2, 0])
            if T2M.size == 0:
                T2M = t2m
            else:
                T2M = np.concatenate((T2M, t2m), axis=2)

            u50m = np.transpose(subset(f['U50M'], param), [1, 2, 0])
            if U50M.size == 0:
                U50M = u50m
            else:
                U50M = np.concatenate((U50M, u50m), axis=2)

            v50m = np.transpose(subset(f['V50M'], param), [1, 2, 0])
            if V50M.size == 0:
                V50M = v50m
            else:
                V50M = np.concatenate((V50M, v50m), axis=2)
        if date.year != tomorrow.year:
            # Create the overall wind speed
            W50M = abs(U50M + (1j * V50M))
            # Calculate the clearness index
            CLEARNESS = np.divide(SWGDN, SWTDN, where=SWTDN != 0)

            sys.stdout.write('\n')
            timecheck('Writing Files: T2M, W50M, CLEARNESS')
            hdf5storage.writes({'T2M': T2M}, paths["T2M"], store_python_metadata=True, matlab_compatible=True)
            create_json(paths["T2M"], param, ["MERRA_coverage", "region_name", "Crd_all", "res_weather"], paths,
                        ["MERRA_IN", "T2M"])
            hdf5storage.writes({'W50M': W50M}, paths["W50M"], store_python_metadata=True, matlab_compatible=True)
            create_json(paths["W50M"], param, ["MERRA_coverage", "region_name", "Crd_all", "res_weather"], paths,
                        ["MERRA_IN", "W50M"])
            hdf5storage.writes({'CLEARNESS': CLEARNESS}, paths["CLEARNESS"], store_python_metadata=True,
                               matlab_compatible=True)
            create_json(paths["CLEARNESS"], param, ["MERRA_coverage", "region_name", "Crd_all", "res_weather"], paths,
                        ["MERRA_IN", "CLEARNESS"])
    timecheck('End')


def clean_weather_data(paths, param):
    """
    This function detects data outliers in the wind input file W50M.mat. An outlier is a data point, for which
    the absolute value of the difference between the yearly average value and the mean of the direct neighbors
    (Moore neighborhood) is higher than a user-defined threshold *MERRA_correction*. It replaces the hourly values
    with the hourly values of the mean of the neighbors, and overwrites the original W50M.mat file.

    :param paths: Dictionary including the path to the file W50M.mat.
    :type paths: dict
    :param param: Dictionary including the threshold value *MERRA_correction*.
    :type param: dict
    :return: The file W50M.mat is overwritten after the correction, along with its metadata in a JSON file.
    :rtype: None
    """

    timecheck('Start')

    # Read Wind Data
    W50M = hdf5storage.read('W50M', paths["W50M"])
    wind = np.mean(W50M, 2)

    # Set convolution mask
    kernel = np.ones((3, 3))
    kernel[1, 1] = 0

    # Compute average Convolution
    averagewind = generic_filter(wind, np.nanmean, footprint=kernel, mode='constant', cval=np.NaN)
    ratio = wind / averagewind

    # Extract over threshold Points
    points = np.where(abs(ratio - np.mean(ratio)) > param["MERRA_correction"])

    # Correct points hourly
    for t in range(W50M.shape[2]):
        W50M[points[0], points[1], t] = W50M[points[0], points[1], t] / ratio[points[0], points[1]]

    # Save corrected Wind
    hdf5storage.writes({'W50M': W50M}, paths["W50M"], store_python_metadata=True, matlab_compatible=True)
    create_json(paths["W50M"], param, ["MERRA_coverage", "region_name", "Crd_all", "res_weather", "MERRA_correction"],
                paths, ["MERRA_IN", "W50M"])
    timecheck('End')


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

    timecheck('Start')
    timecheck('Start Land')
    # Extract land areas
    countries_shp = param["regions_land"]
    Crd_regions_land = param["Crd_regions"][:nRegions_land]
    Ind = ind_merra(Crd_regions_land, Crd_all, res_desired)
    A_land = np.zeros((m_high, n_high))
    status = 0
    for reg in range(0, param["nRegions_land"]):
        # Show status bar
        status = status + 1
        sys.stdout.write('\r')
        sys.stdout.write('Creating A_land ' + '[%-50s] %d%%' % (
            '=' * ((status * 50) // nRegions_land), (status * 100) // nRegions_land))
        sys.stdout.flush()

        # Calculate A_region
        A_region = calc_region(countries_shp.iloc[reg], Crd_regions_land[reg, :], res_desired, GeoRef)

        # Include A_region in A_land
        A_land[(Ind[reg, 2] - 1):Ind[reg, 0], (Ind[reg, 3] - 1):Ind[reg, 1]] = \
            A_land[(Ind[reg, 2] - 1):Ind[reg, 0], (Ind[reg, 3] - 1):Ind[reg, 1]] + A_region
    # Saving file
    array2raster(paths["LAND"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_land)
    create_json(paths["LAND"], param,
                ["region_name", "m_high", "n_high", "Crd_all", "res_desired", "GeoRef", "nRegions_land"], paths,
                ["Countries", "LAND"])
    print('\nfiles saved: ' + paths["LAND"])
    timecheck('Finish Land')

    timecheck('Start Sea')
    # Extract sea areas
    eez_shp = param["regions_sea"]
    Crd_regions_sea = param["Crd_regions"][-nRegions_sea:]
    Ind = ind_merra(Crd_regions_sea, Crd_all, res_desired)
    A_sea = np.zeros((m_high, n_high))
    status = 0
    for reg in range(0, param["nRegions_sea"]):
        # Show status bar
        status = status + 1
        sys.stdout.write('\r')
        sys.stdout.write('Creating A_sea ' + '[%-50s] %d%%' % (
            '=' * ((status * 50) // param["nRegions_sea"]), (status * 100) // param["nRegions_sea"]))
        sys.stdout.flush()

        # Calculate A_region
        A_region = calc_region(eez_shp.iloc[reg], Crd_regions_sea[reg, :], res_desired, GeoRef)

        # Include A_region in A_sea
        A_sea[(Ind[reg, 2] - 1):Ind[reg, 0], (Ind[reg, 3] - 1):Ind[reg, 1]] = \
            A_sea[(Ind[reg, 2] - 1):Ind[reg, 0], (Ind[reg, 3] - 1):Ind[reg, 1]] + A_region

    # Fixing pixels on the borders to avoid duplicates
    A_sea[A_sea > 0] = 1
    A_sea[A_land > 0] = 0
    # Saving file
    array2raster(paths["EEZ"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_sea)
    create_json(paths["EEZ"], param,
                ["region_name", "m_high", "n_high", "Crd_all", "res_desired", "GeoRef", "nRegions_sea"], paths,
                ["EEZ_global", "EEZ"])
    print('\nfiles saved: ' + paths["EEZ"])
    timecheck('Finish Sea')


def generate_subregions(paths, param):
    """
    This function reads the shapefiles of the subregions within the scope, and creates a raster out of it.

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

    timecheck('Start Subregions')
    # Read shapefile of regions
    regions_shp = param["regions_sub"]
    Crd_regions_sub = param["Crd_subregions"]
    Ind = ind_merra(Crd_regions_sub, Crd_all, res_desired)
    A_sub = np.zeros((m_high, n_high))
    status = 0
    for reg in range(0, nRegions_sub):
        # Show status bar
        status = status + 1
        sys.stdout.write('\r')
        sys.stdout.write('Creating A_subregions ' + '[%-50s] %d%%' % (
            '=' * ((status * 50) // nRegions_sub), (status * 100) // nRegions_sub))
        sys.stdout.flush()

        # Calculate A_region
        A_region = calc_region(regions_shp.iloc[reg], Crd_regions_sub[reg, :], res_desired, GeoRef)

        # Include A_region in A_sub
        A_sub[(Ind[reg, 2] - 1):Ind[reg, 0], (Ind[reg, 3] - 1):Ind[reg, 1]] = \
            A_sub[(Ind[reg, 2] - 1):Ind[reg, 0], (Ind[reg, 3] - 1):Ind[reg, 1]] + A_region

    # Fixing pixels on the borders
    with rasterio.open(paths["EEZ"]) as src:
        A_sea = np.flipud(src.read(1)).astype(int)
    with rasterio.open(paths["LAND"]) as src:
        A_land = np.flipud(src.read(1)).astype(int)
    A_sub = A_sub * (A_land + A_sea)

    # Saving file
    array2raster(paths["SUB"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_sub)
    create_json(paths["SUB"], param,
                ["subregions_name", "m_high", "n_high", "Crd_all", "res_desired", "GeoRef", "nRegions_sea"], paths,
                ["subregions", "SUB"])
    print('\nfiles saved: ' + paths["SUB"])
    timecheck('Finish Subregions')

    timecheck('End')


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

    timecheck('Start')
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, res_desired)[0]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["LU_global"]) as src:
        w = src.read(1, window=windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]),
                                                          slice(Ind[3] - 1, Ind[1])))
    w = np.flipud(w)
    array2raster(paths["LU"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], w)
    create_json(paths["LU"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["LU_global", "LU"])
    print("files saved: " + paths["LU"])
    timecheck('End')


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

    timecheck('Start')
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, res_desired)[0]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["Bathym_global"]) as src:
        A_BATH = src.read(1)
    A_BATH = resizem(A_BATH, 180 * 240, 360 * 240)
    A_BATH = np.flipud(A_BATH[Ind[0] - 1: Ind[2], Ind[3] - 1: Ind[1]])
    array2raster(paths['BATH'], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_BATH)
    create_json(paths["BATH"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                ["Bathym_global", "BATH"])
    print("files saved: " + paths["BATH"])
    timecheck('End')


def generate_topography(paths, param):
    """
    This function reads the tiles that make the global map of topography, picks those that lie completely or partially in the scope,
    and creates a raster out of them for the desired scope. The values are in meter.

    :param paths: Dictionary including the paths to the tiles of the global topography raster Topo_tiles and to the output path TOPO.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict
    :return: The tif file for TOPO is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """

    timecheck('Start')
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, res_desired)[0]
    GeoRef = param["GeoRef"]
    Topo = np.zeros((180 * 240, 360 * 240))
    tile_extents = np.zeros((24, 4), dtype=int)
    i = 1
    j = 1
    for letter in char_range('A', 'X'):
        north = (i - 1) * 45 * 240 + 1
        east = j * 60 * 240
        south = i * 45 * 240
        west = (j - 1) * 60 * 240 + 1
        tile_extents[ord(letter) - ord('A'), :] = [north, east, south, west]
        j = j + 1
        if j == 7:
            i = i + 1
            j = 1
    n_min = (Ind[0] // (45 * 240)) * 45 * 240 + 1
    e_max = (Ind[1] // (60 * 240) + 1) * 60 * 240
    s_max = (Ind[2] // (45 * 240) + 1) * 45 * 240
    w_min = (Ind[3] // (60 * 240)) * 60 * 240 + 1

    need = np.logical_and((np.logical_and((tile_extents[:, 0] >= n_min), (tile_extents[:, 1] <= e_max))),
                          np.logical_and((tile_extents[:, 2] <= s_max), (tile_extents[:, 3] >= w_min)))

    status = 0
    for letter in char_range('A', 'X'):
        index = ord(letter) - ord('A')
        if need[index]:
            # Show status bar
            status = status + 1
            sys.stdout.write('\r')
            sys.stdout.write('Generating topography map from tiles ' + '[%-50s] %d%%' % (
                '=' * ((status * 50) // sum(need)), (status * 100) // sum(need)))
            sys.stdout.flush()

            with rasterio.open(paths["Topo_tiles"] + '15-' + letter + '.tif') as src:
                tile = src.read()
            Topo[tile_extents[index, 0] - 1: tile_extents[index, 2],
            tile_extents[index, 3] - 1: tile_extents[index, 1]] = \
                tile[0, 0:-1, 0:-1]

    A_TOPO = np.flipud(Topo[Ind[0] - 1:Ind[2], Ind[3] - 1:Ind[1]])
    array2raster(paths["TOPO"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_TOPO)
    create_json(paths["TOPO"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                ["Topo_tiles", "TOPO"])
    print("\nfiles saved: " + paths["TOPO"])
    timecheck('End')


def generate_slope(paths, param):
    """
    This function reads the topography raster for the scope, and creates a raster of slope out of it. The slope is calculated in
    percentage, although this can be changed easily at the end of the code.

    :param paths: Dictionary including the paths to the topography map of the scope TOPO and to the output path SLOPE.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict
    :return: The tif file for SLOPE is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """

    timecheck('Start')
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, res_desired)[0]
    GeoRef = param["GeoRef"]
    Lat1 = np.arange(-90, 90, 1 / 240)
    Lat2 = np.arange(-90 + 1 / 240, 90 + 1 / 240, 1 / 240)
    latMid = (Lat1 + Lat2) / 2
    deltaLat = abs(Lat1 - Lat2)

    Lat1 = np.arange(-90, 90, 1 / 240)
    Lat2 = np.arange(-90 + 1 / 240, 90 + 1 / 240, 1 / 240)
    latMid_2 = (Lat1 + Lat2) / 2

    Lon1 = np.arange(-180, 180, 1 / 240)
    Lon2 = np.arange(-180 + 1 / 240, 180 + 1 / 240, 1 / 240)
    deltaLon = abs(Lon1 - Lon2)

    m_per_deg_lat = 111132.954 - 559.822 * cos(np.deg2rad(2 * latMid)) + 1.175 * cos(np.deg2rad(4 * latMid))
    m_per_deg_lon = (np.pi / 180) * 6367449 * cos(np.deg2rad(latMid_2))

    x_cell = repmat(deltaLon, 180 * 240, 1) * repmat(m_per_deg_lon, 360 * 240, 1).T
    x_cell = x_cell[Ind[0] - 1:Ind[2], Ind[3] - 1: Ind[1]]
    x_cell = np.flipud(x_cell)

    y_cell = repmat((deltaLat * m_per_deg_lat), 360 * 240, 1).T
    y_cell = y_cell[Ind[0] - 1:Ind[2], Ind[3] - 1: Ind[1]]
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
    create_json(paths["SLOPE"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["TOPO", "SLOPE"])
    print("files saved: " + paths["SLOPE"])
    timecheck('End')


def generate_population(paths, param):
    """
    This function reads the tiles that make the global map of population, picks those that lie completely or partially in the scope,
    and creates a raster out of them for the desired scope. The values are in population by pixel.

    :param paths: Dictionary including the paths to the tiles of the global population raster Pop_tiles and to the output path POP.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict
    :return: The tif file for POP is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck('Start')
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, res_desired)[0]
    GeoRef = param["GeoRef"]
    Pop = np.zeros((180 * 240, 360 * 240))
    tile_extents = np.zeros((24, 4), dtype=int)
    i = 1
    j = 1

    for letter in char_range('A', 'X'):
        north = (i - 1) * 45 * 240 + 1
        east = j * 60 * 240
        south = i * 45 * 240
        west = (j - 1) * 60 * 240 + 1
        tile_extents[ord(letter) - ord('A'), :] = [north, east, south, west]
        j = j + 1
        if j == 7:
            i = i + 1
            j = 1
    n_min = (Ind[0] // (45 * 240)) * 45 * 240 + 1
    e_max = (Ind[1] // (60 * 240) + 1) * 60 * 240
    s_max = (Ind[2] // (45 * 240) + 1) * 45 * 240
    w_min = (Ind[3] // (60 * 240)) * 60 * 240 + 1

    need = np.logical_and((np.logical_and((tile_extents[:, 0] >= n_min), (tile_extents[:, 1] <= e_max))),
                          np.logical_and((tile_extents[:, 2] <= s_max), (tile_extents[:, 3] >= w_min)))

    status = 0
    for letter in char_range('A', 'X'):
        index = ord(letter) - ord('A')
        if need[index]:
            # Show status bar
            status = status + 1
            sys.stdout.write('\r')
            sys.stdout.write('Generating population map from tiles ' + '[%-50s] %d%%' % (
                '=' * ((status * 50) // sum(need)), (status * 100) // sum(need)))
            sys.stdout.flush()

            with rasterio.open(paths["Pop_tiles"] + letter + '.tif') as src:
                tile = src.read()
            Pop[tile_extents[index, 0] - 1: tile_extents[index, 2],
            tile_extents[index, 3] - 1: tile_extents[index, 1]] = \
                tile[0]

    A_POP = np.flipud(Pop[Ind[0] - 1:Ind[2], Ind[3] - 1:Ind[1]])
    array2raster(paths["POP"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_POP)
    create_json(paths["POP"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["Pop_tiles", "POP"])
    print("\nfiles saved: " + paths["POP"])
    timecheck('End')


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
    timecheck('Start')
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
    if not field_exists('Raster', shp_path):
        new_field = ogr.FieldDefn('Raster', ogr.OFTInteger)
        layer.CreateField(new_field)

        for feat in layer:
            pt = feat.GetField('IUCN_CAT')
            feat.SetField('Raster', protection_type[pt])
            layer.SetFeature(feat)
            feat = None

    # Create a second (modified) layer
    outdriver = ogr.GetDriverByName('MEMORY')
    source = outdriver.CreateDataSource('memData')

    # Create the raster dataset
    memory_driver = gdal.GetDriverByName('GTiff')
    out_raster_ds = memory_driver.Create(paths["PA"], ncol, nrow, 1,
                                         gdal.GDT_Byte)

    # Set the ROI image's projection and extent to our input raster's projection and extent
    out_raster_ds.SetProjection(proj)
    out_raster_ds.SetGeoTransform(ext)

    # Fill our output band with the 0 blank, no class label, value
    b = out_raster_ds.GetRasterBand(1)
    b.Fill(0)

    # Rasterize the shapefile layer to our new dataset
    gdal.RasterizeLayer(out_raster_ds,  # output to our new dataset
                        [1],  # output to our new dataset's first band
                        layer,  # rasterize this layer
                        None, None,  # don't worry about transformations since we're in same projection
                        [0],  # burn value 0
                        ['ALL_TOUCHED=FALSE',  # rasterize all pixels touched by polygons
                         'ATTRIBUTE=Raster']  # put raster values according to the 'Raster' field values
                        )
    create_json(paths["PA"], param, ["region_name", "protected_areas", "Crd_all", "res_desired", "GeoRef"], paths,
                ["Protected", "PA"])

    # Close dataset
    out_raster_ds = None
    print("files saved: " + paths["PA"])
    timecheck('End')


def generate_buffered_population(paths, param):
    """
    This function reads the land use raster, identifies urban areas, and exclude pixels around them based on a
    user-defined buffer (buffer_pixel_amount). It creates a masking raster of boolean values (0 or 1) for the scope.
    Zero means the pixel is excluded, one means it is suitable.
    The function is useful in case there is a policy to exclude renewable energy projects next to urban settlements.

    :param paths: Dictionary including the path to the land use raster for the scope, and to the output path BUFFER.
    :type paths: dict
    :param param: Dictionary including the user-defined buffer (buffer_pixel_amount), the urban type within the land use map (type_urban), and the georeference dictionary.
    :type param: dict
    :return: The tif file for BUFFER is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck('Start')
    buffer_pixel_amount = param["WindOn"]["mask"]["buffer_pixel_amount"]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
    A_lu = np.flipud(A_lu).astype(int)
    A_lu = A_lu == param["landuse"]["type_urban"]  # Land use type for Urban and built-up
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_lu_buffered = generic_filter(A_lu, np.nanmean, footprint=kernel, mode='constant', cval=np.NaN)
    A_notPopulated = (~A_lu_buffered).astype(int)

    array2raster(paths["BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"],
                 A_notPopulated)
    create_json(paths["BUFFER"], param, ["region_name", "landuse", "WindOn", "Crd_all", "res_desired", "GeoRef"], paths,
                ["LU", "BUFFER"])
    print("files saved: " + paths["BUFFER"])
    timecheck('End')


def generate_area(paths, param):
    timecheck('Start')
    Crd_all = param["Crd_all"]
    n_high = param["n_high"]
    res_desired = param["res_desired"]

    # Calculate available Area
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

    lowerSliceAreas = np.pi * b ** 2 * ((2 * np.arctanh(e * sin(f_lower))) / (2 * e) +
                                        (sin(f_lower) / (zp_lower * zm_lower)))

    # Upper slice areas
    # Areas between the equator and the upper pixel latitudes circling the globe
    f_upper = np.deg2rad(lat_vec + res_desired[0])

    zm_upper = 1 - (e * sin(f_upper))
    zp_upper = 1 + (e * sin(f_upper))

    upperSliceAreas = np.pi * b ** 2 * ((2 * np.arctanh((e * sin(f_upper)))) / (2 * e) +
                                        (sin(f_upper) / (zp_upper * zm_upper)))

    # Pixel areas
    # Finding the latitudinal pixel-sized globe slice areas then dividing them by the longitudinal pixel size
    area_vec = ((upperSliceAreas - lowerSliceAreas) * res_desired[1] / 360).T
    A_area = np.tile(area_vec, (1, n_high))

    # Save to HDF File
    hdf5storage.writes({'A_area': A_area}, paths["AREA"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths["AREA"])

    timecheck('End')


def generate_wind_correction(paths, param):
    timecheck('Start')
    res_correction_on = param["WindOn"]["resource"]["res_correction"]
    res_correction_off = param["WindOff"]["resource"]["res_correction"]
    topo_correction = param["WindOn"]["resource"]["topo_correction"]
    GeoRef = param["GeoRef"]
    landuse = param["landuse"]
    with rasterio.open(paths["LU"]) as src:
        A_lu = np.flipud(src.read(1)).astype(int)
    A_hellmann = changem(A_lu, landuse["hellmann"], landuse["type"]).astype(float)

    # Onshore resolution correction
    if 'WindOn' in param["technology"]:
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
            A_cf_on = ((turbine_height_on / 50) * turbine_height_on /
                       A_gradient_height) ** A_hellmann / resizem(Sigma, m_high, n_high)
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
                calc_gwa_correction(param, paths)
            gwa_correction = hdf5storage.read('correction_' + param["WindOn"]["resource"]["topo_weight"],
                                              paths["CORR_GWA"])
            A_cf_on = A_cf_on * gwa_correction
        array2raster(paths["CORR_ON"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_cf_on)
        print("\nfiles saved: " + paths["CORR_ON"])

    # Offshore resolution correction
    if 'WindOff' in param["technology"]:
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
            A_cf_off = ((turbine_height_off / 50) * turbine_height_off /
                        A_gradient_height) ** A_hellmann / resizem(Sigma, m_high, n_high)
            del A_gradient_height, Sigma
        else:
            A_cf_off = (turbine_height_off / 50) ** A_hellmann
        del A_hellmann
        with rasterio.open(paths["EEZ"]) as src:
            A_eez = np.flipud(src.read(1)).astype(int)
        A_cf_off = A_cf_off * A_eez

        array2raster(paths["CORR_OFF"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_cf_off)
        print("\nfiles saved: " + paths["CORR_OFF"])
    timecheck('End')


def calculate_FLH(paths, param, tech):
    timecheck('Start')
    print('Region: ' + param["region_name"])

    if tech in ["WindOn", "WindOff"]:
        print("\n" + tech + " - HUB_HEIGHTS: " + str(param[tech]["technical"]["hub_height"]))
    elif tech in ["PV"] and 'orientation' in param["PV"]["technical"].keys():
        print("\n" + tech + " - Orientation: " + str(param[tech]["technical"]["orientation"]))

    nproc = param["nproc"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    CPU_limit = np.full((1, nproc), param["CPU_limit"])

    if tech == "WindOff":
        with rasterio.open(paths["EEZ"]) as src:
            w = src.read(1)
    else:
        res_weather = param["res_weather"]
        Crd_all = param["Crd_all"]
        Ind = ind_merra(Crd_all, Crd_all, res_weather)[0]
        with rasterio.open(paths["LAND"]) as src:
            w = src.read(1)
    param["Ind_nz"] = np.nonzero(np.flipud(w))
    del w

    # Obtain weather and correction matrices
    merraData, rasterData = get_merra_raster_Data(paths, param, tech)

    if tech in ['PV', 'CSP']:

        day_filter = np.nonzero(merraData["CLEARNESS"][Ind[2] - 1:Ind[0], Ind[3] - 1:Ind[1], :].sum(axis=(0, 1)))
        list_hours = np.arange(0, 8760)
        if nproc == 1:
            param["status_bar_limit"] = list_hours[-1]
            results = calc_FLH_solar(list_hours[day_filter], [param, tech, rasterData, merraData])
        else:
            list_hours = np.array_split(list_hours[day_filter], nproc)
            param["status_bar_limit"] = list_hours[0][-1]
            results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(calc_FLH_solar,
                                                                                               product(list_hours,
                                                                                                       [[param,
                                                                                                         tech,
                                                                                                         rasterData,
                                                                                                         merraData]]))
    elif tech in ['WindOn', 'WindOff']:

        list_hours = np.array_split(np.arange(0, 8760), nproc)
        param["status_bar_limit"] = list_hours[0][-1]
        results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(calc_FLH_wind,
                                                                                           product(list_hours,
                                                                                                   [[param,
                                                                                                     tech, rasterData,
                                                                                                     merraData]]))
    # Collecting results
    FLH = np.zeros((m_high, n_high))
    if nproc > 1:
        for p in range(len(results)):
            FLH[param["Ind_nz"]] = FLH[param["Ind_nz"]] + results[p]
    else:
        FLH[param["Ind_nz"]] = results
    FLH[FLH == 0] = np.nan

    hdf5storage.writes({'FLH': FLH}, paths[tech]["FLH"], store_python_metadata=True, matlab_compatible=True)
    print("\nfiles saved: " + paths[tech]["FLH"])

    # Save GEOTIFF files
    if param["savetiff"]:
        GeoRef = param["GeoRef"]
        array2raster(changeExt2tif(paths[tech]["mask"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelHeight"],
                     FLH)
        print("files saved:" + changeExt2tif(paths[tech]["mask"]))

    timecheck('End')


def masking(paths, param, tech):
    timecheck('Start')
    mask = param[tech]["mask"]


    if tech in ['PV', 'CSP']:
        with rasterio.open(paths["PA"]) as src:
            A_protect = src.read(1)
            A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified
        # Exclude protection categories that are not suitable
        A_suitability_pa = changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)
        A_suitability_pa = (A_suitability_pa > 0).astype(int)
        with rasterio.open(paths["LU"]) as src:
            A_lu = src.read(1)
            A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified
        # Exclude landuse types types that are not suitable
        A_suitability_lu = changem(A_lu, mask["lu_suitability"], param["landuse"]["type"]).astype(float)
        A_suitability_lu = (A_suitability_lu > 0).astype(int)
        with rasterio.open(paths["SLOPE"]) as src:
            A_slope = src.read(1)
            A_slope = np.flipud(A_slope)  # Slope in percentage
            A_slope = (A_slope <= mask["slope"]).astype(int)
        # Irrelevant parameters
        A_notPopulated = 1
        A_bathymetry = 1

    if tech == 'WindOn':
        with rasterio.open(paths["PA"]) as src:
            A_protect = src.read(1)
            A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified
        # Exclude protection categories that are not suitable
        A_suitability_pa = changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)
        A_suitability_pa = (A_suitability_pa > 0).astype(int)
        with rasterio.open(paths["LU"]) as src:
            A_lu = src.read(1)
            A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified
        # Exclude landuse types types that are not suitable
        A_suitability_lu = changem(A_lu, mask["lu_suitability"], param["landuse"]["type"]).astype(float)
        A_suitability_lu = (A_suitability_lu > 0).astype(int)
        with rasterio.open(paths["SLOPE"]) as src:
            A_slope = src.read(1)
            A_slope = np.flipud(A_slope)  # Slope in percentage
            A_slope = (A_slope <= mask["slope"]).astype(int)
        with rasterio.open(paths["BUFFER"]) as src:
            A_notPopulated = src.read(1)
            A_notPopulated = (np.flipud(A_notPopulated)).astype(int)  # Is 1 for not populated areas
        # Irrelevant parameters
        A_bathymetry = 1

    if tech == 'WindOff':
        with rasterio.open(paths["EEZ"]) as src:
            A_suitability_lu = src.read(1)
            A_suitability_lu = np.flipud(A_suitability_lu).astype(int)
        with rasterio.open(paths["PA"]) as src:
            A_protect = src.read(1)
            A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified
        # Exclude protection categories that are not suitable
        A_suitability_pa = changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)
        A_suitability_pa = (A_suitability_pa > 0).astype(int)
        with rasterio.open(paths["BATH"]) as src:
            A_bathymetry = src.read(1)
            A_bathymetry = np.flipud(A_bathymetry)  # Bathymetry (depth) in meter
            A_bathymetry = (A_bathymetry >= mask["depth"]).astype(int)  # (boolean)
        # Irrelevant parameters
        A_slope = 1
        A_notPopulated = 1

    # Masking matrix for the suitable sites (pixels)
    A_mask = (A_suitability_pa * A_suitability_lu * A_slope * A_notPopulated * A_bathymetry).astype(float)

    del A_suitability_lu, A_suitability_pa, A_slope, A_notPopulated, A_bathymetry

    # Calculate masked FLH
    FLH = hdf5storage.read('FLH', paths[tech]["FLH"])
    FLH_mask = FLH * A_mask
    FLH_mask[FLH_mask == 0] = np.nan

    # Save HDF5 Files
    hdf5storage.writes({'A_mask': A_mask}, paths[tech]["mask"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths[tech]["mask"])
    hdf5storage.writes({'FLH_mask': FLH_mask}, paths[tech]["FLH_mask"], store_python_metadata=True,
                       matlab_compatible=True)
    print("files saved: " + paths[tech]["FLH_mask"])

    # Save GEOTIFF files
    if param["savetiff"]:
        GeoRef = param["GeoRef"]
        array2raster(changeExt2tif(paths[tech]["mask"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelHeight"],
                     A_mask)
        print("files saved:" + changeExt2tif(paths[tech]["mask"]))

        array2raster(changeExt2tif(paths[tech]["FLH_mask"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelHeight"],
                     FLH_mask)
        print("files saved:" + changeExt2tif(paths[tech]["FLH_mask"]))

    timecheck('End')


def weighting(paths, param, tech):
    timecheck('Start')
    weight = param[tech]["weight"]
    Crd_all = param["Crd_all"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    res_desired = param["res_desired"]
    GeoRef = param["GeoRef"]

    if tech == 'PV':
        # Ground Cover Ratio - defines spacing between PV arrays
        A_GCR = calc_gcr(Crd_all, m_high, n_high, res_desired, weight["GCR"])
    else:
        A_GCR = 1

    with rasterio.open(paths["PA"]) as src:
        A_protect = src.read(1)
        A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified

    # Calculate availability based on protection categories
    A_availability_pa = changem(A_protect, weight["pa_availability"], param["protected_areas"]["type"]).astype(float)

    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
        A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified

    # Calculate availability based on landuse types
    A_availability_lu = changem(A_lu, weight["lu_availability"], param["landuse"]["type"]).astype(float)

    # Calculate availability
    A_availability = np.minimum(A_availability_pa, A_availability_lu)
    del A_availability_pa, A_availability_lu

    # Read available areas
    A_area = hdf5storage.read('A_area', paths["AREA"])

    # Weighting matrix for the power output (technical potential) in MWp
    A_weight = A_area * A_availability * A_GCR * weight["power_density"] * weight["f_performance"]

    # Calculate weighted FLH in MWh
    FLH = hdf5storage.read('FLH', paths[tech]["FLH"])
    FLH_weight = FLH * A_weight

    # Save HDF5 Files
    hdf5storage.writes({'A_weight': A_weight}, paths[tech]["weight"], store_python_metadata=True,
                       matlab_compatible=True)
    print("files saved: " + paths[tech]["weight"])
    hdf5storage.writes({'FLH_weight': FLH_weight}, paths[tech]["FLH_weight"], store_python_metadata=True,
                       matlab_compatible=True)
    print("files saved: " + paths[tech]["FLH_weight"])

    # Save GEOTIFF files
    if param["savetiff"]:
        array2raster(changeExt2tif(paths[tech]["weight"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelHeight"],
                     A_weight)
        print("files saved:" + changeExt2tif(paths[tech]["weight"]))

        array2raster(changeExt2tif(paths[tech]["FLH_weight"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelHeight"],
                     FLH_weight)
        print("files saved:" + changeExt2tif(paths[tech]["FLH_weight"]))
    timecheck('End')


def reporting(paths, param, tech):
    timecheck('Start')
    # read FLH, masking, area, and weighting matrix
    FLH = hdf5storage.read('FLH', paths[tech]["FLH"])
    A_mask = hdf5storage.read('A_mask', paths[tech]["mask"])
    A_weight = hdf5storage.read('A_weight', paths[tech]["weight"])
    A_area = hdf5storage.read('A_area', paths["AREA"])
    density = param[tech]["weight"]["power_density"]

    # Check if land or see
    if tech in ['PV', 'CSP', 'WindOn']:
        location = "land"
    elif tech in ['WindOff']:
        location = "sea"

    # Initialize region masking parameters
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    nRegions = param["nRegions_sub"]
    regions_shp = param["regions_sub"]

    # Initialize regions list of sorted FLH, FLH_M, and FLH_W
    sorted_FLH_list = {}

    # Define sampling for sorted lists
    sampling = param["report_sampling"]

    # Initialize dataframe
    regions = pd.DataFrame(0, index=range(0, nRegions),
                           columns=['Region', 'Available', 'Available_Masked', 'Available_Area_km2', 'FLH_Mean',
                                    'FLH_Median',
                                    'FLH_Max', 'FLH_Min', 'FLH_Mean_Masked', 'FLH_Median_Masked', 'FLH_Max_Masked',
                                    'FLH_Min_Masked', 'FLH_Std_Masked', 'Power_Potential_GW',
                                    'Power_Potential_Weighted_GW',
                                    'Energy_Potential_TWh', 'Energy_Potential_Weighted_TWh',
                                    'Energy_Potential_Weighted_Masked_TWh'])
    status = 0
    # Loop over each region
    for reg in range(0, nRegions):
        # Display Progress
        status += 1
        display_progress('Reporting ', (nRegions, status))

        # Get name of region
        regions.loc[reg, "Region"] = regions_shp.loc[reg]["NAME_SHORT"] + "_" + location

        # Compute region_mask
        A_region_extended = calc_region(regions_shp.loc[reg], Crd_all, res_desired, GeoRef)

        # Sum available : available pixels
        available = np.sum(A_region_extended)
        regions.loc[reg, "Available"] = int(available)

        # Sum availabe_masked : available pixels after masking
        A_masked = A_region_extended * A_mask
        available_masked = np.nansum(A_masked)
        regions.loc[reg, "Available_Masked"] = int(available_masked)

        # Interrupt reporting of region if no available pixels
        if int(available_masked) == 0:
            regions.drop([reg], axis=0, inplace=True)
            continue

        # Interrupt reporting of region already reported (may occur due to discrepancy in borders)
        if regions.loc[reg, "Region"] in regions.loc[:reg - 1, "Region"].to_list():
            ind_prev = regions.loc[regions["Region"] == regions.loc[reg, "Region"]].index[0]
            if regions.loc[ind_prev, "Available_Masked"] > int(available_masked):
                regions.drop([reg], axis=0, inplace=True)
                continue
            else:
                regions.drop([ind_prev], axis=0, inplace=True)

        # Sum area: available area in km2
        A_area_region = A_region_extended * A_area
        Total_area = np.nansum(A_area_region) / (10 ** 6)
        regions.loc[reg, "Available_Area_km2"] = Total_area

        # Stats for FLH
        FLH_region = A_region_extended * FLH
        FLH_region[FLH_region == 0] = np.nan
        regions.loc[reg, "FLH_Mean"] = np.nanmean(FLH_region)
        regions.loc[reg, "FLH_Median"] = np.nanmedian(FLH_region)
        regions.loc[reg, "FLH_Max"] = np.nanmax(FLH_region)
        regions.loc[reg, "FLH_Min"] = np.nanmin(FLH_region)
        regions.loc[reg, "FLH_Std"] = np.nanstd(FLH_region)

        # Stats for FLH_masked
        FLH_region_masked = A_masked * FLH_region
        FLH_region_masked[FLH_region_masked == 0] = np.nan
        if int(np.nansum(FLH_region_masked)) == 0:
            continue
        regions.loc[reg, "FLH_Mean_Masked"] = np.nanmean(FLH_region_masked)
        regions.loc[reg, "FLH_Median_Masked"] = np.nanmedian(FLH_region_masked)
        regions.loc[reg, "FLH_Max_Masked"] = np.nanmax(FLH_region_masked)
        regions.loc[reg, "FLH_Min_Masked"] = np.nanmin(FLH_region_masked)
        regions.loc[reg, "FLH_Std_Masked"] = np.nanstd(FLH_region_masked)

        # Power Potential
        A_P_potential = A_area_region * density
        power_potential = np.nansum(A_P_potential)
        regions.loc[reg, "Power_Potential_GW"] = power_potential / (10 ** 3)

        # Power Potential after weighting
        A_P_W_potential = A_region_extended * A_weight
        power_potential_weighted = np.nansum(A_P_W_potential)
        regions.loc[reg, "Power_Potential_Weighted_GW"] = power_potential_weighted / (10 ** 3)

        # Energy Potential
        A_E_potential = A_P_potential * FLH_region
        energy_potential = np.nansum(A_E_potential)
        regions.loc[reg, "Energy_Potential_TWh"] = energy_potential / (10 ** 6)

        # Energy Potential after weighting
        A_E_W_potential = FLH_region * A_weight
        energy_potential_weighted = np.nansum(A_E_W_potential)
        regions.loc[reg, "Energy_Potential_Weighted_TWh"] = energy_potential_weighted / (10 ** 6)

        # Energy Potential After weighting and masking
        A_E_W_M_potential = A_E_W_potential * A_masked
        energy_potential_weighted_masked = np.nansum(A_E_W_M_potential)
        regions.loc[reg, "Energy_Potential_Weighted_Masked_TWh"] = energy_potential_weighted_masked / (10 ** 6)

        sort = {}
        # Sorted FLH Sampling
        sorted_sampled_FLH = sampled_sorting(FLH_region[~np.isnan(FLH_region)], sampling)
        sort["FLH"] = sorted_sampled_FLH

        # Sorted FLH Sampling after masking
        sorted_sampled_FLH_masked = sampled_sorting(FLH_region_masked[~np.isnan(FLH_region_masked)], sampling)
        sort["FLH_M"] = sorted_sampled_FLH_masked

        # Sorted FLH Sampling after masking and wieghting
        FLH_region_masked_weighted = FLH_region_masked * A_weight
        FLH_region_masked_weighted[FLH_region_masked_weighted == 0] = np.nan

        sorted_sampled_FLH_masked_weighted = \
            sampled_sorting(FLH_region_masked_weighted[~np.isnan(FLH_region_masked_weighted)], sampling)
        sort["FLH_M_W"] = sorted_sampled_FLH_masked_weighted

        sorted_FLH_list[regions.loc[reg, "Region"]] = sort

    # Export the dataframe as CSV
    regions.to_csv(paths[tech]["Region_Stats"], sep=';', decimal=',', index=True)
    print("files saved: " + paths[tech]["Region_Stats"])

    # Save Sorted lists to .mat file
    for reg in sorted_FLH_list.keys():
        hdf5storage.writes({reg + '/FLH': sorted_FLH_list[reg]["FLH"],
                            reg + '/FLH_masked': sorted_FLH_list[reg]["FLH_M"],
                            reg + '/FLH_masked_weighted': sorted_FLH_list[reg]["FLH_M_W"]
                            }, paths[tech]["Sorted_FLH"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths[tech]["Sorted_FLH"])
    timecheck('End')


def find_locations_quantiles(paths, param, tech):
    timecheck('Start')
    FLH_mask = hdf5storage.read('FLH_mask', paths[tech]["FLH_mask"])
    quantiles = param["quantiles"]
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    # Select only indices in the report
    filter = pd.read_csv(paths[tech]["Region_Stats"], sep=';', decimal=',', index_col=0).index
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

        FLH_reg = A_region * FLH_mask[Ind[reg, 2] - 1:Ind[reg, 0], Ind[reg, 3] - 1:Ind[reg, 1]]
        FLH_reg[FLH_reg == 0] = np.nan
        X = FLH_reg.flatten(order='F')
        I_old = np.argsort(X)

        # Escape loop if intersection only yields NaN
        if np.isnan(X).all():
            # do nothing
            continue

        q_rank = 0
        for q in quantiles:
            list_names.append(regions_shp["NAME_SHORT"].loc[reg])
            list_quantiles.append('q' + str(q))
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

    reg_ind = np.reshape(reg_ind, (-1, 2), 'C').astype(int)
    reg_ind = (reg_ind[:, 0], reg_ind[:, 1])
    param[tech]["Ind_points"] = reg_ind
    param[tech]["Crd_points"] = crd_exact_points(reg_ind, Crd_all, res_desired)
    param[tech]["Crd_points"] = (param[tech]["Crd_points"][0], param[tech]["Crd_points"][1], list_names, list_quantiles)

    # Format point locations
    points = [(param[tech]["Crd_points"][1][i], param[tech]["Crd_points"][0][i]) for i in
              range(0, len(param[tech]["Crd_points"][0]))]

    # Create Shapefile
    schema = {'geometry': 'Point', 'properties': {'NAME_SHORT': 'str', 'quantile': 'str'}}
    with fiona.open(paths[tech]["Locations"], 'w', 'ESRI Shapefile', schema) as c:
        c.writerecords([{'geometry': mapping(Point(points[i])),
                         'properties': {'NAME_SHORT': list_names[i], 'quantile': list_quantiles[i]}} for i in
                        range(0, len(points))])
    hdf5storage.writes({'Ind_points': param[tech]["Ind_points"]}, paths[tech]["Locations"][:-4] + '_Ind.mat',
                       store_python_metadata=True, matlab_compatible=True)
    hdf5storage.writes({'Crd_points': param[tech]["Crd_points"]}, paths[tech]["Locations"][:-4] + '_Crd.mat',
                       store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths[tech]["Locations"])
    timecheck('End')


def generate_time_series(paths, param, tech):
    timecheck('Start')


    nproc = param["nproc"]
    CPU_limit = np.full((1, nproc), param["CPU_limit"])
    param[tech]["Crd_points"] = hdf5storage.read('Crd_points', paths[tech]["Locations"][:-4] + '_Crd.mat')
    param[tech]["Ind_points"] = hdf5storage.read('Ind_points', paths[tech]["Locations"][:-4] + '_Ind.mat')
    list_names = param[tech]["Crd_points"][2]
    list_quantiles = param[tech]["Crd_points"][3]
    if tech in ['PV', 'CSP']:
        res_weather = param["res_weather"]
        Crd_all = param["Crd_all"]
        Ind = ind_merra(Crd_all, Crd_all, res_weather)[0]

    list_hours = np.arange(0, 8760)
    param["status_bar_limit"] = list_hours[-1]

    # Obtain weather and correction matrices
    param["Ind_nz"] = param[tech]["Ind_points"]
    merraData, rasterData = get_merra_raster_Data(paths, param, tech)

    if tech in ['PV', 'CSP']:

        day_filter = np.nonzero(merraData["CLEARNESS"][Ind[2] - 1:Ind[0], Ind[3] - 1:Ind[1], :].sum(axis=(0, 1)))
        list_hours = np.arange(0, 8760)
        if nproc == 1:
            param["status_bar_limit"] = list_hours[-1]
            results = calc_TS_solar(list_hours[day_filter], [param, tech, rasterData, merraData])
        else:
            list_hours = np.array_split(list_hours[day_filter], nproc)
            param["status_bar_limit"] = list_hours[0][-1]
            results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
                calc_TS_solar, product(list_hours, [[param, tech, rasterData, merraData]]))

    elif tech in ['WindOn', 'WindOff']:

        list_hours = np.array_split(np.arange(0, 8760), nproc)
        param["status_bar_limit"] = list_hours[0][-1]
        results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(calc_TS_wind,
                                                                                           product(list_hours, [
                                                                                               [param, tech, rasterData,
                                                                                                merraData]]))
    print('\n')

    # Collecting results
    TS = np.zeros((len(param[tech]["Ind_points"][0]), 8760))
    if nproc > 1:
        for p in range(len(results)):
            TS = TS + results[p]
    else:
        TS = results

    # Restructuring results
    tuples = list(zip(list_names, list_quantiles))
    column_names = pd.MultiIndex.from_tuples(tuples, names=['NAME_SHORT', 'Quantile'])
    results = pd.DataFrame(TS.transpose(), columns=column_names)
    results.to_csv(paths[tech]["TS"], sep=';', decimal=',')
    print("files saved: " + paths[tech]["TS"])
    timecheck('End')


def regression_coefficients(paths, param, tech):
    """
    Solves the following optimization problem:

    Express a given model timeseries provided by EMHIRES as a combination timeseries
    for different Hub-Heights/orientations and Quantiles, while constraining the total sum of
    the obtained TS to the FLH given by IRENA

    :param paths: Dictionary including the paths to the shapefile of the globally protected areas, to the landuse raster of the scope, and to the output path PA.
    :type paths: dict
    :param param: Dictionary including the dictionary of regression parameters and year.
    :type param: dict
    :param tech: Name of the technology used for calculations
    :type tech: str
    :return:
        Copy the regression parameters Irena FLH and EMHIRES TS under Regression_Outputs folder,

        Save the regression coefficients, and the result Time-series in a .csv file under Regression_Outputs folder
    :rtype: None
    :raise Missing Data: No Time-series present for technology tech
    :raise Missing Data for Setting: Missing Time-series with desired settings (hub-heights/orientations)
    """

    timecheck('Start')
    year = str(param["year"])
    
    try:
        combinations = combinations_for_regression(paths, param, tech)
    except UserWarning:
        timecheck('End')
        return

    # Display the combinations of settings to be used
    if tech in ['WindOn', 'WindOff']:
        print("Combinations of hub heights to be used for the regression: ", combinations)
    elif tech in ['PV']:
        print("Orientations to be used for the regression: ", combinations)
    
    # Create FLH file for regression
    if not os.path.isfile(paths["FLH_regression"]):
        clean_FLH_regression(param, paths)

    # Create TS file for regression
    if not os.path.isfile(paths[tech]["TS_regression"]):
        clean_TS_regression(param, paths, tech)

    FLH, TS_reg = check_regression_model(paths, tech)

    param["FLH_regression"] = FLH
    param["TS_regression"] = TS_reg

    # Find intersection between FLH and shapefile subregions 
    list_regions = param["regions_sub"]["NAME_SHORT"].values.tolist()
    list_regions = sorted(list(set(list_regions).intersection(set(FLH.index))))

    # Summary Variables
    summary = None
    nodata = ''
    no_sol_high = ''
    no_sol_low = ''
    solution = ''

    # loop over all combinations and regions
    for settings in combinations:
        status = 0
        print("Regions under study : ", list_regions)
        for reg in list_regions:
            # Show progress of the simulation
            status = status + 1
            display_progress('Regression Coefficients ' + tech + ' ' + param["subregions_name"], (len(list_regions), status))
        
            region_data = regmodel_load_data(paths, param, tech, settings, reg)
        
            # Skip regions not present in the generated TS
            if region_data is None:
                nodata = nodata + reg + ', '
                continue
                
            settings_sorted = region_data[None]['s'][None].tolist()
        
            if region_data[None]["IRENA_best_worst"] == (True, True):
        
                # create model instance
                solver = SolverFactory(param["regression"]["solver"])
                model = pyomo_regression_model()
                regression = model.create_instance(region_data)
        
                # solve model and return results
                solver.solve(regression)
        
                # Retrieve results
                r = np.zeros((len(param["quantiles"]), len(settings_sorted)))
                c = 0
                for q in param["quantiles"]:
                    p = 0
                    for s in settings_sorted:
                        r[c, p] = pyo.value(regression.coef[s, q])
                        p += 1
                    c += 1
                r[r < 10 ** (-5)] = 0
                solution = solution + reg + ', '
        
            elif region_data[None]["IRENA_best_worst"] == (False, True):
                # Select best TS (highest height, highest quantile)
                r = np.full((len(param["quantiles"]), len(settings_sorted)), 0)
                r[0, 0] = 1
                no_sol_high = no_sol_high + reg + ', '
        
            elif region_data[None]["IRENA_best_worst"] == (True, False):
                # Select worst TS (lowest height, lowest quantile)
                r = np.full((len(param["quantiles"]), len(settings_sorted)), 0)
                r[-1, -1] = 1
                no_sol_low = no_sol_low + reg + ', '
        
            else:
                r = np.full((len(param["quantiles"]), len(settings_sorted)), 0)
        
            if settings_sorted != [0]:
                result = pd.DataFrame(r, param["quantiles"], (reg + "_" + str(s) for s in settings_sorted))
            else:
                result = pd.DataFrame(r, param["quantiles"], [reg])
        
            if summary is None:
                summary = result
            else:
                summary = pd.concat([summary, result], axis=1)
        
        # Print Regression Summary
        if solution != '':
            print("\nA solution was found for the following regions: " + solution.rstrip(', '))
        if no_sol_low != '':
            print("\nNo Solution was found for the following regions because they are too high: " + no_sol_low.rstrip(', '))
        if no_sol_high != '':
            print(
                "\nNo Solution was found for the following regions because they are too low: " + no_sol_high.rstrip(', '))
        if nodata != '':
            print("\nNo data was available for the following regions: " + nodata.rstrip(', '))
        
        st = ''
        for setting in settings_sorted:
            st = st + str(setting) + '_'
            
        summary.to_csv(paths[tech]["Regression_coefficients"] + st + year + '.csv', sep=';', decimal=',')
        print("\nfiles saved: " + paths[tech]["Regression_coefficients"] + st + year + '.csv')

    timecheck('End')
    
    
def generate_stratified_timeseries(paths, param, tech):
    '''
    description
    '''
    timecheck('Start')
    modes = param["modes"]
    subregions = param["subregions_name"]
    year = str(param["year"])
    
    try:
        settings_existing, inputfiles = combinations_for_stratified_timeseries(paths, param, tech)
    except UserWarning:
        timecheck('End')
        return
    
    # Display the combinations of settings to be used
    if tech in ['WindOn', 'WindOff']:
        print("Combinations of hub heights to be used for the stratified time series: ", settings_existing)
    elif tech in ['PV']:
        print("Orientations to be used for the stratified time series: ", settings_existing)
        
    for tag, combo in param["combo"][tech].items():
        if combo == []:
            combo = list(set([item for sublist in param["combo"][tech].values() for item in sublist]))
        ind = settings_existing.index(sorted(combo))
        coef = pd.read_csv(inputfiles[ind], sep=';', decimal=',', index_col=[0])
        
        # Extract names of regions
        regions = sorted(list(set([col.split('_')[0] for col in coef.columns])))
        
        # Load the TS files
        TS_files = {}
        for setting in combo:
            setting_path = paths["regional_analysis"] + subregions + '_' + tech + '_' + str(setting) + '_TS_' + year + '.csv'
            TS_files[setting] = pd.read_csv(setting_path, sep=';', decimal=',', header=[0,1], index_col=[0])
        
        # Loop over modes and regions
        TS_df = pd.DataFrame(index=range(8760), dtype='float16')
        for mode_tag, quantiles in modes.items():
            for reg in regions:
                col_name = reg + '_' + tech + '_' + tag + '_' + mode_tag
                TS_df[col_name] = np.zeros((8760, 1))
                filter_reg = [col for col in coef if col.startswith(reg)]
                for setting in combo:
                    sum_quantiles = coef.loc[quantiles, filter_reg].sum().sum()
                    for quantile in quantiles:
                        if sum_quantiles:
                            TS_df[col_name] = TS_df[col_name] + TS_files[setting][reg, 'q'+str(quantile)] * coef.loc[quantile, reg+'_'+str(setting)] / sum_quantiles
                        else:
                            TS_df[col_name] = TS_df[col_name] + TS_files[setting][reg, 'q'+str(quantile)] / len(quantiles) / len(combo)
                        
        st = ''
        for setting in combo:
            st = st + str(setting) + '_'
            
        TS_df.to_csv(paths[tech]["Regression_TS"] + st + year + '.csv', sep=';', decimal=',')
        print("File Saved: " + paths[tech]["Regression_TS"] + st + year + '.csv')
    timecheck('End')


if __name__ == '__main__':

    paths, param = initialization()
    generate_weather_files(paths, param)
    clean_weather_data(paths, param)
    generate_landsea(paths, param)  # Land and Sea
    generate_subregions(paths, param)  # Subregions
    generate_area(paths, param)
    generate_landuse(paths, param)  # Landuse
    generate_bathymetry(paths, param)  # Bathymetry
    generate_topography(paths, param)  # Topography
    generate_slope(paths, param)  # Slope
    generate_population(paths, param)  # Population
    generate_protected_areas(paths, param)  # Protected areas
    generate_buffered_population(paths, param)  # Buffered Population
    generate_area(paths, param)
    
    for tech in param["technology"]:
        print("Tech: " + tech)
        calculate_FLH(paths, param, tech)
        masking(paths, param, tech)
        weighting(paths, param, tech)
        reporting(paths, param, tech)
        find_locations_quantiles(paths, param, tech)
        generate_time_series(paths, param, tech)

    for tech in param["technology"]:
        print("Tech: " + tech)
        regression_coefficients(paths, param, tech)
        generate_stratified_timeseries(paths, param, tech)
