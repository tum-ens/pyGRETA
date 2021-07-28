from .correction_functions import clean_weather_data
from .spatial_functions import *
#import cython
import urllib.request
import os

def downloadGWA(paths, param):
    local_file = paths["GWA_global"]
    if not os.path.isfile(local_file):
        print('Download GWA:' + local_file)
        remote_url = 'https://globalwindatlas.info/api/gis/country/' + param["country_code"] + '/wind-speed/50'
        urllib.request.urlretrieve(remote_url, local_file)


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
    generate_array_coordinates(paths, param)
    generate_landsea(paths, param)  # Land and Sea
    generate_subregions(paths, param)  # Subregions
    generate_area(paths, param)  # Area Gradient
    generate_landuse(paths, param)  # Landuse
    # generate_bathymetry(paths, param)  # Bathymetry
    generate_topography(paths, param)  # Topography
    generate_slope(paths, param)  # Slope
    generate_protected_areas(paths, param)  # Protected areas
    # generate_livestock(paths,param)
    # generate_settlements(paths, param)
    # generate_osm(paths, param)
    # generate_population(paths, param)  # Population #not used anywhere?
   

def generate_buffered_maps(paths, param):

    # if "WindOn" in param["technology"]:
       # generate_buffered_population(paths, param)

    generate_buffered_population(paths, param)
    generate_buffered_water(paths, param)
    generate_buffered_wetland(paths, param)
    generate_buffered_snow(paths, param)
    generate_airports(paths, param)
    generate_country_boarders(paths, param)
    generate_buffered_protected_areas(paths, param)
    


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
    for l in range(gwa_cols):
        x_gwa[:,l] = GWA_speed.xy(0,l,offset='center')[0]
    for k in range(gwa_rows):
        y_gwa[k,:] = GWA_speed.xy(k,0,offset='center')[1]
    
    hdf5storage.writes({"GWA_X": x_gwa}, paths["GWA_X"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths["GWA_X"])
    hdf5storage.writes({"GWA_Y": y_gwa}, paths["GWA_Y"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths["GWA_Y"])
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
    # m_high = param["m_high"]
    # n_high = param["n_high"]
    # w_new = np.full((m_high, n_high), np.nan)
    with rasterio.open(paths["LU_global"]) as src:
        w = src.read(1, window=windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]), slice(Ind[3] - 1, Ind[1])))
        w = np.flipud(w)
    print (w.shape)
    # w = adjust_resolution(w, param["res_landuse"], param["res_desired"], "category")
    # w_new = recalc_lu_resolution(w, param["res_landuse"], param["res_desired"], w_new, param["m_low"], param["n_low"])
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
    Ind = ind_global(Crd_all, param["res_topography"])[0]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["Bathym_global"]) as src:
        A_BATH = src.read(1)
    #A_BATH = resizem(A_BATH, 180 * 240, 360 * 240)
    A_BATH = adjust_resolution(A_BATH, param["res_bathymetry"], param["res_topography"], "mean")
    A_BATH = np.flipud(A_BATH[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]])
    print (A_BATH.shape)
    A_BATH = recalc_topo_resolution(A_BATH, param["res_topography"], param["res_desired"])
    
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
            feat.SetField("Raster", int(protection_type[pt]))
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
  
  
def generate_buffered_protected_areas(paths, param):
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
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["PA"]) as src:
        A_pa = src.read(1)
    A_pa = np.flipud(A_pa).astype(int)
    A_pa = (A_pa>0) & (A_pa<6)
    
    buffer_pixel_amount1 = param["PV"]["mask"]["pa_buffer_pixel_amount"]
    kernel1 = np.tri(2 * buffer_pixel_amount1 + 1, 2 * buffer_pixel_amount1 + 1, buffer_pixel_amount1).astype(int)
    kernel1 = kernel1 * kernel1.T * np.flipud(kernel1) * np.fliplr(kernel1)
    A_pa_buffered1 = maximum_filter(A_pa, footprint=kernel1, mode="constant", cval=0)
    A_notProtected1 = (~A_pa_buffered1).astype(int)

    array2raster(paths["PV_PA_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notProtected1)
    print("files saved: " + paths["PV_PA_BUFFER"])
    create_json(paths["PV_PA_BUFFER"], param, ["region_name", "protected_areas", "PV", "Crd_all", "res_desired", "GeoRef"], paths, ["PA", "PV_PA_BUFFER"])
    

    buffer_pixel_amount = param["WindOn"]["mask"]["pa_buffer_pixel_amount"]
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_pa_buffered = maximum_filter(A_pa, footprint=kernel, mode="constant", cval=0)
    A_notProtected = (~A_pa_buffered).astype(int)
    
    array2raster(paths["WINDON_PA_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notProtected)
    print("files saved: " + paths["WINDON_PA_BUFFER"])
    create_json(paths["WINDON_PA_BUFFER"], param, ["region_name", "protected_areas", "WindOn", "Crd_all", "res_desired", "GeoRef"], paths, ["PA", "WINDON_PA_BUFFER"])
    
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
    buffer_pixel_amount = param["WindOn"]["mask"]["urban_buffer_pixel_amount"]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
    A_lu = np.flipud(A_lu).astype(int)
    A_lu = A_lu == 190  # Land use type for Urban and built-up
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_lu_buffered = maximum_filter(A_lu, footprint=kernel, mode="constant", cval=0)
    A_notPopulated = (~A_lu_buffered).astype(int)

    array2raster(paths["POP_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notPopulated)
    print("files saved: " + paths["POP_BUFFER"])
    create_json(paths["POP_BUFFER"], param, ["region_name", "landuse", "WindOn", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "POP_BUFFER"])
    timecheck("End")


def generate_buffered_water(paths, param):
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
    buffer_pixel_amount = param["landuse"]["water_buffer"]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
    A_lu = np.flipud(A_lu).astype(int)
    A_lu = A_lu == 210 # Land use type for water
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_lu_buffered = maximum_filter(A_lu, footprint=kernel, mode="constant", cval=0)
    A_notWater = (~A_lu_buffered).astype(int)

    array2raster(paths["WATER_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notWater)
    print("files saved: " + paths["WATER_BUFFER"])
    create_json(paths["WATER_BUFFER"], param, ["region_name", "landuse", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "WATER_BUFFER"])
    timecheck("End")


def generate_buffered_wetland(paths, param):
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
    buffer_pixel_amount = param["landuse"]["wetland_buffer"]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
    A_lu = np.flipud(A_lu).astype(int)
    A_lu = ((A_lu == 160) | (A_lu == 170) | (A_lu == 180)) # Land use type for wetland
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_lu_buffered = maximum_filter(A_lu, footprint=kernel, mode="constant", cval=0)
    A_notWetland = (~A_lu_buffered).astype(int)

    array2raster(paths["WETLAND_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notWetland)
    print("files saved: " + paths["WETLAND_BUFFER"])
    create_json(paths["WETLAND_BUFFER"], param, ["region_name", "landuse", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "WETLAND_BUFFER"])
    timecheck("End")
   
   
def generate_buffered_snow(paths, param):
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
    buffer_pixel_amount = param["landuse"]["snow_buffer"]
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
    A_lu = np.flipud(A_lu).astype(int)
    A_lu = A_lu == 220 # Land use type for snow
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_lu_buffered = maximum_filter(A_lu, footprint=kernel, mode="constant", cval=0)
    A_notSnow = (~A_lu_buffered).astype(int)

    array2raster(paths["SNOW_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notSnow)
    print("files saved: " + paths["SNOW_BUFFER"])
    create_json(paths["SNOW_BUFFER"], param, ["region_name", "landuse", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "SNOW_BUFFER"])
    timecheck("End")

   
def generate_airports(paths,param):
    timecheck("Start")
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    countries_shp = param["regions_land"]
    nCountries = param["nRegions_land"]
    
     # Load Airports dictionary
    airports_list = pd.read_csv(paths["Airports"], index_col = ["iso_country"],usecols=["iso_country","name","latitude_deg","longitude_deg","type"])
    
     # Load IRENA dictionary
    IRENA_dict = pd.read_csv(paths["IRENA_dict"], sep=";",index_col = ["Countries shapefile"],usecols=["Countries shapefile","Countries Alpha-2 code"])
    
    airports = []
    for reg in range(0, nCountries):
        alpha2code = IRENA_dict["Countries Alpha-2 code"][countries_shp.iloc[reg]["GID_0"]]
        #print (alpha2code)
        airports_filtered = airports_list[airports_list.index==alpha2code]
        airports_filtered = airports_filtered[(airports_filtered["type"]=="small_airport")|(airports_filtered["type"]=="medium_airport")|(airports_filtered["type"]=="large_airport")]
        # print (airports_filtered)
        airports.append(airports_filtered)
    airports = pd.concat(airports)
    print ("Airports are filtered \n")
    
    # Filter points outside spatial scope
    lat_max, lon_max, lat_min, lon_min = param["spatial_scope"][0]
    
    # Points inside the scope bounds
    airports = airports.loc[
        (lat_min <= airports["latitude_deg"]) & (lat_max >= airports["latitude_deg"]) & (lon_min <= airports["longitude_deg"]) & (lon_max >= airports["longitude_deg"])
    ].copy()
    
    with rasterio.open(paths["LAND"]) as src:
        A_land = src.read(1)
    A_land = np.flipud(A_land).astype(int)
    
    if not airports.empty:
        # Prepare input
        crd = (airports["latitude_deg"].to_numpy(), airports["longitude_deg"].to_numpy())
        ind = ind_exact_points(crd, Crd_all, res_desired)

        A_land[tuple(ind)]=100
        airport_raster = A_land == 100

        buffer_pixel_amount = param["WindOn"]["mask"]["airport_buffer_pixel_amount"]
        kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
        kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
        airport_raster = maximum_filter(airport_raster, footprint=kernel, mode="constant", cval=0)
        A_notAirport = (~airport_raster).astype(int)

        array2raster(paths["AIRPORTS"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notAirport)
        print("files saved: " + paths["AIRPORTS"])
        create_json(paths["AIRPORTS"], param, ["region_name", "landuse", "Biomass", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "AIRPORTS"])
    
    timecheck("End") 


def generate_country_boarders(paths,param):
    timecheck("Start")
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    countries_shp = param["regions_land"]
    nCountries = param["nRegions_land"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    
    A_countries_buffered = np.zeros([m_high, n_high]).astype(int)

    buffer_pixel_amount = param["landuse"]["boarder_buffer_pixel_amount"]
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
        
    for reg in range(0, nCountries):
        try:
            A_country_area = calc_region(countries_shp.loc[reg], Crd_all, res_desired, GeoRef)
            A_country_buffered = minimum_filter(A_country_area, footprint=kernel, mode="constant", cval=1)
            A_countries_buffered = A_countries_buffered + A_country_buffered 
        except:
            pass
    print (np.sum(A_countries_buffered))
    A_countries_buffered = A_countries_buffered > 0
    print (np.sum(A_countries_buffered))
    A_notBoarder = (A_countries_buffered).astype(int)
        
    
    array2raster(paths["BOARDERS"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notBoarder)
    print("files saved: " + paths["BOARDERS"])
    create_json(paths["BOARDERS"], param, ["region_name", "landuse", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "BOARDERS"])
    
    timecheck("End") 
    

def generate_roads(paths, param):
    #"fclass" ILIKE 'primary' OR "fclass" ILIKE 'secondary'
    timecheck("Start")
    
    # First we will open our raster image, to understand how we will want to rasterize our vector
    raster_ds = gdal.Open(paths["LU"], gdal.GA_ReadOnly)

    # Fetch number of rows and columns
    ncol = raster_ds.RasterXSize
    nrow = raster_ds.RasterYSize

    # Fetch projection and extent
    proj = raster_ds.GetProjectionRef()
    ext = raster_ds.GetGeoTransform()

    raster_ds = None
    shp_path = paths["OSM"]+"germany-roads.shp"
    # Open the dataset from the file
    dataset = ogr.Open(shp_path, 1)
    layer = dataset.GetLayerByIndex(0)

    # Add a new field
    if not field_exists("Raster", shp_path):
        new_field = ogr.FieldDefn("Raster", ogr.OFTInteger)
        layer.CreateField(new_field)

    for feat in layer:
        feat.SetField("Raster", 1)
        layer.SetFeature(feat)
        feat = None

    # Create a second (modified) layer
    outdriver = ogr.GetDriverByName("MEMORY")
    source = outdriver.CreateDataSource("memData")

    # Create the raster dataset
    memory_driver = gdal.GetDriverByName("GTiff")
    out_raster_ds = memory_driver.Create(paths["ROADS"], ncol, nrow, 1, gdal.GDT_Byte)

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
    create_json(paths["ROADS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["ROADS"])

    # Close dataset
    out_raster_ds = None
    print("files saved: " + paths["ROADS"])    
    timecheck("End")
    

def generate_railways(paths, param):
    #"fclass" ILIKE 'primary' OR "fclass" ILIKE 'secondary'
    timecheck("Start")
    
    # First we will open our raster image, to understand how we will want to rasterize our vector
    raster_ds = gdal.Open(paths["LU"], gdal.GA_ReadOnly)

    # Fetch number of rows and columns
    ncol = raster_ds.RasterXSize
    nrow = raster_ds.RasterYSize

    # Fetch projection and extent
    proj = raster_ds.GetProjectionRef()
    ext = raster_ds.GetGeoTransform()

    raster_ds = None
    shp_path = paths["OSM"]+"germany-railways.shp"
    # Open the dataset from the file
    dataset = ogr.Open(shp_path, 1)
    layer = dataset.GetLayerByIndex(0)

    # Add a new field
    if not field_exists("Raster", shp_path):
        new_field = ogr.FieldDefn("Raster", ogr.OFTInteger)
        layer.CreateField(new_field)

    for feat in layer:
        feat.SetField("Raster", 1)
        layer.SetFeature(feat)
        feat = None

    # Create a second (modified) layer
    outdriver = ogr.GetDriverByName("MEMORY")
    source = outdriver.CreateDataSource("memData")

    # Create the raster dataset
    memory_driver = gdal.GetDriverByName("GTiff")
    out_raster_ds = memory_driver.Create(paths["RAILS"], ncol, nrow, 1, gdal.GDT_Byte)

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
    create_json(paths["RAILS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["RAILS"])

    # Close dataset
    out_raster_ds = None
    print("files saved: " + paths["RAILS"])    
    timecheck("End")
    
    
def generate_osm_areas(paths, param):

    #"fclass" ILIKE 'commercial' OR "fclass" ILIKE 'industrial' OR "fclass" ILIKE 'quarry' OR "fclass" ILIKE 'military' OR "fclass" ILIKE 'park' OR "fclass" ILIKE 'recreation_ground'
    
    timecheck("Start")
    
    osm_areas = param["osm_areas"]
    # set up osm areas dictionary
    osm_type = dict(zip(osm_areas["Category"], osm_areas["type"]))
    
    # First we will open our raster image, to understand how we will want to rasterize our vector
    raster_ds = gdal.Open(paths["LU"], gdal.GA_ReadOnly)

    # Fetch number of rows and columns
    ncol = raster_ds.RasterXSize
    nrow = raster_ds.RasterYSize

    # Fetch projection and extent
    proj = raster_ds.GetProjectionRef()
    ext = raster_ds.GetGeoTransform()

    raster_ds = None
    shp_path = paths["OSM"]+"germany-landuse-all.shp"
    # Open the dataset from the file
    dataset = ogr.Open(shp_path, 1)
    layer = dataset.GetLayerByIndex(0)
    # layer.SetAttributeFilter("fclass" ILIKE 'commercial' OR "fclass" ILIKE 'industrial' OR "fclass" ILIKE 'quarry' OR "fclass" ILIKE 'military' OR "fclass" ILIKE 'park' OR "fclass" ILIKE 'recreation_ground'")

    # Add a new field
    if not field_exists("Raster", shp_path):
        new_field = ogr.FieldDefn("Raster", ogr.OFTInteger)
        layer.CreateField(new_field)
    
    # layer = layer.SetAttributeFilter('fclass = "commercial" OR fclass = "industrial" OR fclass = "quarry" OR fclass = "military" OR fclass = "park" OR fclass = "recreation_ground"')
    
    for feat in layer:
        pt = feat.GetField("fclass")
        if pt == "commercial" or pt == "industrial" or pt == "quarry" or pt == "military" or pt == "park" or pt == "recreation_ground":
            feat.SetField("Raster", int(osm_type[pt]))
        else:
            feat.SetField("Raster", int(0))
        layer.SetFeature(feat)
        feat = None

    # Create a second (modified) layer
    outdriver = ogr.GetDriverByName("MEMORY")
    source = outdriver.CreateDataSource("memData")

    # Create the raster dataset
    memory_driver = gdal.GetDriverByName("GTiff")
    out_raster_ds = memory_driver.Create(paths["OSM_AREAS"], ncol, nrow, 1, gdal.GDT_Byte)

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
    create_json(paths["OSM_AREAS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["OSM","OSM_AREAS"])

    # Close dataset
    out_raster_ds = None
    print("files saved: " + paths["OSM_AREAS"])    
    timecheck("End")
    

def generate_buffered_osm_areas(paths, param):
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
    GeoRef = param["GeoRef"]
    with rasterio.open(paths["OSM_AREAS"]) as src:
        A_osma = src.read(1)
    A_osma = np.flipud(A_osma).astype(int)
    
    #Commercial Areas
    A_com = A_osma == 1
    buffer_pixel_amount = 2
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_com_buffered = maximum_filter(A_com, footprint=kernel, mode="constant", cval=0)
    A_notCommercial = (~A_com_buffered).astype(int)
    array2raster(paths["OSM_COM_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notCommercial)
    print("files saved: " + paths["OSM_COM_BUFFER"])
    create_json(paths["OSM_COM_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["OSM", "OSM_COM_BUFFER"])
    
    #Industrical Areas
    A_ind = A_osma == 2
    buffer_pixel_amount = 1
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_ind_buffered = maximum_filter(A_ind, footprint=kernel, mode="constant", cval=0)
    A_notIndustrial = (~A_ind_buffered).astype(int)
    array2raster(paths["OSM_IND_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notIndustrial)
    print("files saved: " + paths["OSM_IND_BUFFER"])
    create_json(paths["OSM_IND_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["OSM", "OSM_IND_BUFFER"])
    
    #Mining Areas
    A_min = A_osma == 3
    buffer_pixel_amount = 1
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_min_buffered = maximum_filter(A_min, footprint=kernel, mode="constant", cval=0)
    A_notMine = (~A_min_buffered).astype(int)
    array2raster(paths["OSM_MINE_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notMine)
    print("files saved: " + paths["OSM_MINE_BUFFER"])
    create_json(paths["OSM_MINE_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["OSM", "OSM_MINE_BUFFER"])
   
    #Military Areas
    A_mil = A_osma == 4
    buffer_pixel_amount = 2
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_mil_buffered = maximum_filter(A_mil, footprint=kernel, mode="constant", cval=0)
    A_notMilitary = (~A_mil_buffered).astype(int)
    array2raster(paths["OSM_MIL_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notMilitary)
    print("files saved: " + paths["OSM_MIL_BUFFER"])
    create_json(paths["OSM_MIL_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["OSM", "OSM_MIL_BUFFER"])
    
    #Parks
    A_park = A_osma == 5
    #Wind onshore
    buffer_pixel_amount = 2
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_park_buffered = maximum_filter(A_park, footprint=kernel, mode="constant", cval=0)
    A_notPark = (~A_park_buffered).astype(int)
    array2raster(paths["WINDON_OSM_PARK_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notPark)
    print("files saved: " + paths["WINDON_OSM_PARK_BUFFER"])
    create_json(paths["WINDON_OSM_PARK_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["OSM", "WINDON_OSM_PARK_BUFFER"])
    #PV
    buffer_pixel_amount = 1
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_park_buffered = maximum_filter(A_park, footprint=kernel, mode="constant", cval=0)
    A_notPark = (~A_park_buffered).astype(int)
    array2raster(paths["PV_OSM_PARK_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notPark)
    print("files saved: " + paths["PV_OSM_PARK_BUFFER"])
    create_json(paths["PV_OSM_PARK_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["OSM", "PV_OSM_PARK_BUFFER"])
   
    #Recreational Areas
    A_rec = A_osma == 6
    buffer_pixel_amount = 1
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_rec_buffered = maximum_filter(A_rec, footprint=kernel, mode="constant", cval=0)
    A_notRecreation = (~A_rec_buffered).astype(int)
    array2raster(paths["OSM_REC_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notRecreation)
    print("files saved: " + paths["OSM_REC_BUFFER"])
    create_json(paths["OSM_REC_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["OSM", "OSM_REC_BUFFER"])
    
    timecheck("End")
    
    
def generate_settlements(paths, param):
    """
    This function reads the global map of settlements, and creates a raster out of it for the desired scope.
       See :mod:`config.py` for more information on the settlements map.

    :param paths: Dictionary including the paths to the global settlements raster *WSF_global* and to the output path *WSF*.
    :type paths: dict
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
    :type param: dict

    :return: The tif file for *WSF* is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """
    timecheck("Start")
    Crd_all = param["Crd_all"]
    Ind = ind_global(Crd_all, param["res_settlements"])[0]
    GeoRef = param["GeoRef"]

    with rasterio.open(paths["WSF_global"]) as src:
        A_wsf = src.read(1, window=windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]), slice(Ind[3] - 1, Ind[1]))).astype(bool)
        A_wsf = np.flipud(A_wsf)
    array2raster(paths["WSF"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_wsf)
    print("files saved: " + paths["WSF"])
    create_json(paths["WSF"], param, ["region_name", "Crd_all", "res_settlements", "res_desired", "GeoRef"], paths, ["WSF_global", "WSF"])
    
    # PV buffer
    buffer_pixel_amount = 1
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_pv_wsf = maximum_filter(A_wsf, footprint=kernel, mode="constant", cval=0)
    A_NotSettlement = (~A_pv_wsf).astype(int)
     
    array2raster(paths["PV_WSF_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_NotSettlement)
    print("files saved: " + paths["PV_WSF_BUFFER"])
    create_json(paths["PV_WSF_BUFFER"], param, ["region_name", "Crd_all", "res_settlements", "res_desired", "GeoRef"], paths, ["WSF_global", "PV_WSF_BUFFER"])
    
    # Onshore wind buffer
    buffer_pixel_amount = 4
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_wind_wsf = maximum_filter(A_wsf, footprint=kernel, mode="constant", cval=0)
    A_NotSettlement = (~A_wind_wsf).astype(int)
     
    array2raster(paths["WINDON_WSF_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_NotSettlement)
    print("files saved: " + paths["WINDON_WSF_BUFFER"])
    create_json(paths["WINDON_WSF_BUFFER"], param, ["region_name", "Crd_all", "res_settlements", "res_desired", "GeoRef"], paths, ["WSF_global", "WINDON_WSF_BUFFER"])
    
    timecheck("End")
    
    
def generate_HydroLakes(paths, param):
    #"fclass" ILIKE 'primary' OR "fclass" ILIKE 'secondary'
    timecheck("Start")
    
    # First we will open our raster image, to understand how we will want to rasterize our vector
    raster_ds = gdal.Open(paths["LU"], gdal.GA_ReadOnly)

    # Fetch number of rows and columns
    ncol = raster_ds.RasterXSize
    nrow = raster_ds.RasterYSize

    # Fetch projection and extent
    proj = raster_ds.GetProjectionRef()
    ext = raster_ds.GetGeoTransform()

    raster_ds = None
    shp_path = paths["HydroLakes"]
    # Open the dataset from the file
    dataset = ogr.Open(shp_path, 1)
    layer = dataset.GetLayerByIndex(0)

    # Add a new field
    if not field_exists("Raster", shp_path):
        new_field = ogr.FieldDefn("Raster", ogr.OFTInteger)
        layer.CreateField(new_field)

    for feat in layer:
        feat.SetField("Raster", 1)
        layer.SetFeature(feat)
        feat = None

    # Create a second (modified) layer
    outdriver = ogr.GetDriverByName("MEMORY")
    source = outdriver.CreateDataSource("memData")

    # Create the raster dataset
    memory_driver = gdal.GetDriverByName("GTiff")
    out_raster_ds = memory_driver.Create(paths["HYDROLAKES"], ncol, nrow, 1, gdal.GDT_Byte)

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
    create_json(paths["HYDROLAKES"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["HYDROLAKES"])

    # Close dataset
    out_raster_ds = None
    print("files saved: " + paths["HYDROLAKES"])
    
    # Create Buffer
    
    GeoRef = param["GeoRef"]

    with rasterio.open(paths["HYDROLAKES"]) as src:
        A_lake = src.read(1).astype(bool)
        A_lake = np.flipud(A_lake)
    buffer_pixel_amount = 1
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_lake = maximum_filter(A_lake, footprint=kernel, mode="constant", cval=0)
    A_NotLake = (~A_lake).astype(int)
     
    array2raster(paths["HYDROLAKES_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_NotLake)
    print("files saved: " + paths["HYDROLAKES_BUFFER"])
    create_json(paths["HYDROLAKES_BUFFER"], param, ["region_name", "Crd_all", "res_settlements", "res_desired", "GeoRef"], paths, ["HYDROLAKES", "HYDROLAKES_BUFFER"])
    
    timecheck("End")


def generate_HydroRivers(paths, param):
    #"fclass" ILIKE 'primary' OR "fclass" ILIKE 'secondary'
    timecheck("Start")
    
    # First we will open our raster image, to understand how we will want to rasterize our vector
    raster_ds = gdal.Open(paths["LU"], gdal.GA_ReadOnly)

    # Fetch number of rows and columns
    ncol = raster_ds.RasterXSize
    nrow = raster_ds.RasterYSize

    # Fetch projection and extent
    proj = raster_ds.GetProjectionRef()
    ext = raster_ds.GetGeoTransform()

    raster_ds = None
    shp_path = paths["HydroRivers"]
    # Open the dataset from the file
    dataset = ogr.Open(shp_path, 1)
    layer = dataset.GetLayerByIndex(0)

    # Add a new field
    if not field_exists("Raster", shp_path):
        new_field = ogr.FieldDefn("Raster", ogr.OFTInteger)
        layer.CreateField(new_field)

    for feat in layer:
        feat.SetField("Raster", 1)
        layer.SetFeature(feat)
        feat = None

    # Create a second (modified) layer
    outdriver = ogr.GetDriverByName("MEMORY")
    source = outdriver.CreateDataSource("memData")

    # Create the raster dataset
    memory_driver = gdal.GetDriverByName("GTiff")
    out_raster_ds = memory_driver.Create(paths["HYDRORIVERS"], ncol, nrow, 1, gdal.GDT_Byte)

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
    create_json(paths["HYDRORIVERS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["HYDROLAKES"])

    # Close dataset
    out_raster_ds = None
    print("files saved: " + paths["HYDRORIVERS"])
    
    # Create Buffer
    
    GeoRef = param["GeoRef"]

    with rasterio.open(paths["HYDRORIVERS"]) as src:
        A_Riv = src.read(1).astype(bool)
        A_Riv = np.flipud(A_Riv)
    buffer_pixel_amount = 1
    kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_Riv = maximum_filter(A_Riv, footprint=kernel, mode="constant", cval=0)
    A_NotRiver = (~A_Riv).astype(int)
     
    array2raster(paths["HYDRORIVERS_BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_NotRiver)
    print("files saved: " + paths["HYDRORIVERS_BUFFER"])
    create_json(paths["HYDRORIVERS_BUFFER"], param, ["region_name", "Crd_all", "res_settlements", "res_desired", "GeoRef"], paths, ["HYDRORIVERS", "HYDRORIVERS_BUFFER"])
    
    timecheck("End")