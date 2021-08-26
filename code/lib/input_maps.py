from . import correction_functions as cf
from . import spatial_functions as sf
from . import util as ul
from .log import logger
import multiprocessing as mp
import numpy as np
import pandas as pd
import urllib.request
import os
import rasterio
import h5netcdf
import hdf5storage
import datetime
import scipy.ndimage
import sys


def downloadGWA(paths, param):
    """
    This function downloads wind speed data from Global Wind Atlas (www.globalwindatlas.info) if it does not already exist

    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict

    :return: The wind data is saved directly in the desired paths.
    :rtype: None
    """
    if os.path.isfile(paths["GWA_global"]):
        logger.info('Skip')
    else:
        logger.info('Downlad GWA:' + paths["GWA_global"])

        remote_url = 'https://globalwindatlas.info/api/gis/country/' + param["country_code"] + '/wind-speed/50'
        urllib.request.urlretrieve(remote_url, paths["GWA_global"])


def generate_maps_for_scope(paths, param, multiprocessing):
    """
    This function calls the individual functions that generate the maps for the geographic scope.
    
    :param paths: Dictionary including the paths.
    :type paths: dict
    :param param: Dictionary including the user preferences.
    :type param: dict
    :param multiprocessing: Determines if multiprocessing is applied.
    :type param: bool
    
    :return: The maps are saved directly in the desired paths.
    :rtype: None
    """

    generate_land(paths, param)  # Subregions
    generate_weather_files(paths, param)  # MERRA Weather data

    if multiprocessing:
        processes = []
        processes.append(mp.Process(target=generate_sea, args=(paths, param)))
        processes.append(mp.Process(target=generate_area, args=(paths, param)))
        processes.append(mp.Process(target=generate_topography, args=(paths, param)))
        # processes.append(mp.Process(target=generate_bathymetry, args=(paths, param))) # ToDo: not tested
        processes.append(mp.Process(target=generate_landuse, args=(paths, param)))
        processes.append(mp.Process(target=generate_protected_areas, args=(paths, param)))
        processes.append(mp.Process(target=generate_airports, args=(paths, param)))
        processes.append(mp.Process(target=generate_country_boarders, args=(paths, param)))
        processes.append(mp.Process(target=generate_roads, args=(paths, param)))
        processes.append(mp.Process(target=generate_railways, args=(paths, param)))
        processes.append(mp.Process(target=generate_osm_areas, args=(paths, param)))
        processes.append(mp.Process(target=generate_settlements, args=(paths, param)))
        processes.append(mp.Process(target=generate_HydroLakes, args=(paths, param)))
        processes.append(mp.Process(target=generate_HydroRivers, args=(paths, param)))

        logger.debug('Starting processes')
        for p in processes:
            p.start()  # Start all single processes
        logger.debug('All processes started')

        for p in processes:
            p.join()  # Wait until all processes are finished
        logger.info('All processes finished')

    else:
        generate_sea(paths, param)  # Land and Sea
        generate_area(paths, param)  # Area Gradient
        generate_topography(paths, param)  # Topography
        # generate_bathymetry(paths, param)  # Bathymetry # ToDo: not tested
        generate_landuse(paths, param)  # Landuse
        generate_protected_areas(paths,param)
        generate_airports(paths, param)
        generate_country_boarders(paths, param)
        generate_roads(paths, param)
        generate_railways(paths, param)
        generate_osm_areas(paths, param)
        generate_settlements(paths, param)
        generate_HydroLakes(paths, param)
        generate_HydroRivers(paths, param)


def generate_land(paths, param):
    """
    This function reads the shapefile of the subregions within the scope, and creates a raster out of it.

    :param paths: Dictionary including the paths *SUB*, *LAND*, *EEZ*.
    :type paths: dict
    :param param: Dictionary including the geodataframe of the shapefile, the number of features, the coordinates of the bounding box of the spatial scope, and the number of rows and columns.
    :type param: dict

    :return: The tif file for *SUB* is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """

    if os.path.isfile(paths["LAND"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

        m_high = param["m_high"]
        n_high = param["n_high"]
        Crd_all = param["Crd_all"]
        res_desired = param["res_desired"]
        GeoRef = param["GeoRef"]
        nRegions_land = param["nRegions_land"]

        # Read shapefile of regions
        regions_shp = param["regions_land"]
        Crd_regions_land = param["Crd_regions_land"]
        Ind = sf.ind_merra(Crd_regions_land, Crd_all, res_desired)
        A_land = np.zeros((m_high, n_high))
        for reg in range(0, nRegions_land):
            logger.debug('Region: ' + str(reg))

            # Calculate A_region
            A_region = sf.calc_region(regions_shp.iloc[reg], Crd_regions_land[reg, :], res_desired, GeoRef)

            # Include A_region in A_sub
            A_land[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] = (
                A_land[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] + A_region
            )

        # Saving file
        sf.array2raster(paths["LAND"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_land)
        logger.info("files saved: " + paths["LAND"])
        ul.create_json(
            paths["LAND"], param, ["subregions_name", "m_high", "n_high", "Crd_all", "res_desired", "GeoRef", "nRegions_sea"], paths, ["Countries", "LAND"] # FIXME: replaced regions by Countries
        )

        logger.debug("End")


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
    # logger.info(paths["T2M"])
    if os.path.isfile(paths["T2M"]) and os.path.isfile(paths["W50M"]) and os.path.isfile(paths["CLEARNESS"]) and \
            os.path.isfile(paths["MERRA_XMIN"]) and os.path.isfile(paths["MERRA_XMAX"]) and \
            os.path.isfile(paths["MERRA_YMIN"]) and os.path.isfile(paths["MERRA_YMAX"]) and \
            os.path.isfile(paths["GWA_X"]) and os.path.isfile(paths["GWA_Y"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

        SWGDN = np.array([])
        SWTDN = np.array([])
        T2M = np.array([])
        U50M = np.array([])
        V50M = np.array([])

        status = 0
        start = datetime.date(param["year"], 1, 1)
        end = datetime.date(param["year"], 12, 31)
        delta = (end - start).days + 1

        for date in pd.date_range(start, end):
            # Show status bar
            status = status + 1
            sys.stdout.write("\r")
            sys.stdout.write(
                "Reading NetCDF files " + "[%-50s] %d%%" % ("=" * ((status * 50) // delta), (status * 100) // delta))
            sys.stdout.flush()

            # tomorrow = date + pd.Timedelta("1 day")
            if date.day == 29 and date.month == 2:
                continue    # Skip additional day of non leap year

            # Name and path of the NetCDF file to be read
            name = paths["MERRA_IN"] + "MERRA2_400.tavg1_2d_rad_Nx." + date.strftime("%Y%m%d") + ".nc4.nc4"
            name2 = paths["MERRA_IN"] + "MERRA2_400.tavg1_2d_slv_Nx." + date.strftime("%Y%m%d") + ".nc4.nc4"

            # Read NetCDF file, extract hourly tables
            with h5netcdf.File(name, "r") as f:
                # [time, lat 361, lon 576]
                swgdn = np.transpose(sf.subset(f["SWGDN"], param), [1, 2, 0])
                if SWGDN.size == 0:
                    SWGDN = swgdn
                else:
                    SWGDN = np.concatenate((SWGDN, swgdn), axis=2)

                swtdn = np.transpose(sf.subset(f["SWTDN"], param), [1, 2, 0])
                if SWTDN.size == 0:
                    SWTDN = swtdn
                else:
                    SWTDN = np.concatenate((SWTDN, swtdn), axis=2)

            with h5netcdf.File(name2, "r") as f:
                t2m = np.transpose(sf.subset(f["T2M"], param), [1, 2, 0])
                if T2M.size == 0:
                    T2M = t2m
                else:
                    T2M = np.concatenate((T2M, t2m), axis=2)

                u50m = np.transpose(sf.subset(f["U50M"], param), [1, 2, 0])
                if U50M.size == 0:
                    U50M = u50m
                else:
                    U50M = np.concatenate((U50M, u50m), axis=2)

                v50m = np.transpose(sf.subset(f["V50M"], param), [1, 2, 0])
                if V50M.size == 0:
                    V50M = v50m
                else:
                    V50M = np.concatenate((V50M, v50m), axis=2)

            # if date.year != tomorrow.year:
        # Create the overall wind speed
        W50M = abs(U50M + (1j * V50M))
        # Calculate the clearness index
        CLEARNESS = np.divide(SWGDN, SWTDN, out=np.zeros_like(SWGDN), where=SWTDN != 0)

        sys.stdout.write("\n")
        logger.info("Writing Files: T2M, W50M, CLEARNESS")
        hdf5storage.writes({"T2M": T2M}, paths["T2M"], store_python_metadata=True, matlab_compatible=True)
        hdf5storage.writes({"W50M": W50M}, paths["W50M"], store_python_metadata=True, matlab_compatible=True)
        hdf5storage.writes({"CLEARNESS": CLEARNESS}, paths["CLEARNESS"], store_python_metadata=True, matlab_compatible=True)

        if param["MERRA_correction"]:
            cf.clean_weather_data(paths, param)

        ul.create_json(
            paths["T2M"],
            param,
            ["MERRA_coverage", "region_name", "Crd_all", "res_weather", "MERRA_correction", "MERRA_correction_factor"],
            paths,
            ["MERRA_IN", "T2M"],
        )
        ul.create_json(
            paths["W50M"],
            param,
            ["MERRA_coverage", "region_name", "Crd_all", "res_weather", "MERRA_correction", "MERRA_correction_factor"],
            paths,
            ["MERRA_IN", "W50M"],
        )
        ul.create_json(
            paths["CLEARNESS"],
            param,
            ["MERRA_coverage", "region_name", "Crd_all", "res_weather", "MERRA_correction", "MERRA_correction_factor"],
            paths,
            ["MERRA_IN", "CLEARNESS"],
        )

        generate_array_coordinates(paths, param, W50M)

        logger.debug("End")

    
def generate_array_coordinates(paths, param, W50M):
    """
    ToDo: All of this docstring
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

    if os.path.isfile(paths["MERRA_XMIN"]) and os.path.isfile(paths["MERRA_XMAX"]) and \
            os.path.isfile(paths["MERRA_YMIN"]) and os.path.isfile(paths["MERRA_YMAX"]) and \
            os.path.isfile(paths["GWA_X"]) and os.path.isfile(paths["GWA_Y"]):

        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

        ymax, xmax, ymin, xmin = param["Crd_all"]
        res_weather = param["res_weather"]

        # W50M = hdf5storage.read("W50M", paths["W50M"]) # Already given by the function as parameter
        w50m_shape = W50M.shape

        #bounding box coordinates of each pixel in merra
        b_xmin = np.zeros([w50m_shape[0],w50m_shape[1]])
        b_ymin = np.zeros([w50m_shape[0],w50m_shape[1]])
        b_xmax = np.zeros([w50m_shape[0],w50m_shape[1]])
        b_ymax = np.zeros([w50m_shape[0],w50m_shape[1]])
        for row in range(w50m_shape[0]):
            b_ymin[row, :] = ymax - (row + 1) * res_weather[0]
            b_ymax[row, :] = ymax - row * res_weather[0]
        for column in range(w50m_shape[1]):
            b_xmin[:,column] = xmin + column * res_weather[1]
            b_xmax[:,column] = xmin + (column+1) * res_weather[1]

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
            x_gwa[:, l] = GWA_speed.xy(0, l, offset='center')[0]
        for k in range(gwa_rows):
            y_gwa[k, :] = GWA_speed.xy(k, 0, offset='center')[1]

        hdf5storage.writes({"GWA_X": x_gwa}, paths["GWA_X"], store_python_metadata=True, matlab_compatible=True)
        logger.info("files saved: " + paths["GWA_X"])
        hdf5storage.writes({"GWA_Y": y_gwa}, paths["GWA_Y"], store_python_metadata=True, matlab_compatible=True)
        logger.info("files saved: " + paths["GWA_Y"])
        logger.debug("End")


def generate_sea(paths, param):
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

    if os.path.isfile(paths["EEZ"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

        m_high = param["m_high"]
        n_high = param["n_high"]
        Crd_all = param["Crd_all"]
        res_desired = param["res_desired"]
        GeoRef = param["GeoRef"]
        nRegions_sea = param["nRegions_sea"]

        # Extract sea areas
        eez_shp = param["regions_sea"]
        Crd_regions_sea = param["Crd_regions"][-nRegions_sea:]
        Ind = sf.ind_merra(Crd_regions_sea, Crd_all, res_desired)
        A_sea = np.zeros((m_high, n_high))

        for reg in range(0, param["nRegions_sea"]):
            logger.debug('Region: ' + str(reg))
            A_region = sf.calc_region(eez_shp.iloc[reg], Crd_regions_sea[reg, :], res_desired, GeoRef)

            # Include A_region in A_sea
            A_sea[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] = (
                A_sea[(Ind[reg, 2] - 1) : Ind[reg, 0], (Ind[reg, 3] - 1) : Ind[reg, 1]] + A_region
            )

        # Fixing pixels on the borders to avoid duplicates
        A_sea[A_sea > 0] = 1

        # Saving file
        sf.array2raster(paths["EEZ"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_sea)
        logger.info("files saved: " + paths["EEZ"])
        ul.create_json(
            paths["EEZ"], param, ["region_name", "m_high", "n_high", "Crd_all", "res_desired", "GeoRef", "nRegions_sea"], paths, ["EEZ_global", "EEZ"]
        )
        logger.debug("End")


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
    if os.path.isfile(paths["AREA"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

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
        zm_lower = 1 - (e * ul.sin(f_lower))
        zp_lower = 1 + (e * ul.sin(f_lower))

        lowerSliceAreas = np.pi * b ** 2 * ((2 * np.arctanh(e * ul.sin(f_lower))) / (2 * e) + (ul.sin(f_lower) / (zp_lower * zm_lower)))

        # Upper slice areas
        # Areas between the equator and the upper pixel latitudes circling the globe
        f_upper = np.deg2rad(lat_vec + res_desired[0])

        zm_upper = 1 - (e * ul.sin(f_upper))
        zp_upper = 1 + (e * ul.sin(f_upper))

        upperSliceAreas = np.pi * b ** 2 * ((2 * np.arctanh((e * ul.sin(f_upper)))) / (2 * e) + (ul.sin(f_upper) / (zp_upper * zm_upper)))

        # Pixel areas
        # Finding the latitudinal pixel-sized globe slice areas then dividing them by the longitudinal pixel size
        area_vec = ((upperSliceAreas - lowerSliceAreas) * res_desired[1] / 360).T
        A_area = np.tile(area_vec, (1, n_high))

        # Save to HDF File
        hdf5storage.writes({"A_area": A_area}, paths["AREA"], store_python_metadata=True, matlab_compatible=True)
        logger.info("files saved: " + paths["AREA"])
        ul.create_json(paths["AREA"], param, ["Crd_all", "res_desired", "n_high"], paths, [])

        logger.debug("End")


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
    if os.path.isfile(paths["TOPO"]) and os.path.isfile(paths["SLOPE"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

        Crd_all = param["Crd_all"]
        Ind = sf.ind_global(Crd_all, param["res_topography"])[0]
        GeoRef = param["GeoRef"]
        Topo = np.zeros((int(180 / param["res_topography"][0]), int(360 / param["res_topography"][1])))
        tile_extents = np.zeros((24, 4), dtype=int)
        i = 1
        j = 1
        for letter in ul.char_range("A", "X"):
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

        for letter in ul.char_range("A", "X"):
            index = ord(letter) - ord("A")
            if need[index]:
                with rasterio.open(paths["Topo_tiles"] + "15-" + letter + ".tif") as src:
                    tile = src.read()
                Topo[tile_extents[index, 0] - 1 : tile_extents[index, 2], tile_extents[index, 3] - 1 : tile_extents[index, 1]] = tile[0, 0:-1, 0:-1]

        A_TOPO = np.flipud(Topo[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]])
        A_TOPO = sf.adjust_resolution(A_TOPO, param["res_topography"], param["res_desired"], "mean")
        A_TOPO = sf.recalc_topo_resolution(A_TOPO, param["res_topography"], param["res_desired"])
        sf.array2raster(paths["TOPO"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_TOPO)
        logger.info("files saved: " + paths["TOPO"])
        ul.create_json(paths["TOPO"], param, ["region_name", "Crd_all", "res_topography", "res_desired", "GeoRef"], paths, ["Topo_tiles", "TOPO"])

        logger.debug("End")

        generate_slope(paths, param, A_TOPO)


def generate_slope(paths, param, A_TOPO):
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

    if os.path.isfile(paths["SLOPE"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

        # ul.timecheck("Start")
        res_desired = param["res_desired"]
        Crd_all = param["Crd_all"]
        Ind = sf.ind_global(Crd_all, res_desired)[0]
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

        m_per_deg_lat = 111132.954 - 559.822 * ul.cos(np.deg2rad(2 * latMid)) + 1.175 * ul.cos(np.deg2rad(4 * latMid))
        m_per_deg_lon = (np.pi / 180) * 6367449 * ul.cos(np.deg2rad(latMid_2))

        x_cell = ul.repmat(deltaLon[Ind[3] - 1 : Ind[1]], Ind[2]-Ind[0]+1, 1) * ul.repmat(m_per_deg_lon[Ind[0] - 1 : Ind[2]], Ind[1]-Ind[3]+1, 1).T
        x_cell = np.flipud(x_cell)

        y_cell = ul.repmat((deltaLat[Ind[0] - 1 : Ind[2]] * m_per_deg_lat[Ind[0] - 1 : Ind[2]]), Ind[1]-Ind[3]+1, 1).T
        y_cell = np.flipud(y_cell)

        # with rasterio.open(paths["TOPO"]) as src:
        #     A_TOPO = src.read(1)

        kernel = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]) / 8
        dzdx = scipy.ndimage.convolve(A_TOPO, kernel) / x_cell
        kernel = np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]]) / 8
        dzdy = scipy.ndimage.convolve(A_TOPO, kernel) / y_cell

        slope_deg = ul.arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / np.pi
        slope_pc = ul.tan(np.deg2rad(slope_deg)) * 100

        A_SLP = np.flipud(slope_pc)
        sf.array2raster(paths["SLOPE"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_SLP)
        logger.info("files saved: " + paths["SLOPE"])
        ul.create_json(paths["SLOPE"], param, ["region_name", "Crd_all", "res_topography", "res_desired", "GeoRef"], paths, ["TOPO", "SLOPE"])
        logger.debug("End")


def generate_bathymetry(paths, param): #ToDo Adapt to new resolution
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
    if os.path.isfile(paths["BATH"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

        Crd_all = param["Crd_all"]
        Ind = sf.ind_global(Crd_all, param["res_topography"])[0]
        GeoRef = param["GeoRef"]
        with rasterio.open(paths["Bathym_global"]) as src:
            A_BATH = src.read(1)
        #A_BATH = resizem(A_BATH, 180 * 240, 360 * 240)
        A_BATH = sf.adjust_resolution(A_BATH, param["res_bathymetry"], param["res_topography"], "mean")
        A_BATH = np.flipud(A_BATH[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]])
        print (A_BATH.shape)
        A_BATH = sf.recalc_topo_resolution(A_BATH, param["res_topography"], param["res_desired"])

        sf.array2raster(paths["BATH"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_BATH)
        ul.create_json(paths["BATH"], param, ["region_name", "Crd_all", "res_bathymetry", "res_desired", "GeoRef"], paths, ["Bathym_global", "BATH"])
        logger.info("files saved: " + paths["BATH"])

        logger.debug("End")


def generate_landuse(paths, param): #ToDo remove the res_landuse
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

    if os.path.isfile(paths["LU"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")
        Crd_all = param["Crd_all"]
        Ind = sf.ind_global(Crd_all, param["res_desired"])[0]
        GeoRef = param["GeoRef"]
        with rasterio.open(paths["LU_global"]) as src:
            A_lu = src.read(1, window=rasterio.windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]), slice(Ind[3] - 1, Ind[1])))
            A_lu = np.flipud(A_lu).astype(int)
        sf.array2raster(paths["LU"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_lu)
        logger.info("files saved: " + paths["LU"])
        ul.create_json(paths["LU"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["LU_global", "LU"])

    # Buffer map for snow
    if os.path.isfile(paths["SNOW_BUFFER"]):
        logger.info('Skip-SnowBuffer')    # Skip generation if files are already there
    else:
        logger.info("Start-SnowBuffer")
        A_snow = A_lu == 220 # Land use type for snow
        sf.create_buffered_raster(A_snow, param["buffer"]["snow"], param["GeoRef"],
                                  paths["SNOW_BUFFER"])
        logger.info("files saved: " + paths["SNOW_BUFFER"])
        ul.create_json(paths["SNOW_BUFFER"], param, ["region_name", "landuse", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "SNOW_BUFFER"])

    # Buffer map for wetlands
    if os.path.isfile(paths["WETLAND_BUFFER"]):
        logger.info('Skip-WetlandBuffer')  # Skip generation if files are already there
    else:
        logger.info("Start-WetlandBuffer")
        A_wetland = ((A_lu == 160) | (A_lu == 170) | (A_lu == 180))  # Land use type for wetland
        sf.create_buffered_raster(A_wetland, param["buffer"]["wetland"], param["GeoRef"],
                                  paths["WETLAND_BUFFER"])
        logger.info("files saved: " + paths["WETLAND_BUFFER"])
        ul.create_json(paths["WETLAND_BUFFER"], param, ["region_name", "landuse", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "WETLAND_BUFFER"])

    # Buffer map for water
    if os.path.isfile(paths["WATER_BUFFER"]):
        logger.info('Skip-WaterBuffer')    # Skip generation if files are already there
    else:
        logger.info("Start-WaterBuffer")
        A_water = A_lu == 210 # Land use type for water
        sf.create_buffered_raster(A_water, param["buffer"]["water"], param["GeoRef"],
                                  paths["WATER_BUFFER"])
        logger.info("files saved: " + paths["WATER_BUFFER"])
        ul.create_json(paths["WATER_BUFFER"], param, ["region_name", "landuse", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "WATER_BUFFER"])

    logger.debug("End")


def generate_protected_areas(paths, param):
    """
    This function reads the shapefile of the globally protected areas, adds an attribute whose values are based on the dictionary
    of conversion (protected_areas) to identify the protection category, then converts the shapefile into a raster for the scope.
    The values are integers from 0 to 10.

    :param paths: Dictionary including the paths to the shapefile of the globally protected areas, to the land raster of the scope, and to the output path PA.
    :type paths: dict
    :param param: Dictionary including the dictionary of conversion of protection categories (protected_areas).
    :type param: dict
    :return: The tif file for PA is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """

    if os.path.isfile(paths["PA"]):
        logger.info('Skip')    # Skip generation if files are already there
    else:
        logger.info("Start")
        protected_areas = param["protected_areas"]
        # set up protected areas dictionary
        protection_type = dict(zip(protected_areas["IUCN_Category"], protected_areas["type"]))
        sf.shape2raster(paths["Protected"],paths["PA"],"IUCN_CAT",protection_type,paths["LAND"])
        logger.info("files saved: " + paths["PA"])
        ul.create_json(paths["PA"], param, ["region_name", "protected_areas", "Crd_all", "res_desired", "GeoRef"], paths, ["Protected", "PA"])

    # Buffer maps for protected areas
    logger.info("Start-Buffer")
    with rasterio.open(paths["PA"]) as src:
        A_pa = src.read(1)
    A_pa = np.flipud(A_pa).astype(int)
    A_pa = A_pa > 0 #All protected areas pixels

    if os.path.isfile(paths["PV_PA_BUFFER"]):
        logger.info('Skip-PV')  # Skip generation if files are already there
    else:
        logger.info("Start-PV")
        sf.create_buffered_raster(A_pa, param["buffer"]["protected_areas_pv"], param["GeoRef"],
                                  paths["PV_PA_BUFFER"])
        logger.info("files saved: " + paths["PV_PA_BUFFER"])
        ul.create_json(paths["PV_PA_BUFFER"], param,
                       ["region_name", "protected_areas", "Crd_all", "res_desired", "GeoRef"], paths,
                       ["PA", "PV_PA_BUFFER"])

    # if "WindOn" in param["technology"]:
    if os.path.isfile(paths["WINDON_PA_BUFFER"]):
        logger.info('Skip-WindOn')  # Skip generation if files are already there
    else:
        logger.info("Start-WindOn")
        sf.create_buffered_raster(A_pa, param["buffer"]["protected_areas_windon"], param["GeoRef"],
                                  paths["WINDON_PA_BUFFER"])
        logger.info("files saved: " + paths["WINDON_PA_BUFFER"])
        ul.create_json(paths["WINDON_PA_BUFFER"], param,
                       ["region_name", "protected_areas", "WindOn", "Crd_all", "res_desired", "GeoRef"], paths,
                       ["PA", "WINDON_PA_BUFFER"])

    logger.debug("End")


def generate_airports(paths,param):

    if os.path.isfile(paths["AIRPORTS"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")
        Crd_all = param["Crd_all"]
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
            airports_filtered = airports_list[airports_list.index==alpha2code]
            airports_filtered = airports_filtered[(airports_filtered["type"] == "small_airport")
                                                  | (airports_filtered["type"] == "medium_airport")
                                                  | (airports_filtered["type"] == "large_airport")]
            airports.append(airports_filtered)
        airports = pd.concat(airports)
        logger.info("Airports are filtered")

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
            ind = sf.crd2ind(crd, Crd_all, res_desired)

            A_land[tuple(ind)]=100
            airport_raster = A_land == 100
            sf.create_buffered_raster(airport_raster, param["buffer"]["airport_windon"], param["GeoRef"],
                                      paths["AIRPORTS"])
            logger.info("files saved: " + paths["AIRPORTS"])
            ul.create_json(paths["AIRPORTS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["LAND", "AIRPORTS"])

        logger.debug("End")


def generate_country_boarders(paths,param):

    if os.path.isfile(paths["BOARDERS"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")

        Crd_all = param["Crd_all"]
        GeoRef = param["GeoRef"]
        res_desired = param["res_desired"]
        countries_shp = param["regions_land"]
        nCountries = param["nRegions_land"]
        m_high = param["m_high"]
        n_high = param["n_high"]

        A_countries_buffered = np.zeros([m_high, n_high]).astype(int)

        buffer_pixel_amount = param["buffer"]["boarder"]
        kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
        kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)

        for reg in range(0, nCountries):
            try:
                A_country_area = sf.calc_region(countries_shp.loc[reg], Crd_all, res_desired, GeoRef)
                A_country_buffered = scipy.ndimage.minimum_filter(A_country_area, footprint=kernel, mode="constant", cval=1)
                A_countries_buffered = A_countries_buffered + A_country_buffered
            except:
                pass
        A_countries_buffered = A_countries_buffered > 0
        A_notBoarder = (A_countries_buffered).astype(int)

        sf.array2raster(paths["BOARDERS"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_notBoarder)
        logger.info("files saved: " + paths["BOARDERS"])
        ul.create_json(paths["BOARDERS"], param, ["region_name", "landuse", "Crd_all", "res_desired", "GeoRef"], paths, ["LAND", "BOARDERS"])

        logger.debug("End")


def generate_roads(paths, param):
    # "fclass" ILIKE 'primary' OR "fclass" ILIKE 'secondary'
    if os.path.isfile(paths["ROADS"]):
        logger.info('Skip')  # Skip generation if files are already there

    else:
        logger.info("Start")

        sf.shape2raster(paths["OSM_Roads"], paths["ROADS"], "", [], paths["LAND"])
        logger.info("files saved: " + paths["ROADS"])
        ul.create_json(paths["ROADS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["ROADS"])

        logger.debug("End")


def generate_railways(paths, param):
    # "fclass" ILIKE 'primary' OR "fclass" ILIKE 'secondary'
    if os.path.isfile(paths["RAILS"]):
        logger.info('Skip')  # Skip generation if files are already there

    else:
        logger.info("Start")

        sf.shape2raster(paths["OSM_Railways"], paths["RAILS"], "", [], paths["LAND"])
        logger.info("files saved: " + paths["RAILS"])
        ul.create_json(paths["RAILS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["RAILS"])

        logger.debug("End")


def generate_osm_areas(paths, param):
    # "fclass" ILIKE 'commercial' OR "fclass" ILIKE 'industrial' OR "fclass" ILIKE 'quarry' OR "fclass" ILIKE 'military' OR "fclass" ILIKE 'park' OR "fclass" ILIKE 'recreation_ground'

    """
    This function reads the osm land use shapefile, identifies several areas, and excludes pixels around them based on a
    user-defined buffers *buffer_pixel_amount*. It creates a masking raster of boolean values (0 or 1) for the scope.
    Zero means the pixel is excluded, one means it is suitable.
    The function is useful in case there is a policy to exclude renewable energy projects next to certain type of areas.

    :param paths: Dictionary including the path to the osm land-use shapefile, and to the output path BUFFER.
    :type paths: dict
    :param param: Dictionary including the user-defined buffers (buffer_pixel_amount) and the georeference dictionary.
    :type param: dict

    :return: The tif file for BUFFER is saved in its respective path, along with its metadata in a JSON file.
    :rtype: None
    """

    if os.path.isfile(paths["OSM_AREAS"]):
        logger.info('Skip')  # Skip generation if files are already there
    else:
        logger.info("Start")

        osm_areas = param["osm_areas"]
        # set up osm areas dictionary
        osm_type = dict(zip(osm_areas["Category"], osm_areas["type"]))
        sf.shape2raster(paths["OSM_Landuse"], paths["OSM_AREAS"], "fclass", osm_type, paths["LAND"])
        logger.info("files saved: " + paths["OSM_AREAS"])
        ul.create_json(paths["OSM_AREAS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                    ["OSM_AREAS"])

    # Create Buffer osm maps
    logger.info("Start-Buffer")
    with rasterio.open(paths["OSM_AREAS"]) as src:
        A_osma = src.read(1)
    A_osma = np.flipud(A_osma).astype(int)

    # Commercial Areas
    if os.path.isfile(paths["OSM_COM_BUFFER"]):
        logger.info('Skip-Commercial')  # Skip generation if files are already there
    else:
        logger.info('Start-Commercial')
        A_com = A_osma == 1
        sf.create_buffered_raster(A_com, param["buffer"]["commercial_windon"], param["GeoRef"],
                                  paths["OSM_COM_BUFFER"])
        logger.info("files saved: " + paths["OSM_COM_BUFFER"])
        ul.create_json(paths["OSM_COM_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                       ["OSM_COM_BUFFER"])

    # Industrical Areas
    if os.path.isfile(paths["OSM_IND_BUFFER"]):
        logger.info('Skip-Industrial')  # Skip generation if files are already there
    else:
        logger.info('Start-Industrial')
        A_ind = A_osma == 2
        sf.create_buffered_raster(A_ind, param["buffer"]["industrial_windon"], param["GeoRef"],
                                  paths["OSM_IND_BUFFER"])
        logger.info("files saved: " + paths["OSM_IND_BUFFER"])
        ul.create_json(paths["OSM_IND_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                       ["OSM_IND_BUFFER"])

    # Mining Areas
    if os.path.isfile(paths["OSM_MINE_BUFFER"]):
        logger.info('Skip-Mining')  # Skip generation if files are already there
    else:
        logger.info('Start-Mining')
        A_min = A_osma == 3
        sf.create_buffered_raster(A_min, param["buffer"]["mining"], param["GeoRef"],
                                  paths["OSM_MINE_BUFFER"])
        logger.info("files saved: " + paths["OSM_MINE_BUFFER"])
        ul.create_json(paths["OSM_MINE_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                       ["OSM_MINE_BUFFER"])

    # Military Areas
    if os.path.isfile(paths["OSM_MIL_BUFFER"]):
        logger.info('Skip-Military')  # Skip generation if files are already there
    else:
        logger.info('Start-Military')
        A_mil = A_osma == 4
        sf.create_buffered_raster(A_mil, param["buffer"]["military_windon"], param["GeoRef"],
                                  paths["OSM_MIL_BUFFER"])
        logger.info("files saved: " + paths["OSM_MIL_BUFFER"])
        ul.create_json(paths["OSM_MIL_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                       ["OSM_MIL_BUFFER"])

    # Parks
    logger.info('Start-Park')
    A_park = A_osma == 5
    # Wind onshore
    if os.path.isfile(paths["WINDON_OSM_PARK_BUFFER"]):
        logger.info('Skip-WindOn')  # Skip generation if files are already there
    else:
        logger.info('Start-WindOn')
        sf.create_buffered_raster(A_park, param["buffer"]["park_windon"], param["GeoRef"],
                                  paths["WINDON_OSM_PARK_BUFFER"])
        logger.info("files saved: " + paths["WINDON_OSM_PARK_BUFFER"])
        ul.create_json(paths["WINDON_OSM_PARK_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"],
                       paths,
                       ["WINDON_OSM_PARK_BUFFER"])
    # PV
    if os.path.isfile(paths["PV_OSM_PARK_BUFFER"]):
        logger.info('Skip-PV')  # Skip generation if files are already there
    else:
        logger.info('Start-PV')
        sf.create_buffered_raster(A_park, param["buffer"]["park_pv"], param["GeoRef"],
                                  paths["PV_OSM_PARK_BUFFER"])
        logger.info("files saved: " + paths["PV_OSM_PARK_BUFFER"])
        ul.create_json(paths["PV_OSM_PARK_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"],
                       paths,
                       ["PV_OSM_PARK_BUFFER"])

    # Recreational Areas
    if os.path.isfile(paths["OSM_REC_BUFFER"]):
        logger.info('Skip-Recreation')  # Skip generation if files are already there
    else:
        logger.info('Start-Recreation')
        A_rec = A_osma == 6
        sf.create_buffered_raster(A_rec, param["buffer"]["recreation_windon"], param["GeoRef"],
                                  paths["OSM_REC_BUFFER"])
        logger.info("files saved: " + paths["OSM_REC_BUFFER"])
        ul.create_json(paths["OSM_REC_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                       ["OSM_REC_BUFFER"])

    logger.debug("End")


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
    Crd_all = param["Crd_all"]
    Ind = sf.ind_global(Crd_all, param["res_desired"])[0]
    GeoRef = param["GeoRef"]
    logger.info("Start")

    if os.path.isfile(paths["WSF"]):
        logger.info('Skip')  # Skip generation if files are already there
    else:
        with rasterio.open(paths["WSF_global"]) as src:
            A_wsf = src.read(1, window=rasterio.windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]),
                                                                  slice(Ind[3] - 1, Ind[1]))).astype(bool)
            A_wsf = np.flipud(A_wsf)
        sf.array2raster(paths["WSF"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_wsf)
        logger.info("files saved: " + paths["WSF"])
        ul.create_json(paths["WSF"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                    ["WSF_global", "WSF"])

    # PV buffer
    logger.info('Start-Buffer')
    if os.path.isfile(paths["PV_WSF_BUFFER"]):
        logger.info('Skip-PV')  # Skip generation if files are already there
    else:
        logger.info('Start-PV')
        sf.create_buffered_raster(A_wsf, param["buffer"]["settlement_pv"], param["GeoRef"],
                                  paths["PV_WSF_BUFFER"])
        logger.info("files saved: " + paths["PV_WSF_BUFFER"])
        ul.create_json(paths["PV_WSF_BUFFER"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"],
                    paths, ["WSF_global", "PV_WSF_BUFFER"])

    # Onshore wind buffer
    if os.path.isfile(paths["WINDON_WSF_BUFFER"]):
        logger.info('Skip-WindOn')  # Skip generation if files are already there
    else:
        logger.info('Start-WindOn')
        sf.create_buffered_raster(A_wsf, param["buffer"]["settlement_windon"], param["GeoRef"],
                                  paths["WINDON_WSF_BUFFER"])
        logger.info("files saved: " + paths["WINDON_WSF_BUFFER"])
        ul.create_json(paths["WINDON_WSF_BUFFER"], param,
                    ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                    ["WSF_global", "WINDON_WSF_BUFFER"])

    logger.debug("End")


def generate_HydroLakes(paths, param):
    if os.path.isfile(paths["HYDROLAKES"]):
        logger.info('Skip')  # Skip generation if files are already there
    else:
        logger.info('Start')
        sf.shape2raster(paths["HydroLakes"], paths["HYDROLAKES"], "", [], paths["LAND"])
        logger.info("files saved: " + paths["HYDROLAKES"])
        ul.create_json(paths["HYDROLAKES"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["HYDROLAKES"])

    # Create Buffer
    if os.path.isfile(paths["HYDROLAKES_BUFFER"]):
        logger.info('Skip-Buffer')  # Skip generation if files are already there
    else:
        logger.info('Start-Buffer')
        with rasterio.open(paths["HYDROLAKES"]) as src:
            A_lake = src.read(1).astype(bool)
            A_lake = np.flipud(A_lake)
        sf.create_buffered_raster(A_lake, param["buffer"]["hydrolakes"], param["GeoRef"],
                                  paths["HYDROLAKES_BUFFER"])
        logger.info("files saved: " + paths["HYDROLAKES_BUFFER"])
        ul.create_json(paths["HYDROLAKES_BUFFER"], param,
                    ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                    ["HYDROLAKES", "HYDROLAKES_BUFFER"])

    logger.debug("End")


def generate_HydroRivers(paths, param):
    if os.path.isfile(paths["HYDRORIVERS"]):
        logger.info('Skip')  # Skip generation if files are already there
    else:
        logger.info('Start')
        sf.shape2raster(paths["HydroRivers"], paths["HYDRORIVERS"], "", [], paths["LAND"])
        logger.info("files saved: " + paths["HYDRORIVERS"])
        ul.create_json(paths["HYDRORIVERS"], param, ["region_name", "Crd_all", "res_desired", "GeoRef"], paths, ["HYDROLAKES"])

    # Create Buffer
    if os.path.isfile(paths["HYDRORIVERS_BUFFER"]):
        logger.info('Skip-Buffer')  # Skip generation if files are already there
    else:
        logger.info('Start-Buffer')
        with rasterio.open(paths["HYDRORIVERS"]) as src:
            A_Riv = src.read(1).astype(bool)
            A_Riv = np.flipud(A_Riv)
        sf.create_buffered_raster(A_Riv,param["buffer"]["hydrorivers_pv"],param["GeoRef"],paths["HYDRORIVERS_BUFFER"])
        logger.info("files saved: " + paths["HYDRORIVERS_BUFFER"])
        ul.create_json(paths["HYDRORIVERS_BUFFER"], param,
                    ["region_name", "Crd_all", "res_desired", "GeoRef"], paths,
                    ["HYDRORIVERS", "HYDRORIVERS_BUFFER"])

    logger.debug("End")


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
    logger.info("Start")
    res_desired = param["res_desired"]
    Crd_all = param["Crd_all"]
    Ind = sf.ind_global(Crd_all, param["res_livestock"])[0]
    GeoRef = param["GeoRef"]

    A_area = hdf5storage.read("A_area", paths["AREA"])

    for animal in param["Biomass"]["livestock"]["animal"]:
        with rasterio.open(paths["LS_global"] + animal + "_2006.tif") as src:
            A_LS = src.read(1, window=rasterio.windows.Window.from_slices(slice(Ind[0] - 1, Ind[2]),
                                                                          slice(Ind[3] - 1, Ind[1])))
        A_LS = np.flipud(A_LS)
        A_LS = sf.recalc_livestock_resolution(A_LS, param["res_livestock"], param["res_desired"])
        # print (np.size(A_LS))
        A_LS[A_LS < 0] = float(0)
        A_LS = np.multiply(A_LS, A_area) / (10 ** 6)
        sf.array2raster(paths["LS"] + animal + ".tif", GeoRef["RasterOrigin"], GeoRef["pixelWidth"],
                        GeoRef["pixelHeight"], A_LS)
        logger.info("files saved: " + paths["LS"] + animal + ".tif")
        ul.create_json(paths["LS"] + animal + ".tif", param,
                       ["region_name", "Crd_all", "res_livestock", "res_desired", "GeoRef"], paths, ["LS_global", "LS"])

    logger.debug("End")


# def generate_population(paths, param): #ToDo function Not needed
#     """
#     This function reads the global map of population density, resizes it, and creates a raster out of it for the desired scope.
#     The values are in population per pixel.
#
#     :param paths: Dictionary including the paths to the global population raster *Pop_global* and to the output path *POP*.
#     :type paths: dict
#     :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and the georeference dictionary.
#     :type param: dict
#
#     :return: The tif file for *POP* is saved in its respective path, along with its metadata in a JSON file.
#     :rtype: None
#     """
#     ul.timecheck("Start")
#     res_desired = param["res_desired"]
#     Crd_all = param["Crd_all"]
#     Ind = sf.ind_global(Crd_all, param["res_desired"])[0]
#     GeoRef = param["GeoRef"]
#     with rasterio.open(paths["Pop_global"]) as src:
#         A_POP_part = src.read(1)  # map is only between latitudes -60 and 85
#     A_POP = np.zeros((21600, 43200))
#     A_POP[600:18000, :] = A_POP_part
#     #A_POP = adjust_resolution(A_POP, param["res_population"], param["res_desired"], "sum")
#     #A_POP = resizem(A_POP, 180 * 240, 360 * 240) / 4  # density is divided by 4
#     A_POP = np.flipud(A_POP[Ind[0] - 1 : Ind[2], Ind[3] - 1 : Ind[1]])
#     #if "WindOn" in param["technology"]:
#     A_POP = sf.recalc_topo_resolution(A_POP, param["res_landuse"], param["res_desired"])
#     sf.array2raster(paths["POP"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_POP)
#     logger.info("files saved: " + paths["POP"])
#     ul.create_json(paths["POP"], param, ["region_name", "Crd_all", "res_population", "res_desired", "GeoRef"], paths, ["Pop_global", "POP"])
#     ul.timecheck("End")