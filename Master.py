import os
from data_functions import *
from model_functions import *
import numpy as np
import datetime
import geopandas as gpd
import pandas as pd
from rasterio import windows
from shapely.geometry import mapping, Point
import fiona
import hdf5storage
from multiprocessing import Pool
from itertools import product


def initialization():
    """
    Initialization of the script, reads config file and returns the paths and param dictionaries
    Paths: contains all required paths for loading and saving
    param: contains many parameters used throughout the script
    """

    # import param and paths
    from config import paths, param
    res = param["res"]

    # read shapefile of regions
    regions_shp = gpd.read_file(paths["SHP"])
    # Extract onshore and offshore areas separately
    param["regions_land"] = regions_shp.drop(regions_shp[regions_shp["Population"] == 0].index)
    param["regions_eez"] = regions_shp.drop(regions_shp[regions_shp["Population"] != 0].index)
    # Recombine the maps in this order: onshore then offshore
    regions_all = gpd.GeoDataFrame(pd.concat([param["regions_land"], param["regions_eez"]],
                                             ignore_index=True), crs=param["regions_land"].crs)

    param["nRegions_land"] = len(param["regions_land"])
    param["nRegions_eez"] = len(param["regions_eez"])

    nRegions = param["nRegions_land"] + param["nRegions_eez"]
    Crd = np.zeros((nRegions + 1, 4))
    for reg in range(0, nRegions):
        # Box coordinates for MERRA2 data
        r = regions_all.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd[reg, :] = crd_merra_low(box, res)
    Crd[-1, :] = [max(Crd[:-1, 0]), max(Crd[:-1, 1]), min(Crd[:-1, 2]), min(Crd[:-1, 3])]
    param["Crd"] = Crd

    # Indices and matrix dimensions
    Ind = np.zeros((3, nRegions + 1, 4))
    Ind[0] = ind_merra(Crd, res[0, :])  # Range indices for MERRA2 data (centroids)
    Ind[1] = ind_merra(Crd, res[1, :])  # Range indices for high resolution matrices, superposed to MERRA2 data
    Ind[2] = ind_global(Crd, res[1, :])
    param["Ind"] = Ind.astype(int)

    m = Ind[:, :, 0] - Ind[:, :, 2] + 1  # #rows
    m[2, :] = Ind[2, :, 2] - Ind[2, :, 0] + 1  # starts counting from North
    n = Ind[:, :, 1] - Ind[:, :, 3] + 1  # #Cols
    param["m"] = m.astype(int)
    param["n"] = n.astype(int)
    param["GeoRef"] = calc_geotiff(Crd, res)
    return paths, param


def generate_weather_files(paths):
    """
    This function reads the daily NetCDF data (from MERRA) for SWGDN, SWTDN, T2M, U50m, and V50m,
    and saves them in matrices with yearly time series with low spatial resolution.
    This code has to be run only once

    :param paths: paths dictionary containing the input file for NetCDF data
    """

    if not os.path.isfile(paths["U50M"]):
        start = datetime.date(paths["year"], 1, 1)
        end = datetime.date(paths["year"], 12, 31)
        root = paths["MERRA_IN"]

        SWGDN = np.array([])
        SWTDN = np.array([])
        T2M = np.array([])
        U50M = np.array([])
        V50M = np.array([])

        for date in pd.date_range(start, end):
            tomorrow = date + pd.Timedelta('1 day')
            if date.day == 29 and date.month == 2:
                continue

            # Name and path of the NetCDF file to be read
            name = root + 'MERRA2_400.tavg1_2d_rad_Nx.' + date.strftime('%Y%m%d') + '.SUB.nc'
            name2 = root + 'MERRA2_400.tavg1_2d_slv_Nx.' + date.strftime('%Y%m%d') + '.SUB.nc'

            # Read NetCDF file, extract hourly tables
            with h5netcdf.File(name, 'r') as f:
                swgdn = np.transpose(f['SWGDN'], [1, 0, 2])
                if SWGDN.size == 0:
                    SWGDN = swgdn
                else:
                    SWGDN = np.concatenate((SWGDN, swgdn), axis=2)

                swtdn = np.transpose(f['SWTDN'], [1, 0, 2])
                if SWTDN.size == 0:
                    SWTDN = swtdn
                else:
                    SWTDN = np.concatenate((SWTDN, swtdn), axis=2)

            with h5netcdf.File(name2, 'r') as f:
                t2m = np.transpose(f['T2M'], [1, 0, 2])
                if T2M.size == 0:
                    T2M = t2m
                else:
                    T2M = np.concatenate((T2M, t2m), axis=2)

                u50m = np.transpose(f['U50M'], [1, 0, 2])
                if U50M.size == 0:
                    U50M = u50m
                else:
                    U50M = np.concatenate((U50M, u50m), axis=2)

                v50m = np.transpose(f['V50M'], [1, 0, 2])
                if V50M.size == 0:
                    V50M = v50m
                else:
                    V50M = np.concatenate((V50M, v50m), axis=2)
            if date.year != tomorrow.year:
                hdf5storage.writes({'SWGDN': SWGDN}, paths["SWGDN"], store_python_metadata=True, matlab_compatible=True)
                hdf5storage.writes({'SWTDN': SWTDN}, paths["SWTDN"], store_python_metadata=True, matlab_compatible=True)
                hdf5storage.writes({'T2M': T2M}, paths["T2M"], store_python_metadata=True, matlab_compatible=True)
                hdf5storage.writes({'U50M': U50M}, paths["U50M"], store_python_metadata=True, matlab_compatible=True)
                hdf5storage.writes({'V50M': V50M}, paths["V50M"], store_python_metadata=True, matlab_compatible=True)
                # Create the overall wind speed
                W50M = abs(U50M + (1j * V50M))
                hdf5storage.writes({'W50M': W50M}, paths["W50M"], store_python_metadata=True, matlab_compatible=True)
                # Calculate the clearness index
                CLEARNESS = np.divide(SWGDN, SWTDN, where=SWTDN != 0)
                hdf5storage.writes({'CLEARNESS': CLEARNESS}, paths["CLEARNESS"], store_python_metadata=True,
                                   matlab_compatible=True)


def generate_landsea(paths, param):
    m = param["m"]
    n = param["n"]
    res = param["res"]
    GeoRef = param["GeoRef"]
    if not os.path.isfile(paths["LAND"]):
        nRegions = param["nRegions_land"]
        regions_shp = param["regions_land"]
        Crd = param["Crd"][0:nRegions, :]
        Ind = param["Ind"][:, 0:nRegions, :]
        A_land = np.zeros((m[1, -1], n[1, -1]))

        for reg in range(0, nRegions):
            A_region = calc_region(regions_shp.iloc[reg], Crd[reg, :], res, GeoRef)
            A_land[(Ind[1, reg, 2] - 1):Ind[1, reg, 0], (Ind[1, reg, 3] - 1):Ind[1, reg, 1]] = \
                A_land[(Ind[1, reg, 2] - 1):Ind[1, reg, 0], (Ind[1, reg, 3] - 1):Ind[1, reg, 1]] + A_region
        array2raster(paths["LAND"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_land)
        print("files saved: " + paths["LAND"])

    if not os.path.isfile(paths["EEZ"]):
        nRegions = param["nRegions_eez"]
        regions_shp = param["regions_eez"]
        Crd = param["Crd"][-nRegions-1:-1, :]
        Ind = param["Ind"][:, -nRegions-1:-1, :]
        A_eez = np.zeros((m[1, -1], n[1, -1]))

        for reg in range(0, nRegions):
            A_region = calc_region(regions_shp.iloc[reg], Crd[reg, :], res, GeoRef)
            A_eez[(Ind[1, reg, 2] - 1):Ind[1, reg, 0], (Ind[1, reg, 3] - 1):Ind[1, reg, 1]] = \
                A_eez[(Ind[1, reg, 2] - 1):Ind[1, reg, 0], (Ind[1, reg, 3] - 1):Ind[1, reg, 1]] + A_region
        with rasterio.open(paths["LAND"]) as src:
            A_land = np.flipud(src.read(1)).astype(int)
        A_eez = A_eez * (1 - A_land)
        array2raster(paths["EEZ"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_eez)
        print("files saved: " + paths["EEZ"])


def generate_landuse(paths, param):
    if not os.path.isfile(paths['LU']):
        Ind = param["Ind"]
        GeoRef = param["GeoRef"]
        with rasterio.open(paths["LU_global"]) as src:
            w = src.read(1, window=windows.Window.from_slices(slice(Ind[2, -1, 0] - 1, Ind[2, -1, 2]),
                                                              slice(Ind[2, -1, 3] - 1, Ind[2, -1, 1])))
        w = np.flipud(w)
        array2raster(paths["LU"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], w)
        print("files saved: " + paths["LU"])


def generate_bathymetry(paths, param):
    if not os.path.isfile(paths['BATH']):
        Ind = param["Ind"]
        GeoRef = param["GeoRef"]
        with rasterio.open(paths["Bathym_global"]) as src:
            A_BATH = src.read(1)
        A_BATH = resizem(A_BATH, 180 * 240, 360 * 240)
        A_BATH = np.flipud(A_BATH[Ind[2, -1, 0] - 1: Ind[2, -1, 2], Ind[2, -1, 3] - 1: Ind[2, -1, 1]])
        array2raster(paths['BATH'], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_BATH)
        print("files saved: " + paths["BATH"])


def generate_topography(paths, param):
    if not os.path.isfile(paths["TOPO"]):
        Ind = param["Ind"]
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
        n_min = (Ind[2, -1, 0] // (45 * 240)) * 45 * 240 + 1
        e_max = (Ind[2, -1, 1] // (60 * 240) + 1) * 60 * 240
        s_max = (Ind[2, -1, 2] // (45 * 240) + 1) * 45 * 240
        w_min = (Ind[2, -1, 3] // (60 * 240)) * 60 * 240 + 1

        need = np.logical_and((np.logical_and((tile_extents[:, 0] >= n_min), (tile_extents[:, 1] <= e_max))),
                              np.logical_and((tile_extents[:, 2] <= s_max), (tile_extents[:, 3] >= w_min)))

        for letter in char_range('A', 'X'):
            index = ord(letter) - ord('A')
            if need[index]:
                with rasterio.open(paths["Topo_tiles"] + '15-' + letter + '.tif') as src:
                    tile = src.read()
                Topo[tile_extents[index, 0] - 1: tile_extents[index, 2],
                tile_extents[index, 3] - 1: tile_extents[index, 1]] = \
                    tile[0, 0:-1, 0:-1]

        A_TOPO = np.flipud(Topo[Ind[2, -1, 0] - 1:Ind[2, -1, 2], Ind[2, -1, 3] - 1:Ind[2, -1, 1]])
        array2raster(paths["TOPO"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_TOPO)
        print("files saved: " + paths["TOPO"])


def generate_slope(paths, param):
    if not os.path.isfile(paths["SLOPE"]):
        Ind = param["Ind"]
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
        x_cell = x_cell[Ind[2, -1, 0] - 2:Ind[2, -1, 2] + 1, Ind[2, -1, 3] - 2: Ind[2, -1, 1] + 1]
        x_cell = np.flipud(x_cell)

        y_cell = repmat((deltaLat * m_per_deg_lat), 360 * 240, 1).T
        y_cell = y_cell[Ind[2, -1, 0] - 2:Ind[2, -1, 2] + 1, Ind[2, -1, 3] - 2: Ind[2, -1, 1] + 1]
        y_cell = np.flipud(y_cell)

        with rasterio.open(paths["TOPO"]) as src:
            A_TOPO = src.read(1)

        topo_a = np.zeros((Ind[2, -1, 2] - Ind[2, -1, 0] + 3, Ind[2, -1, 1] - Ind[2, -1, 3] + 3))
        topo_a[0:-2, 0:-2] = A_TOPO

        topo_b = np.zeros((Ind[2, -1, 2] - Ind[2, -1, 0] + 3, Ind[2, -1, 1] - Ind[2, -1, 3] + 3))
        topo_b[0:-2, 1:-1] = A_TOPO

        topo_c = np.zeros((Ind[2, -1, 2] - Ind[2, -1, 0] + 3, Ind[2, -1, 1] - Ind[2, -1, 3] + 3))
        topo_c[0:-2, 2:] = A_TOPO

        topo_d = np.zeros((Ind[2, -1, 2] - Ind[2, -1, 0] + 3, Ind[2, -1, 1] - Ind[2, -1, 3] + 3))
        topo_d[1:-1, 0:-2] = A_TOPO

        topo_f = np.zeros((Ind[2, -1, 2] - Ind[2, -1, 0] + 3, Ind[2, -1, 1] - Ind[2, -1, 3] + 3))
        topo_f[1:-1, 2:] = A_TOPO

        topo_g = np.zeros((Ind[2, -1, 2] - Ind[2, -1, 0] + 3, Ind[2, -1, 1] - Ind[2, -1, 3] + 3))
        topo_g[2:, 0:-2] = A_TOPO

        topo_h = np.zeros((Ind[2, -1, 2] - Ind[2, -1, 0] + 3, Ind[2, -1, 1] - Ind[2, -1, 3] + 3))
        topo_h[2:, 1:-1] = A_TOPO

        topo_i = np.zeros((Ind[2, -1, 2] - Ind[2, -1, 0] + 3, Ind[2, -1, 1] - Ind[2, -1, 3] + 3))
        topo_i[2:, 2:] = A_TOPO

        dzdx = ((topo_c + 2 * topo_f + topo_i) - (topo_a + 2 * topo_d + topo_g)) / (8 * x_cell)
        dzdy = ((topo_g + 2 * topo_h + topo_i) - (topo_a + 2 * topo_b + topo_c)) / (8 * y_cell)

        slope_deg = arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / np.pi
        slope_pc = tan(np.deg2rad(slope_deg)) * 100

        A_SLP = np.flipud(slope_pc[1:-1, 1:-1])
        array2raster(paths["SLOPE"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_SLP)
        print("files saved: " + paths["SLOPE"])


def generate_population(paths, param):
    if not os.path.isfile(paths["POP"]):
        Ind = param["Ind"]
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
        n_min = (Ind[2, -1, 0] // (45 * 240)) * 45 * 240 + 1
        e_max = (Ind[2, -1, 1] // (60 * 240) + 1) * 60 * 240
        s_max = (Ind[2, -1, 2] // (45 * 240) + 1) * 45 * 240
        w_min = (Ind[2, -1, 3] // (60 * 240)) * 60 * 240 + 1

        need = np.logical_and((np.logical_and((tile_extents[:, 0] >= n_min), (tile_extents[:, 1] <= e_max))),
                              np.logical_and((tile_extents[:, 2] <= s_max), (tile_extents[:, 3] >= w_min)))

        for letter in char_range('A', 'X'):
            index = ord(letter) - ord('A')
            if need[index]:
                with rasterio.open(paths["Pop_tiles"] + letter + '.tif') as src:
                    tile = src.read()
                Pop[tile_extents[index, 0] - 1: tile_extents[index, 2],
                tile_extents[index, 3] - 1: tile_extents[index, 1]] = \
                    tile[0]
        A_POP = np.flipud(Pop[Ind[2, -1, 0] - 1:Ind[2, -1, 2], Ind[2, -1, 3] - 1:Ind[2, -1, 1]])
        array2raster(paths["POP"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_POP)
        print("files saved: " + paths["POP"])


def generate_protected_areas(paths, param):
    if not os.path.isfile(paths["PA"]):
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

        # Close dataset
        out_raster_ds = None
        print("files saved: " + paths["PA"])


def generate_buffered_population(paths, param):
    buffer_pixel_amount = param["WindOn"]["mask"]["buffer_pixel_amount"]
    if buffer_pixel_amount:
        GeoRef = param["GeoRef"]
        with rasterio.open(paths["LU"]) as src:
            A_lu = src.read(1)
            A_lu = np.flipud(A_lu).astype(int)

        A_lu_buffered = create_buffer(A_lu, buffer_pixel_amount)
        A_notPopulated = (~A_lu_buffered).astype(int)

        array2raster(paths["BUFFER"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"],
                     A_notPopulated)
        print("files saved: " + paths["BUFFER"])


def generate_wind_correction(paths, param):
    res_correction = param["WindOn"]["resource"]["res_correction"]
    topo_correction = param["WindOn"]["resource"]["topo_correction"]
    GeoRef = param["GeoRef"]
    if res_correction:
        landuse = param["landuse"]
        turbine_height_on = param["WindOn"]["technical"]["hub_height"]
        turbine_height_off = param["WindOff"]["technical"]["hub_height"]
        m = param["m"]
        n = param["n"]
        res = param["res"]
        with rasterio.open(paths["LU"]) as src:
            A_lu = np.flipud(src.read(1)).astype(int)
        A_hellmann = changem(A_lu, landuse["hellmann"], landuse["type"]).astype(float) / 100
        A_gradient_height = changem(A_lu.astype(float), landuse["height"], landuse["type"])
        del A_lu
        Sigma = sumnorm_MERRA2((50 / A_gradient_height) ** A_hellmann, m[0, -1], n[0, -1], res)
        A_cf_on = ((turbine_height_on / 50) * turbine_height_on /
                   A_gradient_height) ** A_hellmann / resizem(Sigma, m[1, -1], n[1, -1])
        A_cf_off = ((turbine_height_off / 50) * turbine_height_off /
                    A_gradient_height) ** A_hellmann / resizem(Sigma, m[1, -1], n[1, -1])
        del A_gradient_height, Sigma
    else:
        A_cf_on = (turbine_height_on / 50) ** A_hellmann
        A_cf_off = (turbine_height_on / 50) ** A_hellmann
    del A_hellmann
    with rasterio.open(paths["LAND"]) as src:
        A_land = np.flipud(src.read(1)).astype(int)
    A_cf_on = A_cf_on * A_land
    del A_land
    if topo_correction:
        with rasterio.open(paths["TOPO"]) as src:
            A_topo = np.flipud(src.read(1)).astype(float)
        (a, b) = param["WindOn"]["resource"]["topo_factors"]
        A_cf_on = A_cf_on * np.minimum(np.exp(a * A_topo + b), 3.5)
    with rasterio.open(paths["EEZ"]) as src:
        A_eez = np.flipud(src.read(1)).astype(int)
    A_cf = A_cf_off * A_eez + A_cf_on

    array2raster(paths["CORR"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_cf)
    print("files saved: " + paths["CORR"])


def calculate_FLH(paths, param, tech):
    nproc = param["nproc"]
    year = param["year"]
    res = param["res"]
    GeoRef = param["GeoRef"]
    # H = param["Calculated_hours"]

    if tech == "WindOff":
        regions_shp = param["regions_eez"]
        nRegions = param["nRegions_eez"]
        Crd = param["Crd"][-nRegions-1:-1, :]
        m = param["m"][:, -nRegions-1:-1]
        n = param["n"][:, -nRegions-1:-1]
    else:
        regions_shp = param["regions_land"]
        nRegions = param["nRegions_land"]
        Crd = param["Crd"][0:nRegions, :]
        Ind = param["Ind"][0:nRegions, :]
        m = param["m"][:, 0:nRegions]
        n = param["n"][:, 0:nRegions]

    list_hours = np.array_split(np.arange(0, 8760), nproc)

    for reg in range(0, nRegions):

        region_name = regions_shp["NAME_SHORT"][reg]

        rasterData = {}
        # A_region
        rasterData["A_region"] = calc_region(regions_shp.iloc[reg], Crd[reg, :], res, GeoRef)

        # results = calc_FLH_solar(range(0,20), [reg, paths, param, nRegions, region_name, rasterData, tech])
        if tech in ['PV', 'CSP']:
            CLEARNESS = hdf5storage.read('CLEARNESS', paths["CLEARNESS"])
            day_filter = CLEARNESS[Ind[0, reg, 2] - 1:Ind[0, reg, 0], Ind[0, reg, 3] - 1:Ind[0, reg, 1], :].sum(
                axis=(0, 1)) != 0
            list_hours = np.arange(0, 8760)
            # results = calc_FLH_solar(list_hours[day_filter], [reg, paths, param, nRegions, region_name, rasterData, tech])
            list_hours = np.array_split(list_hours[day_filter], nproc)
            print(len(list_hours[0]))
            param["status_bar_limit"] = list_hours[0][-1]
            results = Pool(processes=nproc).starmap(calc_FLH_solar, product(list_hours, [
                [reg, paths, param, nRegions, region_name, rasterData, tech]]))
        elif tech in ['WindOn', 'WindOff']:
            print(len(list_hours[0]))
            param["status_bar_limit"] = list_hours[0][-1]
            results = Pool(processes=nproc).starmap(calc_FLH_wind, product(list_hours, [
                [reg, paths, param, nRegions, region_name, rasterData, tech]]))

        # Collecting results
        TS = np.zeros((8760, 1))
        FLH = np.zeros((m[1, reg], n[1, reg]))
        for p in range(len(results)):
            FLH = FLH + results[p][0]
            TS = TS + results[p][1]

        FLH[FLH == 0] = np.nan
        TS[np.isnan(TS)] = 0

        path_FLH_reg = paths["OUT"] + tech + '_FLH_' + region_name + '_' + year + '.mat'
        path_TS_mean_reg = paths["OUT"] + tech + '_Mean_Timeseries_' + region_name + '_' + year + '.mat'

        hdf5storage.writes({'FLH': FLH}, path_FLH_reg, store_python_metadata=True, matlab_compatible=True)
        print("files saved: " + path_FLH_reg)

        hdf5storage.writes({'TS_mean': TS}, path_TS_mean_reg, store_python_metadata=True, matlab_compatible=True)
        print("files saved: " + path_TS_mean_reg)

        del FLH, TS


def combine_FLH(paths, param, tech):
    if tech == "WindOff":
        regions_shp = param["regions_eez"]
        nRegions = param["nRegions_eez"]
        Ind = param["Ind"][:, -nRegions-1:-1, :]
    else:
        regions_shp = param["regions_land"]
        nRegions = param["nRegions_land"]
        Ind = param["Ind"][:, 0:nRegions, :]

    m = param["m"][1, -1]
    n = param["n"][1, -1]
    year = param["year"]
    S = ', '.join(list_files(paths["OUT"], tech + '_FLH_' + '*.mat'))
    if len(S):
        FLH_all = np.zeros((m, n))
        for reg in range(0, nRegions):
            region_name = regions_shp["NAME_SHORT"][reg]
            region_path = tech + '_FLH_' + region_name + '_' + year + '.mat'
            if region_path in S:
                FLH = hdf5storage.read('FLH', paths["OUT"] + region_path)
                FLH_all[Ind[1, reg, 2] - 1:Ind[1, reg, 0], Ind[1, reg, 3] - 1:Ind[1, reg, 1]] = \
                    np.fmax(FLH_all[Ind[1, reg, 2] - 1:Ind[1, reg, 0], Ind[1, reg, 3] - 1:Ind[1, reg, 1]], FLH)

        hdf5storage.writes({'FLH': FLH_all}, paths[tech]["FLH"], store_python_metadata=True, matlab_compatible=True)
        print("files saved: " + paths[tech]["FLH"])


def masking(paths, param, tech):
    mask = param[tech]["mask"]
    GeoRef = param["GeoRef"]

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
        A_suitability_lu = 1
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
        Georef = param["Georef"]
        array2raster(changeExt2tif(paths[tech]["mask"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelheight"],
                     A_mask)
        print("files saved:" + changeExt2tif(paths[tech]["mask"]))

        array2raster(changeExt2tif(paths[tech]["FLH_mask"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelheight"],
                     FLH_mask)
        print("files saved:" + changeExt2tif(paths[tech]["FLH_mask"]))


def weighting(paths, param, tech):
    weight = param[tech]["weight"]
    Crd = param["Crd"]
    m = param["m"]
    n = param["n"]
    res = param["res"]
    GeoRef = param["GeoRef"]

    if tech == 'PV':
        # Ground Cover Ratio - defines spacing between PV arrays
        A_GCR = calc_gcr(Crd[-1, :][np.newaxis], m[1, -1], n[1, -1], res, weight["GCR"])
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

    # Calculate available areas
    A_area = calc_areas(Crd, n, res, -1) * A_availability

    # Weighting matrix for the energy output (technical potential) in MWp
    A_weight = A_area * A_GCR * weight["power_density"] * weight["f_performance"]

    # Calculate weighted FLH in MWh
    FLH = hdf5storage.read('FLH', paths[tech]["FLH"])
    FLH_weight = FLH * A_weight

    # Save HDF5 Files
    hdf5storage.writes({'A_area': A_area}, paths[tech]["area"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths[tech]["area"])
    hdf5storage.writes({'A_weight': A_weight}, paths[tech]["weight"], store_python_metadata=True,
                       matlab_compatible=True)
    print("files saved: " + paths[tech]["weight"])
    hdf5storage.writes({'FLH_weight': FLH_weight}, paths[tech]["FLH_weight"], store_python_metadata=True,
                       matlab_compatible=True)
    print("files saved: " + paths[tech]["FLH_weight"])

    # Save GEOTIFF files
    if param["savetiff"]:
        array2raster(changeExt2tif(paths[tech]["area"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelheight"],
                     A_area)
        print("files saved:" + changeExt2tif(paths[tech]["area"]))

        array2raster(changeExt2tif(paths[tech]["weight"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelheight"],
                     A_weight)
        print("files saved:" + changeExt2tif(paths[tech]["weight"]))

        array2raster(changeExt2tif(paths[tech]["FLH_weight"]),
                     GeoRef["RasterOrigin"],
                     GeoRef["pixelWidth"],
                     GeoRef["pixelheight"],
                     FLH_weight)
        print("files saved:" + changeExt2tif(paths[tech]["FLH_weight"]))


# #########
# # Module : reporting
# #########

# # Compute Sums

# A_area_mask = A_area * A_mask
# area_sum_km_2 = np.nansum(A_area_mask) / 1e6
# power_potential_sum_TWp = np.nansum(A_weight) / 1e6
# energy_potential_sum_TWh = np.nansum(A_FLH_weight) / 1e6


def find_locations_quantiles(paths, param, tech):
    FLH_mask = hdf5storage.read('FLH_mask', paths[tech]["FLH_mask"])
    quantiles = param["quantiles"]
    res = param["res"]
    GeoRef = param["GeoRef"]

    if tech == "WindOff":
        regions_shp = param["regions_eez"]
        nRegions = param["nRegions_eez"]
        Crd = param["Crd"][-nRegions - 1:, :]
        Ind = param["Ind"][:, -nRegions - 1:-1, :]
    else:
        regions_shp = param["regions_land"]
        nRegions = param["nRegions_land"]
        Crd = param["Crd"][0:nRegions, :]
        Crd = np.r_[Crd, param["Crd"][-1, :]]
        Ind = param["Ind"][:, 0:nRegions, :]

    Crd_Locations = np.zeros((nRegions, len(quantiles), 4))
    Ind_Locations = np.zeros((2, nRegions, len(quantiles), 4))
    region_names = []
    for reg in range(0, nRegions):
        region_names.append(regions_shp["NAME_SHORT"][reg])

        # A_region
        A_region = calc_region(regions_shp.iloc[reg], Crd[reg, :], res, GeoRef)

        FLH_reg = A_region * FLH_mask[Ind[1, reg, 2] - 1:Ind[1, reg, 0], Ind[1, reg, 3] - 1:Ind[1, reg, 1]]
        FLH_reg[FLH_reg == 0] = np.nan
        X = FLH_reg.flatten(order='F')
        I_old = np.argsort(X)

        # ESCAPE FOR LOOP IF INTERSECTION WITH RASTER ONLY YIELDS NAN
        if sum(np.isnan(X).astype(int)) == len(X):
            # do something
            continue

        for q in range(0, len(quantiles)):
            if quantiles[q] == 100:
                I = I_old[(len(X) - 1) - sum(np.isnan(X).astype(int))]
            elif quantiles[q] == 0:
                I = I_old[0]
            else:
                I = I_old[int(np.round(quantiles[q] / 100 * (len(X) - 1 - sum(np.isnan(X).astype(int)))))]

            # Convert the indices to row-column indices
            I, J = ind2sub(FLH_reg.shape, I)
            Ind_Locations[1, reg, q, :] = [(I + Ind[1, reg, 2]), (J + Ind[1, reg, 3]),
                                           (I + Ind[1, reg, 2]), (J + Ind[1, reg, 3])]
            Crd_Locations[reg, q, :] = crd_exact_high(np.squeeze(Ind_Locations[1, reg, q, :]).T, Crd[-1, :][np.newaxis],
                                                      res)

        Crd_Locations[reg, :, 2:4] = repmat(Crd[-1, 2:4], len(quantiles), 1)
        Ind_Locations[0, reg, :, :] = ind_merra(crd_merra_low(np.squeeze(Crd_Locations[reg, :, :]), res), res[0, :])

    param["Crd_Locations"] = Crd_Locations
    param["Ind_Locations"] = Ind_Locations

    # Format point locations
    points = list(map(tuple, np.reshape(Crd_Locations[:, :, [1, 0]], (-1, 2), 'C')))
    attribute1 = np.reshape(np.tile(region_names, (len(quantiles), 1)), (-1, 1), 'F')
    attribute2 = ['q' + str(i) for i in np.tile(quantiles, nRegions)]
    schema = {'geometry': 'Point', 'properties': {'NAME_SHORT': 'str', 'quantile': 'str'}}

    # Create Shapefile
    with fiona.open(paths[tech]["Locations"], 'w', 'ESRI Shapefile', schema) as c:
        c.writerecords([{'geometry': mapping(Point(points[i])),
                         'properties': {'NAME_SHORT': attribute1[i], 'quantile': attribute2[i]}} for i in
                        range(0, len(points))])
    print("files saved: " + paths[tech]["Locations"])
    return param


def generate_time_series(paths, param, tech):
    return
    # # Calculate A matrices
    # # Landuse classes 0-16, to be reclassified
    # with rasterio.open(paths["LU"]) as src:
    # w = src.read(1, window=rasterio.windows.Window.from_slices(slice((m[1, 0] - Ind[1, reg, 0]),
    # (m[1, 0] - Ind[1, reg, 2] + 1)),
    # slice(Ind[1, reg, 3] - 1,
    # Ind[1, reg, 1])))
    # rasterData["A_lu"] = np.flipud(w)

    # if technology == 'PV' or technology == 'CSP':
    # # Calculate special A matrices
    # with rasterio.open(paths["TOPO"]) as src:
    # w = src.read(1, window=rasterio.windows.Window.from_slices(slice((m[1, 0] - Ind[1, reg, 0]),
    # (m[1, 0] - Ind[1, reg, 2] + 1)),
    # slice(Ind[1, reg, 3] - 1,
    # Ind[1, reg, 1])))
    # rasterData["A_topo"] = np.flipud(w)

    # # Temperature coefficients for heating losses
    # rasterData["A_Ross"] = changem(rasterData["A_lu"], landuse["Ross_coeff"], landuse["type"]).astype(
    # float) / 10000
    # # Reflectivity coefficients
    # rasterData["A_albedo"] = changem(rasterData["A_lu"], landuse["albedo"], landuse["type"]).astype(float) / 100

    # # Calculate CLR_Mean and CLR_MAX
    # _, CLR_MAX = calc_clearness(merraData, reg, Ind)

    # elif technology == 'Wind':
    # rasterData["A_topo"] = rasterData["A_topo"][Locations_ind_high[:,0], Locations_ind_high[:,1]]
    # rasterData["A_cf"] = A_cf[tech][np.ix_(Ind[3:, reg, 0] - 1, Ind[3:, reg, 1] - 1)] * \
    # np.exp(a * rasterData["A_topo"] / 5000 - b)
    # rasterData["A_cf"] = np.diagonal(rasterData["A_cf"])

    # TS = np.zeros((len(quantiles), 8760))

    # for hour in range(0, 8760):
    # # Show progress of the simulation
    # print(str(reg) + '/' + str(nRegions - 1) + ' ' + region_name + ' ' + str(hour + 1))

    # CF = 0
    # if technology == 'PV':
    # CF, _ = calc_CF_solar(hour, reg, Ind, Crd, res, CLR_MAX, merraData, rasterData, pv, Locations)
    # elif technology == 'CSP':
    # _, CF = calc_CF_solar(hour, reg, Ind, Crd, res, CLR_MAX, merraData, rasterData, pv, Locations)
    # elif technology == 'Wind':
    # # Load MERRA data, increase its resolution, and fit it to the extent
    # w50m_h = W50M[:, :, hour]
    # w50m_h = w50m_h[np.ix_(Locations_ind_low[:, 0] - 1, Locations_ind_low[:, 1] - 1)]
    # w50m_h = np.diagonal(w50m_h)

    # # Calculate hourly capacity factor
    # CF = calc_CF_wind(w50m_h, rasterData, turbine, tech)

    # TS[:, hour] = CF

    # TS[np.isnan(TS)] = 0

    # paths["TS_quantiles"] = paths["OUT"] + technology + '_Timeseries_' + region_name + '_' + year + '.mat'
    # paths["Locations"] = paths["OUT"] + technology + '_Locations_' + region_name + '_' + year + '.mat'
    # description["paths"] = paths

    # with h5py.File(paths["TS_quantiles"], 'w') as f:
    # f.create_dataset('TS', data=TS)
    # recursive_dict_save(f, 'description/', description)
    # print("files saved:" + paths["TS_quantiles"])

    # with h5py.File(paths["Locations"], 'w') as f:
    # f.create_dataset('Locations', data=Locations)
    # recursive_dict_save(f, 'description/', description)
    # print("files saved:" + paths["Locations"])
    # del TS


# # ## Timeseries for all regions

# S = np.array(list_files(paths["OUT"], technology + '_Timeseries_' + '*.mat'))
# if S.size != 0:
# timeseries = {}
# for f in range(0, len(S)):
# with h5py.File(paths["OUT"] + S[f], 'r') as src:
# TS = np.array(src["TS"])
# reg = src["description/region"][()] + 1
# suffix = ''

# if technology == 'PV':
# suffix = '.solar'
# elif technology == 'CSP':
# suffix = '.beam'
# elif technology == 'Wind':
# if regions_shp["Population"][reg - 1] == 0:
# suffix = '.WindOff'
# else:
# suffix = '.WindOn'

# timeseries[regions_shp["NAME_SHORT"][reg - 1] + suffix] = {}

# for q in range(0, len(quantiles)):
# newfield = 'q' + str(quantiles[q])
# timeseries[regions_shp["NAME_SHORT"][reg - 1] + suffix][newfield] = TS[q, :].T
# TS = timeseries
# del timeseries

# description["region"] = region
# paths["TS"] = paths["OUT"] + region + '_' + technology + '_Timeseries_' + year + '.mat'
# description["paths"] = paths

# with h5py.File(paths["TS"], 'w') as f:
# recursive_dict_save(f, 'Region/', TS)
# recursive_dict_save(f, 'description/', description)
# print("files saved:" + paths["TS"])


if __name__ == '__main__':
    paths, param = initialization()

    # Check if Merra2 mat files have been generated
    if not (os.path.isfile(paths["GHI"])
            and os.path.isfile((paths["TOA"]))
            and os.path.isfile(paths["T2M"])
            and os.path.isfile(paths["U50M"])
            and os.path.isfile(paths["V50M"])):
        generate_weather_files(paths)

    generate_landsea(paths, param)  # Land and Sea
    generate_landuse(paths, param)  # Landuse
    generate_bathymetry(paths, param)  # Bathymetry
    generate_topography(paths, param)  # Topography
    generate_slope(paths, param)  # Slope
    generate_population(paths, param)  # Population
    generate_protected_areas(paths, param)  # Protected areas
    # generate_buffered_population(paths, param)  # Buffered Population
    # generate_wind_correction(paths, param)  # Correction factors for wind speeds
    for tech in param["technology"]:
        calculate_FLH(paths, param, tech)
        combine_FLH(paths, param, tech)
        masking(paths, param, tech)
        weighting(paths, param, tech)
        # param = find_locations_quantiles(paths, param, tech)
