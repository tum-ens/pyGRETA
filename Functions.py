import os
# os.environ['MKL_NUM_THREADS'] = '8'
import math as m
import numpy as np
import gdal
import osr
import rasterio
from osgeo import ogr
from glob import glob
from os import getcwd, chdir
from rasterio import windows, mask
from numpy.matlib import repmat, reshape, sin, arcsin, cos, arccos, tan, arctan


#################################
# Move these to a separate file #
#################################
def sind(alpha):
    return sin(np.deg2rad(alpha))


def cosd(alpha):
    return cos(np.deg2rad(alpha))


def tand(alpha):
    return tan(np.deg2rad(alpha))


def arcsind(digit):
    return np.rad2deg(arcsin(digit))


def arccosd(digit):
    return np.rad2deg(arccos(digit))


####################################

def calc_ext(regb, ext, res):
    minRow = m.floor(regb["miny"] / res[1, 0]) * res[1, 0]
    maxRow = m.ceil(regb["maxy"] / res[1, 0]) * res[1, 0]
    minCol = m.floor(regb["minx"] / res[1, 1]) * res[1, 1]
    maxCol = m.ceil(regb["maxx"] / res[1, 1]) * res[1, 1]

    return [[min(m.ceil((ext[0, 0] - res[0, 0] / 2) / res[0, 0]) * res[0, 0] + res[0, 0] / 2, maxRow),
             min(m.ceil((ext[0, 1] - res[0, 1] / 2) / res[0, 1]) * res[0, 1] + res[0, 1] / 2, maxCol),
             max(m.ceil((ext[0, 2] - res[0, 0] / 2) / res[0, 0]) * res[0, 0] + res[0, 0] / 2, minRow),
             max(m.ceil((ext[0, 3] - res[0, 1] / 2) / res[0, 1]) * res[0, 1] + res[0, 1] / 2, minCol)]]


def crd_merra_low(Ext, res):
    Crd = np.array([(np.ceil((Ext[:, 0] - res[0, 0] / 2) / res[0, 0]) * res[0, 0] + res[0, 0] / 2),
                    (np.ceil((Ext[:, 1] - res[0, 1] / 2) / res[0, 1]) * res[0, 1] + res[0, 1] / 2),
                    (np.floor((Ext[:, 2] + res[0, 0] / 2) / res[0, 0]) * res[0, 0] - res[0, 0] / 2),
                    (np.floor((Ext[:, 3] + res[0, 1] / 2) / res[0, 1]) * res[0, 1] - res[0, 1] / 2)])
    Crd = Crd.T
    return Crd


def crd_exact_high(Ind, Ext, res):
    Ind = Ind[np.newaxis]

    Crd = [Ind[:, 0] * res[1, 0] + Ext[0, 2],
           Ind[:, 1] * res[1, 1] + Ext[0, 3],
           (Ind[:, 2] - 1) * res[1, 0] + Ext[0, 2],
           (Ind[:, 3] - 1) * res[1, 1] + Ext[0, 3]]
    return Crd


def ind_merra_low(Crd, res):
    ind = np.array([(Crd[:, 0] - Crd[0, 2]) / res[0, 0],
                    (Crd[:, 1] - Crd[0, 3]) / res[0, 1],
                    (Crd[:, 2] - Crd[0, 2]) / res[0, 0] + 1,
                    (Crd[:, 3] - Crd[0, 3]) / res[0, 1] + 1])
    ind = np.transpose(ind)
    ind = ind.astype(int)
    return ind


def ind_merra_high(Crd, res):
    ind = np.array([(Crd[:, 0] - Crd[0, 2]) / res[1, 0],
                    (Crd[:, 1] - Crd[0, 3]) / res[1, 1],
                    (Crd[:, 2] - Crd[0, 2]) / res[1, 0] + 1,
                    (Crd[:, 3] - Crd[0, 3]) / res[1, 1] + 1])
    ind = np.transpose(ind)
    ind = ind.astype(int)
    return ind


def ind_global(Ext_PV, res):
    ind = np.array([np.round((90 - Ext_PV[:, 0]) / res[0]) + 1,
                    np.round((180 + Ext_PV[:, 1]) / res[1]),
                    np.round((90 - Ext_PV[:, 2]) / res[0]),
                    np.round((180 + Ext_PV[:, 3]) / res[1]) + 1])
    ind = np.transpose(ind)
    ind = ind.astype(int)
    return ind


def calc_geotiff(Crd, res):
    R1 = {"RasterOrigin": [Crd[0, 3], Crd[0, 0]],
          "RasterOrigin_alt": [Crd[0, 3], Crd[0, 2]],
          "pixelWidth": res[1, 1],
          "pixelheight": -res[1, 0]}
    return R1


def calc_region(region, Crd, res, R1):

    latlim = Crd[2] - Crd[0]
    lonlim = Crd[3] - Crd[1]
    M = int(m.fabs(latlim) / res[1, 0])
    N = int(m.fabs(lonlim) / res[1, 1])
    A_region = np.ones((M, N))
    origin = [Crd[3], Crd[2]]

    array2raster("test.tiff", origin, R1["pixelWidth"], -R1["pixelheight"], A_region)

    if region.geometry.geom_type == 'MultiPolygon':
        features = [feature for feature in region.geometry]
    else:
        features = [region.geometry]

    with rasterio.open("test.tiff", "r") as src:
        out_image, out_transform = mask.mask(src, features, crop=False, nodata=0, all_touched=False, filled=True)
    A_region = out_image[0]
    os.remove("test.tiff")
    return A_region


def generate_landsea(paths, regions_shp, Crd, Ind, m, n, res, R1):
    nRegions = len(regions_shp) + 1
    if not os.path.isfile(paths["LAND"]):
        A_eez = np.zeros((m[1, 0], n[1, 0]))
        A_land = np.zeros((m[1, 0], n[1, 0]))
        for reg in range(1, nRegions):

            A_region = calc_region(regions_shp.iloc[reg - 1], Crd[reg, :], res, R1)

            if regions_shp['Population'][reg - 1] == 0:
                A_eez[(Ind[1, reg, 2] - 1):Ind[1, reg, 0], (Ind[1, reg, 3] - 1): Ind[1, reg, 1]] = \
                    A_eez[(Ind[1, reg, 2] - 1):Ind[1, reg, 0], (Ind[1, reg, 3] - 1):Ind[1, reg, 1]] + A_region
                array2raster(paths["EEZ"], R1["RasterOrigin"], R1["pixelWidth"], R1["pixelheight"], A_eez)
                print("files saved:" + paths["EEZ"])
            else:
                A_land[(Ind[1, reg, 2] - 1):Ind[1, reg, 0], (Ind[1, reg, 3] - 1):Ind[1, reg, 1]] = \
                    A_land[(Ind[1, reg, 2] - 1):Ind[1, reg, 0], (Ind[1, reg, 3] - 1):Ind[1, reg, 1]] + A_region
                array2raster(paths["LAND"], R1["RasterOrigin"], R1["pixelWidth"], R1["pixelheight"], A_land)
                print("files saved:" + paths["LAND"])


def generate_landuse(paths, Ind, R1):
    if not os.path.isfile(paths['LU']):
        with rasterio.open(paths["LU_global"]) as src:
            w = src.read(1, window=windows.Window.from_slices(slice(Ind[2, 0, 0]-1, Ind[2, 0, 2]),
                                                              slice(Ind[2, 0, 3]-1, Ind[2, 0, 1])))
        w = np.flipud(w)
        array2raster(paths["LU"], R1["RasterOrigin"], R1["pixelWidth"], R1["pixelheight"], w)
        print("files saved:" + paths["LU"])


def generate_bathymetry(paths, Ind, R1):
    if not os.path.isfile(paths['BATH']):
        with rasterio.open(paths["Bathym_global"]) as src:
            A_BATH = src.read(1)
        A_BATH = resizem(A_BATH, 180 * 240, 360 * 240)
        A_BATH = np.flipud(A_BATH[Ind[2, 0, 0] - 1: Ind[2, 0, 2], Ind[2, 0, 3] - 1: Ind[2, 0, 1]])
        array2raster(paths['BATH'], R1["RasterOrigin"], R1["pixelWidth"], R1["pixelheight"], A_BATH)
        print("files saved:" + paths["BATH"])


def generate_topography(paths, Ind, R1):
    if not os.path.isfile(paths["TOPO"]):
        Topo = np.zeros((180 * 240, 360 * 240))
        tile_extents = np.zeros((24, 4), dtype=int)
        i = 1
        j = 1
        for letter in char_range('A', 'X'):
            north = (i - 1) * 45 * 240 + 1
            east = j * 60 * 240
            south = i * 45 * 240
            west = (j - 1) * 60 * 240 + 1
            tile_extents[ord(letter), :] = [north, east, south, west]
            j = j + 1
            if j == 7:
                i = i + 1
                j = 1
        n_min = (Ind[2, 0, 0] // (45 * 240)) * 45 * 240 + 1
        e_max = (Ind[2, 0, 1] // (60 * 240) + 1) * 60 * 240
        s_max = (Ind[2, 0, 2] // (45 * 240) + 1) * 45 * 240
        w_min = (Ind[2, 0, 3] // (60 * 240)) * 60 * 240 + 1

        need = np.logical_and((np.logical_and((tile_extents[:, 0] >= n_min), (tile_extents[:, 1] <= e_max))),
                              np.logical_and((tile_extents[:, 2] <= s_max), (tile_extents[:, 3] >= w_min)))

        for letter in char_range('A', 'X'):
            index = ord(letter) - ord('A')
            if need[index]:
                with rasterio.open(paths["LU_global"] + '15-' + letter + '.tif') as src:
                    tile = src.read()
                Topo[tile_extents[index, 0]-1: tile_extents[index, 2], tile_extents[index, 3]-1: tile_extents[index, 1]] = \
                    tile[0:-1, 0:-1]

        A_TOPO = np.flipud(Topo[Ind[2,0,0]-1:Ind[2,0,2], Ind[2,0,3]-1:Ind[2,0,1]])
        array2raster(paths["TOPO"], R1["RasterOrigin"], R1["pixelWidth"], R1["pixelheight"], A_TOPO)
        print("files saved:" + paths["TOPO"])


def generate_slope(paths, Ind, R1):
    if not os.path.isfile(paths["SLOPE"]):
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
        x_cell = x_cell[Ind[2, 0, 0] - 2:Ind[2, 0, 2]+1, Ind[2, 0, 3] - 2: Ind[3, 0, 1]+1]
        x_cell = np.flipud(x_cell)

        y_cell = repmat((deltaLat * m_per_deg_lat), 360 * 240, 1).T
        y_cell = y_cell[Ind[2, 0, 0] - 2:Ind[2, 0, 2]+1, Ind[2, 0, 3] - 2: Ind[3, 0, 1]+1]
        y_cell = np.flipud(y_cell)

        with rasterio.open(paths["TOPO"]) as src:
            A_TOPO = src.read(1)

        topo_a = np.zeros((Ind[2, 0, 2] - Ind[2, 0, 0] + 3, Ind[2, 0, 1] - Ind[2, 0, 3] + 3))
        topo_a[0: - 2, 0: - 2] = A_TOPO

        topo_b = np.zeros((Ind[2, 0, 2] - Ind[2, 0, 0] + 3, Ind[2, 0, 1] - Ind[2, 0, 3] + 3))
        topo_b[0: - 2, 1: - 1] = A_TOPO

        topo_c = np.zeros((Ind[2, 0, 2] - Ind[2, 0, 0] + 3, Ind[2, 0, 1] - Ind[2, 0, 3] + 3))
        topo_c[0: - 2, 2: - 0] = A_TOPO

        topo_d = np.zeros((Ind[2, 0, 2] - Ind[2, 0, 0] + 3, Ind[2, 0, 1] - Ind[2, 0, 3] + 3))
        topo_d[1: - 1, 0: - 2] = A_TOPO

        topo_f = np.zeros((Ind[2, 0, 2] - Ind[2, 0, 0] + 3, Ind[2, 0, 1] - Ind[2, 0, 3] + 3))
        topo_f[1: - 1, 2: - 0] = A_TOPO

        topo_g = np.zeros((Ind[2, 0, 2] - Ind[2, 0, 0] + 3, Ind[2, 0, 1] - Ind[2, 0, 3] + 3))
        topo_g[2: - 1, 0: - 2] = A_TOPO

        topo_h = np.zeros((Ind[2, 0, 2] - Ind[2, 0, 0] + 3, Ind[2, 0, 1] - Ind[2, 0, 3] + 3))
        topo_h[2: - 1, 1: - 1] = A_TOPO

        topo_i = np.zeros((Ind[2, 0, 2] - Ind[2, 0, 0] + 3, Ind[2, 0, 1] - Ind[2, 0, 3] + 3))
        topo_i[2: - 1, 2: - 0] = A_TOPO

        dzdx = ((topo_c + 2 * topo_f + topo_i) - (topo_a + 2 * topo_d + topo_g)) / (8 * x_cell)
        dzdy = ((topo_g + 2 * topo_h + topo_i) - (topo_a + 2 * topo_b + topo_c)) / (8 * y_cell)

        slope_deg = arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / np.pi
        slope_pc = tan(np.deg2rad(slope_deg)) * 100

        A_SLP = slope_pc[1:-1, 1:-1]
        array2raster(paths["SLOPE"], R1["RasterOrigin"], R1["pixelWidth"], R1["pixelheight"], A_SLP)
        print("files saved:" + paths["SLOPE"])


def generate_population(paths, Ind, R1):
    if not os.path.isfile(paths["POP"]):

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
        n_min = (Ind[2, 0, 0] // (45 * 240)) * 45 * 240 + 1
        e_max = (Ind[2, 0, 1] // (60 * 240) + 1) * 60 * 240
        s_max = (Ind[2, 0, 2] // (45 * 240) + 1) * 45 * 240
        w_min = (Ind[2, 0, 3] // (60 * 240)) * 60 * 240 + 1

        need = np.logical_and((np.logical_and((tile_extents[:, 0] >= n_min), (tile_extents[:, 1] <= e_max))),
                              np.logical_and((tile_extents[:, 2] <= s_max), (tile_extents[:, 3] >= w_min)))

        for letter in char_range('A', 'X'):
            index = ord(letter) - ord('A')
            if need[index]:
                with rasterio.open(paths["Pop_tiles"] + letter + '.tif') as src:
                    tile = src.read()
                Pop[tile_extents[index, 0]-1: tile_extents[index, 2], tile_extents[index, 3]-1: tile_extents[index, 1]] = \
                    tile[0]
        A_POP = np.flipud(Pop[Ind[2, 0, 0]-1:Ind[2, 0, 2], Ind[2, 0, 3]-1:Ind[2, 0, 1]])
        array2raster(paths["POP"], R1["RasterOrigin"], R1["pixelWidth"], R1["pixelheight"], A_POP)
        print("files saved:" + paths["POP"])


def generate_protected_areas(paths, protected_areas, R1):
    if not os.path.isfile(paths["PA"]):
        # print("Protected areas:")
        # set up protected areas dictionary
        protection_type = dict(zip(protected_areas["IUCN_Category"], protected_areas["pa_type"]))

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

        # How many fields are in the shapefile, and what are their names?
        # First we need to capture the layer definition
        defn = layer.GetLayerDefn()

        # How many fields
        field_count = defn.GetFieldCount()
        # print('Layer has {n} fields'.format(n=field_count))

        # What are their names?
        # print('Their names are: ')
        for i in range(field_count):
            field_defn = defn.GetFieldDefn(i)
            print('\t{name} - {datatype}'.format(name=field_defn.GetName(),
                                                 datatype=field_defn.GetTypeName()))

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


def generate_buffered_population(paths, buffer_pixel_amount, R1):
    if not os.path.isfile(paths["PopBuff"]):
        with rasterio.open(paths["LU"]) as src:
            A_lu = src.read(1)
            A_lu = np.flipud(A_lu).astype(int)

        A_lu_buffered = create_buffer(A_lu, buffer_pixel_amount)
        A_notPopulated = (~A_lu_buffered).astype(int)

        array2raster(paths["PopBuff"], R1["RasterOrigin"], R1["pixelWidth"], R1["pixelheight"], A_notPopulated)
        print("files saved:" + paths["PopBuff"])

def sumnorm_MERRA2(A, m, n, res):
    s = np.zeros((m, n))
    row_step = int(res[0, 0] / res[1, 0])
    col_step = int(res[0, 1] / res[1, 1])
    for i in range(0, m):
        for j in range(0, n):
            s[i, j] = np.sum(A[(row_step * i):(row_step * (i + 1)),
                             (col_step * j):(col_step * (j + 1))]) / (row_step * col_step)
    return s


def masking(FLH_all, mask, landuse, paths, technology, windtechnology, buffer_pixel_amount):
    if technology == 'PV' or technology == 'CSP':
        return maskingtemp(FLH_all, mask, landuse, paths, technology, buffer_pixel_amount)
    elif technology == 'Wind':
        on = [0, 0]
        off = [0, 0]
        if 'Onshore' in windtechnology:
            # Onshore landuse variables
            landuse.update(landuse["Onshore"])
            landuse["tech"] = "Onshore"
            on[0], on[1] = maskingtemp(FLH_all, mask, landuse, paths, technology, buffer_pixel_amount)
        if 'Offshore' in windtechnology:
            # Offshore landuse variables
            landuse.update(landuse["Offshore"])
            landuse["tech"] = "Offshore"
            off[0], off[1] = maskingtemp(FLH_all, mask, landuse, paths, technology, buffer_pixel_amount)
        # Summing up the results if both technologies are considered
        return (on[0] + off[0]), (on[0] + off[0]) * FLH_all


def maskingtemp(FLH_all, mask, landuse, paths, technology, buffer_pixel_amount):
    if technology == 'PV' or technology == 'CSP':
        suitability = landuse["suit_s"]
        slope = mask["slope_s"]
    elif technology == 'Wind':
        suitability = landuse["suit_w"]
        slope = mask["slope_w"]
        # Separate land and sea areas
        with rasterio.open(paths["EEZ"]) as src:
            A_eez = src.read(1)
            A_eez = np.flipud(A_eez)
        with rasterio.open(paths["LAND"]) as src:
            A_land = src.read(1)
            A_land = np.flipud(A_land)
        A_eez[A_land.astype(bool)] = 0

    with rasterio.open(paths["PA"]) as src:
        A_protect = src.read(1)
        A_protect = 1 - np.flipud(A_protect)  # Is 1 if non-protected area (boolean)
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
        A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified
    with rasterio.open(paths["SLOPE"]) as src:
        A_slope = src.read(1)
        A_slope = np.flipud(A_slope) <= slope  # Is 1 for slopes <'slope'% (boolean)
    with rasterio.open(paths["PopBuff"]) as src:
        A_notPopulated = src.read(1)
        A_notPopulated = np.flipud(A_notPopulated).astype(bool)  # Is 1 for slopes <'slope'% (boolean)

    # Get buffer around populated areas

    # A_lu_buffered = create_buffer(A_lu, buffer_pixel_amount)
    # A_notPopulated = ~A_lu_buffered

    if technology == 'PV' or technology == 'CSP':
        # Exclude landuse types that are not suitable for PV applications
        A_suitability = changem(A_lu.astype(float), suitability, landuse["type"].astype(float))
        # Masking matrix for the suitable sites (pixels)
        A_mask = A_protect.astype(int) * (A_suitability > 0).astype(int) * A_slope.astype(int) * A_notPopulated.astype(int)
    elif technology == 'Wind':
        if landuse["tech"] == 'Onshore':
            A_suitability = changem(A_lu.astype(float), suitability, landuse["type"].astype(float))

            A_mask = A_notPopulated.astype(int) * A_protect.astype(float) * (A_suitability > 0).astype(float) * \
                     ((A_land * A_slope) > 0).astype(float)
        else:
            with rasterio.open(paths["BATH"]) as src:
                A_bathymetry = src.read(1)
                A_bathymetry = np.flipud(A_bathymetry) >= mask["depth"]  # Is 1 for depths >= 'depth'% (boolean)

            # Masking matrix for the available sites (pixels)
            A_mask = ((A_eez * A_bathymetry) > 0).astype(float)

        A_mask[A_mask > 1] = 1

    FLH_mask = FLH_all * A_mask
    FLH_mask[FLH_mask == 0] = np.nan

    return A_mask, FLH_mask


def weighting(FLH_all, weight, landuse, paths, technology, windtechnology, Crd, n, res, pa_table):
    if technology == 'PV' or technology == 'CSP':
        return weightingtemp(FLH_all, weight, landuse, paths, technology, Crd, n, res, pa_table)
    elif technology == 'Wind':
        on = [0, 0, 0]
        off = [0, 0, 0]
        if 'Onshore' in windtechnology:
            # Onshore landuse and weighting parameters
            w = weight["Onshore"]
            landuse.update(landuse["Onshore"])
            landuse["tech"] = "Onshore"

            on[0], on[1], on[2] = weightingtemp(FLH_all, w, landuse, paths, technology, Crd, n, res, pa_table)
        if 'Offshore' in windtechnology:
            # Offshore landuse and weighting parameters
            w = weight["Offshore"]
            landuse.update(landuse["Offshore"])
            landuse["tech"] = "Offshore"

            off[0], off[1], on[2] = weightingtemp(FLH_all, w, landuse, paths, technology, Crd, n, res, pa_table)
        # Summing up the results if both technologies are considered
        return (on[0] + off[0]), (on[1] + off[1]), (on[2] + off[2])


def weightingtemp(FLH_all, weight, landuse, paths, technology, Crd, n, res, pa_table):

    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
        A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified

    if technology == 'PV':
        availability = landuse["avail_s"]
        f_pd = weight["f_pd_pv"]
        f_performance = weight["f_performance_pv"]

        # Ground Cover Ratio - defines spacing between PV arrays
        A_GCR = calc_gcr(Crd[0, :][np.newaxis], res, weight["GCR"])

    elif technology == 'CSP':
        availability = landuse["avail_s"]
        f_pd = weight["f_pd_csp"]
        f_performance = weight["f_performance_csp"]
        A_GCR = 1  # TO BE CHANGED!
    elif technology == 'Wind':
        availability = landuse["avail_w"]
        f_pd = weight["f_pd_w"]
        f_performance = weight["f_performance_w"]
        A_GCR = 1  # irrelevant for wind, can be a function of the slope, to be changed!!!!

    with rasterio.open(paths["PA"]) as src:
        PA = np.flipud(src.read(1))
    A_notprotect = changem(PA.astype(float), pa_table["pa_availability"], pa_table["pa_type"])

    A_availability = (changem(A_lu, availability, landuse["type"]).astype(float)/100) * A_notprotect
    A_area = calc_areas(Crd, n, res, 1) * A_availability

    # Weighting matrix for the energy output (technical potential) in MWp
    A_weight = A_area * A_GCR * f_pd * f_performance

    # Now in MWh
    FLH_weight = FLH_all * A_weight

    return A_area, A_weight, FLH_weight


def calc_areas(Crd, n, res, reg):
    # WSG84 ellipsoid constants
    a = 6378137  # major axis
    b = 6356752.3142  # minor axis
    reg = reg - 1
    e = np.sqrt(1 - (b / a) ** 2)

    # Lower pixel latitudes
    lat_vec = np.arange(Crd[reg, 2], Crd[reg, 0], res[1, 0])
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
    f_upper = np.deg2rad(lat_vec + res[1, 0])

    zm_upper = 1 - (e * sin(f_upper))
    zp_upper = 1 + (e * sin(f_upper))

    upperSliceAreas = np.pi * b ** 2 * ((2 * np.arctanh((e * sin(f_upper)))) / (2 * e) +
                                        (sin(f_lower) / (zp_upper * zm_upper)))

    # Pixel areas
    # Finding the latitudinal pixel-sized globe slice areas then dividing them by the longitudinal pixel size
    area_vec = ((upperSliceAreas - lowerSliceAreas) * res[1, 1] / 360).T
    A_area = np.tile(area_vec, (1, n[1, reg]))
    return A_area


# ## Miscellaneous Functions

def hourofmonth():
    h = 24 * np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]).astype(int)
    for i in range(1, 12):
        h[i] = h[i] + h[i - 1]
    return h.astype(int)


def resizem(A_in, row_new, col_new):
    row_rep = row_new // np.shape(A_in)[0]
    col_rep = col_new // np.shape(A_in)[1]
    A_inf = (A_in.flatten(order='F')[np.newaxis])
    A_out = reshape(repmat(reshape(reshape(repmat((A_in.flatten(order='F')[np.newaxis]), row_rep, 1), (row_new, -1), order='F').T, (-1, 1), order='F'), 1, col_rep).T, (col_new, row_new), order='F').T

    return A_out


def array2raster(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float64, ['COMPRESS=PACKBITS'])
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(np.flipud(array))
    outband.FlushCache()
    outband = None


def char_range(c1, c2):
    # Generates the characters from `c1` to `c2`, inclusive.
    for c in range(ord(c1), ord(c2) + 1):
        yield chr(c)


def changem(A, newval, oldval):
    Out = np.zeros(A.shape)
    z = np.array((oldval, newval)).T
    for i, j in z:
        np.place(Out, A == i, j)
    return Out


def list_files(directory, format):
    saved = getcwd()
    chdir(directory)
    it = glob(format)
    chdir(saved)
    return it


def forloop():
    return 1


def recursive_dict_save(h5file, path, dic):
    for key, item in dic.items():
        if isinstance(item, (int, float, np.ndarray, np.int64, np.float64, str, bytes)):
            h5file[path + key] = item
        elif isinstance(item, dict):
            recursive_dict_save(h5file, path + key + '/', item)
        else:
            raise ValueError('Cannot save %s type' % type(item))


def sub2ind(array_shape, rows, cols):
    return np.ravel_multi_index((rows, cols), array_shape, order='F')


def ind2sub(array_shape, ind):
    return np.unravel_index(ind, array_shape, order='F')


def create_buffer(A_lu, buffer_pixel_amount):
    # A_lu matrix element values range from 0 to 16:
    # 0   -- Water
    # 1   -- Evergreen needle leaf forest
    # 2   -- Evergreen broad leaf forest
    # 3   -- Deciduous needle leaf forest
    # 4   -- deciduous broad leaf forest
    # 5   -- Mixed forests
    # 6   -- Closed shrublands
    # 7   -- Open shrublands
    # 8   -- Woody savannas
    # 9   -- Grasslands
    # 10  -- Permanent wetland
    # 12  -- Croplands
    # 13  -- URBAN AND BUILT-UP
    # 14  -- Croplands / natural vegetation mosaic
    # 15  -- Snow and ice
    # 16  -- Barren or sparsely vegetated

    # Mark the matrix elements with values 13
    A_lu = A_lu == 13

    # modify
    # create a buffer around the cities
    shifted_A_lu = A_lu

    for p in range(0, buffer_pixel_amount):
        n = 1  # Number of pixel shifts per loop
        shifted_left = superpose_left(shifted_A_lu, n)
        shifted_right = superpose_right(shifted_A_lu, n)
        shifted_up = superpose_up(shifted_A_lu, n)
        shifted_down = superpose_down(shifted_A_lu, n)

        superposed = shifted_left + shifted_right + shifted_up + shifted_down

        superposed = superposed != 0
        shifted_A_lu = superposed

    A_lu_buffered = shifted_A_lu
    return A_lu_buffered


def superpose_left(A_lu, buffer_pixed_amount):
    # shift the matrix to the left
    # shift amount is defined by buffer_pixel amount
    left = np.append(A_lu[:, buffer_pixed_amount:], np.zeros((A_lu.shape[0], buffer_pixed_amount)), axis=1)
    shifted_left = A_lu + left
    shifted_left = shifted_left != 0
    return shifted_left


def superpose_right(A_lu, buffer_pixed_amount):
    # shift the matrix to the right
    # shift amount is defined by buffer_pixel amount
    right = np.append(np.zeros((A_lu.shape[0], buffer_pixed_amount)), A_lu[:, :-buffer_pixed_amount], axis=1)
    shifted_right = A_lu + right
    shifted_right = shifted_right != 0
    return shifted_right


def superpose_up(A_lu, buffer_pixed_amount):
    # shift the matrix to up
    # shift amount is defined by buffer_pixel_amount
    up = np.append(A_lu[buffer_pixed_amount:, :], np.zeros((buffer_pixed_amount, A_lu.shape[1])), axis=0)
    shifted_up = A_lu + up
    shifted_up = shifted_up != 0
    return shifted_up


def superpose_down(A_lu, buffer_pixed_amount):
    # shift the matrix to down
    # shift amount is defined by buffer_pixel_amount
    down = np.append(np.zeros((buffer_pixed_amount, A_lu.shape[1])), A_lu[:-buffer_pixed_amount, :], axis=0)
    shifted_down = A_lu + down
    shifted_down = shifted_down != 0
    return shifted_down


def field_exists(field_name, shp_path):
    shp = ogr.Open(shp_path, 0)
    lyr = shp.GetLayer()
    lyr_dfn = lyr.GetLayerDefn()

    exists = False
    for i in range(lyr_dfn.GetFieldCount()):
        exists = exists or (field_name == lyr_dfn.GetFieldDefn(i).GetName())
    return exists


def changeExt2tif(filepath):

    base = os.path.splitext(filepath)[0]

    return base + '.tif'


def calc_gcr(Crd, res, GCR):
    # This code creates a GCR wieghing matrix for the deisred geographic extent. The sizing of the PV system is
    # conducted on Dec 22 for a shade-free exposure to the Sun during a given number of hours.
    # INPUTS:
    # north_, east_, south_, west_: desired geographic extent
    # res: resolution of MERRA data & desired resolution in lat/lon
    # Shadefree_period: duration of the shade-free period

    # Initialisation
    Ind = ind_merra_high(Crd, res)  # Range indices for high resolution matrices, superposed to MERRA data
    m = Ind[:, 0] - Ind[:, 2] + 1  # number of rows for the high resolution matrix over MERRA
    n = Ind[:, 1] - Ind[:, 3] + 1  # number of cols for the high resolution matrix over MERRA

    # Vector of latitudes between (south) and (north), with resolution (res_should) degrees
    lat = np.arange((Crd[0, 2] + res[1, 0] / 2), (Crd[0, 0] - res[1, 0] / 2), res[1, 0])[np.newaxis].T

    # Solar time where shade-free exposure starts
    omegast = 12 - GCR["shadefree_period"] / 2

    # Calculation
    omega = 15 * (omegast - 12)  # Hour angle
    phi = abs(lat)  # Latitude angle

    beta = np.maximum(phi, 15)  # Tilt angle = latitude, but at least 15 degrees

    if Crd[0, 2] > 0:
        day = GCR["day_north"]
        # Declination angle
        delta = repmat(arcsind(0.3978) * sin(
            day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), int(m), 1)

    if Crd[0, 0] < 0:
        day = GCR["day_south"]
        # Declination angle
        delta = repmat(arcsind(0.3978) * sin(
            day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), int(m), 1)

    if (Crd[0, 2] * Crd[0, 0]) < 0:
        lat_pos = np.sum((lat > 0).astype(int))
        day = GCR["day_north"]
        # Declination angle
        delta_pos = repmat(arcsind(0.3978) * sin(
            day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), lat_pos, 1)

        lat_neg = np.sum((lat < 0).astype(int))
        day = GCR["day_south"]
        # Declination angle
        delta_neg = repmat(arcsind(0.3978) * sin(
            day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), lat_neg, 1)
        delta = np.append(delta_neg, delta_pos, axis=0)

    # Elevation angle
    alpha = arcsind(sind(delta) * sind(phi) + cosd(delta) * cosd(phi) * cosd(omega))

    # Azimuth angle
    azi = arccosd((sind(delta) * cosd(phi) - cosd(delta) * sind(phi) * cosd(omega)) / cosd(alpha))

    # The GCR applies for each line, independently from the longitude
    A_GCR = repmat((1 / (cosd(beta) + np.abs(cosd(azi)) * sind(beta) / tand(alpha))), 1, int(n))

    # Fix too large and too small values of GCR
    A_GCR[A_GCR < 0.2] = 0.2
    A_GCR[A_GCR > 0.9] = 0.9

    return A_GCR


