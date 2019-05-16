import numpy as np
import gdal
import osr
from osgeo import ogr
from numpy.matlib import repmat, reshape, sin, arcsin, cos, arccos, tan, arctan
import os
from os import getcwd, chdir
from glob import glob
import psutil
import datetime
import inspect


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


def hourofmonth():
    h = 24 * np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]).astype(int)
    for i in range(1, 12):
        h[i] = h[i] + h[i - 1]
    return h.astype(int)


def intersection(lst1, lst2):
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3


def resizem(A_in, row_new, col_new):
    """
    Resize regular data grid, by copying and pasting parts of the original array

    :param A_in: Input matrix
    :param row_new: new number of rows
    :param col_new: new number of columns
    :return: resized matrix
    """
    row_rep = row_new // np.shape(A_in)[0]
    col_rep = col_new // np.shape(A_in)[1]
    A_inf = (A_in.flatten(order='F')[np.newaxis])
    A_out = reshape(repmat(
        reshape(reshape(repmat((A_in.flatten(order='F')[np.newaxis]), row_rep, 1), (row_new, -1), order='F').T, (-1, 1),
                order='F'), 1, col_rep).T, (col_new, row_new), order='F').T

    return A_out


def array2raster(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array):
    """
    Save array to geotiff raster format based on EPSG 4326

    :param newRasterfn: desired paths
    :param rasterOrigin: latitude and longitube of the top left corner of the raster
    :param pixelWidth:  pixel width (might be negative)
    :param pixelHeight: pixel height(might be negative)
    :param array: Raster array
    :return:
    """
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
    """
    This function replaces values in a data array.

    :param A: Input Matrix
    :param newval: New values to be set
    :param oldval: Old values to be replaced
    """
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


def sub2ind(array_shape, rows, cols):
    """
    Convert subscripts to linear indices

    :param array_shape: tuple (# of rows, # of columns)
    :param rows: row values
    :param cols: column values
    :return: Index
    """
    return np.ravel_multi_index((rows, cols), array_shape, order='F')


def ind2sub(array_shape, ind):
    """
    Convert linear indices to subscripts
    :param array_shape: tuple (# of rows, # of columns)
    :param ind: Index
    :return: tuple (row values, column values)
    """
    return np.unravel_index(ind, array_shape, order='F')


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


def sumnorm_MERRA2(A, m, n, res_low, res_desired):
    s = np.zeros((m, n))
    row_step = int(res_low[0] / res_desired[0])
    col_step = int(res_low[1] / res_desired[1])
    for i in range(0, m):
        for j in range(0, n):
            s[i, j] = np.sum(A[(row_step * i):(row_step * (i + 1)),
                             (col_step * j):(col_step * (j + 1))]) / (row_step * col_step)
    return s


def limit_cpu(check):
    from sys import platform
    """
    Is called at every process start to set its priority
    :return:
    """
    check = check[0]
    p = psutil.Process(os.getpid())
    if check:
        if platform.startswith('win'):
            # Windows priority
            p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        elif platform.startswith('linux'):
            # Linux priority
            p.nice(1)
    else:
        if platform.startswith('win'):
            # Windows priority
            p.nice(psutil.NORMAL_PRIORITY_CLASS)
        elif platform.startswith('linux'):
            # Linux priority
            p.nice(0)


def timecheck(*args):

    if len(args) == 0:
        print(inspect.stack()[1].function + str(datetime.datetime.now().strftime(": %H:%M:%S:%f")))

    elif len(args) == 1:
        print(inspect.stack()[1].function + ' - ' + str(args[0])
              + str(datetime.datetime.now().strftime(": %H:%M:%S:%f")))

    else:
        raise Exception('Too many arguments have been passed.\nExpected: zero or one \nPassed: ' + format(len(args)))


def reg_miss_folder(paths):

    print("\nNo Regression Folder found, It has been created under: \n" + paths["regression"])
    print("\nPlease populate the input folder with the required data. \n"
          "(look at README.txt in the folder for more information)\n")


def reg_miss_files():

    print("\nPlease populate the input folder with the technology corresponding data.\n"
          "(look at README.txt in the folder for more information)\n")
