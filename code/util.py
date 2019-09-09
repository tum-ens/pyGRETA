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
import sys
import math
import rasterio
from rasterio import windows, mask, MemoryFile
import pandas as pd
import numpy as np
from scipy.ndimage import generic_filter, convolve
import geopandas as gpd
from shapely.geometry import mapping, Point, Polygon
import fiona
import hdf5storage
from multiprocessing import Pool
from itertools import product
import h5netcdf
import cProfile
import pstats
import shutil
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import json

def sind(alpha):
    """
    Calculates the sine of an angle in degrees.
    
    :param alpha: Angle in degrees.
    :type alpha: float
    :return: The sine of the angle.
    :rtype: float
    """
    return sin(np.deg2rad(alpha))


def cosd(alpha):
    """
    Calculates the cosine of an angle in degrees.
    
    :param alpha: Angle in degrees.
    :type alpha: float
    :return: The cosine of the angle.
    :rtype: float
    """
    return cos(np.deg2rad(alpha))


def tand(alpha):
    """
    Calculates the tangent of an angle in degrees.
    
    :param alpha: Angle in degrees.
    :type alpha: float
    :return: The tangent of the angle.
    :rtype: float
    """
    return tan(np.deg2rad(alpha))


def arcsind(digit):
    """
    Calculates the inverse sine of a number.
    
    :param digit: Number between -1 and 1.
    :type digit: float
    :return: The inverse sine of the number in degrees.
    :rtype: float
    """
    return np.rad2deg(arcsin(digit))


def arccosd(digit):
    """
    Calculates the inverse cosine of a number.
    
    :param digit: Number between -1 and 1.
    :type digit: float
    :return: The inverse cosine of the number in degrees.
    :rtype: float
    """
    return np.rad2deg(arccos(digit))


def arctand(digit):
    """
    Calculates the inverse tangent of a number.
    
    :param digit: Number.
    :type digit: float
    :return: The inverse tangent of the number in degrees.
    :rtype: float
    """
    return np.rad2deg(arctan(digit))


def hourofmonth():
    """
    Calculates the rank within a year of the first hour of each month.
    
    :return: The rank of the first hour of each month.
    :rtype: list
    """
    h = 24 * np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]).astype(int)
    for i in range(1, 12):
        h[i] = h[i] + h[i - 1]
    return h.astype(int)


def intersection(lst1, lst2):
    """
    Calculates the intersection between two lists.
    
    :param lst1: First list of elements.
    :type lst1: list
    :param lst2: Second list of elements.
    :type lst2: list
    :return: The unique elements that exist in both lists, without repetition.
    :rtype: list
    """
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3


def resizem(A_in, row_new, col_new):
    """
    Resizes regular data grid, by copying and pasting parts of the original array.

    :param A_in: Input matrix.
    :type A_in: numpy array
    :param row_new: New number of rows.
    :type row_new: integer
    :param col_new: New number of columns.
    :type col_new: integer
    :return: Resized matrix.
    :rtype: numpy array
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
    Saves array to geotiff raster format based on EPSG 4326.

    :param newRasterfn: Output path of the raster.
    :type newRasterfn: string
    :param rasterOrigin: Latitude and longitude of the Northwestern corner of the raster.
    :type rasterOrigin: list of two floats
    :param pixelWidth:  Pixel width (might be negative).
    :type pixelWidth: integer
    :param pixelHeight: Pixel height (might be negative).
    :type pixelHeight: integer
    :param array: Array to be converted into a raster.
    :type array: numpy array
    :return: The raster file will be saved in the desired path *newRasterfn*.
    :rtype: None
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
    """
    Creates a generator to iterate between the characters *c1* and *c2*, including the latter.
    
    :param c1: First character in the iteration.
    :type c1: char
    :param c2: Last character in the iteration (included).
    :type c2: char
    :return: Generator to iterate between the characters *c1* and *c2*.
    :rtype: generator ?????
    """
    for c in range(ord(c1), ord(c2) + 1):
        yield chr(c)


def changem(A, newval, oldval):
    """
    This function replaces existing values *oldval* in a data array *A* by new values *newval*.
    
    *oldval* and *newval* must have the same size.

    :param A: Input matrix.
    :type A: numpy array
    :param newval: Vector of new values to be set.
    :type newval: numpy array
    :param oldval: Vector of old values to be replaced.
    :param oldval: numpy array
    :return: The updated array.
    :rtype: numpy array
    """
    Out = np.zeros(A.shape)
    z = np.array((oldval, newval)).T
    for i, j in z:
        np.place(Out, A == i, j)
    return Out


def list_files(directory, format):
    """
    ?????
    """
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
    """
    Is called at every process start to set its priority
    :return:
    """
    check = check[0]
    p = psutil.Process(os.getpid())
    if check:
        if sys.platform.startswith('win'):
            # Windows priority
            p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        elif sys.platform.startswith('linux'):
            # Linux priority
            p.nice(1)
    else:
        if sys.platform.startswith('win'):
            # Windows priority
            p.nice(psutil.NORMAL_PRIORITY_CLASS)
        elif sys.platform.startswith('linux'):
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

    print("\nNo Regression Folder found, It has been created under: \n" + paths["regression_in"])
    print("\nPlease populate the input folder with the required data. \n"
          "(look at README.txt in the folder for more information)\n")


def reg_miss_files():

    print("\nPlease populate the input folder with the technology corresponding data.\n"
          "(look at README.txt in the folder for more information)\n")


def display_progress(message, progress_stat):
    length = progress_stat[0]
    status = progress_stat[1]
    sys.stdout.write('\r')
    sys.stdout.write(message + ' ' + '[%-50s] %d%%' % ('=' * ((status * 50) // length), (status * 100) // length))
    sys.stdout.flush()
    if status == length:
        print('\n')
