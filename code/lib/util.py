from osgeo import ogr, gdal, osr
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
import shutil
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import json
from warnings import warn


def sind(alpha):
    """
    This function calculates the sine of an angle in degrees.
    
    :param alpha: Angle in degrees.
    :type alpha: float

    :return: The sine of the angle.
    :rtype: float
    """
    return sin(np.deg2rad(alpha))


def cosd(alpha):
    """
    This function calculates the cosine of an angle in degrees.
    
    :param alpha: Angle in degrees.
    :type alpha: float

    :return: The cosine of the angle.
    :rtype: float
    """
    return cos(np.deg2rad(alpha))


def tand(alpha):
    """
    This function calculates the tangent of an angle in degrees.
    
    :param alpha: Angle in degrees.
    :type alpha: float

    :return: The tangent of the angle.
    :rtype: float
    """
    return tan(np.deg2rad(alpha))


def arcsind(digit):
    """
    This function calculates the inverse sine of a number.
    
    :param digit: Number between -1 and 1.
    :type digit: float

    :return: The inverse sine of the number in degrees.
    :rtype: float
    """
    return np.rad2deg(arcsin(digit))


def arccosd(digit):
    """
    This function calculates the inverse cosine of a number.
    
    :param digit: Number between -1 and 1.
    :type digit: float

    :return: The inverse cosine of the number in degrees.
    :rtype: float
    """
    return np.rad2deg(arccos(digit))


def arctand(digit):
    """
    This function calculates the inverse tangent of a number.
    
    :param digit: Number.
    :type digit: float

    :return: The inverse tangent of the number in degrees.
    :rtype: float
    """
    return np.rad2deg(arctan(digit))


def hourofmonth():
    """
    This function calculates the rank within a year of the first hour of each month.
    
    :return: The rank of the first hour of each month.
    :rtype: list
    """
    h = 24 * np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]).astype(int)
    for i in range(1, 12):
        h[i] = h[i] + h[i - 1]
    return h.astype(int)


def intersection(lst1, lst2):
    """
    This function calculates the intersection between two lists.
    
    :param lst1: First list of elements.
    :type lst1: list
    :param lst2: Second list of elements.
    :type lst2: list

    :return lst3: The unique elements that exist in both lists, without repetition.
    :rtype: list
    """
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3


def resizem(A_in, row_new, col_new):
    """
    This function resizes regular data grid, by copying and pasting parts of the original array.

    :param A_in: Input matrix.
    :type A_in: numpy array
    :param row_new: New number of rows.
    :type row_new: integer
    :param col_new: New number of columns.
    :type col_new: integer

    :return A_out: Resized matrix.
    :rtype: numpy array
    """
    row_rep = row_new // np.shape(A_in)[0]
    col_rep = col_new // np.shape(A_in)[1]
    A_inf = A_in.flatten(order="F")[np.newaxis]
    A_out = reshape(
        repmat(
            reshape(reshape(repmat((A_in.flatten(order="F")[np.newaxis]), row_rep, 1), (row_new, -1), order="F").T, (-1, 1), order="F"), 1, col_rep
        ).T,
        (col_new, row_new),
        order="F",
    ).T

    return A_out


def char_range(c1, c2):
    """
    This function creates a generator to iterate between the characters *c1* and *c2*, including the latter.
    
    :param c1: First character in the iteration.
    :type c1: char
    :param c2: Last character in the iteration (included).
    :type c2: char

    :return: Generator to iterate between the characters *c1* and *c2*.
    :rtype: python generator
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
    :type oldval: numpy array

    :return Out: The updated array.
    :rtype: numpy array
    """
    Out = np.zeros(A.shape)
    z = np.array((oldval, newval)).T
    for i, j in z:
        np.place(Out, A == i, j)
    return Out


def ind2sub(array_shape, ind):
    """
    This function converts linear indices to subscripts.
    
    :param array_shape: Dimensions of the array (# of rows, # of columns).
    :type array_shape: tuple (int, int)
    :param ind: Linear index.
    :type index: int

    :return: Tuple of indices in each dimension (row index, column index).
    :rtype: tuple(int, int)
    """
    return np.unravel_index(ind, array_shape, order="F")


def field_exists(field_name, shp_path):
    """
    This function returns whether the specified field exists or not in the shapefile linked by a path.

    :param field_name: Name of the field to be checked for.
    :type field_name: str
    :param shp_path: Path to the shapefile.
    :type shp_path: str

    :return: ``True`` if it exists or ``False`` if it doesn't exist.
    :rtype: bool
    """
    shp = ogr.Open(shp_path, 0)
    lyr = shp.GetLayer()
    lyr_dfn = lyr.GetLayerDefn()

    exists = False
    for i in range(lyr_dfn.GetFieldCount()):
        exists = exists or (field_name == lyr_dfn.GetFieldDefn(i).GetName())
    return exists


def changeExt2tif(filepath):
    """
    This function changes the extension of a file path to .tif.

    :param filepath: Path to the file.
    :type filepath: str

    :return: New path with .tif as extension.
    :rtype: str
    """
    base = os.path.splitext(filepath)[0]
    return base + ".tif"


def sumnorm_MERRA2(A, m, n, res_low, res_desired):
    """
    This function calculates the average of high resolution data if it is aggregated into a lower resolution.

    :param A: High-resolution data.
    :type A: numpy array
    :param m: Number of rows in the low resolution.
    :type m: int
    :param n: Number of columns in the low resolution.
    :type n: int
    :param res_low: Numpy array with with two numbers. The first number is the resolution in the vertical dimension (in degrees of latitude),
        the second is for the horizontal dimension (in degrees of longitude).
    :type res_low: numpy array
    :param res_desired: Numpy array with with two numbers. The first number is the resolution in the vertical dimension (in degrees of latitude),
        the second is for the horizontal dimension (in degrees of longitude).
    :type res_desired: numpy array

    :return s: Aggregated average of *A* on the low resolution.
    :rtype: numpy array
    """
    s = np.zeros((m, n))
    row_step = int(res_low[0] / res_desired[0])
    col_step = int(res_low[1] / res_desired[1])
    for i in range(0, m):
        for j in range(0, n):
            s[i, j] = np.sum(A[(row_step * i) : (row_step * (i + 1)), (col_step * j) : (col_step * (j + 1))]) / (row_step * col_step)
    return s


def limit_cpu(check):
    """
    This functions sets the priority of a process for CPU time and RAM allocation at two levels: average or below average.

    :param check: If ``True``, the process is set a below average priority rating allowing other programs to run undisturbed.
        if ``False``, the process is given the same priority as all other user processes currently running on the machine,
        leading to faster calculation times.
    :type check: boolean

    :return: The priority of the process is set.
    :rtype: None
    """
    check = check[0]
    p = psutil.Process(os.getpid())
    if check:
        if sys.platform.startswith("win"):
            # Windows priority
            p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        elif sys.platform.startswith("linux"):
            # Linux priority
            p.nice(1)
    else:
        if sys.platform.startswith("win"):
            # Windows priority
            p.nice(psutil.NORMAL_PRIORITY_CLASS)
        elif sys.platform.startswith("linux"):
            # Linux priority
            p.nice(0)


def timecheck(*args):
    """
    This function prints information about the progress of the script by displaying the function currently running, and optionally
    an input message, with a corresponding timestamp. If more than one argument is passed to the function, it will raise an exception.

    :param args: Message to be displayed with the function name and the timestamp (optional).
    :type args: string

    :return: The time stamp is printed.
    :rtype: None
    :raise: Too many arguments have been passed to the function, the maximum is only one string.
    """
    if len(args) == 0:
        print(inspect.stack()[1].function + str(datetime.datetime.now().strftime(": %H:%M:%S:%f")) + "\n")

    elif len(args) == 1:
        print(inspect.stack()[1].function + " - " + str(args[0]) + str(datetime.datetime.now().strftime(": %H:%M:%S:%f")) + "\n")

    else:
        raise Exception("Too many arguments have been passed.\nExpected: zero or one \nPassed: " + format(len(args)))


def display_progress(message, progress_stat):
    """
    This function displays a progress bar for long computations. To be used as part of a loop or with multiprocessing.

    :param message: Message to be displayed with the progress bar.
    :type message: string
    :param progress_stat: Tuple containing the total length of the calculation and the current status or progress.
    :type progress_stat: tuple(int, int)

    :return: The status bar is printed.
    :rtype: None
    """
    length = progress_stat[0]
    status = progress_stat[1]
    sys.stdout.write("\r")
    sys.stdout.write(message + " " + "[%-50s] %d%%" % ("=" * ((status * 50) // length), (status * 100) // length))
    sys.stdout.flush()
    if status == length:
        print("\n")


def create_json(filepath, param, param_keys, paths, paths_keys):
    """
    Creates a metadata JSON file containing information about the file in filepath by storing the relevant keys from
    both the param and path dictionaries.

    :param filepath: Path to the file for which the JSON file will be created.
    :type filepath: string
    :param param: Dictionary of dictionaries containing the user input parameters and intermediate outputs.
    :type param: dict
    :param param_keys: Keys of the parameters to be extracted from the *param* dictionary and saved into the JSON file.
    :type param_keys: list of strings
    :param paths: Dictionary of dictionaries containing the paths for all files.
    :type paths: dict
    :param paths_keys: Keys of the paths to be extracted from the *paths* dictionary and saved into the JSON file.
    :type paths_keys: list of strings

    :return: The JSON file will be saved in the desired path *filepath*.
    :rtype: None
    """
    new_file = os.path.splitext(filepath)[0] + ".json"
    new_dict = {}
    # Add standard keys
    param_keys = param_keys + ["author", "comment"]
    for key in param_keys:
        new_dict[key] = param[key]
        if type(param[key]) == np.ndarray:
            new_dict[key] = param[key].tolist()
        if type(param[key]) == tuple:
            param[key] = list(param[key])
            c = 0
            for e in param[key]:
                if type(e) == np.ndarray:
                    new_dict[key][c] = e.tolist()
                c += 1
        if type(param[key]) == dict:
            for k, v in param[key].items():
                if type(v) == np.ndarray:
                    new_dict[key][k] = v.tolist()
                if type(v) == tuple:
                    param[key][k] = list(param[key][k])
                    c = 0
                    for e in param[key][k]:
                        if type(e) == np.ndarray:
                            new_dict[key][k][c] = e.tolist()
                        c += 1
                if type(v) == dict:
                    for k2, v2 in v.items():
                        if type(v2) == np.ndarray:
                            new_dict[key][k][k2] = v2.tolist()
                        if type(v2) == tuple:
                            param[key][k][k2] = list(param[key][k][k2])
                            c = 0
                            for e in param[key][k][k2]:
                                if type(e) == np.ndarray:
                                    new_dict[key][k][k2][c] = e.tolist()
                                c += 1
                        if type(v2) == dict:
                            for k3, v3 in v.items():
                                if type(v3) == np.ndarray:
                                    new_dict[key][k][k2][k3] = v3.tolist()
                                if type(v3) == tuple:
                                    param[key][k][k2][k3] = list(param[key][k][k2][k3])
                                    c = 0
                                    for e in param[key][k][k2][k3]:
                                        if type(e) == np.ndarray:
                                            new_dict[key][k][k2][k3][c] = e.tolist()
                                        c += 1

    for key in paths_keys:
        new_dict[key] = paths[key]
    # Add timestamp
    new_dict["timestamp"] = str(datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
    # Add caller function's name
    new_dict["function"] = inspect.stack()[1][3]
    with open(new_file, "w") as json_file:
        json.dump(new_dict, json_file)
    print("files saved: " + new_file)
