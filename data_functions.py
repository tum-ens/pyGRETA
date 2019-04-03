import os
from util import *
import math as m
import numpy as np
import rasterio
from rasterio import windows, mask, MemoryFile


def calc_ext(regb, ext, res):
    minRow = m.floor(regb["miny"] / res[1, 0]) * res[1, 0]
    maxRow = m.ceil(regb["maxy"] / res[1, 0]) * res[1, 0]
    minCol = m.floor(regb["minx"] / res[1, 1]) * res[1, 1]
    maxCol = m.ceil(regb["maxx"] / res[1, 1]) * res[1, 1]

    return [[min(m.ceil((ext[0, 0] - res[0, 0] / 2) / res[0, 0]) * res[0, 0] + res[0, 0] / 2, maxRow),
             min(m.ceil((ext[0, 1] - res[0, 1] / 2) / res[0, 1]) * res[0, 1] + res[0, 1] / 2, maxCol),
             max(m.ceil((ext[0, 2] - res[0, 0] / 2) / res[0, 0]) * res[0, 0] + res[0, 0] / 2, minRow),
             max(m.ceil((ext[0, 3] - res[0, 1] / 2) / res[0, 1]) * res[0, 1] + res[0, 1] / 2, minCol)]]


def crd_merra(Crd_regions, res_low):
    ''' description '''
    Crd = np.array([(np.ceil((Crd_regions[:, 0] - res_low[0] / 2) / res_low[0]) * res_low[0] + res_low[0] / 2),
                    (np.ceil((Crd_regions[:, 1] - res_low[1] / 2) / res_low[1]) * res_low[1] + res_low[1] / 2),
                    (np.floor((Crd_regions[:, 2] + res_low[0] / 2) / res_low[0]) * res_low[0] - res_low[0] / 2),
                    (np.floor((Crd_regions[:, 3] + res_low[1] / 2) / res_low[1]) * res_low[1] - res_low[1] / 2)])
    Crd = Crd.T
    return Crd


def crd_exact_box(Ind, Crd_all, res_high):
    Ind = Ind[np.newaxis]

    Crd = [Ind[:, 0] * res_high[0] + Crd_all[2],
           Ind[:, 1] * res_high[1] + Crd_all[3],
           (Ind[:, 2] - 1) * res_high[0] + Crd_all[2],
           (Ind[:, 3] - 1) * res_high[1] + Crd_all[3]]
    return Crd

  
def crd_exact_points(Ind_points, Crd_all, res):
    ''' description 
    Ind_points: tuple of indices in the vertical and horizontal axes. '''

    Crd_points = [Ind_points[0] * res[0] + Crd_all[2],
                  Ind_points[1] * res[1] + Crd_all[3]]
    return Crd_points


def ind_merra(Crd, Crd_all, res):
    ''' description '''
    if len(Crd.shape) == 1:
        Crd = Crd[np.newaxis]
    Ind = np.array([(Crd[:, 0] - Crd_all[2]) / res[0],
                    (Crd[:, 1] - Crd_all[3]) / res[1],
                    (Crd[:, 2] - Crd_all[2]) / res[0] + 1,
                    (Crd[:, 3] - Crd_all[3]) / res[1] + 1])
    Ind = np.transpose(Ind.astype(int))
    return Ind


def ind_global(Crd, res_high):
    ''' description '''
    if len(Crd.shape) == 1:
        Crd = Crd[np.newaxis]
    Ind = np.array([np.round((90 - Crd[:, 0]) / res_high[0]) + 1,
                    np.round((180 + Crd[:, 1]) / res_high[1]),
                    np.round((90 - Crd[:, 2]) / res_high[0]),
                    np.round((180 + Crd[:, 3]) / res_high[1]) + 1])
    Ind = np.transpose(Ind.astype(int))
    return Ind


def calc_geotiff(Crd_all, res_high):
    """
    Returns dictionary containing the Georefferencing parameters for geotiff creation,
    based on the desired extent and resolution

    :param Crd: Extent
    :param res: resolution
    """
    GeoRef = {"RasterOrigin": [Crd_all[3], Crd_all[0]],
              "RasterOrigin_alt": [Crd_all[3], Crd_all[2]],
              "pixelWidth": res_high[1],
              "pixelheight": -res_high[0]}
    return GeoRef


def calc_region(region, Crd_reg, res_high, GeoRef):
    ''' description - why is there a minus sign?'''
    latlim = Crd_reg[2] - Crd_reg[0]
    lonlim = Crd_reg[3] - Crd_reg[1]
    M = int(m.fabs(latlim) / res_high[0])
    N = int(m.fabs(lonlim) / res_high[1])
    A_region = np.ones((M, N))
    origin = [Crd_reg[3], Crd_reg[2]]

    if region.geometry.geom_type == 'MultiPolygon':
        features = [feature for feature in region.geometry]
    else:
        features = [region.geometry]
    west = origin[0]
    north = origin[1]
    profile = {'driver': 'GTiff',
               'height': M,
               'width': N,
               'count': 1,
               'dtype': rasterio.float64,
               'crs': 'EPSG:4326',
               'transform': rasterio.transform.from_origin(west, north, GeoRef["pixelWidth"], GeoRef["pixelHeight"])}

    with MemoryFile() as memfile:
        with memfile.open(**profile) as f:
            f.write(A_region, 1)
            out_image, out_transform = mask.mask(f, features, crop=False, nodata=0, all_touched=False, filled=True)
        A_region = out_image[0]

    return A_region


def calc_gcr(Crd, m, n, res, GCR):
    """
    This function creates a GCR weighting matrix for the desired geographic extent
    The sizing of the PV system is conducted on DEC 22 for a shade-free exposure
    to the sun during a given number of hours

    :param Crd: desired geographic extent (north, east, south, west)
    :param m, n, res: desired resolution in lat and lon, and resolution of MERRA data
    :param GCR: Duration of the shade-free period
    """

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


def calc_areas(Crd, n, res, reg):
    # WSG84 ellipsoid constants
    a = 6378137  # major axis
    b = 6356752.3142  # minor axis
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


def create_buffer(A_lu, buffer_pixel_amount):
    """
    This function creates a buffer around urban areas, based on a Von Neumann neighborhood.
    A_lu matrix element values range from 0 to 16:
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

    :param A_lu: Landuse matrix
    :param buffer_pixel_amount: Buffer amount
    """

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
    """
    Used as part of create_buffer()
    Shift and superpose to the left, shift amount is defined by buffer_pixel amount
    """

    left = np.append(A_lu[:, buffer_pixed_amount:], np.zeros((A_lu.shape[0], buffer_pixed_amount)), axis=1)
    shifted_left = A_lu + left
    shifted_left = shifted_left != 0
    return shifted_left


def superpose_right(A_lu, buffer_pixed_amount):
    """
    Used as part of create_buffer()
    Shift and superpose to the right, shift amount is defined by buffer_pixel amount
    """

    right = np.append(np.zeros((A_lu.shape[0], buffer_pixed_amount)), A_lu[:, :-buffer_pixed_amount], axis=1)
    shifted_right = A_lu + right
    shifted_right = shifted_right != 0
    return shifted_right


def superpose_up(A_lu, buffer_pixed_amount):
    """
    Used as part of create_buffer()
    Shift and superpose upward, shift amount is defined by buffer_pixel amount
    """

    up = np.append(A_lu[buffer_pixed_amount:, :], np.zeros((buffer_pixed_amount, A_lu.shape[1])), axis=0)
    shifted_up = A_lu + up
    shifted_up = shifted_up != 0
    return shifted_up


def superpose_down(A_lu, buffer_pixed_amount):
    """
    Used as part of create_buffer()
    Shift and superpose to the downward, shift amount is defined by buffer_pixel amount
    """

    down = np.append(np.zeros((buffer_pixed_amount, A_lu.shape[1])), A_lu[:-buffer_pixed_amount, :], axis=0)
    shifted_down = A_lu + down
    shifted_down = shifted_down != 0
    return shifted_down
