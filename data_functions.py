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


def ind_merra(Crd, res):
    ind = np.array([(Crd[:, 0] - Crd[-1, 2]) / res[0],
                    (Crd[:, 1] - Crd[-1, 3]) / res[1],
                    (Crd[:, 2] - Crd[-1, 2]) / res[0] + 1,
                    (Crd[:, 3] - Crd[-1, 3]) / res[1] + 1])
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
    GeoRef = {"RasterOrigin": [Crd[-1, 3], Crd[-1, 0]],
              "RasterOrigin_alt": [Crd[-1, 3], Crd[-1, 2]],
              "pixelWidth": res[1, 1],
              "pixelHeight": -res[1, 0]}
    return GeoRef


def calc_region(region, Crd, res, GeoRef):

    latlim = Crd[2] - Crd[0]
    lonlim = Crd[3] - Crd[1]
    M = int(m.fabs(latlim) / res[1, 0])
    N = int(m.fabs(lonlim) / res[1, 1])
    A_region = np.ones((M, N))
    origin = [Crd[3], Crd[2]]

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


