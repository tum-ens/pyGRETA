import os
import sys
# os.environ['MKL_NUM_THREADS'] = '8'
from data_functions import *
from util import *
import numpy as np

np.seterr(divide='ignore')  # Repress Invalid value or division by zero error
from numpy.matlib import repmat, sin, cos
import hdf5storage


def calc_CF_solar(hour, reg_ind, param, merraData, rasterData):
    pv = param["PV"]["technical"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    res_high = param["res_high"]
    Crd_all = param["Crd_all"]

    # Load MERRA data, increase its resolution, and fit it to the extent
    CLEARNESS_h = merraData["CLEARNESS"][:, :, hour]
    CLEARNESS_h = resizem(CLEARNESS_h, m_high, n_high)
    reg_ind_h = np.nonzero(CLEARNESS_h)
    # Filter out night hours for every valid point
    filter_lat = np.logical_and(reg_ind[0] >= reg_ind_h[0].min(), reg_ind[0] <= reg_ind_h[0].max())
    filter_lon = np.logical_and(reg_ind[1] >= reg_ind_h[1].min(), reg_ind[1] <= reg_ind_h[1].max())
    filter = np.logical_and(filter_lat, filter_lon)
    reg_ind_h = (reg_ind[0][filter], reg_ind[1][filter])
    if len(reg_ind_h[0]) == 0:
        CF_pv = np.zeros(reg_ind[0].shape)
        CF_csp = np.zeros(reg_ind[0].shape)
        return CF_pv, CF_csp

    CLEARNESS_h = CLEARNESS_h[reg_ind_h]
    TEMP_h = merraData["T2M"][:, :, hour]
    TEMP_h = resizem(TEMP_h, m_high, n_high) - 273.15  # Convert to Celsius
    TEMP_h = TEMP_h[reg_ind_h]
    # Compute the angles
    A_phi, A_omega, A_delta, A_alpha, A_beta, A_azimuth, A_orientation, sunrise, sunset = \
        angles(hour, reg_ind_h, Crd_all, res_high)
    # Other matrices
    A_albedo = rasterData["A_albedo"][reg_ind_h]
    A_Ross = rasterData["A_Ross"][reg_ind_h]

    if pv["tracking"] == 1:
        A_beta = 90 - A_alpha
    elif pv["tracking"] == 2:
        A_beta = 90 - A_alpha
        A_orientation = A_azimuth

    aux = np.maximum(np.minimum((sind(A_delta) * sind(A_phi) * cosd(A_beta)
                                 - sind(A_delta) * cosd(A_phi) * sind(A_beta) * cosd(A_orientation)
                                 + cosd(A_delta) * cosd(A_phi) * cosd(A_beta) * cosd(A_omega)
                                 + cosd(A_delta) * sind(A_phi) * sind(A_beta) * cosd(A_orientation) * cosd(A_omega)
                                 + cosd(A_delta) * sind(A_beta) * sind(A_orientation) * sind(A_omega)),
                                1), -1)
    A_incidence = arccosd(aux)

    # Compute the hourly TOA radiation
    TOA_h = toa_hourly(A_alpha, hour)

    # Compute the ratio of diffuse radiation
    RATIO = global2diff(CLEARNESS_h, A_alpha.shape)

    # Compute the coefficients for the HDKR model
    R_b = cosd(A_incidence) / sind(A_alpha)
    R_b[R_b <= 0] = 0
    A_i = (1 - RATIO) * CLEARNESS_h
    f = (1 - RATIO) ** 0.5

    F_direct, F_diffuse, F_reflected = coefficients(A_beta, RATIO, R_b, A_i, f, hour, sunrise, sunset)

    # Compute the shading losses
    # The following line requires long computation time, replace with SHADING = 0 for fast computation
    SHADING = 0

    F = F_diffuse + F_direct * (1 - SHADING) + F_reflected * A_albedo
    #import pdb; pdb.set_trace()
    F[F > 1] = 1

    # Compute the incident radiation
    GHI_h = TOA_h * CLEARNESS_h
    GHI_h[np.isnan(GHI_h)] = 0
    G_tilt_h = GHI_h * F

    # Compute losses due to heating of the PV cells
    LOSS_TEMP = loss(G_tilt_h, TEMP_h, A_Ross, pv)

    # Compute the hourly capacity factor
    CF_pv = ((A_alpha > 0) * 1) * G_tilt_h / 1000 * (1 - LOSS_TEMP)

    # For CSP: tracking like pv.tracking = 1
    A_beta = 90 - A_alpha
    aux = np.maximum(np.minimum((sind(A_delta) * sind(A_phi) * cosd(A_beta)
                                 - sind(A_delta) * cosd(A_phi) * sind(A_beta) * cosd(A_orientation)
                                 + cosd(A_delta) * cosd(A_phi) * cosd(A_beta) * cosd(A_omega)
                                 + cosd(A_delta) * sind(A_phi) * sind(A_beta) * cosd(A_orientation) * cosd(A_omega)
                                 + cosd(A_delta) * sind(A_beta) * sind(A_orientation) * sind(A_omega)),
                                1), -1)
    A_incidence = arccosd(aux)
    R_b = cosd(A_incidence) / sind(A_alpha)
    F_direct_csp, _, _ = coefficients(90 - A_alpha, RATIO, R_b, A_i, f, hour, sunrise, sunset)
    CF_csp = TOA_h * CLEARNESS_h * F_direct_csp * (1 - SHADING)
    CF_csp[CF_csp > 1] = 1

    # Adjusting the length of the matrices
    aux = np.zeros(len(reg_ind[0]))
    aux[filter] = CF_pv
    CF_pv = aux
    aux = np.zeros(len(reg_ind[0]))
    aux[filter] = CF_csp
    CF_csp = aux
    return CF_pv, CF_csp


def calc_FLH_solar(hours, args):
    # Decomposing the tuple args
    paths = args[0]
    param = args[1]
    tech = args[2]

    landuse = param["landuse"]
    reg_ind = param["Ind_nz"]

    # Obtain weather matrices
    merraData = {}
    # Clearness index - stored variable CLEARNESS
    merraData["CLEARNESS"] = hdf5storage.read('CLEARNESS', paths["CLEARNESS"])
	
    # Temperature 2m above the ground - stored variable T2M
    merraData["T2M"] = hdf5storage.read('T2M', paths["T2M"])

    rasterData = {}
    # Calculate A matrices
    # A_lu
    with rasterio.open(paths["LU"]) as src:
        w = src.read(1)
    rasterData["A_lu"] = np.flipud(w)
    # rasterData["A_lu"] = rasterData["A_lu"][reg_ind]
    # A_Ross (Temperature coefficients for heating losses)
    rasterData["A_Ross"] = changem(rasterData["A_lu"], param["landuse"]["Ross_coeff"], param["landuse"]["type"]).astype(
        float) / 10000
    # A_albedo (Reflectivity coefficients)
    rasterData["A_albedo"] = changem(rasterData["A_lu"], param["landuse"]["albedo"], param["landuse"]["type"]).astype(
        float) / 100

    FLH = np.zeros(len(reg_ind[0]))
    status = 0
    for hour in hours:
        if hour <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            sys.stdout.write('\r')
            sys.stdout.write(tech + ' ' + param["region"] + ' ' + '[%-50s] %d%%' % (
            '=' * ((status * 50) // len(hours)), (status * 100) // len(hours)))
            sys.stdout.flush()

        if tech == 'PV':
            CF, _ = calc_CF_solar(hour, reg_ind, param, merraData, rasterData)
        elif tech == 'CSP':

            _, CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData)

        # Aggregates CF to obtain the yearly FLH
        CF[np.isnan(CF)] = 0
        FLH = FLH + CF

    return FLH


def calc_TS_solar(hours, args):
    # Decomposing the tuple args
    paths = args[0]
    param = args[1]
    tech = args[2]

    landuse = param["landuse"]
    reg_ind = param[tech]["Ind_points"]

    # Obtain weather matrices
    merraData = {}
    # Clearness index - stored variable CLEARNESS
    merraData["CLEARNESS"] = hdf5storage.read('CLEARNESS', paths["CLEARNESS"])
    # Temperature 2m above the ground - stored variable T2M
    merraData["T2M"] = hdf5storage.read('T2M', paths["T2M"])

    rasterData = {}
    # Calculate A matrices
    # A_lu
    with rasterio.open(paths["LU"]) as src:
        w = src.read(1)
    rasterData["A_lu"] = np.flipud(w)
    # A_Ross (Temperature coefficients for heating losses)
    rasterData["A_Ross"] = changem(rasterData["A_lu"], param["landuse"]["Ross_coeff"], param["landuse"]["type"]).astype(
        float) / 10000
    # A_albedo (Reflectivity coefficients)
    rasterData["A_albedo"] = changem(rasterData["A_lu"], param["landuse"]["albedo"], param["landuse"]["type"]).astype(
        float) / 100

    TS = np.zeros((len(reg_ind[0]), 8760))
    status = 0
    for hour in hours:
        if hour <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            sys.stdout.write('\r')
            sys.stdout.write(tech + ' ' + param["region"] + ' ' + '[%-50s] %d%%' % (
            '=' * ((status * 50) // len(hours)), (status * 100) // len(hours)))
            sys.stdout.flush()

        if tech == 'PV':
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData)[0]
        elif tech == 'CSP':
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData)[1]
        # Aggregates CF to obtain the time series
        CF[np.isnan(CF)] = 0
        TS[:, hour] = CF
    return TS


def angles(hour, reg_ind, Crd_all, res_high):
    # This code creates six matrices for the desired extent, that represent the elevation, tilt, azimuth,
    # and orientation angles in addition to the sunrise and sunset hours of every pixel with the desired resolution.True
    # The inputs are:
    # hour: hour rank in a year (from 1 to 8760)
    # north, est, south, west: desired geographic extent
    # res: resolution of MERRA data & desired resolution in lat/lon
    # Ext_in Extent of the global area of interest

    # Initialization
    Crd_points = crd_exact_points(reg_ind, Crd_all, res_high)
    lat = Crd_points[0]
    lon = Crd_points[1]
    k = len(lat)  # number of points
    N = hour // 24 + 1
    hourofday = hour % 24 + 0.5

    # Calculation
    # Latitude angle
    phi = lat

    # Equation of Time (in hours)
    EOT = -0.128 * sind(N * 360 / 365.25 - 2.80) - 0.165 * sind(2 * N * 360 / 365.25 + 19.7)

    # Time Correction Factor (in hours)
    TC = EOT + lon / 15  # no correction factor for differences to GMT, because data is in GMT

    # Local Solar Time (in hours)
    omegast = hourofday + TC

    # Hour angle (in degrees)
    omega = 15 * (omegast - 12)

    # Declination angle
    delta = np.tile(
        arcsind(0.3978 * sin(N * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(N * 2 * np.pi / 365.25 - 0.0489))), k)
    delta[phi < 0] = - delta[phi < 0]

    # Elevation angle (in degrees)
    alpha = arcsind(sind(delta) * sind(phi) + cosd(delta) * cosd(phi) * cosd(omega))

    # Optimal tilt angle (loosely based on Breyer 2010)
    beta = np.minimum(np.abs(phi), 55)  # The tilt angle is preferably equal to the latitude
    range_lat = np.logical_and(np.abs(phi) >= 35, np.abs(phi) < 65)
    beta[range_lat] = (beta[range_lat] - 35) / 65 * 55 + 35  # Tilt angle does not increase very quickly
    range_lat = np.logical_and(lat >= 35, lat < 65)
    range_lon = np.logical_and(lon >= -20, lon < 30)
    # range_lat = repmat(range_lat, range_lon.shape[1], 1).T
    # range_lon = repmat(range_lon, range_lat.shape[0], 1)
    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat,
                                                                      range_lon)] - 35) / 65 * 45 + 35  # Europe
    range_lat = np.logical_and(lat >= 20, lat < 65)
    range_lon = np.logical_and(lon >= 75, lon < 140)
    # range_lat = repmat(range_lat, range_lon.shape[1], 1).T
    # range_lon = repmat(range_lon, range_lat.shape[0], 1)
    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat,
                                                                      range_lon)] - 20) / 65 * 60 + 20  # Asia/China

    # Azimuth angle (in degrees)
    aux = np.maximum(np.minimum(sind(delta) * cosd(phi) - cosd(delta) * sind(phi) * cosd(omega) / cosd(alpha), 1), -1)
    aziam = arccosd(aux)
    azipm = 360 - aziam
    azi = aziam * ((omega < 0) * 1) + azipm * ((omega >= 0) * 1)

    # Orientation (in degrees)
    orientation = np.zeros(alpha.shape)  # Azimuth of the PV panel is zero for the Southern hemisphere
    orientation[phi < 0] = 180  # Azimuth of the PV panel is 180° for the Southern hemisphere

    # Sunrise and sunset hours in GMT
    # aux = np.maximum(np.minimum((-tand(phi)) * tand(delta), 1), -1)
    aux = np.maximum(np.minimum((-tand(phi)) * tand(delta), 1), -1)
    sunrise = 12 - 1 / 15 * arccosd(aux)
    sunset = 12 + 1 / 15 * arccosd(aux)
    coeff_a = (cosd(phi) / tand(beta)) + sind(phi)
    coeff_b = tand(delta) * (cosd(phi) * cosd(orientation) - sind(phi) / tand(beta))
    aux = np.maximum(np.minimum(((coeff_a * coeff_b - sind(orientation)
                                  * (coeff_a ** 2 - coeff_b ** 2 + sind(orientation) ** 2) ** 0.5
                                  / (coeff_a ** 2 + sind(orientation) ** 2))), 1), -1)
    sunrise_tilt = 12 - 1 / 15 * arccosd(aux)
    sunset_tilt = 12 + 1 / 15 * arccosd(aux)

    #import pdb; pdb.set_trace()
    # sunrise = np.maximum(sunrise, sunrise_tilt) - repmat(TC, m, 1)  # Correction Solar time -> GMT
    # sunset = np.minimum(sunset, sunset_tilt) - repmat(TC, m, 1)  # Correction solar time -> GMT
    sunrise = np.maximum(sunrise, sunrise_tilt) - TC  # Correction Solar time -> GMT
    sunset = np.minimum(sunset, sunset_tilt) - TC  # Correction solar time -> GMT

    # if len(args) == 1:
        # phi = np.diag(phi)
        # omega = np.diag(omega)
        # delta = np.diag(delta)
        # alpha = np.diag(alpha)
        # beta = np.diag(beta)
        # azi = np.diag(azi)
        # orientation = np.diag(orientation)
        # sunrise = np.diag(sunrise)
        # sunset = np.diag(sunset)
    # import pdb; pdb.set_trace()
    return phi, omega, delta, alpha, beta, azi, orientation, sunrise, sunset


def toa_hourly(alpha, *args):
    solarconst = 1367  # in W/m^2

    if len(args) >= 2:
        # Calculate day rank in the year
        month = args[0]
        dayofmonth = args[1]
        DoM = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        N = dayofmonth  # initialization
        for ii in range(0, len(month)):
            if month[ii] != 1:
                N[ii] = dayofmonth[ii] + np.sum(DoM[1:month[ii] - 1], 2)

    if len(args) == 1:
        hour = args[0]
        N = np.floor((hour - 1) / 24) + 1

    TOA_h = solarconst * (1 + 0.03344 * cos(N * 2 * np.pi / 365.25 - 0.048869)) * sind(alpha)
    TOA_h = np.maximum(TOA_h, 0)

    return TOA_h


def coefficients(beta, ratio, R_b, A_i, f, *args):
    """
    This function creates three weighting matrices for the desired extent and with the desired resolution,
    that correspond to the gains/losses caused by tilting to each component of the incident radiation
    (direct, diffuse, and reflected

    The Input arguements are:
    - alpha: matrix of elevation angles
    - beta: matrix of tilting angles
    - azi: matrix of azimuth angles
    - orientation: matrix of surface azimuth (PV panel orientation angle)
    - hour, sunrise, and sunset (optional): hour of the year, and the matrix for the sunrise and sunset hours for every
    location on that day
    """

    duration = np.ones(beta.shape)  # Adjusts for any inaccuracy in case the sunset or the sunrise occurs within hour

    if len(args) == 3:
        hour = args[0]
        sunrise = args[1]
        sunset = args[2]
        hourofday = hour % 24 + 0.5

        # Case 1: Sunrise and sunset occur on the same day (GMT)
        critertion = np.logical_and(hourofday < sunrise, np.logical_and(sunrise > 0, sunset < 24))
        duration[critertion] = 0
        critertion = np.logical_and(np.logical_and((hourofday - 1) < sunrise, hourofday > sunrise),
                                    np.logical_and(sunrise > 0, sunset < 24))
        duration[critertion] = hourofday - sunrise[critertion]
        critertion = np.logical_and(np.logical_and((hourofday - 1) < sunset, hourofday > sunset),
                                    np.logical_and(sunrise > 0, sunset < 24))
        duration[critertion] = sunset[critertion] - hourofday + 1
        critertion = np.logical_and((hourofday - 1) > sunset,
                                    np.logical_and(sunrise > 0, sunset < 24))
        duration[critertion] = 0

        # Case 2: Sunrise occurs on the previous day (GMT) -
        # we can assume that it occurs at the same time as the following day
        critertion = np.logical_and(np.logical_and((hourofday - 1) < sunset, hourofday > sunset),
                                    np.logical_and(sunrise < 0, sunset < 24))
        duration[critertion] = sunset[critertion] - hourofday + 1
        critertion = np.logical_and(np.logical_and((hourofday - 1) > sunset, hourofday < sunrise + 24),
                                    np.logical_and(sunrise < 0, sunset < 24))
        duration[critertion] = 0
        critertion = np.logical_and(np.logical_and((hourofday - 1) < (sunrise + 24), hourofday > (sunrise + 24)),
                                    np.logical_and(sunrise < 0, sunset < 24))
        duration[critertion] = hourofday - (sunrise[critertion] + 24)

        # case 3: Sunset occurs on the next day (GMT) - not relevant for ASEAN
        critertion = np.logical_and(np.logical_and((hourofday - 1) < (sunset - 24), hourofday > (sunset - 24)),
                                    np.logical_and(sunrise > 0, sunset > 24))
        duration[critertion] = (sunset[critertion] - 24) - hourofday + 1
        critertion = np.logical_and(np.logical_and((hourofday - 1) > (sunset - 24), hourofday < sunrise),
                                    np.logical_and(sunrise > 0, sunset > 24))
        duration[critertion] = 0
        critertion = np.logical_and(np.logical_and((hourofday - 1) < sunrise, hourofday > sunrise),
                                    np.logical_and(sunrise > 0, sunset > 24))
        duration[critertion] = hourofday - sunrise[critertion]

    F_direct = duration * (1 - ratio + ratio * A_i) * R_b
    F_direct[F_direct < 0] = 0

    F_diffuse = ratio * (1 - A_i) * (1 + cos(np.deg2rad(beta))) / 2 * (1 + f * sin(np.deg2rad(beta / 2)) ** 3)

    F_reflected = (1 - cos(np.deg2rad(beta))) / 2

    return F_direct, F_diffuse, F_reflected


def loss(G_tilt_h, TEMP, A_Ross, pv):
    # This code creates a temperature loss weighting matrix for the desired extent
    # The inputs are:
    # G_tilt_h: incident radiation on the tilted panel
    # TEMP: matrix of ambient temperatures in °C
    # T_r: rated temperature according to STC
    # loss_coeff: % Heat loss coefficient according to the technical characteristics of the Yingli PV module

    T_cell = TEMP + A_Ross * G_tilt_h  # Cell temperature
    LOSS_TEMP = np.maximum(T_cell - pv["T_r"], 0) * pv["loss_coeff"] / 100
    return LOSS_TEMP


def global2diff(k_t, dims):
    # This code estimates the global-to-diffuse ratio.
    # The inputs are:

    k_d = np.zeros(dims)
    criterion = k_t <= 0.22
    k_d[criterion] = 1 - 0.09 * k_t[criterion]

    criterion = np.logical_and(k_t > 0.22, k_t <= 0.8)
    k_d[criterion] = 0.9511 - 0.1604 * k_t[criterion] + 4.388 * k_t[criterion] ** 2 - 16.638 * k_t[criterion] ** 3 \
                     + 12.336 * k_t[criterion] ** 4

    criterion = k_t > 0.8
    k_d[criterion] = 0.165

    A_ratio = np.minimum(k_d, 1)

    return A_ratio

def calc_CF_wind(hour, reg_ind, turbine, m, n, merraData, rasterData):
    ''' This function calculates the capacity factor for a given wind speed at 50m'''

    # Load MERRA data, increase its resolution, and fit it to the extent
    w50m_h = resizem(merraData["W50M"][:, :, hour], m, n)
    w50m_h = w50m_h[reg_ind]

    # Calculate the wind speed a the desired height
    w_new_h = w50m_h * rasterData["A_cf"]
    del w50m_h

    # Calculate the capacity factor

    a = turbine["w_in"] ** 3 / (turbine["w_in"] ** 3 - turbine["w_r"] ** 3)
    b = 1 / (turbine["w_r"] ** 3 - turbine["w_in"] ** 3)

    CF = np.zeros(w_new_h.shape)
    # Case 1 : above the cut-in speed and below the rated speed
    idx1 = np.logical_and(turbine["w_in"] < w_new_h, w_new_h < turbine["w_r"])
    CF[idx1] = a + b * w_new_h[idx1] ** 3
    # Case 2 : above the rated wind speed and below the cut_off speed
    idx2 = np.logical_and(turbine["w_r"] <= w_new_h, w_new_h <= turbine["w_off"])
    CF[idx2] = 1
    # Other cases (below cut-in or above cut-off
    CF[np.logical_not(np.logical_or(idx1, idx2))] = 0

    return CF


def calc_FLH_wind(hours, args):
    # Decomposing the tuple args
    paths = args[0]
    param = args[1]
    tech = args[2]
	
    m_high = param["m_high"]
    n_high = param["n_high"]
    reg_ind = param["Ind_nz"]

    turbine = param[tech]["technical"]

    # Obtain weather matrices
    merraData = {}
    # Wind speed at 50m
    merraData["W50M"] = hdf5storage.read('W50M', paths["W50M"])

    rasterData = {}
    # Calculate A matrices
    # A_cf
    with rasterio.open(paths["CORR"]) as src:
        w = src.read(1)
    rasterData["A_cf"] = np.flipud(w)
    rasterData["A_cf"] = rasterData["A_cf"][reg_ind]
    del w

    FLH = np.zeros(rasterData["A_cf"].shape)
    status = 0
    for hour in hours:
        if hour <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            sys.stdout.write('\r')
            sys.stdout.write(tech + ' ' + param["region"] + ' ' + '[%-50s] %d%%' % (
            '=' * ((status * 50) // len(hours)), (status * 100) // len(hours)))
            sys.stdout.flush()

        # Calculate hourly capacity factor
        CF = calc_CF_wind(hour, reg_ind, turbine, m_high, n_high, merraData, rasterData)

        # Aggregates CF to obtain the yearly FLH
        CF[np.isnan(CF)] = 0
        FLH = FLH + CF
    return FLH


def calc_TS_wind(hours, args):
    # Decomposing the tuple args
    paths = args[0]
    param = args[1]
    tech = args[2]

    m_high = param["m_high"]
    n_high = param["n_high"]
    reg_ind = param[tech]["Ind_points"]

    turbine = param[tech]["technical"]

    # Obtain weather matrices
    merraData = {}
    # Wind speed at 50m
    merraData["W50M"] = hdf5storage.read('W50M', paths["W50M"])

    rasterData = {}
    # Calculate A matrices
    # A_cf
    with rasterio.open(paths["CORR"]) as src:
        w = src.read(1)
    rasterData["A_cf"] = np.flipud(w)
    rasterData["A_cf"] = rasterData["A_cf"][reg_ind]
    del w

    TS = np.zeros((len(reg_ind[0]), 8760))
    status = 0
    for hour in hours:
        if hour <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            sys.stdout.write('\r')
            sys.stdout.write(tech + ' ' + param["region"] + ' ' + '[%-50s] %d%%' % (
            '=' * ((status * 50) // len(hours)), (status * 100) // len(hours)))
            sys.stdout.flush()

        # Calculate hourly capacity factor
        CF = calc_CF_wind(hour, reg_ind, turbine, m_high, n_high, merraData, rasterData)

        # Aggregates CF to obtain the time series
        CF[np.isnan(CF)] = 0
        TS[:, hour] = CF
    return TS
