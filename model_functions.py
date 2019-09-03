from data_functions import *
from util import *

np.seterr(divide='ignore')  # Repress Invalid value or division by zero error


def calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech):
    pv = param["PV"]["technical"]
    csp = param["CSP"]["technical"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    res_desired = param["res_desired"]
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

    # Check orientation parameter
    if 'orientation' in pv.keys():
        orient = pv["orientation"]
    else:
        orient = 0

    # Compute the angles
    A_phi, A_omega, A_delta, A_alpha, A_beta, A_azimuth, A_orientation = \
        angles(hour, reg_ind_h, Crd_all, res_desired, orient)

    # Compute the hourly TOA radiation
    TOA_h = toa_hourly(A_alpha, hour)

    # If all TOA values are zero, return to main function
    if (TOA_h == 0).all():
        CF_pv = np.zeros(reg_ind[0].shape)
        CF_csp = np.zeros(reg_ind[0].shape)
        return CF_pv, CF_csp

    CLEARNESS_h = CLEARNESS_h[reg_ind_h] * param["PV"]["resource"]["clearness_correction"]
    TEMP_h = merraData["T2M"][:, :, hour]
    TEMP_h = resizem(TEMP_h, m_high, n_high) - 273.15  # Convert to Celsius
    TEMP_h = TEMP_h[reg_ind_h]

    # Other matrices
    A_albedo = rasterData["A_albedo"][reg_ind_h]
    A_Ross = rasterData["A_Ross"][reg_ind_h]

    # Compute the ratio of diffuse radiation
    RATIO = global2diff(CLEARNESS_h, A_alpha.shape)
    A_i = (1 - RATIO) * CLEARNESS_h
    f = (1 - RATIO) ** 0.5

    # Compute the shading losses
    # Currently ignored
    SHADING = 0

    if tech == 'PV':
        # Tracking
        if pv["tracking"] == 1:
            A_orientation, A_beta = tracking(1, A_phi, A_alpha, A_beta, A_azimuth)
        elif pv["tracking"] == 2:
            A_orientation, A_beta = tracking(2, A_phi, A_alpha, A_beta, A_azimuth)

        aux = np.maximum(np.minimum((sind(A_delta) * sind(A_phi) * cosd(A_beta)
                                     - sind(A_delta) * cosd(A_phi) * sind(A_beta) * cosd(A_orientation)
                                     + cosd(A_delta) * cosd(A_phi) * cosd(A_beta) * cosd(A_omega)
                                     + cosd(A_delta) * sind(A_phi) * sind(A_beta) * cosd(A_orientation) * cosd(A_omega)
                                     + cosd(A_delta) * sind(A_beta) * sind(A_orientation) * sind(A_omega)),
                                    1), -1)
        A_incidence = arccosd(aux)
        # Compute the coefficients for the HDKR model
        R_b = cosd(A_incidence) / sind(A_alpha)
        R_b[A_alpha <= 5] = cosd(A_incidence[A_alpha <= 5]) / sind(5)
        R_b[A_alpha <= 0] = 0

        F_direct, F_diffuse, F_reflected = coefficients(A_beta, RATIO, R_b, A_i, f)

        F = F_diffuse + F_direct * (1 - SHADING) + F_reflected * A_albedo
        F[F > 1] = 1

        # Compute the incident radiation
        GHI_h = TOA_h * CLEARNESS_h
        GHI_h[np.isnan(GHI_h)] = 0
        G_tilt_h = GHI_h * F

        # Compute losses due to heating of the PV cells
        LOSS_TEMP = loss(G_tilt_h, TEMP_h, A_Ross, pv)

        # Compute the hourly capacity factor
        CF_pv = G_tilt_h * (1 - LOSS_TEMP) / 1000

        CF_pv[A_alpha <= 0] = 0
        # Adjusting the length of the matrices
        aux = np.zeros(len(reg_ind[0]))
        aux[filter] = CF_pv
        CF_pv = aux
    else:
        CF_pv = None

    if tech == 'CSP':
        # Wind Speed Corrected at 2m
        w2m_h = resizem(merraData["W50M"][:, :, hour], m_high, n_high)
        w2m_h = w2m_h[reg_ind] * rasterData["A_WindSpeed_Corr"][reg_ind]

        # Wind Speed cutoff filter:
        windfilter = w2m_h >= csp["Wind_cutoff"]

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
        S = TOA_h * CLEARNESS_h * F_direct_csp * (1 - SHADING)
        Qu = csp["Flow_coeff"] * (S - csp["AbRe_ratio"] * (csp["loss_coeff"] + csp["loss_coeff_wind"] * w2m_h ** 2)
                                  * (csp["T_avg_HTF"] - TEMP_h))
        CF_csp = Qu / 1000
        CF_csp[CF_csp < 0] = 0
        CF_csp[CF_csp > 1] = 1

        if windfilter.any():
            CF_csp[windfilter] = 0

        aux = np.zeros(len(reg_ind[0]))
        aux[filter] = CF_csp
        CF_csp = aux
    else:
        CF_csp = None

    return CF_pv, CF_csp


def calc_FLH_solar(hours, args):
    # Decomposing the tuple args
    paths = args[0]
    param = args[1]
    tech = args[2]
    rasterData = args[3]
    merraData = args[4]
    reg_ind = param["Ind_nz"]

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
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech)[0]
        elif tech == 'CSP':
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech)[1]

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

    rasterData = {}
    # Calculate A matrices
    # A_lu
    with rasterio.open(paths["LU"]) as src:
        w = src.read(1)
    rasterData["A_lu"] = np.flipud(w)
    # A_Ross (Temperature coefficients for heating losses)
    rasterData["A_Ross"] = changem(rasterData["A_lu"], param["landuse"]["Ross_coeff"], param["landuse"]["type"]).astype(
        float)
    # A_albedo (Reflectivity coefficients)
    rasterData["A_albedo"] = changem(rasterData["A_lu"], param["landuse"]["albedo"], param["landuse"]["type"]).astype(
        float)
    # A_WS_Coef wind Speed at 2m above the ground
    A_hellmann = changem(rasterData["A_lu"], landuse["hellmann"], landuse["type"]).astype(float)
    rasterData["A_WindSpeed_Corr"] = (2 / 50) ** A_hellmann
    del A_hellmann

    # Obtain weather matrices
    merraData = {}
    # Clearness index - stored variable CLEARNESS
    merraData["CLEARNESS"] = hdf5storage.read('CLEARNESS', paths["CLEARNESS"])
    # Temperature 2m above the ground - stored variable T2M
    merraData["T2M"] = hdf5storage.read('T2M', paths["T2M"])
    # Wind Speed
    merraData["W50M"] = hdf5storage.read('W50M', paths["W50M"])

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
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech)[0]
        elif tech == 'CSP':
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech)[1]

        # Aggregates CF to obtain the time series
        CF[np.isnan(CF)] = 0
        TS[:, hour] = CF
    return TS


def angles(hour, reg_ind, Crd_all, res_desired, orient):
    """
    This function creates six matrices for the desired extent, that represent the elevation, tilt, azimuth and
    oerientation angles, in addition to the sunrise and sunset hours of every pixel with the desired resolution

    :param hour: hour rank in a year (from 1 to 8760
    :param args:
        set1
        north, east, south, west: desired geographic extent
        res: resolution of MERRA data & desired resolution in lat/lon
        Ext_in: extent of the global are of interest
        set2
        (lat, lon): latitude and longitude of all points of interest
    """

    # Initialization
    Crd_points = crd_exact_points(reg_ind, Crd_all, res_desired)
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

    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat,
                                                                      range_lon)] - 35) / 65 * 45 + 35  # Europe
    range_lat = np.logical_and(lat >= 20, lat < 65)
    range_lon = np.logical_and(lon >= 75, lon < 140)

    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat,
                                                                      range_lon)] - 20) / 65 * 60 + 20  # Asia/China

    # Azimuth angle (in degrees)
    aux = np.maximum(np.minimum(sind(delta) * cosd(phi) - cosd(delta) * sind(phi) * cosd(omega) / cosd(alpha), 1), -1)
    aziam = arccosd(aux)
    azipm = 360 - aziam
    azi = aziam * ((omega < 0) * 1) + azipm * ((omega >= 0) * 1)

    # Orientation (in degrees)
    orientation = np.full(alpha.shape, orient)  # Azimuth of the PV panel is zero for the Northern hemisphere
    orientation[phi < 0] = 180 - orient  # Azimuth of the PV panel is 180° for the Southern hemisphere

    return phi, omega, delta, alpha, beta, azi, orientation


def tracking(axis, A_phi, A_alpha, A_beta, A_azimuth):
    # One axis Tracking
    if axis == 1:
        A_orientation = np.zeros(A_alpha.shape)
        A_orientation[A_phi < 0] = 180

        x = (cosd(A_alpha) * (-sind(A_azimuth - A_orientation))) \
            / (cosd(A_alpha) * (-cosd(A_azimuth - A_orientation)) * sind(A_beta) + sind(A_alpha) * cosd(A_beta))

        # The difference should be the angular displacement
        Az_dif = A_azimuth - (A_orientation + 180)
        Az_dif[Az_dif < -180] = Az_dif[Az_dif < -180] + 360
        Az_dif[Az_dif > 180] = Az_dif[Az_dif > 180] - 360

        # y is used to locate R on the right quadrant
        y = np.zeros(x.shape())
        crit = np.logical_or(x == 0,
                             np.logical_or(np.logical_and(x > 0, Az_dif > 0), np.logical_and(x < 0, Az_dif < 0)))
        y[crit] = 0
        crit = np.logical_and(x < 0, Az_dif > 0)
        y[crit] = 180
        crit = np.logical_and(x > 0, Az_dif < 0)
        y[crit] = -180

        # Tracking angle
        R = arctand(x) + y

        # New beta
        beta = arccosd(cosd(R) * cosd(A_beta))

        # New orientation
        orientation = np.zeros(A_beta.shape())
        crit = np.logical_and(A_beta != 0, np.logical_and(-90 <= R, R <= 90))
        orientation[crit] = A_orientation[crit] + 180 + arcsind(sind(R[crit]) / sind(beta[crit]))
        crit = np.logical_and(-180 <= R, R < -90)
        orientation[crit] = A_orientation[crit] - arcsind(sind(R[crit]) / sind(beta[crit]))
        crit = np.logical_and(90 < R, R <= 180)
        orientation[crit] = A_orientation[crit] + 360 - arcsind(sind(R[crit]) / sind(beta[crit]))

        A_orientation = orientation - 180
        A_beta = beta

    # Two Axis Tracking
    elif axis == 2:
        A_beta = 90 - A_alpha
        A_orientation = A_azimuth - 180

    return A_orientation, A_beta


def toa_hourly(alpha, hour):
    solarconst = 1367  # in W/m^2
    N = hour // 24 + 1

    TOA_h = solarconst * (1 + 0.03344 * cos(N * 2 * np.pi / 365.25 - 0.048869)) * sind(alpha)
    TOA_h = np.maximum(TOA_h, 0)

    return TOA_h


def coefficients(beta, ratio, R_b, A_i, f):
    """
    This function creates three weighting matrices for the desired extent and width the desired resolution,
    that correspond to the gains/losses caused by tilting to each component of the incident radiation
    (direct, diffuse, and reflected

    The Input arguements are:
    - alpha: matrix of elevation angles
    - beta: matrix of tilting angles
    - azi: matrix of azimuth angles
    - orientation: matrix of surface azimuth (PV panel orientation angle)
    """

    F_direct = (1 - ratio + ratio * A_i) * R_b
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
    rasterData = args[3]
    merraData = args[4]
    m_high = param["m_high"]
    n_high = param["n_high"]
    reg_ind = param["Ind_nz"]

    turbine = param[tech]["technical"]

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
    if tech == 'WindOn':
        with rasterio.open(paths["CORR_ON"]) as src:
            w = src.read(1)
    if tech == 'WindOff':
        with rasterio.open(paths["CORR_OFF"]) as src:
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
