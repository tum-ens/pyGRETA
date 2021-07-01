from .spatial_functions import crd_exact_points
from .util import *

np.seterr(divide="ignore")  # Repress invalid value or division by zero error


def calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech):
    """
    This function computes the hourly capacity factor for PV and CSP technologies for all valid pixels within
    the spatial scope for a given hour.

    :param hour: Hour within the year (from 0 to 8759).
    :type hour: integer
    :param reg_ind: indices of valid pixels within the spatial scope (pixels on land).
    :type reg_ind: tuple of arrays
    :param param: Dictionary including the desired resolution, the coordinates of the bounding box of the spatial scope, and technology parameters.
    :type param: dict
    :param merraData: Dictionary of numpy arrays containing the weather data for every point in *reg_ind*.
    :type merraData: dict
    :param rasterData: Dictionary of numpy arrays containing land use types, Ross coefficients, albedo coefficients, and wind speed correction for every point in *reg_ind*.
    :type rasterData: dict
    :param tech: Name of the technology (``'PV'`` or ``'CSP'``).
    :type tech: str

    :return (CF_pv, CF_csp): the capacity factors for all the points during that hour for PV and CSP.
    :rtype: tuple (numpy array, numpy array)
    """
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
    if "orientation" in pv.keys():
        orient = pv["orientation"]
    else:
        orient = 0

    # Compute the angles
    A_phi, A_omega, A_delta, A_alpha, A_beta, A_azimuth, A_orientation = angles(hour, reg_ind_h, Crd_all, res_desired, orient)

    # Compute the hourly TOA radiation
    TOA_h = toa_hourly(A_alpha, hour)

    # If all TOA values are zero, return to main function
    if (TOA_h == 0).all():
        CF_pv = np.zeros(reg_ind[0].shape)
        CF_csp = np.zeros(reg_ind[0].shape)
        return CF_pv, CF_csp

    CLEARNESS_h = CLEARNESS_h[reg_ind_h] * param[tech]["resource"]["clearness_correction"]
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

    if tech == "PV":
        # Tracking
        if pv["tracking"] == 1:
            A_orientation, A_beta = tracking(1, A_phi, A_alpha, A_beta, A_azimuth)
        elif pv["tracking"] == 2:
            A_orientation, A_beta = tracking(2, A_phi, A_alpha, A_beta, A_azimuth)

        aux = np.maximum(
            np.minimum(
                (
                    sind(A_delta) * sind(A_phi) * cosd(A_beta)
                    - sind(A_delta) * cosd(A_phi) * sind(A_beta) * cosd(A_orientation)
                    + cosd(A_delta) * cosd(A_phi) * cosd(A_beta) * cosd(A_omega)
                    + cosd(A_delta) * sind(A_phi) * sind(A_beta) * cosd(A_orientation) * cosd(A_omega)
                    + cosd(A_delta) * sind(A_beta) * sind(A_orientation) * sind(A_omega)
                ),
                1,
            ),
            -1,
        )
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

    if tech == "CSP":
        # Wind Speed Corrected at 2m
        w2m_h = resizem(merraData["W50M"][:, :, hour], m_high, n_high)
        w2m_h = w2m_h[reg_ind] * rasterData["A_WindSpeed_Corr"][reg_ind]

        # Wind Speed cutoff filter:
        windfilter = w2m_h >= csp["Wind_cutoff"]

        # For CSP: tracking like pv.tracking = 1
        A_beta = 90 - A_alpha
        aux = np.maximum(
            np.minimum(
                (
                    sind(A_delta) * sind(A_phi) * cosd(A_beta)
                    - sind(A_delta) * cosd(A_phi) * sind(A_beta) * cosd(A_orientation)
                    + cosd(A_delta) * cosd(A_phi) * cosd(A_beta) * cosd(A_omega)
                    + cosd(A_delta) * sind(A_phi) * sind(A_beta) * cosd(A_orientation) * cosd(A_omega)
                    + cosd(A_delta) * sind(A_beta) * sind(A_orientation) * sind(A_omega)
                ),
                1,
            ),
            -1,
        )
        A_incidence = arccosd(aux)
        R_b = cosd(A_incidence) / sind(A_alpha)
        F_direct_csp, _, _ = coefficients(90 - A_alpha, RATIO, R_b, A_i, f)
        S = TOA_h * CLEARNESS_h * F_direct_csp * (1 - SHADING)
        Qu = csp["Flow_coeff"] * (S - csp["AbRe_ratio"] * (csp["loss_coeff"] + csp["loss_coeff_wind"] * w2m_h ** 2) * (csp["T_avg_HTF"] - TEMP_h))
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


def angles(hour, reg_ind, Crd_all, res_desired, orient):
    """
    This function creates multiple matrices for the whole scope, that represent the incidence, hour angles, declination,
    elevation, tilt, azimuth and orientation angles of every pixel with the desired resolution.

    :param hour: Hour rank in a year (from 0 to 8759).
    :type hour: int
    :param reg_ind: indices of valid pixels within the spatial scope (pixels on land).
    :type reg_ind: tuple of arrays
    :param Crd_all: Coordinates of the bounding box of the spatial scope.
    :type Crd_all: list
    :param res_desired: Desired high resolution in degrees.
    :type res_desired: list
    :param orient: Azimuth orientation of the module in degrees.
    :type orient: int

    :return (phi, omega, delta, alpha, beta, azi, orientation): Rasters of latitude, hour, declination, elevation,
        tilt, azimuth and orientation angles.
    :rtype: tuple of arrays
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
    delta = np.tile(arcsind(0.3978 * sin(N * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(N * 2 * np.pi / 365.25 - 0.0489))), k)
    delta[phi < 0] = -delta[phi < 0]

    # Elevation angle (in degrees)
    alpha = arcsind(sind(delta) * sind(phi) + cosd(delta) * cosd(phi) * cosd(omega))

    # Optimal tilt angle (loosely based on Breyer 2010)
    beta = np.minimum(np.abs(phi), 55)  # The tilt angle is preferably equal to the latitude
    range_lat = np.logical_and(np.abs(phi) >= 35, np.abs(phi) < 65)
    beta[range_lat] = (beta[range_lat] - 35) / 65 * 55 + 35  # Tilt angle does not increase very quickly
    range_lat = np.logical_and(lat >= 35, lat < 65)
    range_lon = np.logical_and(lon >= -20, lon < 30)

    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat, range_lon)] - 35) / 65 * 45 + 35  # Europe
    range_lat = np.logical_and(lat >= 20, lat < 65)
    range_lon = np.logical_and(lon >= 75, lon < 140)

    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat, range_lon)] - 20) / 65 * 60 + 20  # Asia/China

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
    """
    This function computes the tilt angle and orientation based on the type of tracking, incidence, elevation tilt
    and azimuth angles.

    :param axis: Number of tracking axes (0, 1, 2). The value ``0`` means no tracking (fixed rack),
        ``1`` means single-axis tracking in east-west dimension, and ``2`` means double-axis tracking.

    :type axis: int
    :param A_phi: Raster of latitude angle.
    :type A_phi: numpy array
    :param A_alpha: Raster of elevation angle.
    :type A_alpha: numpy array
    :param A_beta: Raster of tilt angle.
    :type A_beta: numpy array
    :param A_azimuth: Raster of azimuth angle.
    :type A_azimuth: numpy array

    :return (A_orient, A_beta): Tuple of rasters for orientationa and tilt angles for specified tracking type.
    :rtype: tuple of arrays
    """
    # Single-axis tracking
    if axis == 1:
        A_orientation = np.zeros(A_alpha.shape)
        A_orientation[A_phi < 0] = 180

        x = (cosd(A_alpha) * (-sind(A_azimuth - A_orientation))) / (
            cosd(A_alpha) * (-cosd(A_azimuth - A_orientation)) * sind(A_beta) + sind(A_alpha) * cosd(A_beta)
        )

        # The difference should be the angular displacement
        Az_dif = A_azimuth - (A_orientation + 180)
        Az_dif[Az_dif < -180] = Az_dif[Az_dif < -180] + 360
        Az_dif[Az_dif > 180] = Az_dif[Az_dif > 180] - 360

        # y is used to locate R on the right quadrant
        y = np.zeros(x.shape)
        crit = np.logical_or(x == 0, np.logical_or(np.logical_and(x > 0, Az_dif > 0), np.logical_and(x < 0, Az_dif < 0)))
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
        orientation = np.zeros(A_beta.shape)
        crit = np.logical_and(A_beta != 0, np.logical_and(-90 <= R, R <= 90))
        orientation[crit] = A_orientation[crit] + 180 + arcsind(sind(R[crit]) / sind(beta[crit]))
        crit = np.logical_and(-180 <= R, R < -90)
        orientation[crit] = A_orientation[crit] - arcsind(sind(R[crit]) / sind(beta[crit]))
        crit = np.logical_and(90 < R, R <= 180)
        orientation[crit] = A_orientation[crit] + 360 - arcsind(sind(R[crit]) / sind(beta[crit]))

        A_orientation = orientation - 180
        A_beta = beta

    # Double-axis tracking
    elif axis == 2:
        A_beta = 90 - A_alpha
        A_orientation = A_azimuth - 180

    return A_orientation, A_beta


def toa_hourly(alpha, hour):
    """
    This function returns the top of the atmosphere normal irradiance based on the solar constant, hour rank, and incidence angle.

    :param alpha: Raster of elevation angles.
    :type alpha: numpy array
    :param hour: Hour rank of the year (from 0 to 8759).
    :type hour: int

    :return TOA_h: Raster of the normal top of the atmosphere irradiance.
    :rtype: numpy array
    """
    solarconst = 1367  # in W/m^2
    N = hour // 24 + 1

    TOA_h = solarconst * (1 + 0.03344 * cos(N * 2 * np.pi / 365.25 - 0.048869)) * sind(alpha)
    TOA_h = np.maximum(TOA_h, 0)

    return TOA_h


def coefficients(beta, ratio, R_b, A_i, f):
    """
    This function creates three weighting matrices for the spatial scope with the desired resolution,
    that correspond to the gains/losses caused by tilting to each component of the incident irradiance
    (direct, diffuse, and reflected).

    :param beta: Raster of tilt angles.
    :type beta: numpy array
    :param ratio: Diffuse fraction of global horizontal solar radiation using the Erbs model.
    :type ratio: numpy array
    :param R_b: Ratio of incident beam to horizontal beam in the HDKR model.
    :type R_b: numpy array
    :param A_i: Anisotropy index for forward scattering circumsolar diffuse irradiance in the HDKR model.
    :type A_i: numpy array
    :param f: Modulating factor for horizontal brightening correction.
    :type f: numpy array

    :return (F_direct, F_diffuse, F_reflected): Rasters of direct, diffuse and reflected ratios of irradiance.
    :rtype: tuple of arrays
    """

    F_direct = (1 - ratio + ratio * A_i) * R_b
    F_direct[F_direct < 0] = 0

    F_diffuse = ratio * (1 - A_i) * (1 + cos(np.deg2rad(beta))) / 2 * (1 + f * sin(np.deg2rad(beta / 2)) ** 3)

    F_reflected = (1 - cos(np.deg2rad(beta))) / 2

    return F_direct, F_diffuse, F_reflected


def loss(G_tilt_h, TEMP, A_Ross, pv):
    """
    This function creates a temperature loss weighting matrix for the spatial scope.

    :param G_tilt_h: Raster of incident irradiance on the tilted panel.
    :type G_tilt_h: numpy array
    :param TEMP: Raster of ambient temperatures in °C
    :type TEMP: numpy array
    :param A_Ross: Raster of Ross coefficients for temperature sensitivity.
    :type A_Ross: numpy array
    :param pv: Dictionary containing PV-specific parameters for loss coefficient and rated temperature.
    :type pv: dict

    :return LOSS_TEMP: raster of weighting temperature loss.
    :rtype: numpy array
    """

    T_cell = TEMP + A_Ross * G_tilt_h  # Cell temperature
    LOSS_TEMP = np.maximum(T_cell - pv["T_r"], 0) * pv["loss_coeff"] / 100
    return LOSS_TEMP


def global2diff(k_t, dims):
    """
    This function estimates the global-to-diffuse irradiance ratio using the Erb model.

    :param k_t: Raster of clearness indices.
    :type k_t: numpy array
    :param dims: Dimensions of the output (similar to the dimension of the angles).
    :type dims: tuple

    :return A_ratio: Raster of global-to-diffuse irradiance ratios.
    :rtype: numpy array
    """
    k_d = np.zeros(dims)
    criterion = k_t <= 0.22
    k_d[criterion] = 1 - 0.09 * k_t[criterion]

    criterion = np.logical_and(k_t > 0.22, k_t <= 0.8)
    k_d[criterion] = 0.9511 - 0.1604 * k_t[criterion] + 4.388 * k_t[criterion] ** 2 - 16.638 * k_t[criterion] ** 3 + 12.336 * k_t[criterion] ** 4

    criterion = k_t > 0.8
    k_d[criterion] = 0.165

    A_ratio = np.minimum(k_d, 1)

    return A_ratio


def calc_CF_wind(hour, reg_ind, turbine, m, n, merraData, rasterData):
    """
    This function computes the hourly capacity factor for onshore and offshore wind for all valid pixels within
    the spatial scope for a given hour.

    :param hour: Hour within the year (from 0 to 8759).
    :type hour: integer
    :param reg_ind: indices of valid pixels within the spatial scope (pixels on land for onshore wind, on sea for offshore wind).
    :type reg_ind: tuple of arrays
    :param turbine: Dictionary including the turbine parameters (cut-in, cut-off and rated wind speed).
    :type turbine: dict
    :param m: number of rows.
    :type m: int
    :param n: number of columns.
    :type n: int
    :param merraData: Dictionary of numpy arrays containing the weather data for every point in *reg_ind*.
    :type merraData: dict
    :param rasterData: Dictionary of numpy arrays containing the wind speed correction for every point in *reg_ind*.
    :type rasterData: dict

    :return CF: Capacity factors for all the valid points during that hour.
    :rtype: numpy array
    """

    # Load MERRA data, increase its resolution, and fit it to the extent
    w50m_h = resizem(merraData["W50M"][:, :, hour], m, n)
    w50m_h = w50m_h[tuple(reg_ind)]

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
