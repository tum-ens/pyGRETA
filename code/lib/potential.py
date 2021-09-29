from . import physical_models as pm
from . import spatial_functions as sf
from . import util as ul
from .log import logger
import multiprocessing as mp
import itertools as it
import pandas as pd
import numpy as np
import traceback
import rasterio
import hdf5storage

def calculate_full_load_hours(paths, param, tech, multiprocessing):
    """
    This function calculates the yearly FLH for a technology for all valid pixels in a spatial scope. Valid pixels are land pixels
    for WindOn, PV and CSP, and sea pixels for WindOff. The FLH values are calculated by summing up hourly capacity factors.

    :param paths: Dictionary of dictionaries containing the paths to the input weather data, land, sea and land use rasters, and correction rasters.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the spatial scope, and technology and computation parameters.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str
    :param multiprocessing: Determines if the computation uses multiprocessing (True/False)
    :type multiprocessing: bool
    :return: The raster of FLH potential is saved as mat and tif files, along with the json metadata file.
    :rtype: None
    """
    logger.info("Start - Region: " + param["region_name"])

    if tech in ["WindOn", "WindOff"]:
        logger.info(tech + " - HUB_HEIGHTS: " + str(param[tech]["technical"]["hub_height"]))
    elif tech in ["OpenFieldPV","RoofTopPV"]:
        if "orientation" in param[tech]["technical"].keys():
            logger.info(tech + " - Orientation: " + str(param[tech]["technical"]["orientation"]))

    nproc = param["nproc"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    m_low = param["m_low"]
    n_low = param["n_low"]
    CPU_limit = np.full((1, nproc), param["CPU_limit"])
    res_weather = param["res_weather"]
    Crd_all = param["Crd_all"]
    Ind = sf.ind_merra(Crd_all, Crd_all, res_weather)[0]
    
    if tech == "WindOff":
        w = np.flipud(hdf5storage.read("EEZ", paths["EEZ"]))
        # with rasterio.open(paths["EEZ"]) as src:
        #     w = src.read(1)
    else:
        with rasterio.open(paths["LAND"]) as src:
            w = src.read(1)
    
    param["Ind_nz"] = np.nonzero(np.flipud(w))

    del w    
    
    # Obtain weather and correction matrices
    merraData, rasterData = get_merra_raster_data(paths, param, tech)
        
    if tech in ["OpenFieldPV", "RoofTopPV", "CSP"]:

        day_filter = np.nonzero(merraData["CLEARNESS"][Ind[2] - 1: Ind[0], Ind[3] - 1: Ind[1], :].sum(axis=(0, 1)))
        list_hours = np.arange(0, 8760)
        # if nproc == 1:
        #     param["status_bar_limit"] = list_hours[-1]
        #     results = calc_FLH_solar(list_hours[day_filter], [param, tech, rasterData, merraData])
        # else:
        #     list_hours = np.array_split(list_hours[day_filter], nproc)
        #     param["status_bar_limit"] = list_hours[0][-1]
        #     results = mp.Pool(processes=nproc, initializer=ul.limit_cpu, initargs=CPU_limit).starmap(
        #         calc_FLH_solar, it.product(list_hours, [[param, tech, rasterData, merraData]])
        #     )
        #  # Collecting results
        # FLH_low = np.zeros((m_low, n_low))
        #
        # if nproc > 1:
        #     for p in range(len(results)):
        #         FLH_low = FLH_low + results[p]
        # else:
        #     FLH_low = results

        param["status_bar_limit"] = list_hours[-1]
        FLH_low = calc_FLH_solar(list_hours[day_filter], [param, tech, rasterData, merraData])
        
        FLH_high = ul.resizem(FLH_low, m_high, n_high)
        FLH = np.full((m_high, n_high),np.nan)

        FLH[param["Ind_nz"]]=FLH_high[param["Ind_nz"]]
        hdf5storage.writes({"FLH": FLH}, paths[tech]["FLH"], store_python_metadata=True, matlab_compatible=True)
        ul.create_json(
            paths[tech]["FLH"],
            param,
            ["author", "comment", tech, "region_name", "subregions_name", "year", "res_desired", "res_weather"],
            paths,
            ["subregions"],
        )
        logger.info("files saved: " + paths[tech]["FLH"])

        # Save GEOTIFF files
        if param["savetiff_potentials"]:
            GeoRef = param["GeoRef"]
            sf.array2raster(ul.changeExt2tif(paths[tech]["FLH"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], FLH)
            logger.info("files saved:" + ul.changeExt2tif(paths[tech]["FLH"]))
    
    elif tech in ["WindOff"]:
        list_hours = np.array_split(np.arange(0, 8760), nproc)
        param["status_bar_limit"] = list_hours[0][-1]
        results = mp.Pool(processes=nproc, initializer=ul.limit_cpu, initargs=CPU_limit).starmap(
            calc_FLH_windoff, it.product(list_hours, [[param, tech, rasterData, merraData]])
        )
        # Collecting results
        FLH = np.full((m_high, n_high), np.nan)
        FLH[param["Ind_nz"]] = 0
        if nproc > 1:
            for p in range(len(results)):
                FLH[param["Ind_nz"]] = FLH[param["Ind_nz"]] + results[p]
        else:
            FLH[param["Ind_nz"]] = results

        hdf5storage.writes({"FLH": FLH}, paths[tech]["FLH"], store_python_metadata=True, matlab_compatible=True)
        ul.create_json(
            paths[tech]["FLH"],
            param,
            ["author", "comment", tech, "region_name", "subregions_name", "year", "res_desired", "res_weather"],
            paths,
            ["subregions"],
        )
        logger.info("files saved: " + paths[tech]["FLH"])

        # Save GEOTIFF files
        if param["savetiff_potentials"]:
            GeoRef = param["GeoRef"]
            sf.array2raster(ul.changeExt2tif(paths[tech]["FLH"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], FLH)
            logger.info("files saved:" + ul.changeExt2tif(paths[tech]["FLH"]))
    
    elif tech in ["WindOn"]:
        merraData = merraData["W50M"][::-1, :, :]

        rasterData = rasterData["A_cf"]
        #rasterData = rasterData.astype(dtype=np.float32)

        with rasterio.open(paths["GWA_global"]) as src:
            GWA_array = src.read(1)
        GWA_array[np.isnan(GWA_array)] = 0

        #merraData = merraData.astype(dtype=np.float32)
        GWA_array = GWA_array.astype(dtype=np.float16)
        #rasterData = rasterData.astype(dtype=np.float32)

        # param["status_bar_limit"] = list_rows_splitted[0][-1] # ToDo: Wo we need this?

        # -------------------------------------------------------------------------------------
        # Multiprocessing by Patrick 20210615
        # FIXME Import what?
        b_xmin = hdf5storage.read("MERRA_XMIN", paths[
            "MERRA_XMIN"])  # ToDo: move into calculate_fullload_hours? # ToDo: Put outside of function?
        b_xmax = hdf5storage.read("MERRA_XMAX", paths["MERRA_XMAX"])
        b_ymin = hdf5storage.read("MERRA_YMIN", paths["MERRA_YMIN"])
        b_ymax = hdf5storage.read("MERRA_YMAX", paths["MERRA_YMAX"])

        # FIXME Import What?
        x_gwa = hdf5storage.read("GWA_X", paths["GWA_X"])
        y_gwa = hdf5storage.read("GWA_Y", paths["GWA_Y"])

        processes = []  # Store all single process of multiprocessing
        list_results = mp.RawArray('f', m_high * n_high)
        FLH = np.frombuffer(list_results, dtype=np.float32).reshape(m_high, n_high)
        FLH[:] = np.nan
        list_pixles = np.arange(n_low * m_low)  # All pixles within MERRA data

        # multiprocessing = True  # debuging
        if multiprocessing:
            list_pixles_splitted = np.array_split(list_pixles,nproc)  # Splitted list acording to the number of parrallel processes
            logger.info('# of processes: ' + str(len(list_pixles_splitted)))
            logger.debug(list_pixles_splitted)

            for pixles in list_pixles_splitted:  # Run the 'calc_FLH_windon' for each of the splitted rows
                p = mp.Process(target=calc_FLH_windon, args=(
                    param, tech, rasterData, merraData, GWA_array, b_xmin, b_xmax, b_ymin, b_ymax, x_gwa, y_gwa,
                    pixles, list_results))
                processes.append(p)

            logger.debug('Starting processes for wind computation')
            for p in processes:
                p.start()  # Start all single processes
            logger.info('All processes started')

            for p in processes:
                p.join()  # Wait until all processes are finished
            logger.info('All processes finished')
        else:
            # list_pixles = np.arange(100,102)  # debuging
            # list_pixles = [np.ravel_multi_index((1, 2), (m_low, n_low))]  # debuging
            logger.debug(list_pixles)
            calc_FLH_windon(param, tech, rasterData, merraData, GWA_array, b_xmin, b_xmax, b_ymin, b_ymax, x_gwa,
                            y_gwa, list_pixles, list_results)
            logger.info('Calculations for all pixels done')

        # ------------------------------------------------------------------------------------

        # FLH[A_country_area==0] = float("nan")python
        FLH = np.flipud(FLH)
        FLH_scope = np.full((m_high, n_high), np.nan)
        FLH_scope[param["Ind_nz"]] = FLH[param["Ind_nz"]]

        hdf5storage.writes({"FLH": FLH}, paths[tech]["FLH"], store_python_metadata=True, matlab_compatible=True)
        logger.info("\nfiles saved: " + paths[tech]["FLH"])

        if param["savetiff_potentials"]:
            GeoRef = param["GeoRef"]
            sf.array2raster(ul.changeExt2tif(paths[tech]["FLH"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"],
                         GeoRef["pixelHeight"], FLH)
            logger.info("files saved:" +ul.changeExt2tif(paths[tech]["FLH"]))

    logger.debug("End")


def get_merra_raster_data(paths, param, tech): #ToDo clean up unnecessary things
    """
    This function returns a tuple of two dictionaries containing weather and correction rasters for specified technology.

    :param paths: Dictionary of dictionaries containing the paths to the input weather and raster data.
    :type paths: dict
    :param param: Dictionary of dictionaries containing land use, Ross coefficients, albedo, and Hellmann coefficients.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return (merraData, rasterData): Dictionaries for the weather data and for the correction data.
    :rtype: tuple (dict, dict)
    """
    landuse = param["landuse"]
    merraData = {}
    rasterData = {}
    # Wind Speed Data
    merraData["W50M"] = hdf5storage.read("W50M", paths["W50M"])
    if tech in ["OpenFieldPV", "RoofTopPV", "CSP"]:

        # Other weather Data
        # Clearness index - stored variable CLEARNESS
        merraData["CLEARNESS"] = hdf5storage.read("CLEARNESS", paths["CLEARNESS"])
        # Temperature 2m above the ground - stored variable T2M
        merraData["T2M"] = hdf5storage.read("T2M", paths["T2M"])

        # Calculate A matrices correction
        # A_lu
        w = np.flipud(hdf5storage.read("LU", paths["LU"]))
        # with rasterio.open(ul.changeExt2tif(paths["LU"])) as src:
        #     w = src.read(1)
        rasterData["A_lu"] = np.flipud(w)
        # A_Ross (Temperature coefficients for heating losses)
        rasterData["A_Ross"] = ul.changem(rasterData["A_lu"], param["landuse"]["Ross_coeff"],
                                       param["landuse"]["type"]).astype("float16")
        # A_albedo (Reflectivity coefficients)
        rasterData["A_albedo"] = ul.changem(rasterData["A_lu"], param["landuse"]["albedo"],
                                         param["landuse"]["type"]).astype("float16")
        # A_WS_Coef wind Speed at 2m above the ground
        A_hellmann = ul.changem(rasterData["A_lu"], landuse["hellmann"], landuse["type"])
        rasterData["A_WindSpeed_Corr"] = ((2 / 50) ** A_hellmann).astype("float16")
        del A_hellmann

    elif tech in ["WindOn", "WindOff"]:
        reg_ind = param["Ind_nz"]
        # A_cf
        if tech == "WindOn":
            paths_corr = paths["CORR_ON"]
        else:
            paths_corr = paths["CORR_OFF"]
        with rasterio.open(paths_corr) as src:
            w = src.read(1)
        rasterData["A_cf"] = np.flipud(w).astype("float16")
        if tech == "WindOff":
            rasterData["A_cf"] = rasterData["A_cf"][tuple(reg_ind)]
        del w
    return merraData, rasterData


def calc_FLH_solar(hours, args):
    """
    This function computes the full-load hours for all valid pixels specified in *ind_nz* in *param*.
    Due to parallel processing, most of the inputs are collected in the list *args*.

    :param hours: Filtered day hour ranks in a year (from 0 to 8759).
    :type hours: numpy array

    :param args: List of arguments:
        * *param* (dict): Dictionary including multiple parameters such as the status bar limit, the name of the region,
        and others for calculating the hourly capacity factors.
        * *tech* (str): Name of the technology.
        * *rasterData* (dict): Dictionary of numpy arrays containing land use types, Ross coefficients, albedo coefficients,
        and wind speed correction for every point in *reg_ind*.
        * *merraData* (dict): Dictionary of numpy arrays containing the weather data for every point in *reg_ind*.

    :type args: list
    :return FLH: Full-load hours over the year for the technology.
    :rtype: numpy array
    """
    # Decomposing the list args
    param = args[0]
    tech = args[1]
    rasterData = args[2]
    merraData = args[3]
    m_low = param["m_low"]
    n_low = param["n_low"]

    x = np.ones((m_low,n_low))
    ind = np.nonzero(x)

    FLH = np.zeros((m_low, n_low))
    status = 0
    for hour in hours:
        if hour <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            ul.display_progress(tech + " " + param["region_name"], [len(hours), status])

        if tech == "OpenFieldPV":
            CF = pm.calc_CF_solar(hour, ind, param, merraData, rasterData, tech)[0]
        elif tech == "RoofTopPV":
            CF = pm.calc_CF_solar(hour, ind, param, merraData, rasterData, tech)[1]
        elif tech == "CSP":
            CF = pm.calc_CF_solar(hour, ind, param, merraData, rasterData, tech)[2]

        # Aggregates CF to obtain the yearly FLH
        CF[np.isnan(CF)] = 0
        
        FLH[ind] = FLH[ind] + CF
    return FLH


def calc_FLH_windoff(hours, args):
    """
    This function computes the full-load hours for all valid pixels specified in *ind_nz* in *param*. Due to parallel processing,
    most of the inputs are collected in the list *args*.

    :param hours: Hour ranks in a year (from 0 to 8759).
    :type hours: numpy array
    :param args: List of arguments:
        * *param* (dict): Dictionary including multiple parameters such as the status bar limit, the name of the region, and
        others for calculating the hourly capacity factors.
        * *tech* (str): Name of the technology.
        * *rasterData* (dict): Dictionary of numpy arrays containing land use types, Ross coefficients, albedo coefficients,
        and wind speed correction for every point in *reg_ind*.
        * *merraData* (dict): Dictionary of numpy arrays containing the weather data for every point in *reg_ind*.

    :type args: list
    :return FLH: Full-load hours over the year for the technology.
    :rtype: numpy array
    """
    # Decomposing the tuple args
    param = args[0]
    tech = args[1]
    rasterData = args[2]
    merraData = args[3]
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
            ul.display_progress(tech + " " + param["region_name"], [len(hours), status])

        # Calculate hourly capacity factor
        CF = pm.calc_CF_windoff(hour, reg_ind, turbine, m_high, n_high, merraData, rasterData)

        # Aggregates CF to obtain the yearly FLH
        CF[np.isnan(CF)] = 0
        FLH = FLH + CF
    return FLH


# def calc_FLH_windon(row, args):
def calc_FLH_windon(param, tech, rasterData, merraData, GWA_array, b_xmin, b_xmax, b_ymin, b_ymax, x_gwa, y_gwa,
                    pixles, list_results):
    """
    This function computes the full-load hours for all valid pixels specified in *ind_nz* in *param*. Due to parallel processing,
    most of the inputs are collected in the list *args*.

    :param hours: Hour ranks in a year (from 0 to 8759).
    :type hours: numpy array
    :param args: List of arguments:
        * *param* (dict): Dictionary including multiple parameters such as the status bar limit, the name of the region, and
        others for calculating the hourly capacity factors.
        * *tech* (str): Name of the technology.
        * *rasterData* (dict): Dictionary of numpy arrays containing land use types, Ross coefficients, albedo coefficients,
        and wind speed correction for every point in *reg_ind*.
        * *merraData* (dict): Dictionary of numpy arrays containing the weather data for every point in *reg_ind*.

    :type args: list
    :return FLH: Full-load hours over the year for the technology.
    :rtype: numpy array
    """

    logger.debug('Process started - Pixles: (' + str(pixles) + ')')

    m_low = param["m_low"]
    n_low = param["n_low"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    turbine = param[tech]["technical"]

    reRaster = np.flipud(rasterData)  # ToDo: out of loop? # Why doing this?
    FLH_np = np.frombuffer(list_results, dtype=np.float32).reshape(m_high, n_high)

    for pixle in pixles:
        row, column = np.unravel_index(pixle, (m_low, n_low))  # Generate row and column out of the numbered position
        logger.debug('Pixel (' + str(row) + ',' + str(column) + ')')  # Print the recent computing pixle
        try:
            reMerra = redistribution_array(param, merraData[row, column, :], row, column, b_xmin[row, column], b_xmax[row, column], b_ymin[row, column], b_ymax[row, column], GWA_array, x_gwa, y_gwa)
            logger.debug('reMerra (' + str(row) + ',' + str(column) + ')')

            if np.sum(reMerra):  # FIXME: Why is reMerra zero?
                FLH_part = np.zeros([200, 250])
                for hour in np.arange(8760):
                    FLH_part += pm.calc_CF_windon(hour, turbine, reMerra,reRaster[row * 200:((row + 1) * 200), column * 250:((column + 1) * 250)])
                # # Does the same as above but needs more virtual memory (RAM)
                # hours = np.arange(8760)
                # CF = calc_CF_windon(hours, turbine, reMerra,
                #                     reRaster[row * 200:((row + 1) * 200), column * 250:((column + 1) * 250)])
                # FLH_part = np.nansum(CF, axis=2)
            else:
                FLH_part = np.zeros([200, 250])
                logger.debug('zeros(' + str(row) + ',' + str(column) + ')')

            rows_higherResolution = np.arange(row * 200, (row + 1) * 200)  # Extend the rows due to the higher resolution after redistribution Todo: use 'mm/n_high'?
            columns_higherResolution = np.arange(column * 250, (column + 1) * 250)  # same for columns
            logger.debug('Write on FLH_np (' + str(row) + ',' + str(column) + ')')
            FLH_np[np.ix_(rows_higherResolution,
                          columns_higherResolution)] = FLH_part  # Assign the computed FLH to the big FLH array
            logger.debug('Wrote on FLH_np (' + str(row) + ',' + str(column) + ')')
        except:
            traceback.print_exc()
            logger.error('Error on pixel (' + str(row) + ',' + str(column) + ')')

    logger.info('Done: ' + str(pixles))


def redistribution_array(param, merraData, i, j, xmin, xmax, ymin, ymax, GWA_array, x_gwa, y_gwa):
    """
    What does this function do?

    :param param:
    :type param:
    :param merraData:
    :type merraData:
    :param i:
    :type i:
    :param j:
    :type j:
    :param xmin:
    :type xmin:
    :param xmax:
    :type xmax:
    :param ymin:
    :type ymin:
    :param ymax:
    :type ymax:
    :param GWA_array:
    :type GWA_array:
    :param x_gwa:
    :type x_gwa:
    :param y_gwa:
    :type y_gwa:

    :return reMerra:
    :rtype: numpy array

    Aim:
        Increase the resolution of the MERRA wind data by using Global Wind Atlas data.
        For this reason, the low resolution MERRA data is redistributed by the energy distribution of the higher resolution Global Wind Atlas data

    Algorithm:
        1) Import wind data from Global Wind Atlas
        2) Select pixels that are within one MERRA pixle
        3) Convert from wind speed to energy
        4) Redistribute MERRA data

    ToDo: Where do the borders/limits come from? -> not shape file ?!?
    """

    logger.debug('start redistribution')

    # 1) Selection
    GWA_array_copy = GWA_array.astype(dtype=np.float32).copy()  # Create copy so that the origin array doesnt get changed
    # GWA_array_copy = GWA_array.astype(dtype=np.float16).copy()
    selection_index = (xmin <= x_gwa) & (x_gwa < xmax) & (ymin <= y_gwa) & (
            y_gwa < ymax)  # Determine the pixels that are inbetween the range
    GWA_array_copy[np.invert(selection_index)] = 0  # Set pixel not covered by the shapfile to zero
    value_num_cells = np.count_nonzero(GWA_array_copy)  # Number of non-zero pixels

    coordinates_nonzero_pixels_x, coordinates_nonzero_pixels_y = selection_index.nonzero()
    # Determine the coordinates of non-zero pixels in order to determine the outer rows and columns. Tuple(x,y)
    K_first = coordinates_nonzero_pixels_x[0]  # First x-coordinate of non-zero pixels
    L_first = coordinates_nonzero_pixels_y[0]  # First y-coordinate of non-zero pixels
    K_last = coordinates_nonzero_pixels_x[-1]  # Last x-coordinate of non-zero pixels
    L_last = coordinates_nonzero_pixels_y[-1]  # Last y-coordinate of non-zero pixels

    # 2) Redistribute MERRA data
    #reMerra = np.zeros([200, 250, 8760], dtype=np.float16)
    reMerra = np.zeros([200, 250, 8760], dtype = np.float32)
    if np.sum(GWA_array_copy):
        i_offset = 0
        j_offset = 0

        if K_last + 1 - K_first != 200:
            if i < param["m_low"] / 2:
                i_offset = 200 - (K_last + 1 - K_first)
        if L_last + 1 - L_first != 250:
            if j < param["n_low"] / 2:
                j_offset = 250 - (L_last + 1 - L_first)

        gwa_cut_wind = GWA_array_copy[K_first:K_last + 1, L_first:L_last + 1]  # Computation only for the nonzero values

        gwa_cut_energy = np.power(gwa_cut_wind, 3)
        merra_cut_energy_weighting = gwa_cut_energy / np.sum(gwa_cut_energy)  # Compute the weighting matrix how the energy is distributed within Global Wind Atlas data
        merra_cut_speed_weighting = np.cbrt(merra_cut_energy_weighting)  # Convert the weighting matrix from energy to wind speeds
        merra_cut_speed_redistributed = np.repeat(merra_cut_speed_weighting[..., None], 8760,
                                                  axis=2) * merraData * np.cbrt(
            value_num_cells)  # Expand the Merra time series of this pixle weighted by energy based distribution of Global Wind Atlas

        reMerra[i_offset:i_offset + K_last + 1 - K_first, j_offset:j_offset + L_last + 1 - L_first,:] = merra_cut_speed_redistributed

    return reMerra


def mask_potential_maps(paths, param, tech): #ToDo optimize no. of lines
    """
    This function first reads the rasters for land use, slope, bathymetry, and protected areas for the scope. Based on user-defined assumptions on
    their suitabilities, it generates a masking raster to exclude the unsuitable pixels. Both the mask itself
    and the masked potential rasters can be saved.

    :param paths: Dictionary of dictionaries containing user-defined parameters for masking, protected areas, and landuse.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the paths to the land use, protected areas, slope and bathymetry, in addition to output paths.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: The files for the mask and the masked FLH are saved as tif and mat files, along with their metadata json files.
    :rtype: None
    """
    logger.info("Start")
    mask = param[tech]["mask"]

    if tech in ["OpenFieldPV", "RoofTopPV", "CSP", "WindOn"]:
        with rasterio.open(paths["PA"]) as src:
            A_protect = src.read(1)
            A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified
        # Exclude protection categories that are not suitable
        A_suitability_pa = ul.changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)
        A_suitability_pa = (A_suitability_pa > 0).astype(int)
        A_lu = hdf5storage.read("LU", paths["LU"]).astype(int) # Landuse classes 0-16, to be reclassified
        # Exclude landuse types types that are not suitable
        A_suitability_lu = ul.changem(A_lu, mask["lu_suitability"], param["landuse"]["type"]).astype(float)
        A_suitability_lu = (A_suitability_lu > 0).astype(int)
        A_notWater = hdf5storage.read("BUFFER", paths["WATER_BUFFER"]).astype(int)
        A_notWetland = hdf5storage.read("BUFFER", paths["WETLAND_BUFFER"]).astype(int)
        A_notSnow = hdf5storage.read("BUFFER", paths["SNOW_BUFFER"]).astype(int)
        A_notBoarder = hdf5storage.read("BOARDERS", paths["BOARDERS"]).astype(int)
        with rasterio.open(paths["ROADS"]) as src:
            A_Roads = src.read(1)
            A_notRoads = (np.flipud(~A_Roads.astype(bool))).astype(int)
        with rasterio.open(paths["RAILS"]) as src:
            A_Rails = src.read(1)
            A_notRails = (np.flipud(~A_Rails.astype(bool))).astype(int)
        A_notMine = hdf5storage.read("BUFFER", paths["OSM_MINE_BUFFER"]).astype(int)
        A_NotLake = hdf5storage.read("BUFFER", paths["HYDROLAKES_BUFFER"]).astype(int)
        # Irrelevant parameters
        A_bathymetry = 1

    if tech in ["OpenFieldPV", "RoofTopPV", "CSP"]:
        A_notProtected = hdf5storage.read("BUFFER", paths["PV_PA_BUFFER"]).astype(int)
        with rasterio.open(paths["OSM_AREAS"]) as src:
            A_osma = src.read(1)
            A_com = A_osma == 1
            A_notCommercial = (np.flipud(~A_com.astype(bool))).astype(int)
            A_ind = A_osma == 2
            A_notIndustrial = (np.flipud(~A_ind.astype(bool))).astype(int)
            A_mil = A_osma == 4
            A_notMilitary = (np.flipud(~A_mil.astype(bool))).astype(int)
            A_rec = A_osma == 6
            A_notRecreation = (np.flipud(~A_rec.astype(bool))).astype(int)
        A_notPark = hdf5storage.read("BUFFER", paths["PV_OSM_PARK_BUFFER"]).astype(int)
        A_notSettlement = hdf5storage.read("BUFFER", paths["PV_WSF_BUFFER"]).astype(int)
        A_NotRiver = hdf5storage.read("BUFFER", paths["HYDRORIVERS_BUFFER"]).astype(int)
        A_notAirport = 1

    if tech in ["OpenFieldPV", "CSP"]:
        A_slope = hdf5storage.read("SLOPE", paths["SLOPE"])
        A_slope = (A_slope <= mask["slope"]).astype(int)

    if tech == "RoofTopPV":
        A_Settlement = hdf5storage.read("WSF", paths["WSF"]).astype(bool)
        A_Settlement = (np.flipud(A_Settlement)).astype(int)
        A_notSettlement = A_Settlement  # Just for using the same mask equation.
        # Irrelevant parameters
        A_slope = 1

    if tech == "WindOn":
        A_slope = hdf5storage.read("SLOPE", paths["SLOPE"])
        A_slope = (A_slope <= mask["slope"]).astype(int)
        A_notProtected = hdf5storage.read("BUFFER", paths["WINDON_PA_BUFFER"]).astype(int)
        A_notAirport = hdf5storage.read("BUFFER", paths["AIRPORTS"]).astype(int)
        A_notCommercial = hdf5storage.read("BUFFER", paths["OSM_COM_BUFFER"]).astype(int)
        A_notIndustrial = hdf5storage.read("BUFFER", paths["OSM_IND_BUFFER"]).astype(int)
        A_notMilitary = hdf5storage.read("BUFFER", paths["OSM_MIL_BUFFER"]).astype(int)
        A_notPark = hdf5storage.read("BUFFER", paths["WINDON_OSM_PARK_BUFFER"]).astype(int)
        A_notRecreation = hdf5storage.read("BUFFER", paths["OSM_REC_BUFFER"]).astype(int)
        A_notSettlement = hdf5storage.read("BUFFER", paths["WINDON_WSF_BUFFER"]).astype(int)
        with rasterio.open(paths["HYDRORIVERS"]) as src:
            A_Riv = src.read(1)
            A_NotRiver = (np.flipud(~A_Riv.astype(bool))).astype(int)
        # Irrelevant parameters
        A_bathymetry = 1

    if tech == "WindOff":
        A_suitability_lu = hdf5storage.read("EEZ", paths["EEZ"]).astype(int)
        with rasterio.open(paths["PA"]) as src:
            A_protect = src.read(1)
            A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified
        # Exclude protection categories that are not suitable
        A_suitability_pa = ul.changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)
        A_suitability_pa = (A_suitability_pa > 0).astype(int)
        A_bathymetry = hdf5storage.read("BATH", paths["BATH"]) # Bathymetry (depth) in meter
        A_bathymetry = (A_bathymetry >= mask["depth"]).astype(int) # (boolean)
        # Irrelevant parameters
        A_slope = 1
        A_notSettlement= 1
        A_notWater = 1
        A_notWetland = 1
        A_notSnow = 1
        A_notProtected = 1
        A_notAirport = 1
        A_notBoarder = 1
        A_notRoads = 1
        A_notCommercial = 1
        A_notIndustrial = 1
        A_notMine = 1
        A_notMilitary = 1
        A_notPark = 1
        A_notRecreation = 1

    # Masking matrix for the suitable sites (pixels)
    A_mask = (
            A_suitability_pa * A_suitability_lu * A_slope * A_bathymetry
            * A_notProtected * A_notAirport * A_notBoarder
            * A_notWater * A_notWetland * A_notSnow
            * A_notRoads * A_notRails * A_notCommercial * A_notIndustrial
            * A_notMine * A_notMilitary * A_notPark * A_notRecreation
            * A_notSettlement * A_NotLake * A_NotRiver
    ).astype(float)

    del A_suitability_lu, A_suitability_pa, A_slope, A_bathymetry
    del A_notProtected, A_notAirport, A_notBoarder
    del A_notWater, A_notWetland, A_notSnow
    del A_notRoads, A_notRails, A_notCommercial, A_notIndustrial
    del A_notMine, A_notMilitary, A_notPark, A_notRecreation
    del A_notSettlement, A_NotLake, A_NotRiver

    # Calculate masked FLH
    FLH = hdf5storage.read("FLH", paths[tech]["FLH"])
    FLH_mask = FLH * A_mask
    FLH_mask[FLH_mask == 0] = np.nan

    # Save HDF5 Files
    hdf5storage.writes({"A_mask": A_mask}, paths[tech]["mask"], store_python_metadata=True, matlab_compatible=True)
    logger.info("files saved: " + paths[tech]["mask"])
    hdf5storage.writes({"FLH_mask": FLH_mask}, paths[tech]["FLH_mask"], store_python_metadata=True,
                       matlab_compatible=True)
    logger.info("files saved: " + paths[tech]["FLH_mask"])

    ul.create_json(
        paths[tech]["mask"],
        param,
        ["author", "comment", tech, "region_name", "year", "GeoRef", "landuse", "protected_areas"],
        paths,
        ["subregions", "PA", "LU", "SLOPE", "BATH"],
    )

    # Save GEOTIFF files
    if param["savetiff_potentials"]:
        GeoRef = param["GeoRef"]
        sf.array2raster(ul.changeExt2tif(paths[tech]["mask"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"],
                     GeoRef["pixelHeight"], A_mask)
        logger.info("files saved: " + ul.changeExt2tif(paths[tech]["mask"]))

        # sf.array2raster(ul.changeExt2tif(paths[tech]["FLH_mask"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"],
        #              GeoRef["pixelHeight"], FLH_mask)
        # logger.info("files saved: " + ul.changeExt2tif(paths[tech]["FLH_mask"]))

    logger.debug("End")


def calc_gcr(Crd_all, m_high, n_high, res_desired, GCR):
    """
    This function creates a GCR weighting matrix for the desired geographic extent.
    The sizing of the PV system is conducted on a user-defined day for a shade-free exposure
    to the sun during a given number of hours.

    :param Crd_all: Desired geographic extent of the whole region (north, east, south, west).
    :type Crd_all: list
    :param m_high: Number of rows.
    :type m_high: int
    :param n_high: Number of columns.
    :type n_high: int
    :param res_desired: Resolution of the high resolution map.
    :type res_desired: list
    :param GCR: Dictionary that includes the user-defined day and the duration of the shade-free period.
    :type GCR: dict

    :return A_GCR: GCR weighting raster.
    :rtype: numpy array
    """
    # Vector of latitudes between (south) and (north), with resolution (res_should) degrees
    lat = np.arange((Crd_all[2] + res_desired[0] / 2), Crd_all[0], res_desired[0])[np.newaxis]
    lon = np.arange((Crd_all[3] + res_desired[1] / 2), Crd_all[1], res_desired[1])[np.newaxis]

    # Repeating for all longitudes/latitudes
    lat = ul.repmat(lat.transpose(), 1, n_high)
    lon = ul.repmat(lon, m_high, 1)

    # Solar time where shade-free exposure starts
    omegast = 12 - GCR["shadefree_period"] / 2

    # Calculation
    omega = 15 * (omegast - 12)  # Hour angle
    phi = abs(lat)  # Latitude angle

    beta = np.maximum(phi, 15)  # Tilt angle = latitude, but at least 15 degrees
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

    if Crd_all[2] > 0:
        day = GCR["day_north"]
        # Declination angle
        delta = ul.repmat(
            ul.arcsind(0.3978) * ul.sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * ul.sin(day * 2 * np.pi / 365.25 - 0.0489)),
            m_high, 1)

    if Crd_all[0] < 0:
        day = GCR["day_south"]
        # Declination angle
        delta = ul.repmat(
            ul.arcsind(0.3978) * ul.sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * ul.sin(day * 2 * np.pi / 365.25 - 0.0489)),
            m_high, 1)

    if (Crd_all[2] * Crd_all[0]) < 0:
        lat_pos = int(np.sum(lat >= 0, axis=0)[0])
        day = GCR["day_north"]
        # Declination angle
        delta_pos = ul.repmat(
            ul.arcsind(0.3978) * ul.sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * ul.sin(day * 2 * np.pi / 365.25 - 0.0489)),
            lat_pos, 1)

        lat_neg = int(np.sum(lat < 0, axis=0)[0])
        day = GCR["day_south"]
        # Declination angle
        delta_neg = ul.repmat(
            ul.arcsind(0.3978) * ul.sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * ul.sin(day * 2 * np.pi / 365.25 - 0.0489)),
            lat_neg, 1)
        delta = np.append(delta_neg, delta_pos, axis=0)

    # Elevation angle
    alpha = ul.arcsind(ul.sind(delta) * ul.sind(phi) + ul.cosd(delta) * ul.cosd(phi) * ul.cosd(omega))

    # Azimuth angle
    azi = ul.arccosd((ul.sind(delta) * ul.cosd(phi) - ul.cosd(delta) * ul.sind(phi) * ul.cosd(omega)) / ul.cosd(alpha))

    # The GCR
    A_GCR = 1 / (ul.cosd(beta) + np.abs(ul.cosd(azi)) * ul.sind(beta) / ul.tand(alpha))

    # Fix too large and too small values of GCR
    A_GCR[A_GCR < 0.2] = 0.2
    A_GCR[A_GCR > 0.9] = 0.9

    return A_GCR


def weight_potential_maps(paths, param, tech): #ToDo change variable names
    """
    This function weights the power potential by including assumptions on the power density and the available area.
    Therefore, it reads the rasters for land use and protected areas for the scope. Based on user-defined assumptions on
    their availabilities, it generates a weighting raster to exclude the unsuitable pixels. Both the weight itself
    and the weighted potential rasters can be saved.

    :param paths: Dictionary of dictionaries containing user-defined parameters for weighting, protected areas, and landuse.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the paths to the land use, protected areas, area, in addition to output paths.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: The files for the weight and the weighted FLH are saved as tif and mat files, along with their metadata json files.
    :rtype: None
    """
    logger.info("Start")
    weight = param[tech]["weight"]
    mask = param[tech]["mask"]
    Crd_all = param["Crd_all"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    res_desired = param["res_desired"]
    A_mask = hdf5storage.read("A_mask", paths[tech]["mask"])
    GeoRef = param["GeoRef"]

    if tech in ["OpenFieldPV", "RoofTopPV"]:
        # Ground Cover Ratio - defines spacing between PV arrays
        A_GCR = calc_gcr(Crd_all, m_high, n_high, res_desired, weight["GCR"])
    else:
        A_GCR = 1

    with rasterio.open(paths["PA"]) as src:
        A_protect = src.read(1)
        A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified

    # Calculate availability based on protection categories
    A_availability_pa = ul.changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)

    A_lu = hdf5storage.read("LU", paths["LU"]).astype(int)
    # with rasterio.open(paths["LU"]) as src:
    #     A_lu = src.read(1)
    #     A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified

    # Calculate availability based on landuse types
    A_availability_lu = ul.changem(A_lu, mask["lu_suitability"], param["landuse"]["type"]).astype(float)

    # Calculate availability
    A_availability = np.minimum(A_availability_pa, A_availability_lu)
    del A_availability_pa, A_availability_lu

    # Read available areas
    A_area = hdf5storage.read("A_area", paths["AREA"])

    # Weighting matrix for the power output (technical potential) in MWp
    A_weight = A_area * A_mask * A_GCR * weight["power_density"] * weight["f_performance"]
    if tech == "RoofTopPV":
        A_weight = A_weight * weight["suitable_roofs"]

    # Calculate weighted FLH in MWh
    FLH = hdf5storage.read("FLH", paths[tech]["FLH"])
    FLH_weight = FLH * A_weight

    # Save HDF5 Files
    hdf5storage.writes({"A_weight": A_weight}, paths[tech]["weight"], store_python_metadata=True,
                       matlab_compatible=True)
    logger.info("files saved: " + paths[tech]["weight"])
    hdf5storage.writes({"FLH_weight": FLH_weight}, paths[tech]["FLH_weight"], store_python_metadata=True,
                       matlab_compatible=True)
    logger.info("files saved: " + paths[tech]["FLH_weight"])
    ul.create_json(
        paths[tech]["weight"],
        param,
        ["author", "comment", tech, "region_name", "year", "GeoRef", "landuse", "protected_areas"],
        paths,
        ["subregions", "PA", "LU", "AREA"],
    )

    # Save GEOTIFF files
    # if param["savetiff_potentials"]:
    #     sf.array2raster(ul.changeExt2tif(paths[tech]["weight"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"],
    #                  GeoRef["pixelHeight"], A_weight)
    #     logger.info("files saved: " + ul.changeExt2tif(paths[tech]["weight"]))
    #
    #     sf.array2raster(ul.changeExt2tif(paths[tech]["FLH_weight"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"],
    #                  GeoRef["pixelHeight"], FLH_weight)
    #     logger.info("files saved: " + ul.changeExt2tif(paths[tech]["FLH_weight"]))
    logger.debug("End")


def sampled_sorting(Raster, sampling):
    """
    This function returns a list with a defined length of sorted values sampled from a numpy array.

    :param Raster: Input raster to be sorted.
    :type Raster: numpy array
    :param sampling: Number of values to be sampled from the raster, defines length of output list.
    :type sampling: int

    :return s: List of sorted values sampled from *Raster*.
    :rtype: list
    """
    # Flatten the raster and sort raster from highest to lowest
    Sorted_FLH = np.sort(Raster.flatten(order="F"))
    Sorted_FLH = np.flipud(Sorted_FLH)

    # Loop over list with sampling increment
    s = Sorted_FLH[0]  # Highest value
    for n in np.arange(sampling, len(Sorted_FLH), sampling):
        s = np.append(s, Sorted_FLH[n])
    s = np.append(s, Sorted_FLH[-1])  # Lowest value

    return s


def report_potentials(paths, param, tech):
    """
    This function reads the FLH files and the subregion shapefile, and creates a CSV file containing various statistics:

    * Available number of pixels, before and after masking
    * Available area in in km²
    * FLH mean, median, max, min values, before and after masking
    * FLH standard deviation after masking
    * Power Potential in GW, before and after weighting
    * Energy Potential in TWh in total, after weighting, and after masking and weighting
    * Sorted sample of FLH values for each region

    :param paths: Dictionary of dictionaries containing the paths to FLH, Masking, Weighting, and Area rasters.
    :type paths: dict
    :param param: Dictionary of dictionaries containing technology parameters and sampling parameters.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: The CSV files with the report and the sorted FLH are saved directly in the desired paths, along with the corresponding metadata in JSON files.
    :rtype: None
    """
    logger.info("Start")
    # read FLH, masking, area, and weighting matrix
    FLH = hdf5storage.read("FLH", paths[tech]["FLH"])
    A_mask = hdf5storage.read("A_mask", paths[tech]["mask"])
    A_weight = hdf5storage.read("A_weight", paths[tech]["weight"])
    A_area = hdf5storage.read("A_area", paths["AREA"])
    density = param[tech]["weight"]["power_density"]

    # Check if land or see
    if tech in ["OpenFieldPV", "RoofTopPV", "CSP", "WindOn"]:
        location = "land"
    elif tech in ["WindOff"]:
        location = "sea"

    # Initialize region masking parameters
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    nRegions = param["nRegions_land"]
    regions_shp = param["regions_land"]

    # Initialize regions list of sorted FLH, FLH_M, and FLH_W
    sorted_FLH_list = {}

    # Define sampling for sorted lists
    sampling = param["report_sampling"]

    # Initialize dataframe
    regions = pd.DataFrame(
        0,
        index=range(0, nRegions),
        columns=[
            "Region",
            "Available",
            "Available_Masked",
            "Available_Area_km2",
            "FLH_Mean",
            "FLH_Median",
            "FLH_Max",
            "FLH_Min",
            "FLH_Mean_Masked",
            "FLH_Median_Masked",
            "FLH_Max_Masked",
            "FLH_Min_Masked",
            "FLH_Std_Masked",
            "Power_Potential_GW",
            "Power_Potential_Weighted_GW",
            # "Power_Potential_Weighted_Masked_GW",
            "Energy_Potential_TWh",
            "Energy_Potential_Weighted_TWh",
            # "Energy_Potential_Weighted_Masked_TWh",
        ],
    )
    # Loop over each region
    # Display Progress
    status = 0
    ul.display_progress("Reporting ", (nRegions, status))
    for reg in range(0, nRegions):
        # Get name of region
        # regions.loc[reg, "Region"] = regions_shp.loc[reg]["NAME_SHORT"] + "_" + location
        regions.loc[reg, "Region"] = regions_shp.loc[reg]["GID_0"] + "_" + location

        # Compute region_mask
        A_region_extended = sf.calc_region(regions_shp.loc[reg], Crd_all, res_desired, GeoRef)

        # Sum available : available pixels
        available = np.sum(A_region_extended)
        regions.loc[reg, "Available"] = int(available)

        # Sum availabe_masked : available pixels after masking
        A_masked = A_region_extended * A_mask
        available_masked = np.nansum(A_masked)
        regions.loc[reg, "Available_Masked"] = int(available_masked)

        # Interrupt reporting of region if no available pixels
        if int(available_masked) == 0:
            regions.drop([reg], axis=0, inplace=True)
            continue

        # Interrupt reporting of region already reported (may occur due to discrepancy in borders)
        if regions.loc[reg, "Region"] in regions.loc[: reg - 1, "Region"].to_list():
            ind_prev = regions.loc[regions["Region"] == regions.loc[reg, "Region"]].index[0]
            if regions.loc[ind_prev, "Available_Masked"] > int(available_masked):
                regions.drop([reg], axis=0, inplace=True)
                continue
            else:
                regions.drop([ind_prev], axis=0, inplace=True)

        # Sum area: available area in km2
        A_area_region = A_region_extended * A_area
        Total_area = np.nansum(A_area_region) / (10 ** 6)
        regions.loc[reg, "Available_Area_km2"] = Total_area

        # Stats for FLH
        FLH_region = A_region_extended * FLH
        FLH_region[FLH_region == 0] = np.nan
        regions.loc[reg, "FLH_Mean"] = np.nanmean(FLH_region)
        regions.loc[reg, "FLH_Median"] = np.nanmedian(FLH_region)
        regions.loc[reg, "FLH_Max"] = np.nanmax(FLH_region)
        regions.loc[reg, "FLH_Min"] = np.nanmin(FLH_region)
        regions.loc[reg, "FLH_Std"] = np.nanstd(FLH_region)

        # Stats for FLH_masked
        FLH_region_masked = A_masked * FLH_region
        FLH_region_masked[FLH_region_masked == 0] = np.nan
        if int(np.nansum(FLH_region_masked)) == 0:
            continue
        regions.loc[reg, "FLH_Mean_Masked"] = np.nanmean(FLH_region_masked)
        regions.loc[reg, "FLH_Median_Masked"] = np.nanmedian(FLH_region_masked)
        regions.loc[reg, "FLH_Max_Masked"] = np.nanmax(FLH_region_masked)
        regions.loc[reg, "FLH_Min_Masked"] = np.nanmin(FLH_region_masked)
        regions.loc[reg, "FLH_Std_Masked"] = np.nanstd(FLH_region_masked)

        # Power Potential
        A_P_potential = A_area_region * density
        if tech == "RoofTopPV":
            A_P_potential = A_P_potential * param[tech]["weight"]["suitable_roofs"]
        power_potential = np.nansum(A_P_potential)
        regions.loc[reg, "Power_Potential_GW"] = power_potential / (10 ** 3)

        # Power Potential after weighting and masking
        A_P_W_potential = A_region_extended * A_weight
        power_potential_weighted = np.nansum(A_P_W_potential)
        regions.loc[reg, "Power_Potential_Weighted_GW"] = power_potential_weighted / (10 ** 3)

        # Power Potential after weighting and masking
        # A_P_W_M_potential = A_P_W_potential * A_masked
        # power_potential_weighted_masked = np.nansum(A_P_W_M_potential)
        # regions.loc[reg, "Power_Potential_Weighted_Masked_GW"] = power_potential_weighted_masked / (10 ** 3)

        # Energy Potential
        A_E_potential = A_P_potential * FLH_region
        energy_potential = np.nansum(A_E_potential)
        regions.loc[reg, "Energy_Potential_TWh"] = energy_potential / (10 ** 6)

        # Energy Potential after weighting and masking
        A_E_W_potential = FLH_region * A_weight
        energy_potential_weighted = np.nansum(A_E_W_potential)
        regions.loc[reg, "Energy_Potential_Weighted_TWh"] = energy_potential_weighted / (10 ** 6)

        # Energy Potential After weighting and masking
        # A_E_W_M_potential = A_E_W_potential * A_masked
        # energy_potential_weighted_masked = np.nansum(A_E_W_M_potential)
        # regions.loc[reg, "Energy_Potential_Weighted_Masked_TWh"] = energy_potential_weighted_masked / (10 ** 6)

        sort = {}
        # Sorted FLH Sampling
        sorted_sampled_FLH = sampled_sorting(FLH_region[~np.isnan(FLH_region)], sampling)
        sort["FLH"] = sorted_sampled_FLH

        # Sorted FLH Sampling after masking
        sorted_sampled_FLH_masked = sampled_sorting(FLH_region_masked[~np.isnan(FLH_region_masked)], sampling)
        sort["FLH_M"] = sorted_sampled_FLH_masked

        # Sorted FLH Sampling after masking and wieghting
        FLH_region_masked_weighted = FLH_region_masked * A_weight
        FLH_region_masked_weighted[FLH_region_masked_weighted == 0] = np.nan

        sorted_sampled_FLH_masked_weighted = sampled_sorting(
            FLH_region_masked_weighted[~np.isnan(FLH_region_masked_weighted)], sampling)
        sort["FLH_M_W"] = sorted_sampled_FLH_masked_weighted

        sorted_FLH_list[regions.loc[reg, "Region"]] = sort
        # Display Progress
        status += 1
        ul.display_progress("Reporting ", (nRegions, status))

    # Export the dataframe as CSV
    regions.to_csv(paths[tech]["Region_Stats"], sep=";", decimal=",", index=True)
    ul.create_json(
        paths[tech]["Region_Stats"],
        param,
        ["author", "comment", tech, "region_name", "subregions_name", "year", "res_desired", "Crd_all", "GeoRef"],
        paths,
        ["subregions", "subregions", "AREA", tech],
    )
    logger.info("files saved: " + paths[tech]["Region_Stats"])

    # Save sorted lists to mat file
    for reg in sorted_FLH_list.keys():
        hdf5storage.writes(
            {
                reg + "/FLH": sorted_FLH_list[reg]["FLH"],
                reg + "/FLH_masked": sorted_FLH_list[reg]["FLH_M"],
                reg + "/FLH_masked_weighted": sorted_FLH_list[reg]["FLH_M_W"],
            },
            paths[tech]["Sorted_FLH"],
            store_python_metadata=True,
            matlab_compatible=True,
        )
    ul.create_json(
        paths[tech]["Sorted_FLH"],
        param,
        ["author", "comment", tech, "region_name", "subregions_name", "year", "res_desired", "Crd_all", "GeoRef",
         "report_sampling"],
        paths,
        ["subregions", "subregions", "AREA", tech],
    )
    logger.info("files saved: " + paths[tech]["Sorted_FLH"])
    logger.debug("End")


def generate_biomass_production(paths, param, tech): #ToDo update to .mat files
    logger.info("Start")
    # Crd_all = param["Crd_all"]
    # GeoRef = param["GeoRef"]
    # res_desired = param["res_desired"]
    nRegions_land = param["nRegions_land"]
    # countries_shp = param["regions_land"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    biomass = param["Biomass"]

    # open land use raster to read the land type
    A_lu = hdf5storage.read("LU", paths["LU"]).astype(int)
    # with rasterio.open(paths["LU"]) as src:
    #     A_lu = src.read(1)
    # A_lu = np.flipud(A_lu).astype(int)
    A_lu_crop = (A_lu == 10) | (A_lu == 11) | (A_lu == 12) | (A_lu == 20) | (A_lu == 30)  # Agriculture pixels
    A_lu_veg = (A_lu == 40)  # Natural vegetation mosaic pixels
    A_lu_forest = (A_lu >= 50) | (A_lu <= 90)  # Forest pixels

    # open the protected areas raster to read the protected area type
    with rasterio.open(paths["PA"]) as src:
        A_pa = src.read(1)
    A_pa = np.flipud(A_pa).astype(int)
    A_Notprotected = A_pa == 0 # Unprotected pixels

    # Read the country codes
    IRENA_dict = pd.read_csv(paths["IRENA_dict"], sep=";", index_col=["Countries shapefile"],
                             usecols=["Countries shapefile", "IRENA"])

    # Define the result arrays for bioenergy and bioco2 for whole scope
    A_Bioenergy = np.zeros(A_lu.shape)
    A_Bioco2 = np.zeros(A_lu.shape)

    # for each country in the scope
    # for country in range(0, nRegions_land):
    with rasterio.open(paths["LAND"]) as src:
        A_country_area = src.read(1)
    A_country_area = np.flipud(A_country_area).astype(int)
        # A_country_area = sf.calc_region(countries_shp.loc[country], Crd_all, res_desired, GeoRef) # Country pixels
    country_name = IRENA_dict["IRENA"][param["country_code"]]

    # ==========Crop Residues Biomass potential==========#
    logger.info("Crop Residues Start")
    A_crop_country = np.multiply(A_lu_crop, A_country_area)
    n_crop_country = np.sum(A_crop_country)
    A_veg_country = np.multiply(A_lu_veg, A_country_area)
    n_veg_country = np.sum(A_veg_country)

        # if there are agriculture or vegetation mosaic pixels within country, read the annual crop production values from FAO database
    if n_crop_country or n_veg_country:
        production = pd.read_csv(paths["Biomass_Crops"], index_col=["Area"], usecols=["Area", "Item", "Value"])
        bio_energy = 0
        bio_co2 = 0
        if country_name in production.index:
            production_country = production[production.index == country_name]
            for crop in biomass["agriculture"]["crops"]:
                production_country_crop = production_country[production_country["Item"] == crop]
                if not production_country_crop.empty:
                    for residue in biomass["agriculture"]["residue"][crop]:
                        bio_energy = bio_energy \
                                         + (float(production_country_crop["Value"].iloc[0])
                                            * biomass["agriculture"]["rpr"][crop][residue]
                                            * biomass["agriculture"]["af"][crop][residue]
                                            * biomass["agriculture"]["lhv"][crop][residue])
                        bio_co2 = bio_co2 \
                                      + (float(production_country_crop["Value"].iloc[0])
                                         * biomass["agriculture"]["rpr"][crop][residue]
                                         * biomass["agriculture"]["af"][crop][residue]
                                         * biomass["agriculture"]["emission factor"])
        A_Bioenergy = A_Bioenergy\
                          + (A_crop_country * 2 * bio_energy / (2 * n_crop_country + n_veg_country))\
                          + (A_veg_country * bio_energy / (2 * n_crop_country + n_veg_country))
        A_Bioco2 = A_Bioco2\
                       + (A_crop_country * 2 * bio_co2 / (2 * n_crop_country + n_veg_country))\
                       + (A_veg_country * bio_co2 / (2 * n_crop_country + n_veg_country))
    logger.info("Crop Residues End")

        # ==========Forest wood Biomass potential==========#
    logger.info("Forest Wood Start")
    A_forest_country = np.multiply(A_lu_forest, A_country_area)
    A_forest_country = np.multiply(A_forest_country, A_Notprotected) # Only unprotected forest pixels
    n_forest_country = np.sum(A_forest_country)

    # if there are forest pixels within country, read the annual wood production values from FAO database
    if n_forest_country:
        production = pd.read_csv(paths["Biomass_Forestry"], index_col=["Area"], usecols=["Area", "Item", "Value"])
        if country_name in production.index:
            production_country = production[production.index == country_name]
            for wood in biomass["forest"]["woods"]:
                production_country_wood = production_country[production_country["Item"] == wood]
                if not production_country_wood.empty:
                    A_Bioenergy = A_Bioenergy \
                                      + (A_forest_country / n_forest_country
                                         * float(production_country_wood["Value"].iloc[0]) * biomass["forest"]["density"][wood]
                                         * biomass["forest"]["rpr"] * biomass["forest"]["af"] * biomass["forest"]["lhv"])

                    A_Bioco2 = A_Bioco2 \
                                   + (A_forest_country / n_forest_country
                                      * float(production_country_wood["Value"].iloc[0]) * biomass["forest"]["density"][wood]
                                      * biomass["forest"]["rpr"] * biomass["forest"]["af"] * biomass["forest"]["emission factor"])
    logger.info("Forest Wood End")

    # ==========Livestock Biomass potential==========#
    logger.info("Livestock Start")
    # Simple process for whole scope irrespective of number of countries in the scope
    # open the Land raster to read the valid pixels within scope
    with rasterio.open(paths["LAND"]) as src:
        A_scope_area = src.read(1)
    A_scope_area = np.flipud(A_scope_area).astype(int)
    param["Ind_nz"] = np.nonzero(A_scope_area)

    n_animal = 0  # number of animals
    for animal in param["Biomass"]["livestock"]["animal"]:
        # Extract Livestock density numbers
        # with rasterio.open(paths["LS"] + animal + ".tif") as src:
        #     A_LS_animal = src.read(1)
        # A_LS_animal = np.flipud(A_LS_animal)
        A_LS_animal = hdf5storage.read("LS", paths["LS"] + animal + ".mat")
        A_LS_animal = np.multiply(A_LS_animal, A_Notprotected)
        A_LS_animal = np.multiply(A_LS_animal, A_scope_area)

        A_Bioenergy = A_Bioenergy\
                      + (A_LS_animal
                         * biomass["livestock"]["rpr"][n_animal]
                         * biomass["livestock"]["af"][n_animal]
                         * biomass["livestock"]["lhv"][n_animal])
        A_Bioco2 = A_Bioco2\
                   + (A_LS_animal
                     * biomass["livestock"]["rpr"][n_animal]
                     * biomass["livestock"]["af"][n_animal]
                     * biomass["livestock"]["emission factor"])
        n_animal = n_animal + 1  # call for next animal type
    logger.info("Livestock End")

    energy_map = np.full((m_high, n_high), np.nan)
    energy_map[param["Ind_nz"]] = A_Bioenergy[param["Ind_nz"]]
    co2_map = np.full((m_high, n_high), np.nan)
    co2_map[param["Ind_nz"]] = A_Bioco2[param["Ind_nz"]]
    hdf5storage.writes({"BIOMASS_ENERGY": energy_map}, paths[tech]["BIOMASS_ENERGY"], store_python_metadata=True,
                       matlab_compatible=True)
    ul.create_json(
        paths[tech]["BIOMASS_ENERGY"],
        param,
        ["author", "comment", tech, "region_name", "year", "res_desired"],
        paths,
        ["LS","LU","PA"],
    )
    logger.info("\nfiles saved: " + paths[tech]["BIOMASS_ENERGY"])

    # Save GEOTIFF files
    if param["savetiff_potentials"]:
        GeoRef = param["GeoRef"]
        sf.array2raster(ul.changeExt2tif(paths[tech]["BIOMASS_ENERGY"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"],
                     GeoRef["pixelHeight"], energy_map)
        logger.info("files saved:" + ul.changeExt2tif(paths[tech]["BIOMASS_ENERGY"]))

    hdf5storage.writes({"BIOMASS_CO2": co2_map}, paths[tech]["BIOMASS_CO2"], store_python_metadata=True,
                       matlab_compatible=True)
    ul.create_json(
        paths[tech]["BIOMASS_CO2"],
        param,
        ["author", "comment", tech, "region_name", "year", "res_desired"],
        paths,
        ["LS","LU","PA"],
    )
    logger.info("files saved: " + paths[tech]["BIOMASS_CO2"])

    # Save GEOTIFF files
    # if param["savetiff_potentials"]:
    #     GeoRef = param["GeoRef"]
    #     sf.array2raster(ul.changeExt2tif(paths[tech]["BIOMASS_CO2"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"],
    #                  GeoRef["pixelHeight"], co2_map)
    #     logger.info("files saved:" + ul.changeExt2tif(paths[tech]["BIOMASS_CO2"]))

    logger.debug("End")


def report_biomass_potentials(paths, param, tech):
    logger.info("Start")
    # read FLH, masking, area, and weighting matrix
    A_Bioenergy = hdf5storage.read("BIOMASS_ENERGY", paths[tech]["BIOMASS_ENERGY"])
    A_Bioco2 = hdf5storage.read("BIOMASS_CO2", paths[tech]["BIOMASS_CO2"])

    # Initialize region masking parameters
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    nRegions_land = param["nRegions_land"]
    countries_shp = param["regions_land"]

    # Initialize dataframe
    countries = pd.DataFrame(
        0,
        index=range(0, nRegions_land),
        columns=[
            "Region",
            "Bio_Energy_Potential_TWh",
            "Bio_CO2_emissions_million_tons"
        ],
    )
    # Loop over each region
    # Display Progress
    status = 0
    ul.display_progress("Reporting ", (nRegions_land, status))
    for reg in range(0, nRegions_land):
        # Get name of Country
        # countries.loc[reg, "Region"] = countries_shp.loc[reg]["NAME_SHORT"]
        countries.loc[reg, "Region"] = countries_shp.loc[reg]["GID_0"]

        # Compute region_mask
        A_country_extended = sf.calc_region(countries_shp.loc[reg], Crd_all, res_desired, GeoRef)

        # Energy Potential
        A_E_potential = A_Bioenergy * A_country_extended
        energy_potential = np.nansum(A_E_potential)
        countries.loc[reg, "Bio_Energy_Potential_TWh"] = energy_potential / (10 ** 6)

        # Interrupt reporting of region if no available energy
        if int(energy_potential) == 0:
            countries.drop([reg], axis=0, inplace=True)
            continue
        # Interrupt reporting of region already reported (may occur due to discrepancy in borders)
        if countries.loc[reg, "Region"] in countries.loc[: reg - 1, "Region"].to_list():
            ind_prev = countries.loc[countries["Region"] == countries.loc[reg, "Region"]].index[0]
            if countries.loc[ind_prev, "Bio_Energy_Potential_TWh"] > int(energy_potential):
                countries.drop([reg], axis=0, inplace=True)
                continue
            else:
                countries.drop([ind_prev], axis=0, inplace=True)

        A_CO2_emissions = A_Bioco2 * A_country_extended
        CO2_emissions = np.nansum(A_CO2_emissions)
        countries.loc[reg, "Bio_CO2_emissions_million_tons"] = CO2_emissions / (10 ** 9)

        # Display Progress
        status += 1
        ul.display_progress("Reporting ", (nRegions_land, status))

    # Export the dataframe as CSV
    countries.to_csv(paths[tech]["Region_Stats"], sep=";", decimal=",", index=True)
    ul.create_json(
        paths[tech]["Region_Stats"],
        param,
        ["author", "comment", tech, "region_name", "year", "res_desired", "Crd_all", "GeoRef"],
        paths,
        ["AREA", tech],
    )
    logger.info("files saved: " + paths[tech]["Region_Stats"])
    logger.debug("End")