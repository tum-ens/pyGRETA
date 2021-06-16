from lib.physical_models import calc_CF_solar, calc_CF_windon, calc_CF_windoff
from lib.spatial_functions import *
import multiprocessing as mp


def calculate_full_load_hours(paths, param, tech):
    """
    This function calculates the yearly FLH for a technology for all valid pixels in a spatial scope. Valid pixels are land pixels
    for WindOn, PV and CSP, and sea pixels for WindOff. The FLH values are calculated by summing up hourly capacity factors.

    :param paths: Dictionary of dictionaries containing the paths to the input weather data, land, sea and land use rasters, and correction rasters.
    :type paths: dict
    :param param: Dictionary of dictionaries containing the spatial scope, and technology and computation parameters.
    :type param: dict
    :param tech: Technology under study.
    :type tech: str

    :return: The raster of FLH potential is saved as mat and tif files, along with the json metadata file.
    :rtype: None
    """
    timecheck("Start")
    print("Region: " + param["region_name"])

    if tech in ["WindOn", "WindOff"]:
        print("\n" + tech + " - HUB_HEIGHTS: " + str(param[tech]["technical"]["hub_height"]))
    elif tech in ["PV"] and "orientation" in param["PV"]["technical"].keys():
        print("\n" + tech + " - Orientation: " + str(param[tech]["technical"]["orientation"]))

    nproc = param["nproc"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    m_low = param["m_low"]
    n_low = param["n_low"]
    CPU_limit = np.full((1, nproc), param["CPU_limit"])
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    countries_shp = param["regions_land"]
    nRegions = param["nRegions_sub"]
    regions_shp = param["regions_sub"]
    res_weather = param["res_weather"]
    Crd_all = param["Crd_all"]
    Ind = ind_merra(Crd_all, Crd_all, res_weather)[0]
    
    if tech == "WindOff":
        with rasterio.open(paths["EEZ"]) as src:
            w = src.read(1)
    else:
        with rasterio.open(paths["LAND"]) as src:
            w = src.read(1)
    
    param["Ind_nz"] = np.nonzero(np.flipud(w))
  
    A_country_area = np.zeros(w.shape)
    for reg in range(0, nRegions):
        r = calc_region(regions_shp.loc[reg], Crd_all, res_desired, GeoRef)  
        A_country_area = A_country_area + r
    A_country_area = np.flipud(A_country_area)
    del w    
    
    # Obtain weather and correction matrices
    merraData, rasterData = get_merra_raster_data(paths, param, tech)
        
    if tech in ["PV", "CSP"]:

        day_filter = np.nonzero(merraData["CLEARNESS"][Ind[2] - 1 : Ind[0], Ind[3] - 1 : Ind[1], :].sum(axis=(0, 1)))
        list_hours = np.arange(0, 8760)
        if nproc == 1:
            param["status_bar_limit"] = list_hours[-1]
            results = calc_FLH_solar(list_hours[day_filter], [param, tech, rasterData, merraData])
        else:
            list_hours = np.array_split(list_hours[day_filter], nproc)
            param["status_bar_limit"] = list_hours[0][-1]
            results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
                calc_FLH_solar, product(list_hours, [[param, tech, rasterData, merraData]])
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
        create_json(
            paths[tech]["FLH"],
            param,
            ["author", "comment", tech, "region_name", "subregions_name", "year", "res_desired", "res_weather"],
            paths,
            ["spatial_scope"],
        )
        print("\nfiles saved: " + paths[tech]["FLH"])

        # Save GEOTIFF files
        if param["savetiff"]:
            GeoRef = param["GeoRef"]
            array2raster(changeExt2tif(paths[tech]["FLH"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], FLH)
            print("files saved:" + changeExt2tif(paths[tech]["FLH"]))
    
    elif tech in ["WindOff"]:
        list_hours = np.array_split(np.arange(0, 8760), nproc)
        param["status_bar_limit"] = list_hours[0][-1]
        results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
            calc_FLH_windoff, product(list_hours, [[param, tech, rasterData, merraData]])
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
        create_json(
            paths[tech]["FLH"],
            param,
            ["author", "comment", tech, "region_name", "subregions_name", "year", "res_desired", "res_weather"],
            paths,
            ["spatial_scope"],
        )
        print("\nfiles saved: " + paths[tech]["FLH"])

        # Save GEOTIFF files
        if param["savetiff"]:
            GeoRef = param["GeoRef"]
            array2raster(changeExt2tif(paths[tech]["FLH"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], FLH)
            print("files saved:" + changeExt2tif(paths[tech]["FLH"]))
    
    elif tech in ["WindOn"]:

        merraData = merraData["W50M"][::-1,:,:]
        rasterData = rasterData["A_cf"]
        
        with rasterio.open(paths["GWA_global"]) as src:
            GWA_array = src.read(1)
        GWA_array = np.power(GWA_array, 3)
        where_are_NaNs = np.isnan(GWA_array)
        GWA_array[where_are_NaNs] = 0

        # param["status_bar_limit"] = list_rows_splitted[0][-1] # ToDo: Wo we need this?

        # -------------------------------------------------------------------------------------
        # Multiprocessing by Patrick 20210615
        # FIXME Import what?
        b_xmin = hdf5storage.read("MERRA_XMIN", paths["MERRA_XMIN"])  # ToDo: move into calculate_fullload_hours? # ToDo: Put outside of function?
        b_xmax = hdf5storage.read("MERRA_XMAX", paths["MERRA_XMAX"])
        b_ymin = hdf5storage.read("MERRA_YMIN", paths["MERRA_YMIN"])
        b_ymax = hdf5storage.read("MERRA_YMAX", paths["MERRA_YMAX"])

        # FIXME Import What?
        x_gwa = hdf5storage.read("GWA_X", paths["GWA_X"])
        y_gwa = hdf5storage.read("GWA_Y", paths["GWA_Y"])


        list_pixles = np.arange(n_low*m_low)     # All pixles within MERRA data
        # list_pixles = np.arange(100, 110) # debuging

        list_pixles_splitted = np.array_split(list_pixles, nproc)   # Splitted list acording to the number of parrallel processes
        print('# of processes: ' + str(len(list_pixles_splitted)))
        # print(list_rows_splitted)

        processes = []  # Store all single process of multiprocessing
        list_results = mp.Array('f', m_high*n_high, lock=False)
        FLH = np.frombuffer(list_results, dtype=np.float32).reshape(m_high, n_high)
        # np.copyto(FLH, np.full((m_high, n_high), np.nan))
        FLH[:] = np.nan


        for pixles in list_pixles_splitted:  # Run the 'calc_FLH_windon' for each of the splitted rows
            p = mp.Process(target=calc_FLH_windon, args=(param, tech, paths, rasterData, merraData, GWA_array, b_xmin, b_xmax, b_ymin, b_ymax, x_gwa, y_gwa, pixles, list_results))
            processes.append(p)

        for p in processes:
            p.start()   # Start all single processes

        for p in processes:
            p.join()    # Wait until all processes are finished
        print('All processes finished')

        # ------------------------------------------------------------------------------------

        # FLH[A_country_area==0] = float("nan")python
        FLH = np.flipud(FLH)
        hdf5storage.writes({"FLH": FLH}, paths[tech]["FLH"], store_python_metadata=True, matlab_compatible=True)

        print("\nfiles saved: " + paths[tech]["FLH"])

        GeoRef = param["GeoRef"]
        array2raster(changeExt2tif(paths[tech]["FLH"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], FLH)
        print("files saved:" + changeExt2tif(paths[tech]["FLH"]))

    timecheck("End")


def get_merra_raster_data(paths, param, tech):
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
    if tech in ["PV", "CSP"]:

        # Other weather Data
        # Clearness index - stored variable CLEARNESS
        merraData["CLEARNESS"] = hdf5storage.read("CLEARNESS", paths["CLEARNESS"])
        # Temperature 2m above the ground - stored variable T2M
        merraData["T2M"] = hdf5storage.read("T2M", paths["T2M"])

        # Calculate A matrices correction
        # A_lu
        with rasterio.open(paths["LU"]) as src:
            w = src.read(1)
        rasterData["A_lu"] = np.flipud(w)
        # A_Ross (Temperature coefficients for heating losses)
        rasterData["A_Ross"] = changem(rasterData["A_lu"], param["landuse"]["Ross_coeff"], param["landuse"]["type"]).astype("float16")
        # A_albedo (Reflectivity coefficients)
        rasterData["A_albedo"] = changem(rasterData["A_lu"], param["landuse"]["albedo"], param["landuse"]["type"]).astype("float16")
        # A_WS_Coef wind Speed at 2m above the ground
        A_hellmann = changem(rasterData["A_lu"], landuse["hellmann"], landuse["type"])
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
    reg_ind = param["Ind_nz"]

    FLH = np.zeros(len(reg_ind[0]))
    status = 0
    for hour in hours:
        if hour <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            display_progress(tech + " " + param["region_name"], [len(hours), status])

        if tech == "PV":
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech)[0]
        elif tech == "CSP":
            CF = calc_CF_solar(hour, reg_ind, param, merraData, rasterData, tech)[1]

        # Aggregates CF to obtain the yearly FLH
        CF[np.isnan(CF)] = 0
        FLH = FLH + CF
        
        #print(str(hour)+"_"+str(sum(FLH)))

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
            display_progress(tech + " " + param["region_name"], [len(hours), status])

        # Calculate hourly capacity factor
        CF = calc_CF_windoff(hour, reg_ind, turbine, m_high, n_high, merraData, rasterData)

        # Aggregates CF to obtain the yearly FLH
        CF[np.isnan(CF)] = 0
        FLH = FLH + CF
    return FLH


# def calc_FLH_windon(row, args):
def calc_FLH_windon(param, tech, paths, rasterData, merraData, GWA_array, b_xmin, b_xmax, b_ymin, b_ymax, x_gwa, y_gwa, pixles, list_results):
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
    # param = args[0]
    # tech = args[1]
    # paths = args[2]
    # rasterData = args[3]
    # merraData = args[4]
    # GWA_array = args[5]
    # pixles = args[6]


    
    # reg_ind = param["Ind_nz"]
    m_low = param["m_low"]
    n_low = param["n_low"]
    # res_weather = param["res_weather"]
    # res_desired = param["res_desired"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    turbine = param[tech]["technical"]

    # FLH = np.zeros([len(rows)*200, n_high])
    # print(type(FLH))
    FLH =[]
    rows_higherResolution = []
    columns_higherResolution = []
    reRaster = np.flipud(rasterData)  # ToDo: out of loop? # Why doing this?

    # rows = [10]   # In case for debuging
    # olumns = [10] # In case for debuging
    # columns = range(n_low)      # Loop over all columns

    FLH_np = np.frombuffer(list_results, dtype=np.float32).reshape(m_high, n_high)

    for pixle in pixles:
        row, column = np.unravel_index(pixle, (m_low, n_low))   # Generate row and column out of the numbered position
        print('Pixle (' + str(row) + ',' + str(column) + ')')   # Print the recent computing pixle

        FLH_part = np.zeros([200, 250])
        reMerra = redistribution_array(param, paths, merraData[row,column,:],row,column,b_xmin[row,column],b_xmax[row,column],b_ymin[row,column],b_ymax[row,column],GWA_array,x_gwa,y_gwa)

        if np.sum(reMerra):
            for hour in range(8760):
                CF = calc_CF_windon(hour, turbine, reMerra, reRaster[row*200:((row+1)*200),column*250:((column+1)*250)])
                CF[np.isnan(CF)] = 0
                FLH_part = FLH_part + CF    # Aggregates CF to obtain the yearly FLH

            # rows_higherResolution = np.arange(row * 200, (row + 1) * 200)  # Extend the rows due to the higher resolution after redistribution Todo: use 'mm/n_high'?
            # columns_higherResolution = np.arange(column * 250, (column + 1) * 250) # Extend the columns due to the higher resolution after redistribution

            # FLH_np[rows_higherResolution, columns_higherResolution] = FLH_part
            FLH_np[row * 200:((row + 1) * 200), column * 250:((column + 1) * 250)] = FLH_part



def redistribution_array(param, paths, merraData, i, j, xmin, xmax, ymin, ymax, GWA_array, x_gwa, y_gwa):

    """
    What does this function do?

    :param param:
    :type param:
    :param paths:
    :type paths:
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
    :param Num_pix:
    :type Num_pix:

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

    # 1) Selection
    GWA_array_copy = GWA_array.copy()  # Create copy so that the origin array doesnt get changed
    selection_index = (xmin <= x_gwa) & (x_gwa < xmax) & (ymin <= y_gwa) & (y_gwa < ymax)   # Determine the pixels that are inbetween the range
    GWA_array_copy[np.invert(selection_index)] = 0  # Set pixel not covered by the shapfile to zero
    value_num_cells = np.count_nonzero(GWA_array_copy)       # Number of non-zero pixels

    coordinates_nonzero_pixels_x, coordinates_nonzero_pixels_y = selection_index.nonzero() # Determine the coordinates of non-zero pixels in order to determine the outer rows and columns. Tuple(x,y)
    K_first = coordinates_nonzero_pixels_x[0]      # First x-coordinate of non-zero pixels
    L_first = coordinates_nonzero_pixels_y[0]      # First y-coordinate of non-zero pixels
    K_last = coordinates_nonzero_pixels_x[-1]      # Last x-coordinate of non-zero pixels
    L_last = coordinates_nonzero_pixels_y[-1]      # Last y-coordinate of non-zero pixels

    # 2) Redistribute MERRA data
    reMerra = np.zeros([200, 250, 8760])

    if np.sum(GWA_array_copy):
        i_offset = 0
        j_offset = 0

        if K_last+1-K_first != 200:
            if i < param["m_low"]/2:
                i_offset = 200 - (K_last+1-K_first)
        if L_last+1-L_first != 250:
            if j < param["n_low"]/2:
                j_offset = 250 - (L_last+1-L_first)

        gwa_cut_energy = GWA_array_copy[K_first:K_last+1, L_first:L_last+1]     # Computation only for the nonzero values
        merra_cut_energy_weighting = gwa_cut_energy / np.sum(gwa_cut_energy)    # Compute the weighting matrix how the energy is distributed within Global Wind Atlas data
        merra_cut_speed_weighting = np.cbrt(merra_cut_energy_weighting)     # Convert the weighting matrix from energy to wind speeds
        merra_cut_speed_redistributed = np.repeat(merra_cut_speed_weighting[..., None], 8760, axis=2) * merraData * np.cbrt(value_num_cells)    # Expand the Merra time series of this pixle weighted by energy based distribution of Global Wind Atlas
        reMerra[i_offset:i_offset + K_last + 1 - K_first, j_offset:j_offset + L_last + 1 - L_first, :] = merra_cut_speed_redistributed

    return reMerra


def mask_potential_maps(paths, param, tech):
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
    timecheck("Start")
    mask = param[tech]["mask"]

    if tech in ["PV", "CSP"]:
        with rasterio.open(paths["PA"]) as src:
            A_protect = src.read(1)
            A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified
        # Exclude protection categories that are not suitable
        A_suitability_pa = changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)
        A_suitability_pa = (A_suitability_pa > 0).astype(int)
        with rasterio.open(paths["LU"]) as src:
            A_lu = src.read(1)
            A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified
        # Exclude landuse types types that are not suitable
        A_suitability_lu = changem(A_lu, mask["lu_suitability"], param["landuse"]["type"]).astype(float)
        A_suitability_lu = (A_suitability_lu > 0).astype(int)
        with rasterio.open(paths["SLOPE"]) as src:
            A_slope = src.read(1)
            A_slope = np.flipud(A_slope)  # Slope in percentage
            A_slope = (A_slope <= mask["slope"]).astype(int)
        with rasterio.open(paths["WATER_BUFFER"]) as src:
            A_notWater = src.read(1)
            A_notWater = (np.flipud(A_notWater)).astype(int)
        with rasterio.open(paths["WETLAND_BUFFER"]) as src:
            A_notWetland = src.read(1)
            A_notWetland = (np.flipud(A_notWetland)).astype(int)
        with rasterio.open(paths["SNOW_BUFFER"]) as src:
            A_notSnow = src.read(1)
            A_notSnow = (np.flipud(A_notSnow)).astype(int)
        with rasterio.open(paths["PV_PA_BUFFER"]) as src:
            A_notProtected = src.read(1)
            A_notProtected = (np.flipud(A_notProtected)).astype(int)
        with rasterio.open(paths["BOARDERS"]) as src:
            A_notBoarder = src.read(1)
            A_notBoarder = (np.flipud(A_notBoarder)).astype(int)
        # Irrelevant parameters
        A_notPopulated = 1
        A_bathymetry = 1
        A_notAirport = 1


    if tech == "WindOn":
        with rasterio.open(paths["PA"]) as src:
            A_protect = src.read(1)
            A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified
        # Exclude protection categories that are not suitable
        A_suitability_pa = changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)
        A_suitability_pa = (A_suitability_pa > 0).astype(int)
        with rasterio.open(paths["LU"]) as src:
            A_lu = src.read(1)
            A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified
        # Exclude landuse types types that are not suitable
        A_suitability_lu = changem(A_lu, mask["lu_suitability"], param["landuse"]["type"]).astype(float)
        A_suitability_lu = (A_suitability_lu > 0).astype(int)
        with rasterio.open(paths["SLOPE"]) as src:
            A_slope = src.read(1)
            A_slope = np.flipud(A_slope)  # Slope in percentage
            A_slope = (A_slope <= mask["slope"]).astype(int)
        with rasterio.open(paths["POP_BUFFER"]) as src:
            A_notPopulated = src.read(1)
            A_notPopulated = (np.flipud(A_notPopulated)).astype(int)  # Is 1 for not populated areas
        with rasterio.open(paths["WATER_BUFFER"]) as src:
            A_notWater = src.read(1)
            A_notWater = (np.flipud(A_notWater)).astype(int)
        with rasterio.open(paths["WETLAND_BUFFER"]) as src:
            A_notWetland = src.read(1)
            A_notWetland = (np.flipud(A_notWetland)).astype(int)
        with rasterio.open(paths["SNOW_BUFFER"]) as src:
            A_notSnow = src.read(1)
            A_notSnow = (np.flipud(A_notSnow)).astype(int)
        with rasterio.open(paths["WINDON_PA_BUFFER"]) as src:
            A_notProtected = src.read(1)
            A_notProtected = (np.flipud(A_notProtected)).astype(int)
        with rasterio.open(paths["AIRPORTS"]) as src:
            A_notAirport = src.read(1)
            A_notAirport = (np.flipud(A_notAirport)).astype(int)
        with rasterio.open(paths["BOARDERS"]) as src:
            A_notBoarder = src.read(1)
            A_notBoarder = (np.flipud(A_notBoarder)).astype(int)
        # Irrelevant parameters
        A_bathymetry = 1

    if tech == "WindOff":
        with rasterio.open(paths["EEZ"]) as src:
            A_suitability_lu = src.read(1)
            A_suitability_lu = np.flipud(A_suitability_lu).astype(int)
        with rasterio.open(paths["PA"]) as src:
            A_protect = src.read(1)
            A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified
        # Exclude protection categories that are not suitable
        A_suitability_pa = changem(A_protect, mask["pa_suitability"], param["protected_areas"]["type"]).astype(float)
        A_suitability_pa = (A_suitability_pa > 0).astype(int)
        with rasterio.open(paths["BATH"]) as src:
            A_bathymetry = src.read(1)
            A_bathymetry = np.flipud(A_bathymetry)  # Bathymetry (depth) in meter
            A_bathymetry = (A_bathymetry >= mask["depth"]).astype(int)  # (boolean)
        # Irrelevant parameters
        A_slope = 1
        A_notPopulated = 1
        A_notWater = 1
        A_notWetland = 1 
        A_notSnow = 1
        A_notProtected = 1
        A_notAirport = 1
        A_notBoarder = 1

    # Masking matrix for the suitable sites (pixels)
    A_mask = (A_suitability_pa * A_suitability_lu * A_slope * A_notPopulated * A_notWater * A_notWetland * A_notSnow * A_notProtected * A_notAirport * A_notBoarder * A_bathymetry).astype(float)

    del A_suitability_lu, A_suitability_pa, A_slope, A_notPopulated, A_notWater, A_notWetland, A_notSnow, A_bathymetry

    # Calculate masked FLH
    FLH = hdf5storage.read("FLH", paths[tech]["FLH"])
    FLH_mask = FLH * A_mask
    FLH_mask[FLH_mask == 0] = np.nan

    # Save HDF5 Files
    hdf5storage.writes({"A_mask": A_mask}, paths[tech]["mask"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths[tech]["mask"])
    hdf5storage.writes({"FLH_mask": FLH_mask}, paths[tech]["FLH_mask"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths[tech]["FLH_mask"])

    create_json(
        paths[tech]["mask"],
        param,
        ["author", "comment", tech, "region_name", "year", "GeoRef", "landuse", "protected_areas"],
        paths,
        ["spatial_scope", "PA", "LU", "SLOPE", "BATH"],
    )

    # Save GEOTIFF files
    if param["savetiff"]:
        GeoRef = param["GeoRef"]
        array2raster(changeExt2tif(paths[tech]["mask"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_mask)
        print("files saved: " + changeExt2tif(paths[tech]["mask"]))

        array2raster(changeExt2tif(paths[tech]["FLH_mask"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], FLH_mask)
        print("files saved: " + changeExt2tif(paths[tech]["FLH_mask"]))

    timecheck("End")


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
    lat = repmat(lat.transpose(), 1, n_high)
    lon = repmat(lon, m_high, 1)

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
    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat, range_lon)] - 35) / 65 * 45 + 35  # Europe
    range_lat = np.logical_and(lat >= 20, lat < 65)
    range_lon = np.logical_and(lon >= 75, lon < 140)
    beta[np.logical_and(range_lat, range_lon)] = (beta[np.logical_and(range_lat, range_lon)] - 20) / 65 * 60 + 20  # Asia/China

    if Crd_all[2] > 0:
        day = GCR["day_north"]
        # Declination angle
        delta = repmat(arcsind(0.3978) * sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), m_high, 1)

    if Crd_all[0] < 0:
        day = GCR["day_south"]
        # Declination angle
        delta = repmat(arcsind(0.3978) * sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), m_high, 1)

    if (Crd_all[2] * Crd_all[0]) < 0:
        lat_pos = int(np.sum(lat >= 0, axis=0)[0])
        day = GCR["day_north"]
        # Declination angle
        delta_pos = repmat(arcsind(0.3978) * sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), lat_pos, 1)

        lat_neg = int(np.sum(lat < 0, axis=0)[0])
        day = GCR["day_south"]
        # Declination angle
        delta_neg = repmat(arcsind(0.3978) * sin(day * 2 * np.pi / 365.25 - 1.400 + 0.0355 * sin(day * 2 * np.pi / 365.25 - 0.0489)), lat_neg, 1)
        delta = np.append(delta_neg, delta_pos, axis=0)

    # Elevation angle
    alpha = arcsind(sind(delta) * sind(phi) + cosd(delta) * cosd(phi) * cosd(omega))

    # Azimuth angle
    azi = arccosd((sind(delta) * cosd(phi) - cosd(delta) * sind(phi) * cosd(omega)) / cosd(alpha))

    # The GCR
    A_GCR = 1 / (cosd(beta) + np.abs(cosd(azi)) * sind(beta) / tand(alpha))

    # Fix too large and too small values of GCR
    A_GCR[A_GCR < 0.2] = 0.2
    A_GCR[A_GCR > 0.9] = 0.9

    return A_GCR


def weight_potential_maps(paths, param, tech):
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
    timecheck("Start")
    weight = param[tech]["weight"]
    Crd_all = param["Crd_all"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    res_desired = param["res_desired"]
    GeoRef = param["GeoRef"]

    if tech == "PV":
        # Ground Cover Ratio - defines spacing between PV arrays
        A_GCR = calc_gcr(Crd_all, m_high, n_high, res_desired, weight["GCR"])
    else:
        A_GCR = 1

    with rasterio.open(paths["PA"]) as src:
        A_protect = src.read(1)
        A_protect = np.flipud(A_protect).astype(int)  # Protection categories 0-10, to be classified

    # Calculate availability based on protection categories
    A_availability_pa = changem(A_protect, weight["pa_availability"], param["protected_areas"]["type"]).astype(float)

    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
        A_lu = np.flipud(A_lu).astype(int)  # Landuse classes 0-16, to be reclassified

    # Calculate availability based on landuse types
    A_availability_lu = changem(A_lu, weight["lu_availability"], param["landuse"]["type"]).astype(float)

    # Calculate availability
    A_availability = np.minimum(A_availability_pa, A_availability_lu)
    del A_availability_pa, A_availability_lu

    # Read available areas
    A_area = hdf5storage.read("A_area", paths["AREA"])

    # Weighting matrix for the power output (technical potential) in MWp
    A_weight = A_area * A_availability * A_GCR * weight["power_density"] * weight["f_performance"]
    
    
    # Calculate weighted FLH in MWh
    FLH = hdf5storage.read("FLH", paths[tech]["FLH"])
    FLH_weight = FLH * A_weight
    


    # Save HDF5 Files
    hdf5storage.writes({"A_weight": A_weight}, paths[tech]["weight"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths[tech]["weight"])
    hdf5storage.writes({"FLH_weight": FLH_weight}, paths[tech]["FLH_weight"], store_python_metadata=True, matlab_compatible=True)
    print("files saved: " + paths[tech]["FLH_weight"])
    create_json(
        paths[tech]["weight"],
        param,
        ["author", "comment", tech, "region_name", "year", "GeoRef", "landuse", "protected_areas"],
        paths,
        ["spatial_scope", "PA", "LU", "AREA"],
    )

    # Save GEOTIFF files
    if param["savetiff"]:
        array2raster(changeExt2tif(paths[tech]["weight"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_weight)
        print("files saved: " + changeExt2tif(paths[tech]["weight"]))

        array2raster(changeExt2tif(paths[tech]["FLH_weight"]), GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], FLH_weight)
        print("files saved: " + changeExt2tif(paths[tech]["FLH_weight"]))
    timecheck("End")


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
    timecheck("Start")
    # read FLH, masking, area, and weighting matrix
    FLH = hdf5storage.read("FLH", paths[tech]["FLH"])
    A_mask = hdf5storage.read("A_mask", paths[tech]["mask"])
    A_weight = hdf5storage.read("A_weight", paths[tech]["weight"])
    A_area = hdf5storage.read("A_area", paths["AREA"])
    density = param[tech]["weight"]["power_density"]

    # Check if land or see
    if tech in ["PV", "CSP", "WindOn"]:
        location = "land"
    elif tech in ["WindOff"]:
        location = "sea"

    # Initialize region masking parameters
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    nRegions = param["nRegions_sub"]
    regions_shp = param["regions_sub"]

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
            "Power_Potential_Weighted_Masked_GW",
            "Energy_Potential_TWh",
            "Energy_Potential_Weighted_TWh",
            "Energy_Potential_Weighted_Masked_TWh",
        ],
    )
    # Loop over each region
    # Display Progress
    status = 0
    display_progress("Reporting ", (nRegions, status))
    for reg in range(0, nRegions):
        # Get name of region
        regions.loc[reg, "Region"] = regions_shp.loc[reg]["NAME_SHORT"] + "_" + location

        # Compute region_mask
        A_region_extended = calc_region(regions_shp.loc[reg], Crd_all, res_desired, GeoRef)

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
        power_potential = np.nansum(A_P_potential)
        regions.loc[reg, "Power_Potential_GW"] = power_potential / (10 ** 3)

        # Power Potential after weighting
        A_P_W_potential = A_region_extended * A_weight
        power_potential_weighted = np.nansum(A_P_W_potential)
        regions.loc[reg, "Power_Potential_Weighted_GW"] = power_potential_weighted / (10 ** 3)
        
        # Power Potential after weighting and masking
        A_P_W_M_potential = A_P_W_potential * A_masked
        power_potential_weighted_masked = np.nansum(A_P_W_M_potential)
        regions.loc[reg, "Power_Potential_Weighted_Masked_GW"] = power_potential_weighted_masked / (10 ** 3)

        # Energy Potential
        A_E_potential = A_P_potential * FLH_region
        energy_potential = np.nansum(A_E_potential)
        regions.loc[reg, "Energy_Potential_TWh"] = energy_potential / (10 ** 6)

        # Energy Potential after weighting
        A_E_W_potential = FLH_region * A_weight
        energy_potential_weighted = np.nansum(A_E_W_potential)
        regions.loc[reg, "Energy_Potential_Weighted_TWh"] = energy_potential_weighted / (10 ** 6)

        # Energy Potential After weighting and masking
        A_E_W_M_potential = A_E_W_potential * A_masked
        energy_potential_weighted_masked = np.nansum(A_E_W_M_potential)
        regions.loc[reg, "Energy_Potential_Weighted_Masked_TWh"] = energy_potential_weighted_masked / (10 ** 6)

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

        sorted_sampled_FLH_masked_weighted = sampled_sorting(FLH_region_masked_weighted[~np.isnan(FLH_region_masked_weighted)], sampling)
        sort["FLH_M_W"] = sorted_sampled_FLH_masked_weighted

        sorted_FLH_list[regions.loc[reg, "Region"]] = sort
        # Display Progress
        status += 1
        display_progress("Reporting ", (nRegions, status))

    # Export the dataframe as CSV
    regions.to_csv(paths[tech]["Region_Stats"], sep=";", decimal=",", index=True)
    create_json(
        paths[tech]["Region_Stats"],
        param,
        ["author", "comment", tech, "region_name", "subregions_name", "year", "res_desired", "Crd_all", "GeoRef"],
        paths,
        ["spatial_scope", "subregions", "AREA", tech],
    )
    print("files saved: " + paths[tech]["Region_Stats"])

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
    create_json(
        paths[tech]["Sorted_FLH"],
        param,
        ["author", "comment", tech, "region_name", "subregions_name", "year", "res_desired", "Crd_all", "GeoRef", "report_sampling"],
        paths,
        ["spatial_scope", "subregions", "AREA", tech],
    )
    print("files saved: " + paths[tech]["Sorted_FLH"])
    timecheck("End")


def generate_biomass_production(paths, param):
    
    timecheck("Start")
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    nRegions = param["nRegions_sub"]
    regions_shp = param["regions_sub"]
    nproc = param["nproc"]
    CPU_limit = np.full((1, nproc), param["CPU_limit"])
    
    with rasterio.open(paths["LU"]) as src:
        A_lu = src.read(1)
    A_lu = np.flipud(A_lu).astype(int)
    
    with rasterio.open(paths["PA"]) as src:
        A_protect = src.read(1)
    A_protect = np.flipud(A_protect).astype(int)
    
    with rasterio.open(paths["SUB"]) as src:
        A_country_area = src.read(1)
    A_country_area = np.flipud(A_country_area).astype(int)
    
    #Pixels within the shapefile or country!
    # A_country_area = np.zeros(A_lu.shape)
    # list_regions = np.array_split(np.arange(0, nRegions), nproc)
    # param["status_bar_limit"] = list_regions[0][-1]
    # results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
        # country_region, product(list_regions,[[param, paths, regions_shp, Crd_all, res_desired, GeoRef]])
    # )
    # # Collecting results
    # if nproc > 1:
        # for p in range(len(results)):
            # A_country_area = A_country_area + results[p]
    # else:
        # A_country_area = results[p]
   
    #Extract un-protected areas
    A_Notprotected = np.zeros(A_protect.shape) 
    val_include = [0,5,6,7,8,9,10]
    for i in val_include:
        A_i = A_protect == i
        A_Notprotected = A_Notprotected + A_i
    
    A_Bioenergy = np.zeros(A_lu.shape)
    A_Bioco2 = np.zeros(A_lu.shape)
    
    #Crop Residues potential
    list_regions = np.array_split(np.arange(0, nRegions), nproc)
    param["status_bar_limit"] = list_regions[0][-1]
    results = Pool(processes=nproc, initializer=limit_cpu, initargs=CPU_limit).starmap(
        crop_residues_potential, product(list_regions,[[param, paths, A_lu]])
    )
    # Collecting results
    if nproc > 1:
        for p in range(len(results)):
            A_Bioenergy = A_Bioenergy + results[p]["A_Bioenergy"]
            A_Bioco2 = A_Bioco2 + results[p]["A_Bioco2"]
    else:
        A_Bioenergy = results[p]["A_Bioenergy"]
        A_Bioco2 = results[p]["A_Bioco2"]
    
    # #Extract forest areas
    # A_lu_forest = np.zeros(A_lu.shape)
    # val_forest = [1,2,3,4,5]
    # for j in val_forest:
        # A_j = A_lu == j
        # A_lu_forest = A_lu_forest + A_j
    # A_lu_forest = np.multiply(A_lu_forest, A_country_area)
    # A_lu_forest = np.multiply(A_lu_forest, A_Notprotected)
    # n_lu_forest = np.sum(A_lu_forest)
    # print (n_lu_forest)    
    
    # IRENA_dict = pd.read_csv(paths["IRENA_dict"], sep=";",index_col = ["Countries shapefile"],usecols=["Countries shapefile","IRENA"])
    # country_name = IRENA_dict["IRENA"][param["country_code"]]
    # print (country_name)
    # #Add Bio potential from forest wood
    # forest = param["Biomass"]["forest"]
    # if n_lu_forest:   
        # production_wood = pd.read_csv(paths["Biomass_Forestry"], index_col = ["Area"],usecols=["Area","Item","Value"])
        # production_wood = production_wood[production_wood.index==country_name].to_numpy()
        # print (production_wood)
        # for row in range(len(production_wood[:,0])):
            # if production_wood[row,0] == "Wood fuel, coniferous":
                # wood_coniferous = production_wood[row,1]
                # print (wood_coniferous)
            # else:
                # wood_coniferous = 0
            # if production_wood[row,0] == "Wood fuel, non-coniferous":
                # wood_nonconiferous = production_wood[row,1]
                # print (wood_nonconiferous)
            # else:
                # wood_nonconiferous = 0
                    
            # #Forest Bio Energy
            # energy_lu_forest = A_lu_forest * (wood_coniferous*forest["density, coniferous"] + wood_nonconiferous*forest["density, non-coniferous"])*forest["rpr"]*forest["af"]*forest["lhv"] / n_lu_forest
            # A_Bioenergy = A_Bioenergy + energy_lu_forest
        
            # #Forest Bio CO2
            # co2_lu_forest = A_lu_forest * (wood_coniferous*forest["density, coniferous"] + wood_nonconiferous*forest["density, non-coniferous"])*forest["rpr"]*forest["af"]*forest["emission factor"] / n_lu_forest
            # A_Bioco2 = A_Bioco2 + co2_lu_forest
    
    
    # #Add Bio potential from Livestock
    # n_animal = 0
    # for animal in param["Biomass"]["livestock"]["animal"]:
        # #Extract Livestock numbers
        # with rasterio.open(paths["LS"]+animal+".tif") as src:
            # A_LS_animal = src.read(1)
        # A_LS_animal = np.flipud(A_LS_animal)
        # A_LS_animal = np.multiply(A_LS_animal,A_Notprotected)
        # A_LS_animal = np.multiply(A_LS_animal,A_country_area)
        
        # #Livestock Bio Energy
        # energy_ls_animal = A_LS_animal * param["Biomass"]["livestock"]["rpr"][n_animal] * param["Biomass"]["livestock"]["af"][n_animal] * param["Biomass"]["livestock"]["lhv"][n_animal]
        # A_Bioenergy = A_Bioenergy + energy_ls_animal
        
        # #Livestock Bio CO2
        # co2_ls_animal = A_LS_animal * param["Biomass"]["livestock"]["rpr"][n_animal] * param["Biomass"]["livestock"]["af"][n_animal] * param["Biomass"]["livestock"]["emission factor"]
        # A_Bioco2 = A_Bioco2 + co2_ls_animal
        
        # n_animal = n_animal+1
        
    array2raster(paths["BIOMASS_ENERGY"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_Bioenergy)
    print("files saved: " + paths["BIOMASS_ENERGY"])
    create_json(paths["BIOMASS_ENERGY"], param, ["region_name", "landuse", "Biomass", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "BIOMASS_ENERGY"])
    
    array2raster(paths["BIOMASS_CO2"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_Bioco2)
    print("files saved: " + paths["BIOMASS_CO2"])
    create_json(paths["BIOMASS_CO2"], param, ["region_name", "landuse", "Biomass", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "BIOMASS_CO2"])
    
    timecheck("End") 

def country_region(regions, args):
    param = args[0]
    paths = args[1]
    regions_shp = args[2]
    Crd_all = args[3]
    res_desired = args[4]
    GeoRef = args[5]
    
    m_high = param["m_high"]
    n_high = param["n_high"]
    
    A_country_area = np.zeros([m_high,n_high])
    status = 0
    for reg in regions:
        if reg <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            display_progress("Country Area Pixels " + param["region_name"], [len(regions), status])
        
        A_region_extended = calc_region(regions_shp.loc[reg], Crd_all, res_desired, GeoRef)
        A_country_area = A_country_area + A_region_extended
    
    return A_country_area

def crop_residues_potential(regions, args):    
    
    param = args[0]
    paths = args[1]
    A_lu = args[2]
    
    
    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    nRegions = param["nRegions_sub"]
    regions_shp = param["regions_sub"]
    m_high = param["m_high"]
    n_high = param["n_high"]
    
    
    A_Bioenergy = np.zeros([m_high,n_high])
    A_Bioco2 = np.zeros([m_high,n_high])
    A_country_area = np.zeros([m_high,n_high])
    
    status = 0
    for reg in regions:
        if reg <= param["status_bar_limit"]:
            # Show progress of the simulation
            status = status + 1
            display_progress("Crop Residues " + param["region_name"], [len(regions), status])
        
        A_region_extended = calc_region(regions_shp.loc[reg], Crd_all, res_desired, GeoRef)
        A_country_area = A_country_area + A_region_extended    
        
        A_lu_crop = A_lu == param["landuse"]["type_croplands"]
        A_lu_veg = A_lu == param["landuse"]["type_vegetation"]
        
        A_lu_crop = np.multiply(A_lu_crop,A_region_extended)
        n_lu_crop = np.sum(A_lu_crop)
        #print (n_lu_crop)
        A_lu_veg = np.multiply(A_lu_veg,A_region_extended)
        n_lu_veg = np.sum(A_lu_veg)
        #print (n_lu_veg)

        if n_lu_crop or n_lu_veg:        
            region_code = regions_shp.loc[reg]["GID_2"][:-2]
            if len(region_code)==11:
                region_code = region_code[0:3]+region_code[4:7]+region_code[8:]
            elif len(region_code) == 10:
                if region_code[7] == ".":
                    region_code = region_code[0:3]+region_code[4:7]+"0"+region_code[8:]
                else:
                    region_code = region_code[0:3]+"0"+region_code[4:6]+region_code[7:]
            elif len(region_code) == 9:
                if region_code[7] == ".":
                    region_code = region_code[0:3]+region_code[4:7]+"00"+region_code[8:]
                elif region_code[6] == ".":
                    region_code = region_code[0:3]+"0"+region_code[4:6]+"0"+region_code[7:]
                else:
                    region_code = region_code[0:3]+"00"+region_code[4]+region_code[6:]
            elif len(region_code) == 8:
                if region_code[6] == ".":
                    region_code = region_code[0:3]+"0"+region_code[4:6]+"00"+region_code[7:]
                else:
                    region_code = region_code[0:3]+"00"+region_code[4]+"0"+region_code[6:]
            else:
                region_code = region_code[0:3]+"00"+region_code[4]+"00"+region_code[6:]
            #print (region_code)
            bio_energy = 0
            bio_co2 = 0
            for crop in param["Biomass"]["agriculture"]["crops"]:
                production_crop = pd.read_csv(paths["Biomass_Crops"]+crop+".csv", index_col = ["agromap_area_2"],usecols=["agromap_area_2","Production"])
                production_crop = production_crop.dropna()
                if region_code in production_crop.index:
                    production_crop = production_crop[production_crop.index==region_code]
                    for residue in param["Biomass"]["agriculture"]["residue"][crop]:
                        if isinstance(production_crop["Production"].iloc[0], str):
                            bio_energy = bio_energy + float(production_crop["Production"].iloc[0].replace(',',''))* param["Biomass"]["agro_rpr"][crop][residue]*param["Biomass"]["agro_af"][crop][residue]*param["Biomass"]["agro_lhv"][crop][residue]
                            bio_co2 = bio_co2 + float(production_crop["Production"].iloc[0].replace(',',''))* param["Biomass"]["agro_rpr"][crop][residue]*param["Biomass"]["agro_af"][crop][residue]*param["Biomass"]["agro_emission factor"]
                        else:
                            bio_energy = bio_energy + float(production_crop["Production"].iloc[0])* param["Biomass"]["agro_rpr"][crop][residue]*param["Biomass"]["agro_af"][crop][residue]*param["Biomass"]["agro_lhv"][crop][residue]
                            bio_co2 = bio_co2 + float(production_crop["Production"].iloc[0])* param["Biomass"]["agro_rpr"][crop][residue]*param["Biomass"]["agro_af"][crop][residue]*param["Biomass"]["agro_emission factor"]

            energy_lu_crop = A_lu_crop * 2 * bio_energy / (2*n_lu_crop+n_lu_veg)
            energy_lu_veg = A_lu_veg * bio_energy / (2*n_lu_crop+n_lu_veg)
            
            co2_lu_crop = A_lu_crop * 2 * bio_co2 / (2*n_lu_crop+n_lu_veg)
            co2_lu_veg = A_lu_veg * bio_co2 / (2*n_lu_crop+n_lu_veg)
        
            energy_reg = energy_lu_crop + energy_lu_veg
            A_Bioenergy = A_Bioenergy + energy_reg

            co2_reg = co2_lu_crop + co2_lu_veg
            A_Bioco2 = A_Bioco2 + co2_reg
    Bioresults = dict()
    Bioresults["A_Bioenergy"] = A_Bioenergy
    Bioresults["A_Bioco2"] = A_Bioco2
    
    return Bioresults
    
def club_biomass(paths,param):

    Crd_all = param["Crd_all"]
    GeoRef = param["GeoRef"]
    res_desired = param["res_desired"]
    
    Bioenergy = rasterio.open(paths["BIOMASS_ENERGY"])
    A_Bioenergy = Bioenergy.read(1)
    A_Bioenergy_rows, A_Bioenergy_cols = A_Bioenergy.shape

    West = rasterio.open(paths["West"])
    A_West = West.read(1)
    A_West_rows, A_West_cols = A_West.shape
    ind_West = ind_exact_points([West.xy(0,0,offset='ul')[1],West.xy(0,0,offset='ul')[0]], Crd_all, res_desired)
    A_Bioenergy[A_Bioenergy_rows-ind_West[0]:A_Bioenergy_rows-ind_West[0]+A_West_rows, ind_West[1]:ind_West[1]+A_West_cols] = A_Bioenergy[A_Bioenergy_rows-ind_West[0]:A_Bioenergy_rows-ind_West[0]+A_West_rows, ind_West[1]:ind_West[1]+A_West_cols] + A_West
    
    North = rasterio.open(paths["North"])
    A_North = North.read(1)
    A_North_rows, A_North_cols = A_North.shape
    ind_North = ind_exact_points([North.xy(0,0,offset='ul')[1],North.xy(0,0,offset='ul')[0]], Crd_all, res_desired)
    A_Bioenergy[A_Bioenergy_rows-ind_North[0]:A_Bioenergy_rows-ind_North[0]+A_North_rows, ind_North[1]:ind_North[1]+A_North_cols] = A_Bioenergy[A_Bioenergy_rows-ind_North[0]:A_Bioenergy_rows-ind_North[0]+A_North_rows, ind_North[1]:ind_North[1]+A_North_cols] + A_North
    
    East = rasterio.open(paths["East"])
    A_East = East.read(1)
    A_East_rows, A_East_cols = A_East.shape
    ind_East = ind_exact_points([East.xy(0,0,offset='ul')[1],East.xy(0,0,offset='ul')[0]], Crd_all, res_desired)
    A_Bioenergy[A_Bioenergy_rows-ind_East[0]:A_Bioenergy_rows-ind_East[0]+A_East_rows, ind_East[1]:ind_East[1]+A_East_cols] = A_Bioenergy[A_Bioenergy_rows-ind_East[0]:A_Bioenergy_rows-ind_East[0]+A_East_rows, ind_East[1]:ind_East[1]+A_East_cols] + A_East
    
    South = rasterio.open(paths["South"])
    A_South = South.read(1)
    A_South_rows, A_South_cols = A_South.shape
    ind_South = ind_exact_points([South.xy(0,0,offset='ul')[1],South.xy(0,0,offset='ul')[0]], Crd_all, res_desired)
    A_Bioenergy[A_Bioenergy_rows-ind_South[0]:A_Bioenergy_rows-ind_South[0]+A_South_rows, ind_South[1]:ind_South[1]+A_South_cols] = A_Bioenergy[A_Bioenergy_rows-ind_South[0]:A_Bioenergy_rows-ind_South[0]+A_South_rows, ind_South[1]:ind_South[1]+A_South_cols] + A_South
    
    A_Bioenergy = np.flipud(A_Bioenergy)    
    array2raster(paths["CLUB_CO2"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_Bioenergy)
    print("files saved: " + paths["CLUB_CO2"])
    create_json(paths["CLUB_CO2"], param, ["region_name", "landuse", "Biomass", "Crd_all", "res_desired", "GeoRef"], paths, ["LU", "BIOMASS_ENERGY"])
    