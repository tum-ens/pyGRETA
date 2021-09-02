from . import spatial_functions as sf
from . import util as ul
from .log import logger
import numpy as np
import pandas as pd
import hdf5storage
import os
import rasterio
from scipy.ndimage import generic_filter


def clean_weather_data(paths, param):
    """
    This function detects data outliers in the weather input .mat files. An outlier is a data point, for which
    the absolute value of the difference between the yearly average value and the mean of the direct neighbors
    (Moore neighborhood) is higher than a user-defined threshold *MERRA_correction_factor*. It replaces the hourly values
    with the hourly values of the mean of the neighbors, and overwrites the original .mat file.

    :param paths: Dictionary including the path to the weather .mat files.
    :type paths: dict
    :param param: Dictionary including the threshold value *MERRA_correction_factor*.
    :type param: dict

    :return: The file weather .mat files are overwritten after the correction.
    :rtype: None
    """
    logger.info('Start')
    for p in ["W50M", "CLEARNESS", "T2M"]:

        # Read Weather Data
        weather = hdf5storage.read(p, paths[p])
        mean = np.mean(weather, 2)

        # Set convolution mask
        kernel = np.ones((3, 3))
        kernel[1, 1] = 0

        # Compute average Convolution
        neighbors = generic_filter(mean, np.nanmean, footprint=kernel, mode="constant", cval=np.NaN)
        ratio = mean / neighbors

        # Extract over threshold Points
        points = np.where(abs(ratio - np.mean(ratio)) > param["MERRA_correction_factor"][p])

        # Correct points hourly
        for t in range(weather.shape[2]):
            weather[points[0], points[1], t] = weather[points[0], points[1], t] / ratio[points[0], points[1]]

        # Save corrected Wind
        hdf5storage.writes({p: weather}, paths[p], store_python_metadata=True, matlab_compatible=True)
    logger.debug("End")


def generate_wind_correction(paths, param):
    """
    This function creates a matrix of correction factors for onshore and/or offshore wind.
    There are different types of correction:

    * Gradient correction: Adjusts for the hub height of the wind turbines, based on the Hellmann coefficients of each land use type.
      This correction applies always.
    * Resolution correction: Performs a redistribution of wind speed when increasing the resolution based on land use types, while ensuring that
      the average of each MERRA-2 cell at 50m is still the same. This correction is optional, and is activated if *res_correction* is 1.
      If not activated, the same value from the low resolution is repeated.
    * Topographic/Orographic correction: Takes into account the elevation of the terrain, because MERRA-2 usually underestimates
      the wind speed in mountains. This correction is optional, uses data from the Global Wind Atlas for all countries in the scope,
      and is activated only for onshore wind if *topo_correction* is 1

    :param paths: Dictionary of dictionaries containing the paths to the land, land use, and topography rasters, and to the output files CORR_ON and CORR_OFF.
    :type paths: dict
    :param param: Dictionary of dictionaries containing user-preferences regarding the wind correction, landuse, hub height, weather and desired resolutions.
    :type param: dict

    :return: The rasters for wind correction CORR_ON and/or CORR_OFF are saved directly in the user-defined paths, along with their metadata in JSON files.
    :rtype: None
    """
    if os.path.isfile(paths["CORR_ON"]) and os.path.isfile(paths["CORR_OFF"]):
        logger.info('Skip')    # Skip generation if files are already there

    else:
        logger.info("Start")
        GeoRef = param["GeoRef"]
        landuse = param["landuse"]
        A_lu = hdf5storage.read("LU", paths["LU"]).astype(int)
        # with rasterio.open(paths["LU"]) as src:
        #     A_lu = np.flipud(src.read(1)).astype(int)
        A_hellmann = ul.changem(A_lu, landuse["hellmann"], landuse["type"]).astype(float)

        # Onshore height correction
        if "WindOn" in param["technology"]:
            turbine_height_on = param["WindOn"]["technical"]["hub_height"]
            A_cf_on = (turbine_height_on / 50) ** A_hellmann
            with rasterio.open(paths["LAND"]) as src:
                A_land = np.flipud(src.read(1)).astype(int)
            A_cf_on = A_cf_on * A_land
            del A_land
            sf.array2raster(paths["CORR_ON"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_cf_on)
            ul.create_json(paths["CORR_ON"], param, ["region_name", "year", "WindOn", "landuse", "res_weather", "res_desired"], paths, ["LAND"])
            logger.info("files saved: " + paths["CORR_ON"])

        # Offshore height correction
        if "WindOff" in param["technology"]:
            turbine_height_off = param["WindOff"]["technical"]["hub_height"]
            A_cf_off = (turbine_height_off / 50) ** A_hellmann
            del A_hellmann
            A_eez = hdf5storage.read("EEZ", paths["EEZ"]).astype(int)
            # with rasterio.open(paths["EEZ"]) as src:
            #     A_eez = np.flipud(src.read(1)).astype(int)
            A_cf_off = A_cf_off * A_eez

            sf.array2raster(paths["CORR_OFF"], GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"], A_cf_off)
            ul.create_json(paths["CORR_OFF"], param, ["region_name", "year", "WindOff", "landuse", "res_weather", "res_desired"], paths, ["EEZ"])
            logger.info("files saved: " + paths["CORR_OFF"])
        logger.debug("End")


def clean_IRENA_summary(paths, param):
    """
    This function reads the IRENA database, format the output for selected regions and computes the FLH based on the
    installed capacity and yearly energy production. The results are saved in CSV file.

    :param param: Dictionary of dictionaries containing list of subregions, and year.
    :type param: dict
    :param paths: Dictionary of dictionaries containing the paths to the IRENA country name dictionary, and IRENA database.
    :type paths: dict

    :return: The CSV file containing the summary of IRENA data for the countries within the scope is saved directly in the desired path, along with the corresponding metadata in a JSON file.
    :rtype: None
    """
    year = str(param["year"])
    filter_countries = param["regions_land"]["GID_0"].to_list()
    IRENA_dict = pd.read_csv(paths["IRENA_dict"], sep=";", index_col=0)
    IRENA_dict = IRENA_dict["Countries shapefile"].to_dict()
    IRENA = pd.read_csv(paths["IRENA"], skiprows=7, sep=";", index_col=False, usecols=[0, 1, 2, 3])
    for i in IRENA.index:
        if pd.isnull(IRENA.loc[i, "Country/area"]):
            IRENA.loc[i, "Country/area"] = IRENA.loc[i - 1, "Country/area"]
        if pd.isnull(IRENA.loc[i, "Technology"]):
            IRENA.loc[i, "Technology"] = IRENA.loc[i - 1, "Technology"]

    for c in IRENA["Country/area"].unique():
        IRENA.loc[IRENA["Country/area"] == c, "Country/area"] = IRENA_dict[c]

    IRENA = IRENA.set_index(["Country/area", "Technology"])

    IRENA = IRENA.fillna(0).sort_index()

    for (c, t) in IRENA.index.unique():
        sub_df = IRENA.loc[(c, t), :]
        inst_cap = sub_df.loc[sub_df["Indicator"] == "Electricity capacity (MW)", year][0]
        if isinstance(inst_cap, str):
            inst_cap = int(inst_cap.replace(" ", ""))
            IRENA.loc[(IRENA.index.isin([(c, t)])) & (IRENA["Indicator"] == "Electricity capacity (MW)"), year] = inst_cap
        gen_prod = sub_df.loc[sub_df["Indicator"] == "Electricity generation (GWh)", year][0]
        if isinstance(gen_prod, str):
            gen_prod = 1000 * int(gen_prod.replace(" ", ""))
            IRENA.loc[(IRENA.index.isin([(c, t)])) & (IRENA["Indicator"] == "Electricity generation (GWh)"), year] = gen_prod
        if inst_cap == 0:
            FLH = 0
        else:
            FLH = gen_prod / inst_cap
        IRENA = IRENA.append(pd.DataFrame([["FLH (h)", FLH]], index=[(c, t)], columns=["Indicator", year])).sort_index()

    # Filter countries
    IRENA = IRENA.reset_index()
    IRENA = IRENA.set_index(["Country/area"]).sort_index()
    IRENA = IRENA.loc[IRENA.index.isin(filter_countries)]
    # Reshape
    IRENA = IRENA.reset_index()
    IRENA = IRENA.set_index(["Country/area", "Technology"])
    IRENA = IRENA.pivot(columns="Indicator")[year].rename(
        columns={"Electricity capacity (MW)": "inst-cap (MW)", "Electricity generation (GWh)": "prod (MWh)"}
    )
    IRENA = IRENA.astype(float)
    IRENA.to_csv(paths["IRENA_summary"], sep=";", decimal=",", index=True)
    ul.create_json(paths["IRENA_summary"], param, ["author", "comment", "region_name", "year"], paths, ["Countries", "IRENA", "IRENA_dict"])
    logger.info("files saved: " + paths["IRENA_summary"])
