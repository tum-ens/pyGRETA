from config import configuration
from lib import spatial_functions as sf
from lib.log import logger
import geopandas as gpd
import numpy as np
import urllib.request
import os
import zipfile
import shapely
# import matplotlib.pyplot as plt


def initialization(config_file):
    """
    This function reads the user-defined parameters and paths from :mod:`config.py`, then adds additional parameters related
    to the shapefiles. 
    First, it saves the spatial scope of the problem.
    Then, it distinguishes between countries, exclusive economic zones and subregions. For each one of them, 
    it saves the geodataframes, the number of features, and the coordinates of the bounding boxes of each feature.
    Finally, it saves the number of rows and columns in the low and righ resolution, and a georeference dictionary
    used for saving tif files.

   
   :return: The updated dictionaries param and paths.
    :rtype: tuple(dict, dict)
    """
    logger.info("Start")

    paths, param = configuration(config_file)   # Load configurations depended on the recent config_file
    param = read_shapefiles(paths, param)       # Import land and eez regions based on the given shapefile

    # Indices and matrix dimensions
    Crd_all = param["Crd_all"]
    Ind_all_low = sf.ind_merra(Crd_all, Crd_all, param["res_weather"])
    Ind_all_high = sf.ind_merra(Crd_all, Crd_all, param["res_desired"])

    param["m_high"] = int((Ind_all_high[:, 0] - Ind_all_high[:, 2] + 1)[0])  # number of rows
    param["n_high"] = int((Ind_all_high[:, 1] - Ind_all_high[:, 3] + 1)[0])  # number of columns
    param["m_low"] = int((Ind_all_low[:, 0] - Ind_all_low[:, 2] + 1)[0])  # number of rows
    param["n_low"] = int((Ind_all_low[:, 1] - Ind_all_low[:, 3] + 1)[0])  # number of columns
    param["GeoRef"] = sf.calc_geotiff(Crd_all, param["res_desired"])

    # Display initial information
    logger.info("Region: " + param["region_name"] + " - Year: " + str(param["year"]))
    logger.info("Folder Path: " + paths["region"])

    return paths, param

def read_shapefiles(paths, param):
    logger.info("Start")
    scope_shp = gpd.read_file(paths["subregions"])      # Import shapefile of regions
    scope_shp = scope_shp.to_crs(epsg=4326)             # In case it is not already in this format

    param["spatial_scope"] = sf.define_spatial_scope(scope_shp)
    param["Crd_all"] = sf.crd_merra(param["spatial_scope"], param["res_weather"])   # rectangle coordinates

    param = read_regions(param, scope_shp)
    param = read_EEZ(paths, param, scope_shp)

    logger.debug("End")
    return param

def read_regions(param, regions_shp):
    logger.info("Start")
    # regions_shp.sort_values(by=["NAME_SHORT"], inplace=True)
    regions_shp.sort_values(by=[param["gid"]], inplace=True)
    regions_shp.reset_index(inplace=True)
    param["regions_land"] = regions_shp
    param["nRegions_land"] = len(regions_shp)

    Crd_regions_land = np.zeros((param["nRegions_land"], 4))
    for reg in range(0, param["nRegions_land"]):
        r = regions_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])    # Box coordinates for MERRA2 data
        Crd_regions_land[reg, :] = sf.crd_merra(box, param["res_weather"])
    param["Crd_regions_land"] = Crd_regions_land
    logger.debug('End')

    return param

def read_EEZ(paths, param, scope_shp):
    logger.info("Start")

    eez_shp = gpd.read_file(paths["EEZ_global"])    # Import all exclusive economic zones worldwide
    eez_shp = eez_shp.to_crs(epsg=4326)             # In case it is not already in this format

    countries = scope_shp['ID_0'].drop_duplicates()        # Obtain the countries out of the regions shapefile
    eez_shp = eez_shp[eez_shp['ISO_Ter1'].isin(countries)].reset_index()    # Remain only eez areas of countries involved

    param["offshore_scope"] = sf.define_spatial_scope(eez_shp)
    param["Crd_offshore"] = sf.crd_merra(param["offshore_scope"], param["res_weather"])  # rectangle coordinates

    Crd_offshore = param["Crd_offshore"]
    Ind_offshore_low = sf.ind_merra(Crd_offshore, Crd_offshore, param["res_weather"])
    Ind_offshore_high = sf.ind_merra(Crd_offshore, Crd_offshore, param["res_desired"])

    param["m_high_offshore"] = int((Ind_offshore_high[:, 0] - Ind_offshore_high[:, 2] + 1)[0])  # number of rows
    param["n_high_offshore"] = int((Ind_offshore_high[:, 1] - Ind_offshore_high[:, 3] + 1)[0])  # number of columns
    param["m_low_offshore"] = int((Ind_offshore_low[:, 0] - Ind_offshore_low[:, 2] + 1)[0])  # number of rows
    param["n_low_offshore"] = int((Ind_offshore_low[:, 1] - Ind_offshore_low[:, 3] + 1)[0])  # number of columns
    param["GeoRef_offshore"] = sf.calc_geotiff(Crd_offshore, param["res_desired"])

    param["regions_sea"] = eez_shp
    param["nRegions_sea"] = len(eez_shp)

    Crd_regions_sea = np.zeros((param["nRegions_sea"], 4))
    for reg in range(0, param["nRegions_sea"]):
        r = eez_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])    # Box coordinates for MERRA2 data
        Crd_regions_sea[reg, :] = sf.crd_merra(box, param["res_weather"])
    param['Crd_regions_sea'] = Crd_regions_sea

    return param
