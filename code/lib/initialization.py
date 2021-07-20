from config import configuration
from .spatial_functions import *
import numpy as np
from .log import logger
import logging
#import matplotlib.pyplot as plt


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

    param = read_shapefiles(paths, param)    # Load

    # Saving parameters
    param["Crd_regions"] = np.concatenate((param['Crd_regions_land'], param['Crd_regions_sea']), axis=0)

    # Indices and matrix dimensions
    Crd_all = param["Crd_all"]
    Ind_all_low = ind_merra(Crd_all, Crd_all, param["res_weather"])
    Ind_all_high = ind_merra(Crd_all, Crd_all, param["res_desired"])

    param["m_high"] = int((Ind_all_high[:, 0] - Ind_all_high[:, 2] + 1)[0])  # number of rows
    param["n_high"] = int((Ind_all_high[:, 1] - Ind_all_high[:, 3] + 1)[0])  # number of columns
    param["m_low"] = int((Ind_all_low[:, 0] - Ind_all_low[:, 2] + 1)[0])  # number of rows
    param["n_low"] = int((Ind_all_low[:, 1] - Ind_all_low[:, 3] + 1)[0])  # number of columns
    param["GeoRef"] = calc_geotiff(Crd_all, param["res_desired"])
    # timecheck("End")

    # Display initial information
    logger.info("Region: " + param["region_name"] + " - Year: " + str(param["year"]))
    logger.info("Folder Path: " + paths["region"])

    return paths, param

# def read_shapfile_countries(paths, param, scope_shp, bounds_box):
#     #logger.setLevel(logging.DEBUG)
#     logger.info("Read shapefile of countries")
#     # Extract land areas
#     countries_shp = gpd.read_file(paths["Countries"], bbox=scope_shp)
#     countries_shp.plot()
#     logger.debug('File read')
#     countries_shp_epsg = countries_shp.to_crs({"init": "epsg:4326"})
#     countries_shp_epsg.plot()
#     logger.debug('shp')
#     # Crop all polygons and take the part inside the bounding box
#     countries_shp_buffer = countries_shp_epsg["geometry"].buffer(0)
#     countries_shp_buffer.plot()
#     logging.debug('buffer')
#     countries_shp["geometry"] = countries_shp_buffer.intersection(bounds_box)
#     countries_shp.plot()
#     plt.show()
#     logger.debug('intersection')
#     countries_shp = countries_shp[countries_shp.geometry.area > 0]
#     param["regions_land"] = countries_shp
#     param["nRegions_land"] = len(param["regions_land"])
#
#     logger.debug('cropped')
#     Crd_regions_land = np.zeros((param["nRegions_land"], 4))
#     for reg in range(0, param["nRegions_land"]):
#         # Box coordinates for MERRA2 data
#         r = countries_shp.bounds.iloc[reg]
#         box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
#         Crd_regions_land[reg, :] = crd_merra(box, param["res_weather"])
#     param['Crd_regions_land'] = Crd_regions_land
#     logger.debug('done')
#
#     return param

def read_shapefiles(paths, param):

    scope_shp = gpd.read_file(paths["subregions"])      # Read shapefile of region
    scope_shp = scope_shp.to_crs({"init": "epsg:4326"})     # If it is not already in this format
    param["spatial_scope"] = define_spatial_scope(scope_shp)

    param["Crd_all"] = crd_merra(param["spatial_scope"], param["res_weather"])[0]    # rectangle coordinates

    param = read_shapefile_regions(param, scope_shp)
    param = read_shapefile_EEZ(paths, param, scope_shp)

    return param


def read_shapefile_regions(param, regions_shp):
    logger.info("Read shapefile of regions")

    regions_shp = regions_shp[regions_shp.geometry.area > 0] # ToDo: Keep it?
    # regions_shp.sort_values(by=["NAME_SHORT"], inplace=True)
    regions_shp.sort_values(by=["GID_0"], inplace=True)     # ToDo: Replace 'GID_0' by variable
    regions_shp.reset_index(inplace=True)
    param["regions_land"] = regions_shp
    param["nRegions_land"] = len(regions_shp)

    Crd_regions_land = np.zeros((param["nRegions_land"], 4))
    for reg in range(0, param["nRegions_land"]):
        # Box coordinates for MERRA2 data
        r = regions_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd_regions_land[reg, :] = crd_merra(box, param["res_weather"])
    param["Crd_regions_land"] = Crd_regions_land
    logger.debug('Finished')

    return param


def read_shapefile_EEZ(paths, param, scope_shp):
    logger.info("Read shapefile of EEZ")

    ymax, xmax, ymin, xmin = param['Crd_all']
    bounds_box = Polygon([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)])
    eez_shp = gpd.read_file(paths["EEZ_global"], bbox=scope_shp)    # Extract sea areas

    eez_shp = eez_shp.to_crs({"init": "epsg:4326"})
    eez_shp['geometry'] = eez_shp['geometry'].intersection(bounds_box)  # Crop all polygons and take the part inside the bounding box
    eez_shp = eez_shp[eez_shp.geometry.area > 0]
    param["regions_sea"] = eez_shp
    param["nRegions_sea"] = len(eez_shp)

    Crd_regions_sea = np.zeros((param["nRegions_sea"], 4))
    for reg in range(0, param["nRegions_sea"]):
        # Box coordinates for MERRA2 data
        r = eez_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd_regions_sea[reg, :] = crd_merra(box, param["res_weather"])
    param['Crd_regions_sea'] = Crd_regions_sea

    return param


def test(paths, param, scope_shp, bounds_box):
    regions_shp = gpd.read_file(paths["subregions"], bbox=scope_shp)
    regions_shp.plot()
    regions_shp_2 = gpd.read_file(paths["subregions"])
    regions_shp_2.plot()
    countries_shp = gpd.read_file(paths["Countries"], bbox=scope_shp)
    countries_shp.plot()
    countries_shp = countries_shp.intersection(bounds_box)  # Crop all polygons and take the part inside the bounding box
    countries_shp.plot()
    countries_shp = countries_shp[countries_shp.geometry.area > 0]
    countries_shp.plot()
    plt.show()