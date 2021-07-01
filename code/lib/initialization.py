from config import configuration
from .spatial_functions import *


def initialization():
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
    timecheck("Start")
    # import param and paths
    paths, param = configuration()

    # Read shapefile of scope
    scope_shp = gpd.read_file(paths["spatial_scope"])
    param["spatial_scope"] = define_spatial_scope(scope_shp)

    res_weather = param["res_weather"]
    res_desired = param["res_desired"]
    Crd_all = crd_merra(param["spatial_scope"], res_weather)[0]
    param["Crd_all"] = Crd_all
    ymax, xmax, ymin, xmin = Crd_all
    bounds_box = Polygon([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)])

    timecheck("Read shapefile of countries")
    # Extract land areas
    countries_shp = gpd.read_file(paths["Countries"], bbox=scope_shp)
    countries_shp = countries_shp.to_crs({"init": "epsg:4326"})

    # Crop all polygons and take the part inside the bounding box
    countries_shp["geometry"] = countries_shp["geometry"].buffer(0)
    countries_shp["geometry"] = countries_shp["geometry"].intersection(bounds_box)
    countries_shp = countries_shp[countries_shp.geometry.area > 0]
    param["regions_land"] = countries_shp
    param["nRegions_land"] = len(param["regions_land"])
    Crd_regions_land = np.zeros((param["nRegions_land"], 4))

    for reg in range(0, param["nRegions_land"]):
        # Box coordinates for MERRA2 data
        r = countries_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd_regions_land[reg, :] = crd_merra(box, res_weather)

    timecheck("Read shapefile of EEZ")
    # For the tutorial, a try-except was added. With this it is not necessary to download all raw data.
    try:
        # Extract sea areas
        eez_shp = gpd.read_file(paths["EEZ_global"], bbox=scope_shp)
        eez_shp = eez_shp.to_crs({"init": "epsg:4326"})

        # Crop all polygons and take the part inside the bounding box
        eez_shp["geometry"] = eez_shp["geometry"].buffer(0)
        eez_shp["geometry"] = eez_shp["geometry"].intersection(bounds_box)
        eez_shp = eez_shp[eez_shp.geometry.area > 0]

    except:
        # If the file is missing, create an empty dataframe
        eez_shp = gpd.GeoDataFrame()
        print("#############################################################")
        print("EEZ_global shapefile is missing in the database.")
        print("If you want to run this step, make sure, it is downloaded and at the defined path (check config.py).")
        print("For the tutorial, you can ignore this.")
        print("#############################################################")

    param["regions_sea"] = eez_shp
    param["nRegions_sea"] = len(param["regions_sea"])
    Crd_regions_sea = np.zeros((param["nRegions_sea"], 4))

    for reg in range(0, param["nRegions_sea"]):
        # Box coordinates for MERRA2 data
        r = eez_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd_regions_sea[reg, :] = crd_merra(box, res_weather)


    timecheck("Read shapefile of subregions")
    # Read shapefile of regions
    regions_shp = gpd.read_file(paths["subregions"], bbox=scope_shp)
    regions_shp = regions_shp.to_crs({"init": "epsg:4326"})

    # Crop all polygons and take the part inside the bounding box
    regions_shp["geometry"] = regions_shp["geometry"].intersection(bounds_box)
    regions_shp = regions_shp[regions_shp.geometry.area > 0]
    regions_shp.sort_values(by=["NAME_SHORT"], inplace=True)
    regions_shp.reset_index(inplace=True)
    param["regions_sub"] = regions_shp
    param["nRegions_sub"] = len(param["regions_sub"])
    Crd_regions_sub = np.zeros((param["nRegions_sub"], 4))

    for reg in range(0, param["nRegions_sub"]):
        # Box coordinates for MERRA2 data
        r = regions_shp.bounds.iloc[reg]
        box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
        Crd_regions_sub[reg, :] = crd_merra(box, res_weather)

    # Saving parameters
    param["Crd_subregions"] = Crd_regions_sub
    param["Crd_regions"] = np.concatenate((Crd_regions_land, Crd_regions_sea), axis=0)

    # Indices and matrix dimensions
    Ind_all_low = ind_merra(Crd_all, Crd_all, res_weather)
    Ind_all_high = ind_merra(Crd_all, Crd_all, res_desired)

    param["m_high"] = int((Ind_all_high[:, 0] - Ind_all_high[:, 2] + 1)[0])  # number of rows
    param["n_high"] = int((Ind_all_high[:, 1] - Ind_all_high[:, 3] + 1)[0])  # number of columns
    param["m_low"] = int((Ind_all_low[:, 0] - Ind_all_low[:, 2] + 1)[0])  # number of rows
    param["n_low"] = int((Ind_all_low[:, 1] - Ind_all_low[:, 3] + 1)[0])  # number of columns
    param["GeoRef"] = calc_geotiff(Crd_all, res_desired)
    timecheck("End")

    # Display initial information
    print("\nRegion: " + param["region_name"] + " - Year: " + str(param["year"]))
    print("Folder Path: " + paths["region"] + "\n")

    return paths, param
