from . import util as ul
from osgeo import gdal, osr, ogr
#from rasterio import MemoryFile
import rasterio.mask
import numpy as np
import math
import rasterio
import scipy.ndimage
import os
import hdf5storage


def define_spatial_scope(scope_shp):
    """
    This function reads the spatial scope shapefile and returns its bounding box.

    :param scope_shp: Spatial scope shapefile.
    :type scope_shp: Geopandas dataframe

    :return box: List of the bounding box coordinates.
    :rtype: list
    """

    lon_min, lat_min, lon_max, lat_max = scope_shp.total_bounds

    return np.array([lat_max, lon_max, lat_min, lon_min])


def crd_merra(Crd_regions, res_weather):
    """
    This function calculates coordinates of the bounding box covering MERRA-2 data.

    :param Crd_regions: Coordinates of the bounding boxes of the regions.
    :type Crd_regions: numpy array
    :param res_weather: Weather data resolution.
    :type res_weather: list

    :return Crd: Coordinates of the bounding box covering MERRA-2 data for each region.
    :rtype: numpy array
    """
    Crd = np.array(
        [
            np.ceil((Crd_regions[0] + res_weather[0] / 2) / res_weather[0]) * res_weather[0] - res_weather[0] / 2,
            np.ceil(Crd_regions[1] / res_weather[1]) * res_weather[1],
            #np.ceil((Crd_regions[:, 1] + res_weather[1] / 2) / res_weather[1]) * res_weather[1] - res_weather[1] / 2,
            np.floor((Crd_regions[2] + res_weather[0] / 2) / res_weather[0]) * res_weather[0] - res_weather[0] / 2,
            np.floor(Crd_regions[3] / res_weather[1]) * res_weather[1],
            #np.floor((Crd_regions[:, 3] + res_weather[1] / 2) / res_weather[1]) * res_weather[1] - res_weather[1] / 2,
        ]
    )
    Crd = Crd.T

    # Todo: New approach because coordinates are not correct (Patrick)
    # resolution_latitude, resolution_longitude = res_weather
    # lat_max, lon_max, lat_min, lon_min = bounding_box
    #
    #
    # Crd = np.array(
    #     [
    #         # |-x-|....|-x-|,     'x': center coordinate of raster, '|': raster
    #         resolution_latitude * np.ceil(lat_max / resolution_latitude) - resolution_latitude / 2,
    #         resolution_longitude * np.ceil(lon_max / resolution_longitude),# - resolution_longitude / 2,
    #         resolution_latitude * np.floor(lat_min / resolution_latitude) - resolution_latitude / 2,
    #         resolution_longitude * np.floor(lon_min / resolution_longitude),# + resolution_longitude / 2
    #     ]
    # )

    return Crd


def ind2crd(Ind_points, Crd_all, resolution):
    """
    This function converts indices of points in high resolution rasters into longitude and latitude coordinates.

    :param Ind_points: Tuple of arrays of indices in the vertical and horizontal axes.
    :type Ind_points: tuple of arrays
    :param Crd_all: Array of coordinates of the bounding box of the spatial scope.
    :type Crd_all: numpy array
    :param resolution: Data resolution in the vertical and horizontal dimensions.
    :type resolution: list
    
    :return Crd_points: Coordinates of the points in the vertical and horizontal dimensions.
    :rtype: list of arrays
    """
    Crd_points = [Ind_points[0] * resolution[0] + Crd_all[2] + resolution[0] / 2, Ind_points[1] * resolution[1] + Crd_all[3] + resolution[1] / 2]
    return Crd_points


def crd2ind(Crd_points, Crd_all, resolution):
    """
    This function converts latitude and longitude of points in high resolution rasters into indices.

    :param Crd_points: Coordinates of the points in the vertical and horizontal dimensions.
    :type Crd_points: tuple of arrays
    :param Crd_all: Array of coordinates of the bounding box of the spatial scope.
    :type Crd_all: numpy array
    :param resolution: Data resolution in the vertical and horizontal dimensions.
    :type resolution: list

    :return Ind_points: Tuple of arrays of indices in the vertical and horizontal axes.
    :rtype: list of arrays
    """

    Ind_points = [np.around((Crd_points[0] - Crd_all[2]) / resolution[0]).astype(int),
                  np.around((Crd_points[1] - Crd_all[3]) / resolution[1]).astype(int)]
    return Ind_points


def subset(A, param):
    """
    This function retrieves a subset of the global MERRA-2 coverage based on weather resolution and the
    bounding box coordinates of the spatial scope.

    :param A: Weather data on a global scale.
    :type A: numpy array
    :param param: Dictionary of parameters containing MERRA-2 coverage and the name of the region.
    :type param: dict

    :return subset: The subset of the weather data contained in the bounding box of *spatial_scope*.
    :rtype: numpy array
    """
    if param["MERRA_coverage"] == "World" and param["region_name"] != "World":
        crd = param["Crd_all"]
        res = param["res_weather"]
        southlim = int(math.floor((crd[2] + res[0] / 10 + 90 + res[0] / 2) / res[0]))
        northlim = int(math.ceil((crd[0] - res[0] / 10 + 90 + res[0] / 2) / res[0]))
        westlim = int(math.floor((crd[3] + res[1] / 10 + 180 + res[1] / 2) / res[1]))
        eastlim = int(math.ceil((crd[1] - res[1] / 10 + 180) / res[1]))
        subset = A[:, southlim:northlim, westlim:eastlim]
    else:
        subset = A
    return subset


def ind_merra(Crd, Crd_all, res):
    """
    This function converts longitude and latitude coordinates into indices within the spatial scope of MERRA-2 data.

    :param Crd: Coordinates to be converted into indices.
    :type Crd: numpy array
    :param Crd_all: Coordinates of the bounding box of the spatial scope.
    :type Crd_all: numpy array
    :param res: Resolution of the data, for which the indices are produced.
    :type res: list
    
    :return Ind: Indices within the spatial scope of MERRA-2 data.
    :rtype: numpy array
    """
    if len(Crd.shape) == 1:
        Crd = Crd[np.newaxis]
    Ind = np.array(
        [
            (Crd[:, 0] - Crd_all[2]) / res[0],
            (Crd[:, 1] - Crd_all[3]) / res[1],
            (Crd[:, 2] - Crd_all[2]) / res[0] + 1,
            (Crd[:, 3] - Crd_all[3]) / res[1] + 1,
        ]
    )
    Ind = np.transpose(Ind).astype(int)
    return Ind


def ind_global(Crd, res_desired):
    """
    This function converts longitude and latitude coordinates into indices on a global data scope, where the origin is at (-90, -180).

    :param Crd: Coordinates to be converted into indices.
    :type Crd: numpy array
    :param res_desired: Desired resolution in the vertical and horizontal dimensions.
    :type res_desired: list
    
    :return Ind: Indices on a global data scope.
    :rtype: numpy array
    """
    if len(Crd.shape) == 1:
        Crd = Crd[np.newaxis]
    Ind = np.array(
        [
            np.round((90 - Crd[:, 0]) / res_desired[0]) + 1,
            np.round((180 + Crd[:, 1]) / res_desired[1]),
            np.round((90 - Crd[:, 2]) / res_desired[0]),
            np.round((180 + Crd[:, 3]) / res_desired[1]) + 1,
        ]
    )
    Ind = np.transpose(Ind.astype(int))
    return Ind


def calc_geotiff(Crd_all, res_desired):
    """
    This function returns a dictionary containing the georeferencing parameters for geotiff creation,
    based on the desired extent and resolution.

    :param Crd_all: Coordinates of the bounding box of the spatial scope.
    :type Crd_all: numpy array
    :param res_desired: Desired data resolution in the vertical and horizontal dimensions.
    :type res_desired: list

    :return GeoRef: Georeference dictionary containing *RasterOrigin*, *RasterOrigin_alt*, *pixelWidth*, and *pixelHeight*.
    :rtype: dict
    """
    GeoRef = {
        "RasterOrigin": [Crd_all[3], Crd_all[0]],
        "RasterOrigin_alt": [Crd_all[3], Crd_all[2]],
        "pixelWidth": res_desired[1],
        "pixelHeight": -res_desired[0],
    }
    return GeoRef


def calc_region(region, Crd_reg, res_desired, GeoRef):
    """
    This function reads the region geometry, and returns a masking raster equal to 1 for pixels within and 0 outside of
    the region.

    :param region: Region geometry
    :type region: Geopandas series
    :param Crd_reg: Coordinates of the region
    :type Crd_reg: list
    :param res_desired: Desired high resolution of the output raster
    :type res_desired: list
    :param GeoRef: Georeference dictionary containing *RasterOrigin*, *RasterOrigin_alt*, *pixelWidth*, and *pixelHeight*.
    :type GeoRef: dict

    :return A_region: Masking raster of the region.
    :rtype: numpy array
    """
    latlim = Crd_reg[2] - Crd_reg[0]
    lonlim = Crd_reg[3] - Crd_reg[1]
    M = int(math.fabs(latlim) / res_desired[0])
    N = int(math.fabs(lonlim) / res_desired[1])
    A_region = np.ones((M, N))
    origin = [Crd_reg[3], Crd_reg[2]]

    if region['geometry'].geom_type == "MultiPolygon":
        features = [feature for feature in region['geometry']]
    else:
        features = [region['geometry']]
    west = origin[0]
    south = origin[1]
    profile = {
        "driver": "GTiff",
        "height": M,
        "width": N,
        "count": 1,
        "dtype": rasterio.float64,
        "crs": "EPSG:4326",
        "transform": rasterio.transform.from_origin(west, south, GeoRef["pixelWidth"], GeoRef["pixelHeight"]),
    }

    with rasterio.MemoryFile() as memfile:
        with memfile.open(**profile) as f:
            f.write(A_region, 1)
            out_image, out_transform = rasterio.mask.mask(f, features, crop=False, nodata=0, all_touched=False, filled=True)
        A_region = out_image[0]

    return A_region


def array2raster(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array):
    """
    This function saves array to geotiff raster format based on EPSG 4326.

    :param newRasterfn: Output path of the raster.
    :type newRasterfn: string
    :param rasterOrigin: Latitude and longitude of the Northwestern corner of the raster.
    :type rasterOrigin: list of two floats
    :param pixelWidth:  Pixel width (might be negative).
    :type pixelWidth: integer
    :param pixelHeight: Pixel height (might be negative).
    :type pixelHeight: integer
    :param array: Array to be converted into a raster.
    :type array: numpy array

    :return: The raster file will be saved in the desired path *newRasterfn*.
    :rtype: None
    """
    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName("GTiff")
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float64, ["COMPRESS=PACKBITS"])
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(np.flipud(array))
    outband.FlushCache()
    outband = None


def aggregate_x_dim(array, res_data, res_desired, aggfun):
    """
    description
    """
    if aggfun == "category":
        array2 = array.reshape(-1, int(res_desired[1] / res_data[1]))
        array2 = np.array([np.bincount(b).argmax() for b in array2]).reshape(array.shape[0], int(array.shape[1]/(res_desired[1] / res_data[1])))
    else:
        array2 = array.reshape(array.shape[0], int(array.shape[1]/(res_desired[1] / res_data[1])), -1)
        if aggfun == "mean":
            array2 = np.mean(array2, 2)
        elif aggfun == "sum":
            array2 = np.sum(array2, 2)
    return array2


def aggregate_y_dim(array, res_data, res_desired, aggfun):
    """
    description
    """
    if aggfun == "category":
        array2 = array.transpose().reshape(-1, int(res_desired[0] / res_data[0]))
        array2 = np.array([np.bincount(b).argmax() for b in array2]).reshape(int(array.shape[0]/(res_desired[0] / res_data[0])), array.shape[1], order="F")
    else:
        array2 = array.transpose().reshape(array.shape[1], int(array.shape[0]/(res_desired[0] / res_data[0])), -1)
        if aggfun == "mean":
            array2 = np.mean(array2, 2).transpose()
        elif aggfun == "sum":
            array2 = np.sum(array2, 2).transpose()
    return array2


def adjust_resolution(array, res_data, res_desired, aggfun=None):
    """
    description
    """
    if ((res_data[1] % res_desired[1] < 1e-10) and (res_data[1] > res_desired[1])): # data is coarse on x dimension (columns)
        array = ul.resizem(array, array.shape[0], int(array.shape[1]*(res_data[1] / res_desired[1])))
        if aggfun == "sum":
            array = array / (res_data[1] % res_desired[1])
    if ((res_desired[1] % res_data[1] < 1e-10) and (res_desired[1] > res_data[1])): # data is too detailed on x dimension
        array = aggregate_x_dim(array, res_data, res_desired, aggfun)
    if ((res_data[0] % res_desired[0] < 1e-10) and (res_data[0] > res_desired[0])): # data is coarse on y dimension (rows)
        array = ul.resizem(array, int(array.shape[0]*(res_data[0] / res_desired[0])), array.shape[1])
        if aggfun == "sum":
            array = array / (res_data[0] % res_desired[0])
    if ((res_desired[0] % res_data[0] < 1e-10) and (res_desired[0] > res_data[0])): # data is too detailed on y dimension
        array = aggregate_y_dim(array, res_data, res_desired, aggfun)
    return array

def recalc_topo_resolution(array, res_data, res_desired):
    
    array1 = np.zeros([int(array.shape[0]*res_data[0]/res_desired[0]),int(array.shape[1]*res_data[1]/res_desired[1])])
    array1[:] = np.NaN
    for i in range(0,array.shape[0],3):
        for j in range(0,array.shape[1],3):
            array1[int(i*5/3),int(j*5/3)] = array[i,j] #first row, first column
            array1[int(i*5/3),int(j*5/3)+2] = array[i,j+1] 
            array1[int(i*5/3),int(j*5/3)+4] = array[i,j+2]
            
            array1[int(i*5/3)+2,int(j*5/3)] = array[i+1,j] #third row, first column
            array1[int(i*5/3)+2,int(j*5/3)+2] = array[i+1,j+1] 
            array1[int(i*5/3)+2,int(j*5/3)+4] = array[i+1,j+2]
            
            array1[int(i*5/3)+4,int(j*5/3)] = array[i+2,j] #fifth row, first column
            array1[int(i*5/3)+4,int(j*5/3)+2] = array[i+2,j+1] 
            array1[int(i*5/3)+4,int(j*5/3)+4] = array[i+2,j+2]

            array1[int(i*5/3),int(j*5/3)+1] = np.mean([array1[int(i*5/3),int(j*5/3)],array1[int(i*5/3),int(j*5/3)+2]],dtype=np.float64)
            array1[int(i*5/3),int(j*5/3)+3] = np.mean([array1[int(i*5/3),int(j*5/3)+2],array1[int(i*5/3),int(j*5/3)+4]],dtype=np.float64)
                            
            array1[int(i*5/3)+1,int(j*5/3)] = np.mean([array1[int(i*5/3),int(j*5/3)],array1[int(i*5/3)+2,int(j*5/3)]],dtype=np.float64)
            array1[int(i*5/3)+1,int(j*5/3)+2] = np.mean([array1[int(i*5/3),int(j*5/3)+2],array1[int(i*5/3)+2,int(j*5/3)+2]],dtype=np.float64)
            array1[int(i*5/3)+1,int(j*5/3)+4] = np.mean([array1[int(i*5/3),int(j*5/3)+4],array1[int(i*5/3)+2,int(j*5/3)+4]],dtype=np.float64)
                          
            array1[int(i*5/3)+2,int(j*5/3)+1] = np.mean([array1[int(i*5/3)+2,int(j*5/3)],array1[int(i*5/3)+2,int(j*5/3)+2]],dtype=np.float64)
            array1[int(i*5/3)+2,int(j*5/3)+3] = np.mean([array1[int(i*5/3)+2,int(j*5/3)+2],array1[int(i*5/3)+2,int(j*5/3)+4]],dtype=np.float64)
                            
            array1[int(i*5/3)+3,int(j*5/3)] = np.mean([array1[int(i*5/3)+2,int(j*5/3)],array1[int(i*5/3)+4,int(j*5/3)]],dtype=np.float64)
            array1[int(i*5/3)+3,int(j*5/3)+2] = np.mean([array1[int(i*5/3)+2,int(j*5/3)+2],array1[int(i*5/3)+4,int(j*5/3)+2]],dtype=np.float64)
            array1[int(i*5/3)+3,int(j*5/3)+4] = np.mean([array1[int(i*5/3)+2,int(j*5/3)+4],array1[int(i*5/3)+4,int(j*5/3)+4]],dtype=np.float64)
                            
            array1[int(i*5/3)+4,int(j*5/3)+1] = np.mean([array1[int(i*5/3)+4,int(j*5/3)],array1[int(i*5/3)+4,int(j*5/3)+2]],dtype=np.float64)
            array1[int(i*5/3)+4,int(j*5/3)+3] = np.mean([array1[int(i*5/3)+4,int(j*5/3)+2],array1[int(i*5/3)+4,int(j*5/3)+4]],dtype=np.float64)
            
            array1[int(i*5/3)+1,int(j*5/3)+1] = np.mean([array1[int(i*5/3)+1,int(j*5/3)],array1[int(i*5/3)+1,int(j*5/3)+2],array1[int(i*5/3),int(j*5/3)+1],array1[int(i*5/3)+2,int(j*5/3)+1]],dtype=np.float64)
            array1[int(i*5/3)+1,int(j*5/3)+3] = np.mean([array1[int(i*5/3)+1,int(j*5/3)+2],array1[int(i*5/3)+1,int(j*5/3)+4],array1[int(i*5/3),int(j*5/3)+3],array1[int(i*5/3)+2,int(j*5/3)+3]],dtype=np.float64)
            
            array1[int(i*5/3)+3,int(j*5/3)+1] = np.mean([array1[int(i*5/3)+3,int(j*5/3)],array1[int(i*5/3)+3,int(j*5/3)+2],array1[int(i*5/3)+2,int(j*5/3)+1],array1[int(i*5/3)+4,int(j*5/3)+1]],dtype=np.float64)
            array1[int(i*5/3)+3,int(j*5/3)+3] = np.mean([array1[int(i*5/3)+3,int(j*5/3)+2],array1[int(i*5/3)+3,int(j*5/3)+4],array1[int(i*5/3)+2,int(j*5/3)+3],array1[int(i*5/3)+4,int(j*5/3)+3]],dtype=np.float64)
                
    return array1
    

def recalc_bath_resolution(array, res_data, res_desired):
    
    array1 = np.zeros([int(array.shape[0]*res_data[0]/res_desired[0]),int(array.shape[1]*res_data[1]/res_desired[1])])
    array1[:] = np.NaN
    for i in range(0,array.shape[0],3):
        for j in range(0,array.shape[1],3):
            #print (str(i)+"_"+str(j))
            array1[int(i*20/3):int(i*20/3)+6,int(j*20/3):int(j*20/3)+6] = array[i,j] #first row, first column
            array1[int(i*20/3):int(i*20/3)+6,int(j*20/3)+7:int(j*20/3)+13] = array[i,j+1] 
            array1[int(i*20/3):int(i*20/3)+6,int(j*20/3)+14:int(j*20/3)+20] = array[i,j+2]
            
            array1[int(i*20/3)+7:int(i*20/3)+13,int(j*20/3):int(j*20/3)+6] = array[i+1,j] #third row, first column
            array1[int(i*20/3)+7:int(i*20/3)+13,int(j*20/3)+7:int(j*20/3)+13] = array[i+1,j+1] 
            array1[int(i*20/3)+7:int(i*20/3)+13,int(j*20/3)+14:int(j*20/3)+20] = array[i+1,j+2]
            
            array1[int(i*20/3)+14:int(i*20/3)+20,int(j*20/3):int(j*20/3)+6] = array[i+2,j] #fifth row, first column
            array1[int(i*20/3)+14:int(i*20/3)+20,int(j*20/3)+7:int(j*20/3)+13] = array[i+2,j+1] 
            array1[int(i*20/3)+14:int(i*20/3)+20,int(j*20/3)+14:int(j*20/3)+20] = array[i+2,j+2]

            array1[int(i*20/3):int(i*20/3)+6,int(j*20/3)+6] = np.mean([array1[int(i*20/3),int(j*20/3)],array1[int(i*20/3),int(j*20/3)+7]],dtype=np.float64)
            array1[int(i*20/3):int(i*20/3)+6,int(j*20/3)+13] = np.mean([array1[int(i*20/3),int(j*20/3)+7],array1[int(i*20/3),int(j*20/3)+14]],dtype=np.float64)
                            
            array1[int(i*20/3)+6,int(j*20/3):int(j*20/3)+6] = np.mean([array1[int(i*20/3),int(j*20/3)],array1[int(i*20/3)+7,int(j*20/3)]],dtype=np.float64)
            array1[int(i*20/3)+6,int(j*20/3)+7:int(j*20/3)+13] = np.mean([array1[int(i*20/3),int(j*20/3)+7],array1[int(i*20/3)+7,int(j*20/3)+7]],dtype=np.float64)
            array1[int(i*20/3)+6,int(j*20/3)+14:int(j*20/3)+20] = np.mean([array1[int(i*20/3),int(j*20/3)+14],array1[int(i*20/3)+7,int(j*20/3)+14]],dtype=np.float64)
                          
            array1[int(i*20/3)+7:int(i*20/3)+13,int(j*20/3)+6] = np.mean([array1[int(i*20/3)+7,int(j*20/3)],array1[int(i*20/3)+7,int(j*20/3)+7]],dtype=np.float64)
            array1[int(i*20/3)+7:int(i*20/3)+13,int(j*20/3)+13] = np.mean([array1[int(i*20/3)+7,int(j*20/3)+7],array1[int(i*20/3)+7,int(j*20/3)+14]],dtype=np.float64)
                            
            array1[int(i*20/3)+13,int(j*20/3):int(j*20/3)+6] = np.mean([array1[int(i*20/3)+7,int(j*20/3)],array1[int(i*20/3)+14,int(j*20/3)]],dtype=np.float64)
            array1[int(i*20/3)+13,int(j*20/3)+7:int(j*20/3)+13] = np.mean([array1[int(i*20/3)+7,int(j*20/3)+7],array1[int(i*20/3)+14,int(j*20/3)+7]],dtype=np.float64)
            array1[int(i*20/3)+13,int(j*20/3)+14:int(j*20/3)+20] = np.mean([array1[int(i*20/3)+7,int(j*20/3)+14],array1[int(i*20/3)+14,int(j*20/3)+14]],dtype=np.float64)
                            
            array1[int(i*20/3)+14:int(i*20/3)+20,int(j*20/3)+6] = np.mean([array1[int(i*20/3)+14,int(j*20/3)],array1[int(i*20/3)+14,int(j*20/3)+7]],dtype=np.float64)
            array1[int(i*20/3)+14:int(i*20/3)+20,int(j*20/3)+13] = np.mean([array1[int(i*20/3)+14,int(j*20/3)+7],array1[int(i*20/3)+14,int(j*20/3)+14]],dtype=np.float64)
            
            array1[int(i*20/3)+6,int(j*20/3)+6] = np.mean([array1[int(i*20/3)+6,int(j*20/3)],array1[int(i*20/3)+6,int(j*20/3)+7],array1[int(i*20/3),int(j*20/3)+6],array1[int(i*20/3)+7,int(j*20/3)+6]],dtype=np.float64)
            array1[int(i*20/3)+6,int(j*20/3)+13] = np.mean([array1[int(i*20/3)+6,int(j*20/3)+7],array1[int(i*20/3)+6,int(j*20/3)+14],array1[int(i*20/3),int(j*20/3)+13],array1[int(i*20/3)+7,int(j*20/3)+13]],dtype=np.float64)
            
            array1[int(i*20/3)+13,int(j*20/3)+6] = np.mean([array1[int(i*20/3)+13,int(j*20/3)],array1[int(i*20/3)+13,int(j*20/3)+7],array1[int(i*20/3)+7,int(j*20/3)+6],array1[int(i*20/3)+14,int(j*20/3)+6]],dtype=np.float64)
            array1[int(i*20/3)+13,int(j*20/3)+13] = np.mean([array1[int(i*20/3)+13,int(j*20/3)+7],array1[int(i*20/3)+13,int(j*20/3)+14],array1[int(i*20/3)+7,int(j*20/3)+13],array1[int(i*20/3)+14,int(j*20/3)+13]],dtype=np.float64)
                
    return array1

    
def recalc_livestock_resolution(array, res_data, res_desired):
    
    #The livestock density (number of animals per sq.km) will be adjusted to higher resolution
    #The new pixel just takes the value from the old pixel with which it has more area
    #It is an approximation
    
    array1 = np.zeros([int(array.shape[0]*res_data[0]/res_desired[0]),int(array.shape[1]*res_data[1]/res_desired[1])])
    array1[:] = np.NaN
    for i in range(0, array.shape[0],3):
        for j in range(0,array.shape[1],3):
            array1[int(i*10/3):int(i*10/3)+3,int(j*10/3):int(j*10/3)+3] = array[i,j]
            array1[int(i*10/3):int(i*10/2)+3,int(j*10/3)+3:int(j*10/3)+7] = array[i,j+1] 
            array1[int(i*10/3):int(i*10/3)+3,int(j*10/3)+7:int(j*10/3)+10] = array[i,j+2]
            
            array1[int(i*10/3)+3:int(i*10/3)+7,int(j*10/3):int(j*10/3)+3] = array[i+1,j]
            array1[int(i*10/3)+3:int(i*10/2)+7,int(j*10/3)+3:int(j*10/3)+7] = array[i+1,j+1] 
            array1[int(i*10/3)+3:int(i*10/3)+7,int(j*10/3)+7:int(j*10/3)+10] = array[i+1,j+2]
            
            array1[int(i*10/3)+7:int(i*10/3)+10,int(j*10/3):int(j*10/3)+3] = array[i+2,j]
            array1[int(i*10/3)+7:int(i*10/2)+10,int(j*10/3)+3:int(j*10/3)+7] = array[i+2,j+1] 
            array1[int(i*10/3)+7:int(i*10/3)+10,int(j*10/3)+7:int(j*10/3)+10] = array[i+2,j+2]
    
    return array1


def create_buffer(param, array, buffer_pixel_amount, GeoRef, Outraster):
    kernel = np.zeros((2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1))
    y, x = np.ogrid[-buffer_pixel_amount:buffer_pixel_amount + 1, -buffer_pixel_amount:buffer_pixel_amount + 1]
    mask = x ** 2 + y ** 2 <= buffer_pixel_amount ** 2
    kernel[mask] = 1
    # kernel = np.tri(2 * buffer_pixel_amount + 1, 2 * buffer_pixel_amount + 1, buffer_pixel_amount).astype(int)
    # kernel = kernel * kernel.T * np.flipud(kernel) * np.fliplr(kernel)
    A_array = scipy.ndimage.maximum_filter(array, footprint=kernel, mode="constant", cval=0)
    A_NotArray = (~A_array).astype(int)

    #saving file
    hdf5storage.writes({"BUFFER": A_NotArray}, Outraster, store_python_metadata=True, matlab_compatible=True)
    if param["savetiff_inputmaps"]:
        array2raster(Outraster, GeoRef["RasterOrigin"], GeoRef["pixelWidth"], GeoRef["pixelHeight"],
                    A_NotArray)


def shape2raster(fileinput, fileoutput, fieldname, dictinput, rastertyp):

    if not os.path.isfile(fileinput):
        raise FileNotFoundError(fileinput)  # Raise error if input file doesn't exist
    else:

        # First we will open our raster image, to understand how we will want to rasterize our vector
        raster_ds = gdal.Open(rastertyp, gdal.GA_ReadOnly)

        # Fetch number of rows and columns
        ncol = raster_ds.RasterXSize
        nrow = raster_ds.RasterYSize

        # Fetch projection and extent
        proj = raster_ds.GetProjectionRef()
        ext = raster_ds.GetGeoTransform()

        raster_ds = None
        shp_path = fileinput
        # Open the dataset from the file
        dataset = ogr.Open(shp_path, 1)
        layer = dataset.GetLayerByIndex(0)

        # Add a new field
        if not ul.field_exists("Raster", shp_path):
            new_field = ogr.FieldDefn("Raster", ogr.OFTInteger)
            layer.CreateField(new_field)

            for feat in layer:
               if dictinput == []:
                    feat.SetField("Raster", 1)
                    layer.SetFeature(feat)
                    feat = None
               else:
                    pt = feat.GetField(fieldname)
                    feat.SetField("Raster", int(dictinput[pt]))
                    layer.SetFeature(feat)
                    feat = None

        # Create a second (modified) layer
        # outdriver = ogr.GetDriverByName("MEMORY")
        # source = outdriver.CreateDataSource("memData")

        # Create the raster dataset
        memory_driver = gdal.GetDriverByName("GTiff")
        out_raster_ds = memory_driver.Create(fileoutput, ncol, nrow, 1, gdal.GDT_Byte)

        # Set the ROI image's projection and extent to our input raster's projection and extent
        out_raster_ds.SetProjection(proj)
        out_raster_ds.SetGeoTransform(ext)

        # Fill our output band with the 0 blank, no class label, value
        b = out_raster_ds.GetRasterBand(1)
        b.Fill(0)

        # Rasterize the shapefile layer to our new dataset
        gdal.RasterizeLayer(
            out_raster_ds,  # output to our new dataset
            [1],  # output to our new dataset's first band
            layer,  # rasterize this layer
            None,
            None,  # don't worry about transformations ul.since we're in same projection
            [0],  # burn value 0
            [
                "ALL_TOUCHED=FALSE",  # rasterize all pixels touched by polygons
                "ATTRIBUTE=Raster",
            ],  # put raster values according to the 'Raster' field values
        )

        # Close dataset
        out_raster_ds = None