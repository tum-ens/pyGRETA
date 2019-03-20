import os
# os.environ['MKL_NUM_THREADS'] = '8'
from SolarFunc import *
from WindFunc import *
from config import *
import numpy as np
import geopandas as gpd
from rasterio import windows
from shapely.geometry import mapping, MultiPoint
import fiona
import h5py
import hdf5storage
from multiprocessing import Pool
from itertools import product
#import pickle
#import sys

# Initialization - Enter Location Info
# read shapefile of regions
regions_shp = gpd.read_file(paths["SHP"])

if technology == 'PV' or technology == 'CSP' or (technology == 'Wind' and 'Offshore' not in windtechnology):
    # Remove Sea regions from regions_shp
    regions_shp = regions_shp.drop(regions_shp[regions_shp["Population"] == 0].index)
elif technology == 'Wind' and 'Onshore' not in windtechnology:
    # Remove Land regions from regions_shp
    regions_shp = regions_shp.drop(regions_shp[regions_shp["Population"] != 0].index)

nRegions = len(regions_shp) + 1

Crd = np.zeros((nRegions, 4))
for reg in range(1, nRegions):
    # Box coordinates for MERRA2 data
    r = regions_shp.bounds.iloc[reg - 1]
    box = np.array([r["maxy"], r["maxx"], r["miny"], r["minx"]])[np.newaxis]
    Crd[reg, :] = crd_merra_low(box, res)
Crd[0, :] = [max(Crd[1:, 0]), max(Crd[1:, 1]), min(Crd[1:, 2]), min(Crd[1:, 3])]

# Indices and matrix dimensions
Ind = np.zeros((3 + len(quantiles), nRegions, 4))
Ind[0] = ind_merra_low(Crd, res)  # Range indices for MERRA2 data (centroids)
Ind[1] = ind_merra_high(Crd, res)  # Range indices for high resolution matrices, superposed to MERRA2 data
Ind[2] = ind_global(Crd, res[1, :])
Ind = Ind.astype(int)

m = Ind[:, :, 0] - Ind[:, :, 2] + 1  # #rows
m[2, :] = Ind[2, :, 2] - Ind[2, :, 0] + 1  # starts counting from North
n = Ind[:, :, 1] - Ind[:, :, 3] + 1  # #Cols
m = m.astype(int)
n = n.astype(int)

R1 = calc_geotiff(Crd, res)

generate_landsea(paths, regions_shp, Crd, Ind, m, n, res, R1)  # Land and Sea
generate_landuse(paths, Ind, R1)  # Landuse
generate_bathymetry(paths, Ind, R1)  # Bathymetry
generate_topography(paths, Ind, R1)  # Topography
generate_slope(paths, Ind, R1)  # Slope
generate_population(paths, Ind, R1)  # Population
generate_protected_areas(paths, protected_areas, R1)  # Protected areas
generate_buffered_population(paths, buffer_pixel_amount, R1)  # Buffered Population

# Access to data files
merraData = {}
if technology == 'PV' or technology == 'CSP':
    # Downward shortwave radiation on the ground - stored variable SWGDN
    with h5py.File(paths["GHI"]) as f:
        merraData["file_pv"] = {"SWGDN": np.array(f["SWGDN"][:]).T}
    # Downward shortwave radiation at the top of the atmosphere SWTDN
    with h5py.File(paths["TOA"]) as f:
        merraData["file_toa"] = {"SWTDN": np.array(f["SWTDN"][:]).T}
    # Temperature 2m above the ground - stored variable T2M
    with h5py.File(paths["T2M"]) as f:
        merraData["file_t"] = {"T2M": np.array(f["T2M"][:]).T}
elif technology == 'Wind':
    with h5py.File(paths["U50M"]) as f:
        merraData["file_u50m"] = {"U50M": np.array(f["U50M"][:]).T}
    with h5py.File(paths["V50M"]) as f:
        merraData["file_v50m"] = {"V50M": np.array(f["V50M"][:]).T}
    W50M = abs(merraData["file_u50m"]["U50M"] + (1j * merraData["file_v50m"]["V50M"]))

# Read land use map for wind adjustment

if technology == 'Wind':
    with rasterio.open(paths["LU"]) as src:
        A_lu = np.flipud(src.read(1)).astype(int)
    A_hellmann = changem(A_lu, landuse["hellmann"], landuse["type"]).astype(float) / 100
    A_cf = {}
    if correction:
        A_gradient_height = changem(A_lu.astype(float), landuse["height"], landuse["type"])
        del A_lu
        Sigma = sumnorm_MERRA2((50 / A_gradient_height) ** A_hellmann, m[0, 0], n[0, 0], res)
        if 'Onshore' in windtechnology:
            A_cf["Onshore"] = ((turbine["Onshore"]["hub_height"] / 50) * turbine["Onshore"]["hub_height"] /
                               A_gradient_height) ** A_hellmann / resizem(Sigma, m[1, 0], n[1, 0])
        if 'Offshore' in windtechnology:
            A_cf["Offshore"] = ((turbine["Offshore"]["hub_height"] / 50) * turbine["Offshore"]["hub_height"] /
                                A_gradient_height) ** A_hellmann / resizem(Sigma, m[1, 0], n[1, 0])
        del A_hellmann, A_gradient_height, Sigma

    else:
        if 'Onshore' in windtechnology:
            A_cf = (turbine["Onshore"]["hub_height"] / 50) ** A_hellmann
        if 'Offshore' in windtechnology:
            A_cf = (turbine["Offshore"]["hub_height"] / 50) ** A_hellmann
        del A_lu, A_hellmann

# Per Region

for reg in range(1, nRegions):
    description["region"] = reg - 1
    suffix = ''
    if technology == 'Wind' and regions_shp["Population"][reg - 1] == 0:
        tech = 'Offshore'
        suffix = '_offshore'
    elif technology == 'Wind':
        tech = 'Onshore'
    region_name = regions_shp["NAME_SHORT"][reg - 1] + suffix

    # Calculate A matrices
    rasterData = {"A_region": calc_region(regions_shp.iloc[reg - 1], Crd[reg, :], res, R1)}

    with rasterio.open(paths["LU"]) as src:
        w = src.read(1, window=rasterio.windows.Window.from_slices(slice((m[1, 0] - Ind[1, reg, 0]),
                                                                         (m[1, 0] - Ind[1, reg, 2] + 1)),
                                                                   slice(Ind[1, reg, 3] - 1,
                                                                         Ind[1, reg, 1])))
    rasterData["A_lu"] = np.flipud(w)

    # Calculate special A matrices
    # Topographic map with 15 arcsec resolution
    with rasterio.open(paths["TOPO"]) as src:
        w = src.read(1, window=rasterio.windows.Window.from_slices(slice((m[1, 0] - Ind[1, reg, 0]),
                                                                         (m[1, 0] - Ind[1, reg, 2] + 1)),
                                                                   slice(Ind[1, reg, 3] - 1,
                                                                         Ind[1, reg, 1])))
        rasterData["A_topo"] = np.flipud(w)

    if technology == 'PV' or technology == 'CSP':
        # Temperature coefficients for heating losses
        rasterData["A_Ross"] = changem(rasterData["A_lu"], landuse["Ross_coeff"], landuse["type"]).astype(float) / 10000
        # Reflectivity coefficients
        rasterData["A_albedo"] = changem(rasterData["A_lu"], landuse["albedo"], landuse["type"]).astype(float) / 100

        # Calculate CLR_Mean and CLR_MAX
        _, CLR_MAX = calc_clearness(merraData, reg, Ind)

    elif technology == 'Wind':

        rasterData["A_cf"] = A_cf[tech][Ind[1, reg, 2] - 1:Ind[1, reg, 0], Ind[1, reg, 3] - 1:Ind[1, reg, 1]] * \
                             np.exp(a * rasterData["A_topo"] / 5000 - b)
	
    list_hours = []
    chunk = 20 // nproc
    for i in range(0, nproc-1):
        list_hours.append(list(range(chunk*i, chunk*(i+1))))
    list_hours.append(list(range(chunk*(nproc-1), 20)))
	
    if technology == 'PV' or technology == 'CSP':
        args = [reg, nRegions, region_name, Ind, Crd, res, merraData, rasterData, m, n, technology, CLR_MAX, pv]
        myfun = calc_FLH_solar
        #test = pickle.dumps(args)
        #print(region_name, sys.getsizeof(test))
    elif technology == 'Wind':
        args = [reg, nRegions, region_name, Ind, Crd, res, merraData, rasterData, m, n, W50M, turbine, tech]
        myfun = calc_FLH_wind
    results = Pool(processes=nproc).starmap(myfun, product(list_hours, [args]))
	
    # Collecting results
    TS = np.zeros((8760, 1))
    FLH = np.zeros((m[1, reg], n[1, reg]))
    for p in range(len(results)):
        FLH = FLH + results[p][0]
        TS = TS + results[p][1]

    FLH[FLH == 0] = np.nan
    TS[np.isnan(TS)] = 0

    paths["FLH"] = paths["OUT"] + technology + '_FLH_' + region_name + '_' + year + '.mat'
    paths["TS_mean"] = paths["OUT"] + technology + '_Mean_Timeseries_' + region_name + '_' + year + '.mat'
    description["paths"] = paths
    if not os.path.isdir(paths["OUT"]):
        os.mkdir(paths["OUT"])

    hdf5storage.writes({'FLH':FLH, 'description':description}, paths["FLH"], store_python_metadata=True, matlab_compatible=True)
    print("files saved:" + paths["FLH"])

    hdf5storage.writes({'TS_mean':TS}, paths["TS_mean"], store_python_metadata=True, matlab_compatible=True)
    print("files saved:" + paths["TS_mean"])

    del FLH, TS

# New!  Aggregate regions in one single file
S = np.array(list_files(paths["OUT"], technology + '_FLH_' + '*.mat'))
if S.size != 0:
    FLH_all = np.zeros((m[1, 0], n[1, 0]))
    for pd in S:
        with h5py.File(paths["OUT"] + pd, 'r') as f:
            FLH = np.array(f["FLH"])
            reg = f["description/region"][()] + 1
            FLH_all[Ind[1, reg, 2] - 1:Ind[1, reg, 0], Ind[1, reg, 3] - 1:Ind[1, reg, 1]] = \
                np.fmax(FLH_all[Ind[1, reg, 2] - 1:Ind[1, reg, 0], Ind[1, reg, 3] - 1:Ind[1, reg, 1]], FLH)

    paths["FLH"] = paths["OUT"] + region + '_' + technology + '_FLH_' + year + '.mat'
    description["region"] = 'all'
    description["region"] = paths

    with h5py.File(paths["FLH"], 'w') as f:
        f.create_dataset('FLH_all', data=FLH_all)
        recursive_dict_save(f, 'description/', description)
    print("files saved:" + paths["FLH"])
    del FLH, FLH_all

# Mask
description["region"] = region
with h5py.File(paths["FLH"], 'r') as f:
    FLH_all = np.array(f["FLH_all"])

A_mask, A_FLH_mask = masking(FLH_all, mask, landuse, paths, technology, windtechnology, buffer_pixel_amount)

# SAVE HDF5 Files and GEOTIFF

with h5py.File(paths["mask"], 'w') as f:
    f.create_dataset('A_mask', data=A_mask)
    recursive_dict_save(f, 'description/', description)
print("files saved:" + paths["mask"])

with h5py.File(paths["FLH_mask"], 'w') as f:
    f.create_dataset('A_FLH_mask', data=A_FLH_mask)
    recursive_dict_save(f, 'description/', description)
print("files saved:" + paths["FLH_mask"])

if savetiff:
    array2raster(changeExt2tif(paths["mask"]),
                 R1["RasterOrigin"],
                 R1["pixelWidth"],
                 R1["pixelheight"],
                 A_mask)
    print("files saved:" + changeExt2tif(paths["mask"]))

    array2raster(changeExt2tif(paths["FLH_mask"]),
                 R1["RasterOrigin"],
                 R1["pixelWidth"],
                 R1["pixelheight"],
                 A_FLH_mask)
    print("files saved:" + changeExt2tif(paths["FLH_mask"]))

del FLH_all, A_mask, A_FLH_mask

# Weight
description["region"] = region
with h5py.File(paths["FLH_mask"], 'r') as f:
    A_FLH_mask = np.array(f["A_FLH_mask"])

with h5py.File(paths["mask"], 'r') as f:
    A_mask = np.array(f["A_mask"])

A_area, A_weight, A_FLH_weight = weighting(A_FLH_mask, weight, landuse, paths, technology, windtechnology, Crd, n, res, protected_areas)

# Compute Sums

A_area_mask = A_area * A_mask
area_sum_km_2 = np.nansum(A_area_mask) / 1e6
power_potential_sum_TWp = np.nansum(A_weight) / 1e6
energy_potential_sum_TWh = np.nansum(A_FLH_weight) / 1e6

# SAVE HDF5 Files and GEOTIFF

with h5py.File(paths["area"], 'w') as f:
    f.create_dataset('A_area', data=A_area)
    f.create_dataset('area_sum_km_2', data=area_sum_km_2)
    recursive_dict_save(f, 'description/', description)
print("files saved:" + paths["area"])

with h5py.File(paths["weight"], 'w') as f:
    f.create_dataset('A_weight', data=A_weight)
    f.create_dataset('power_potential_sum_TWp', data=power_potential_sum_TWp)
    recursive_dict_save(f, 'description/', description)
print("files saved:" + paths["weight"])

with h5py.File(paths["FLH_weight"], 'w') as f:
    f.create_dataset('A_FLH_weight', data=A_FLH_weight)
    f.create_dataset('energy_potential_sum_TWh', data=energy_potential_sum_TWh)
    recursive_dict_save(f, 'description/', description)
print("files saved:" + paths["FLH_weight"])

if savetiff:
    array2raster(changeExt2tif(paths["FLH_weight"]),
                 R1["RasterOrigin"],
                 R1["pixelWidth"],
                 R1["pixelheight"],
                 A_FLH_weight)
    print("files saved:" + changeExt2tif(paths["FLH_weight"]))

    array2raster(changeExt2tif(paths["area"]),
                 R1["RasterOrigin"],
                 R1["pixelWidth"],
                 R1["pixelheight"],
                 A_area)
    print("files saved:" + changeExt2tif(paths["area"]))

    array2raster(changeExt2tif(paths["weight"]),
                 R1["RasterOrigin"],
                 R1["pixelWidth"],
                 R1["pixelheight"],
                 A_weight)
    print("files saved:" + changeExt2tif(paths["weight"]))

# ## Cost (To be done later)

# Time series per region

with h5py.File(paths["FLH_mask"], 'r') as f:
    A_FLH_mask = np.array(f["A_FLH_mask"])

for reg in range(1, nRegions):
    description["region"] = reg - 1
    suffix = ''
    if technology == 'Wind' and regions_shp["Population"][reg - 1] == 0:
        tech = 'Offshore'
        suffix = '_offshore'
    elif technology == 'Wind':
        tech = 'Onshore'
    if technology == 'PV' and regions_shp["Population"][reg - 1] == 0:
        continue
    region_name = regions_shp["NAME_SHORT"][reg - 1] + suffix

    # Calculate A matrices
    rasterData["A_region"] = calc_region(regions_shp.iloc[reg - 1], Crd[reg, :], res, R1)

    # Find position of quantiles
    Locations = np.zeros((len(quantiles), 4))
    FLH_reg = rasterData["A_region"] * A_FLH_mask[Ind[1, reg, 2] - 1:Ind[1, reg, 0], Ind[1, reg, 3] - 1:Ind[1, reg, 1]]
    FLH_reg[FLH_reg == 0] = np.nan
    X = FLH_reg.flatten(order='F')
    I_old = np.argsort(X)

    # ESCAPE FOR LOOP IF INTERSECTION WITH RASTER ONLY YIELDS NAN
    for q in range(0, len(quantiles)):
        if quantiles[q] == 100:
            I = I_old[(len(X) - 1) - sum(np.isnan(X).astype(int))]
        elif quantiles[q] == 0:
            I = I_old[0]
        else:
            I = I_old[int(np.round(quantiles[q] / 100 * (len(X) - 1 - sum(np.isnan(X).astype(int)))))]

        I, J = ind2sub(FLH_reg.shape, I)  # Convert the indices to row-column indices
        I = I + 1
        J = J + 1
        Ind[3 + q, reg, :] = [(I + Ind[1, reg, 2] - 1), (J + Ind[1, reg, 3] - 1),
                              (I + Ind[1, reg, 2] - 1), (J + Ind[1, reg, 3] - 1)]
        Locations[q, :] = crd_exact_high(np.squeeze(Ind[3 + q, reg, :]).T, Crd, res)

    Locations[:, 2:4] = repmat(Crd[0, 2:4], Locations.shape[0], 1)
    Locations_ind_low = ind_merra_low(crd_merra_low(Locations, res), res)
    Locations_ind_high = np.squeeze(Ind[3:, reg, :]) - \
                         repmat(np.squeeze(Ind[1, reg, 2:4]).T, Locations.shape[0], 2) + 1

    # Calculate A matrices
    if Ind[1, reg, 0] == m[1, 0]:
        # Landuse classes 0-16, to be reclassified
        with rasterio.open(paths["LU"]) as src:
            w = src.read(1, window=rasterio.windows.Window.from_slices(slice((m[1, 0] - Ind[1, reg, 0]),
                                                                             (m[1, 0] - Ind[1, reg, 2] + 1)),
                                                                       slice(Ind[1, reg, 3] - 1,
                                                                             Ind[1, reg, 1])))
            rasterData["A_lu"] = np.flipud(w)
    else:
        with rasterio.open(paths["LU"]) as src:
            w = src.read(1, window=rasterio.windows.Window.from_slices(slice((m[1, 0] - Ind[1, reg, 0] - 1),
                                                                             (m[1, 0] - Ind[1, reg, 2])),
                                                                       slice(Ind[1, reg, 3] - 1,
                                                                             Ind[1, reg, 1])))
            rasterData["A_lu"] = np.flipud(w)

    if technology == 'PV' or technology == 'CSP':
        # Calculate special A matrices
        with rasterio.open(paths["TOPO"]) as src:
            w = src.read(1, window=rasterio.windows.Window.from_slices(slice((m[1, 0] - Ind[1, reg, 0] - 1),
                                                                             (m[1, 0] - Ind[1, reg, 2])),
                                                                       slice(Ind[1, reg, 3] - 1,
                                                                             Ind[1, reg, 1])))
            rasterData["A_topo"] = np.flipud(w)

        # Temperature coefficients for heating losses
        rasterData["A_Ross"] = changem(rasterData["A_lu"], landuse["Ross_coeff"], landuse["type"]).astype(
            float) / 10000
        # Reflectivity coefficients
        rasterData["A_albedo"] = changem(rasterData["A_lu"], landuse["albedo"], landuse["type"]).astype(float) / 100

        # Calculate CLR_Mean and CLR_MAX
        _, CLR_MAX = calc_clearness(merraData, reg, Ind)

    elif technology == 'Wind':

        rasterData["A_cf"] = A_cf[tech][np.ix_(Ind[3:, reg, 0] - 1, Ind[3:, reg, 1] - 1)]
        rasterData["A_cf"] = np.diagonal(rasterData["A_cf"])

    TS = np.zeros((len(quantiles), 8760))

    for hour in range(0, 8760):
        # Show progress of the simulation
        print(str(reg) + '/' + str(nRegions - 1) + ' ' + region_name + ' ' + str(hour + 1))

        CF = 0
        if technology == 'PV':
            CF, _ = calc_CF_solar(hour, reg, Ind, Crd, res, CLR_MAX, merraData, rasterData, pv, Locations)
        elif technology == 'CSP':
            _, CF = calc_CF_solar(hour, reg, Ind, Crd, res, CLR_MAX, merraData, rasterData, pv, Locations)
        elif technology == 'Wind':
            # Load MERRA data, increase its resolution, and fit it to the extent
            Locations_ind_low = ind_merra_low(crd_merra_low(Locations, res), res)
            w50m_h = W50M[:, :, hour]
            w50m_h = w50m_h[np.ix_(Locations_ind_low[:, 0] - 1, Locations_ind_low[:, 1] - 1)]
            w50m_h = np.diagonal(w50m_h)

            # Calculate hourly capacity factor
            CF = calc_CF_wind(w50m_h, rasterData, turbine, tech)

        TS[:, hour] = CF

    TS[np.isnan(TS)] = 0

    paths["TS_quantiles"] = paths["OUT"] + technology + '_Timeseries_' + region_name + '_' + year + '.mat'
    paths["Locations"] = paths["OUT"] + technology + '_Locations_' + region_name + '_' + year + '.mat'
    description["paths"] = paths

    with h5py.File(paths["TS_quantiles"], 'w') as f:
        f.create_dataset('TS', data=TS)
        recursive_dict_save(f, 'description/', description)
    print("files saved:" + paths["TS_quantiles"])

    with h5py.File(paths["Locations"], 'w') as f:
        f.create_dataset('Locations', data=Locations)
        recursive_dict_save(f, 'description/', description)
    print("files saved:" + paths["Locations"])
    del TS

# ## Timeseries for all regions

S = np.array(list_files(paths["OUT"], technology + '_Timeseries_' + '*.mat'))
if S.size != 0:
    timeseries = {}
    for f in range(0, len(S)):
        with h5py.File(paths["OUT"] + S[f], 'r') as src:
            TS = np.array(src["TS"])
            reg = src["description/region"][()] + 1
        suffix = ''

        if technology == 'PV':
            suffix = '.solar'
        elif technology == 'CSP':
            suffix = '.beam'
        elif technology == 'Wind':
            if regions_shp["Population"][reg - 1] == 0:
                suffix = '.WindOff'
            else:
                suffix = '.WindOn'

        timeseries[regions_shp["NAME_SHORT"][reg - 1] + suffix] = {}

        for q in range(0, len(quantiles)):
            newfield = 'q' + str(quantiles[q])
            timeseries[regions_shp["NAME_SHORT"][reg - 1] + suffix][newfield] = TS[q, :].T
TS = timeseries
del timeseries

description["region"] = region
paths["TS"] = paths["OUT"] + region + '_' + technology + '_Timeseries_' + year + '.mat'
description["paths"] = paths

with h5py.File(paths["TS"], 'w') as f:
    recursive_dict_save(f, 'Region/', TS)
    recursive_dict_save(f, 'description/', description)
print("files saved:" + paths["TS"])

# Locations for all regions
S = np.array(list_files(paths["OUT"], technology + '_Locations_' + '*.mat'))
if S.size != 0:
    Locations_all = np.zeros((len(S), len(quantiles), 4))
    for f in range(0, len(S)):
        with h5py.File(paths["OUT"] + S[f], 'r') as src:
            Locations = np.array(src["Locations"])
            Locations_all[f, :, :] = Locations[0:len(quantiles), :]
    for q in range(0, len(quantiles)):

        # format point locations
        points = np.zeros((1, len(Locations_all[:, q, 0]))).astype(tuple)
        for r in range(0, len(Locations_all[:, q, 0])):
            points[0, r] = (Locations_all[r, q, 1], Locations_all[r, q, 0])

        # Create MultiPoint object
        P_q = MultiPoint(list(points[0, :]))
        schema = {
            'geometry': 'MultiPoint',
            'properties': {'quantile': 'str'}
        }

        filepath = paths["OUT"] + region + '_Locations_' + str(quantiles[q]) + '.shp'
        # Create Shapefile
        with fiona.open(filepath, 'w', 'ESRI Shapefile', schema) as c:
            c.write({
                'geometry': mapping(P_q),
                'properties': {'quantile': 'q' + str(quantiles[q])}
            })
        print("files saved:" + filepath)

# if __name__ == '__main__':
#
#     # Hard coded Output Folder Timestamp
