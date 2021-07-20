import lib.correction_functions as cf    # old: generate_wind_correction
import lib.initialization as ii          # old: import initialization
import lib.input_maps as im              # old: generate_maps_for_scope, generate_buffered_maps
import lib.potential as pl               # old: import calculate_full_load_hours, mask_potential_maps, weight_potential_maps, report_potentials, generate_biomass_production, club_biomass
from lib.regression import get_regression_coefficients
from lib.time_series import (
    find_representative_locations,
    generate_time_series_for_representative_locations,
    generate_time_series_for_regions,
    generate_time_series_for_specific_locations
)
import os
from lib.log import logger
import logging
import shutil
import psutil

if __name__ == "__main__":

    print(psutil.virtual_memory().available)
    if psutil.virtual_memory().available > 50*10**9:
        multiprocessing = True
    else:
        multiprocessing = False
    logger.info('Multiprocessing: ' + str(multiprocessing))
    #logger.setLevel(logging.DEBUG)    # Comment out to get more information on the console

    configs = os.listdir('../configs')     # Iterate over all config files for each country in folder 'configs'
    for config in configs:      # Loop over all countries

        logger.info('Started: ' + str(config))
        # Initialize for each country with the corresponding config defined in folder 'configs'
        paths, param = ii.initialization(config)

        im.downloadGWA(paths, param)    # Download wind speed data from Global Wind Atlas

        im.generate_maps_for_scope(paths, param, multiprocessing)    # Generate input raster maps
        im.generate_buffered_maps(paths, param, multiprocessing)     # Generate buffer maps

        cf.generate_wind_correction(paths, param)  # TODO: Into calculate full load hours?

        # if "Biomass" in param["technology"]:
        #     pl.generate_biomass_production(paths, param)
        #     # club_biomass(paths,param)

        for tech in param["technology"]:
            # print("Tech: " + tech)
            logger.info("Tech: " + tech)

            # Generate potential maps and reports
            pl.calculate_full_load_hours(paths, param, tech, multiprocessing)
            pl.mask_potential_maps(paths, param, tech)
            pl.weight_potential_maps(paths, param, tech)
            pl.report_potentials(paths, param, tech)

            # Generate time series
            #find_representative_locations(paths, param, tech)
            #generate_time_series_for_representative_locations(paths, param, tech)
            #generate_time_series_for_specific_locations(paths, param, tech)

        #for tech in param["technology"]:
         #   print("Tech: " + tech)

            # Generate regression coefficients for FLH and TS model matching
            #get_regression_coefficients(paths, param, tech)

            # Generate times series for combinations of technologies and locations
            #generate_time_series_for_regions(paths, param, tech)