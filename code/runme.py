import lib.correction_functions as cf
import lib.input_maps as im
import lib.potential as pl
from lib.log import logger
import initialization as ii
import os
import psutil

if __name__ == "__main__":

    # logger.setLevel(logging.DEBUG)    # Comment out to get more information on the console

    if psutil.virtual_memory().available > 50*10**9:    # Check if memory size is large enough for multiprocessing
        multiprocessing = True
    else:
        multiprocessing = False
    logger.info('Multiprocessing: ' + str(multiprocessing))

    configs = sorted(os.listdir('../configs'))
    for config in configs:       # Iterate over all config files for each country in folder 'configs'

        logger.info('Started: ' + str(config))
        paths, param = ii.initialization(config)    # Initialize for each country with the corresponding config defined in folder 'configs'

        im.downloadGWA(paths, param)    # Download wind speed data from Global Wind Atlas
        im.generate_maps_for_scope(paths, param, multiprocessing)    # Generate input raster maps

        cf.generate_wind_correction(paths, param)


        for tech in param["technology"]:
            logger.info("Tech: " + tech)
            if tech == "Biomass":
                pl.generate_biomass_production(paths, param, tech)
                pl.report_biomass_potentials(paths, param, tech)
                # pl.club_biomass(paths,param)
            else:
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