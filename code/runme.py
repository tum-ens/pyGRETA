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

if __name__ == "__main__":

    configs = os.listdir('configs')     # Look for config files of each country in folder 'configs'
    for config in configs:      # Loop over all countries
        print('Started: ' + str(config))

        paths, param = ii.initialization(config)   # Initialize with the corresponding config for each country defined in folder 'configs'

        # # Generate input raster maps
        im.generate_maps_for_scope(paths, param)

        # # Generate buffer maps
        im.generate_buffered_maps(paths, param)

         # #Wind speed correction
        if "WindOn" in param["technology"] or "WindOff" in param["technology"]:
            cf.generate_wind_correction(paths, param)

        if "Biomass" in param["technology"]:
            pl.generate_biomass_production(paths, param)
            # club_biomass(paths,param)

        for tech in param["technology"]:
            print("Tech: " + tech)

            # Generate potential maps and reports
            pl.calculate_full_load_hours(paths, param, tech)
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