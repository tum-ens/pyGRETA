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

        # Initialize with the corresponding config for each country defined in folder 'configs'
        paths, param = ii.initialization(config)
              
        # Download data from Global Wind Atlas
        # im.downloadGWA(paths, param)
        
        # im.generate_landuse(paths, param)
        
        # im.generate_roads(paths, param)
        # im.generate_railways(paths, param)
        # im.generate_osm_areas(paths, param)
        # im.generate_buffered_osm_areas(paths, param)
        # im.generate_buffered_protected_areas(paths, param)
        # im.generate_airports(paths,param)
        # im.generate_subregions(paths, param)
        # im.generate_landuse(paths, param)
        # im.generate_buffered_population(paths, param)
        # im.generate_buffered_water(paths, param)
        # im.generate_buffered_wetland(paths, param)
        # im.generate_buffered_snow(paths, param)
        # im.generate_settlements(paths, param)
        # im.generate_HydroLakes(paths, param)
        # im.generate_HydroRivers(paths, param)
        
        
        # Generate input raster maps
        if os.path.isfile(paths["PA"]):
            print('Skip generate_maps_for_scope')
        else:
            im.generate_maps_for_scope(paths, param)
        # Generate buffer maps
        if os.path.isfile(paths["WINDON_PA_BUFFER"]):   # Check if there are already maps
            print('Skip generate_buffered_maps')
        else:
            im.generate_buffered_maps(paths, param)

         # Wind speed correction
        if ("WindOn" in param["technology"] or "WindOff" in param["technology"]):
            if os.path.isfile(paths["CORR_ON"]):
                print('Skip generate_wind_correction')
            else:
                cf.generate_wind_correction(paths, param)  # TODO: Into calculate full load hours?
        
        if "Biomass" in param["technology"]:
            pl.generate_biomass_production(paths, param, "Biomass")
            pl.report_biomass_potentials(paths, param, "Biomass")
            # club_biomass(paths,param)

        for tech in param["technology"]:
            print("Tech: " + tech)

            # Generate potential maps and reports
            # pl.calculate_full_load_hours(paths, param, tech)
            pl.mask_potential_maps(paths, param, tech)
            pl.weight_potential_maps(paths, param, tech)
            pl.report_potentials(paths, param, tech)

            # # Generate time series
            # find_representative_locations(paths, param, tech)
            # generate_time_series_for_representative_locations(paths, param, tech)
            #generate_time_series_for_specific_locations(paths, param, tech)

        #for tech in param["technology"]:
         #   print("Tech: " + tech)

            # Generate regression coefficients for FLH and TS model matching
            #get_regression_coefficients(paths, param, tech)

            # Generate times series for combinations of technologies and locations
            #generate_time_series_for_regions(paths, param, tech)