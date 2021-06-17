from lib.correction_functions import generate_wind_correction
from lib.initialization import initialization
from lib.input_maps import generate_maps_for_scope, generate_buffered_maps
from lib.potential import calculate_full_load_hours, mask_potential_maps, weight_potential_maps, report_potentials, generate_biomass_production, club_biomass
from lib.regression import get_regression_coefficients
from lib.time_series import (
    find_representative_locations,
    generate_time_series_for_representative_locations,
    generate_time_series_for_regions,
    generate_time_series_for_specific_locations,
    generate_time_series_for_full_scope
)

config_file = 'sample.txt'

if __name__ == "__main__":

    paths, param = initialization(config_file)

    # Generate input raster maps
    # generate_maps_for_scope(paths, param)
    
    # Generate buffer maps
    # generate_buffered_maps(paths,param)
    
    

     #Wind speed correction
    #if "WindOn" in param["technology"] or "WindOff" in param["technology"]:
        #generate_wind_correction(paths, param)
        
    if "Biomass" in param["technology"]:
        #generate_biomass_production(paths, param)
        club_biomass(paths,param)

    for tech in param["technology"]:
        print("Tech: " + tech)

        # Generate potential maps and reports
        calculate_full_load_hours(paths, param, tech)
        #mask_potential_maps(paths, param, tech)
        #weight_potential_maps(paths, param, tech)
        #report_potentials(paths, param, tech)

        # Generate time series
        #generate_time_series_for_full_scope(paths, param, tech)
        #find_representative_locations(paths, param, tech)
        #generate_time_series_for_representative_locations(paths, param, tech)
        #generate_time_series_for_specific_locations(paths, param, tech)

    for tech in param["technology"]:
        print("Tech: " + tech)

        # Generate regression coefficients for FLH and TS model matching
        #get_regression_coefficients(paths, param, tech)

        # Generate times series for combinations of technologies and locations
        #generate_time_series_for_regions(paths, param, tech)
