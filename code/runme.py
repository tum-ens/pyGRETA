from correction_functions import generate_wind_correction
from initialization import initialization
from input_maps import *
from potential import calculate_FLH, masking, weighting, reporting
from regression import regression_coefficients
from time_series import find_locations_quantiles, generate_time_series, generate_stratified_timeseries

if __name__ == "__main__":

    paths, param = initialization()

    # Generate Input raster maps
    #############################
    generate_input_maps(paths, param)

    # Wind speed correction
    if "WindOn" in param["technology"] or "WindOff" in param["technology"]:
        generate_wind_correction(paths, param)

    for tech in param["technology"]:
        print("Tech: " + tech)

        # Generate Potential Maps and Reports
        #####################################
        calculate_FLH(paths, param, tech)
        masking(paths, param, tech)
        weighting(paths, param, tech)
        reporting(paths, param, tech)

        # Generate Timeseries
        #####################
        find_locations_quantiles(paths, param, tech)
        generate_time_series(paths, param, tech)

    for tech in param["technology"]:
        print("Tech: " + tech)

        # Generate regression coefficients for FLH and TS model matching
        ################################################################
        regression_coefficients(paths, param, tech)

        # Generate Stratified timeseries
        ################################
        generate_stratified_timeseries(paths, param, tech)
