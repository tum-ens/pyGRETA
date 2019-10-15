from initialization import initialization
from input_maps import *
from correction_functions import generate_wind_correction
from potential import calculate_FLH, masking, weighting, reporting
from time_series import find_locations_quantiles, generate_time_series, generate_stratified_timeseries
from regression import regression_coefficients

if __name__ == "__main__":

    paths, param = initialization()
    generate_input_maps(paths, param)

    # Wind speed correction
    if "WindOn" in param["technology"] or "WindOff" in param["technology"]:
        generate_wind_correction(paths, param)

    for tech in param["technology"]:
        print("Tech: " + tech)
        calculate_FLH(paths, param, tech)
        masking(paths, param, tech)
        weighting(paths, param, tech)
        reporting(paths, param, tech)
        find_locations_quantiles(paths, param, tech)
        generate_time_series(paths, param, tech)

    for tech in param["technology"]:
        print("Tech: " + tech)
        regression_coefficients(paths, param, tech)
        generate_stratified_timeseries(paths, param, tech)
