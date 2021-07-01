from lib.correction_functions import generate_wind_correction
from lib.initialization import initialization
from lib.input_maps import generate_maps_for_scope
from lib.potential import calculate_full_load_hours, mask_potential_maps, weight_potential_maps, report_potentials
from lib.regression import get_regression_coefficients
from lib.time_series import (
    find_representative_locations,
    generate_time_series_for_representative_locations,
    generate_time_series_for_regions,
    generate_time_series_for_specific_locations,
)

if __name__ == "__main__":

    paths, param = initialization()

    # Generate input raster maps
    # For the tutorial, a try-except was added. With this it is not necessary to download all raw data.
    try:
        generate_maps_for_scope(paths, param)
    except:
        print("#############################################################")
        print("Raw input data is missing to run this step.")
        print("If you want to run this step, make sure, all the raw input data is downloaded and at the right place.")
        print("If not, you can also skip this step and continue with Step 2.")
        print("#############################################################")

    # Wind speed correction
    if "WindOn" in param["technology"] or "WindOff" in param["technology"]:
        generate_wind_correction(paths, param)

    for tech in param["technology"]:
        print("Tech: " + tech)

        # Generate potential maps and reports
        calculate_full_load_hours(paths, param, tech)
        mask_potential_maps(paths, param, tech)
        weight_potential_maps(paths, param, tech)
        #report_potentials(paths, param, tech)

        # Generate time series
        #find_representative_locations(paths, param, tech)
        #generate_time_series_for_representative_locations(paths, param, tech)
        #generate_time_series_for_specific_locations(paths, param, tech)

    for tech in param["technology"]:
        print("Tech: " + tech)

        # Generate regression coefficients for FLH and TS model matching
        # get_regression_coefficients(paths, param, tech)

        # Generate times series for combinations of technologies and locations
        #generate_time_series_for_regions(paths, param, tech)
