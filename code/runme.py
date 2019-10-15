from model_functions import *

if __name__ == "__main__":

    paths, param = initialization()
    generate_weather_files(paths, param)  # MERRA Weather data
    clean_weather_data(paths, param)  # Outlier smoothing
    generate_landsea(paths, param)  # Land and Sea
    generate_subregions(paths, param)  # Subregions
    generate_area(paths, param)  # Area Gradient
    generate_landuse(paths, param)  # Landuse
    generate_bathymetry(paths, param)  # Bathymetry
    generate_topography(paths, param)  # Topography
    generate_slope(paths, param)  # Slope
    generate_population(paths, param)  # Population
    generate_protected_areas(paths, param)  # Protected areas
    generate_buffered_population(paths, param)  # Buffered Population

    # Wind Speed correction for hub_height and ????
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
