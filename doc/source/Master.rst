Master.py
==========

Initialization
---------------

.. automodule:: Master
	:members: initialization
   
Weather files for the scope
----------------------------

.. automodule:: Master
   :noindex:
   :members: generate_weather_files, clean_weather_data
   
Rasters for the scope
----------------------

.. automodule:: Master
   :noindex:
   :members: generate_landsea, generate_subregions, generate_landuse, generate_bathymetry, generate_topography, generate_slope, generate_population, generate_protected_areas, generate_buffered_population, generate_wind_correction
   
Potential maps for the scope
-----------------------------

.. automodule:: Master
   :noindex:
   :members: calculate_FLH, masking, weighting
   
Summary report for subregions
------------------------------

.. automodule:: Master
   :noindex:
   :members: reporting
   
Time series for subregions
---------------------------

.. automodule:: Master
   :noindex:
   :members: find_locations_quantiles, generate_time_series, regression_coefficients