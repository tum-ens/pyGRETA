Master.py
==========

Initialization
---------------

.. automodule:: Master
   :members: initialization
   
Weather files for the scope
----------------------------

.. autofunction:: Master.generate_weather_files
.. autofunction:: Master.clean_weather_data
   
Rasters for the scope
----------------------

.. autofunction:: Master.generate_landsea
.. autofunction:: Master.generate_subregions
.. autofunction:: Master.generate_landuse
.. autofunction:: Master.generate_bathymetry
.. autofunction:: Master.generate_topography
.. autofunction:: Master.generate_slope
.. autofunction:: Master.generate_area
.. autofunction:: Master.generate_population
.. autofunction:: Master.generate_protected_areas
.. autofunction:: Master.generate_buffered_population
.. autofunction:: Master.generate_wind_correction
   
Potential maps for the scope
-----------------------------

.. autofunction:: Master.calculate_FLH
.. autofunction:: Master.masking
.. autofunction:: Master.weighting
   
Summary report for subregions
------------------------------

.. autofunction:: Master.reporting
   
Time series for subregions
---------------------------

.. autofunction:: Master.find_locations_quantiles
.. autofunction:: Master.generate_time_series
.. autofunction:: Master.regression_coefficients
