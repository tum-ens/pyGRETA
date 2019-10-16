Main configuration function
---------------------------

.. automodule:: config
   :members: configuration, general_settings
   
.. NOTE::
   Both *param* and *paths* will be updated in the code after running the function :mod:`config.configuration`.

.. NOTE::
   ``root`` points to the directory that contains all the inputs and outputs.
   All the paths will be defined relatively to the root, which is located in a relative position to the current folder.
   
The code differentiates between the geographic scope and the subregions of interest.
You can run the first part of the script ``runme.py`` once and save results for the whole scope, and then repeat the second part using different subregions within the scope.

.. automodule:: config
   :noindex:
   :members: scope_paths_and_parameters

.. NOTE::
   We recommend using a name tag that describes the scope of the bounding box of the regions of interest.
   For example, ``'Europe'`` and ``'Europe_without_Switzerland'`` will actually lead to the same output for the first part of the code.
   
.. NOTE::
   As of version |version|, it is possible to use different technologies in the same run, but not the same technology with different settings.

.. WARNING::
   If you intend to use the wind correction feature relying on the `Global Wind Atlas <https://globalwindatlas.info/>`_,
   it is recommended that *spatial_scope* covers **all** the countries that you are interested in,
   because the correction is done on a country-level. Also, you have to download the data from the Global Wind Atlas
   for each country that lies within the scope, even partially, and put it in the corresponding location.


User preferences
----------------

.. automodule:: config
   :noindex:
   :members: computation_parameters, resolution_parameters
   
.. NOTE::
   As of version |version|, these settings should not be changed. Only MERRA-2 data can be used in the tool.
   Its spatial resolution is 0.5° of latitudes and 0.625° of longitudes. The high resolution is 15 arcsec in both directions.

.. automodule:: config
   :noindex:
   :members: weather_data_parameters, file_saving_options, time_series_parameters, landuse_parameters, protected_areas_parameters, pv_parameters, csp_parameters, onshore_wind_parameters, offshore_wind_paramters


Paths
------

.. _path-to-gwa:

.. automodule:: config
   :noindex:
   :members: weather_input_folder, global_maps_input_paths, output_folders, weather_output_paths, local_maps_paths, irena_paths, regression_paths, emhires_input_paths, potential_output_paths, regional_analysis_output_paths