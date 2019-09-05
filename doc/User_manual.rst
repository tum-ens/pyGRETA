User manual
===========

Installation
------------

.. NOTE:: We assume that you are familiar with *git* and *conda*.

First, clone the git repository in a directory of your choice using a Command Prompt window::

	$ ~\directory-of-my-choice> git clone https://github.com/tum-ens/renewable-timeseries.git

We recommend using conda and installing the environment from the file ``PyMTS37.yml`` that you can find in the repository. In the Command Prompt window, type::

	$ cd renewable-timeseries/
	$ conda env create -f PyMTS37.yml

Then activate the environment::

	$ conda activate PyMTS37

In the folder ``code``, you will find multiple files:

+-----------------------+---------------------------------------------------------------------------+
| File                  | Description                                                               |
+=======================+===========================================================================+
| config.py             | used for configuration, see below.                                        |
+-----------------------+---------------------------------------------------------------------------+
| Master.py             | main file, which will be run later using ``python Master.py``.            |
+-----------------------+---------------------------------------------------------------------------+
| data_functions.py     | contains helping functions for data formatting, filtering, reshaping, etc.|
+-----------------------+---------------------------------------------------------------------------+
| model_functions.py    | contains helping functions for the modeling.                              |
+-----------------------+---------------------------------------------------------------------------+
| util.py               | contains minor helping functions and the necessary python libraries to be |
|                       | imported.                                                                 |
+-----------------------+---------------------------------------------------------------------------+

config.py                                                                                           
---------
This file contains the user preferences, the links to the input files, and the paths where the outputs should be saved.
The paths are initialized in a way that follows a particular folder hierarchy. However, you can change the hierarchy as you wish.

User preferences
^^^^^^^^^^^^^^^^^

General parameters
++++++++++++++++++
The first lines of ``config.py`` start with initializing some general parameters::

	$ param = {}
	$ paths = {}
	$ fs = os.path.sep
	$ current_folder = os.path.dirname(os.path.abspath(__file__))
	$ ...
	$ root = str(Path(current_folder).parent.parent) + "Database_KS" + fs

where ``param`` is the dictionary that will include all the user preferences, and ``paths`` the dictionary that will include all the paths for inputs and outputs.
Both will be updated in the code after running the function :func: initialization() in ``Master.py`` to include other parameters as well.
``fs`` stands for the system-specific file separator, and ``root`` points out to the directory that contains all the inputs and outputs.
All the paths will be defined relatively to the root, which is located in relative position to ``current_folder``.	

The code will generate some general outputs, that only depend on the geographic scope and year, and others that are user-specific.
For the latter, they will be saved in a separate subfolder, which is named using a timestamp, unless you overwrite the name::

	$ # Custom timestamp
	$ timestamp = str(datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
	$ timestamp = 'custom-name' # Custom name, can be commented out
	
Scope
++++++

Now come the most important settings, which are saved in ``param``. The key *region* is a name tag for the geographic scope of the output.
We recommend using a name tag that describes the scope of the bounding box of the regions of interest.
For example, 'Europe' and 'Europe_without_Switzerland' will actually lead to the same output for the first part of the code.
This is because of the flexibility offered by the tool. It differentiates between the geographic scope and the subregions of interest.
You can run the first part of the script once and save results for the whole scope, and then repeat the second part using different subregions within the scope.
This distinction will be highlighted again when we introduce the paths of the regions.
*year* defines the year of the weather data, and *technology* the list of technologies that you are interested in.
Currently, four technologies are defined: onshore wind ``WindOn``, offshore wind ``WindOff``, photovoltaics ``PV``, concentrated solar power ``CSP``.
As of version |version|, it is possible to use different technologies in the same run, but not the same technology with different settings::

	$ param["region"] = 'Australia' # Name of the spatial scope, define path to shapefile below!
	$ param["year"] = 2015
	$ param["technology"] = ['WindOn', 'PV']  # ['PV', 'CSP', 'WindOn', 'WindOff']
	
Computation
++++++++++++

Some modules in ``Master.py`` allow parallel processing. The key *nproc*, which takes an integer as a value, limits the number of parallel processes.
*CPU_limit* is a boolean parameter ?????::

	$ param["nproc"] = 10
	$ param["CPU_limit"] = True

Data resolution
++++++++++++++++

As of version |version|, only MERRA-2 data can be used in the tool. Its spatial resolution is defined in *res_weather* with an array of numbers.
The first number is the resolution in the vertical dimension (0.5° of latitudes), and the second is for the horizontal dimension (0.625° of longitudes).
*res_desired* is a similar array characterizing the high resolution of the results, which is 15 arcsec in both directions::

	$ param["res_weather"] = np.array([1 / 2, 5 / 8])
	$ param["res_desired"] = np.array([1 / 240, 1 / 240])
	
.. NOTE:: As of version |version|, these settings should not be changed.

Weather
++++++++

In *MERRA_coverage*, you can set the spatial coverage of the MERRA-2 files that you would like to use.
In case you have downloaded MERRA-2 data for the whole world for ::

	$ param["MERRA_coverage"] = 'World'
	$ param["MERRA_correction"] = 0.35
	
	
Mask / Weight
++++++++++++++
	$ param["savetiff"] = 1  # Save geotiff files of mask and weight rasters

Reporting
++++++++++

	$ param["report_sampling"] = 100

Time series
++++++++++++

	$ param["quantiles"] = np.array([100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0])

Regression
+++++++++++
	$ regression = {
	$ 	 "solver": 'gurobi',
	$ 	 "hub_heights": [80],
	$ 	 "orientations": []}
	$ param["regression"] = regression
	


Paths
^^^^^^
Inputs
+++++++
The first paths to be defined are those of the *spatial_scope* and of the *subregions*. Both should be shapefiles of polygons or multipolygons.
For the *spatial_scope*, only the bounding box around all the features matters.
In case of Europe, whether a shapefile of Europe as one multipolygon, or as a set of multiple features (countries, states, etc.) is used, does not make a difference.
Potential maps (theoretical and technical) are generated for the whole scope of the bounding box.
For the *subregions*, the shapes of the individual features matter, but not their scope.
For each individual feature that lies within the scope, you can generate a summary report and time series.
The shapefile of *subregions* does not have to have the same bounding box as *spatial_scope*.
In case it is larger, features that lie completely outside the scope will be ignored, whereas those that lie partly inside it will be cropped using the bounding box
of *spatial_scope*. In case it is smaller, all features are used with no modification, such as in the following example, where ``"Western_Australia.shp"`` is 
completely within ``"Australia.shp"``::

	$ # Shapefiles
	$ PathTemp = root + "02 Shapefiles for regions" + fs + "User-defined" + fs
	$ paths["spatial_scope"] = PathTemp + "Australia.shp"
	$ paths["subregions"] = PathTemp + "Western_Australia.shp"
	
.. NOTE:: If you intend to use the wind correction feature relying on the `Global Wind Atlas <https://globalwindatlas.info/>`_, it is recommended that the *spatial_scope* covers **all** the countries that you are interested in, because the correction is done on a country-level. Also, you have to download the data from the Global Wind Atlas for each country that lies within the scope, even partially, and put it in the corresponding location. More on that when describing the :ref:`path to GWA <path-to-gwa>`.

.. _path-to-gwa:

test

Links to the main inputs

Outputs
++++++++

Expected output explanation

Recommended input sources
^^^^^^^^^^^^^^^^^^^^^^^^^^

Recommended workflow
^^^^^^^^^^^^^^^^^^^^^
The script is designed to be modular, yet there is a recommended work flow to follow for your first run...
