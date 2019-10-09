User manual
===========

Installation
------------

.. NOTE:: We assume that you are familiar with `git <https://git-scm.com/downloads>`_ and `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html>`_.

First, clone the git repository in a directory of your choice using a Command Prompt window::

	$ ~\directory-of-my-choice> git clone https://github.com/tum-ens/renewable-timeseries.git

We recommend using conda and installing the environment from the file ``PyMTS37.yml`` that you can find in the repository. In the Command Prompt window, type::

	$ cd renewable-timeseries\env\
	$ conda env create -f PyMTS37.yml

Then activate the environment::

	$ conda activate PyMTS37

In the folder ``code``, you will find multiple files:

.. tabularcolumns:: |p{3cm}|p{6cm}|

.. Needs to be filled up.... 
	.. list-table::
	   :header-rows: 1

* - File
  - Description
  
* - config.py
  - used for configuration, see below.

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

.. toctree::
   :maxdepth: 3
   
   source/config


Recommended input sources
-------------------------
Weather data from MERRA-2
^^^^^^^^^^^^^^^^^^^^^^^^^^
The most important  inputs within this model are the weather time series.
These are taken from the Modern-Era Retrospective Analysis for Research and Applications, version 2 (MERRA-2),
which is the latest atmospheric reanalysis of the modern satellite era produced by NASA's Global Modeling and
Assimilation Office (GMAO) :cite:`Gelaro.2017`. The parameters taken from MERRA-2 are:

* Global Horizontal Irradiance (*GHI*): Downward shortwave radiation received by a surface horizontal to the ground
  (*SWGDN* in MERRA-2 nomenclature).
* Top of the Atmosphere Irradiance (*TOA*): Downward shortwave radiation at the top of the atmosphere
  (*SWTDN* in MERRA-2 nomenclature).
* Air temperature 2 meters above the ground (*T2M*).
* Northward wind velocity at 50 meters (*V50M*).
* Eastward wind velocity at 50 meters (*U50M*).

The *GHI* and *TOA* data are time-averaged hourly values given in W/m while *T2M* data are instantaneous
values in Kelvin. *V50M* and *U50M* are instantaneous hourly values given in m/s.

The spatial arrangement of the data consists of a global horizontal grid structure with a resolution of 576 points in
the longitudinal direction and 361 points in the latitudinal direction, resulting in pixels of 5/8° longitude 
and 1/2° latitude :cite:`MERRA2.`.

It is possible to download MERRA-2 dataset for the whole globe or just for a subset of your region of interest.
Depending on the "MERRA_coverage" parameter in config.py, the script can accept both datasets. Note that downloading 
the coverage for the whole globe is easier but will require a significant amount of space on your drive (coverage 
of the whole globe requires 13.6 Gb).

In both cases, please following these instructions to download the MERRA-2 dataset:

1. In order to download MERRA-2 data using the FTP server, you first need to create an Eathdata account (more on that on their `website <https://disc.gsfc.nasa.gov/data-access>`_).
2. Navigate to the link for the FTP sever `here <https://disc.gsfc.nasa.gov/daac-bin/FTPSubset2.pl>`_.
3. In *Data Product*, choose :math:`\texttt{tavg1\_2d\_slv\_NX}` and select the *Parameters* T2M, U50M, V50M to downaload the temperature and the wind speed datasets.
4. In *Spatial Search*, enter the coordinates of the bounding box around your region of interest or leave the default values for the whole globe. 
   To avoid problems at the edge of the MERRA-2 cells, use the following set of formulas:

   .. math::
   	   \begin{align*}
           minLat &= \left\lfloor\dfrac{s + 90 + 3/10}{1/2}\right\rfloor  \\
           maxLat &= \left\lceil\dfrac{n + 90 + 1/5}{1/2}\right\rceil  \\
           minLon &= \left\lfloor\dfrac{w + 180 + 1/16}{5/8}\right\rfloor  \\
           maxLon &= \left\lceil\dfrac{e + 180 - 1/16}{5/8}\right\rceil 
       \end{align*}
	
   where *[s n w e]* are the southern, northern, western, and eastern bounds of
   the region of interest, which you can read from the shapefile properties in
   a GIS software.
	
5. In *Temporal Order Option*, choose the year(s) of interest.
6. Leave the other fields unchanged (no time subsets, no regridding, and NetCDF4 for the output file format).
7. Repeat the steps 4-6 for the *Data Product* :math:`\texttt{tavg1\_2d\_rad\_Nx}`, for which you select the *Parameters* SWGDN and SWTDN, the surface incoming shortwave flux and the top of the atmosphere incoming shortwave flux.
8. Follow the instructions in the `website <https://disc.gsfc.nasa.gov/data-access>`_ to actually download the NetCDF4 files from the urls listed in the text files you obtained. 

If you follow these steps to download the data for the year 2015, you will obtain 730 NetCDF files, one for each day of the year and for each data product. 

Land use Map
^^^^^^^^^^^^
Another important input for this model is the land use type. 
A Land use map is useful in the sense that other parameters can be associated with different landuse types, namely:

* Urban areas
* Ross Coefficient
* Hellmann Coefficient
* Albedo
* Suitability
* Availability
* Installation Cost
* ...

By assigning a value of those parameters to each land use type, a land use raster geographically allocates 
such parameters, which play an important role in the calculations for solar power and wind speed correction. 
The spatial resolution of the land use raster, and therefore the other geographic intermediate rasters,
used for this model are defined by pixels of 1/240° longitude and 1/240° latitude.
This resolution will be used for further calculations and it will be addressed as
high resolution.

Expected outputs
----------------

Recommended workflow
--------------------
The script is designed to be modular, yet there is a recommended work flow to follow for your first run...
