renewable-timeseries
=====================
GIS-based model for renewable energy potential and time series generation
-------------------------------------------------------------------------

:Code developers: Kais Siala, Houssame Houmy
:Documentation authors:	Kais Siala, Sergio Alejandro Huezo Rodriguez, Houssame Houmy		
:Maintainers: Kais Siala <kais.siala@tum.de>
:Organization: `Chair of Renewable and Sustainable Energy Systems`_, Technical University of Munich
:Version: |version|
:Date: |today|
:License:
 The model code is licensed under the `GNU General Public License 3.0  <http://www.gnu.org/licenses/gpl-3.0>`_.  
 This documentation is licensed under a `Creative Commons Attribution 4.0 International <http://creativecommons.org/licenses/by/4.0/>`_ license. 

Features
--------
* Generation of potential maps and time series for user-defined regions within the globe
* Modeled technologies: onshore wind, offshore wind, PV, CSP (user-defined technology characteristics)
* Use of MERRA-2 reanalysis data, with the option to detect and correct outliers
* High resolution potential taking into account the land use suitability/availability, topography, bathymetry, slope, distance to urban areas, etc.
* Statistical reports with summaries (available area, maximum capacity, maximum energy output, etc.) for each user-defined region
* Generation of several time series for each technology and region, based on user's preferences
* Possibility to combine the time series into one using linear regression to match given full-load hours and temporal fluctuations

Applications
-------------
This code is useful if:

* You want to estimate the theoretical and/or technical potential of an area, which you can define through a shapefile
* You want to obtain high resolution maps
* You want to define your own technology characteristics
* You want to generate time series for an area after excluding parts of it that are not suitable for renewable power plants
* You want to generate multiple time series for the same area (best site, upper 10%, median, lower 25%, etc.)
* You want to match historical capacity factors of countries from the IRENA database

You do not need to use the code (*but you can*) if:

* You do not need to exclude unsuitable areas - use the `Global Solar Atlas <https://globalsolaratlas.info/>`_ or `Global Wind Atlas <https://globalwindatlas.info/>`_
* You only need time series for specific points - use other webtools such as `Renewables.ninja <https://www.renewables.ninja/>`_
* You only need time series for administrative divisions (countries, NUTS-2, etc.), for which such data is readily available - see `Renewables.ninja <https://www.renewables.ninja/>`_ or `EMHIRES <https://ec.europa.eu/jrc/en/scientific-tool/emhires>`_

Changes
--------
version 1.0.0
^^^^^^^^^^^^^^
This is the initial version.


Contents
--------
User manual
^^^^^^^^^^^^^

These documents give a general overview and help you getting started from the installation to you first running model.

.. the following section contains the links to the other parts of the documentation, the 'maxdepth' component define how many sub sections should be displayed

.. toctree::
   :maxdepth: 3
   
   User_manual
	

Mathematical documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Continue here if you want to understand the theoretical conception of the model,
the logic behind the equations, and the structure of the features.

.. toctree::
   :maxdepth: 2

   theoretical

Technical documentation
^^^^^^^^^^^^^^^^^^^^^^^

Continue here if you want to understand in detail the model implementation.

.. toctree::
   :maxdepth: 2
   
   implementation


Dependencies
------------


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
