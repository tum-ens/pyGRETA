=======
pyGRETA
=======
.. figure:: img/pyGRETA_logo.png
   :width: 20%
   :align: left
   :alt: pyGRETA_logo

.. |br| raw:: html

   <br />

|br|

|br|
--------------------------------------------------
python Generator of REnewable Time series and mAps 
--------------------------------------------------

:Code developers: Kais Siala, Houssame Houmy
:Documentation authors:	Kais Siala, Houssame Houmy, Sergio Alejandro Huezo Rodriguez		
:Maintainers: Kais Siala <kais.siala@tum.de>
:Organization: `Chair of Renewable and Sustainable Energy Systems <http://www.ens.ei.tum.de/en/homepage/>`_, Technical University of Munich
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
version 1.1.0
^^^^^^^^^^^^^^

* Introducing flexibility for resolutions.
* Adding logo for the tool.
* Minor fixes and improvements.

version 1.0.1
^^^^^^^^^^^^^^

* Fixed the syntax of the code in the PV module tracking (lib.physical_models).
* Edited the formatting of the PDF documentation.
* Edited the list of references.

version 1.0.0
^^^^^^^^^^^^^^
This is the initial version.


Contents
--------
User manual
^^^^^^^^^^^^^

These documents give a general overview and help you get started from the installation to your first running model.

.. the following section contains the links to the other parts of the documentation, the 'maxdepth' component define how many sub sections should be displayed

.. toctree::
   :maxdepth: 3
   
   user_manual
	

Mathematical documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Continue here if you want to understand the concept of the model.

.. toctree::
   :maxdepth: 2

   theory

Technical documentation
^^^^^^^^^^^^^^^^^^^^^^^

Continue here if you want to understand in detail the model implementation.

.. toctree::
   :maxdepth: 2
   
   implementation


Dependencies
------------
A list of the used libraries is available in the environment file:

.. literalinclude:: ../env/ren_ts.yml

Bibliography	
------------

.. toctree::
   :maxdepth: 1
   
   zref

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
