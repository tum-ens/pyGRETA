# renewable-timeseries

This is a python script that generates high-resolution potential maps and time series for user-defined regions within the globe.

[![Documentation Status](https://readthedocs.org/projects/renewable-timeseries/badge/?version=latest)](http://renewable-timeseries.readthedocs.io/en/latest/?badge=latest)

## Features

* Generation of potential maps and time series for user-defined regions within the globe
* Modeled technologies: onshore wind, offshore wind, PV, CSP (user-defined technology characteristics)
* Use of MERRA-2 reanalysis data, with the option to detect and correct outliers
* High resolution potential taking into account the land use suitability/availability, topography, bathymetry, slope, distance to urban areas, etc.
* Statistical reports with summaries (available area, maximum capacity, maximum energy output, etc.) for each user-defined region
* Generation of several time series for each technology and region, based on user's preferences
* Possibility to combine the time series into one using linear regression to match given full-load hours and temporal fluctuations

## Applications

This code is useful if:

* You want to estimate the theoretical and/or technical potential of an area, which you can define through a shapefile
* You want to obtain high resolution maps
* You want to define your own technology characteristics
* You want to generate time series for an area after excluding parts of it that are not suitable for renewable power plants
* You want to generate multiple time series for the same area (best site, upper 10%, median, lower 25%, etc.)
* You want to match historical capacity factors of countries from the IRENA database

You do not need to use the code (*but you can*) if:

* You do not need to exclude unsuitable areas - use the [Global Solar Atlas](https://globalsolaratlas.info/) or [Global Wind Atlas](https://globalwindatlas.info/)
* You only need time series for specific points - use other webtools such as [Renewables.ninja](https://www.renewables.ninja/)
* You only need time series for administrative divisions (countries, NUTS-2, etc.), for which such data is readily available - see [Renewables.ninja](https://www.renewables.ninja/) or [EMHIRES](https://ec.europa.eu/jrc/en/scientific-tool/emhires)

