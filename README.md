<div align="left">
<img src="doc\img\pyGRETA_logo.png" alt="pyGRETA_logo" width="160px">
</div>

[![Documentation Status](https://readthedocs.org/projects/pygreta/badge/?version=latest)](http://pyGRETA.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/174577484.svg)](https://zenodo.org/badge/latestdoi/174577484)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![All Contributors](https://img.shields.io/badge/all_contributors-8-orange.svg?style=flat-square)](#contributors)

**py**thon **G**enerator of **RE**newable **T**ime series and m**A**ps: a tool that generates high-resolution potential maps and time series for user-defined regions within the globe.

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

## Outputs

Potential maps for solar PV and onshore wind in Australia, using weather data for 2015:
<div align="center">
<img src="doc\img\FLH_solar.png" alt="FLH_solar_PV_Australia_2015" width="400px"><img src="doc\img\FLH_wind.png" alt="FLH_wind_onshore_Australia_2015" width="400px">
</div>
<div align="center">
<img src="doc\img\AustraliaQ50WindvsSolar.png" alt="Australia_PV_wo_quant" width="100%">
</div>
 
## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/kais-siala"><img src="https://avatars2.githubusercontent.com/u/21306297?v=4?s=100" width="100px;" alt=""/><br /><sub><b>kais-siala</b></sub></a><br /><a href="#question-kais-siala" title="Answering Questions">💬</a> <a href="https://github.com/tum-ens/pyGRETA/issues?q=author%3Akais-siala" title="Bug reports">🐛</a> <a href="https://github.com/tum-ens/pyGRETA/commits?author=kais-siala" title="Code">💻</a> <a href="https://github.com/tum-ens/pyGRETA/commits?author=kais-siala" title="Documentation">📖</a> <a href="#ideas-kais-siala" title="Ideas, Planning, & Feedback">🤔</a> <a href="#maintenance-kais-siala" title="Maintenance">🚧</a> <a href="https://github.com/tum-ens/pyGRETA/pulls?q=is%3Apr+reviewed-by%3Akais-siala" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/tum-ens/pyGRETA/commits?author=kais-siala" title="Tests">⚠️</a> <a href="#talk-kais-siala" title="Talks">📢</a></td>
    <td align="center"><a href="https://github.com/HoussameH"><img src="https://avatars2.githubusercontent.com/u/48953960?v=4?s=100" width="100px;" alt=""/><br /><sub><b>HoussameH</b></sub></a><br /><a href="#question-HoussameH" title="Answering Questions">💬</a> <a href="https://github.com/tum-ens/pyGRETA/commits?author=HoussameH" title="Code">💻</a> <a href="https://github.com/tum-ens/pyGRETA/commits?author=HoussameH" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/pgrimaud"><img src="https://avatars1.githubusercontent.com/u/1866496?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Pierre Grimaud</b></sub></a><br /><a href="https://github.com/tum-ens/pyGRETA/issues?q=author%3Apgrimaud" title="Bug reports">🐛</a></td>
    <td align="center"><a href="https://github.com/thushara2020"><img src="https://avatars.githubusercontent.com/u/75068469?v=4?s=100" width="100px;" alt=""/><br /><sub><b>thushara2020</b></sub></a><br /><a href="https://github.com/tum-ens/pyGRETA/pulls?q=is%3Apr+reviewed-by%3Athushara2020" title="Reviewed Pull Requests">👀</a></td>
    <td align="center"><a href="https://github.com/lodersky"><img src="https://avatars.githubusercontent.com/u/36160124?v=4?s=100" width="100px;" alt=""/><br /><sub><b>lodersky</b></sub></a><br /><a href="https://github.com/tum-ens/pyGRETA/commits?author=lodersky" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/sonercandas"><img src="https://avatars.githubusercontent.com/u/17904824?v=4?s=100" width="100px;" alt=""/><br /><sub><b>sonercandas</b></sub></a><br /><a href="https://github.com/tum-ens/pyGRETA/commits?author=sonercandas" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/patrick-buchenberg"><img src="https://avatars.githubusercontent.com/u/84960180?v=4?s=100" width="100px;" alt=""/><br /><sub><b>patrick-buchenberg</b></sub></a><br /><a href="#platform-patrick-buchenberg" title="Packaging/porting to new platform">📦</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/molarana"><img src="https://avatars.githubusercontent.com/u/19924540?v=4?s=100" width="100px;" alt=""/><br /><sub><b>molarana</b></sub></a><br /><a href="#design-molarana" title="Design">🎨</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!

## Please cite as:

Kais Siala, & Houssame Houmy. (2020, June 1). tum-ens/pyGRETA: python Generator of REnewable Time series and mAps (Version v1.1.0). Zenodo. http://doi.org/10.5281/zenodo.3872068
