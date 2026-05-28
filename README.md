# CLIMBER-X Earth System Model

[![Documentation](https://img.shields.io/badge/docs-CLIMBER--X-blue?logo=quarto&logoColor=white)](https://cxesmc.github.io/climber-x/)

CLIMBER-X is a comprehensive fast Earth System Model, designed to simulate the evolution of the Earth system on time scales ranging from decades to glacial-interglacial cycles. It includes the following components:
- semi-empirical statistical-dynamical Atmosphere model (SESAM)
- 3-D frictional-geostrophic Ocean model (GOLDSTEIN): [Edwards and Marsh 2005](https://link.springer.com/article/10.1007/s00382-004-0508-8)
- Ocean Carbon Cycle model (HAMOCC): [Ilyina et al. 2013](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012MS000178)
- Sea Ice model (SISIM)
- Land model (PALADYN): [Willeit and Ganopolski 2016](https://gmd.copernicus.org/articles/9/3817/2016/)
- Ice Sheet model (Yelmo): [Yelmo](https://github.com/palma-ice/yelmo)
- Ice Sheet model (SICOPOLIS): [SICOPOLIS](http://www.sicopolis.net/)
- Viscoelastic Lithosphere and Mantle model (VILMA)

The model is described in detail in the following articles:
- Willeit, M., Ganopolski, A., Robinson, A., and Edwards, N. R.: The Earth system model CLIMBER-X v1.0 – Part 1: 
Climate model description and validation, Geosci. Model Dev., 15, 5905–5948, [https://doi.org/10.5194/gmd-15-5905-2022](https://doi.org/10.5194/gmd-15-5905-2022), 2022. 
- Willeit, M., Ilyina, T., Liu, B., Heinze, C., Perrette, M., Heinemann, M., Dalmonech, D., Brovkin, V., Munhoven, G., Börker, J., Hartmann, J., Romero-Mujalli, G., and Ganopolski, A.: The Earth system model CLIMBER-X v1.0 – Part 2: The global carbon cycle, Geosci. Model Dev., 16, 3501–3534, [https://doi.org/10.5194/gmd-16-3501-2023](https://doi.org/10.5194/gmd-15-5905-2022), 2023.
- Willeit, M., Calov, R., Talento, S., Greve, R., Bernales, J., Klemann, V., Bagge, M., and Ganopolski, A.: Glacial inception through rapid ice area increase driven by albedo and vegetation feedbacks, Clim. Past, 20, 597–623, [https://doi.org/10.5194/cp-20-597-2024](https://doi.org/10.5194/cp-20-597-2024), 2024.

While most components of CLIMBER-X are available as open-source code, the access to some parts (ocean biogeochemistry model HAMOCC and solid Earth model Vilma) are still restricted, but the goal is to eventually publish the whole model as open source.

## COPYRIGHT

Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research (PIK)

## AUTHORS

Matteo Willeit (willeit@pik-potsdam.de)
Andrey Ganopolski (andrey@pik-potsdam.de)
Neil R. Edwards (neil.edwards@open.ac.uk)
Alexander Robinson (alexander.robinson@awi.de)
Reinhard Calov (calov@pik-potsdam.de)
Ralf Greve (greve@lowtem.hokudai.ac.jp)

## LICENSE

CLIMBER-X is free software: you can redistribute it and/or modify
it under the terms of the **GNU General Public License** as published by
the Free Software Foundation, **version 3** of the License or later. 
CLIMBER-X is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this
program. See the LICENSE file in the root directory. If not, see
https://www.gnu.org/licenses/gpl.txt

## NOTES

Outside joint collaborative agreements, there is absolutely no support in model download, 
setup, development, application or similar. 

When using a modified version of **CLIMBER-X** which is not identical to versions
in the official main repository at [https://github.com/cxesmc/climber-x](https://github.com/cxesmc/climber-x) add a suffix
to the name to allow distinguishing versions (format **CLIMBER-X-suffix**).

## Getting started

See the [CLIMBER-X documentation](https://cxesmc.github.io/climber-x/) for further details on how to get started using **CLIMBER-X**.


