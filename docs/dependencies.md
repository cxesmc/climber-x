# Dependencies

**CLIMBER-X** is dependent on the following libraries:

- NetCDF: [NetCDF library](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)
- FFTW (ver. 3.9+)

- Python 3.x, which is only needed for automatic configuration of the Makefile
and the use of the `runme` script for job preparation and submission.
- runme: ['runme' Python library (fesmc version)](https://github.com/fesmc/runme), used to stage, run, and submit single simulations and ensembles.

Needed only if running with coupled ice sheets:

- Yelmo: [Yelmo ice sheet model](https://github.com/palma-ice/yelmo). The library will have to be compiled from the original source code.
- LIS: [Library of Iterative Solvers for Linear Systems](http://www.ssisc.org/lis/). The library will have to be compiled from the original source code.

Installation tips for each dependency can be found below.

## Installing NetCDF (preferably version 4.0 or higher)

The NetCDF library is typically available with different distributions (Linux, Mac, etc).
Along with installing `libnetcdf`, it will be necessary to install the package `libnetcdf-dev`.
Installing the NetCDF viewing program `ncview` is also recommended.

If you want to install NetCDF from source, then you must install both the
`netcdf-c` and subsequently `netcdf-fortran` libraries. The source code and
installation instructions are available from the Unidata website:

[https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)

## Install LIS, FFTW and utils

These packages could be installed individually and linked into the main directories of climber-x and yelmo. However, to ensure the right versions are used, etc., we have now made a separate repository for managing the installation of LIS and FFTW from the versions available in that repository, as well as some custom modules. This repository is managed as part of the Fast Earth System Model Community (FESMC).

Please download the code from this repository and see the README for installation instructions:
[https://github.com/fesmc/fesm-utils](https://github.com/fesmc/fesm-utils)

## Installing runme

1. Install `runme` to your system's Python installation via `pip`.

```bash
pip install https://github.com/fesmc/runme
```

That's it! Now check that system command `runme` is available by running `runme -h`.
If the command is not found, it means that the Python bin directory is not available in your `PATH`. To add it, typically something like this is needed in your .profile or .bashrc file:

```bash
PATH=${PATH}:${HOME}/.local/bin
export PATH
```
