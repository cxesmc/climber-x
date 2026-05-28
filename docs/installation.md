# Installing CLIMBER-X

Here you can find the steps needed to install and build **CLIMBER-X**.

There are currently four different flavors of **CLIMBER-X** that can be set up:

- `climber-clim`: minimal climate model configuration with atmosphere, ocean, sea ice and land (including dynamic vegetation)
- `climber-clim-bgc`: coupled climate-carbon cycle model configuration; clim plus ocean biogeochemistry
- `climber-clim-ice`: coupled climate-ice sheet model configuration; clim plus ice sheets
- `climber-clim-bgc-ice`: fully coupled model configuration; clim plus ocean biogeochemistry and ice sheets

## System dependencies

CLIMBER-X is configured and built with [`configme`](https://github.com/fesmc/configme),
which clones, configures, links, and builds the whole stack for you (see
[Quick start](#quick-start) below). The only dependencies you must install
yourself are:

- **netCDF** (C and Fortran libraries) — see [Dependencies](dependencies.md#installing-netcdf-preferably-version-4.0-or-higher) for installation tips.
- **CDO** ([Climate Data Operators](https://code.mpimet.mpg.de/projects/cdo/)) — only needed to (re)generate maps that transform between coordinate grids.

Everything else — `coordinates`, the `fesm-utils` libraries (LIS + FFTW + utils),
`yelmo`, and `runme` — is managed by `configme`. For the full dependency list see
[Dependencies](dependencies.md).

## Install configme (one time)

`configme` is a small Python tool that detects your netCDF installation,
configures every package in the stack for your machine and compiler, and
clones/links/builds the whole thing with one command. It is installed once,
globally, and provides the `configme` command on your `PATH`:

```bash
pip install git+https://github.com/fesmc/configme
```

To upgrade it later, add `--upgrade` to the same command. If the `configme`
command is not found afterwards, your Python user bin directory is probably not
on your `PATH`; add it in your `~/.bashrc` / `~/.zshrc`:

```bash
export PATH="${PATH}:${HOME}/.local/bin"
```

## Quick start

With `configme` installed, build the whole CLIMBER-X stack with a single command
from the directory where you want the checkout to live:

```bash
configme install climber-x
```

This clones CLIMBER-X and its component repositories (`fesm-utils`,
`coordinates`, and `yelmo` on its `climber-x` branch), configures each for your
machine and compiler, links them into the CLIMBER-X directory, builds
`fesm-utils` (LIS + FFTW + utils, which can take 10-30 min), installs `runme`,
creates a `.runme_config`, and clones the large input-data repository into
`input/`. If `configme` can detect your machine from the hostname it does so,
otherwise it prompts you.

It also *attempts* to clone the private `bgc` and `vilma` components (needed for
the `bgc` and `ice` flavors). If you do not have access these are skipped with a
note rather than failing the install — see [Private components](#private-components-bgc-and-vilma).

Common options:

```bash
configme install climber-x -m pik_hpc2024 -c ifx   # pick the machine + compiler explicitly
configme install climber-x -d https                 # clone over HTTPS (no GitHub SSH key needed)
configme install climber-x --dir ~/models/climber-x # put the checkout here instead of ./climber-x
configme install climber-x --overwrite              # re-clone over an existing checkout
configme install climber-x --build-deps             # rebuild dependency packages without prompting
```

Run `configme list` for the supported machines and compilers, and
`configme --help` for the full command surface. The exact
clone/configure/link/build commands `configme install climber-x` runs for you
are recorded in a `.install.sh` script in the checkout, and are also shown for
context on the [configme install details](configme-install-details.md) page —
these are for reference only; `configme install` is the recommended path and you
do not need to run them by hand.

## Private components (bgc and vilma)

Two components are needed only for the carbon-cycle and ice-sheet flavors and
live in private repositories. `configme install climber-x` attempts to clone
them on every run, but skips them with a summary note if you do not yet have
access.

- **`bgc`** (HAMOCC ocean biogeochemistry — required for `climber-clim-bgc`).
  Since the HAMOCC model code is not open source, the `bgc` repository is private
  and you need to be given permission to access it. HAMOCC is covered by the Max
  Planck Institute for Meteorology software licence agreement as part of the
  MPI-ESM ([MPI-ESM_SLA_v3.4.pdf](https://code.mpimet.mpg.de/attachments/download/26986/MPI-ESM_SLA_v3.4.pdf)).
  A pre-requisite to access the `bgc` repository is therefore that you agree to
  the MPI-ESM license by following the steps outlined here:
  [https://code.mpimet.mpg.de/projects/mpi-esm-license](https://code.mpimet.mpg.de/projects/mpi-esm-license).
  Once you have done so, send an email to
  [Matteo Willeit](mailto:matteo.willeit@gmail.com?subject=[GitHub]%20bgc%20source%20code)
  and you will be granted permission.

- **`vilma`** (VILMA solid-Earth model — required for `climber-clim-ice`).
  Since the VILMA model code is not open source, the `vilma` repository is
  private and you need to be given permission to access it. Please send an email
  to [Matteo Willeit and Volker Klemann](mailto:matteo.willeit@gmail.com,volker.klemann@gfz-potsdam.de?subject=[GitHub]%20VILMA%20access)
  and you will be granted permission.

After being granted access, re-run `configme install climber-x` (optionally with
`--overwrite`) to clone the components you can now reach.

## Compile and run

`configme install climber-x` configures the Makefiles but leaves the choice of
flavor to you. From the `climber-x` directory, compile the executable you want
and run a test simulation with `runme`. The executable `climber.x` is produced
in the main directory.

### Climate model

The climate-only version `climber-clim` corresponds to the version described by
Willeit et al. (2022). It requires neither the private `bgc`/`vilma` code nor the
LIS library.

```bash
make clean
make climber-clim

# Run a pre-industrial equilibrium climate-only test simulation
runme -rs -q short --omp 32 -o output/clim
```

To compile with debug flags enabled use `make climber-clim debug=1`. By default
the model is compiled with OpenMP; `make climber-clim openmp=0` builds the serial
version (which should typically not be used).

### Climate and carbon cycle model

Requires the private `bgc` component (see [Private components](#private-components-bgc-and-vilma)).

```bash
make clean
make climber-clim-bgc

# Run a pre-industrial equilibrium simulation with ocean biogeochemistry
runme -rs -q short --omp 16 -o output/clim-bgc -p ctl.flag_bgc=T
```

### Climate and ice sheet model

Uses the `yelmo` ice-sheet model (cloned by `configme`) and requires the private
`vilma` solid-Earth component (see [Private components](#private-components-bgc-and-vilma)).

```bash
make clean
make climber-clim-ice

# Run a pre-industrial equilibrium simulation with an interactive Greenland ice sheet
runme -rs -q short --omp 16 -o output/clim-ice -p ctl.flag_ice=T ctl.flag_geo=T ctl.flag_smb=T ctl.flag_imo=T ctl.ice_model_name=yelmo ctl.ice_domain_name=GRL-16KM
```

### Fully coupled configuration

With both `bgc` and `vilma` available you can build the fully coupled model:

```bash
make clean
make climber-clim-bgc-ice  # or equivalently: make climber

# Run a pre-industrial equilibrium simulation with biogeochemistry and an interactive Greenland ice sheet
runme -s -q short --omp 16 -o output/clim-bgc-ice -p ctl.flag_bgc=T ctl.flag_ice=T ctl.flag_geo=T ctl.flag_smb=T ctl.flag_imo=T ctl.ice_model_name=yelmo ctl.ice_domain_name=GRL-16KM
```

See the [runme notes](runme-notes.md) to learn how `runme` stages and submits
simulations, or use the benchmark script `run_bench.sh` in the main repository.

## Notes for specific systems

On HPC systems, load the netCDF and compiler modules before running `configme
install` and `make`. For convenience you can add these commands to your
`.profile`/`.bashrc`.

### Running at PIK on HPC2024 (foote)

```bash
module purge
module use /p/system/modulefiles/compiler \
           /p/system/modulefiles/gpu \
           /p/system/modulefiles/libraries \
           /p/system/modulefiles/parallel \
           /p/system/modulefiles/tools

module load intel/oneAPI/2024.0.0
module load netcdf-c/4.9.2
module load netcdf-fortran-intel/4.6.1
module load udunits/2.2.28
module load ncview/2.1.10
module load cdo/2.4.2
```

Then install with `configme install climber-x -m pik_hpc2024 -c ifx`.

### Running at AWI on albedo

```bash
module load intel-oneapi-compilers/2024.0.0
module load netcdf-c/4.8.1-openmpi4.1.3-oneapi2022.1.0
module load netcdf-fortran/4.5.4-oneapi2022.1.0
module load udunits/2.2.28
module load ncview/2.1.8
module load cdo/2.2.0
module load python/3.10.4
```

Then install with `configme install climber-x -m awi_albedo -c ifx`.

### Running at DKRZ on levante

```bash
# Tools
module load cdo/2.4.0-gcc-11.2.0
module load esmvaltool/2.5.0
module load ncview/2.1.8-gcc-11.2.0
module load git/2.43.3-gcc-11.2.0
module load python3/2023.01-gcc-11.2.0

# Compilers and libs
module load intel-oneapi-compilers/2023.2.1-gcc-11.2.0
module load netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0
```

Then install with `configme install climber-x -m dkrz_levante -c ifx`.
