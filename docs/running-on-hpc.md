# Running on an HPC

## Notes for specific systems

### Running at PIK on HPC2024 (foote)

The following modules have to be loaded in order to compile and run the model.
For convenience you can also add those commands to your `.profile` file in your home directory.

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

When installing `fesm-utils` (see [Dependencies](dependencies.md)) use the `pik` script:

```bash
./install_pik.sh ifx
```

### Running at AWI on albedo

Load the following modules in your `.bashrc` file in your home directory.

```bash
    module load intel-oneapi-compilers/2024.0.0
    module load netcdf-c/4.8.1-openmpi4.1.3-oneapi2022.1.0
    module load netcdf-fortran/4.5.4-oneapi2022.1.0
    module load udunits/2.2.28
    module load ncview/2.1.8
    module load cdo/2.2.0
    module load python/3.10.4
```

When installing `fesm-utils` (see [Dependencies](dependencies.md)) use the `awi` script (which is actually a link to the `dkrz` script since they work the same way):

```bash
./install_awi.sh ifx
```

### Running at DKRZ on levante

Load the following modules in your `.bashrc` file in your home directory.

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

When installing `fesm-utils` (see [Dependencies](dependencies.md)) use the `dkrz` script:

```bash
./install_dkrz.sh ifx
```

## Get the code

### CLIMBER-X climate model

```bash

### Download the CLIMBER-X code ###

# Clone repository
git clone https://github.com/cxesmc/climber-x.git
git clone git@github.com:cxesmc/climber-x.git # via ssh

# Enter directory 
cd climber-x

# Run configuration script
python config.py config/pik_hpc2024_ifx   # Or config file for your system

# Clone input file directory
git clone https://gitlab.pik-potsdam.de/cxesmc/climber-x-input.git input
git clone git@gitlab.pik-potsdam.de:cxesmc/climber-x-input.git input    # via ssh

# Step back and clone and install external libraries repository
cd ..
git clone git@github.com:fesmc/fesm-utils.git
cd fesm-utils
./install_pik.sh ifx   # Use install_dkrz.sh as needed
FESMUSRC=$PWD
cd ../climber-x/
ln -s $FESMUSRC ./src/utils/

# Download and configure coordinates
cd src/utils/
git clone git@github.com:fesmc/coordinates.git
cd coordinates
python3 config.py config/pik_hpc2024_ifx   # Or config file for your system
cd ../../..   # Return to climber-x parent directory

# Install other external utils library
cd src/utils/fesm-utils/utils
python config.py config/pik_hpc2024_ifx  # replace with config file for your system
make clean
make fesmutils-static openmp=0	# serial version	
make fesmutils-static openmp=1	# parallel version
cd ../../../.. # Return to climber-x parent directory

### Compile and run ###

# Compile the climate model 
make cleanall
make climber-clim

# Set up your `runme` config file for your system
cp .runme/runme_config .runme_config
# - Edit hpc and account name to match your settings

# Make sure to install the `runme` package too
pip install https://github.com/fesmc/runme 

# Run a pre-industrial equilibrium climate-only test simulation
runme -rs -q short --omp 32 -o output/clim
```

### CLIMBER-X climate and carbon cycle model

If you would also like to run CLIMBER-X with an interactive carbon cycle, then the **HAMOCC**
ocean biogeochemistry (`bgc`) code must also be downloaded:

```bash
# bgc
cd src/
git clone git@github.com:cxesmc/bgc.git
cd bgc/
git submodule update --init --recursive     # for submodule M4AGO 
cd ../..
```

Since the HAMOCC model code is not open source, the `bgc` repository is private at the moment and
you need to be given permission in order to access it. HAMOCC is covered by the Max Planck Institute for
Meteorology software licence agreement as part of the MPI-ESM ([https://code.mpimet.mpg.de/attachments/download/26986/MPI-ESM_SLA_v3.4.pdf](https://code.mpimet.mpg.de/attachments/download/26986/MPI-ESM_SLA_v3.4.pdf)).
A pre-requisite to access the `bgc` repository is therefore that you agree to the MPI-ESM license
by following the steps outlined here: [https://code.mpimet.mpg.de/projects/mpi-esm-license](https://code.mpimet.mpg.de/projects/mpi-esm-license).
Once you have done so, send an email to [Matteo Willeit](mailto:matteo.willeit@gmail.com?subject=[GitHub]%20bgc%20source%20code) and you will be granted permission to access the `bgc` repository.

```bash
# Compile the climate and carbon cycle model 
make clean
make climber-clim-bgc

# Run a pre-industrial equilibrium simulation with ocean biogeochemistry
runme -rs -q short --omp 16 -o output/clim-bgc -p ctl.flag_bgc=T
```

### CLIMBER-X climate and ice sheet model

If you would also like to run with an interactive ice sheet, the **Yelmo** ice-sheet code
must be downloaded and configured and the solid Earth model **VILMA** libraries must be
downloaded before compiling:

```bash
# yelmo
cd src
git clone git@github.com:palma-ice/yelmo.git
cd yelmo
git checkout climber-x                     # Get climber-x branch
python3 config.py config/pik_hpc2024_ifx   # Or config file for your system
ln -s $FESMUSRC .              # Link absolute path
cd ../..            # Return to climber-x parent directory

# vilma
cd src/
git clone git@github.com:cxesmc/vilma.git  # private repository, premission needed
cd ..
```

Since the VILMA model code is not open source, the `vilma` repository is private at the moment and you need to be given permission in order to access it. Please send an email to [Matteo Willeit and Volker Klemann](mailto:matteo.willeit@gmail.com,volker.klemann@gfz.de?subject=[GitHub]%20VILMA%20access) and you will be granted permission to access the `vilma` repository.

```bash
# Compile the climate and ice sheet model
make clean
make climber-clim-ice

# Run pre-industrial equilibrium simulation with interactive Greenland ice sheet
runme -rs -q short --omp 16 -o output/clim-ice -p ctl.flag_ice=T ctl.flag_geo=T ctl.flag_smb=T ctl.flag_imo=T ctl.ice_model_name=yelmo ctl.ice_domain_name=GRL-16KM
```

### Fully coupled CLIMBER-X configuration

If you have followed all steps above you will also be ready to run fully coupled simulations:

```bash
# Compile the fully coupled model
make clean
make climber-clim-bgc-ice  # or equivalently make climber

# Run pre-industrial equilibrium simulation with ocean biogeochemistry and interactive Greenland ice sheet
runme -s -q short --omp 16 -o output/clim-bgc-ice -p ctl.flag_bgc=T ctl.flag_ice=T ctl.flag_geo=T ctl.flag_smb=T ctl.flag_imo=T ctl.ice_model_name=yelmo ctl.ice_domain_name=GRL-16KM
```
