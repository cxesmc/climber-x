# Pre-industrial spinup

Instructions to perform a pre-industrial spinup simulation.

## Climate only

The default parameters in the control namelist (`nml/control.nml`) are set to pre-industrial conditions. 
Running a climate-only pre-industrial equilibrium spinup therefore doesn't require any specific changes: 

```bash
runme -rs -q short -w 24:00:00 --omp 32 -o output/pi_spinup \
-p ctl.nyears=5000 ctl.i_write_restart=1
```

This simulation corresponds to the standard CMIP experiment **piControl**.

The length of the simulation is set to 5000 years (which is also the default), which is generally enough to reach an equilibrium of the climate system.
`ctl.i_write_restart=1` enforces a regular writing of the restart files every 1000 years. By default the restart files are written out only at the end of the simulation.

The pre-industrial equilibrium state of a given model release is generally included in the `restart/pi` directory, which is also the default directory from where restart files are read (`ctl.restart_in_dir`).

## With carbon cycle (*closed*)

To spin up the model including the carbon cycle, the ocean biogeochemistry module needs to be enabled (`ctl.flag_bgc=T`):

```bash
runme -rs -q medium -w 60:00:00 --omp 32 -o output/pi_spinup_cc -p \
ctl.nyears=10000 ctl.flag_bgc=T ctl.bgc_restart=F ctl.i_write_restart=1
```

This simulation should generally be run a bit longer than the climate-only spinup, in particular to allow equilibration of the slow soil carbon pools. The `ctl.bgc_restart=F` specifies that the biogeochemistry module is not initialized from a previously computed state, but uses present-day observation of the tracer concentrations instead.
The land carbon cycle is active also in the climate-only setup, and whenever `ctl.flag_lnd=T`.

This simulation corresponds to the standard CMIP experiment **esm-piControl**.

The default carbon cycle configuration in the model is the so-called *closed* setup, in which marine sediments and weathering fluxes from continental weathering are disabled. This setup is expected to work reasonably well for time scales up to a few millennial. For longer-term simulations the `open` carbon cycle setup described below should be used.

The pre-industrial equilibrium state of a given model release, including biogeochemistry, is generally included in the `restart/pi` directory, which is also the default directory from where restart files are read (`ctl.restart_in_dir`).

## With carbon cycle (*open*)

The spinup procedure for the *open* carbon cycle spinup is a bit more complicated and requires a simulation of at least 100,000 years in order to reach an approximate equilibrium of the carbon cycle system. Chemical weathering on land and marine sediments need to be enabled with `ctl.l_weathering=T bgc.l_sediments=T` and various carbon cycle spinup flags need to be set (`ctl.l_spinup_cc=T ctl.nyears_spinup_bgc=5000 ctl.year_start_offline=1000000 bgc.l_spinup_bgc=T bgc.l_spinup_sed=T`):

```bash
runme -rs -q long -w 200:00:00 --omp 32 -o output/pi_spinup_cc_open \
-p ctl.nyears=100000 ctl.flag_bgc=T ctl.bgc_restart=F ctl.i_write_restart=1 \
ctl.l_spinup_cc=T ctl.nyears_spinup_bgc=5000 ctl.year_start_offline=1000000 bgc.l_spinup_bgc=T bgc.l_spinup_sed=T \
ctl.l_weathering=T bgc.l_sediments=T
```

The pre-industrial equilibrium state of a given model release, including biogeochemistry in an *open* carbon cycle setup, is generally included in the `restart/pi_cc_open`. 

## With ice sheets

```bash
runme -rs -q medium --omp 32 -o output/pi_ice_nh -p ctl.nyears=100000 ctl.n_accel=10 ctl.ice_domain_name=NH-16KM ctl.ice_model_name=yelmo \
ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T
```


