# Last Glacial Maximum (LGM)

Instructions to perform simulations of the LGM.

## Climate only

This simulation corresponds to the standard PMIP experiment **lgm**.

The length of the simulation is set to 10000 years, which is generally enough to reach an equilibrium of the climate system.

```bash
runme -rs -q short -w 24:00:00 --omp 32 -o output/lgm -p ctl.nyears=10000 \
ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 \
ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200 \
ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_geo_ref_file=input/geo_ice_tarasov_lgc_0ka.nc geo.geo_ref_file=input/RTopo-2.0.1_0.125deg_DRThydrocorr.nc \
ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc 
```

The simulation starts from a pre-industrial equilibrium state of the climate and switches instantaneoulsy to LGM boundary conditions for orbital parameters, GHGs and continental ice sheets. 

The ice sheet reconstruction used here is GLAC1D from Lev Tarasov. Prescribing ice sheets different from present-day requires to separately prescribe the bedrock elevation (`ctl.fake_geo_const_file`) and ice thickness (`ctl.fake_ice_const_file`). From that, land-sea mask, surface elevation and runoff routing directions are derived internally in the model.
For a realistic dynamic derivation of the runoff routing directions a very high resolution topography is needed, while paleo-reconstructions of topography are usually relatively coarse. That is why bedrock elevation is computed from **anomalies** relative to present day, set as the difference between `ctl.fake_geo_const_file` and `ctl.fake_geo_ref_file` (they must be on the same resolution!). These anomalies are then interpolated and added on top of the high-resolution reference present-day bedrock elevation defined by `geo.geo_ref_file`. The `ctl.fake_geo_ref_file` and `geo.geo_ref_file` in the command above are actually the default ones, and are included only for clarity.
The ice thickness is read from the file specified by `ctl.fake_ice_const_file` and will be simply bilinearly interpolated onto the high-resolution of `geo_ref_file`.

The change in topography and ice sheets also affect the ocean volume. In the process of switching to LGM conditions, the total salt content in the ocean is conserved, which leads to an increase in volume-averaged salinty by about 1 psu as a result of the ~100 m sea level drop.

## With interactive ice sheets

The following command runs an LGM simulation with interactive ice sheets in the NH, starting from the LGM ice sheet reconstructions of GLAC1D:

```bash
runme -rs -q medium --omp 32 -o output/lgm_ice_nh -p ctl.nyears=100000 ctl.n_accel=10 \
ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 \
ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200 \
ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_geo_ref_file=input/geo_ice_tarasov_lgc_0ka.nc geo.geo_ref_file=input/RTopo-2.0.1_0.125deg_DRThydrocorr.nc \
ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc \
ctl.ice_domain_name=NH-16KM ctl.ice_model_name=yelmo \
ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T \
ctl.restart_in_dir=restart/lgm
```

Before running an LGM simulation with interactive ice sheets it is recommended to run a climate-only snapshot to get the climate into equilibrium with the LGM boundary conditions.
Then the simulation with interactive ice sheets is initialized using the climate state from that simulation (`ctl.restart_in_dir=restart/lgm`, needs to be created first, i.e.: `cp -r output/lgm/restart_out/year_10000 restart/lgm`) and the same GLAC1D ice sheets.

Since ice sheets have a long equilibration time, the simulation is run for 100,000 years, but using an acceleration of the climate of a factor 10 (climate is updated only every 10 years).

