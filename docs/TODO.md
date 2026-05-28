# TO DO #

- Update docs with latest compilation instructions.
- Delete compiled Yelmo version.
- Add `real(...,sp)` wrapper to writing yelmo output in `ice_model::yelmo_write_step_2D()`. 
- Get climber-clim version compiling with `gfortran`. 
- Delete compiled `coordinates` version. Add instructions to download `coordinates` and configure it.

- Consider removing dummy geo and bgc, and rather use preprocessor statements. At least for VILMA, this probably makes a lot of sense. For bgc, perhaps see below to make dummy usage easier.
- Where possible, make _def.f90 files for model components, so that the derived types are defined separately from the model itself. Then remaining dummy files will be much smaller and easier to maintain. 

## Compilation

- Make different compilation commands for different sections of the code (could help with compiling during development phases, to avoid repeatedly compiling full code)
- Make default model version the minimal version without major dependencies.

- Find out how to check which Makefile 'target' was called, to be able to properly set the LDFLAGS (LFLAGS+FFLAGS)

```
climber-clim: minimal configuration with ocn,atm,lnd,sic
climber-clim-bgc: plus with bgc
climber-clim-ice
climber-clim-bgc-ice
climber: climber-clim-bgc-ice

atm
bgc
bnd
ch4
co2
geo
ice
ice_sico
imo
lnd
lndvc
main
ocn
sic
smb
utils
vilma
yelmo: yelmo-static

possible later options:
ice-yelmo
ice-sico
```
