# Only ocean simulations

Instructions to perform simulations of the ocean with imposed boundary conditions.


## Notes on parameter choices

In nml/ocn_par.nml, this section relates to options for imposing fixed boundary conditions for the ocean.

```
 !-------------------------------------------------------------------------------------------
 ! fixed boundary conditions settings
 !-------------------------------------------------------------------------------------------
 l_ocn_fix_wind = F           ! used prescribed wind stress daily climatology?
 l_ocn_fix_fw   = F           ! used prescribed surface freshwater flux daily climatology?
 l_ocn_fix_flx  = F           ! used prescribed surface heat flux daily climatology?
 i_ocn_input_fix = 1          ! what fixed ocean input to use, only if any l_ocn_fix_*==T :
                              ! 1 = average of first 30 years of simulation
                              ! 2 = read from file ocn_input_fix_file
 l_ocn_input_fix_write = F    ! write average daily ocean input variables from last 30 years of simulation to file?
 ocn_input_fix_file = "input/ocn_input_fix_pi.nc"     ! file with daily climatology of ocean input fields to be used for simulations with any l_ocn_fix_*==T
 ```

So you would run a PI equilibrium simulation with ocn.l_ocn_input_fix_write=T, which writes an ocn_input_fix.nc file with the daily climatology of the last 30 years of all ocean input fields.

Here is the command from the pi-spinup case, with the additional argument added.

```bash
runme -rs -q short -w 24:00:00 --omp 32 -o output/pi_spinup \
-p ctl.nyears=5000 ctl.i_write_restart=1 ocn.l_ocn_input_fix_write=T
```