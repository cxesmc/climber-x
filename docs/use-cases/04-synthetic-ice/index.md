# Synthetic ice-sheet geometry

CLIMBER-X can build an ice-sheet geometry from a prescribed **target extent** instead of, or alongside, a dynamic ice sheet model. The synthetic geometry is a stateless function of the target mask, the current bedrock and sea level: a signed distance to the mask boundary sets a perfect-plastic (or linear) surface profile on top of the bedrock, and the grounded ice thickness follows from that. It lets the whole climate (atmosphere, land, ocean, SMB) see an ice sheet that follows a reconstruction through time, while a real ice sheet model can run in parallel and receive the resulting surface mass balance.

## Ice-sheet configurations

| Mode | `flag_ice` | `ice_model_name` | Climate topography | Solid Earth and sea level |
|---|---|---|---|---|
| Imposed thickness | `F` (+ `ifake_ice=1` for transient) | – | `fake_ice_*_file` | prescribed ice |
| Interactive | `T` | `yelmo`, `sico` or **`syn`** | ice model | ice model |
| Shadow | `T`, `l_ice_syn=T`, `i_ice_topo_clim=1` | `yelmo` or `sico` | **synthetic** | ice sheet model |

In the **shadow** mode the synthetic geometry is computed every `n_year_ice` on the same regional grid, bedrock and sea level as the ice sheet model. `ice_to_geo` passes the synthetic thickness to the geography module (`geo%hires%h_ice`, from which surface elevation, masks, lakes and runoff routing are derived), while the ice sheet model's thickness is passed separately as the solid-Earth load (`geo%hires%h_ice_load`) and sets the sea-level change. `ice_to_smb` passes the synthetic surface and ice mask to the SMB grid, so SEMI (or the simple SMB scheme) computes the mass balance on the synthetic ice sheet, and that SMB is applied to the ice sheet model. The ice sheet model therefore tends to grow toward the reconstructed geometry without its own biases feeding back into the climate.

With `l_ice_syn=T` and `i_ice_topo_clim=0` the synthetic geometry is only diagnosed (written to `ice_syn_<domain>.nc`) and the ice sheet model drives everything, which is useful for comparison.

## Target extent

The target extent is read from a lat/lon netCDF file set in `nml/ice_syn_par.nml` (group `&ice_syn`, runme alias `icesyn`):

```fortran
mask_file  = "input/Batchelor2019_LGM_icemask.nc"
mask_var   = "mask"    ! mask_var(lon,lat) or mask_var(lon,lat,time)
var_thresh = 0.5       ! cells with mask_var > var_thresh are target ice
```

A file with a `time` dimension (model years, same convention as `fake_ice_var_file`) is linearly interpolated between slices and clamped outside its range. An ice-thickness reconstruction can be used as the target by setting `var_thresh` to a thickness (e.g. `10.`). The mask is mapped conservatively to the ice grid and thresholded at a fraction of 0.5.

## Geometry parameters

```fortran
use_plastic = T        ! T => plastic Nye/Vialov profile; F => linear wedge
tau0        = 36.0e3   ! Pa, plastic basal yield stress
z_max_in    = 2500.0   ! m, cap on the profile height above the margin
h_grd_min   = 10.0     ! m, minimum thickness above flotation inside the mask
```

Inside the target mask the thickness is at least the flotation thickness plus `h_grd_min`, so all synthetic ice is grounded; outside the mask there is none. The profile kernels live in `src/ice/ice_syn_topo.f90` and are shared with the simple SMB scheme.

## Simple SMB scheme

The simple SMB scheme (`smb.i_smb=3`) builds its own synthetic elevation from a static target mask by default. When the synthetic geometry comes from the ice component, set

```fortran
&smb_simple
  l_z_syn_external = T
```

so the scheme uses the delivered surface elevation and ice mask directly (target mask refreshed every year).

## Examples

LGM with the synthetic ice sheet as the ice model (mode 2):

```bash
runme -rs -q short --omp 32 -o output/lgm_syn -p ctl.nyears=1000 \
ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 \
ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200 \
ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc \
ctl.ice_domain_name=NH-32KM ctl.ice_model_name=syn \
ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=F \
ctl.restart_in_dir=restart/lgm
```

Yelmo driven by the SMB computed on the synthetic ice sheet (mode 3):

```bash
runme -rs -q medium --omp 32 -o output/lgm_yelmo_syn -p ctl.nyears=20000 ctl.n_accel=10 \
ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 \
ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200 \
ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc \
ctl.ice_domain_name=NH-32KM ctl.ice_model_name=yelmo ctl.l_ice_syn=T ctl.i_ice_topo_clim=1 \
ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T \
ctl.restart_in_dir=restart/lgm
```

## Output

- `ice_<domain>.nc` (mode 2) or `ice_syn_<domain>.nc` (shadow): `f_ice_target`, `mask_target`, `dist`, `z_syn`, `H_ice`, `z_srf`, `z_bed`, `z_sl` and the received `smb` (diagnostic only for the synthetic model).
- `geo_hires.nc` and the geo restart: `h_ice` (climate) and `h_ice_load` (solid Earth), identical unless the shadow mode is active.
