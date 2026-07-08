# Changelog

All notable changes to CLIMBER-X are documented here. 
Entries marked **Results** may change model output relative to the previous version.

## [1.5] - 2026-07-08

Major release: new solid-earth backend, restructured build/dependencies, and a
substantial overhaul of atmosphere, ocean, sediment, and conservation physics.
Expect model output to differ from v1.4.x.

### Added
- **FastEarth3D solid-earth backend** (`i_geo=3`) as a swappable alternative to VILMA, with its own sea-level diagnostic; VILMA viscosity-uncertainty perturbation (`f_visc_sd`).
- **Sediment hypsometry and depth classes** (sub-grid sediment hypsometry now the default); treatment of sediments during shelf drying/re-flooding and upward sediment shift.
- **Passive water tracers in the land model** (first step toward O18 isotope tracing); energy- and water-conservation check diagnostics (`water_check`) for atmosphere, ocean, and land.
- Ocean time-step ramp-up option (`l_ocn_dt_ramp`) for initial stability; Pa/Th tracers and sites.
- Freshwater-hosing domains for standard experiments and optional volume compensation of hosing.
- Land parameter for critical soil temperature suppressing respiration (enzyme denaturation); option to transfer slow→fast carbon on seasonal permafrost thaw.
- `runme` `--config`/`--list` options.

### Changed / Results
- **Atmosphere:** removed `sam2`, simplified precipitation/cloud parameterisations; uniformized lapse rates across surface types with increased moisture diffusivity; reduced ocean moisture exchange coefficient (lowers global precipitation); common minimum Coriolis parameter (1e-5) and retuned tropical sea-level pressure.
- **Ocean stability:** implicit diffusion, implicit EKE solution, and a Shapiro 2DX filter to suppress checkerboard patterns; **fast time step doubled**; ocean and sea-ice OpenMP/code optimisation; islands ordered by decreasing area; refined smoothing/extrapolation of initial T/S fields.
- **Biogeochemistry:** new simplified brine parameterisation (+coastal penetration-depth parameter); `q10=1.5` for marine NPP now default; weathering-flux tuning to balance ocean alkalinity; retuned CH4/N2O land emissions, LUC soil-carbon emissions, and permafrost thaw priming.
- **Reference density:** `rho0` updated from 1000 to 1025 kg/m³, with `rho_w` used for freshwater where appropriate.
- **SMB:** replaced TG24 simple SMB (`i_smb==3`) with a synthetic-elevation scheme.
- More dynamic topography (removed `k1_pot`); limited CO2 radiative-forcing increase at very high CO2.

### Build / Infrastructure
- **Migrated grid/mapping to the `fesm-utils` coords API and dropped the standalone `coordinates` dependency**; `coords` is now the default remapping backend (CDO not needed anymore).
- Added a `configme` manifest to track sub-packages; updated for **Yelmo v2.2** (`yhyd`/FastHydrology, nested under `yelmo/`).
- Migrated documentation to a Quarto docs site; installation reworked around `configme install climber-x`.

### Fixed
- Segfault in `ocn_to_cmn` from uninitialised `k_bmb/nk` on non-start-of-year ocean steps; ocean tracer conservation when ocean volume changes.
- Water-conservation fixes over lakes, in the sea-ice model, and in the atmosphere (`wcon` now prognostic instead of `qam`); stomatal resistance causing transpiration at `wilt==0`.
- Geothermal heat flux at the ocean bottom; `ice_model` 2D output (`H_w` → `W_til`); netCDF output file dimensions; duplicate `mo_m4ago_*` symbols at link (macOS).

## [1.4.3] - 2026-04-22

### Added
- Alternative fire disturbance scheme based on cumulative water deficit; direct CO2 emissions from fires (Da Nian).
- Time-dependent `scale_dhdt_ice` and 2D ocean temperature bias correction for basal mass balance.
- Option to scale marine NPP and POC sinking speed with global temperature anomaly (Christine); smooth CO2 fertilisation limitation.
- CFCs added to atmospheric output; freshwater-hosing contribution to relative sea level (`rsl_mass`) and refined steric diagnostics.
- Resolution-independent sub-grid bedrock std-dev handling (`l_use_z_bed_std_lowres`); island minimum-area namelist parameter.
- Separate `openmp_yelmo` build flag to control OpenMP in Yelmo independently.

### Changed / Results
- Restructured hosing module for multiple hosing domains; harmonized ice-mask definition (thickness-based) across components.
- CH4 radiative forcing corrected for very high concentrations (>5000 ppb).
- Tracer conservation reworked when ocean volume changes after a geo update.

### Fixed
- `q_geo` used instead of `q_geo_ice` during ice-sheet initialisation; sub-daily start/end-of-day detection; fractional hosing mask.
- Reversed i/j loop in land restart read (very long Intel oneAPI 2025 compile times).

## [1.4.2] - 2025-11-09

### Added
- Pangea ice grids; `bmb_grnd` to Yelmo output; shelf-area diagnostic and option to scale BGC tracers by shelf area.

### Changed / Results
- No temperature bias correction in SEMIX by default; `grounded_melt = F` by default in Yelmo (had been erroneously true).
- Reduce ocean time step instead of Coriolis-term reduction on CFL instability; broad code cleanup and neighbor-index optimisation.

### Fixed
- Sea-level change altering ocean layer count in top 1000 m, causing out-of-bounds access and segfaults in bmb fields.
- Uninitialised `ssh` in ocean model; orbital parameters out of bounds in insolation; time-dependent `fake_atm` forcing; multiple Yelmo domains; FFTW issues.

## [1.4.1] - 2025-09-16

### Added
- Interactive atmospheric N2O with land and ocean fluxes (constant lifetime); preformed tracers in ocean biogeochemistry.
- Separated equilibrium from reference bedrock topography with iterative-spinup options; separate geothermal heat flux for ice sheets vs. ocean (`dz_bed_ini`).
- Temperature-dependent climate sensitivity option; git version written to output.

### Changed / Results
- Elevation correction of precipitation uses `gamma` instead of a fixed 6.5 K/km lapse rate.
- Refined sea-level diagnostics (separated steric and mass contributions).

### Fixed
- `i_geo=0`: bedrock topography was not updated for sea-level changes.

## [1.4.0] - 2025-06-19

### Added
- Options for separate tropical/boreal deforestation and separate tree/grass interception efficiency.
- Variable atmospheric CH4 lifetime (Gang Liu); potential-energy-from-convection and deep-ocean salinity diagnostics.
- GRL-PAL-16KM / GRL-PAL-32KM grid definitions.

### Changed / Results
- **Parameter retuning**, including sea-level pressure at the equator and over topography (improved ITCZ/monsoon).
- Modified and retuned soil evaporation (only top ~20 cm accessible); removed fuel factor from fire disturbance.
- Timer reworked for sub-daily component timesteps (modules called on first timestep of the day). *Note: shorter BGC timestep affects climate — cause under investigation.*

### Fixed
- Corrected sign of `fw_dhdt_ice` in coupler; Yelmo restart fixes.

## [1.3.1] - 2025-02-07

### Added
- Yelmo set as the default ice-sheet model; option to treat ice-melt freshwater flux into the ocean separately.
- Multiple brine-rejection parameterisations (disabled by default); volume-averaged ideal age and additional ocean diagnostics.

### Changed / Results
- Migrated to `fesm-utils` static library (replaces `exlib`); CDO SCRIP now always used for mapping (removed `i_map==1`).
- Renamed `imo` → `bmb` (basal mass balance) and `mb_applied` → `mb_net`; smooth phenology/autotrophic-respiration transition made default.
- Increased default `tau_dmb` to 30 years (Yelmo and SICOPOLIS).

### Fixed
- `fovn` diagnostic bug; missing trailing `/` in SICO namelist; `ico2_emis==3` bug; Yelmo `H_calv_ref_thin` → `Hc_ref_thin`.

[1.5.0]: https://github.com/cxesmc/climber-x/compare/v1.4.3...v1.5
[1.4.3]: https://github.com/cxesmc/climber-x/compare/v1.4.2...v1.4.3
[1.4.2]: https://github.com/cxesmc/climber-x/compare/v1.4.1...v1.4.2
[1.4.1]: https://github.com/cxesmc/climber-x/compare/v1.4.0...v1.4.1
[1.4.0]: https://github.com/cxesmc/climber-x/compare/v1.3.1...v1.4.0
[1.3.1]: https://github.com/cxesmc/climber-x/compare/v1.3.0...v1.3.1
