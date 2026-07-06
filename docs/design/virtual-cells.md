# Virtual-cell framework for land and surface mass balance

**Status:** design draft (not published — lives under `docs/design/`, excluded from the Quarto site)
**Scope:** a subgrid "virtual cell" framework that lets the land (`lnd`) and surface-mass-balance (`smb`) models be computed at many sub-grid samples (elevation bands × surface class) per coarse cell, then aggregated for coupling with CLIMBER-X.

---

## 1. Motivation

Processes strongly controlled by local topography — snow, melt, surface mass balance, near-surface temperature, vegetation limits — are poorly represented when a coarse coupler cell (~5°×5°) carries a single mean elevation. Today two separate pieces of code already compute overlapping physics at different resolutions:

- the **land model's ice tile** (`ebal_ice`, `ice_temp`, snow reservoirs) — snow + ice-surface energy balance over the ice fraction of a coarse cell, at cell-mean elevation;
- the **SMB model** (SEMI / PDD / simple) — the *same kind* of physics, on a dedicated high-resolution ice-sheet grid, downscaling coarse forcing to each fine cell's true elevation.

The goal is a single framework in which a coarse cell is decomposed into **virtual cells** — sub-grid samples each carrying their own elevation and surface class — the shared physics is run per virtual cell, and the results are aggregated back to the coarse cell for coupling. This:

1. gives the coarse land model a genuine topographic sub-grid (elevation-resolved snow, melt, temperature, vegetation);
2. unifies the land ice-tile and the SMB model onto one physics core, removing duplication;
3. exposes a tunable accuracy/cost knob (`n_vc`) instead of a fixed fine mesh — the speed-up opportunity.

**Design intent for the port:** `lnd` and `smb` source is *not* rewritten. Their physics routines are ported into the framework as single-column modules and called by the new driver. The reference model must be reproducible as the `n_vc = 1` special case.

---

## 2. Current architecture (summary)

- **Land** (`src/lnd/`) runs on the coarse coupler grid `ni×nj`. State is `lnd_class%l2d(i,j)` (one `lnd_2d_class` per cell) + `id_map`/`ij_1d` compression. `lnd_update_wrapper` OMP-parallelizes `lnd_update` over active land cells; each call handles one cell. Sub-grid tiling already exists — by **surface type** `nsurf=8` (5 PFTs + bare + lake + ice), plus vertical layers (`nl`, `nl_l`, `nlc`) and snow reservoirs (`nsoil=3`). Physics submodules loop internally over `nsurf` and aggregate by `frac_surf(n)` weights. **All tiles share one cell-mean elevation** — the limitation being removed.
- **SMB** (`src/smb/`) runs on separate high-res ice-sheet grids, called ~every 10 yr. `semi(i,j,…)` is invoked per cell in an embarrassingly-parallel loop; it takes **scalar** inputs (including pre-computed slopes `dz_dx_sur`, `dz_dy_sur`) and does its own forcing downscaling (`downscaling.f90`: T/q/LW/SW/precip lapse-rate + orographic factors). All horizontal work — topo filtering, slope maps, `map_field` regridding, aggregation to the coupler grid — is pre/post-processing. **SEMI is already single-column.**
- **Snow is always represented in the land model**, independent of `flag_smb`, via `nsoil` reservoirs (`is_veg`/`is_ice`/`is_lake`) with prognostic `w_snow`/`h_snow`. With `flag_smb=F` the SMB model is a no-op stub; what is lost is only the elevation/spatially-resolved SMB field needed to drive ice-sheet dynamics.
- **`src/lndvc/`** is a dormant skeleton (`vc_class{surf,zsrf,f_land,f_ice}`, `veg_vc(:,:,n_vc)`, dispatch on `surf`, empty `aggregate_cell`). It is compiled by a `make lndvc` target but never referenced in `climber.f90`.

---

## 3. Design principles / invariants

1. **The virtual cell is the fundamental unit of computation.** Everything is a loop over virtual cells.
2. **Physics is single-column and grid-agnostic.** A process module takes one virtual cell's state block + its downscaled forcing and returns tendencies/fluxes. It knows nothing about how many virtual cells exist, how they were enumerated, or how they aggregate. **All horizontal coupling lives in the geometry/decomposition service, never in the column physics.** (SEMI already satisfies this; the land submodules must be refactored to satisfy it by lifting the `nsurf` loop out into the driver.)
3. **Decomposition is a composable pipeline** of refine stages ending in leaf columns. Horizontal refinement (regrid) and vertical refinement (hypsometric split) are two stages of the same pipeline.
4. **Coupling with CLIMBER stays at the coarse cell.** The aggregated cell state presents the same interface the coupler consumes today, bounding the blast radius on `coupler.f90`.
5. **Conservation is enforced at the aggregation boundary** (energy, water, carbon) and is unit-testable in isolation.
6. **Reproduce-first.** `n_vc = 1` at cell-mean elevation must reproduce the reference land model. Every phase keeps a reproducing baseline.

---

## 4. Core abstractions

### 4.1 Virtual cell (leaf column)

A virtual cell is a *patch sample* of the sub-grid:

```
vc = { z        : elevation (m)
       dz       : band width (m)
       w        : area weight relative to the coarse cell (Σ w = 1)
       slope    : (dz_dx, dz_dy, |grad z|) for orographic downscaling
       class    : {land, lake, ice, shelf} }
```

Elevation and surface class are *different* sub-grid coordinates: elevation is continuous (band + area weight from hypsometry); class is categorical. A virtual cell is **one class at one elevation band** — a mosaic tile, not a fraction-of-everything tile. This lets an ice vc at 2000 m and a veg vc at 500 m coexist in one coarse cell.

**PFTs stay inside a land vc.** A land vc holds `pft_frac(npft)+bare`; competition (`dyn_veg`) is a local within-patch process, so PFTs are a sub-classification of a land vc, not their own virtual cells. This matches the skeleton's `class ∈ {land,lake,ice}` dispatch.

### 4.2 Single-column physics + per-class state blocks

State is decomposed into reusable blocks; a vc instantiates only the blocks its class needs:

```
land vc  → veg(pft) + soil_column + snowpack + soil_carbon
ice  vc  → snowpack + ice/firn_column + subglacial_carbon      (SMB is a diagnostic here)
lake vc  → lake_column + snowpack
shelf vc → shelf_column
+ every vc emits one shared surface_flux block (SH, LH, albedo, T_skin, runoff, evap)
```

The **snowpack is the unification linchpin**: one snow state block and one snow module, shared by every class that can hold snow. Then **SMB is not a separate model** — over an ice vc it is the mass budget of the shared snow+ice column:

```
smb = snowfall + refreezing − melt − sublimation − runoff
```

The shared `surface_flux` block is the coupling currency that aggregation reduces.

**Class blocks are allocated per class** (Fortran allocatable derived-type components), so a vc carries only what it needs. This makes reduced configurations fall out for free: an **SMB-only / no-vegetation** run instantiates only ice (and snow-bearing) vcs — the `veg`/`soil_carbon` blocks are never allocated and the veg path never dispatched; a **land-only** run omits the ice blocks.

**Permafrost.** Permafrost is not a separate model — it is emergent from the soil thermal column (`t_soil`, `theta_w/theta_i` freezing/phase) plus soil-carbon thaw dynamics (`frozen_years`, `thaw_timer`, `k_slow_to_fast`, active-layer thickness `alt`). It lives entirely in a land vc's `soil_col` + `soil_carbon` blocks (thermal and permafrost state kept together; `alt` is diagnosed from the profile). The key structural point: **depth is never a virtual dimension — only the surface (elevation × class) is.** Each vc owns a full-depth soil column with the standard shared depth grid; what varies per vc is the *surface boundary condition* (elevation-downscaled temperature and snow cover). So active-layer thickness, thaw, and permafrost carbon become **elevation-resolved within a coarse cell**. The deep column state is necessarily per-vc (thermal history integrates the per-vc surface BC); resolving only the active layer per-vc and sharing deep thermal is a possible future optimization, not the base design. On a land→ice class change the soil column becomes subglacial and is handed over by the conservative `remap_state`.

### 4.3 Decomposition pipeline

```
coarse coupler cell (5×5, atmosphere + coupling)
   │  ① regrid forcing (bilinear) ▼            ← optional horizontal refine stage
   mid-res cell (2×2 / 1×1, own mean elevation from high-res reference)
      │  ② hypsometric split ▼                 ← vertical refine stage
      virtual cell (elevation band × class)
         │  ③ downscale forcing to band elevation ▼
         single-column physics  →  surface_flux + state
   ▲ aggregate all leaf columns (area-weighted, conserving) back to 5×5 for coupling
```

The pipeline is configuration, not special-case code. Backends are pipeline configs over one physics core:

| backend | pipeline | purpose |
|---|---|---|
| virtual-cell | [hypsometric split] | coarse land + climate coupling, cheap |
| mid-res + vc | [bilinear regrid] → [hypsometric split] | horizontal + vertical refinement |
| high-res ice | [regrid to ice grid] (→ optional split) | ice-sheet dynamics, high fidelity |

---

## 5. Shared services (the cell ↔ vc boundary)

- **Downscaling** — promote `smb/downscaling.f90` to shared infrastructure. Input: coarse (or mid-res) forcing + the reference elevation it is valid at + a vc's `{z, slope}`. Output: vc-level T/q/LW/SW/precip/wind. Used by *every* class, so land vcs finally get proper per-elevation forcing. **Interpolate forcing and its reference elevation together** (SMB's `z_sur_i` vs `z_sur` pattern) so the downscaling delta `z_vc − z_ref` is correct.
- **Snowpack + surface energy balance** — one shared module set, used by all snow-bearing classes; the mechanism that collapses the land ice-tile and SMB duplication.
- **Aggregation / conservation** — vc ensemble + area weights + per-variable reduction rule (area-weighted mean / sum / flux-conserving) → coarse-cell state for coupling. Conservation (energy/water/carbon) enforced and unit-tested here. **Conservative fields (precip above all) must be renormalized after any non-conservative regrid so the coarse-cell mean is preserved**, and the leaf-column reduction must be strictly area-conservative — otherwise the coarse budget CLIMBER sees drifts.
- **Decomposition remap** — when ice advances/retreats or topography shifts, band membership/class changes → a named, tested, conservative state-transfer operation. Replaces the ad-hoc "initialize newly vegetated/ice/lake cell" branches scattered through `lnd_update`.

---

## 6. Keeping (and unifying with) the high-res SMB

The high-res ice grid is **not a compromise** against this design — it is the design's high-fidelity mode. Elevation-classes and the fine ice mesh are two ways to enumerate columns over one physics core; the physics (snow, energy balance, SMB budget) is identical, so there is no duplicated *physics*, only a second **decomposition/mapping frontend**.

Elevation-class SMB alone may be insufficient for ice sheets because hypsometric binning discards horizontal structure that matters most at margins/ablation zones (aspect, slope, distance-to-coast, wind exposure, the sharply nonlinear SMB–ELA feedback, and the 1:1 map to ice-grid cells the elevation feedback needs). **All of these are geometry/forcing concerns**, so they belong to the high-res backend's decomposition + downscaling layer — a mapping/regridding module you would write for the virtual cells anyway — not a second snow/energy-balance model.

Payoff: a **tiered/adaptive SMB** from one codebase — cheap elevation-class SMB globally, full high-res SMB over ice domains — **consistent by construction**: where both exist, the coarse aggregate is the conservative reduction of the high-res field. Today's inconsistency (land ice-tile formulation ≠ SMB-model formulation) collapses to a mere resolution difference of one formulation.

**Invariant to protect:** the moment a genuinely neighbor-coupled process appears (e.g. lateral snowdrift/wind redistribution), it must go in the decomposition/geometry layer, never inside the column physics, or the two backends diverge.

---

## 7. Refinement knobs and the horizontal/vertical trade-off

Two complementary refinement axes, both configurable, both with an identity baseline:

- **Horizontal — mid-res zoom.** Bilinearly interpolate coarse forcing to a mid-res grid (2×2 / 1×1) before the hypsometric split, using the existing `coords` `map_field`/`map_init` machinery (`method="bil"`). Identity baseline: mid-res = coarse.
- **Vertical — `n_vc`.** Number of elevation bands per (mid-res) cell. Identity baseline: `n_vc = 1` at cell-mean elevation ≡ reference land model.

Cost ≈ `n_midres_per_coarse × n_vc_per_midres × ncells`. The axes substitute: finer horizontal resolution shrinks the elevation spread within each cell, so fewer bands are needed. Reallocate a refinement budget by field character:

- smooth synoptic gradients (T with latitude, precip with distance-to-coast) → the **horizontal/bilinear** axis is efficient;
- terrain-driven structure (SMB near the ELA, snow line) → the **elevation-band** axis is efficient.

**Honesty caveat:** bilinear interpolation adds no new sub-grid information — it resolves the coarse field's existing gradient smoothly and removes the blocky cell-boundary artifact, and pairs forcing better with high-res elevation before downscaling. The terrain-driven variance still comes entirely from the hypsometric split + downscaling. The two are complementary; the zoom sharpens the baseline the elevation classes then perturb.

Parallelism is over the flat `(cell, midres, vc)` active leaf list (generalizing today's `ij_1d`) — embarrassingly parallel.

---

## 8. Data layout

```
framework container
  ├─ decomposition: per coarse cell → list of leaf vc descriptors {z,dz,w,slope,class}
  │                 (compressed active-leaf list for OMP, like ij_1d over (cell,vc))
  ├─ per-vc state: class-specific blocks (veg/soil/snow/ice/lake/carbon) + surface_flux
  ├─ aggregated cell state: coarse-cell means presented to the coupler (same interface as today)
  └─ global 0-D state (carbon totals, isotope ratios, …), as lnd_0d_class today
```

This supersedes the skeleton's single 150-field `lndvc_veg_class` used for both `veg_vc` and `veg`: state is split into per-class blocks so a vc allocates only what it needs, and the aggregated cell state is a distinct (thin) structure.

---

## 9. Relationship to the existing `src/lndvc/` skeleton

The skeleton already reaches for this: `vc_class{surf,zsrf,f_land,f_ice}`, `veg_vc(:,:,n_vc)`, `select case(surf)` dispatch, an `aggregate_cell` reducer, and a `lndvc_grid` copy of the land levels. What it is missing and this design adds:

- shared snow / energy-balance / downscaling services (skeleton has an empty `smb_vc` copy instead of SMB-as-diagnostic);
- per-class decomposed state blocks instead of one mega-type;
- the decomposition pipeline (hypsometry + mid-res regrid) and conservation-checked aggregation;
- the scalar-per-field choice is correct **provided** PFTs live inside a land vc.

The framework is built out in `src/lndvc/`, continuing the skeleton.

---

## 10. Porting strategy and phased plan

**`lnd` and `smb` are not rewritten.** Their physics routines are ported into the framework as single-column modules and called by the new driver. Development happens on a feature branch (and a worktree if experiments are running). Each phase keeps a reproducing baseline and is independently testable and mergeable.

- **Phase 1 — scaffold + identity check.** Build the container, decomposition (`set_vc` from hypsometry), dispatch (`update_vc`), and conservation-checked `aggregate_cell`; wire into `climber.f90` behind `flag_lndvc`. Run with `n_vc = 1` at cell-mean elevation → reproduce the reference land model. Validates the whole scaffold with zero physics change.
- **Phase 2 — per-vc forcing + single-column physics.** Port the land submodules to single-column, single-class signatures (lift the `nsurf` loop into the driver); add per-vc forcing downscaling (shared `downscaling`). Still `n_vc = 1`, still reproducing.
- **Phase 3 — multiple land bands.** `n_vc > 1`, real hypsometry from the high-res reference, downscale per band, two-stage aggregation. Validate topographic signal + strict conservation.
- **Phase 4 — mid-res zoom.** Add the optional bilinear horizontal refine stage with precip renormalization + conservative aggregation. Tune the horizontal/vertical split.
- **Phase 5 — SMB as virtual cells.** Route ice-surface SMB/snow/melt through the ice-vc path; SMB becomes a mass-budget diagnostic. Keep the high-res ice backend as the high-fidelity decomposition; make the coarse aggregate the conservative reduction of the high-res field. Benchmark speed vs. the current high-res SMB.
- **Phase 6 — consolidation.** Retire the duplicated land ice-tile / SMB snow-energy code paths in favour of the shared services; tune; document.

---

## 11. Decisions

**Settled:**
- **Band scheme:** fixed elevation edges. A band always means the same elevation everywhere and forever, so a vc's prognostic column stays bound to its band — no continuous vertical remapping. Only the band *area weight* changes over time (cheap; used in aggregation), and discrete occupied↔empty / class-change events fire the conservative `remap_state`. Fixed edges also make cross-cell comparison and NetCDF output trivial.
- **Type layout:** per-class blocks, built fresh in `lndvc_def.f90`, functionality kept as close as possible to the reference `lnd`/`smb`. Permafrost thermal+state kept together in `soil_col`. Single shared `snowpack` block across land/ice/lake.
- **Naming:** keep `lndvc` for now; a neutral `vcell`/`vc` rename is a mechanical module-prefix pass later if warranted.
- **Reduced configs:** SMB-only (no vegetation) and land-only are supported via class-allocated blocks.

**Open:**
1. **Whether to fully unify the snowpack** in Phase 5/6 (retire the land ice-tile snow path) or leave both during a transition.
2. **Mid-res default** target resolution and per-region `(mid-res, n_vc)` budgets.
3. **Aggregation currency:** confirm the exact set of coarse-cell fields the coupler requires so `vc_cell_t` is a drop-in for `lnd_to_cmn` / `smb_to_cmn`.

---

## 12. Implementation status (branch `lndvc-framework`)

**Done:**
- **Phase 1 scaffold + integration** — `flag_lndvc` (control + `nml/control.nml`, default off) wired into `climber.f90` (init/update/end); `obj_lndvc` linked into all four executable variants; `make climber-clim` builds and links. Data model rebuilt as per-class blocks in `lndvc_def.f90` (`vc_desc`, `forcing`, `surface_flux`, `snowpack`, `soil_col` incl. permafrost, `veg`, `soil_carbon`, `ice_col` incl. SMB diagnostic, `lake_col`; thin `vc_cell_t`; `veg_glob`), isotopes included. Services: `lndvc_decomp` (identity-baseline `decompose`: one leaf per present class at cell-mean elevation, fed from the land model's own tiling), `lndvc_aggregate` (fraction roundtrip + sanity check), `lndvc_downscale`, `lndvc_model` driver (alloc/update/dispatch).
- **SMB single-column physics ported** (pattern A = bodies unchanged, B = params brought in, so `lndvc` is standalone; only shared utils `precision`/`tridiag`/`nml` reused). Files under `src/lndvc/`: `constants.f90` (`const_m`), `thermo.f90` (`thermo_m`), `smb/params.f90` (`smb_par_m`), `smb/grid.f90`, `smb/snow.f90`, `smb/downscaling.f90`, `smb/temp.f90`, `smb/ebal.f90`, `smb/surface_par.f90`, `smb/semi.f90`. The full SEMI chain (downscaling → albedo → resistance → ebal → temp → update_tskin → snow) compiles and links. Distinct object names (`lndvc_smb_*.o`) avoid the flat-`obj/` clash with `src/smb`.
- **(a) Ice/SMB path is live end-to-end** (compile + link only; run validation gated on the cluster). Four steps:
  1. **`forcing_t` expanded** with the reference-level (`_i`) forcing SEMI needs (`tam_i`, `t2m_bias_i`, `gam_i`, `tstd_i`, `ram_i`, `u700_i/v700_i/wind_i`, `prc_i`, `prc_bias_i`, SW components + `dswd_*` sensitivities, `swd_toa_i/min`, `cld_i`, `lwdown_i`, `gam_lw_i`, `dust_i`, `coszm_i`) + `z_sur_i`/`dTvar`/`f_ice`/`alb_ice`; vc geometry (`z_sur_std`, `dz_sur`, `f_ele`) on `vc_desc_t`. New carry-over homes: `refreezing_sum`/`dt_snowfree` (snowpack), `t_prof`/`t_prof_old`/`runoff` (ice_col). Ice-vc inner arrays sized + initialized once in `set_leaf`/`alloc_ice_blocks` (lean: only fields SEMI reads/writes as carry-over; diagnostics are call-local).
  2. **`cmn_to_lndvc`** (coupler) populates each ice vc's `_i` forcing from the coarse cell — analogue of `cmn_to_smb`'s coarse block, `z_sur_i` = coarse-cell mean; wired into `climber.f90` before `lndvc_update`. Identity baseline: no spatial filter, no bias correction, `dTvar=0`, `tstd_i=2`.
  3. **`semi` called in `lndvc_update_ice`** (pattern A unpack of `forc`/`snow`/`ice`/`flx`/`desc`); barometric `pressure` + `f_ele` computed framework-side from vc geometry, `f_ice=1`/`alb_ice=firn`. SMB mass budget written into `ice_col`: `smb = snow − evp − snowmelt − icemelt + refreezing`; `melt`, `runoff`.
  4. **`aggregate_cell`** reduces the ice-surface coupling currency (`t_skin`, albedo, `flx_sh/lh/g`, evap) as an area-weighted mean over live ice vcs and the mass budget (`smb`/`melt`/`runoff`) as a conserving cell-area-weighted sum (`smb`/`melt` added to `vc_cell_t`).

**Next step — (b) validate + widen:**
1. **Cluster run** the identity ice path (`flag_lndvc=T`) and sanity-check aggregated SMB against the reference `smb` model over an ice domain.
2. **Port the land + lake single-column paths** (`lndvc_update_land`/`_lake` are still no-ops), lifting the `nsurf` loop into the driver so the full coupling currency aggregates across all classes (not just ice).
3. **Phase 3 — multiple land bands**: real hypsometry + per-band downscaling (populate the currently-zeroed `z_sur_std`/`dz_sur`/slopes/`f_ele`).

**Notes / blockers:**
- Local runs cannot currently reach the time loop: a geo/map-generation **segfault** during init (reproduces with `flag_lndvc=F` and in `--gen-maps`), independent of `lndvc` — likely a local **memory limit** on macOS. Validate on the cluster (PI config: `run_bench.sh`, `runme -r ... -o output/... -p ctl.nyears=... ctl.flag_lndvc=T`).
- **Build switch `LNDVC` (default 0 = off).** The framework is gated behind a compile flag so it never forces into standard builds: `make climber-clim` builds *without* it (no `src/lndvc/` objects compiled or linked; the `use`/`call` sites in `climber.f90`/`coupler.f90` are `#ifdef LNDVC`-removed), while `make LNDVC=1 climber-clim` adds `-DLNDVC` and links `obj_lndvc`. Setting `flag_lndvc=T` in an executable built with `LNDVC=0` stops at init with a rebuild hint. Because compiler-flag changes don't bump timestamps, switching `LNDVC` mode needs the main objects removed first (`rm obj/coupler_nobgc.o obj/climber_clim.o`). The former `lndvc_def.mod` dependency gap is fixed: `coupler*.o` now depend on `$(dep_lndvc_def)` and the `climber*.o` rules on `$(dep_lndvc_main)` (both empty under `LNDVC=0`).
- The reference `lnd`/`smb` code is untouched; consolidation (retiring duplicated paths) is the later Phase 6.
