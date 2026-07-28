---
title: "The CLIMBER-X Atmosphere (SESAM)"
subtitle: "Structural Analysis and a Roadmap Toward GCM-like Generality at Low Cost"
author: "Prepared for A. Robinson"
date: "28 July 2026"
toc: true
toc-depth: 2
numbersections: true
geometry: "margin=2.4cm"
fontsize: 11pt
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
mainfont: "Helvetica Neue"
monofont: "Menlo"
header-includes:
  - \usepackage{booktabs}
  - \usepackage{array}
  - \renewcommand{\arraystretch}{1.25}
  - \usepackage{sectsty}
  - \definecolor{accent}{RGB}{20,70,130}
  - \allsectionsfont{\color{accent}}
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[L]{\small SESAM Analysis}
  - \fancyhead[R]{\small \thepage}
  - \fancyfoot[C]{}
---

# Executive summary

CLIMBER-X's atmosphere ("SESAM") is a *statistical–dynamical* model on a 5°×5° grid — a higher-resolution refinement of the CLIMBER-2 atmosphere. Its speed comes from a radical reduction of the prognostic state: the **entire time-integrated atmospheric state is two-dimensional** — five fields stepped on a 3-hour inner loop — while everything with vertical structure, momentum, clouds, or eddy transport is **diagnosed each step by imposing a prescribed shape** onto those 2-D fields. Given boundary conditions (SST, sea ice, land state) and radiative forcing, the atmosphere has a unique, reproducible, weather-free solution.

This design is the source of both its efficiency and its three principal generality gaps relative to a GCM:

1. **Vertical realism** — free-tropospheric temperature and humidity are algebraic functions of a single surface anchor and diagnosed lapse rates; upper-air structure cannot decouple from the surface, and climate sensitivity (ECS) is a tuning knob rather than an emergent property.
2. **Circulation and longitudinal structure** — the mean meridional circulation is an imposed three-cell template, and stationary waves are computed by a thin, barotropic, orography-only, latitude-by-latitude solve. Regional structure is largely absent.
3. **Internal variability** — there is no entropy source anywhere in the atmosphere. It cannot produce NAO/annular-mode, blocking, storm-track, or ENSO-seeding variability; tellingly, the coupled model retrofits stochastic forcing *downstream* in the ocean because the physically correct source — atmospheric eddies and winds — is silent.

The good news is that the codebase is unusually well-suited to targeted upgrades. Layer-resolved transport already exists, an unconditionally-stable implicit-diffusion pattern is established, the synoptic-energy solver is a ready hook for stochastic and richer-eddy physics, and the model is already instrumented with radiative kernels and a partial-radiative-perturbation feedback analysis to *verify* that any change produces GCM-like feedbacks. A staged roadmap can move SESAM substantially toward GCM-like generality — including genuine internal variability — while preserving its defining speed.

---

# How SESAM works today

## The 2.5-D design: what is prognostic

The grid is 72×36 (5°×5°) with `km=13` prescribed pressure levels. The slow atmospheric step is one day; a fast inner loop runs `nstep_fast = 8` sub-steps of 3 h (`tstep = 10800 s`). The **complete prognostic state is five 2-D fields**:

| Prognostic field | Meaning | Updated in |
|---|---|---|
| `tam(im,jm)` | near-surface / column-anchor air temperature | `time_step.f90:301` |
| `wcon(im,jm)` | vertically-integrated column water | `time_step.f90:192` |
| `dam(im,jm)` | column dust | `time_step.f90:272` |
| `cam(im,jm)` | column CO~2~ | `time_step.f90:288` |
| `sam(im,jm)` | synoptic (eddy) kinetic energy | `synop.f90:225` |

Everything with a vertical dimension — `t3, q3, tp, d3` and the 3-D winds — is **reconstructed every sub-step** by imposing prescribed profile shapes. The 13 levels are quadrature points for these profiles, **not prognostic degrees of freedom**. There is no prognostic momentum, vorticity, divergence, or vertical velocity.

## Vertical structure (`vesta.f90`)

Temperature is a piecewise profile with **three diagnosed lapse rates** — surface-layer `gams`, lower-troposphere `gamb`, upper-troposphere `gamt` — up to a prognostic tropopause `htrop`, isothermal above. Humidity is constant `ram` through the boundary layer, then an exponential decay with a single scale height `hrm`, then a stratospheric value. The column-water integral is split into a tropospheric slope and stratospheric intercept so that `wcon = ram·A_trop + W_strat`; the prognostic `wcon` is inverted to recover `ram`, keeping the water budget conservative while the profile itself stays diagnostic. The structure is thus **shape-prescribed, amplitude-diagnosed**.

## Winds — fully balanced and diagnostic

There is no momentum equation. The 3-D wind (`u3d.f90`) is a sum of three balanced components:

- **Barotropic geostrophic** from the sea-level-pressure gradient, with Coriolis floored near the equator;
- **Thermal wind** integrated from horizontal gradients of the reconstructed temperature — the source of the jets — explicitly damped at equator and poles to mask the geostrophic breakdown;
- **Ageostrophic Ekman flow** in the boundary layer, given a prescribed vertical shape plus a compensating return flow so its vertical-integral mass transport vanishes.

Mass continuity is enforced *diagnostically*: each column's net convergence is compensated in the upper troposphere, and the **vertical velocity is the residual of column mass convergence** — never a prognostic equation.

## Mean meridional circulation and stationary waves (`slp.f90`)

The Hadley/Ferrel/Polar structure is an **analytic three-cell template** (`zslp`); its only responsive freedoms are amplitude (from zonal-mean temperature gradients), ITCZ latitude, and Hadley width. It cannot produce a double ITCZ, a different cell count in extreme climates, or the correct nonlinear Hadley widening.

Stationary waves (`azslp`) are a genuine physics-based but thin calculation: a **per-latitude, barotropic, topography-only** linearized stationary Rossby-wave solve by FFT in longitude, plus a thermal term from sea-level-temperature anomalies. Each latitude is solved independently (no meridional coupling), and there is no stationary-wave response to zonally-asymmetric diabatic heating.

## Eddies (`synop.f90`)

Synoptic activity is a **single prognostic scalar** `sam` (eddy kinetic energy): production from the Eady maximum baroclinic growth rate (Hoskins–Valdes style), dissipation `∝ sam^{3/2}`, advection by the steering wind, all advanced by an unconditionally-stable ADI solver. From `sqrt(sam)` the scheme diagnoses **purely down-gradient** heat/moisture diffusivities, the synoptic surface wind and stress, and a synoptic vertical velocity. There is **no random term** — the closure is a smooth relaxation toward the baroclinicity-implied equilibrium, so it generates transport but no variability. (An ERA-Interim monthly EKE climatology is read at init and carried for diagnostics but is not used dynamically — a useful calibration target, currently a dead input.)

## Column physics: clouds, moisture, radiation

- **Clouds** (`clouds.f90`) are fully diagnostic: a relative-humidity/vertical-velocity fit for large-scale cloud, a near-surface RH-gradient fit for inversion cloud, combined by random overlap and heavily time-relaxed (0.1 new / 0.9 old). There is **no prognostic cloud water or ice, no microphysics, a single cloud deck**; fraction, height, and optical thickness are independent empirical functions.
- **Moisture and precipitation** (`time_step.f90`) use a **bulk statistical closure** — no CAPE, no convective plume, no separate large-scale condensation: precipitation is moisture convergence + evaporation + a residence-time term, scaled by relative humidity, with tuned ocean/land residence times. Rain/snow split by a linear 2 m-temperature ramp.
- **Longwave radiation** (`lwr.f90`) is a broadband grey/emissivity **transmission-function** scheme (not a band model, not correlated-k). Non-CO~2~ greenhouse gases are folded into an *effective* CO~2~. The CO~2~ transmission needs an ad-hoc corrective factor to behave at high concentrations, and **ECS is set by an explicit tuning multiplier** (`ecs_scale`, with a state-dependent adjustment) rather than derived.
- **Shortwave radiation** (`swr.f90`) is a two-band integral-transmission-function column model. A genuine strength: it analytically provides surface-albedo and elevation derivatives of downward SW, so the column scheme is **differentiable by construction** — used by the surface downscaling coupler.
- **Surface fluxes** are computed by the surface components, not the atmosphere; SESAM supplies a **neutral (stability-independent) drag** and a fixed-depth boundary layer, with an Ekman cross-isobar turning angle for the surface wind.

## Temporal structure and coupling

The atmosphere integrates on two clocks: a **daily** "slow" step for radiation, pressure, winds, and statistics (all diagnosed from the current mean state), and the **3-hourly** fast loop that steps only the column budgets. It couples to ocean, sea ice, and land — all daily by default — by exchanging **daily-mean deterministic fields** (SLP, winds, stress, precipitation, cloud, temperature, humidity, surface radiation). There is intra-seasonal memory in the prognostic columns, but **no intrinsic variability generator**: any interannual atmospheric variability is inherited entirely from the slowly-evolving boundary and the prescribed seasonal cycle.

---

# The three structural generality gaps

## Gap 1 — Vertical realism

Mid- and upper-tropospheric temperature and humidity are algebraic functions of the surface anchor `tam` and the diagnosed lapse rates. Consequences: upper-level temperature can never decouple from the surface (no independent response of static stability to remote or elevated heating), and the vertical shear — hence the jets — is only as good as the prescribed lapse-rate diagnosis. Lapse-rate and water-vapour feedbacks are therefore semi-prescribed, and **ECS is a knob** (`ecs_scale`) rather than an emergent property — a fundamental generality gap for paleoclimate and high-CO~2~ states. The broadband radiation compounds this: spectral overlap and the CO~2~ logarithmic behaviour are approximated, requiring patches at high concentrations.

## Gap 2 — Circulation and longitudinal structure

The mean meridional circulation is imposed, not solved: cell number, position, and the qualitative response to very different climates are fixed by the template. Longitudinal structure is thin — the stationary-wave solve is barotropic, orography-only, and meridionally uncoupled, with no response to land–sea diabatic-heating contrasts. Regional features (realistic monsoon lows, blocking, heating-driven stationary waves) are largely absent. The scalar eddy closure carries no information about anisotropy, tilt, or the separation of heat and momentum fluxes, so it cannot represent the up-gradient eddy-momentum transport that positions the midlatitude jet and drives the Ferrel cell.

## Gap 3 — Internal variability

There is no stochastic term anywhere in the atmospheric source. SESAM cannot produce NAO/annular-mode, PNA, blocking, or storm-track variability, and cannot seed ENSO through stochastic wind bursts. The **Hasselmann mechanism** — atmospheric white-noise fluxes integrated by the ocean mixed layer into red-noise SST and low-frequency climate variability — is absent, so mixed-layer and sea-ice variance and much interannual–decadal variability are structurally underestimated. With no entropy source, initial-condition ensembles are degenerate (no separation of forced signal from internal variability), and there are no extremes. Symptomatically, the coupled model injects noise *downstream* in the ocean surface fluxes to obtain AMOC/convection variability — the correct physical location (atmospheric eddies and winds) is silent, inviting double-counting.

---

# Roadmap toward GCM-like generality

Proposals are ranked by generality-gain ÷ cost/risk. Several are cheap because the architecture already provides the needed infrastructure: layer-resolved transport fluxes (`adifa.f90`), an unconditionally-stable implicit-diffusion pattern (`diffuse_impl.f90`), the implicit synoptic-energy solver (`synop.f90`) as a hook, and existing conservation guards (`check_energy`, `check_water`) plus radiative kernels and a partial-radiative-perturbation feedback analysis (`feedbacks.f90`, `rad_kernels.f90`) to **verify** GCM-like feedbacks at every step.

| # | Change | Unlocks | Cost / Risk |
|---|--------|---------|-------------|
| **1** | **Stochastic `sam` (SPPT on eddy energy):** AR(1) multiplicative noise on the synoptic production term, calibrated to the *variance* of the already-loaded ERA-Interim EKE climatology. One perturbation propagates coherently into all eddy transports, winds, stress, and precipitation. | Internal variability; extremes; Hasselmann mechanism; retires the ocean-noise crutch | Very low / low |
| **2** | **Second prognostic thermal level:** a mid/upper-tropospheric temperature with its own dry-static-energy budget; lapse rates become diagnosed *from two anchors*. | Static stability, jets, tropopause decouple from surface; better baroclinicity and lapse-rate feedback | Moderate / moderate |
| **3** | **Heating-forced 2-D stationary waves:** generalize `azslp` to a meridionally-coupled, diabatic-heating-forced linear solve (spectral in longitude, banded in latitude). | Longitudinal/regional skill: troughs, ridges, monsoon lows | Moderate / moderate |
| **4** | **Solved mean meridional circulation:** replace the three-cell template with a zonally-averaged QG-ω balance forced by diabatic + eddy heat flux. | Emergent cell number/position, double-ITCZ, physical Hadley widening, eddy-driven Ferrel cell | Moderate / mod-high |
| **5** | **Prognostic cloud condensate:** derive optical thickness from liquid/ice water path via the existing supersaturation overflow, phase by temperature. | Emergent (not tuned) cloud-optical feedback and its state-dependence | Moderate / moderate |
| **6** | **Two-layer clouds:** separate low (SW-dominated) and high cirrus (LW-dominated) decks, each with its own height and optical depth. | Correct cloud-radiative-effect partition and feedback sign | Moderate / moderate |
| **7** | **Richer eddy closure:** anisotropic diffusivity plus a parameterized eddy-momentum-flux convergence acting on the mean jet. | Eddy-driven jet; correct (up-gradient) momentum transport | Moderate / moderate |
| **8** | **Few-band radiation with explicit CO~2~/H~2~O**, calibrated offline to line-by-line / RRTMG references. | ECS *derived*; removes the high-CO~2~ patch and the `ecs_scale` crutch; generality across climate states | Higher / higher |
| **9** | **Lower-cost cleanups:** virtual-temperature hypsometry (consistent density/pressure); Monin–Obukhov boundary-layer stability; stochastic wind/stress with sub-daily (or variance-carrying) coupling; per-member RNG seeding for reproducible ensembles. | Consistency, coupling fidelity, ensemble capability | Low / low |

## Recommended sequencing

Two mutually reinforcing tracks give the most GCM-like generality per unit effort while protecting speed.

**Variability track (#1, plus the noise items in #9).** The cheapest and highest-novelty work, directly calibratable against the EKE climatology already in the code. It converts SESAM from weather-free to variability-capable and lets the ocean-noise workaround be replaced by a physically-sited source. Two design constraints matter: apply the perturbation on the right-hand side of the implicit `sam` solve so stability is preserved and clip to non-negative energy; and **precompute the daily correlated-noise field on a single thread** so OpenMP results remain reproducible for ensembles. Use the existing `check_energy`/`check_water` diagnostics as guardrails, with a global fixer and a taper near the humidity caps for the SPPT-on-tendencies variant.

**Structure track (#2 $\rightarrow$ #3 $\rightarrow$ #4).** Add the second thermal level first — it improves the shear that the stationary-wave and MMC upgrades depend on — then heating-forced 2-D stationary waves for regional skill, then the solved MMC for regime flexibility. These three are mutually reinforcing and together move the dynamical core substantially toward GCM-like generality while retaining the diagnostic-balance speed.

**Clouds and radiation (#5, #6, #8)** dominate feedback realism but carry the heaviest retuning burden; sequence them after the structure track and validate every step with the existing kernel/feedback machinery against GCM kernels.

## Why the speed is preserved

Every proposal keeps the core architecture: a small number of 2-D prognostics advanced by unconditionally-stable implicit solvers, with balanced/diagnostic reconstruction of the vertical and momentum fields. The added prognostics (one thermal level, one eddy-noise field, cloud condensate) are each a single 2-D budget reusing existing transport machinery; the new elliptic solves (stationary waves, MMC) are small banded/tridiagonal systems; and correlated-noise schemes add negligible cost. The result targets GCM-like *generality and variability* without adopting a GCM's prognostic 3-D dynamical core.

---

*Prepared from a three-part source study of `src/atm/` (dynamical core; column physics; variability and coupling). All file and line references point to the CLIMBER-X `dev` branch.*
