module lndvc_def
    !
    ! Data model for the virtual-cell framework (see docs/design/virtual-cells.md).
    !
    ! A virtual cell (vc) is a sub-grid patch sample: one surface class at one
    ! fixed elevation band, with an area weight relative to its coarse cell.
    ! State is decomposed into per-class blocks; a vc allocates only the blocks
    ! its class needs (so an ice/smb-only or land-only configuration falls out
    ! naturally). Column depth (soil/ice/lake layers) is internal to a block and
    ! is NOT a virtual dimension — only the surface (elevation x class) is.
    !
    ! Fields are ported faithfully from lnd_def::lnd_2d_class and
    ! lndvc_smb_class to keep functionality close to the reference lnd/smb.
    ! Per-tile arrays remain allocatable(:) and are sized per class: a land vc
    ! sizes them over the vegetation tiles (nveg), an ice/lake vc over 1.
    ! NOTE: nwiso water-isotope arrays are deferred to the physics-port phase.

    use precision, only : sp, dp, wp

    implicit none

! ============================================================================
!  Virtual-cell descriptor (geometry) — fixed elevation edges
! ============================================================================

    type vc_desc_t
        integer  :: class           ! surface class: 1 land, 2 lake, 3 ice, 4 shelf/ocean
        real(wp) :: z               ! band elevation [m] (fixed)
        real(wp) :: dz              ! band width [m]
        real(wp) :: w               ! area weight relative to coarse cell (Sum w = 1)
        real(wp) :: dz_dx, dz_dy    ! surface slope components (for orographic downscaling)
        real(wp) :: grad            ! |grad z|
        real(wp) :: z_sur_std       ! sub-grid surface elevation std [m] (snow/albedo/orog)
        real(wp) :: dz_sur          ! sub-grid elevation range for precip downscaling [m]
        real(wp) :: f_ele           ! elevation-desertification factor for precip [/]
    end type

! ============================================================================
!  Downscaled per-vc atmospheric forcing (scalar at the vc elevation)
! ============================================================================

    type forcing_t
        ! downscaled per-vc forcing (SEMI computes these internally, pattern A)
        real(wp) :: coszm, daylength
        real(wp) :: tatm, t2m, qatm, q2m
        real(wp) :: lwdown, swnet, swnet_min
        real(wp) :: rain, snow
        real(wp) :: wind, pressure
        real(wp) :: disturbance
        ! reference elevation the _i forcing below is valid at [m]
        ! (coarse-cell mean; SEMI downscales from here to vc%desc%z)
        real(wp) :: z_sur_i
        ! reference-level (_i) forcing SEMI needs (populated by cmn_to_lndvc).
        ! Ported 1:1 from the reference smb% *_i fields (see src/smb/semi.f90).
        real(wp) :: tam_i, t2m_bias_i, gam_i, tstd_i, ram_i
        real(wp) :: u700_i, v700_i, wind_i, prc_i, prc_bias_i
        real(wp) :: alb_vis_dir_i, alb_nir_dir_i, alb_vis_dif_i, alb_nir_dif_i
        real(wp) :: swd_sur_vis_dir_i, swd_sur_nir_dir_i, swd_sur_vis_dif_i, swd_sur_nir_dif_i
        real(wp) :: dswd_dalb_vis_dir_i, dswd_dalb_nir_dir_i, dswd_dalb_vis_dif_i, dswd_dalb_nir_dif_i
        real(wp) :: dswd_dz_nir_dir_i, dswd_dz_nir_dif_i
        real(wp) :: dust_i, coszm_i, swd_toa_i, swd_toa_min_i, cld_i
        real(wp) :: lwdown_i, gam_lw_i
        real(wp) :: dTvar                       ! artificial interannual T variability [K] (0 at identity)
        real(wp) :: f_ice, alb_ice              ! sub-grid ice fraction + ice background albedo
    end type

! ============================================================================
!  Shared surface flux block — the coupling currency
!  allocatable(:) sized per class (nveg for land, 1 for ice/lake)
! ============================================================================

    type surface_flux_t
        real(wp), allocatable, dimension(:) :: rough_m, rough_h, Ch, z0m, Ri
        real(wp), allocatable, dimension(:) :: r_a, r_s, beta_s
        real(wp), allocatable, dimension(:) :: r_a_can, r_s_can, beta_s_can
        real(wp), allocatable, dimension(:) :: albedo, alb_vis_dir, alb_vis_dif, alb_nir_dir, alb_nir_dif
        real(wp), allocatable, dimension(:) :: flx_sh, flx_lh, flx_g, dflxg_dT, flx_melt, flx_lwu, lwnet
        real(wp), allocatable, dimension(:) :: t_skin, t_skin_old, t_skin_amp
        real(wp), allocatable, dimension(:) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
        real(wp), allocatable, dimension(:) :: f_sh, f_e, f_t, f_le, f_lt, f_lw, lh_ecan, qsat_e, dqsatdT_e, qsat_t, dqsatdT_t
        real(wp), allocatable, dimension(:) :: transpiration, evap_surface, et
        real(wp), allocatable, dimension(:) :: rain_ground, evap_can, snow_ground, subl_can
        real(wp), allocatable, dimension(:) :: w_can, w_can_old, s_can, s_can_old, f_wat_can, f_snow_can
        real(wp), allocatable, dimension(:) :: frac_surf     ! sub-tile fractions within this vc
        ! water isotopes (sub-tile, nwiso)
        real(wp), allocatable, dimension(:,:) :: rain_iso, snow_iso
        real(wp), allocatable, dimension(:,:) :: rain_ground_iso, snow_ground_iso
        real(wp), allocatable, dimension(:,:) :: evap_can_iso, subl_can_iso
        real(wp), allocatable, dimension(:,:) :: transpiration_iso, evap_surface_iso, et_iso
        real(wp), allocatable, dimension(:,:) :: w_can_iso, w_can_iso_old, s_can_iso, s_can_iso_old
    end type

! ============================================================================
!  Shared snowpack block — one snow model for every snow-bearing class
! ============================================================================

    type snowpack_t
        integer  :: mask_snow
        real(wp) :: f_snow, h_snow
        real(wp) :: w_snow, w_snow_max, w_snow_old
        real(wp) :: snowmelt, icemelt, icesub
        real(wp) :: refreezing
        real(wp) :: refreezing_sum      ! refreezing-capacity accumulator (carry-over)
        real(wp) :: dt_snowfree         ! time since snow-free [s] (carry-over)
        real(wp) :: snow_grain, dust_con
        real(wp) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
        ! water isotopes (nwiso) — one snowpack per vc
        real(wp), allocatable, dimension(:) :: w_snow_iso, w_snow_iso_old
        real(wp), allocatable, dimension(:) :: snowmelt_iso, icemelt_iso, icesub_iso
    end type

! ============================================================================
!  Soil column (land vc): thermal + permafrost + hydrology + water balance
!  Permafrost lives here (alt diagnosed from the thermal profile).
! ============================================================================

    type soil_col_t
        ! thermal
        real(wp), allocatable, dimension(:) :: t_soil, t_soil_old, t_soil_max
        real(wp), allocatable, dimension(:) :: theta_w, theta_i, theta, w_w, w_i, w_w_old, w_i_old, w_w_phase, w_i_phase
        real(wp), allocatable, dimension(:) :: lambda_soil, lambda_int_soil, cap_soil
        real(wp), allocatable, dimension(:) :: lambda_s, lambda_dry, kappa_int
        real(wp), allocatable, dimension(:) :: theta_sat, k_sat, psi_sat, theta_field, theta_wilt
        real(wp), allocatable, dimension(:) :: t_soil_cum, theta_w_cum, theta_i_cum
        real(wp), allocatable, dimension(:) :: psi
        integer,  allocatable, dimension(:) :: k_exp, psi_exp
        ! permafrost state
        real(wp) :: alt                                        ! active layer thickness [m]
        real(wp), allocatable, dimension(:) :: frozen_years, thaw_timer
        ! hydrology / water balance
        real(wp) :: infiltration, w_table, w_table_peat
        real(wp) :: w_table_cum, w_table_min
        real(wp) :: f_wet, f_wet_cum, f_wet_max, f_wetland, cti_lim
        real(wp) :: f_wet_mon, w_table_mon, f_wet_long
        real(wp), allocatable, dimension(:) :: runoff, runoff_sur, calving, drainage, water_cons
        real(wp) :: runoff_ann
        real(wp) :: pet, mcwd, mcwd_clim
        real(wp), allocatable, dimension(:) :: wilt, root_frac
        ! water isotopes
        real(wp), allocatable, dimension(:,:) :: w_w_iso, w_i_iso, w_w_iso_old, w_i_iso_old   ! (nl,nwiso)
        real(wp), allocatable, dimension(:)   :: infiltration_iso                              ! (nwiso)
        real(wp), allocatable, dimension(:)   :: runoff_iso, runoff_sur_iso, drainage_iso, calving_iso  ! (nwiso)
        real(wp), allocatable, dimension(:)   :: water_iso_cons                                ! (nwiso)
    end type

! ============================================================================
!  Vegetation (land vc): PFT tiling stays inside the land vc
! ============================================================================

    type veg_t
        real(wp), allocatable, dimension(:) :: ci, g_can, gpp, npp, npp13, npp14, aresp
        real(wp), allocatable, dimension(:) :: discrimination
        real(wp), allocatable, dimension(:) :: lai, sai, phen, phen_acc, gdd, gamma_leaf, lambda, lai_bal
        real(wp), allocatable, dimension(:) :: npp_cum, npp13_cum, npp14_cum
        real(wp), allocatable, dimension(:) :: npp_ann, npp13_ann, npp14_ann
        real(wp), allocatable, dimension(:) :: veg_c, veg_h, pft_frac, seed_frac
        real(wp), allocatable, dimension(:) :: veg_c_below, veg_c13_below, veg_c14_below
        real(wp), allocatable, dimension(:) :: leaf_c, stem_c, root_c, veg_c13, veg_c14
        real(wp), allocatable, dimension(:) :: fire_c_flux_pft, fire_c13_flux_pft, fire_c14_flux_pft
        real(wp), allocatable, dimension(:) :: gamma_fire, gamma_fire_cum, gamma_luc, gamma_ice, gamma_dist, gamma_dist_cum
        real(wp) :: gdd5, gdd5_temp, npp_real, npp13_real, npp14_real
        real(wp) :: veg_c_above, veg_c13_above, veg_c14_above
        real(wp) :: theta_fire_cum, fuel, f_fire_fuel, f_fire_cwd
        real(wp) :: t2m_min_mon, t2m_ann_mean
        real(wp) :: f_crop, f_pasture, df_crop, df_pasture
    end type

! ============================================================================
!  Soil carbon (land vc): mineral/peat pools + decomposition rates
!  (shelf/ice/lake carbon pools travel with the ice/lake blocks or here as
!   the class dictates; kept together for now to mirror lnd_2d_class).
! ============================================================================

    type soil_carbon_t
        real(wp), allocatable, dimension(:) :: ftemp, fmoist, fdepth
        real(wp), allocatable, dimension(:) :: k_litter, k_fast, k_slow, k_litter_wet, k_fast_wet, k_slow_wet, diff_soilc, adv_soilc
        real(wp), allocatable, dimension(:) :: k_slow_to_fast          ! permafrost-thaw slow->fast conversion
        real(wp), allocatable, dimension(:) :: k_cato, ch4_frac_wet, ch4_frac_peat, ch4_frac_shelf, ch4_frac_lake
        real(wp), allocatable, dimension(:) :: frac_soc
        real(wp), allocatable, dimension(:) :: litter_c, fast_c, slow_c, litter_c13, fast_c13, slow_c13, litter_c14, fast_c14, slow_c14
        real(wp), allocatable, dimension(:) :: cato_c, cato_c13, cato_c14
        real(wp), allocatable, dimension(:) :: soil_c_tot, soil_resp, soil_c13_tot, soil_resp13, soil_c14_tot, soil_resp14
        real(wp), allocatable, dimension(:) :: litterfall, litterfall13, litterfall14
        real(wp), allocatable, dimension(:) :: litter_in_frac
        real(wp), allocatable, dimension(:) :: soil_resp_l
        ! peat
        real(wp) :: k_litter_peat, k_acro, k_litter_peat_anox, k_acro_anox, f_oxic_peat
        real(wp) :: litter_c_peat, acro_c, litter_c13_peat, acro_c13, litter_c14_peat, acro_c14
        real(wp) :: f_peat, f_peat_pot, acro_h, cato_h, peat_c_ini_year, dCpeat_dt
        ! emissions / weathering / dust (per-vc contributions; aggregated to cell)
        real(wp) :: ch4_emis_wetland, ch4_emis_shelf, ch4_emis_peat, ch4_emis_lake
        real(wp) :: c13h4_emis_wetland, c13h4_emis_shelf, c13h4_emis_peat, c13h4_emis_lake
        real(wp) :: n2o_emis
        real(wp) :: dust_emis_d, dust_emis_g, dust_emis_s, dust_emis, dust_dep
        real(wp) :: f_carb
        real(wp) :: weath_carb, weath_sil, weath_loess, weath13_carb, weath13_sil, weath14_carb, weath14_sil
        real(wp) :: poc_export, poc13_export, poc14_export, doc_export, doc13_export, doc14_export
        real(wp) :: lithology_gemco2, lithology_uhh, lithology_shelf_uhh
        ! conservation
        real(wp) :: carbon_cons_soil, carbon13_cons_soil, carbon14_cons_soil
    end type

! ============================================================================
!  Ice/firn column (ice vc): thermal + snow-on-ice + SMB mass-budget diagnostic
! ============================================================================

    type ice_col_t
        ! snow+ice/firn thermal profile SEMI integrates (0:nl; index 0 = skin layer)
        real(wp), allocatable, dimension(:) :: t_prof, t_prof_old
        real(wp), allocatable, dimension(:) :: t_ice, t_ice_old
        real(wp), allocatable, dimension(:) :: lambda_ice, lambda_int_ice, cap_ice
        ! subglacial carbon pools (mirror lnd ice/shelf carbon)
        real(wp), allocatable, dimension(:) :: litter_c_ice, fast_c_ice, slow_c_ice
        real(wp), allocatable, dimension(:) :: litter_c_shelf, fast_c_shelf, slow_c_shelf
        real(wp), allocatable, dimension(:) :: cap_shelf, lambda_int_shelf
        real(wp), allocatable, dimension(:) :: t_shelf, t_shelf_old, t_shelf_max
        real(wp), allocatable, dimension(:) :: theta_w_shelf, theta_i_shelf, w_w_shelf, w_i_shelf
        ! SMB diagnostic (mass budget of the shared snow+ice column)
        real(wp) :: smb                 ! surface mass balance [kg/m2/s]
        real(wp) :: melt                ! snow+ice melt [kg/m2/s]
        real(wp) :: runoff              ! snowmelt + icemelt + rain - refreezing [kg/m2/s]
        real(wp) :: f_rfz_to_snow
        real(wp) :: energy_cons_ice, energy_cons_shelf
    end type

! ============================================================================
!  Lake column (lake vc)
! ============================================================================

    type lake_col_t
        real(wp), allocatable, dimension(:) :: t_lake, t_lake_old
        real(wp), allocatable, dimension(:) :: lambda_lake, lambda_int_lake, cap_lake
        real(wp), allocatable, dimension(:) :: t_sublake, lambda_sublake, lambda_int_sublake, cap_sublake
        real(wp), allocatable, dimension(:) :: theta_w_sublake, theta_i_sublake, w_w_sublake, w_i_sublake
        real(wp), allocatable, dimension(:) :: w_w_lake, w_i_lake, f_i_lake
        real(wp), allocatable, dimension(:) :: t_sublake_cum, theta_w_sublake_cum, theta_i_sublake_cum
        real(wp), allocatable, dimension(:) :: t_shelf_cum, theta_w_shelf_cum, theta_i_shelf_cum
        real(wp) :: h_lake, h_lake_conv, h_lake_mix, f_lake_ice, lake_water_tendency
        real(wp) :: energy_cons_lake
        ! water isotopes (nl_l,nwiso)
        real(wp), allocatable, dimension(:,:) :: w_w_lake_iso, w_i_lake_iso
    end type

! ============================================================================
!  One virtual cell — blocks allocated by class
! ============================================================================

    type vc_t
        type(vc_desc_t)     :: desc
        type(forcing_t)     :: forc
        type(surface_flux_t):: flx
        type(snowpack_t),   allocatable :: snow
        type(soil_col_t),   allocatable :: soil
        type(veg_t),        allocatable :: veg
        type(soil_carbon_t),allocatable :: carb
        type(ice_col_t),    allocatable :: ice
        type(lake_col_t),   allocatable :: lake
        ! per-vc conservation residuals
        real(wp) :: energy_cons_surf1, energy_cons_surf2, energy_cons_soil
    end type

! ============================================================================
!  Aggregated coarse-cell state — thin, the drop-in for coupling
! ============================================================================

    type vc_cell_t
        integer  :: mask_lnd
        ! class fractions (derived from summing vc area weights by class)
        real(wp) :: f_land, f_land0, f_ice, f_ice_old, f_ice_grd, f_ice_grd_old, f_ice_nbr
        real(wp) :: f_shelf, f_shelf_old, f_lake, f_lake_old, f_veg, f_veg_old
        real(wp) :: z_veg, z_veg_std, z_veg_min, z_veg_max
        ! aggregated coupling currency (area-weighted means of surface_flux_t)
        real(wp) :: t_skin, albedo, flx_sh, flx_lh, flx_g, runoff, evap, et
        ! aggregated ice-surface mass budget (area-weighted cell-mean fluxes)
        real(wp) :: smb, melt
        ! cell-level carbon fluxes / emissions
        real(wp) :: Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd
        real(wp) :: fire_c_flux, fire_c13_flux, fire_c14_flux
        real(wp) :: ch4_emis, n2o_emis
    end type

! ============================================================================
!  Global 0-D state (mirror lnd_0d_class)
! ============================================================================

    type veg_glob_t
        real(wp) :: co2, c13_c12_atm, c14_c_atm
        real(wp) :: Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd
        real(wp) :: Cflx_avg
        real(wp) :: Cflx_burial, C13flx_burial, C14flx_burial
        real(wp) :: ch4_emis, n2o_emis
        real(wp) :: landc, landc13, landc14
        real(wp) :: burc, burc13, burc14
        real(wp) :: weath_scale
        real(wp) :: weath_carb_avg, weath_sil_avg
    end type

! ============================================================================
!  Framework container
! ============================================================================

    type lndvc_class
        ! grid / index maps
        integer               :: ncells
        integer,  allocatable :: id_map(:,:)
        integer,  allocatable :: ij_1d(:,:)
        ! decomposition metadata (fixed elevation edges)
        integer               :: n_vc              ! max leaf virtual cells per coarse cell
        real(wp), allocatable :: z_edges(:)        ! fixed band edges [m]
        integer,  allocatable :: leaf_1d(:,:)      ! compressed active-leaf list (cell,vc) for OMP
        real(wp), allocatable :: z_lake(:)
        ! leaf virtual cells: (nx, ny, n_vc)
        type(vc_t),      allocatable :: vc(:,:,:)
        ! aggregated coarse-cell state for coupling: (nx, ny)
        type(vc_cell_t), allocatable :: cell(:,:)
        ! global 0-D
        type(veg_glob_t) :: glob
    end type

    public

contains

end module lndvc_def
