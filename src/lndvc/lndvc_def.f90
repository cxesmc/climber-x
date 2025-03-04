module lndvc_def

    use precision, only : sp, dp, wp

    implicit none

! === VIRTUAL CELL CLASSES =======

    type lndvc_vc_class
        ! Definition of characteristics of a given virtual cell
        integer  :: surf            ! Surface type (0: ocean, 1: land, 2: lake, 3: ice)
        real(wp) :: zsrf            ! Surface elevation (m)
        real(wp) :: f_land, f_ice
    end type

! === VEGETATION MODEL VARIABLES ====

    type lndvc_veg_global_class
        real(wp) :: co2, c13_c12_atm, c14_c_atm
        real(wp) :: Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd
        real(wp) :: Cflx_avg
        real(wp) :: Cflx_burial, C13flx_burial, C14flx_burial
        real(wp) :: ch4_emis
        real(wp) :: landc, landc13, landc14
        real(wp) :: burc, burc13, burc14
        real(wp) :: weath_scale
        real(wp) :: weath_carb_avg, weath_sil_avg
    end type

    type lndvc_veg_class
        
        integer  :: mask_lnd
        real(wp) :: Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd
        real(wp) :: f_land, f_land0, f_ice, f_ice_old, f_ice_grd, f_ice_grd_old, f_ice_nbr, f_shelf, f_shelf_old, f_lake, f_lake_old, f_veg, f_veg_old
        real(wp) :: z_veg_std, z_veg, z_veg_min, z_veg_max
        real(wp) :: f_crop, f_pasture
        real(wp) :: t2m_min_mon
        real(wp) :: t2m_ann_mean
        real(wp) :: infiltration, w_table, w_table_peat
        real(wp) :: f_wet, f_wet_cum, f_wet_max, f_wetland, w_table_cum, w_table_min
        real(wp) :: t_skin_veg, flx_g_veg, dflxg_dT_veg, flx_melt_veg, t_2m
        real(wp) :: f_lake_ice
        real(wp) :: h_lake
        real(wp) :: h_lake_conv
        real(wp) :: h_lake_mix
        real(wp) :: lake_water_tendency
        real(wp) :: gdd5, gdd5_temp, npp_real, npp13_real, npp14_real
        real(wp) :: veg_c_above, veg_c13_above, veg_c14_above
        real(wp) :: theta_fire_cum
        real(wp) :: alt
        real(wp) :: k_litter_peat, k_acro, k_litter_peat_anox, k_acro_anox, f_oxic_peat
        real(wp) :: litter_c_peat, acro_c, litter_c13_peat, acro_c13, litter_c14_peat, acro_c14
        real(wp) :: f_peat, f_peat_pot, acro_h, cato_h, peat_c_ini_year, dCpeat_dt
        real(wp) :: ch4_emis_wetland, ch4_emis_shelf, ch4_emis_peat, ch4_emis_lake
        real(wp) :: c13h4_emis_wetland, c13h4_emis_shelf, c13h4_emis_peat, c13h4_emis_lake
        real(wp) :: dust_emis_d, dust_emis_g, dust_emis_s, dust_emis
        real(wp) :: dust_dep
        real(wp) :: runoff_ann
        real(wp) :: f_carb
        real(wp) :: weath_carb, weath_sil, weath_loess
        real(wp) :: weath13_carb, weath13_sil
        real(wp) :: weath14_carb, weath14_sil
        real(wp) :: poc_export, poc13_export, poc14_export
        real(wp) :: doc_export, doc13_export, doc14_export
        real(wp) :: energy_cons_soil, energy_cons_ice, energy_cons_shelf, energy_cons_lake
        real(wp) :: carbon_cons_veg, carbon13_cons_veg, carbon14_cons_veg

        ! Variables that are calculated as virtual cell variables:
        ! (used to be defined with "allocatable, dimension(:)" here)

        real(wp) :: frac_surf
        real(wp) :: disturbance
        real(wp) :: coszm, daylength
        real(wp) :: tatm, t2m, qatm, q2m, lwdown, swnet, swnet_min
        real(wp) :: rain, snow
        real(wp) :: wind
        real(wp) :: pressure
        real(wp) :: rough_m, rough_h, Ch, z0m, Ri
        real(wp) :: r_a, r_s, beta_s
        real(wp) :: r_a_can, r_s_can, beta_s_can
        real(wp) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
        real(wp) :: albedo, alb_vis_dir, alb_vis_dif, alb_nir_dir, alb_nir_dif
        integer,  allocatable :: mask_snow(:)
        real(wp) :: f_snow, h_snow, w_snow, w_snow_max, w_snow_old, snowmelt, icemelt, snow_grain, dust_con
        real(wp) :: runoff, runoff_sur, calving, drainage, water_cons
        real(wp) :: rain_ground, evap_can, snow_ground, subl_can
        real(wp) :: w_can, w_can_old, s_can, s_can_old, f_snow_can
        real(wp) :: transpiration, evap_surface, et
        real(wp) :: f_wet_mon, w_table_mon
        real(wp) :: f_wet_long
        real(wp) :: flx_sh, flx_lh, flx_g, dflxg_dT, flx_melt, flx_lwu, lwnet
        real(wp) :: t_skin, t_skin_old, t_skin_amp
        real(wp) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
        real(wp) :: f_sh, f_e, f_t, f_le, f_lt, f_lw, lh_ecan, qsat_e, dqsatdT_e, qsat_t, dqsatdT_t
        real(wp) :: ci, g_can, gpp, npp, npp13, npp14, aresp
        real(wp) :: discrimination
        real(wp) :: lai, sai, phen, phen_acc, gdd, gamma_leaf, lambda, lai_bal
        real(wp) :: npp_cum, npp13_cum, npp14_cum
        real(wp) :: npp_ann, npp13_ann, npp14_ann
        real(wp) :: veg_c, veg_h, pft_frac, seed_frac
        real(wp) :: veg_c_below, veg_c13_below, veg_c14_below
        real(wp) :: leaf_c, stem_c, root_c, veg_c13, veg_c14
        real(wp) :: gamma_luc, gamma_ice, gamma_dist, gamma_dist_cum
        real(wp) :: lambda_soil, lambda_int_soil, cap_soil
        real(wp) :: lambda_ice, lambda_int_ice, cap_ice
        real(wp) :: lambda_lake, lambda_int_lake, cap_lake 
        real(wp) :: lambda_sublake, lambda_int_sublake, cap_sublake
        real(wp) :: cap_shelf
        real(wp) :: lambda_int_shelf
        real(wp) :: lambda_s, lambda_dry
        real(wp) :: kappa_int
        real(wp) :: theta_sat, k_sat, psi_sat, theta_field, theta_wilt
        real(wp) :: t_soil, t_soil_old, t_soil_max
        real(wp) :: t_ice, t_ice_old
        real(wp) :: t_shelf, t_shelf_old, t_shelf_max
        real(wp) :: t_lake, t_lake_old
        real(wp) :: t_sublake
        real(wp) :: theta_w, theta_i, theta, w_w, w_i, w_w_old, w_i_old, w_w_phase, w_i_phase
        real(wp) :: theta_w_shelf, theta_i_shelf, w_w_shelf, w_i_shelf
        real(wp) :: theta_w_sublake, theta_i_sublake, w_w_sublake, w_i_sublake
        real(wp) :: w_w_lake, w_i_lake, f_i_lake
        real(wp) :: t_soil_cum, theta_w_cum, theta_i_cum
        real(wp) :: t_shelf_cum, theta_w_shelf_cum, theta_i_shelf_cum
        real(wp) :: t_sublake_cum, theta_w_sublake_cum, theta_i_sublake_cum
        real(wp) :: psi
        integer  :: k_exp, psi_exp
        real(wp) :: ftemp, fmoist, fdepth
        real(wp) :: k_litter, k_fast, k_slow, k_litter_wet, k_fast_wet, k_slow_wet, diff_soilc, adv_soilc
        real(wp) :: k_litter_shelf, k_fast_shelf, k_slow_shelf, diff_shelfc, adv_shelfc
        real(wp) :: k_litter_lake, k_fast_lake, k_slow_lake, diff_lakec, adv_lakec
        real(wp) :: k_litter_ice, k_fast_ice, k_slow_ice, diff_icec, adv_icec
        real(wp) :: k_cato, ch4_frac_wet, ch4_frac_peat, ch4_frac_shelf, ch4_frac_lake
        real(wp) :: frac_soc
        real(wp) :: litter_c, fast_c, slow_c, litter_c13, fast_c13, slow_c13, litter_c14, fast_c14, slow_c14
        real(wp) :: litter_c_shelf, fast_c_shelf, slow_c_shelf, litter_c13_shelf, fast_c13_shelf, slow_c13_shelf, litter_c14_shelf, fast_c14_shelf, slow_c14_shelf
        real(wp) :: litter_c_lake, fast_c_lake, slow_c_lake, litter_c13_lake, fast_c13_lake, slow_c13_lake, litter_c14_lake, fast_c14_lake, slow_c14_lake
        real(wp) :: litter_c_ice, fast_c_ice, slow_c_ice, litter_c13_ice, fast_c13_ice, slow_c13_ice, litter_c14_ice, fast_c14_ice, slow_c14_ice
        real(wp) :: cato_c, cato_c13, cato_c14
        real(wp) :: soil_c_tot, soil_resp, soil_c13_tot, soil_resp13, soil_c14_tot, soil_resp14
        real(wp) :: litter_in_frac
        real(wp), allocatable :: litterfall(:), litterfall13(:), litterfall14(:)
        real(wp), allocatable :: wilt(:)
        real(wp), allocatable :: root_frac(:)
        real(wp), allocatable :: soil_resp_l(:)

        ! === end virtual cell variables ====

        real(wp) :: lithology_gemco2
        real(wp) :: lithology_uhh
        real(wp) :: lithology_shelf_uhh

        real(wp) :: energy_cons_surf1, energy_cons_surf2
        real(wp) :: carbon_cons_soil, carbon13_cons_soil, carbon14_cons_soil
        
    end type

    type lndvc_veg_param_class
        integer :: nx
        integer :: ny 
        integer :: n_vc

        integer :: ni, nj, lat
        integer :: nl_l

    end type


! ==== SMB =======

    type lndvc_smb_param_class
        integer :: par1

    end type 


    type lndvc_smb_class

        ! ajr: these will be removed potentially, or some will be moved to lndvc_vc_class...
        ! ! input variables from ice sheet (matches smb grid)
        ! integer :: mask_ice   !! ice mask on ice sheet grid []
        ! integer :: mask_ice_old   !! old ice mask on ice sheet grid []
        ! integer :: mask_smb   !! mask where smb is computed []
        ! integer :: mask_smb_tmp   !! mask where smb is computed []
        ! integer :: mask_maxice  !! mask of maximum ice extent []
        ! integer :: mask_margin    !! ice margin mask on ice sheet grid [/]
        ! real(wp) :: z_sur      !! surface elevation on ice sheet grid [m]
        ! real(wp) :: z_sur_eff     !! effective surface elevation on ice sheet grid [m]
        ! real(wp) :: z_sur_fil     !! filtered surface elevation on ice sheet grid [m]
        ! real(wp) :: z_sur_std     !! sub-grid standard deviation of surface elevation [m]
        ! real(wp) :: z_bed_std     !! sub-grid standard deviation of bedrock elevation [m]
        ! real(wp) :: h_ice      !! ice thickness on ice sheet grid [m]

        ! boundary variables

        real(wp) :: co2         ! atmospheric CO2 concentration [ppm]
        real(wp) :: Smax65N     ! maximum summer insolation at 65N [W/m2]

        ! input variables on coupler grid
        real(wp) :: z_sur  !! grid cell mean elevation [m]
        real(wp) :: t2m    !! surface air temperature at grid cell mean elevation [K]
        real(wp), allocatable :: t2m_bias(:)    !! surface air temperature bias at present-day [K]
        real(wp), allocatable :: prc_bias(:)    !! precipitation bias at present-day (ratio model/obs) [1]
        real(wp) :: tam    !! atmospheric temperature [K] 
        real(wp) :: ram    !! atmospheric relative humidity [1] 
        real(wp) :: gam    !! atmospheric lapse rate at the surface [K/m] 
        real(wp) :: tstd   !! standard deviation on daily 2m temperature [K] 
        real(wp) :: prc    !! total precipitation rate [kg/m2/s] 
        real(wp) :: u700   !! zonal wind component at 700 hPa [m/s]
        real(wp) :: v700   !! meridional wind component at 700 hPa [m/s]
        real(wp) :: wind   !! surface wind speed [m/s]
        real(wp) :: cld    !! cloud cover fraction [/]
        real(wp) :: dust   !! dust deposition rate [kg/m2/s]
        real(wp) :: swd_toa   !! downward shortwave radiation at TOA [W/m2]
        real(wp) :: swd_toa_min   !! minimum diurnal downward shortwave radiation at TOA [W/m2]
        real(wp) :: swd_sur_vis_dir   !! downward shortwave visible radiation at surface, clear sky [W/m2]
        real(wp) :: swd_sur_nir_dir   !! downward shortwave near-infrared radiation at surface, clear sky [W/m2]
        real(wp) :: swd_sur_vis_dif   !! downward shortwave visible radiation at surface, cloudy [W/m2]
        real(wp) :: swd_sur_nir_dif   !! downward shortwave near-infrared radiation at surface, cloudy [W/m2]
        real(wp) :: dswd_dalb_vis_dir
        real(wp) :: dswd_dalb_nir_dir
        real(wp) :: dswd_dalb_vis_dif
        real(wp) :: dswd_dalb_nir_dif
        real(wp) :: dswd_dz_nir_dir
        real(wp) :: dswd_dz_nir_dif
        real(wp) :: alb_vis_dir   !! surface albedo for visible radiation at surface, clear sky []
        real(wp) :: alb_nir_dir   !! surface albedo for near-infrared radiation at surface, clear sky []
        real(wp) :: alb_vis_dif   !! surface albedo for visible radiation at surface, cloudy []
        real(wp) :: alb_nir_dif   !! surface albedo for near-infrared radiation at surface, cloudy []
        real(wp) :: lwdown   !! downward surface longwave radiation at grid-cell mean elevation [W/m2]
        real(wp) :: gam_lw   !! downward surface longwave radiation 'lapse rate' [W/m2]
        real(wp) :: t_ground  !! annual mean ground temperature (soil over land and bottom water over ocean) [degC]
        real(wp) :: coszm  !! daily mean cosine of solar zenith angle []

        ! internal variables
        real(wp) :: dz_dx_sur  !! ice surface x-slope on ice sheet grid [/]
        real(wp) :: dz_dy_sur  !! ice surface y-slope on ice sheet grid [/]
        real(wp) :: dz_sur     !! ice surface slope on ice sheet grid [/]
        real(wp) :: t_skin   !! skin temperature [K]
        real(wp) :: t_skin_amp   !! amplitude of diurnal cycle of skin temperature [K]
        real(wp) :: t_skin_old   !! old skin temperature [K]
        real(wp), allocatable :: t_prof(:)  !! snow/ice temperature [K]
        real(wp), allocatable :: t_prof_old(:)  !! old snow/ice temperature [K]
        real(wp) :: q2m    !! surface air specific humidity [kg/kg]
        real(wp) :: pressure     !! surface pressure [Pa]
        integer  :: mask_snow    !! snow mask [/]
        real(wp) :: f_snow !! snow cover fraction [1]
        real(wp) :: h_snow !! snow thickness [m]
        real(wp) :: w_snow !! snow water equivalent [kg/m2]
        real(wp) :: w_snow_old !! old snow water equivalent [kg/m2]
        real(wp) :: w_snow_max !! max snow water equivalent [kg/m2]
        real(wp) :: snowmelt !! snow melt [kg/m2/s]
        real(wp) :: icemelt !! ice melt [kg/m2/s]
        real(wp) :: melt !! ice + snow melt [kg/m2/s]
        real(wp) :: refreezing !! refreezing [kg/m2/s]
        real(wp) :: refreezing_sum !! refreezing summed up over the year [kg/m2/s]
        real(wp) :: f_rfz_to_snow !! fraction of refreezing going to snow [1]
        real(wp) :: runoff  !! runoff [kg/m2/s]
        real(wp) :: cod  !! cloud optical depth [/]
        real(wp) :: albedo  !! allwave surface albedo [/]
        real(wp) :: alb_bg  !! allwave background albedo [/]
        real(wp) :: f_ice !! ice cover fraction [1]
        real(wp) :: f_ice_old !! old ice cover fraction [1]
        real(wp) :: dt_snowfree !! cumulated time over the year when snowfree [s]
        real(wp) :: alb_ice  !! allwave ice albedo [/]
        real(wp) :: alb_snow_vis_dir   !! snow_albedo for visible radiation for clear sky interpolated to ice sheet grid [W/m2]
        real(wp) :: alb_snow_nir_dir   !! snow_albedo for near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
        real(wp) :: alb_snow_vis_dif   !! snow_albedo for visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
        real(wp) :: alb_snow_nir_dif   !! snow_albedo for near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
        real(wp) :: snow_grain  !! snow grain size [um]
        real(wp) :: dust_con  !! dust concentration in snow [kg/kg]
        real(wp) :: swdown  !! downward surface shortwave radiation [W/m2]
        real(wp) :: swnet  !! net surface shortwave radiation [W/m2]
        real(wp) :: swnet_min  !! minimum diurnal net surface shortwave radiation [W/m2]
        real(wp) :: r_a  !! aerodynamic resistance [s/m]
        real(wp) :: flx_g  !! ground heat flux [W/m2]
        real(wp) :: dflxg_dT  !! derivative of ground heat flux wrt top layer temperature [W/m2/K]
        real(wp) :: flx_melt !! heat flux used to melt snow [W/m2]
        real(wp) :: flx_sh  !! sensible heat flux [W/m2]
        real(wp) :: flx_lwu  !! upwelling longwave radiation [W/m2]
        real(wp) :: flx_lh  !! latent heat flux [W/m2]
        real(wp) :: evp   !! evaporation rate [kg/m2/s]
        real(wp) :: f_ele  !! precipitation rate, downscaled [kg/m2/s]
        real(wp) :: f_wind  !! precipitation rate, downscaled [kg/m2/s]
        real(wp) :: rain   !! rainfall rate, downscaled [kg/m2/s]
        real(wp) :: snow   !! snowfall rate, downscaled [kg/m2/s]

        real(wp) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
        real(wp) :: f_sh, f_e, f_lh, f_lw, qsat, dqsatdT

    end type

    type lndvc_class

        ! Define all parameters here
        
        ! Define all geometry related information here
        ! (topographies, masks, maps, etc)
        type(lndvc_vc_class), allocatable :: vc(:,:,:)

        integer :: ncells
        integer,  allocatable :: id_map(:,:)
        integer,  allocatable :: ij_1d(:,:)
        real(wp), allocatable :: z_lake(:)      ! par%nl_l

        ! === VEGETATION MODEL VARIABLES ===

        type(lndvc_veg_param_class) :: veg_par
        type(lndvc_veg_global_class) :: veg_glob
        type(lndvc_veg_class), allocatable :: veg_vc(:,:,:)
        type(lndvc_veg_class), allocatable :: veg(:,:)
        
        ! === SOIL MODEL ===========



        ! === SMB MODEL (SEMIX) VARIABLES ===

        type(lndvc_smb_class), allocatable :: smb_vc(:,:,:)     ! cmn grid plus vc dimension

        type(lndvc_smb_class), allocatable :: smb(:,:)          ! smb on the cmn grid, 5x5
        type(lndvc_smb_class), allocatable :: smb_nh(:,:)       ! High-res projected grid, NH
        type(lndvc_smb_class), allocatable :: smb_sh(:,:)       ! High-res projected grid, SH

    end type


    public

contains


end module lndvc_def