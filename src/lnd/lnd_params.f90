!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l n d _ p a r a m s
!
!  Purpose : land model parameters
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
!
! This file is part of CLIMBER-X.
!
! CLIMBER-X is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! CLIMBER-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with CLIMBER-X.  If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module lnd_params

  use precision, only : wp
  use ncio
  use timer, only : dt_lnd, nmon_year, nday_year, sec_year, sec_day
  use control, only : out_dir
  use lnd_grid, only : nx, ny, npft, nl, nsurf, z, z_int, lat
  use lnd_grid, only : flag_tree, flag_grass, flag_shrub

  implicit none

  integer :: nstep_v
  logical :: time_call_carb, time_call_veg, time_call_carb_p
  integer :: dt_day_carb
  integer :: dt_day_veg 

  real(wp) :: dt, dt_day, rdt
  integer :: dt_day_c, dt_day_v
  real(wp) :: dt_c, dt_v

  character (len=256) :: lnd_par_file
  logical :: write_surf, write_surf_n, write_carbon, write_soil, write_soil_par, write_lake, write_cons, l_daily_output

  integer :: i_init_veg
  logical :: l_dynveg
  character (len=256) :: pft_fix_file
  logical :: l_fixlai
  character (len=256) :: lai_fix_file
  logical :: l_co2_fert_lim
  real(wp) :: co2_fert_lim_min
  real(wp) :: co2_fert_lim_max
  logical :: l_diurnal_cycle
  logical :: l_neutral
  real(wp) :: sw_par_frac    ! fraction of SW down radiation which is PAR
  integer :: i_ci
  integer :: i_beta
  integer :: i_vcmax
  integer :: i_disc
  real(wp) :: z_sfl
  integer :: i_racan
  real(wp) :: p_cdense
  integer :: i_roots

  character (len=256) :: lithology_uhh_file
  character (len=256) :: lithology_gemco2_file
  integer :: i_weathering
  integer :: i_weath_sc
  real(wp) :: kweath_scale
  logical :: l_river_export
  real(wp) :: kexport
  logical :: l_match_alk_weath

  real(wp) :: K_eddy_lake_bg 
  real(wp) :: K_eddy_lake_max
  real(wp) :: z_mix_lake_min
  real(wp) :: dz_lake_ice

  logical :: l_ice_albedo_semi

  ! hydrology parameters
  type hydro_par_type
    integer :: i_frz_imp
    integer :: i_cond_theta
    integer :: i_hydro_par
    integer :: i_wtab
    integer :: i_fwet
    logical :: l_dew 
    logical :: l_prc_intercept
    integer :: i_runoff
    logical :: l_runoff_icemelt
    integer :: i_evp_soil
    real(wp) :: theta_crit_evp
    real(wp) :: dz_evp
    real(wp) :: frac_inf_2
    real(wp) :: p_psi_min
    real(wp) :: p_psi_max
    real(wp) :: kappa_max = 10._wp   ! kg/m2/day
    real(wp) :: theta_min = 0.01_wp ! m3/m3
    real(wp), dimension(npft) :: alpha_int_w ! interception factor for water, CLM, TUNABLE
    real(wp) :: alpha_int_s ! interception factor for snow
    real(wp) :: can_max_w      ! kg/m2, canopy interception capacity parameter Verseghy 1991
    real(wp) :: can_max_s      ! kg/m2, canopy interception capacity parameter for snow
    real(wp) :: tau_w       ! s, removal time scale of water from canopy
    real(wp) :: tau_s       ! s, removal time scale of snow from canopy
    real(wp) :: f_wtab
    real(wp) :: cti_min
    real(wp) :: wtab_scale
    real(wp) :: fmax_crit
    real(wp) :: cti_mean_crit
  end type
  type(hydro_par_type) :: hydro_par

  ! TOPMODEL parameters
  type topmodel_type
    real(wp) :: cti_mean
    real(wp) :: cti_cdf(15)
  end type
  type(topmodel_type), dimension(nx,ny) :: topmodel

  ! dyptop parameters
  type dyptop_type
    real(wp) :: k
    real(wp) :: xm
    real(wp) :: v
    real(wp) :: fmax
  end type
  type(dyptop_type), dimension(nx,ny) :: dyptop

  ! snow parameters
  type snow_par_type
    integer :: i_snow_albedo
    logical :: f_snow_sub
    logical :: l_snow_aging
    logical :: l_snow_dust
    real(wp) :: alb_snow_vis_dif_new = 0.99_wp
    real(wp) :: alb_snow_nir_dif_new = 0.65_wp
    real(wp) :: snow_grain_fresh
    real(wp) :: snow_grain_old
    real(wp) :: f_age_t 
    real(wp) :: snow_0
    real(wp) :: snow_1
    real(wp) :: rho_snow
    real(wp) :: lambda_snow
    real(wp) :: w_snow_crit     ! kg/m2, minimum snow water equivalent for explicit snow layer
    real(wp) :: w_snow_off   ! kg/m2, maximum snow water equivalent
    logical :: l_fsnow_orog
    real(wp) :: c_fsnow
    real(wp) :: c_fsnow_orog
    real(wp) :: z_veg_std_min
  end type
  type(snow_par_type) :: snow_par

  ! surface parameters
  type surf_par_type
    real(wp) :: alb_vis_dir_ice  
    real(wp) :: alb_vis_dif_ice 
    real(wp) :: alb_nir_dir_ice 
    real(wp) :: alb_nir_dif_ice 
    ! CHECK
    real(wp) :: alb_vis_dir_water = 0.1
    real(wp) :: alb_vis_dif_water = 0.1
    real(wp) :: alb_nir_dir_water = 0.1
    real(wp) :: alb_nir_dif_water = 0.1

    real(wp) :: d_alb_age_vis = 0.05_wp
    real(wp) :: d_alb_age_nir = 0.25_wp

    real(wp) :: z0m_ice
    real(wp) :: z0m_lake
    real(wp) :: z0m_lake_ice
    real(wp) :: z0m_bare 
    real(wp) :: z0m_snow = 0.0024_wp

    integer :: i_z0h
    real(wp) :: zm_to_zh_const != exp(-2._wp)  ! LandLad, Garratt 1992

    real(wp) :: f_Ri_unstab
    real(wp) :: f_Ri_stab

    ! longwave emissivity, Jin 2006, Walters 2014           BL     NL     C3     C4     SH    BARE  WATER  ICE
    real(wp), dimension(nsurf) :: emissivity = (/ 0.96_wp,0.96_wp,0.96_wp,0.96_wp,0.96_wp,0.96_wp,0.98_wp,0.99_wp /)
    !real(wp), dimension(nsurf) :: emissivity = (/ 1._wp,1._wp,1._wp,1._wp,1._wp,1._wp,1._wp,1._wp /)
    real(wp) :: emissivity_snow = 0.99_wp
    real(wp), dimension(nx,ny) :: alb_bare_vis, alb_bare_nir

  end type
  type(surf_par_type) :: surf_par

  ! dust parameters
  type dust_par_type
    real(wp) :: b0  !! dust emission flux constant [kg s^2 m^-5]
    real(wp) :: qd  !! desert dust source scaling constant
    real(wp) :: qg  !! grass dust source scaling constant
    real(wp) :: u0  !! wind velocity threshold [m/s]
    real(wp) :: wind_gust_fac   !! factor to scale wind to account for topographic roughness-dependent wind gusts [1/1000 m]
    integer :: i_fsnow  !! snow fraction used for dust emission limitation []
    integer :: i_theta  !! soil moisture used for dust emission limitation []
    real(wp) :: sm_t   !! soil moisture of transition [m3/m3]
    real(wp) :: sm_n   !! normalisation constant
    real(wp) :: lai_t  !! grass leaf area index of transition [m2/m2]
    real(wp) :: lai_n  !! normalization constant [m2/m2]
    real(wp) :: rho_dust = 2500._wp  !! dust density [kg/m3]
    logical :: l_dust_stab
    logical :: l_dust_topo
    integer :: topo_exp
  end type
  type(dust_par_type) :: dust_par

  ! soil parameters
  type soil_par_type
    integer :: soil_texture
    logical :: uniform_porosity, uniform_soil_par_therm, uniform_soil_par_hydro
    logical :: constant_porosity, constant_soil_par_therm, constant_soil_par_hydro
    real(wp) :: cap_s = 2.3d6   ! J/m3/K, volumetric heat capacity of soil
    real(wp) :: theta_sat_u
    real(wp) :: theta_field_u
    real(wp) :: theta_wilt_u
    real(wp) :: lambda_s_u
    real(wp) :: lambda_dry_u
    real(wp) :: k_sat_u
    real(wp) :: psi_sat_u
    integer :: b_u
  end type
  type(soil_par_type) :: soil_par
  ! organic soil properties
  type organic_par_type ! values from Lawrence 2008, Letts 2000
    real(wp) :: theta_sat  = 0.7_wp    ! m3/m3
    real(wp) :: k_sat      = 0.01_wp   ! mm/s, kg/m2/s, Avis 2011
    real(wp) :: psi_sat    = -10.3d-3 ! m, -0.1 for mineral soil
    real(wp) :: Bi         = 4._wp     ! between 2.7 and 12 in Letts 2010, 2.7 used by lawrence 2008 
    real(wp) :: lambda_s   = 0.25_wp   ! W/m/K 
    real(wp) :: lambda_dry = 0.05_wp   ! W/m/K 
  end type
  type(organic_par_type) :: organic
  ! mineral soil properties
  type mineral_vars
    real(wp) :: theta_sat
    real(wp) :: psi_sat
    real(wp) :: k_sat
    real(wp) :: Bi
    real(wp) :: lambda_s
    real(wp) :: lambda_dry
    real(wp) :: heatcap
  end type
  type(mineral_vars), dimension(nx,ny) :: mineral

  ! vegetation parameters
  type veg_par_type
    real(wp) :: ext_coef  = 0.5_wp  !! extinction coefficient for radiation through vegetation []
    real(wp) :: sai_scale   !! factor to scale balanced leaf area index to stem area index []
    real(wp) :: veg_h_min   !! minimum vegetation height [m]
    real(wp) :: leaf_alb  = 0.17_wp !! leaf albedo in PAR range []
    real(wp) :: gamma_down  !! photosynthesis downregulation parameter for high CO2 [] 
    real(wp) :: c_lambda_co2    !! factor for lai dependence of fraction of NPP used for spatial expansion []
    real(wp) :: lambda_exp  !! exponent of lai dependence of fraction of NPP used for spatial expansion []
    real(wp) :: lambda_min  !! minimum fraction of NPP used for spatial expansion []
    real(wp) :: lambda_max  !! maximum fraction of NPP used for spatial expansion []
    real(wp) :: delta_lai_conv  !! critical LAI difference for convergence in iterative solver
    logical :: lroot_frac   !! adjust vertical root distribution so that all roots in active layer?
    real(wp) :: gamma_phen  !! leaf senescenece rate [1/s]
    real(wp) :: gamma_luc   !! disturbance rate for land use change [1/s]
    real(wp) :: gamma_ice   !! disturbance rate from ice sheets [1/s]
    real(wp) :: z_veg_std_crit !! critical sub-grid topography standard deviation for disturbance rate from ice sheets [m]
    real(wp) :: theta_fire_crit !! critical soil moisture for fire, Thonicke 2001 ~0.4 for relative soil moisture [m3/m3]
    real(wp) :: cveg_fire_low  = 0.2_wp   !! value of aboveground biomass below which no fire [kgC/m2]
    real(wp) :: cveg_fire_high = 1._wp    !! value of aboveground biomass above which no fire limitation by fuel [kgC/m2]
    integer :: iseed    !! seeds parameterisation flag
    real(wp) :: seed_pft_min    !! minimum pft fraction needed to seed neighbouring cells
    real(wp) :: seed_fraction = 0.001_wp  !! PFT seed fraction []
    real(wp) :: seed_lim   = 0.05_wp   !! 
    real(wp) :: f_veg_crit = 0.2_wp    !! critical vegetation fraction for buried litter []
    real(wp) :: f_lit_to_ice    !! fraction of vegetation buried when ice sheets are expanding []
    integer :: i_deforest
  end type
  type(veg_par_type) :: veg_par

  ! pft specific parameters
  type pft_par_type                                                   !  BL    NL    C3    C4    SH
    ! root distribution parameters
    real(wp), dimension(npft) :: root_jules            = (/ 3._wp , 1._wp , 0.5_wp , 0.5_wp , 0.5_wp /)  ! JULES, Clark 2011
    real(wp), dimension(npft) :: root_jackson          = (/ 0.965_wp , 0.95_wp , 0.95_wp , 0.95_wp , 0.97_wp /) ! beta factor in root distribution from Jackson 1996
    real(wp), dimension(npft) :: root_clm1             = (/ 6.5_wp, 7._wp, 11._wp, 11._wp, 7._wp /) ! CLM root parameters from Zeng 2001
    real(wp), dimension(npft) :: root_clm2             = (/ 1.5_wp, 2._wp, 2._wp , 2._wp , 1.5_wp/) ! CLM root parameters from Zeng 2001
    real(wp), dimension(nl,npft) :: root_frac
    ! aboveground litter input vertical distribution
    real(wp), dimension(nl) :: litter_in_frac
    ! PFT albedo values from Houldcroft 2009, derived from MODIS for JULES, modified
    real(wp), dimension(npft) :: alb_can_vis_dir = (/ 0.04_wp, 0.04_wp, 0.04_wp, 0.03_wp, 0.04_wp/) 
    real(wp), dimension(npft) :: alb_can_vis_dif = (/ 0.05_wp, 0.05_wp, 0.05_wp, 0.04_wp, 0.05_wp/)
    real(wp), dimension(npft) :: alb_can_nir_dir = (/ 0.26_wp, 0.23_wp, 0.27_wp, 0.24_wp, 0.23_wp/)
    real(wp), dimension(npft) :: alb_can_nir_dif = (/ 0.28_wp, 0.25_wp, 0.31_wp, 0.28_wp, 0.25_wp/)
    ! snow-covered PFT albedo values for diffuse radiation (used also for direct beam) from MODIS, Moody 2007, modified
    real(wp), dimension(npft) :: alb_can_vis_dir_snow = (/ 0.45_wp, 0.4_wp, 0.70_wp, 0.70_wp, 0.6_wp/)
    real(wp), dimension(npft) :: alb_can_vis_dif_snow = (/ 0.45_wp, 0.4_wp, 0.70_wp, 0.70_wp, 0.6_wp/)
    real(wp), dimension(npft) :: alb_can_nir_dir_snow = (/ 0.35_wp, 0.3_wp, 0.48_wp, 0.48_wp, 0.4_wp/)
    real(wp), dimension(npft) :: alb_can_nir_dif_snow = (/ 0.35_wp, 0.3_wp, 0.48_wp, 0.48_wp, 0.4_wp/)
    ! PFT-specific background albedo (representing litter etc...), e.g. Vamborg 2010, Stärz 2016
    real(wp), dimension(npft) :: alb_bg_vis = (/ 0.1_wp, 0.1_wp, 0.1_wp, 0.1_wp, 0.1_wp/)
    real(wp), dimension(npft) :: alb_bg_nir = (/ 0.15_wp, 0.15_wp, 0.15_wp, 0.15_wp, 0.15_wp/)
    ! phenology parameters
    real(wp), dimension(npft) :: gdd5_phen   ! K, threshold growing degree above 5 degC for decidous
    real(wp), dimension(npft) :: t_base_phen != (/ 5.0_wp , 2.0_wp , 0.0_wp , 5.0_wp , 2.0_wp /) ! degC,base temperature for phenology
    real(wp), dimension(npft) :: t_cmon_phen != (/ 5.0_wp ,-30._wp, 0.0_wp , 0.0_wp , -35._wp /) ! degC, threshold temperature for decidous
    real(wp), dimension(npft) :: ramp        = (/ 300._wp , 200._wp , 50._wp , 100._wp , 100._wp /) ! LPJ, gdd when full leaf
    ! bioclimatic limits
    real(wp), dimension(npft) :: t_cmon_min != (/ -17._wp ,-999._wp, -999._wp , 15.5_wp , -999.0_wp /) ! °C, minimum coldest month temperature, LPJ
    real(wp), dimension(npft) :: t_cmon_max != (/ 999._wp ,-5._wp, 15.5_wp , 999._wp , 999.0_wp /) ! °C, maximum coldest month temperature, LJ
    real(wp), dimension(npft) :: gdd5_min   = (/ 1200._wp , 350._wp , 0._wp , 0._wp , 0._wp /) ! minimum gdd(5°C) for establishment, LJ
    ! photosynthesis parameters
    integer, dimension(npft)          :: i_c4                  = (/ 0   , 0   , 0   , 1   , 0   /) 
    real(wp), dimension(npft) :: t_co2_low             = (/-4.0_wp , -4.0_wp, -4.0_wp,  6.0_wp, -4.0_wp/)
    real(wp), dimension(npft) :: t_co2_high            = (/55.0_wp , 42.0_wp, 45.0_wp, 55.0_wp, 40.0_wp/)
    real(wp), dimension(npft) :: t_photos_low          = (/20.0_wp , 15.0_wp, 10.0_wp, 20.0_wp, 10.0_wp/)
    real(wp), dimension(npft) :: t_photos_high         = (/30.0_wp , 30.0_wp, 30.0_wp, 45.0_wp, 30.0_wp/) !! high temperature for photosynthesis
    real(wp), dimension(npft) :: g_min = (/0.5_wp,0.3_wp, 0.5_wp,0.5_wp,0.5_wp/) !! minimum stomatal conductance [mm/s]
    real(wp), dimension(npft) :: g1    = (/4._wp,2.3_wp,3._wp,1.6_wp,4.2_wp/)  !! parameter in optimal stomatal conductance model, Lin2015 []
    !real(wp), dimension(npft) :: beta  = (/200._wp,140._wp,180._wp,100._wp,220._wp/)  !! parameter in least cost optimality model
    real(wp), dimension(npft) :: beta  = (/150._wp,150._wp,150._wp,150._wp,150._wp/)  !! parameter in least cost optimality model, Lavergne 2020 []
    real(wp), dimension(npft) :: photos_k1, photos_k2, photos_k3
    real(wp), dimension(npft) :: flnr = (/0.08,0.06,0.15,0.09,0.08/)  !!  CLM 4.5, Table 8.1 [gNinRub/gN]
    ! vegetation dynamics parameters
    real(wp), dimension(npft) :: gamma_leaf = (/0.7_wp,0.3_wp,1._wp,2._wp,0.5_wp/)/sec_year !! leaf turnover rate [1/s]
    real(wp), dimension(npft) :: gamma_root = (/0.15_wp, 0.15_wp, 0.3_wp, 0.5_wp, 0.3_wp/)/sec_year !! root turnover rate from Gill 2000 [1/s]
    real(wp), dimension(npft) :: gamma_stem = (/0.005_wp,0.005_wp,0.1_wp,0.1_wp,0.1_wp/)/sec_year !! stem turnover rate [1/s]
    real(wp), dimension(npft) :: gamma_dist_min  !! minimum disturbance rate [1/s]
    real(wp), dimension(npft) :: tau_fire   !! fire disturbance time scale [s]
    real(wp), dimension(npft) :: lai_min !! minimum leaf area index [m2/m2]
    real(wp), dimension(npft) :: lai_max !! maximum leaf area index [m2/m2]
    real(wp), dimension(npft) :: sla    !! specific leaf area, TRY 2011 [m2/(kg dry leaf) * (kg dry leaf)/kgC = m2/kgC]
    !real(wp), dimension(npft) :: awl = (/2._wp,2._wp,0.01_wp,0.01_wp,0.5_wp/) !! allometric coefficient, modified from Cox 2001 to account for bwl=1 [kgC/m2,]
    real(wp), dimension(npft) :: awl = (/0.65_wp,0.65_wp,0.005_wp,0.005_wp,0.3_wp/) !! allometric coefficient, Cox 2001 for bwl=5/3 [kgC/m2,]
    real(wp) :: bwl = 5._wp/3._wp !! allometric coefficient, Cox 2001, bwl=5/3 
    real(wp), dimension(npft) :: aws = (/10._wp,10._wp,1._wp,1._wp,5._wp/) !! ratio of total to respiring stem carbon, Cox 2001 []
    real(wp), dimension(npft) :: awh = (/3.5_wp,6._wp,0.15_wp,0.17_wp,1._wp/) !! factor relating plant height to balanced LAI, derived from fitting TRY plant heights to LAI []
    real(wp), dimension(npft) :: hveg_z0_scale

  end type
  type(pft_par_type) :: pft_par

  ! peatland parameters
  type peat_pars
    logical :: peat_carb
    logical :: peat_area
    integer :: nmonwet_peat
    integer :: nyearwet_peat = 30
    real(wp) :: f_peat_min  = 0.001_wp !! minimum peatland area fraction []
    real(wp) :: peat_ch_rate   !! rate of change of peat fraction [1/s]
    real(wp) :: dCpeat_dt_min  !! minimum peat change rate for peat survival [kgC/m2/s]
    real(wp) :: k10_acro   !! !! k10 for acrotelm decomposition [1/s]
    real(wp) :: k10_cato   !! k10 for catotelm decomposition [1/s]
    real(wp) :: acroc_crit    !! critical acrotelm carbon content for catotelm formation [kgC/m2]
    real(wp) :: k_acro_to_cato    !! acrotelm to catotelm transfer rate [1/s] 
    real(wp) :: fmoist_peat = 0.2_wp   !! moisture limitation factor for anoxic acrotelm soil decomposition, Wania 2009, Koven 2013 []
    real(wp) :: Cpeat_min  !! minimum peat carbon for peat survival, Stocker 2014 [kgC/m2]
    real(wp) :: rho_acro    = 20._wp !! acrotelm carbon density, Wania 2009, Kleinen 2012, Spahni 2013 [kgC/m3]
    real(wp) :: rho_cato    = 50._wp !! catotelm carbon density, Kleinen 2012, Lawrence 2008 [kgC/m3]
  end type
  type(peat_pars) :: peat_par
  integer, parameter :: nmonwet = nmon_year*30    !! 

  ! soil carbon parameters
  type soilc_par_type
    logical :: l_burial !! allow for carbon burial?
    real(wp) :: f_resp_litter  !! fraction of decomposed litter going to atmosphere
    real(wp) :: f_litter_to_fast  !! fraction of decomposed litter carbon going to fast pool
    real(wp) :: f_litter_to_slow  !! fraction of decomposed litter carbon going to slow pool
    real(wp) :: z_litter_in   !! e-folding depth for vertical profile of litter input
    real(wp) :: k_min  !! minimum carbon turnover rate [1/s]
    real(wp) :: k10_litter  !! litter carbon turnover rate at 10 degC [1/s]
    real(wp) :: k10_fast    !! fast soil carbon turnover rate at 10 degC [1/s] 
    real(wp) :: k10_slow    !! slow soil carbon turnover rate at 10 degC [1/s]
    real(wp) :: z_tau   !! e-folding depth for soil carbon decompositions reduction
    real(wp) :: heat_decomp_spec = 40.d6  !! specific heat of soil carbon decomposition [J/kgC], 40 MJ/kgC in Khvorostyanov 2008, range 25-40
    real(wp) :: n_alt
    real(wp) :: k_ice   !! turnover rate of carbon below ice sheets [1/s]
    logical :: heat_soc
    integer :: iresp_temp
    integer :: iresp_moist
    real(wp) :: theta_crit_fmoist
    real(wp) :: psi_min, psi_max
    real(wp) :: Ea
    real(wp) :: q10_c
    real(wp) :: diff_bio, diff_cryo, diff_shelf, diff_ice, diff_lake
    real(wp) :: diff_min = 0._wp/sec_year   !! minimum soil carbon vertical diffusivity (for stability) [m2/s]
    real(wp) :: z_diff   !! e-folding depth for diffusivity reduction with depth
    integer :: iadv_soilc  !! parameterisation of vertical advection of soil carbon []
    real(wp) :: rho_loess = 1400._wp  !! loess bulk density [kg/m3]
    real(wp) :: adv_soil    !! advection soil carbon [m/s]
    real(wp) :: adv_shelf   !! advection shelf carbon [m/s] 
    real(wp) :: adv_ice     !! advection ice carbon [m/s]
    real(wp) :: rho_soc_max = 50._wp !! maximum soil organic carbon density, be consistent with rho_cato [kgC/m3]
  end type
  type(soilc_par_type) :: soilc_par

  ! methane parameters
  type ch4_par_type
    integer :: ich4_ftemp
    real(wp) :: q10_ch4
    real(wp) :: Ea_ch4
    real(wp) :: ch4_frac_wet
    real(wp) :: ch4_frac_peat
    real(wp) :: ch4_frac_shelf
    real(wp) :: ch4_frac_lake
    real(wp) :: c_ch4_conv = 16._wp/12._wp ! conversion from C to CH4, mass ratio of CH4/C
  end type
  type(ch4_par_type) :: ch4_par

  ! weathering parameters
  integer, parameter :: nlit_gemco2 = 6
  type weath_gemco2_type
    integer :: nlit = nlit_gemco2
    ! GEM-CO2 model, values from Table 2 of Colbourn 2013
    real(wp), dimension(nlit_gemco2) :: frac_carb= [0.93_wp,    0.39_wp,  0.48_wp,  0._wp,    0._wp,    0._wp]   !! fraction to weather as carbonate rocks [1]
    real(wp), dimension(nlit_gemco2) :: kweath   = [1.586_wp,   0.626_wp, 0.152_wp, 0.479_wp, 0.095_wp, 0.222_wp]*1.e-3_wp    !! weathering rate [mol C / kg water]
    real(wp), dimension(nx,ny,nlit_gemco2) :: lit_map
  end type weath_gemco2_type
  type(weath_gemco2_type) :: weath_gemco2_par

  integer, parameter :: nlit_uhh = 14
  type weath_uhh_type
    integer :: nlit = nlit_uhh
    integer :: i_loess = 1
    integer :: i_ev = 2
    integer :: i_mt = 3
    integer :: i_pa = 4
    integer :: i_pb = 5
    integer :: i_pi = 6
    integer :: i_py = 7
    integer :: i_sc = 8
    integer :: i_sm = 9
    integer :: i_ss = 10
    integer :: i_su = 11
    integer :: i_va = 12
    integer :: i_vb = 13
    integer :: i_vi = 14
    real(wp), dimension(nlit_uhh) :: b         = 1._wp/12._wp * &      !! weathering rate [mol C / kg water]
      [0._wp,0._wp,0.007626_wp,0.005095_wp,0.007015_wp,0.007015_wp,0.0061_wp,0._wp,0.012481_wp,0.005341_wp,0.003364_wp,0.002455_wp,0.007015_wp,0.007015_wp]
    real(wp), dimension(nlit_uhh) :: frac_carb = &      !! fraction to weather as carbonate rocks [1]
      [1._wp,0._wp,0.75_wp,    0.42_wp,    0._wp,      0.42_wp,    0._wp,    1._wp,0.76_wp,    0.36_wp,    0._wp,      0._wp,      0._wp,      0._wp]
    real(wp), dimension(nlit_uhh) :: sa        = &      !! activation energy silicates
      [0._wp,0._wp,60._wp,     60._wp,     50._wp,     60._wp,     46._wp,   0._wp,60._wp,     60._wp,     60._wp,     60._wp,     50._wp,     50._wp]
    real(wp) :: ca = 14._wp     !! activation energy carbonates
    real(wp), dimension(nx,ny,nlit_uhh) :: lit_map
    real(wp), dimension(nx,ny,nlit_uhh) :: lit_map_lgm
    real(wp), dimension(nx,ny,nlit_uhh) :: lit_map_shelf
  end type weath_uhh_type
  type(weath_uhh_type) :: weath_uhh_par

contains

  subroutine lnd_params_init

  implicit none

  integer :: i, j, n, k
  real(wp), dimension(nx,ny) :: frac_sand, frac_clay
  real(wp), dimension(nx,ny) :: bulk_den


  call lnd_par_load

  ! compute all timesteps based on namelist input
  dt = dt_lnd
  dt_day = dt/sec_day
  rdt = 1._wp/dt
  dt_day_c = dt_day_carb
  dt_day_v = dt_day_veg
  dt_c = real(dt_day_carb,wp)*sec_day
  dt_v = real(dt_day_veg,wp)*sec_day
  nstep_v = nint(dt_day_v/dt_day)

  ! define standard vertical root distribution in the absence of permafrost
  do n=1,npft
    do k=1,nl
      if( i_roots .eq. 1 ) then
        ! JULES
        pft_par%root_frac(k,n) = ( exp(-2._wp*z_int(k-1)/pft_par%root_jules(n)) - exp(-2._wp*z_int(k)/pft_par%root_jules(n)) ) &
        / ( 1._wp - exp(-2._wp*z_int(nl)/pft_par%root_jules(n)) )
      else if( i_roots .eq. 2 ) then
        ! Jackson 1996
        pft_par%root_frac(k,n) = (pft_par%root_jackson(n)**(z_int(k-1)*100._wp) - pft_par%root_jackson(n)**(z_int(k)*100._wp)) &
        / (1._wp - pft_par%root_jackson(n)**(z_int(nl)*100._wp))
      else if( i_roots .eq. 3 ) then
        ! CLM, Zeng 2001
        if( k .lt. nl ) then
          pft_par%root_frac(k,n) = 0.5_wp * (exp(-pft_par%root_clm1(n)*z_int(k-1)) + exp(-pft_par%root_clm2(n)*z_int(k-1)) &
          - exp(-pft_par%root_clm1(n)*z_int(k)) - exp(-pft_par%root_clm2(n)*z_int(k)))
        else
          pft_par%root_frac(k,n) = 0.5_wp * (exp(-pft_par%root_clm1(n)*z_int(k)) + exp(-pft_par%root_clm2(n)*z_int(k)))
        endif
      endif
    enddo
    pft_par%root_frac(:,n) = pft_par%root_frac(:,n) / sum(pft_par%root_frac(:,n))
     !print *
     !print *,n,pft_par%root_frac(:,n)
     !print *,sum(root_frac(:,n))
  enddo

  ! define standard vertical litter input distribution in the absence of permafrost
  do k=1,nl
    pft_par%litter_in_frac(k) = ( exp(-z_int(k-1)/soilc_par%z_litter_in) - exp(-z_int(k)/soilc_par%z_litter_in) ) &
      / ( 1._wp - exp(-z_int(nl)/soilc_par%z_litter_in) )
  enddo
  pft_par%litter_in_frac = pft_par%litter_in_frac / sum(pft_par%litter_in_frac)
  !print *
  !print *,pft_par%litter_in_frac
  !print *,sum(pft_par%litter_in_frac)

  ! photosynthesis temperature stress factors
  pft_par%photos_k1 = 2._wp * log(1._wp/0.99_wp-1._wp) / (pft_par%t_co2_low - pft_par%t_photos_low)
  pft_par%photos_k2 = 0.5_wp * (pft_par%t_co2_low + pft_par%t_photos_low)
  pft_par%photos_k3 = log(0.99_wp/0.01_wp) / (pft_par%t_co2_high - pft_par%t_photos_high)

  ! bare soil/desert albedo
  call nc_read( trim(lnd_par_file), "a_vis",surf_par%alb_bare_vis )
  call nc_read( trim(lnd_par_file), "a_nir",surf_par%alb_bare_nir )
  ! set to constant value in mid- to high latitudes (poleward of 50°)
  do i=1,nx
    do j=1,ny
      if (abs(lat(j)).gt.50._wp) then
        surf_par%alb_bare_vis(i,j) = 0.2_wp
        surf_par%alb_bare_nir(i,j) = 0.2_wp
      endif
    enddo
  enddo

  ! soil texture
  if( soil_par%soil_texture .eq. 1 ) then
    call nc_read( trim(lnd_par_file), "fclay_hwsd",frac_clay )
    call nc_read( trim(lnd_par_file), "fsand_hwsd",frac_sand )
    elseif( soil_par%soil_texture .eq. 2 ) then
    call nc_read( trim(lnd_par_file), "fclay_gsde",frac_clay )
    call nc_read( trim(lnd_par_file), "fsand_gsde",frac_sand )
  endif
  
  ! mineral soil characteristics
  mineral%theta_sat = 0.489_wp-0.00126_wp*frac_sand
  mineral%lambda_s = (8.80_wp*frac_sand+2.92_wp*frac_clay)/(frac_sand+frac_clay)
  bulk_den = 2700.0_wp*(1.0_wp-mineral%theta_sat)
  mineral%lambda_dry = (0.135_wp*bulk_den+64.7_wp)/(2700.0_wp-0.947_wp*bulk_den)
  mineral%heatcap = (2.128_wp*frac_sand+2.385_wp*frac_clay)/(frac_sand+frac_clay)*1.0d6
  mineral%psi_sat = (-10.0_wp*10**(1.88_wp-0.0131_wp*frac_sand)) *1.d-3  ! m
  mineral%k_sat = 0.0070556_wp*10**(-0.884_wp+0.0153_wp*frac_sand)  ! kg/m2/s
  mineral%Bi = 2.91_wp+0.159_wp*frac_clay

  ! read topmodel parameters
  call nc_read('input/cti_5x5.nc', "cti_mean",topmodel%cti_mean)
  do i=1,15
    call nc_read('input/cti_5x5.nc', "cti_cdf",topmodel%cti_cdf(i),start=[1,1,i],count=[nx,ny,1])
  enddo

  !call nc_read( "input/fmax_5x5.nc", "fmax",clm_fmax )
  ! read dyptop parameters
  call nc_read( lnd_par_file, "k",dyptop%k )
  call nc_read( lnd_par_file, "xm",dyptop%xm )
  call nc_read( lnd_par_file, "v",dyptop%v )
  call nc_read( lnd_par_file, "fmax",dyptop%fmax )

  ! lithologies

  if (i_weathering.eq.1) then

    call nc_read( lithology_gemco2_file, "lit_map",weath_gemco2_par%lit_map )
    do i=1,nx
      do j=1,ny
        if (sum(weath_gemco2_par%lit_map(i,j,:)).eq.0._wp) then
          ! equal share of all lithologies if map not available
          weath_gemco2_par%lit_map(i,j,:) = 1._wp/real(weath_gemco2_par%nlit,wp)
        endif
      enddo
    enddo

  endif
 
  if (i_weathering.eq.2) then

    call nc_read( lithology_uhh_file, "loess",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_loess) )
    call nc_read( lithology_uhh_file, "ev",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_ev) )
    call nc_read( lithology_uhh_file, "mt",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_mt) )
    call nc_read( lithology_uhh_file, "pa",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_pa) )
    call nc_read( lithology_uhh_file, "pb",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_pb) )
    call nc_read( lithology_uhh_file, "pi",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_pi) )
    call nc_read( lithology_uhh_file, "py",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_py) )
    call nc_read( lithology_uhh_file, "sc",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_sc) )
    call nc_read( lithology_uhh_file, "sm",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_sm) )
    call nc_read( lithology_uhh_file, "ss",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_ss) )
    call nc_read( lithology_uhh_file, "su",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_su) )
    call nc_read( lithology_uhh_file, "va",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_va) )
    call nc_read( lithology_uhh_file, "vb",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_vb) )
    call nc_read( lithology_uhh_file, "vi",weath_uhh_par%lit_map(:,:,weath_uhh_par%i_vi) )
    do i=1,nx
      do j=1,ny
        if (sum(weath_uhh_par%lit_map(i,j,:)).eq.0._wp) then
          ! equal share of all lithologies if map not available
          weath_uhh_par%lit_map(i,j,:) = 1._wp/real(weath_uhh_par%nlit,wp)
        endif
      enddo
    enddo

    ! lgm
    call nc_read( "input/Lithology_lgm_UHH.nc", "loess",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_loess) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "ev",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_ev) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "mt",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_mt) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "pa",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_pa) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "pb",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_pb) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "pi",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_pi) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "py",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_py) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "sc",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_sc) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "sm",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_sm) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "ss",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_ss) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "su",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_su) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "va",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_va) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "vb",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_vb) )
    call nc_read( "input/Lithology_lgm_UHH.nc", "vi",weath_uhh_par%lit_map_lgm(:,:,weath_uhh_par%i_vi) )
    do i=1,nx
      do j=1,ny
        if (sum(weath_uhh_par%lit_map_lgm(i,j,:)).eq.0._wp) then
          ! equal share of all lithologies if map not available
          weath_uhh_par%lit_map_lgm(i,j,:) = 1._wp/real(weath_uhh_par%nlit,wp)
        endif
      enddo
    enddo

    call nc_read( 'input/Lithology_shelves_UHH.nc', "loess",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_loess) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "ev",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_ev) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "mt",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_mt) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "pa",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_pa) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "pb",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_pb) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "pi",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_pi) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "py",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_py) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "sc",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_sc) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "sm",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_sm) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "ss",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_ss) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "su",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_su) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "va",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_va) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "vb",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_vb) )
    call nc_read( 'input/Lithology_shelves_UHH.nc', "vi",weath_uhh_par%lit_map_shelf(:,:,weath_uhh_par%i_vi) )
    do i=1,nx
      do j=1,ny
        if (sum(weath_uhh_par%lit_map_shelf(i,j,:)).eq.0._wp) then
          ! equal share of all lithologies if map not available
          weath_uhh_par%lit_map_shelf(i,j,:) = 1._wp/real(weath_uhh_par%nlit,wp)
        endif
      enddo
    enddo

  endif


  return

  end subroutine lnd_params_init


subroutine lnd_par_load

    use nml

    implicit none 

    integer :: k, n
    real(wp) :: sla_nl, sla_bl, sla_c3, sla_c4, sla_sh
    real(wp) :: gamma_dist_tree, gamma_dist_grass, gamma_dist_shrub
    real(wp) :: tau_fire
    real(wp) :: alb_ice
    real(wp) :: alpha_int_w_tree, alpha_int_w_grass

    character (len=256) :: filename


    filename = trim(out_dir)//"/lnd_par.nml"

    ! read parameter namelist
    call nml_read(filename,"lnd_par","i_init_veg",i_init_veg)
    call nml_read(filename,"lnd_par","l_dynveg",l_dynveg)
    call nml_read(filename,"lnd_par","pft_fix_file",pft_fix_file)
    call nml_read(filename,"lnd_par","l_fixlai",l_fixlai)
    call nml_read(filename,"lnd_par","lai_fix_file",lai_fix_file)
    call nml_read(filename,"lnd_par","l_co2_fert_lim",l_co2_fert_lim)
    call nml_read(filename,"lnd_par","co2_fert_lim_min",co2_fert_lim_min)
    call nml_read(filename,"lnd_par","co2_fert_lim_max",co2_fert_lim_max)
    call nml_read(filename,"lnd_par","dt_day_veg",dt_day_veg)
    call nml_read(filename,"lnd_par","dt_day_carb",dt_day_carb)
    call nml_read(filename,"lnd_par","lnd_par_file",lnd_par_file)

    call nml_read(filename,"lnd_par","l_diurnal_cycle" ,l_diurnal_cycle )
    call nml_read(filename,"lnd_par","l_neutral",l_neutral)
    call nml_read(filename,"lnd_par","z_sfl",z_sfl)
    call nml_read(filename,"lnd_par","i_racan",i_racan)
    call nml_read(filename,"lnd_par","p_cdense",p_cdense)
    call nml_read(filename,"lnd_par","i_roots",i_roots)
    call nml_read(filename,"lnd_par","z0m_bare",surf_par%z0m_bare)
    call nml_read(filename,"lnd_par","z0m_ice",surf_par%z0m_ice)
    call nml_read(filename,"lnd_par","z0m_lake",surf_par%z0m_lake)
    call nml_read(filename,"lnd_par","z0m_lake_ice",surf_par%z0m_lake_ice)
    call nml_read(filename,"lnd_par","i_z0h",surf_par%i_z0h)
    call nml_read(filename,"lnd_par","zm_to_zh_const",surf_par%zm_to_zh_const)
    call nml_read(filename,"lnd_par","f_Ri_unstab",surf_par%f_Ri_unstab)
    call nml_read(filename,"lnd_par","f_Ri_stab",surf_par%f_Ri_stab)

    call nml_read(filename,"lnd_par","K_eddy_lake_bg",K_eddy_lake_bg)
    call nml_read(filename,"lnd_par","K_eddy_lake_max",K_eddy_lake_max)
    call nml_read(filename,"lnd_par","z_mix_lake_min",z_mix_lake_min)
    call nml_read(filename,"lnd_par","dz_lake_ice",dz_lake_ice)

    call nml_read(filename,"lnd_par","i_snow_albedo",snow_par%i_snow_albedo)
    call nml_read(filename,"lnd_par","f_snow_sub",snow_par%f_snow_sub)
    call nml_read(filename,"lnd_par","l_snow_aging",snow_par%l_snow_aging)
    call nml_read(filename,"lnd_par","l_snow_dust",snow_par%l_snow_dust)
    call nml_read(filename,"lnd_par","snow_grain_fresh" ,snow_par%snow_grain_fresh )
    call nml_read(filename,"lnd_par","snow_grain_old" ,snow_par%snow_grain_old )
    call nml_read(filename,"lnd_par","f_age_t" ,snow_par%f_age_t)
    call nml_read(filename,"lnd_par","snow_0" ,snow_par%snow_0 )
    snow_par%snow_0 = snow_par%snow_0/sec_day
    call nml_read(filename,"lnd_par","snow_1" ,snow_par%snow_1 )
    call nml_read(filename,"lnd_par","rho_snow",snow_par%rho_snow)
    call nml_read(filename,"lnd_par","lambda_snow",snow_par%lambda_snow)
    call nml_read(filename,"lnd_par","w_snow_crit",snow_par%w_snow_crit)
    call nml_read(filename,"lnd_par","w_snow_off",snow_par%w_snow_off)
    call nml_read(filename,"lnd_par","l_fsnow_orog" ,snow_par%l_fsnow_orog)
    call nml_read(filename,"lnd_par","c_fsnow" ,snow_par%c_fsnow)
    call nml_read(filename,"lnd_par","c_fsnow_orog" ,snow_par%c_fsnow_orog)
    call nml_read(filename,"lnd_par","z_veg_std_min" ,snow_par%z_veg_std_min)

    call nml_read(filename,"lnd_par","l_ice_albedo_semi",l_ice_albedo_semi)
    call nml_read(filename,"lnd_par","alb_ice",alb_ice)
    surf_par%alb_vis_dir_ice = alb_ice 
    surf_par%alb_vis_dif_ice = alb_ice
    surf_par%alb_nir_dir_ice = alb_ice
    surf_par%alb_nir_dif_ice = alb_ice

    call nml_read(filename,"lnd_par","i_runoff",hydro_par%i_runoff)
    call nml_read(filename,"lnd_par","l_runoff_icemelt",hydro_par%l_runoff_icemelt)
    call nml_read(filename,"lnd_par","l_dew",hydro_par%l_dew)
    call nml_read(filename,"lnd_par","l_prc_intercept",hydro_par%l_prc_intercept)
    call nml_read(filename,"lnd_par","i_cond_theta",hydro_par%i_cond_theta)
    call nml_read(filename,"lnd_par","i_wtab",hydro_par%i_wtab)
    call nml_read(filename,"lnd_par","i_fwet",hydro_par%i_fwet)
    call nml_read(filename,"lnd_par","i_evp_soil",hydro_par%i_evp_soil)
    call nml_read(filename,"lnd_par","theta_crit_evp",hydro_par%theta_crit_evp)
    call nml_read(filename,"lnd_par","dz_evp",hydro_par%dz_evp)
    call nml_read(filename,"lnd_par","frac_inf_2",hydro_par%frac_inf_2)
    call nml_read(filename,"lnd_par","p_psi_min",hydro_par%p_psi_min)
    call nml_read(filename,"lnd_par","p_psi_max",hydro_par%p_psi_max)
    call nml_read(filename,"lnd_par","wtab_scale",hydro_par%wtab_scale)
    call nml_read(filename,"lnd_par","f_wtab",hydro_par%f_wtab)
    call nml_read(filename,"lnd_par","cti_min",hydro_par%cti_min)
    call nml_read(filename,"lnd_par","fmax_crit",hydro_par%fmax_crit)
    call nml_read(filename,"lnd_par","cti_mean_crit",hydro_par%cti_mean_crit)
    call nml_read(filename,"lnd_par","alpha_int_w_tree",alpha_int_w_tree)
    call nml_read(filename,"lnd_par","alpha_int_w_grass",alpha_int_w_grass)
    do n=1,npft
      if (flag_tree(n) .or. flag_shrub(n)) then
        hydro_par%alpha_int_w(n) = alpha_int_w_tree
      else
        hydro_par%alpha_int_w(n) = alpha_int_w_grass
      endif
    enddo
    call nml_read(filename,"lnd_par","alpha_int_s",hydro_par%alpha_int_s)
    call nml_read(filename,"lnd_par","can_max_w",hydro_par%can_max_w  )
    call nml_read(filename,"lnd_par","can_max_s",hydro_par%can_max_s  )
    call nml_read(filename,"lnd_par","tau_w",hydro_par%tau_w)
    call nml_read(filename,"lnd_par","tau_s",hydro_par%tau_s)
    hydro_par%tau_w = hydro_par%tau_w*sec_day
    hydro_par%tau_s = hydro_par%tau_s*sec_day

    call nml_read(filename,"lnd_par","b0",dust_par%b0)
    call nml_read(filename,"lnd_par","qd",dust_par%qd)
    call nml_read(filename,"lnd_par","qg",dust_par%qg)
    call nml_read(filename,"lnd_par","u0",dust_par%u0)
    call nml_read(filename,"lnd_par","wind_gust_fac",dust_par%wind_gust_fac)
    call nml_read(filename,"lnd_par","i_fsnow",dust_par%i_fsnow)
    call nml_read(filename,"lnd_par","i_theta",dust_par%i_theta)
    call nml_read(filename,"lnd_par","sm_t",dust_par%sm_t)
    call nml_read(filename,"lnd_par","sm_n",dust_par%sm_n)
    call nml_read(filename,"lnd_par","lai_t",dust_par%lai_t)
    call nml_read(filename,"lnd_par","lai_n",dust_par%lai_n)
    call nml_read(filename,"lnd_par","l_dust_stab",dust_par%l_dust_stab)
    call nml_read(filename,"lnd_par","l_dust_topo",dust_par%l_dust_topo)
    call nml_read(filename,"lnd_par","topo_exp",dust_par%topo_exp)

    call nml_read(filename,"lnd_par","sw_par_frac",sw_par_frac)

    call nml_read(filename,"lnd_par","i_ci",i_ci)
    call nml_read(filename,"lnd_par","i_beta",i_beta)
    call nml_read(filename,"lnd_par","i_vcmax",i_vcmax)
    call nml_read(filename,"lnd_par","i_disc",i_disc)

    call nml_read(filename,"lnd_par","sla_bl",sla_bl)
    call nml_read(filename,"lnd_par","sla_nl",sla_nl)
    call nml_read(filename,"lnd_par","sla_c3",sla_c3)
    call nml_read(filename,"lnd_par","sla_c4",sla_c4)
    call nml_read(filename,"lnd_par","sla_sh",sla_sh)
    ! convert from m2/(kg dry leaf) to m2/kgC
    pft_par%sla(1) = sla_bl*2._wp
    pft_par%sla(2) = sla_nl*2._wp
    pft_par%sla(3) = sla_c3*2._wp
    pft_par%sla(4) = sla_c4*2._wp
    pft_par%sla(5) = sla_sh*2._wp
    call nml_read(filename,"lnd_par","gamma_dist_tree",gamma_dist_tree)
    call nml_read(filename,"lnd_par","gamma_dist_grass",gamma_dist_grass)
    call nml_read(filename,"lnd_par","gamma_dist_shrub",gamma_dist_shrub)
    ! convert from 1/yr to 1/s
    pft_par%gamma_dist_min(1) = gamma_dist_tree/sec_year 
    pft_par%gamma_dist_min(2) = gamma_dist_tree/sec_year
    pft_par%gamma_dist_min(3) = gamma_dist_grass/sec_year
    pft_par%gamma_dist_min(4) = gamma_dist_grass/sec_year
    pft_par%gamma_dist_min(5) = gamma_dist_shrub/sec_year
    call nml_read(filename,"lnd_par","gamma_luc",veg_par%gamma_luc)
    ! convert from 1/yr to 1/s
    veg_par%gamma_luc = veg_par%gamma_luc/sec_year
    call nml_read(filename,"lnd_par","gamma_ice",veg_par%gamma_ice)
    ! convert from 1/yr to 1/s
    veg_par%gamma_ice = veg_par%gamma_ice/sec_year
    call nml_read(filename,"lnd_par","z_veg_std_crit",veg_par%z_veg_std_crit)

    call nml_read(filename,"lnd_par","iseed",veg_par%iseed)
    call nml_read(filename,"lnd_par","seed_pft_min",veg_par%seed_pft_min)
    call nml_read(filename,"lnd_par","lai_min",pft_par%lai_min)
    call nml_read(filename,"lnd_par","lai_max",pft_par%lai_max)
    call nml_read(filename,"lnd_par","c_lambda_co2",veg_par%c_lambda_co2)
    call nml_read(filename,"lnd_par","lambda_exp",veg_par%lambda_exp)
    call nml_read(filename,"lnd_par","lambda_min",veg_par%lambda_min)
    call nml_read(filename,"lnd_par","lambda_max",veg_par%lambda_max)
    call nml_read(filename,"lnd_par","delta_lai_conv",veg_par%delta_lai_conv)
    call nml_read(filename,"lnd_par","gamma_down",veg_par%gamma_down)
    call nml_read(filename,"lnd_par","lroot_frac",veg_par%lroot_frac)
    call nml_read(filename,"lnd_par","gamma_phen",veg_par%gamma_phen)
    veg_par%gamma_phen = 1._wp/veg_par%gamma_phen / sec_day
    call nml_read(filename,"lnd_par","gdd5_phen",pft_par%gdd5_phen)
    call nml_read(filename,"lnd_par","t_base_phen",pft_par%t_base_phen)
    call nml_read(filename,"lnd_par","t_cmon_phen",pft_par%t_cmon_phen)
    call nml_read(filename,"lnd_par","t_cmon_min",pft_par%t_cmon_min)
    call nml_read(filename,"lnd_par","t_cmon_max",pft_par%t_cmon_max)
    call nml_read(filename,"lnd_par","hveg_z0_scale",pft_par%hveg_z0_scale)
    call nml_read(filename,"lnd_par","sai_scale",veg_par%sai_scale)
    call nml_read(filename,"lnd_par","veg_h_min",veg_par%veg_h_min)
    call nml_read(filename,"lnd_par","f_lit_to_ice",veg_par%f_lit_to_ice)
    call nml_read(filename,"lnd_par","tau_fire",tau_fire)
    pft_par%tau_fire = tau_fire*sec_year
    call nml_read(filename,"lnd_par","theta_fire_crit",veg_par%theta_fire_crit)
    call nml_read(filename,"lnd_par","i_deforest",veg_par%i_deforest)

    call nml_read(filename,"lnd_par","soil_texture",soil_par%soil_texture)
    call nml_read(filename,"lnd_par","uniform_porosity",soil_par%uniform_porosity)
    call nml_read(filename,"lnd_par","uniform_soil_par_hydro",soil_par%uniform_soil_par_hydro)
    call nml_read(filename,"lnd_par","uniform_soil_par_therm",soil_par%uniform_soil_par_therm)
    call nml_read(filename,"lnd_par","constant_porosity",soil_par%constant_porosity)
    call nml_read(filename,"lnd_par","constant_soil_par_hydro",soil_par%constant_soil_par_hydro)
    call nml_read(filename,"lnd_par","constant_soil_par_therm",soil_par%constant_soil_par_therm)
    call nml_read(filename,"lnd_par","k_sat_u",soil_par%k_sat_u)
    soil_par%k_sat_u = soil_par%k_sat_u/sec_day
    call nml_read(filename,"lnd_par","psi_sat_u",soil_par%psi_sat_u)
    call nml_read(filename,"lnd_par","b_u",soil_par%b_u)
    call nml_read(filename,"lnd_par","theta_sat_u",soil_par%theta_sat_u)
    call nml_read(filename,"lnd_par","theta_field_u",soil_par%theta_field_u)
    call nml_read(filename,"lnd_par","theta_wilt_u",soil_par%theta_wilt_u)

    call nml_read(filename,"lnd_par","lambda_s_u",soil_par%lambda_s_u)
    call nml_read(filename,"lnd_par","lambda_dry_u",soil_par%lambda_dry_u)

    call nml_read(filename,"lnd_par","l_burial",soilc_par%l_burial)
    call nml_read(filename,"lnd_par","heat_soc",soilc_par%heat_soc)
    call nml_read(filename,"lnd_par","f_resp_litter",soilc_par%f_resp_litter)
    call nml_read(filename,"lnd_par","f_litter_to_fast",soilc_par%f_litter_to_fast)
    soilc_par%f_litter_to_slow = 1._wp - soilc_par%f_litter_to_fast
    call nml_read(filename,"lnd_par","k_min",soilc_par%k_min)
    call nml_read(filename,"lnd_par","k10_litter",soilc_par%k10_litter)
    call nml_read(filename,"lnd_par","k10_fast",soilc_par%k10_fast)
    call nml_read(filename,"lnd_par","k10_slow",soilc_par%k10_slow)
    soilc_par%k_min      = 1._wp/soilc_par%k_min/sec_year   ! 1/s
    soilc_par%k10_litter = 1._wp/soilc_par%k10_litter/sec_year   ! 1/s
    soilc_par%k10_fast   = 1._wp/soilc_par%k10_fast/sec_year   ! 1/s
    soilc_par%k10_slow   = 1._wp/soilc_par%k10_slow/sec_year  ! 1/s
    call nml_read(filename,"lnd_par","z_tau",soilc_par%z_tau)
    call nml_read(filename,"lnd_par","z_litter_in",soilc_par%z_litter_in)
    call nml_read(filename,"lnd_par","iresp_temp",soilc_par%iresp_temp)
    call nml_read(filename,"lnd_par","iresp_moist",soilc_par%iresp_moist)
    call nml_read(filename,"lnd_par","theta_crit_fmoist",soilc_par%theta_crit_fmoist)
    call nml_read(filename,"lnd_par","psi_min",soilc_par%psi_min)
    call nml_read(filename,"lnd_par","psi_max",soilc_par%psi_max)
    call nml_read(filename,"lnd_par","Ea",soilc_par%Ea)
    call nml_read(filename,"lnd_par","q10_c",soilc_par%q10_c)
    call nml_read(filename,"lnd_par","n_alt",soilc_par%n_alt)
    call nml_read(filename,"lnd_par","k_ice",soilc_par%k_ice)
    soilc_par%k_ice = 1._wp/soilc_par%k_ice/sec_year ! 1/s

    call nml_read(filename,"lnd_par","diff_bio",soilc_par%diff_bio)
    call nml_read(filename,"lnd_par","diff_cryo",soilc_par%diff_cryo)
    call nml_read(filename,"lnd_par","diff_ice",soilc_par%diff_ice)
    call nml_read(filename,"lnd_par","diff_shelf",soilc_par%diff_shelf)
    call nml_read(filename,"lnd_par","diff_lake",soilc_par%diff_lake)
    soilc_par%diff_bio   = soilc_par%diff_bio/sec_year ! m2/s
    soilc_par%diff_cryo  = soilc_par%diff_cryo/sec_year ! m2/s
    soilc_par%diff_ice   = soilc_par%diff_ice/sec_year ! m2/s
    soilc_par%diff_shelf = soilc_par%diff_shelf/sec_year ! m2/s
    soilc_par%diff_lake  = soilc_par%diff_lake/sec_year ! m2/s
    call nml_read(filename,"lnd_par","z_diff",soilc_par%z_diff)
    call nml_read(filename,"lnd_par","iadv_soilc",soilc_par%iadv_soilc)
    call nml_read(filename,"lnd_par","adv_soil",soilc_par%adv_soil)
    call nml_read(filename,"lnd_par","adv_ice",soilc_par%adv_ice)
    call nml_read(filename,"lnd_par","adv_shelf",soilc_par%adv_shelf)
    soilc_par%adv_soil  = soilc_par%adv_soil/sec_year ! m/s
    soilc_par%adv_ice   = soilc_par%adv_ice/sec_year ! m/s
    soilc_par%adv_shelf = soilc_par%adv_shelf/sec_year ! m/s

    call nml_read(filename,"lnd_par","peat_carb",peat_par%peat_carb)
    call nml_read(filename,"lnd_par","peat_area",peat_par%peat_area)
    call nml_read(filename,"lnd_par","nmonwet_peat",peat_par%nmonwet_peat)
    call nml_read(filename,"lnd_par","acroc_crit",peat_par%acroc_crit)
    call nml_read(filename,"lnd_par","k10_acro",peat_par%k10_acro)
    call nml_read(filename,"lnd_par","k10_cato",peat_par%k10_cato)
    call nml_read(filename,"lnd_par","k_acro_to_cato",peat_par%k_acro_to_cato)
    call nml_read(filename,"lnd_par","peat_ch_rate",peat_par%peat_ch_rate)
    call nml_read(filename,"lnd_par","dCpeat_dt_min",peat_par%dCpeat_dt_min)
    call nml_read(filename,"lnd_par","Cpeat_min",peat_par%Cpeat_min)
    peat_par%k_acro_to_cato = peat_par%k_acro_to_cato / sec_year ! 1/s
    peat_par%peat_ch_rate   = peat_par%peat_ch_rate / sec_year ! 1/s
    peat_par%dCpeat_dt_min  = peat_par%dCpeat_dt_min / sec_year ! kgC/m2/s
    peat_par%k10_acro   = 1._wp/peat_par%k10_acro/sec_year    ! 1/s
    peat_par%k10_cato   = 1._wp/peat_par%k10_cato/sec_year    ! 1/s

    call nml_read(filename,"lnd_par","ich4_ftemp",ch4_par%ich4_ftemp)
    call nml_read(filename,"lnd_par","q10_ch4",ch4_par%q10_ch4)
    call nml_read(filename,"lnd_par","Ea_ch4",ch4_par%Ea_ch4)
    call nml_read(filename,"lnd_par","ch4_frac_wet",ch4_par%ch4_frac_wet)
    call nml_read(filename,"lnd_par","ch4_frac_peat",ch4_par%ch4_frac_peat)
    call nml_read(filename,"lnd_par","ch4_frac_shelf",ch4_par%ch4_frac_shelf)
    call nml_read(filename,"lnd_par","ch4_frac_lake",ch4_par%ch4_frac_lake)

    call nml_read(filename,"lnd_par","lithology_uhh_file",lithology_uhh_file)
    call nml_read(filename,"lnd_par","lithology_gemco2_file",lithology_gemco2_file)
    call nml_read(filename,"lnd_par","i_weathering",i_weathering)
    call nml_read(filename,"lnd_par","i_weath_sc",i_weath_sc)
    call nml_read(filename,"lnd_par","kweath_scale",kweath_scale)
    call nml_read(filename,"lnd_par","l_river_export",l_river_export)
    call nml_read(filename,"lnd_par","kexport",kexport)
    call nml_read(filename,"lnd_par","l_match_alk_weath",l_match_alk_weath)

    call nml_read(filename,"lnd_par","write_surf",write_surf)
    call nml_read(filename,"lnd_par","write_surf_n",write_surf_n)
    call nml_read(filename,"lnd_par","write_carbon",write_carbon)
    call nml_read(filename,"lnd_par","write_soil",write_soil)
    call nml_read(filename,"lnd_par","write_soil_par",write_soil_par)
    call nml_read(filename,"lnd_par","write_lake",write_lake)
    call nml_read(filename,"lnd_par","write_cons",write_cons)
    call nml_read(filename,"lnd_par","l_daily_output",l_daily_output)

    return 

end subroutine lnd_par_load 

end module lnd_params

