!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o u p l e r
!
!  Purpose : coupling and exchange between all model components 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit, Andrey Ganopolski and
!                         Alexander Robinson
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
module coupler

    use precision, only : wp, dp
    use timer, only : doy, mon, year, sec_day, sec_year, nday_year, nmon_year, n_year_geo, nstep_mon_atm
    use timer, only : time_soy_atm, time_eoy_atm, time_soy_ocn, time_soy_sic, time_soy_lnd, time_eoy_lnd
    use timer, only : time_soy_bnd, time_soy_bgc, time_eoy_bgc, time_soy_smb, time_eoy_smb, time_eom_smb
    use timer, only : dt_atm
    use timer, only : time_call_daily_input_save, time_use_daily_input_save, nyear_avg_offline, nyears_spinup_bgc, nstep_year_ocn, nstep_year_bgc, n_accel
    use timer, only : monthly2daily
    use climber_grid, only : ni, nj, area
    use climber_grid, only : lon, lat, basin_mask, i_atlantic, i_pacific
    use control, only : ico2_rad, ich4_rad, id13c, iD14c, prc_forcing
    use control, only : ocn_restore_temp, ocn_restore_sal, atm_fix_tau
    use control, only : flag_co2, flag_ch4, flag_atm, flag_ocn, flag_bgc, flag_sic, flag_lnd, flag_dust, flag_smb, flag_bmb, flag_ice, flag_geo, flag_lakes
    use control, only : geo_restart
    use control, only : ifake_ice
    use control, only : l_spinup_cc
    use control, only : l_weathering
    use control, only : restart_in_dir
    use control, only : l_aqua_slab
    use constants, only : fqsat, q_sat_w, q_sat_i, Le, Lf, frac_vu
    use constants, only : rho_w, rho_sw, cap_w, rho_i, T0, c13_c12_std, c14_c_std, pi, ppm_to_PgC
    use constants, only : sigma, cap_a
    use coord, only : grid_class, grid_init
    use coord, only : map_scrip_class, map_scrip_init, map_scrip_field
    use ncio
    use filter, only : filter1d, smooth2

    use atm_def, only : atm_class
    use ocn_def, only : ocn_class
    use bgc_def, only : bgc_class, BGC_NTRA
    use sic_def, only : sic_class
    use lnd_def, only : lnd_class
    use ice_def, only : ice_class
    use co2_def, only : co2_class
    use ch4_def, only : ch4_class
    use smb_def, only : smb_in_class, smb_class
    use bmb_def, only : bmb_class
    use bnd_mod, only : bnd_class
    use geo_def, only : geo_class
    use lake_mod, only : lake_type
    use atm_params, only : i_zoro, p0
    use atm_grid, only : hatm, k1000, k900, k850, k700, k500
    use ocn_params, only : n_tracers_tot, n_tracers_ocn, tau_sst, tau_sss
    use ocn_params, only : l_ocn_input_fix, l_ocn_fix_wind, l_ocn_fix_fw, l_ocn_fix_flx, i_ocn_input_fix
    use ocn_params, only : n_cells_dist_runoff, n_cells_dist_calving, n_cells_dist_weath
    use ocn_params, only : relax_run, relax_calv, relax_bmelt, scale_runoff_ice
    use ocn_grid, only : maxi, maxj, maxk, k1_shelf, ocn_area_tot
    use lnd_params, only : dt_lnd => dt, l_ice_albedo_semi
    use lnd_grid, only : is_veg, is_ice, nl
    USE bgc_params, ONLY : l_sediments, l_spinup_bgc, i_compensate, l_conserve_phos, l_conserve_sil, l_conserve_alk, i_bgc_fw
    use bgc_params, only : iatmco2, iatmo2, iatmn2, iatmc13, iatmc14, isssc12, issssil, rcar
    use geo_params, only : h_ice_min

    !$use omp_lib

    implicit none 

    private
    public :: cmn_class
    public :: cmn_init
    public :: bnd_to_cmn
    public :: cmn_to_atm, atm_to_cmn
    public :: cmn_to_ocn, ocn_to_cmn
    public :: cmn_to_bgc, bgc_to_cmn
    public :: cmn_to_sic, sic_to_cmn
    public :: cmn_to_lnd, lnd_to_cmn
    public :: cmn_to_smb, smb_to_cmn
    public :: ice_to_smb, smb_to_ice 
    public :: cmn_to_bmb, bmb_to_cmn
    public :: ice_to_bmb, bmb_to_ice 
    public :: geo_to_bmb
    public :: ice_to_cmn
    public :: ice_to_geo, geo_to_ice
    public :: geo_to_smb, geo_to_cmn
    public :: bnd_to_geo, cmn_to_geo
    public :: cmn_to_co2, co2_to_cmn
    public :: cmn_to_ch4, ch4_to_cmn
    public :: lakes_update, runoff_merge, runoff_to_ocn
    public :: aquaplanet, aqua_init, aqua_end
    public :: nsurf

    integer, parameter :: nsurf = 10  !! total number of surface types in coupler [ocn,sic,5pfts,bare,lake,ice]
    integer, parameter :: i_surf_ocn = 1, i_surf_sic = 2, i_surf_lake = 9, i_surf_ice = 10  !! index of surface types in coupler
    integer, parameter :: nsurf_macro = 5       !! number of surface types in atmosphere model
    integer, parameter :: i_surf_macro_ocn = 1, i_surf_macro_sic = 2, i_surf_macro_lnd = 3, i_surf_macro_ice = 4, i_surf_macro_lake = 5  !! index of surface types in atmosphere model
    !integer, dimension(nsurf) :: st2ast = [1,2,3,3,3,3,3,3,3,4]    !! mapping from surface types to macro surface types for atmosphere (ocean, sea ice, land+lakes, ice sheets)
    integer, dimension(nsurf) :: st2ast = [1,2,3,3,3,3,3,3,5,4]    !! mapping from surface types to macro surface types for atmosphere (ocean, sea ice, land+lakes, ice sheets)
    integer, parameter :: nsurf_lnd = 8 !! number of surface types over land [5pfts,bare,lake,ice]
    integer, dimension(nsurf_lnd) :: i_surf_lnd = [3,4,5,6,7,8,9,10] !! index of surface types over land
    integer, parameter :: i_surf_lnd_ice = 8  !! index of ice surface type in land model
    integer, parameter :: i_surf_lnd_lake = 7  !! index of lake surface type in land model

    integer, dimension(nday_year) :: m0, m1
    real(wp), dimension(nday_year) :: wtm0, wtm1

    type cmn_class
      type(grid_class) :: grid
      integer :: n_lakes    !! number of lakes
      real(wp) :: co2   !! atmopheric CO2 concentration [ppm]
      real(wp) :: co2_rad   !! atmopheric CO2 concentration for radiation [ppm]
      real(wp) :: c13_c12_atm   !! atmospheric C13/C12 ratio []
      real(wp) :: c14_c_atm !! atmospheric C14/C ratio []
      real(wp) :: delta_C_lnd   !! net air-land carbon flux, positive = land uptake[kgC/yr]
      real(wp) :: delta_C13_lnd   !! net air-land carbon 13 flux, positive = land uptake[kgC13/yr]
      real(wp) :: delta_C14_lnd   !! net air-land carbon 14 flux, positive = land uptake[kgC14/yr]
      real(wp) :: delta_C_ocn   !! net air-sea carbon flux, positive = ocean uptake[kgC/yr]
      real(wp) :: delta_C13_ocn   !! net air-sea carbon 13 flux, positive = ocean uptake[kgC13/yr]
      real(wp) :: delta_C14_ocn   !! net air-sea carbon 14 flux, positive = ocean uptake[kgC14/yr]
      real(wp), dimension(:,:), allocatable :: delta_C_lnd_2d   !! net air-land carbon flux, positive = land uptake[kgC/s]
      real(wp), dimension(:,:), allocatable :: delta_C13_lnd_2d   !! net air-land carbon 13 flux, positive = land uptake[kgC13/s]
      real(wp), dimension(:,:), allocatable :: delta_C14_lnd_2d   !! net air-land carbon 14 flux, positive = land uptake[kgC14/s]
      real(wp), dimension(:,:), allocatable :: delta_C_ocn_2d   !! net air-sea carbon flux, positive = ocean uptake[kgC/s]
      real(wp), dimension(:,:), allocatable :: delta_C13_ocn_2d   !! net air-sea carbon 13 flux, positive = ocean uptake[kgC13/s]
      real(wp), dimension(:,:), allocatable :: delta_C14_ocn_2d   !! net air-sea carbon 14 flux, positive = ocean uptake[kgC14/s]
      real(wp) :: Cflx_lnd_avg  !! average net air-land carbon flux, positive = land uptake [kgC/yr]
      real(wp) :: Cflx_ocn_avg  !! average net air-sea carbon flux, positive = land uptake [kgC/yr]
      real(wp) :: weath_carb_avg !! average carbonate weathering carbon flux [kgC/yr]
      real(wp) :: weath_sil_avg  !! average silicate weathering carbon flux [kgC/yr] 
      real(wp) :: ch4   !! atmopheric CH4 concentration [ppb]
      real(wp) :: ch4_rad   !! atmopheric CH4 concentration for radiation [ppb]
      real(wp) :: ch4_flx_lnd   !! land CH4 emissions [kgCH4/yr]
      real(wp) :: ch4_emis  !! anthropogenic CH4 emissions [kgCH4/yr]
      real(wp) :: f_ch4emis_agro  !! fraction of anthropogenic CH4 emissions originating from agriculture [1]
      real(wp) :: n2o   !! atmopheric N2O concentration [ppb]
      real(wp) :: cfc11   !! atmopheric CFC11 concentration [ppt]
      real(wp) :: cfc12   !! atmopheric CFC12 concentration [ppt]
      real(wp), dimension(:,:), allocatable :: so4   !! atmopheric SO4 load [kg/m2]
      real(wp), dimension(:), allocatable :: o3_pl    !! vertical 'pressure' levels for atmospheric O3 concentration
      real(wp), dimension(:,:,:), allocatable :: o3    !! atmopheric O3 concentration [mol/mol]
      real(wp) :: sea_level !! sea level relative to present [m]
      real(wp) :: A_bering !! Bering Strait cross-sectional area [m2]
      real(wp) :: ocn_vol_tot !! actual ocean volume [m3]
      real(wp) :: buoy_sic_NA(6)   !! buoyancy flux from sea ice export over North Atlantic [N]
      real(wp) :: t2m_glob_ann  !! annual mean global temperature [K]
      real(wp) :: alk_to_ocn_glob    !! global alkalinity flux to ocean [molC/yr]
      real(wp) :: alk_from_lnd_glob    !! global alkalinity input flux from land from weathering [molC/yr]
      integer, dimension(:,:),   allocatable :: mask_ocn  !! ocean mask []
      integer, dimension(:,:),   allocatable :: mask_lnd  !! land mask []
      integer, dimension(:,:),   allocatable :: mask_smb  !! mask of coupler grid points covered by smb domain(s) []
      integer, dimension(:,:),   allocatable :: mask_ice   !! mask of coupler grid points covered by ice sheet domain(s) []
      integer, dimension(:,:),   allocatable :: mask_coast   !! mask of coastal cells []
      integer, dimension(:,:),   allocatable :: mask_coast2   !! mask of coastal cells []
      real(wp), dimension(:,:),   allocatable :: f_ocn  !! ocean fraction of grid cell, including floating ice []
      real(wp), dimension(:,:),   allocatable :: f_ocn2  !! ocean fraction of grid cell, excluding floating ice []
      real(wp), dimension(:,:),   allocatable :: f_sic  !! sea ice fraction of ocean fraction (f_ocn2) []
      real(wp), dimension(:,:),   allocatable :: f_lnd0 !! land fraction of grid cell for present day sea level []
      real(wp), dimension(:,:),   allocatable :: f_lnd  !! land fraction of grid cell []
      real(wp), dimension(:,:),   allocatable :: f_ice  !! ice sheet/glacier fraction of grid cell []
      real(wp), dimension(:,:),   allocatable :: f_ice_grd  !! grounded ice sheet/glacier fraction of grid cell []
      real(wp), dimension(:,:),   allocatable :: f_ice_flt  !! flating ice fraction of grid cell []
      real(wp), dimension(:,:),   allocatable :: f_lake !! lake fraction of grid cell []
      real(wp), dimension(:,:,:), allocatable :: f_lake_n !! individual lake fraction in each grid cell []
      real(wp), dimension(:,:),   allocatable :: f_crop !! crop fraction of land fraction []
      real(wp), dimension(:,:),   allocatable :: f_pasture !! pasture fraction of land fraction []
      real(wp), dimension(:,:,:), allocatable :: disturbance !! vegetation disturbance rate [1/s]
      real(wp), dimension(:,:,:), allocatable :: f_stp  !! surface types fraction []
      real(wp), dimension(:,:,:), allocatable :: f_astp  !! surface types fraction in atmosphere model []
      real(wp), dimension(:,:),   allocatable :: z_sur   !! mean surface elevation of grid cell [m]
      real(wp), dimension(:,:,:), allocatable :: z_sur_n   !! mean surface elevation of each surface type in atmosphere model [m]
      real(wp), dimension(:,:),   allocatable :: z_ocn_max   !! max elevation of ocean, Q10 [m]
      real(wp), dimension(:,:),   allocatable :: z_veg   !! mean surface elevation of ice free land [m]
      real(wp), dimension(:,:),   allocatable :: z_veg_min   !! min surface elevation of ice free land [m]
      real(wp), dimension(:,:),   allocatable :: z_veg_max   !! max surface elevation of ice free land [m]
      real(wp), dimension(:,:),   allocatable :: z_ice   !! mean surface elevation of ice sheet [m]
      real(wp), dimension(:,:),   allocatable :: z_lake  !! mean surface elevation of lakes [m]
      real(wp), dimension(:,:),   allocatable :: z_sur_smooth_std   !! standard deviation of smoothed surface elevation of grid cell [m]
      real(wp), dimension(:,:),   allocatable :: z_veg_std   !! standard deviation of surface elevation of grid cell, ice-free land only [m]
      real(wp), dimension(:,:,:), allocatable :: coral_f_area  !! ocean hypsometry for corals []
      real(wp), dimension(:,:,:), allocatable :: coral_f_topo  !! topography factor for corals []
      integer, dimension(:,:),    allocatable :: coast_nbr     !! number of neighbors of coastal cells []
      integer, dimension(:,:,:),  allocatable :: i_coast_nbr     !! i index of neighbors of coastal cells []
      integer, dimension(:,:,:),  allocatable :: j_coast_nbr     !! j index of neighbors of coastal cells []
      real(wp), dimension(:,:),   allocatable :: q_geo   !! geothermal heat flux [W/m2]
      real(wp), dimension(:,:),   allocatable :: sst    !! sea surface temperature [degC]
      real(wp), dimension(:,:),   allocatable :: sss    !! sea surface salinity [psu]
      real(wp), dimension(:,:),   allocatable :: sss_dat    !! sea surface salinity for ocean surface restoring setup [psu]
      real(wp), dimension(:,:),   allocatable :: sst_dat    !! sea surface temperature for ocean surface restoring setup [degC]
      real(wp), dimension(:,:),   allocatable :: flx_ocn    !! heat flux into the ocean [W/m2]
      real(wp), dimension(:,:),   allocatable :: qflux      !! q-flux for slab ocean [W/m2]
      real(wp), dimension(:,:),   allocatable :: p_e_sic_ocn !! freshwater flux into the ocean (excluding runoff) [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: fw_brines !! freshwater flux related to brine rejection [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: t_shelf    !! temperature on ocean shelf surface [K]
      real(wp), dimension(:),     allocatable :: z_ocn_bmb    !! depth of ocean layers for basal melt [m]
      integer,  dimension(:,:,:), allocatable :: mask_ocn_bmb    !! ocean mask for basal melt [1]
      real(wp), dimension(:,:,:), allocatable :: t_ocn_bmb    !! ocean temperature for basal melt [degC]
      real(wp), dimension(:,:,:), allocatable :: s_ocn_bmb    !! ocean salinity for basal melt [psu]
      real(wp), dimension(:,:,:), allocatable :: taux   !! zonal surface wind stress [N/m2]
      real(wp), dimension(:,:,:), allocatable :: tauy   !! meridional surface wind stress [N/m2] 
      real(wp), dimension(:,:),   allocatable :: tauxo  !! zonal sea ice stress on ocean [N/m2]
      real(wp), dimension(:,:),   allocatable :: tauyo  !! meridional sea ice stress on ocean [N/m2] 
      real(wp), dimension(:,:),   allocatable :: rbatm_dat   !! radiative balance of atmosphere for fix rad setup [W/m2]
      real(wp), dimension(:,:),   allocatable :: cld_dat   !! cloud fraction for fix rad setup [/]
      real(wp), dimension(:,:),   allocatable :: cld_day_dat   !! daytime cloud fraction for fix rad setup [/]
      real(wp), dimension(:,:),   allocatable :: slp_dat   !! sea level pressure [Pa]
      real(wp), dimension(:,:),   allocatable :: tsl_dat   !! sea level temperature [K]
      real(wp), dimension(:,:),   allocatable :: htrop_dat   !! tropopause height [m]
      real(wp), dimension(:,:),   allocatable :: taux_dat   !! zonal surface wind stress for fix tau setup [N/m2]
      real(wp), dimension(:,:),   allocatable :: tauy_dat   !! meridional surface wind stress for fix tau setup [N/m2] 
      real(wp), dimension(:,:),   allocatable :: uo1    !! zonal surface ocean velocity on u-grid [m/s]
      real(wp), dimension(:,:),   allocatable :: vo1    !! meridional surface ocean velocity on v-grid [m/s]
      real(wp), dimension(:,:),   allocatable :: ssh !! elevation of the ocean free surface [m]
      real(wp), dimension(:,:,:), allocatable :: t_skin !! skin temperature over each surface type [K]
      real(wp), dimension(:,:,:), allocatable :: t2m   !! surface air temperature over each surface type [K]
      real(wp), dimension(:,:,:), allocatable :: q2m   !! surface air specific humidity over each surface type [kg/kg]
      real(wp), dimension(:,:,:), allocatable :: t2m_mon  !! monthly mean temperature [K]
      real(wp), dimension(:,:,:), allocatable :: t2m_mon_lnd  !! monthly mean temperature over icefree land [K]
      real(wp), dimension(:,:),   allocatable :: t2m_min_mon  !! minimum monthly temperature over icefree land [K]
      real(wp), dimension(:,:),   allocatable :: tam   !! atmospheric air temperature [K]
      real(wp), dimension(:,:),   allocatable :: ram   !! atmospheric air relative humidity [1]
      real(wp), dimension(:,:),   allocatable :: gams   !! atmospheric temperature lapse rate at the surface [K/m]
      real(wp), dimension(:,:),   allocatable :: t_soil   !! mean annual temperature of bottom layer of land model [K]
      real(wp), dimension(:,:),   allocatable :: prc    !! total precipitation rate [kg/m2/s] 
      real(wp), dimension(:,:,:),   allocatable :: rain   !! rainfall rate [kg/m2/s] 
      real(wp), dimension(:,:,:),   allocatable :: snow   !! snowfall rate [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: cld   !! cloud fraction [/]
      real(wp), dimension(:,:),   allocatable :: cld_day   !! daytime cloud fraction [/]
      real(wp), dimension(:,:,:), allocatable :: ps    !! surface pressure [Pa]
      real(wp), dimension(:,:),   allocatable :: slp   !! sea level pressure [Pa]
      real(wp), dimension(:,:),   allocatable :: usur   !! surface zonal wind speed [m/s]
      real(wp), dimension(:,:),   allocatable :: vsur   !! surface meridional wind speed [m/s]
      real(wp), dimension(:,:,:), allocatable :: wind   !! surface wind speed [m/s]
      real(wp), dimension(:,:),   allocatable :: u700   !! zonal wind speed at 700 hPa [m/s]
      real(wp), dimension(:,:),   allocatable :: v700   !! meridional wind speed at 700 hPa [m/s]
      real(wp), dimension(:,:,:), allocatable :: z0m    !! surface roughness length for momentum [m]
      real(wp), dimension(:,:,:), allocatable :: swnet  !! net surface shortwave radiation over each surface type [W/m2]
      real(wp), dimension(:,:),   allocatable :: swd    !! surface downward flux of shortwave radiation [W/m2]
      real(wp), dimension(:,:),   allocatable :: swd_vis_dir    !! surface downward flux of shortwave visible radiation, clear sky [W/m2]
      real(wp), dimension(:,:),   allocatable :: swd_nir_dir    !! surface downward flux of shortwave near-infrared radiation, clear sky [W/m2]
      real(wp), dimension(:,:),   allocatable :: swd_vis_dif    !! surface downward flux of shortwave visible radiation, cloudy sky [W/m2]
      real(wp), dimension(:,:),   allocatable :: swd_nir_dif    !! surface downward flux of shortwave near-infrared radiation, cloudy sky [W/m2]
      real(wp), dimension(:,:),   allocatable :: dswd_dalb_vis_dir
      real(wp), dimension(:,:),   allocatable :: dswd_dalb_nir_dir
      real(wp), dimension(:,:),   allocatable :: dswd_dalb_vis_dif
      real(wp), dimension(:,:),   allocatable :: dswd_dalb_nir_dif
      real(wp), dimension(:,:),   allocatable :: dswd_dz_nir_dir
      real(wp), dimension(:,:),   allocatable :: dswd_dz_nir_dif
      real(wp), dimension(:,:,:), allocatable :: alb_vis_dir    !! visible direct beam surface albedo of each surface type []
      real(wp), dimension(:,:,:), allocatable :: alb_vis_dif    !! visible diffuse surface albedo of each surface type []
      real(wp), dimension(:,:,:), allocatable :: alb_nir_dir    !! near infrared direct beam surface albedo of each surface type []
      real(wp), dimension(:,:,:), allocatable :: alb_nir_dif    !! near infrared diffuse surface albedo of each surface type []
      real(wp), dimension(:,:,:), allocatable :: alb_vis_dir_ice_semi    !! ice sheet visible direct beam albedo from SEMI interpolate to cmn grid []
      real(wp), dimension(:,:,:), allocatable :: alb_vis_dif_ice_semi    !! ice sheet visible diffuse albedo from SEMI interpolate to cmn grid[]
      real(wp), dimension(:,:,:), allocatable :: alb_nir_dir_ice_semi    !! ice sheet near infrared direct beam albedo from SEMI interpolate to cmn grid []
      real(wp), dimension(:,:,:), allocatable :: alb_nir_dif_ice_semi    !! ice sheet near infrared diffuse albedo from SEMI interpolate to cmn grid []
      real(wp), dimension(:,:,:), allocatable :: lwd    !! surface downward flux of longwave radiation [W/m2]
      real(wp), dimension(:,:,:), allocatable :: lwd_cs    !! surface downward flux of longwave radiation, clear sky [W/m2]
      real(wp), dimension(:,:,:), allocatable :: lwd_cld    !! surface downward flux of longwave radiation, cloudy sky [W/m2]
      real(wp), dimension(:,:,:), allocatable :: lwu    !! upward surface longwave radiation over each surface type [W/m2]
      real(wp), dimension(:,:,:), allocatable :: sh !! sensible heat flux over each surface type [W/m2]
      real(wp), dimension(:,:,:), allocatable :: lh !! latent heat flux over each surface type [W/m2]
      real(wp), dimension(:,:,:), allocatable :: evp    !! evaporation from each surface type [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: runoff_veg  !! runoff from vegetated grid cell parts [kg/s]
      real(wp), dimension(:,:),   allocatable :: calving_veg !! calving from vegetated grid cell parts [kg/s]
      real(wp), dimension(:,:),   allocatable :: runoff_ice  !! runoff from ice sheets [kg/s]
      real(wp), dimension(:,:),   allocatable :: melt_ice  !! ice sheet melt [kg/s]
      real(wp), dimension(:,:),   allocatable :: calving_ice !! calving from ice sheets [kg/s]
      real(wp), dimension(:,:),   allocatable :: bmelt_grd !! basal melt from grounded ice sheets [kg/s]
      real(wp), dimension(:,:),   allocatable :: bmelt_flt !! basal melt from floating ice shelfs [kg/s]
      real(wp), dimension(:,:),   allocatable :: runoff_ice_l  !! runoff from ice sheets computed by land model [kg/s]
      real(wp), dimension(:,:),   allocatable :: calving_ice_l !! calving from ice sheets computed by land model [kg/s]
      real(wp), dimension(:,:,:), allocatable :: melt_ice_i_mon  !! monthly ice sheet melt [kg/s]
      real(wp), dimension(:,:,:), allocatable :: acc_ice_i_mon  !! monthly ice sheet net accumulation from prescribed bnd ice thickness changes [kg/s]
      real(wp), dimension(:,:,:), allocatable :: runoff_ice_i_mon  !! monthly runoff from ice sheets [kg/s]
      real(wp), dimension(:,:),   allocatable :: calving_ice_i !! calving from ice sheets computed by ice model [kg/s]
      real(wp), dimension(:,:),   allocatable :: bmelt_grd_i !! basal melt from grounded ice sheets computed by ice model [kg/s]
      real(wp), dimension(:,:),   allocatable :: bmelt_flt_i !! basal melt from floating ice shelfs computed by ice model [kg/s]
      real(wp), dimension(:,:),   allocatable :: runoff_o  !! runoff on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: runoff_o_ann  !! annual runoff on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: calving_o !! calving on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: bmelt_o !! basal melt from ice sheets on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: bmelt_grd_o !! basal melt from grounded ice sheets on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: bmelt_flt_o !!  basal melt from floating ice shelfs on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: runoff_veg_o  !! runoff from vegetated grid cell parts on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: calving_veg_o !! calving from vegetated grid cell parts on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: runoff_ice_o  !! runoff from ice sheets on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: melt_ice_o  !! freshwater flux to the ocean from ice sheet melt [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: calving_ice_o !! calving from ice sheets on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: runoff_lake_o  !! runoff from lakes on ocean domain [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: calving_lake_o !! calving from lakes on ocean domain [kg/m2/s]
      integer, dimension(:,:), allocatable :: i_runoff, j_runoff!! i and j indexes of land runoff destination gridcell 
      integer, dimension(:,:), allocatable :: i_runoff_veg, j_runoff_veg  !! i and j indexes of land runoff destination gridcell 
      integer, dimension(:,:), allocatable :: i_runoff_ice, j_runoff_ice  !! i and j indexes of ice sheet runoff destination gridcell 
      real(wp), dimension(:,:,:),   allocatable :: f_drain_veg !! drainage fractions from land to lakes and ocean []
      real(wp), dimension(:,:,:),   allocatable :: f_drain_ice !! drainage fractions from lice sheets to lakes and ocean []
      type(lake_type), dimension(:), allocatable :: lake  !! lake derived type
      real(wp), dimension(:,:),   allocatable :: lake_p_e  !! P-E for lake [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: h_lake  !! average lake(s) depth on coupler grid [m]
      real(wp), dimension(:,:),   allocatable :: f_ice_lake  !! ice fraction over lake [1]
      real(wp), dimension(:),     allocatable :: z_lake_bmb    !! depth of lake layers for basal melt [m]
      integer,  dimension(:,:), allocatable :: mask_lake_bmb    !! lake mask for basal melt [1]
      real(wp), dimension(:,:,:), allocatable :: t_lake_bmb    !! lake temperature for basal melt [degC]
      real(wp), dimension(:,:,:), allocatable :: s_lake_bmb    !! lake salinity for basal melt [psu]
      integer, dimension(:), allocatable :: idivide_pac_atl 
      integer, dimension(:), allocatable :: idivide_atl_indpac
      real(wp) :: eccentricity
      real(wp) :: precession
      real(wp) :: obliquity
      real(wp), dimension(:,:,:),   allocatable :: cosz  !! cosine of solar zenith angle, with diurnal cycle []
      real(wp), dimension(:,:),   allocatable :: coszm  !! radiation weighted mean cosine of solar zenith angle []
      real(wp), dimension(:,:,:),   allocatable :: solar  !! insolation at top of the atmosphere, with diurnal cycle [W/m2]
      real(wp), dimension(:,:),   allocatable :: solarm     !! mean insolation at top of the atmosphere [W/m2]
      real(wp), dimension(:,:),   allocatable :: solarmin     !! minimum diurnal insolation at top of the atmosphere [W/m2]
      real(wp), dimension(:,:),   allocatable :: solarmax     !! maximum diurnal insolation at top of the atmosphere [W/m2]
      real(wp), dimension(:,:),   allocatable :: daylength  !! length of day [hours]
      real(wp), dimension(:,:),   allocatable :: dust_emis  !! dust emission [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: dust_dep   !! dust deposition [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: poc_export    !! particulate organic carbon ocean input from rivers [kgC/m2/yr]
      real(wp), dimension(:,:),   allocatable :: poc13_export  !! particulate organic carbon 13 ocean input from rivers [kgC/m2/yr]
      real(wp), dimension(:,:),   allocatable :: poc14_export  !! particulate organic carbon 14 ocean input from rivers [kgC/m2/yr]
      real(wp), dimension(:,:),   allocatable :: doc_export    !! dissolved organic carbon ocean input from rivers [kgC/m2/yr]
      real(wp), dimension(:,:),   allocatable :: doc13_export  !! dissolved organic carbon 13 ocean input from rivers [kgC/m2/yr]
      real(wp), dimension(:,:),   allocatable :: doc14_export  !! dissolved organic carbon 14 ocean input from rivers [kgC/m2/yr]
      real(wp), dimension(:,:),   allocatable :: weath_carb   !! calcium ocean input from carbonate weathering [mol C/m2/yr]
      real(wp), dimension(:,:),   allocatable :: weath_sil    !! silica ocean input from silicate weathering [mol C/m2/yr]
      real(wp), dimension(:,:),   allocatable :: weath13_carb   !! C13 ocean input from carbonate weathering [mol C/m2/yr]
      real(wp), dimension(:,:),   allocatable :: weath13_sil    !! C13 ocean input from silicate weathering [mol C/m2/yr]
      real(wp), dimension(:,:),   allocatable :: weath14_carb   !! C14 ocean input from carbonate weathering [mol C/m2/yr]
      real(wp), dimension(:,:),   allocatable :: weath14_sil    !! C14 ocean input from silicate weathering [mol C/m2/yr]
      real(wp), dimension(:,:),   allocatable :: f_carb   !! carbonate fraction on shelf (for weathering) [1]
    end type cmn_class

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ a t m
  ! Purpose    :  from common grid to atmosphere
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_atm(cmn,atm)

    implicit none

    type(cmn_class) :: cmn
    type(atm_class) :: atm

    integer :: i, j, k, kk, n, ja
    real(wp) :: zr_oro
    integer, allocatable, dimension(:), save  :: ki1, ki2
    real(wp), allocatable, dimension(:), save :: w1, w2

    logical, save :: firstcall=.true.


    if (firstcall) then
      allocate(ki1(size(atm%pl)))
      allocate(ki2(size(atm%pl)))
      allocate(w1(size(atm%pl)))
      allocate(w2(size(atm%pl)))
      do k=1,size(atm%pl)
        kk = 1
        do while (cmn%o3_pl(kk).ge.atm%pl(k) .and. (kk+1).le.size(cmn%o3_pl))
          kk=kk+1
        enddo
        ki1(k) = max(1,kk-1)
        ki2(k) = min(size(cmn%o3_pl),kk)
        w1(k)  = 1._wp - (atm%pl(k)-cmn%o3_pl(ki1(k))) / (cmn%o3_pl(ki2(k))-cmn%o3_pl(ki1(k)))
        w2(k)  = 1._wp-w1(k)
      enddo
      firstcall = .false.
    endif

    if (time_soy_atm) then
      ! orbital parameters
      atm%eccentricity = cmn%eccentricity
      atm%precession   = cmn%precession
      atm%obliquity    = cmn%obliquity
      ! well mixed greenhouse gas concentrations
      atm%co2 = cmn%co2_rad
      atm%ch4 = cmn%ch4_rad
      atm%n2o = cmn%n2o
      atm%cfc11 = cmn%cfc11
      atm%cfc12 = cmn%cfc12
    endif

    !$omp parallel do collapse(2) private(i,j,ja,k,n,zr_oro)
    do j=1,nj
      do i=1,ni

        ! NOTE: atmospheric grid is reversed in north-south direction compared to the coupler!
        ja = nj-j+1

        if (time_soy_atm) then

          atm%frlnd(i,j) = cmn%f_lnd(i,ja) 
          atm%frocn(i,j) = cmn%f_ocn2(i,ja)
          atm%sigoro(i,j) = cmn%z_sur_smooth_std(i,ja)

          atm%idivide_pac_atl(j) = cmn%idivide_pac_atl(ja)
          atm%idivide_atl_indpac(j) = cmn%idivide_atl_indpac(ja)

          atm%zs(i,j,i_surf_macro_ocn) = 0._wp     ! ocean surface elevation
          atm%zs(i,j,i_surf_macro_sic) = 0._wp     ! sea ice surface elevation
          atm%zs(i,j,i_surf_macro_lnd) = cmn%z_veg(i,ja)  ! ice-free land surface elevation 
          atm%zs(i,j,i_surf_macro_ice) = cmn%z_ice(i,ja)  ! ice-sheet surface elevation
          atm%zs(i,j,i_surf_macro_lake) = cmn%z_lake(i,ja)  ! lakes surface elevation

          atm%zsa(i,j) = cmn%z_sur(i,ja)   ! grid cell mean surface elevation

          ! insolation related variables
          atm%solar(:,:,j) = cmn%solar(:,:,ja)
          atm%solarm(:,j)  = cmn%solarm(:,ja)
          atm%cosz(:,:,j)  = cmn%cosz(:,:,ja)
          atm%coszm(:,j)   = cmn%coszm(:,ja)


          ! SO4 load
          atm%so4(i,j) = cmn%so4(i,ja)

          ! O3 concentration, interpolated to atmospheric levels
          do k=1,size(atm%pl)
            atm%o3(i,j,k) = w1(k)*cmn%o3(i,ja,ki1(k)) + w2(k)*cmn%o3(i,ja,ki2(k))
          enddo

        endif

        ! TODO, for open cc add weathering and degassing fluxes!
        atm%co2flx(i,j)   = -(cmn%delta_C_ocn_2d(i,ja) + cmn%delta_C_lnd_2d(i,ja)) * 44.0095_wp/12._wp / area(i,ja)    ! kgC/s * kgCO2/kgC / m2 = kgCO2/m2/s
        ! TODO isotopes?
        !atm%C13flx(i,j) = cmn%delta_C13_ocn_2d(i,ja) + cmn%delta_C13_lnd_2d(i,ja)  
        !atm%C14flx(i,j) = cmn%delta_C14_ocn_2d(i,ja) + cmn%delta_C14_lnd_2d(i,ja)  

        if (abs(sum(cmn%f_stp(i,ja,:))-1._wp)>1.e-5) then
          print *
          print *,'sum(f_stp).ne.1!!!',sum(cmn%f_stp(i,ja,:))
          print *,i,ja
          print *,'f_ocn2,f_lnd,f_ocn',cmn%f_ocn2(i,ja), cmn%f_lnd(i,ja), cmn%f_ocn(i,ja)
          print *,'f_ice,f_ice_flt,f_ice_grd',cmn%f_ice(i,ja),cmn%f_ice_flt(i,ja),cmn%f_ice_grd(i,ja)
          print *,'f_stp',cmn%f_stp(i,ja,:)
        endif

        ! 4 surface types are distinguished for radiation in the atmophere model (ocean, sea ice, land, ice sheets)
        ! the 4 surface types generally differ in surface properties and elevation (important for LW radiation!)
        do n=1,nsurf_macro

          ! surface type fraction conversion 
          atm%frst(i,j,n) = sum(cmn%f_stp(i,ja,:), mask=(st2ast==n))

          if (atm%frst(i,j,n).gt.0._wp) then
            ! average albedoes of macro surface types
            atm%alb_vu_s(i,j,n) = sum(cmn%alb_vis_dir(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==n)) / atm%frst(i,j,n) 
            atm%alb_vu_c(i,j,n) = sum(cmn%alb_vis_dif(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==n)) / atm%frst(i,j,n)
            atm%alb_ir_s(i,j,n) = sum(cmn%alb_nir_dir(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==n)) / atm%frst(i,j,n)
            atm%alb_ir_c(i,j,n) = sum(cmn%alb_nir_dif(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==n)) / atm%frst(i,j,n)
            ! upward longwave radiation at the surface
            atm%flwr_up_sur(i,j,n) = sum(cmn%lwu(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==n)) / atm%frst(i,j,n)
            ! skin temperature 
            atm%tskin(i,j,n) = sum(cmn%t_skin(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==n)) / atm%frst(i,j,n)
          else
            atm%alb_vu_s(i,j,n) = 0._wp 
            atm%alb_vu_c(i,j,n) = 0._wp
            atm%alb_ir_s(i,j,n) = 0._wp
            atm%alb_ir_c(i,j,n) = 0._wp
            atm%flwr_up_sur(i,j,n) = 0._wp
            atm%tskin(i,j,n) = 0._wp 
          endif

        enddo

        ! incoming SW radiation at TOA
        atm%swr_dw_top(i,j) = cmn%solarm(doy,ja)

        ! grid average sensible heat flux
        atm%sha(i,j)  = sum(cmn%sh(i,ja,:)*cmn%f_stp(i,ja,:))
        ! grid average latent heat flux
        atm%lha(i,j)  = sum(cmn%lh(i,ja,:)*cmn%f_stp(i,ja,:))
        ! grid average evaporation
        atm%evpa(i,j) = sum(cmn%evp(i,ja,:)*cmn%f_stp(i,ja,:))

        ! orographic roughness
        if (i_zoro.eq.0) then
          zr_oro = 0._wp 
        else if (i_zoro.eq.1) then
          ! orographic roughness length, Hansen 1983
          zr_oro = 0.041*cmn%z_sur_smooth_std(i,ja)**0.71 
        else if (i_zoro.eq.2) then
          ! orographic roughness length, Hansen 1983, mod
          zr_oro = 0.041*(0.5_wp*cmn%z_sur_smooth_std(i,ja))**0.71 
        else if (i_zoro.eq.3) then
          ! orographic roughness length, linear
          zr_oro = 0.004*cmn%z_sur_smooth_std(i,ja)
        else if (i_zoro.eq.4) then
          ! orographic roughness length, linear
          zr_oro = 0.004*cmn%z_sur_smooth_std(i,ja)*(1._wp-cmn%f_ice(i,ja))
        endif
        atm%zoro(i,j) = zr_oro

        ! surface roughness length for momentum and heat/moisture and drag coefficient for heat/moisture
        ! ocean
        atm%z0m(i,j,i_surf_macro_ocn) = cmn%z0m(i,ja,i_surf_macro_ocn)
        ! sea ice
        atm%z0m(i,j,i_surf_macro_sic) = cmn%z0m(i,ja,i_surf_macro_sic)
        ! land 
        if (atm%frst(i,j,i_surf_macro_lnd).gt.0._wp) then
          atm%z0m(i,j,i_surf_macro_lnd) = sum(cmn%z0m(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==i_surf_macro_lnd)) / atm%frst(i,j,i_surf_macro_lnd)    ! ice-free land
        else
          atm%z0m(i,j,i_surf_macro_lnd)  = 1.e-3_wp
        endif
        ! ice sheets
        if (atm%frst(i,j,i_surf_macro_ice).gt.0._wp) then
          atm%z0m(i,j,i_surf_macro_ice) = sum(cmn%z0m(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==i_surf_macro_ice)) / atm%frst(i,j,i_surf_macro_ice) 
        else
          atm%z0m(i,j,i_surf_macro_ice)  = 1.e-3_wp
        endif
        ! lakes
        if (atm%frst(i,j,i_surf_macro_lake).gt.0._wp) then
          atm%z0m(i,j,i_surf_macro_lake) = sum(cmn%z0m(i,ja,:)*cmn%f_stp(i,ja,:), mask=(st2ast==i_surf_macro_lake)) / atm%frst(i,j,i_surf_macro_lake)    
        else
          atm%z0m(i,j,i_surf_macro_lake)  = 1.e-3_wp
        endif


        if (flag_dust) then
          ! dust emissions 
          atm%dust_emis(i,j) = cmn%dust_emis(i,ja) 
        else
          ! no emissions and pass dust deposition for diagnostics only
          atm%dust_emis(i,j) = 0._wp
          atm%dust_dep(i,j) = cmn%dust_dep(i,ja) 
        endif

        ! lake ice fraction
        atm%f_ice_lake(i,j) = cmn%f_ice_lake(i,ja) 
        
      enddo
    enddo
    !$omp end parallel do

   return

  end subroutine cmn_to_atm


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  a t m _ t o _ c m n
  ! Purpose    :  from atmosphere to common grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_to_cmn(atm,cmn)

    implicit none

    type(cmn_class) :: cmn
    type(atm_class) :: atm

    integer :: i, j, n, ja


    !$omp parallel do collapse(2) private(i,j,n,ja)
    do j=1,nj
      do i=1,ni
        ja = nj-j+1 ! flip N-S

        cmn%f_astp(i,j,:)  = atm%frst(i,ja,:)

        cmn%slp(i,j)  = atm%slp(i,ja)
        cmn%prc(i,j)  = atm%prc(i,ja)
        cmn%cld(i,j)  = atm%cld(i,ja)   
        cmn%tam(i,j)  = atm%tam(i,ja) 
        cmn%ram(i,j)  = atm%ram(i,ja) 
        cmn%gams(i,j) = atm%gams(i,ja) 

        cmn%ps(i,j,i_surf_macro_ocn) = atm%ps(i,ja,i_surf_macro_ocn)
        cmn%ps(i,j,i_surf_macro_sic) = atm%ps(i,ja,i_surf_macro_sic)
        do n=1,nsurf_lnd
          if (n.ne.i_surf_lnd_ice .and. n.ne.i_surf_lnd_lake) then
            cmn%ps(i,j,i_surf_lnd(n)) = atm%ps(i,ja,i_surf_macro_lnd)
          endif
        enddo
        cmn%ps(i,j,i_surf_ice) = atm%ps(i,ja,i_surf_macro_ice)
        cmn%ps(i,j,i_surf_lake) = atm%ps(i,ja,i_surf_macro_lake)

        cmn%rain(i,j,i_surf_macro_ocn) = atm%prcw(i,ja,i_surf_macro_ocn)
        cmn%rain(i,j,i_surf_macro_sic) = atm%prcw(i,ja,i_surf_macro_sic)
        do n=1,nsurf_lnd
          if (n.ne.i_surf_lnd_ice .and. n.ne.i_surf_lnd_lake) then
            cmn%rain(i,j,i_surf_lnd(n)) = atm%prcw(i,ja,i_surf_macro_lnd)
          endif
        enddo
        cmn%rain(i,j,i_surf_ice) = atm%prcw(i,ja,i_surf_macro_ice)
        cmn%rain(i,j,i_surf_lake) = atm%prcw(i,ja,i_surf_macro_lake)

        cmn%snow(i,j,i_surf_macro_ocn) = atm%prcs(i,ja,i_surf_macro_ocn)
        cmn%snow(i,j,i_surf_macro_sic) = atm%prcs(i,ja,i_surf_macro_sic)
        do n=1,nsurf_lnd
          if (n.ne.i_surf_lnd_ice .and. n.ne.i_surf_lnd_lake) then
            cmn%snow(i,j,i_surf_lnd(n)) = atm%prcs(i,ja,i_surf_macro_lnd)
          endif
        enddo
        cmn%snow(i,j,i_surf_ice) = atm%prcs(i,ja,i_surf_macro_ice)
        cmn%snow(i,j,i_surf_lake) = atm%prcs(i,ja,i_surf_macro_lake)

        cmn%lwd(i,j,i_surf_macro_ocn)   = atm%flwr_dw_sur(i,ja,i_surf_macro_ocn)
        cmn%lwd(i,j,i_surf_macro_sic)   = atm%flwr_dw_sur(i,ja,i_surf_macro_sic)
        do n=1,nsurf_lnd
          if (n.ne.i_surf_lnd_ice .and. n.ne.i_surf_lnd_lake) then
            cmn%lwd(i,j,i_surf_lnd(n))   = atm%flwr_dw_sur(i,ja,i_surf_macro_lnd)
          endif
        enddo
        cmn%lwd(i,j,i_surf_ice)   = atm%flwr_dw_sur(i,ja,i_surf_macro_ice)
        cmn%lwd(i,j,i_surf_lake)   = atm%flwr_dw_sur(i,ja,i_surf_macro_lake)

        cmn%lwd_cs(i,j,:)  = atm%flwr_dw_sur_cs(i,ja,:)
        cmn%lwd_cld(i,j,:) = atm%flwr_dw_sur_cld(i,ja,:)

        cmn%swnet(i,j,i_surf_macro_ocn) = atm%fswr_sur(i,ja,i_surf_macro_ocn)
        cmn%swnet(i,j,i_surf_macro_sic) = atm%fswr_sur(i,ja,i_surf_macro_sic)
        do n=1,nsurf_lnd
          if (n.ne.i_surf_lnd_ice .and. n.ne.i_surf_lnd_lake) then
            cmn%swnet(i,j,i_surf_lnd(n)) = atm%fswr_sur(i,ja,i_surf_macro_lnd)
          endif
        enddo
        cmn%swnet(i,j,i_surf_ice)   = atm%fswr_sur(i,ja,i_surf_macro_ice)
        cmn%swnet(i,j,i_surf_lake)   = atm%fswr_sur(i,ja,i_surf_macro_lake)

        cmn%swd_vis_dir(i,j) = atm%swr_dw_sur_vis_cs(i,ja)
        cmn%swd_nir_dir(i,j) = atm%swr_dw_sur_nir_cs(i,ja)
        cmn%swd_vis_dif(i,j) = atm%swr_dw_sur_vis_cld(i,ja)
        cmn%swd_nir_dif(i,j) = atm%swr_dw_sur_nir_cld(i,ja)

        cmn%dswd_dalb_vis_dir(i,j) = atm%dswd_dalb_vu_cs(i,ja)  
        cmn%dswd_dalb_nir_dir(i,j) = atm%dswd_dalb_ir_cs(i,ja)
        cmn%dswd_dalb_vis_dif(i,j) = atm%dswd_dalb_vu_cld(i,ja)
        cmn%dswd_dalb_nir_dif(i,j) = atm%dswd_dalb_ir_cld(i,ja)
        cmn%dswd_dz_nir_dir(i,j)   = atm%dswd_dz_ir_cs(i,ja)  
        cmn%dswd_dz_nir_dif(i,j)   = atm%dswd_dz_ir_cld(i,ja)  

        cmn%t2m(i,j,i_surf_macro_ocn) = atm%t2(i,ja,i_surf_macro_ocn)
        cmn%t2m(i,j,i_surf_macro_sic) = atm%t2(i,ja,i_surf_macro_sic)
        do n=1,nsurf_lnd
          if (n.ne.i_surf_lnd_ice .and. n.ne.i_surf_lnd_lake) then
            cmn%t2m(i,j,i_surf_lnd(n)) = atm%t2(i,ja,i_surf_macro_lnd)
          endif
        enddo
        cmn%t2m(i,j,i_surf_ice) = atm%t2(i,ja,i_surf_macro_ice)
        cmn%t2m(i,j,i_surf_lake) = atm%t2(i,ja,i_surf_macro_lake)
        cmn%q2m(i,j,i_surf_macro_ocn) = atm%q2(i,ja,i_surf_macro_ocn)
        cmn%q2m(i,j,i_surf_macro_sic) = atm%q2(i,ja,i_surf_macro_sic)
        do n=1,nsurf_lnd
          if (n.ne.i_surf_lnd_ice .and. n.ne.i_surf_lnd_lake) then
            cmn%q2m(i,j,i_surf_lnd(n)) = atm%q2(i,ja,i_surf_macro_lnd)
          endif
        enddo
        cmn%q2m(i,j,i_surf_ice) = atm%q2(i,ja,i_surf_macro_ice)
        cmn%q2m(i,j,i_surf_lake) = atm%q2(i,ja,i_surf_macro_lake)

        cmn%wind(i,j,i_surf_macro_ocn) = atm%wind(i,ja,i_surf_macro_ocn)
        cmn%wind(i,j,i_surf_macro_sic) = atm%wind(i,ja,i_surf_macro_sic)
        do n=1,nsurf_lnd
          if (n.ne.i_surf_lnd_ice .and. n.ne.i_surf_lnd_lake) then
            cmn%wind(i,j,i_surf_lnd(n)) = atm%wind(i,ja,i_surf_macro_lnd)
          endif
        enddo
        cmn%wind(i,j,i_surf_ice) = atm%wind(i,ja,i_surf_macro_ice)
        cmn%wind(i,j,i_surf_lake) = atm%wind(i,ja,i_surf_macro_lake)
        cmn%usur(i,j)   = sum(atm%us(i,ja,:)*atm%frst(i,ja,:))
        cmn%vsur(i,j)   = sum(atm%vs(i,ja,:)*atm%frst(i,ja,:))
        cmn%u700(i,j)   = atm%u3(i,ja,k700)
        cmn%v700(i,j)   = atm%v3(i,ja,k700)

        cmn%taux(i,j,:) = atm%taux(i,ja,:)
        cmn%tauy(i,j,:) = atm%tauy(i,ja,:)

        if (flag_dust) then
          cmn%dust_dep(i,j) = atm%dust_dep(i,ja)
        endif

        if (time_soy_atm) then         
          cmn%t2m_mon(i,j,:) = 0._wp
          cmn%t2m_mon_lnd(i,j,:) = 0._wp
        endif
        if (doy.gt.0) then
          cmn%t2m_mon(i,j,mon) = cmn%t2m_mon(i,j,mon) &
            + sum(cmn%t2m(i,j,:)*cmn%f_stp(i,j,:)) / nstep_mon_atm
          if (atm%frst(i,ja,i_surf_macro_lnd).gt.0._wp) then
            cmn%t2m_mon_lnd(i,j,mon) = cmn%t2m_mon_lnd(i,j,mon) &
              + sum(cmn%t2m(i,j,:)*cmn%f_stp(i,j,:), mask=(st2ast==i_surf_macro_lnd))/atm%frst(i,ja,i_surf_macro_lnd) / nstep_mon_atm
          endif
        endif
        if (time_eoy_atm) then
          cmn%t2m_min_mon(i,j) = minval(cmn%t2m_mon_lnd(i,j,:))
          cmn%t2m_glob_ann = sum( sum(cmn%t2m_mon(:,:,:),3)/nmon_year * area ) / sum(area)
        endif

      enddo
    enddo
    !$omp end parallel do


   return

  end subroutine atm_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ o c n
  ! Purpose    :  from common grid to ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_ocn(cmn,ocn)

    implicit none

    type(cmn_class) :: cmn
    type(ocn_class) :: ocn

    integer :: i, j, ip1, jp1

    real(wp) :: avg
    real(wp), dimension(:,:), allocatable :: taux, tauy


    ! update masks and sea level at beginning of year
    if (time_soy_ocn) then
      ocn%grid%ocn_vol_tot_real = cmn%ocn_vol_tot
      ocn%A_bering = cmn%A_bering
      ! save old ocean fraction and update
      ocn%f_ocn_old = ocn%f_ocn
      ocn%grid%ocn_area_old = ocn%grid%ocn_area
      ocn%grid%ocn_vol_old  = ocn%grid%ocn_vol
      ocn%z_ocn_max = cmn%z_ocn_max
      ocn%f_ocn     = cmn%f_ocn
      ocn%f_ocn2    = cmn%f_ocn2
      ocn%mask_coast= cmn%mask_coast2
      ocn%grid%mask_ocn = cmn%mask_ocn
      ocn%cfc11_atm = cmn%cfc11
      ocn%cfc12_atm = cmn%cfc12
      ocn%q_geo = cmn%q_geo
    endif

    ocn%buoy_sic_NA(:) = cmn%buoy_sic_NA(:)

    ! weighted average of wind stress over sea ice-free ocean and sea ice stress over sea ice covered fraction
    allocate(taux(ni,nj))
    allocate(tauy(ni,nj))

    if (flag_atm .and. atm_fix_tau) then

      taux = cmn%taux_dat
      tauy = cmn%tauy_dat

    else

      do j=1,nj
        do i=1,ni
          if (cmn%f_ocn(i,j).gt.0._wp) then
            taux(i,j) = cmn%f_ocn2(i,j)/cmn%f_ocn(i,j) * (cmn%f_sic(i,j)*(cmn%tauxo(i,j)) + (1._wp-cmn%f_sic(i,j))*cmn%taux(i,j,i_surf_macro_ocn))
            tauy(i,j) = cmn%f_ocn2(i,j)/cmn%f_ocn(i,j) * (cmn%f_sic(i,j)*(cmn%tauyo(i,j)) + (1._wp-cmn%f_sic(i,j))*cmn%tauy(i,j,i_surf_macro_ocn))
          else
            taux(i,j) = 0._wp
            tauy(i,j) = 0._wp
          endif
        enddo
      enddo

    endif

    ! interpolate stresses on u- and v-grid
    do j=1,nj
      do i=1,ni
        ip1 = i+1
        if (ip1.eq.ni+1) ip1=1
        jp1 = j+1
        if (jp1.eq.nj+1) jp1=nj
        ocn%stressxu(i,j) = 0.5_wp*(taux(i,j)+taux(ip1,j))  ! zonal wind stress on u-grid
        ocn%stressyv(i,j) = 0.5_wp*(tauy(i,j)+tauy(i,jp1))  ! meridional wind stress on v-grid
        ocn%stressxv(i,j) = 0.5_wp*(taux(i,j)+taux(i,jp1))  ! zonal wind stress on v-grid
        ocn%stressyu(i,j) = 0.5_wp*(tauy(i,j)+tauy(ip1,j))  ! meridional wind stress on u-grid
!        if (lat(j)>-20._wp .and. lat(j)<20._wp) then
!          ocn%stressxu(i,j) = min(0.1_wp,ocn%stressxu(i,j)) 
!          ocn%stressyv(i,j) = min(0.05_wp,ocn%stressyv(i,j)) 
!          ocn%stressxv(i,j) = min(0.1_wp,ocn%stressxv(i,j)) 
!          ocn%stressyu(i,j) = min(0.05_wp,ocn%stressyu(i,j)) 
!          ocn%stressxu(i,j) = max(-0.1_wp,ocn%stressxu(i,j)) 
!          ocn%stressxv(i,j) = max(-0.1_wp,ocn%stressxv(i,j)) 
!        endif
      enddo
    enddo
    deallocate(taux,tauy)

    ocn%wind(:,:) = cmn%wind(:,:,i_surf_macro_ocn)

    where (cmn%mask_ocn.eq.1)
      ocn%f_sic = cmn%f_sic
      ocn%slp   = cmn%slp
    endwhere

    if (ocn_restore_sal) then
      ! compute restoring salinity flux
      do j=1,nj
        do i=1,ni
          if (cmn%mask_ocn(i,j).eq.1) then
            ! restore sea surface salinity to climatology, 1 psu = 1 g/kg
            ocn%p_e_sic(i,j)  = -rho_w*50._wp * 1.e-3_wp*(cmn%sss_dat(i,j) - (ocn%ts(i,j,ubound(ocn%ts,3),2))) / (tau_sss*sec_day) ! kg/m2/s
          endif
        enddo
      enddo
      ocn%runoff  = 0._wp
      ocn%calving = 0._wp
      ocn%bmelt_grd = 0._wp
      ocn%bmelt_flt = 0._wp
    else
      where (cmn%mask_ocn.eq.1)
        ! pass P-E, runoff and calving separately
        ! scale from f_ocn2 (for which surface fluxes are computed) to f_ocn (for which the fluxes are applied)
        ocn%p_e_sic   = cmn%p_e_sic_ocn * ocn%f_ocn2/ocn%f_ocn 
        ocn%fw_brines = cmn%fw_brines   * ocn%f_ocn2/ocn%f_ocn 
        ocn%runoff    = cmn%runoff_o
        ocn%runoff_veg    = cmn%runoff_veg_o
        ocn%runoff_ice    = cmn%runoff_ice_o
        ocn%runoff_lake   = cmn%runoff_lake_o
        ocn%melt_ice  = cmn%melt_ice_o
        ocn%calving   = cmn%calving_o
        ocn%bmelt_grd = cmn%bmelt_grd_o
        ocn%bmelt_flt = cmn%bmelt_flt_o
        ocn%bmelt     = ocn%bmelt_grd+ocn%bmelt_flt
      endwhere
    endif

    if (ocn_restore_temp) then
      ! compute restoring heat flux
      do j=1,nj
        do i=1,ni
          if (cmn%mask_ocn(i,j).eq.1) then
            ! restore sea surface temperature to climatology
            ocn%flx(i,j) = rho_w*cap_w*50._wp * (cmn%sst_dat(i,j) - ocn%ts(i,j,ubound(ocn%ts,3),1)) / (tau_sst*sec_day)  ! W/m2
          endif
        enddo
      enddo
    else
      where (cmn%mask_ocn.eq.1)
        ! net heat flux into the ocean
        ocn%flx = cmn%flx_ocn * ocn%f_ocn2/ocn%f_ocn ! scale from f_ocn2 (for which surface fluxes are computed) to f_ocn (for which the fluxes are applied)
      endwhere
    endif

    !------------------------------------------------------------
    ! save daily input fields for offline simulation, if required
    if (time_call_daily_input_save) then
      ! average over nyear_avg_offline years to remove possible noise
      avg = real(nday_year,wp)/real(nstep_year_ocn,wp) / real(nyear_avg_offline,wp)
      ocn%daily_input_save%stressxu(:,:,doy)    = ocn%daily_input_save%stressxu(:,:,doy)    + ocn%stressxu(:,:)     * avg 
      ocn%daily_input_save%stressyv(:,:,doy)    = ocn%daily_input_save%stressyv(:,:,doy)    + ocn%stressyv(:,:)     * avg 
      ocn%daily_input_save%stressxv(:,:,doy)    = ocn%daily_input_save%stressxv(:,:,doy)    + ocn%stressxv(:,:)     * avg 
      ocn%daily_input_save%stressyu(:,:,doy)    = ocn%daily_input_save%stressyu(:,:,doy)    + ocn%stressyu(:,:)     * avg 
      ocn%daily_input_save%wind(:,:,doy)        = ocn%daily_input_save%wind(:,:,doy)        + ocn%wind(:,:)         * avg
      ocn%daily_input_save%f_sic(:,:,doy)       = ocn%daily_input_save%f_sic(:,:,doy)       + ocn%f_sic(:,:)        * avg
      ocn%daily_input_save%slp(:,:,doy)         = ocn%daily_input_save%slp(:,:,doy)         + ocn%slp(:,:)          * avg
      ocn%daily_input_save%p_e_sic(:,:,doy)     = ocn%daily_input_save%p_e_sic(:,:,doy)     + ocn%p_e_sic(:,:)      * avg  
      ocn%daily_input_save%fw_brines(:,:,doy)   = ocn%daily_input_save%fw_brines(:,:,doy)   + ocn%fw_brines(:,:)    * avg
      ocn%daily_input_save%runoff(:,:,doy)      = ocn%daily_input_save%runoff(:,:,doy)      + ocn%runoff(:,:)       * avg
      ocn%daily_input_save%runoff_veg(:,:,doy)  = ocn%daily_input_save%runoff_veg(:,:,doy)  + ocn%runoff_veg(:,:)   * avg 
      ocn%daily_input_save%runoff_ice(:,:,doy)  = ocn%daily_input_save%runoff_ice(:,:,doy)  + ocn%runoff_ice(:,:)   * avg 
      ocn%daily_input_save%runoff_lake(:,:,doy) = ocn%daily_input_save%runoff_lake(:,:,doy) + ocn%runoff_lake(:,:)  * avg 
      ocn%daily_input_save%calving(:,:,doy)     = ocn%daily_input_save%calving(:,:,doy)     + ocn%calving(:,:)      * avg
      ocn%daily_input_save%bmelt_grd(:,:,doy)   = ocn%daily_input_save%bmelt_grd(:,:,doy)   + ocn%bmelt_grd(:,:)    * avg
      ocn%daily_input_save%bmelt_flt(:,:,doy)   = ocn%daily_input_save%bmelt_flt(:,:,doy)   + ocn%bmelt_flt(:,:)    * avg
      ocn%daily_input_save%bmelt(:,:,doy)       = ocn%daily_input_save%bmelt(:,:,doy)       + ocn%bmelt(:,:)        * avg
      ocn%daily_input_save%flx(:,:,doy)         = ocn%daily_input_save%flx(:,:,doy)         + ocn%flx(:,:)          * avg      
    endif
    !------------------------------------------------------------
    ! use saved daily input fields now
    if (time_use_daily_input_save) then
      ocn%stressxu(:,:)    = ocn%daily_input_save%stressxu(:,:,doy)    
      ocn%stressyv(:,:)    = ocn%daily_input_save%stressyv(:,:,doy)    
      ocn%stressxv(:,:)    = ocn%daily_input_save%stressxv(:,:,doy)    
      ocn%stressyu(:,:)    = ocn%daily_input_save%stressyu(:,:,doy)    
      ocn%wind(:,:)        = ocn%daily_input_save%wind(:,:,doy)        
      ocn%f_sic(:,:)       = ocn%daily_input_save%f_sic(:,:,doy)       
      ocn%slp(:,:)         = ocn%daily_input_save%slp(:,:,doy)         
      ocn%p_e_sic(:,:)     = ocn%daily_input_save%p_e_sic(:,:,doy)      
      ocn%fw_brines(:,:)   = ocn%daily_input_save%fw_brines(:,:,doy)   
      ocn%runoff(:,:)      = ocn%daily_input_save%runoff(:,:,doy)      
      ocn%runoff_veg(:,:)  = ocn%daily_input_save%runoff_veg(:,:,doy)  
      ocn%runoff_ice(:,:)  = ocn%daily_input_save%runoff_ice(:,:,doy)  
      ocn%runoff_lake(:,:) = ocn%daily_input_save%runoff_lake(:,:,doy) 
      ocn%calving(:,:)     = ocn%daily_input_save%calving(:,:,doy)     
      ocn%bmelt_grd(:,:)   = ocn%daily_input_save%bmelt_grd(:,:,doy)   
      ocn%bmelt_flt(:,:)   = ocn%daily_input_save%bmelt_flt(:,:,doy)   
      ocn%bmelt(:,:)       = ocn%daily_input_save%bmelt(:,:,doy)       
      ocn%flx(:,:)         = ocn%daily_input_save%flx(:,:,doy)              
    endif

    !------------------------------------------------------------
    ! save daily input fields for simulations with (selected) fixed boundary conditions 
    if (l_ocn_input_fix .and. i_ocn_input_fix.eq.1 .and. year.le.30) then
      ! average over 30 years to remove possible noise
      avg = real(nday_year,wp)/real(nstep_year_ocn,wp) / 30._wp
      ocn%daily_input_save%stressxu(:,:,doy)    = ocn%daily_input_save%stressxu(:,:,doy)    + ocn%stressxu(:,:)     * avg 
      ocn%daily_input_save%stressyv(:,:,doy)    = ocn%daily_input_save%stressyv(:,:,doy)    + ocn%stressyv(:,:)     * avg 
      ocn%daily_input_save%stressxv(:,:,doy)    = ocn%daily_input_save%stressxv(:,:,doy)    + ocn%stressxv(:,:)     * avg 
      ocn%daily_input_save%stressyu(:,:,doy)    = ocn%daily_input_save%stressyu(:,:,doy)    + ocn%stressyu(:,:)     * avg 
      ocn%daily_input_save%wind(:,:,doy)        = ocn%daily_input_save%wind(:,:,doy)        + ocn%wind(:,:)         * avg
      ocn%daily_input_save%f_sic(:,:,doy)       = ocn%daily_input_save%f_sic(:,:,doy)       + ocn%f_sic(:,:)        * avg
      ocn%daily_input_save%slp(:,:,doy)         = ocn%daily_input_save%slp(:,:,doy)         + ocn%slp(:,:)          * avg
      ocn%daily_input_save%p_e_sic(:,:,doy)     = ocn%daily_input_save%p_e_sic(:,:,doy)     + ocn%p_e_sic(:,:)      * avg  
      ocn%daily_input_save%fw_brines(:,:,doy)   = ocn%daily_input_save%fw_brines(:,:,doy)   + ocn%fw_brines(:,:)    * avg
      ocn%daily_input_save%runoff(:,:,doy)      = ocn%daily_input_save%runoff(:,:,doy)      + ocn%runoff(:,:)       * avg
      ocn%daily_input_save%runoff_veg(:,:,doy)  = ocn%daily_input_save%runoff_veg(:,:,doy)  + ocn%runoff_veg(:,:)   * avg 
      ocn%daily_input_save%runoff_ice(:,:,doy)  = ocn%daily_input_save%runoff_ice(:,:,doy)  + ocn%runoff_ice(:,:)   * avg 
      ocn%daily_input_save%runoff_lake(:,:,doy) = ocn%daily_input_save%runoff_lake(:,:,doy) + ocn%runoff_lake(:,:)  * avg 
      ocn%daily_input_save%calving(:,:,doy)     = ocn%daily_input_save%calving(:,:,doy)     + ocn%calving(:,:)      * avg
      ocn%daily_input_save%bmelt_grd(:,:,doy)   = ocn%daily_input_save%bmelt_grd(:,:,doy)   + ocn%bmelt_grd(:,:)    * avg
      ocn%daily_input_save%bmelt_flt(:,:,doy)   = ocn%daily_input_save%bmelt_flt(:,:,doy)   + ocn%bmelt_flt(:,:)    * avg
      ocn%daily_input_save%bmelt(:,:,doy)       = ocn%daily_input_save%bmelt(:,:,doy)       + ocn%bmelt(:,:)        * avg
      ocn%daily_input_save%flx(:,:,doy)         = ocn%daily_input_save%flx(:,:,doy)         + ocn%flx(:,:)          * avg      
    endif
    !------------------------------------------------------------
    ! use saved daily input fields now
    if (l_ocn_input_fix .and. (i_ocn_input_fix.eq.2 .or. (i_ocn_input_fix.eq.1 .and. year.gt.30))) then
      if (l_ocn_fix_wind) then
        ocn%stressxu(:,:)    = ocn%daily_input_save%stressxu(:,:,doy)    
        ocn%stressyv(:,:)    = ocn%daily_input_save%stressyv(:,:,doy)    
        ocn%stressxv(:,:)    = ocn%daily_input_save%stressxv(:,:,doy)    
        ocn%stressyu(:,:)    = ocn%daily_input_save%stressyu(:,:,doy)    
        ocn%wind(:,:)        = ocn%daily_input_save%wind(:,:,doy)        
      endif
      !ocn%f_sic(:,:)       = ocn%daily_input_save%f_sic(:,:,doy)       
      !ocn%slp(:,:)         = ocn%daily_input_save%slp(:,:,doy)         
      if (l_ocn_fix_fw) then
        ocn%p_e_sic(:,:)     = ocn%daily_input_save%p_e_sic(:,:,doy)      
        ocn%fw_brines(:,:)   = ocn%daily_input_save%fw_brines(:,:,doy)   
        ocn%runoff(:,:)      = ocn%daily_input_save%runoff(:,:,doy)      
        ocn%runoff_veg(:,:)  = ocn%daily_input_save%runoff_veg(:,:,doy)  
        ocn%runoff_ice(:,:)  = ocn%daily_input_save%runoff_ice(:,:,doy)  
        ocn%runoff_lake(:,:) = ocn%daily_input_save%runoff_lake(:,:,doy) 
        ocn%calving(:,:)     = ocn%daily_input_save%calving(:,:,doy)     
        ocn%bmelt_grd(:,:)   = ocn%daily_input_save%bmelt_grd(:,:,doy)   
        ocn%bmelt_flt(:,:)   = ocn%daily_input_save%bmelt_flt(:,:,doy)   
        ocn%bmelt(:,:)       = ocn%daily_input_save%bmelt(:,:,doy)       
      endif
      if (l_ocn_fix_flx) then
        ocn%flx(:,:)         = ocn%daily_input_save%flx(:,:,doy)              
      endif
    endif


   return

  end subroutine cmn_to_ocn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  o c n _ t o _ c m n
  ! Purpose    :  from ocean to common grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_to_cmn(ocn,cmn)

    implicit none

    type(cmn_class) :: cmn
    type(ocn_class) :: ocn

    integer :: i, j, k, ii, jj, iii, jjj, n, nk, kshelf, k_bmb


    ! find index of top 1 km ocean layers for basal melt
    if (flag_bmb) then
      k_bmb = minloc(abs(ocn%grid%zro+1000._wp),1)
      nk = ocn%grid%nk - k_bmb + 1
      if (.not.allocated(cmn%z_ocn_bmb)) allocate(cmn%z_ocn_bmb(nk))
      if (.not.allocated(cmn%mask_ocn_bmb)) allocate(cmn%mask_ocn_bmb(ni,nj,nk))
      if (.not.allocated(cmn%t_ocn_bmb)) allocate(cmn%t_ocn_bmb(ni,nj,nk))
      if (.not.allocated(cmn%s_ocn_bmb)) allocate(cmn%s_ocn_bmb(ni,nj,nk))
      cmn%z_ocn_bmb = -ocn%grid%zro(ocn%grid%nk:k_bmb:-1)
    endif

    !$omp parallel do collapse(2) private(i,j,kshelf,ii,jj,iii,jjj,n)
    do j=1,nj
      do i=1,ni   
        cmn%uo1(i,j) = ocn%u(1,i,j,ubound(ocn%ts,3))  ! m/s, zonal surface ocean velocity on u-grid
        cmn%vo1(i,j) = ocn%u(2,i,j,ubound(ocn%ts,3))  ! m/s, meridional surface ocean velocity on v-grid

        cmn%ssh(i,j) = ocn%ssh(i,j)

        if (cmn%mask_ocn(i,j).eq.1) then
          cmn%sst(i,j) = ocn%ts(i,j,ubound(ocn%ts,3),1)   ! degC, top ocean temperature
          cmn%sss(i,j) = ocn%ts(i,j,ubound(ocn%ts,3),2)   ! psu, top ocean salinity
          kshelf = max(ocn%grid%k1(i,j),k1_shelf)  
          cmn%t_shelf(i,j) = ocn%ts(i,j,kshelf,1) + T0  ! K
        endif
        if (flag_bmb) then
          ! temperature and salinity for ice shelf basal melt
          cmn%mask_ocn_bmb(i,j,1:nk) = ocn%grid%mask_c(i,j,ocn%grid%nk:k_bmb:-1)
          cmn%t_ocn_bmb(i,j,1:nk) = ocn%ts(i,j,ocn%grid%nk:k_bmb:-1,1)  ! degC
          cmn%s_ocn_bmb(i,j,1:nk) = ocn%ts(i,j,ocn%grid%nk:k_bmb:-1,2)  ! psu
        endif
      enddo
    enddo
    !$omp end parallel do

   return

  end subroutine ocn_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ s i c
  ! Purpose    :  from common grid to sea ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_sic(cmn,sic)

    implicit none

    type(cmn_class) :: cmn
    type(sic_class) :: sic

    integer :: i, j, ip1


    ! update ocean fraction at beginning of year
    if (time_soy_sic) then
      sic%f_ocn = cmn%f_ocn2
      sic%coszm = cmn%coszm
    endif

    !$omp parallel do collapse(2) private(i,j,ip1)
    do j=1,nj
      do i=1,ni
        ip1 = i+1
        if (ip1.eq.ni+1) ip1=1
        if (cmn%mask_ocn(i,j).eq.1) then
          sic%rain(i,j)     = sic%f_sic(i,j)*cmn%rain(i,j,i_surf_macro_sic) + (1._wp-sic%f_sic(i,j))*cmn%rain(i,j,i_surf_macro_ocn)
          sic%snow(i,j)     = sic%f_sic(i,j)*cmn%snow(i,j,i_surf_macro_sic) + (1._wp-sic%f_sic(i,j))*cmn%snow(i,j,i_surf_macro_ocn)
          sic%pressure(i,j) = cmn%ps(i,j,i_surf_macro_ocn)
          sic%flx_lwd_sic(i,j)  = cmn%lwd(i,j,i_surf_macro_sic)
          sic%flx_lwd_ocn(i,j)  = cmn%lwd(i,j,i_surf_macro_ocn)
          sic%wind(i,j)      = sic%f_sic(i,j)*cmn%wind(i,j,i_surf_macro_sic) + (1._wp-sic%f_sic(i,j))*cmn%wind(i,j,i_surf_macro_ocn)
          sic%t_air_ocn(i,j) = cmn%t2m(i,j,i_surf_macro_ocn)
          sic%t_air_sic(i,j) = cmn%t2m(i,j,i_surf_macro_sic)
          sic%q_air_ocn(i,j) = cmn%q2m(i,j,i_surf_macro_ocn)
          sic%q_air_sic(i,j) = cmn%q2m(i,j,i_surf_macro_sic)
          sic%sst(i,j)       = cmn%sst(i,j) + T0  ! K
          sic%sss(i,j)       = cmn%sss(i,j)
          if (flag_atm) then
            sic%flx_swnet_sic(i,j) = cmn%swnet(i,j,i_surf_macro_sic)
            sic%flx_swnet_ocn(i,j) = cmn%swnet(i,j,i_surf_macro_ocn)
          else
            sic%flx_swnet_sic(i,j) = cmn%swd(i,j) * (1._wp-sic%albedo_sic(i,j))
            sic%flx_swnet_ocn(i,j) = cmn%swd(i,j) * (1._wp-sic%albedo_ocn(i,j))
          endif
        endif
        sic%ssh(i,j) = cmn%ssh(i,j)
        sic%uo(i,j) = cmn%uo1(i,j)
        sic%vo(i,j) = cmn%vo1(i,j)
        ! zonal wind stress on u-grid
        sic%tauxa(i,j)     = 0.5_wp * (cmn%taux(i,j,i_surf_macro_sic) + cmn%taux(ip1,j,i_surf_macro_sic))
        ! meridional wind stress on v-grid
        if (j.lt.nj) then
          sic%tauya(i,j)   = 0.5_wp * (cmn%tauy(i,j,i_surf_macro_sic) + cmn%tauy(i,j+1,i_surf_macro_sic))
        else
          sic%tauya(i,j)   = 0._wp
        endif
        ! dust deposition
        sic%dust_dep(i,j) = cmn%dust_dep(i,j)
      enddo
    enddo
    !$omp end parallel do


   return

  end subroutine cmn_to_sic


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s i c _ t o _ c m n
  ! Purpose    :  from sea ice to common grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_to_cmn(sic,cmn)

    implicit none

    type(cmn_class) :: cmn
    type(sic_class) :: sic
    
    integer :: i, j, im1


    cmn%buoy_sic_NA(:) = sic%buoy_sic_NA(:)

    !$omp parallel do collapse(2) private(i,j,im1)
    do j=1,nj
      do i=1,ni
        cmn%f_stp(i,j,i_surf_macro_sic)       = sic%f_sic(i,j) * cmn%f_ocn2(i,j)
        cmn%f_stp(i,j,i_surf_macro_ocn)       = (1._wp-sic%f_sic(i,j)) * cmn%f_ocn2(i,j)
        if (cmn%mask_ocn(i,j).eq.1) then
          ! if atmosphere model active
          if (flag_atm) then
            cmn%alb_vis_dir(i,j,i_surf_macro_sic) = sic%alb_sic_vis_dir(i,j)
            cmn%alb_vis_dif(i,j,i_surf_macro_sic) = sic%alb_sic_vis_dif(i,j)
            cmn%alb_nir_dir(i,j,i_surf_macro_sic) = sic%alb_sic_nir_dir(i,j)
            cmn%alb_nir_dif(i,j,i_surf_macro_sic) = sic%alb_sic_nir_dif(i,j)
            cmn%alb_vis_dir(i,j,i_surf_macro_ocn) = sic%alb_ocn_vis_dir(i,j)
            cmn%alb_vis_dif(i,j,i_surf_macro_ocn) = sic%alb_ocn_vis_dif(i,j)
            cmn%alb_nir_dir(i,j,i_surf_macro_ocn) = sic%alb_ocn_nir_dir(i,j)
            cmn%alb_nir_dif(i,j,i_surf_macro_ocn) = sic%alb_ocn_nir_dif(i,j)
            cmn%lwu(i,j,i_surf_macro_sic)         = sic%flx_lwu_sic(i,j)
            cmn%lwu(i,j,i_surf_macro_ocn)         = sic%flx_lwu_ocn(i,j)
            cmn%sh(i,j,i_surf_macro_sic)          = sic%flx_sh_sic(i,j)
            cmn%sh(i,j,i_surf_macro_ocn)          = sic%flx_sh_ocn(i,j)
            cmn%lh(i,j,i_surf_macro_sic)          = sic%flx_lh_sic(i,j)
            cmn%lh(i,j,i_surf_macro_ocn)          = sic%flx_lh_ocn(i,j)
            cmn%evp(i,j,i_surf_macro_sic)         = sic%evp_sic(i,j)
            cmn%evp(i,j,i_surf_macro_ocn)         = sic%evp_ocn(i,j)
            cmn%t_skin(i,j,i_surf_macro_sic)      = sic%t_skin_sic(i,j)
            cmn%t_skin(i,j,i_surf_macro_ocn)      = sic%t_skin_ocn(i,j)
            cmn%z0m(i,j,i_surf_macro_sic)         = sic%rough_m_sic(i,j)
            cmn%z0m(i,j,i_surf_macro_ocn)         = sic%rough_m_ocn(i,j)
          endif
          ! if ocean model active
          if (flag_ocn) then
            cmn%flx_ocn(i,j) = sic%flx_ocn(i,j)
            cmn%p_e_sic_ocn(i,j) = sic%fw_ocn(i,j)
            cmn%fw_brines(i,j) = sic%fw_brines(i,j)
            cmn%f_sic(i,j)   = sic%f_sic(i,j)
          endif
        endif
        if (flag_ocn) then
          im1 = i-1
          if (im1.eq.0) im1=ni
          ! intrerpolate zonal sea ice stress on ocean from u-points to c-points
          cmn%tauxo(i,j)   = 0.5_wp * (sic%tauxo(i,j) + sic%tauxo(im1,j))   ! on c-points
          ! intrerpolate meridional sea ice stress on ocean from v-points to c-points
          if (j.gt.1) then
            cmn%tauyo(i,j)   = 0.5_wp * (sic%tauyo(i,j) + sic%tauyo(i,j-1))   ! on c-points
          else
            cmn%tauyo(i,j)   = 0._wp 
          endif
        endif
      enddo
    enddo
    !$omp end parallel do


   return

  end subroutine sic_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ l n d
  ! Purpose    :  from common grid to land
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_lnd(cmn,lnd)

    implicit none

    type(cmn_class) :: cmn
    type(lnd_class) :: lnd

    integer :: i, j, ii, jj, iii, jjj, n, k


    if (time_soy_lnd) then
      ! get at start of year

      lnd%l0d%co2         = cmn%co2
      lnd%l0d%c13_c12_atm = cmn%c13_c12_atm
      lnd%l0d%c14_c_atm   = cmn%c14_c_atm

      ! save fractions before updating
      lnd%l2d%f_veg_old   = lnd%l2d%f_veg
      lnd%l2d%f_ice_old   = lnd%l2d%f_ice
      lnd%l2d%f_ice_grd_old   = lnd%l2d%f_ice_grd
      lnd%l2d%f_shelf_old = lnd%l2d%f_shelf
      lnd%l2d%f_lake_old  = lnd%l2d%f_lake

      ! land fraction 
      lnd%l2d%f_land0 = cmn%f_lnd0
      lnd%l2d%f_land  = cmn%f_lnd
      ! ice fraction
      lnd%l2d%f_ice   = cmn%f_ice
      ! ice fraction, grounded ice only 
      lnd%l2d%f_ice_grd   = cmn%f_ice_grd
      ! for land model, grounded ice fraction can not be larger than land fraction 
      where (lnd%l2d%f_ice_grd.gt.lnd%l2d%f_land) 
        lnd%l2d%f_ice_grd = lnd%l2d%f_land
      endwhere
      ! sum of ice fraction in neighboring cells
      do j=1,nj
        do i=1,ni
          lnd%l2d(i,j)%f_ice_nbr = 0._wp
          if(cmn%mask_lnd(i,j).eq.1) then
            do ii=i-1,i+1
              do jj=j-1,j+1
                iii = ii
                jjj = jj
                if (iii.eq.ni+1) iii = 1
                if (iii.eq.0) iii = ni
                jjj = min(jjj+1,nj)
                jjj = max(jjj-1,1)
                lnd%l2d(i,j)%f_ice_nbr = lnd%l2d(i,j)%f_ice_nbr + lnd%l2d(iii,jjj)%f_ice
              enddo
            enddo
          endif
        enddo
      enddo
      ! water-covered shelf fraction
      lnd%l2d%f_shelf = cmn%f_ocn 
      ! lake fraction
      lnd%l2d%f_lake  = cmn%f_lake

      ! land cover state
      lnd%l2d%f_crop    = cmn%f_crop
      lnd%l2d%f_pasture = cmn%f_pasture

      ! disturbance rate
      do j=1,nj
        do i=1,ni
          do n=1,5
            lnd%l2d(i,j)%disturbance(n) = cmn%disturbance(n,i,j)
          enddo
        enddo
      enddo

      ! potentially vegetated fraction
      where (lnd%l2d%mask_lnd.eq.1) 
        lnd%l2d%f_veg = lnd%l2d%f_land - lnd%l2d%f_ice_grd - lnd%l2d%f_lake
      elsewhere
        lnd%l2d%f_veg = 0._wp
      endwhere
      ! f_veg can be negative for numerical reasons
      where (lnd%l2d%f_veg.lt.0._wp)
        lnd%l2d%f_veg = 0._wp
      endwhere
      where (lnd%l2d%f_veg.lt.1.e-10_wp)
        lnd%l2d%f_veg = 0._wp
      endwhere

      ! standard deviation of topography, ice-free land only
      lnd%l2d%z_veg_std = cmn%z_veg_std

      ! mean grid cell elevation
      lnd%l2d%z_veg     = cmn%z_veg
      ! min and max grid cell elevations
      lnd%l2d%z_veg_min = cmn%z_veg_min
      lnd%l2d%z_veg_max = cmn%z_veg_max
      ! average over 15x15deg (3x3 cells)
      do j=2,nj-1
        do i=2,ni-1
          lnd%l2d(i,j)%z_veg_min = minval(cmn%z_veg_min(i-1:i+1,j-1:j+1))
          lnd%l2d(i,j)%z_veg_max = maxval(cmn%z_veg_max(i-1:i+1,j-1:j+1))
        enddo
      enddo

      ! carbonate fraction for lithology used for weathering
      lnd%l2d%f_carb = cmn%f_carb
      
      ! get cosine of zenith angle and daylength at start of year
      do k=1,nday_year
        do j=1,nj
          do i=1,ni
            lnd%l2d(i,j)%coszm(k)     = cmn%coszm(k,j)
            lnd%l2d(i,j)%daylength(k) = cmn%daylength(k,j)
          enddo
        enddo
      enddo

      ! lake depth
      if (flag_lakes) then
        do j=1,nj
          do i=1,ni
            if (cmn%f_lake(i,j).gt.0._wp) then
              ! compute average depth of lake(s) in grid cell
              lnd%l2d(i,j)%h_lake = 0._wp
              do n=1,cmn%n_lakes
                lnd%l2d(i,j)%h_lake = lnd%l2d(i,j)%h_lake + cmn%lake(n)%depth*cmn%f_lake_n(n,i,j)/cmn%f_lake(i,j)
              enddo
            else
              ! no lakes in grid cell
              lnd%l2d(i,j)%h_lake = 0._wp  ! dummy value
            endif
          enddo
        enddo
      else  
        ! dynamic lakes not active, but lake depth still needed for lake thermodynamics
        lnd%l2d(:,:)%h_lake = 100._wp    
      endif

    endif

    if (time_eoy_lnd) then
      lnd%l2d%t2m_min_mon = cmn%t2m_min_mon
      lnd%l2d%t2m_ann_mean = sum(cmn%t2m_mon_lnd,3)/real(nmon_year,wp)
    endif

    !$omp parallel do collapse(2) private(i,j,n)
    do j=1,nj
      do i=1,ni
        if(cmn%mask_lnd(i,j).eq.1) then
          ! get current surface climate and fluxes down from the atmosphere
          do n=1,nsurf_lnd
            lnd%l2d(i,j)%pressure(n) = cmn%ps(i,j,i_surf_lnd(n))
            lnd%l2d(i,j)%t2m(n)   = cmn%t2m(i,j,i_surf_lnd(n))
            lnd%l2d(i,j)%q2m(n)   = cmn%q2m(i,j,i_surf_lnd(n))
            lnd%l2d(i,j)%tatm(n)  = cmn%t2m(i,j,i_surf_lnd(n))
            lnd%l2d(i,j)%qatm(n)  = cmn%q2m(i,j,i_surf_lnd(n))
            lnd%l2d(i,j)%wind(n)  = cmn%wind(i,j,i_surf_lnd(n))
            lnd%l2d(i,j)%rain(n)  = cmn%rain(i,j,i_surf_lnd(n))
            lnd%l2d(i,j)%snow(n)  = cmn%snow(i,j,i_surf_lnd(n))
            if (flag_atm) then
              lnd%l2d(i,j)%swnet(n)  = cmn%swnet(i,j,i_surf_lnd(n))
            else
              lnd%l2d(i,j)%swnet(n)  = cmn%swd(i,j)*(1._wp-lnd%l2d(i,j)%albedo(n)) ! net shortwave at the surface, approx (using old albedo)
            endif
            if (cmn%solarmin(doy,j).gt.0._wp) then
              lnd%l2d(i,j)%swnet_min(n)  = lnd%l2d(i,j)%swnet(n) * cmn%solarmin(doy,j)/cmn%solarm(doy,j) ! daily minimum shortwave at the surface, approx
            else
              lnd%l2d(i,j)%swnet_min(n) = 0._wp
            endif
            lnd%l2d(i,j)%lwdown(n)   = cmn%lwd(i,j,i_surf_lnd(n))       
          enddo
          ! dust deposition
          lnd%l2d(i,j)%dust_dep = cmn%dust_dep(i,j)
        endif
        ! shelf water temperature
        if (cmn%mask_ocn(i,j).eq.1) then
            lnd%l2d(i,j)%t_shelf(0) = cmn%t_shelf(i,j)
        endif
      enddo
    enddo
    !$omp end parallel do


   return

  end subroutine cmn_to_lnd


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  l n d _ t o _ c m n
  ! Purpose    :  from land to common grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_to_cmn(lnd,cmn)

    implicit none

    type(cmn_class), intent(inout) :: cmn
    type(lnd_class), intent(inout) :: lnd

    integer :: i, j, ii, jj, ir, jr, n, nbr, nk
    real(wp) :: lnd_cell_area, area_nbr, avg


    ! find index of top 1 km ocean layers for basal melt
    if (flag_bmb) then
      nk = size(lnd%z_lake)
      if (.not.allocated(cmn%z_lake_bmb)) allocate(cmn%z_lake_bmb(nk))
      if (.not.allocated(cmn%mask_lake_bmb)) allocate(cmn%mask_lake_bmb(ni,nj))
      if (.not.allocated(cmn%t_lake_bmb)) allocate(cmn%t_lake_bmb(ni,nj,nk))
      if (.not.allocated(cmn%s_lake_bmb)) allocate(cmn%s_lake_bmb(ni,nj,nk))
      cmn%z_lake_bmb = lnd%z_lake
    endif

    !$omp parallel do collapse(2) private(i,j,n,ir,jr)
    do j=1,nj
      do i=1,ni

        if (flag_atm) then
          cmn%delta_C_lnd_2d(i,j)   = lnd%l2d(i,j)%Cflx_atm_lnd ! kgC/s
          cmn%delta_C13_lnd_2d(i,j) = lnd%l2d(i,j)%C13flx_atm_lnd ! kgC13/s
          cmn%delta_C14_lnd_2d(i,j) = lnd%l2d(i,j)%C14flx_atm_lnd ! kgC14/s
          do n=1,nsurf_lnd
            cmn%f_stp(i,j,i_surf_lnd(n))       = lnd%l2d(i,j)%frac_surf(n)
            cmn%alb_vis_dir(i,j,i_surf_lnd(n)) = lnd%l2d(i,j)%alb_vis_dir(n)
            cmn%alb_vis_dif(i,j,i_surf_lnd(n)) = lnd%l2d(i,j)%alb_vis_dif(n)
            cmn%alb_nir_dir(i,j,i_surf_lnd(n)) = lnd%l2d(i,j)%alb_nir_dir(n)
            cmn%alb_nir_dif(i,j,i_surf_lnd(n)) = lnd%l2d(i,j)%alb_nir_dif(n)
            cmn%lwu(i,j,i_surf_lnd(n))     = lnd%l2d(i,j)%flx_lwu(n)
            cmn%sh(i,j,i_surf_lnd(n))      = lnd%l2d(i,j)%flx_sh(n)
            cmn%lh(i,j,i_surf_lnd(n))      = lnd%l2d(i,j)%flx_lh(n)
            cmn%evp(i,j,i_surf_lnd(n))         = lnd%l2d(i,j)%et(n)
            cmn%z0m(i,j,i_surf_lnd(n))         = lnd%l2d(i,j)%rough_m(n)
            cmn%t_skin(i,j,i_surf_lnd(n))      = lnd%l2d(i,j)%t_skin(n)
          enddo
          cmn%dust_emis(i,j) = lnd%l2d(i,j)%dust_emis
          ! overwrite ice sheet albedo with albedo computed in SEMI if required
          if (flag_smb .and. l_ice_albedo_semi) then
            if (cmn%mask_smb(i,j).eq.1) then
              cmn%alb_vis_dir(i,j,i_surf_ice) = cmn%alb_vis_dir_ice_semi(i,j,doy)
              cmn%alb_vis_dif(i,j,i_surf_ice) = cmn%alb_vis_dif_ice_semi(i,j,doy)
              cmn%alb_nir_dir(i,j,i_surf_ice) = cmn%alb_nir_dir_ice_semi(i,j,doy)
              cmn%alb_nir_dif(i,j,i_surf_ice) = cmn%alb_nir_dif_ice_semi(i,j,doy)
            endif
          endif
        endif

        ! bottom soil temperature for smb (then passed to ice sheet model for top boundary condition)
        if (lnd%l2d(i,j)%f_veg.gt.0._wp) then
          cmn%t_soil(i,j) = lnd%l2d(i,j)%t_soil(nl) 
        else
          cmn%t_soil(i,j) = lnd%l2d(i,j)%t_ice(nl) 
        endif

        ! runoff and 'calving' 
        cmn%runoff_veg(i,j)    = lnd%l2d(i,j)%runoff(is_veg)  * area(i,j)*(cmn%f_lnd(i,j)-cmn%f_lake(i,j)-cmn%f_ice_grd(i,j)) ! kg/m2/s -> kg/s
        cmn%calving_veg(i,j)   = lnd%l2d(i,j)%calving(is_veg) * area(i,j)*(cmn%f_lnd(i,j)-cmn%f_lake(i,j)-cmn%f_ice_grd(i,j)) ! kg/m2/s -> kg/s
        cmn%runoff_ice_l(i,j)  = lnd%l2d(i,j)%runoff(is_ice)  * area(i,j)*cmn%f_ice(i,j)  ! kg/m2/s -> kg/s
        cmn%calving_ice_l(i,j) = lnd%l2d(i,j)%calving(is_ice) * area(i,j)*cmn%f_ice(i,j)  ! kg/m2/s -> kg/s

        if (flag_bmb) then
          ! temperature and salinity for ice shelf basal melt
          if (cmn%f_lake(i,j).gt.0._wp) then
            cmn%mask_lake_bmb(i,j) = 1
          else
            cmn%mask_lake_bmb(i,j) = 0
          endif
          cmn%t_lake_bmb(i,j,1:nk) = lnd%l2d(i,j)%t_lake(1:nk)-T0  ! degC
          cmn%s_lake_bmb(i,j,1:nk) = 0._wp  ! psu, assume fresh lake
        endif

        ! lake ice fraction
        cmn%f_ice_lake(i,j) = lnd%l2d(i,j)%f_lake_ice

        ! lake P-E
        cmn%lake_p_e(i,j) = lnd%l2d(i,j)%lake_water_tendency   ! kg/m2/s

      enddo
    enddo
    !$omp end parallel do

    if (time_eoy_lnd) then
      ! transfer weathering runoff fluxes from land to ocean domain at end of year
      cmn%doc_export = 0._wp
      cmn%doc13_export = 0._wp
      cmn%doc14_export = 0._wp
      cmn%poc_export = 0._wp
      cmn%poc13_export = 0._wp
      cmn%poc14_export = 0._wp
      cmn%weath_carb = 0._wp
      cmn%weath_sil = 0._wp
      cmn%weath13_carb = 0._wp
      cmn%weath13_sil = 0._wp
      cmn%weath14_carb = 0._wp
      cmn%weath14_sil = 0._wp
      do j=1,nj
        do i=1,ni
          if (lnd%l2d(i,j)%f_veg.gt.0._wp) then
            ! index of runoff destination cells (i,j) -> (ir,jr)
            ir = cmn%i_runoff(i,j)
            jr = cmn%j_runoff(i,j)
            ! compute total area of neighbors where runoff is distributed
            area_nbr = 0._wp
            nbr = min(cmn%coast_nbr(ir,jr),n_cells_dist_weath)
            do n=1,nbr
              ii = cmn%i_coast_nbr(ir,jr,n)
              jj = cmn%j_coast_nbr(ir,jr,n)
              area_nbr = area_nbr + area(ii,jj)*cmn%f_ocn(ii,jj)  ! m2
            enddo
            ! transfer weathering to ocean domain
            lnd_cell_area = area(i,j)*lnd%l2d(i,j)%f_veg 
            do n=1,nbr
              ii = cmn%i_coast_nbr(ir,jr,n)
              jj = cmn%j_coast_nbr(ir,jr,n)
              ! organic carbon fluxes (DOC and POC)
              cmn%doc_export(ii,jj) = cmn%doc_export(ii,jj) &   ! kgC/m2/yr
                + lnd%l2d(i,j)%doc_export * lnd_cell_area/area_nbr  ! adjust for area of grid cells
              cmn%doc13_export(ii,jj) = cmn%doc13_export(ii,jj) &   ! kgC/m2/yr
                + lnd%l2d(i,j)%doc13_export * lnd_cell_area/area_nbr  ! adjust for area of grid cells
              cmn%doc14_export(ii,jj) = cmn%doc14_export(ii,jj) &   ! kgC/m2/yr
                + lnd%l2d(i,j)%doc14_export * lnd_cell_area/area_nbr  ! adjust for area of grid cells
              cmn%poc_export(ii,jj) = cmn%poc_export(ii,jj) &   ! kgC/m2/yr
                + lnd%l2d(i,j)%poc_export * lnd_cell_area/area_nbr  ! adjust for area of grid cells
              cmn%poc13_export(ii,jj) = cmn%poc13_export(ii,jj) &   ! kgC/m2/yr
                + lnd%l2d(i,j)%poc13_export * lnd_cell_area/area_nbr  ! adjust for area of grid cells
              cmn%poc14_export(ii,jj) = cmn%poc14_export(ii,jj) &   ! kgC/m2/yr
                + lnd%l2d(i,j)%poc14_export * lnd_cell_area/area_nbr  ! adjust for area of grid cells
              ! carbonate and silicate weathering
              if (l_weathering) then
                cmn%weath_carb(ii,jj) = cmn%weath_carb(ii,jj) &   ! mol C/m2/yr
                  + lnd%l2d(i,j)%weath_carb * lnd_cell_area/area_nbr  ! adjust for area of grid cells
                cmn%weath_sil(ii,jj) = cmn%weath_sil(ii,jj) &     ! mol C/m2/yr
                  + lnd%l2d(i,j)%weath_sil * lnd_cell_area/area_nbr  ! adjust for area of grid cells
                cmn%weath13_carb(ii,jj) = cmn%weath13_carb(ii,jj) &   ! mol C/m2/yr
                  + lnd%l2d(i,j)%weath13_carb * lnd_cell_area/area_nbr  ! adjust for area of grid cells
                cmn%weath13_sil(ii,jj) = cmn%weath13_sil(ii,jj) &     ! mol C/m2/yr
                  + lnd%l2d(i,j)%weath13_sil * lnd_cell_area/area_nbr  ! adjust for area of grid cells
                cmn%weath14_carb(ii,jj) = cmn%weath14_carb(ii,jj) &   ! mol C/m2/yr
                  + lnd%l2d(i,j)%weath14_carb * lnd_cell_area/area_nbr  ! adjust for area of grid cells
                cmn%weath14_sil(ii,jj) = cmn%weath14_sil(ii,jj) &     ! mol C/m2/yr
                  + lnd%l2d(i,j)%weath14_sil * lnd_cell_area/area_nbr  ! adjust for area of grid cells
              endif
            enddo
          endif
        enddo
      enddo

      !------------------------------------------------------------
      ! compute globally integrated alkalinity flux from land due to weathering
      if (l_spinup_cc) then
        if (year.le.(nyears_spinup_bgc-2) .and. year.gt.(nyears_spinup_bgc-10-2)) then
          ! average over 10 years to remove possible noise
          avg = 1._wp/10._wp
          do j=1,nj
            do i=1,ni
              if (cmn%f_ocn(i,j).gt.0._wp) then
                cmn%alk_from_lnd_glob = cmn%alk_from_lnd_glob + (cmn%weath_carb(i,j)+cmn%weath_sil(i,j))*area(i,j)*cmn%f_ocn(i,j) * avg
              endif
            enddo
          enddo
        endif
        if (year.eq.(nyears_spinup_bgc-2)) then
          ! compute weathering scaling needed in order to match alkalinity flux
          ! from weathering with alkalinity flux needed to keep alkalinity stable in the ocean
          lnd%l0d%weath_scale = cmn%alk_to_ocn_glob/cmn%alk_from_lnd_glob
          print *,'weath_scale',lnd%l0d%weath_scale
        endif
      endif

    endif

    if (time_eoy_lnd) then
      ! net land carbon flux to atmosphere
      cmn%delta_C_lnd   = lnd%l0d%Cflx_atm_lnd ! kgC/yr
      cmn%delta_C13_lnd = lnd%l0d%C13flx_atm_lnd ! kgC13/yr
      cmn%delta_C14_lnd = lnd%l0d%C14flx_atm_lnd ! kgC14/yr
      ! annual methane emissions
      cmn%ch4_flx_lnd = lnd%l0d%ch4_emis        ! kgCH4/yr
    endif

    cmn%Cflx_lnd_avg = lnd%l0d%Cflx_avg
    cmn%weath_carb_avg = lnd%l0d%weath_carb_avg
    cmn%weath_sil_avg  = lnd%l0d%weath_sil_avg

   return

  end subroutine lnd_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ b g c
  ! Purpose    :  from common grid (and ocean) to bgc
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_bgc(cmn,ocn,bgc)
    type(cmn_class), intent(inout) :: cmn
    type(ocn_class), intent(in) :: ocn
    type(bgc_class), intent(inout) :: bgc

    integer :: i, j, n, l
    integer :: ocn_k1, bgc_k1
    real(wp) :: avg, albedo
    real(wp) :: runoff_tot
    real(wp) :: runoff_weight


    ! get updated grid from ocean
    if (time_soy_bgc) then
      bgc%grid%ocn_area_tot          = ocn%grid%ocn_area_tot
      bgc%grid%ocn_area              = ocn%grid%ocn_area
      bgc%grid%ocn_area_old          = ocn%grid%ocn_area_old
      bgc%grid%mask_ocn              = ocn%grid%mask_ocn
      bgc%grid%kbo                   = -ocn%grid%k1+maxk+1
      bgc%grid%ocn_vol_tot           = ocn%grid%ocn_vol_tot
      bgc%grid%ocn_vol(:,:,:)        = ocn%grid%ocn_vol(:,:,maxk:1:-1)
      bgc%grid%mask3d(:,:,:)         = ocn%grid%mask_c(:,:,maxk:1:-1)
      bgc%grid%layer_thk(1:maxk)     = ocn%grid%dz(maxk:1:-1)
      bgc%grid%level_depth(1:maxk+1) = -ocn%grid%zw(maxk:0:-1)
      bgc%grid%layer_depth(1:maxk)   = -ocn%grid%zro(maxk:1:-1)
      bgc%grid%coral_f_area = cmn%coral_f_area(:,:,nint(cmn%sea_level)-50:nint(cmn%sea_level)-1) ! top 50 m of ocean considering sea level
      ! todo
      bgc%grid%coral_f_topo = cmn%coral_f_topo(:,:,nint(cmn%sea_level)-50:nint(cmn%sea_level)-1) ! top 50 m of ocean considering sea level
    endif

    ! atmospheric CO2
    bgc%atm(iatmco2) = cmn%co2  ! ppm
    bgc%atm(iatmc13) = cmn%c13_c12_atm*cmn%co2   ! ppm
    bgc%atm(iatmc14) = cmn%c14_c_atm*cmn%co2    ! ppm

    ! total runoff
    if (time_soy_bgc) then
      runoff_tot = sum(cmn%runoff_o_ann * area*cmn%f_ocn)  ! kg/s
    endif

    !$omp parallel do collapse(2) private(i,j,n,l,bgc_k1,ocn_k1,albedo,runoff_weight)
    do j=1,maxj
      do i=1,maxi
        ocn_k1 = ocn%grid%k1(i,j)
        if (ocn_k1.le.ocn%grid%nk) then
          n = bgc%grid%id_map(i,j)

          bgc_k1 = bgc%grid%kbo(i,j)

          ! climate input fields
          bgc%clim%temp(i,j,1:bgc_k1)    = ocn%ts(i,j,maxk:ocn_k1:-1,1)  ! temp
          !bgc%clim%temp(i,j,1:bgc_k1)    = ocn%ts(i,j,maxk:ocn_k1:-1,1)-1._wp  ! temp
          bgc%clim%sal(i,j,1:bgc_k1)     = ocn%ts(i,j,maxk:ocn_k1:-1,2)  ! salt
          if (i_bgc_fw.eq.1) then
            bgc%clim%fw(i,j)        = ocn%p_e_sic(i,j) + cmn%runoff_o(i,j) + cmn%calving_o(i,j) + cmn%bmelt_grd_o(i,j) + cmn%bmelt_flt_o(i,j) ! net surface freshwater flux
          else if (i_bgc_fw.eq.2) then
            bgc%clim%fw(i,j)        = ocn%fw(i,j)-ocn%fw_noise(i,j)-ocn%fw_hosing(i,j) ! net surface freshwater flux
          else if (i_bgc_fw.eq.3) then
            bgc%clim%fw(i,j)        = ocn%fw(i,j)-ocn%fw_noise(i,j) ! net surface freshwater flux
          else if (i_bgc_fw.eq.4) then
            bgc%clim%fw(i,j)        = ocn%fw(i,j) ! net surface freshwater flux
          else if (i_bgc_fw.eq.5) then
            bgc%clim%fw(i,j)        = ocn%fw_corr(i,j) ! net surface freshwater flux
          endif
          bgc%clim%swr_toa_24h(j,:) = cmn%solar(doy,:,j)  ! SW at TOA with diurnal cycle, temporal resolution 1h
          bgc%clim%daylength(j)   = cmn%daylength(doy,j) ! daylength [h]
          albedo = (1._wp-cmn%cld(i,j))*(frac_vu*cmn%alb_vis_dir(i,j,i_surf_macro_ocn)+(1._wp-frac_vu)*cmn%alb_nir_dir(i,j,i_surf_macro_ocn)) &
                 + cmn%cld(i,j)*(frac_vu*cmn%alb_vis_dif(i,j,i_surf_macro_ocn)+(1._wp-frac_vu)*cmn%alb_nir_dif(i,j,i_surf_macro_ocn))
          bgc%clim%swr_sur_dw(i,j) = cmn%swnet(i,j,i_surf_macro_ocn)/(1._wp-albedo)  ! SW down radiation flux at the surface
          bgc%clim%ice_frac(i,j)  = (cmn%f_sic(i,j)*cmn%f_ocn2(i,j) + cmn%f_ice_flt(i,j)) / cmn%f_ocn(i,j) ! ice fraction (sea ice + floating shelf ice)
          !if (cmn%grid%lat(1,j)<-45._wp) then
          !  bgc%clim%ice_frac(i,j) = 1._wp
          !endif
          bgc%clim%pressure(i,j)  = cmn%slp(i,j) ! sea level pressure
          bgc%clim%wind10(i,j)    = cmn%wind(i,j,i_surf_macro_ocn) ! 10 m wind speed
          if (time_soy_bgc) then
            bgc%clim%sst_min(i,j) = ocn%sst_min(i,j)    ! minimum annual surface layer temperature [degC]
            bgc%clim%sst_max(i,j) = ocn%sst_max(i,j)    ! maximum annual surface layer temperature [degC]
          endif

          ! dust deposition
          bgc%bgc_1d(n)%flux%dustdep = cmn%dust_dep(i,j)    ! kg/m2/s?
          !if (cmn%grid%lat(1,j)<-40._wp) then
         !   bgc%bgc_1d(n)%flux%dustdep = 3._wp*bgc%bgc_1d(n)%flux%dustdep
          !endif

          ! input from continental weathering 
          if (time_soy_bgc .and. year.gt.1) then

            if (.not.l_sediments) then
              ! sediment model not active

              bgc%bgc_1d(n)%flux%inp_doc   = 0._wp 
              bgc%bgc_1d(n)%flux%inp_doc13 = 0._wp 
              bgc%bgc_1d(n)%flux%inp_doc14 = 0._wp 
              bgc%bgc_1d(n)%flux%inp_poc   = 0._wp  
              bgc%bgc_1d(n)%flux%inp_poc13 = 0._wp 
              bgc%bgc_1d(n)%flux%inp_poc14 = 0._wp 
              bgc%bgc_1d(n)%flux%inp_dic   = 0._wp 
              bgc%bgc_1d(n)%flux%inp_dic13 = 0._wp 
              bgc%bgc_1d(n)%flux%inp_dic14 = 0._wp 
              bgc%bgc_1d(n)%flux%inp_alk   = 0._wp 
              bgc%bgc_1d(n)%flux%inp_sil   = 0._wp 

            else
              ! sediment model active

              runoff_weight = cmn%runoff_o_ann(i,j)/runoff_tot ! kg/m2/s / kg *s = m-2
              cmn%runoff_o_ann(i,j) = 0._wp

              if (l_spinup_bgc) then
                ! balance net sediment fluxes during bgc spinup to conserve water colum inventories

                ! organic carbon input to balance sedimentation flux
                if (i_compensate.eq.1) then
                  bgc%bgc_1d(n)%flux%inp_doc   = bgc%flx_sed%poc(i,j)   / sec_year ! kmol P m-2 s-1 
                  bgc%bgc_1d(n)%flux%inp_doc13 = bgc%flx_sed%poc13(i,j) / sec_year ! kmol C m-2 s-1 
                  bgc%bgc_1d(n)%flux%inp_doc14 = bgc%flx_sed%poc14(i,j) / sec_year ! kmol C m-2 s-1 
                else if (i_compensate.eq.2) then
                  bgc%bgc_1d(n)%flux%inp_doc   = bgc%flx_sed%poc_tot   * runoff_weight / sec_year ! kmol P m-2 s-1 
                  bgc%bgc_1d(n)%flux%inp_doc13 = bgc%flx_sed%poc13_tot * runoff_weight / sec_year ! kmol C m-2 s-1 
                  bgc%bgc_1d(n)%flux%inp_doc14 = bgc%flx_sed%poc14_tot * runoff_weight / sec_year ! kmol C m-2 s-1 
                endif
                bgc%bgc_1d(n)%flux%inp_poc   = 0._wp 
                bgc%bgc_1d(n)%flux%inp_poc13 = 0._wp
                bgc%bgc_1d(n)%flux%inp_poc14 = 0._wp
                ! CaCO3 input to balance sedimentation flux
                if (i_compensate.eq.1) then
                  bgc%bgc_1d(n)%flux%inp_dic   = bgc%flx_sed%caco3(i,j)   / sec_year ! kmol C m-2 s-1
                  bgc%bgc_1d(n)%flux%inp_dic13 = bgc%flx_sed%caco313(i,j) / sec_year ! kmol C m-2 s-1
                  bgc%bgc_1d(n)%flux%inp_dic14 = bgc%flx_sed%caco314(i,j) / sec_year ! kmol C m-2 s-1
                else if (i_compensate.eq.2) then
                  bgc%bgc_1d(n)%flux%inp_dic   = bgc%flx_sed%caco3_tot   * runoff_weight / sec_year ! kmol C m-2 s-1
                  bgc%bgc_1d(n)%flux%inp_dic13 = bgc%flx_sed%caco313_tot * runoff_weight / sec_year ! kmol C m-2 s-1
                  bgc%bgc_1d(n)%flux%inp_dic14 = bgc%flx_sed%caco314_tot * runoff_weight / sec_year ! kmol C m-2 s-1
                endif
                ! Alkalinity input to balance sedimentation flux
                if (i_compensate.eq.1) then
                  bgc%bgc_1d(n)%flux%inp_alk   = 2._wp*bgc%flx_sed%caco3(i,j) / sec_year ! kmol m-2 s-1
                else if (i_compensate.eq.2) then
                  bgc%bgc_1d(n)%flux%inp_alk   = 2._wp*bgc%flx_sed%caco3_tot * runoff_weight / sec_year ! kmol m-2 s-1
                endif
                ! Silica input to balance sedimentation flux
                if (i_compensate.eq.1) then
                  bgc%bgc_1d(n)%flux%inp_sil   = bgc%flx_sed%opal(i,j) / sec_year ! kmol Si m-2 s-1
                else if (i_compensate.eq.2) then
                  bgc%bgc_1d(n)%flux%inp_sil   = bgc%flx_sed%opal_tot * runoff_weight / sec_year ! kmol Si m-2 s-1
                endif

              else
                ! appy weathering from land or compensate for sediment burial, depending on the settings
                    
                if (l_conserve_phos) then
                  ! organic carbon input to balance sediment burial flux
                  if (i_compensate.eq.1) then
                    bgc%bgc_1d(n)%flux%inp_doc   = bgc%flx_sed%poc(i,j)   / sec_year ! kmol P m-2 s-1 
                    bgc%bgc_1d(n)%flux%inp_doc13 = bgc%flx_sed%poc13(i,j) / sec_year ! kmol C m-2 s-1 
                    bgc%bgc_1d(n)%flux%inp_doc14 = bgc%flx_sed%poc14(i,j) / sec_year ! kmol C m-2 s-1 
                  else if (i_compensate.eq.2) then
                    bgc%bgc_1d(n)%flux%inp_doc   = bgc%flx_sed%poc_tot   * runoff_weight / sec_year ! kmol P m-2 s-1 
                    bgc%bgc_1d(n)%flux%inp_doc13 = bgc%flx_sed%poc13_tot * runoff_weight / sec_year ! kmol C m-2 s-1 
                    bgc%bgc_1d(n)%flux%inp_doc14 = bgc%flx_sed%poc14_tot * runoff_weight / sec_year ! kmol C m-2 s-1 
                  endif
                  bgc%bgc_1d(n)%flux%inp_poc   = 0._wp 
                  bgc%bgc_1d(n)%flux%inp_poc13 = 0._wp
                  bgc%bgc_1d(n)%flux%inp_poc14 = 0._wp
                else
                  ! dissolved organic carbon input from rivers, assuming same Redfield ratio on land and in the ocean (supported by Cleveland 2007)
                  bgc%bgc_1d(n)%flux%inp_doc   = cmn%doc_export(i,j)   /rcar/12._wp / sec_year    ! kmol P m-2 s-1 
                  bgc%bgc_1d(n)%flux%inp_doc13 = cmn%doc13_export(i,j) /12._wp / sec_year    ! kmol C m-2 s-1 
                  bgc%bgc_1d(n)%flux%inp_doc14 = cmn%doc14_export(i,j) /12._wp / sec_year    ! kmol C m-2 s-1 
                  ! particulate organic carbon input from rivers, assuming same Redfield ratio on land and in the ocean (supported by Cleveland 2007)
                  bgc%bgc_1d(n)%flux%inp_poc   = cmn%poc_export(i,j)   /rcar/12._wp / sec_year    ! kmol P m-2 s-1 
                  bgc%bgc_1d(n)%flux%inp_poc13 = cmn%poc13_export(i,j) /12._wp / sec_year    ! kmol C m-2 s-1 
                  bgc%bgc_1d(n)%flux%inp_poc14 = cmn%poc14_export(i,j) /12._wp / sec_year    ! kmol C m-2 s-1 
                endif
                ! dissolved inorganic carbon input from weathering (NO short-circuiting of atmosphere)
                bgc%bgc_1d(n)%flux%inp_dic   = (cmn%weath_sil(i,j)+cmn%weath_carb(i,j))*1.e-3_wp / sec_year  ! kmol C m-2 s-1
                bgc%bgc_1d(n)%flux%inp_dic13 = (cmn%weath13_sil(i,j)+cmn%weath13_carb(i,j))*1.e-3_wp / sec_year  ! kmol C m-2 s-1
                bgc%bgc_1d(n)%flux%inp_dic14 = (cmn%weath14_sil(i,j)+cmn%weath14_carb(i,j))*1.e-3_wp / sec_year  ! kmol C m-2 s-1
                if (l_conserve_alk) then
                  ! Alkalinity input to balance sediment burial flux
                  if (i_compensate.eq.1) then
                    bgc%bgc_1d(n)%flux%inp_alk   = 2._wp*bgc%flx_sed%caco3(i,j) / sec_year ! kmol m-2 s-1
                  else if (i_compensate.eq.2) then
                    bgc%bgc_1d(n)%flux%inp_alk   = 2._wp*bgc%flx_sed%caco3_tot * runoff_weight / sec_year ! kmol m-2 s-1
                  endif
                else
                  ! alkalinity input from weathering
                  bgc%bgc_1d(n)%flux%inp_alk = (cmn%weath_sil(i,j)+cmn%weath_carb(i,j))*1.e-3_wp / sec_year  ! kmol m-2 s-1
                endif
                if (l_conserve_sil) then
                  ! Silica input to balance sediment burial flux
                  if (i_compensate.eq.1) then
                    bgc%bgc_1d(n)%flux%inp_sil   = bgc%flx_sed%opal(i,j) / sec_year ! kmol Si m-2 s-1
                  else if (i_compensate.eq.2) then
                    bgc%bgc_1d(n)%flux%inp_sil   = bgc%flx_sed%opal_tot * runoff_weight / sec_year ! kmol Si m-2 s-1
                endif
                else
                  ! Silica input from silicate weathering 
                  ! fixme todo: n_bic/n_sil=1.75 approx. after Munhoven 1997 PhD thesis, page 78
                  !bgc%bgc_1d(n)%flux%inp_sil = 1._wp/1.75_wp*cmn%weath_sil(i,j)*1.e-3_wp / sec_year   ! kmol Si m-2 s-1
                  bgc%bgc_1d(n)%flux%inp_sil = cmn%weath_sil(i,j)*1.e-3_wp / sec_year   ! kmol Si m-2 s-1
                endif
              endif
            endif
          endif

          ! transfer bgc tracers
          do l=1,BGC_NTRA
            !bgc%bgc_1d(n)%tra%ocetra(:, l) = ocn%ts( i, j, maxk:1:-1, n_tracers_ocn+l)
            bgc%bgc_1d(n)%tra%ocetra(1:bgc_k1, l) = ocn%ts(i, j, maxk:ocn_k1:-1,n_tracers_ocn+l )
          enddo

          !------------------------------------------------------------
          ! save daily input fields for offline simulation, if required
          if (time_call_daily_input_save) then
            ! average over nyear_avg_offline years to remove possible noise
            avg = real(nday_year,wp)/real(nstep_year_bgc,wp) / real(nyear_avg_offline,wp)
            bgc%daily_input_save%temp(i,j,1:bgc_k1,doy) = bgc%daily_input_save%temp(i,j,1:bgc_k1,doy) + bgc%clim%temp(i,j,1:bgc_k1)  * avg     
            bgc%daily_input_save%sal(i,j,1:bgc_k1,doy)  = bgc%daily_input_save%sal(i,j,1:bgc_k1,doy)  + bgc%clim%sal(i,j,1:bgc_k1)   * avg
            bgc%daily_input_save%fw(i,j,doy)            = bgc%daily_input_save%fw(i,j,doy)            + bgc%clim%fw(i,j)             * avg 
            bgc%daily_input_save%swr_sur_dw(i,j,doy)    = bgc%daily_input_save%swr_sur_dw(i,j,doy)    + bgc%clim%swr_sur_dw(i,j)     * avg
            bgc%daily_input_save%ice_frac(i,j,doy)      = bgc%daily_input_save%ice_frac(i,j,doy)      + bgc%clim%ice_frac(i,j)       * avg  
            bgc%daily_input_save%pressure(i,j,doy)      = bgc%daily_input_save%pressure(i,j,doy)      + bgc%clim%pressure(i,j)       * avg
            bgc%daily_input_save%wind10(i,j,doy)        = bgc%daily_input_save%wind10(i,j,doy)        + bgc%clim%wind10(i,j)         * avg
            bgc%daily_input_save%dust_dep(i,j,doy)      = bgc%daily_input_save%dust_dep(i,j,doy)      + bgc%bgc_1d(n)%flux%dustdep   * avg 
          endif
          !------------------------------------------------------------
          ! use saved daily input fields now
          if (time_use_daily_input_save) then
            bgc%clim%temp(i,j,1:bgc_k1) = bgc%daily_input_save%temp(i,j,1:bgc_k1,doy)        
            bgc%clim%sal(i,j,1:bgc_k1)  = bgc%daily_input_save%sal(i,j,1:bgc_k1,doy)    
            bgc%clim%fw(i,j)            = bgc%daily_input_save%fw(i,j,doy)               
            bgc%clim%swr_sur_dw(i,j)    = bgc%daily_input_save%swr_sur_dw(i,j,doy)       
            bgc%clim%ice_frac(i,j)      = bgc%daily_input_save%ice_frac(i,j,doy)          
            bgc%clim%pressure(i,j)      = bgc%daily_input_save%pressure(i,j,doy)        
            bgc%clim%wind10(i,j)        = bgc%daily_input_save%wind10(i,j,doy)          
            bgc%bgc_1d(n)%flux%dustdep  = bgc%daily_input_save%dust_dep(i,j,doy)         
          endif
          !------------------------------------------------------------

        endif
      enddo
    enddo
    !$omp end parallel do


  end subroutine cmn_to_bgc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  b g c _ t o _ c m n
  ! Purpose    :  from bgc to common grid (and ocean)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bgc_to_cmn(bgc,cmn,ocn)

    type(bgc_class), intent(in) :: bgc
    type(cmn_class), intent(inout) :: cmn
    type(ocn_class), intent(inout) :: ocn

    integer :: i, j, n, l
    real(wp) :: cal, sil, tot
    real(wp) :: avg

    
    !$omp parallel do collapse(2) private(i,j,n,l,cal,sil,tot)
    do j = 1,maxj
      do i = 1,maxi
        n = bgc%grid%id_map(i,j)
        do l = 1,BGC_NTRA
          ocn%ts(i, j, maxk:1:-1, n_tracers_ocn+l) = bgc%bgc_1d(n)%tra%ocetra(1:maxk, l)
        enddo
        cmn%delta_C_ocn_2d(i,j)   = bgc%bgc_1d(n)%delta_C ! kgC/s
        cmn%delta_C13_ocn_2d(i,j) = bgc%bgc_1d(n)%delta_C13 ! kgC13/s
        cmn%delta_C14_ocn_2d(i,j) = bgc%bgc_1d(n)%delta_C14 ! kgC14/s
        if (time_eoy_bgc) then
          ! carbonate fraction in buried sediments needed for weathering on shelf
          cal = bgc%bgc_1d(n)%sed%burial(isssc12) * 3.85e-2_wp  ! calcium carbonate
          sil = bgc%bgc_1d(n)%sed%burial(issssil) * 2.73e-2_wp  ! silicate
          !org = bgc%bgc_1d(n)%sed%burial(issso12) * 3.e-2_wp    ! organic carbon
          !cla = bgc%bgc_1d(n)%sed%burial(issster) * 3.85e-4_wp  ! clay
          tot = cal + sil ! total weight
          if (tot.gt.0._wp) then
            cmn%f_carb(i,j) = cal/tot
          else
            cmn%f_carb(i,j) = 0._wp
          endif
        endif
      enddo
    enddo
    !$omp end parallel do

    if (time_eoy_bgc) then
      cmn%delta_C_ocn   = bgc%delta_C ! kgC/yr
      cmn%delta_C13_ocn = bgc%delta_C13 ! kgC13/yr
      cmn%delta_C14_ocn = bgc%delta_C14 ! kgC14/yr
    endif

    cmn%Cflx_ocn_avg = bgc%Cflx_avg ! kgC/yr

    !------------------------------------------------------------
    ! compute globally integrated alkalinity input to ocean 
    if (l_spinup_cc .and. time_soy_bgc .and. year.le.(nyears_spinup_bgc-2) .and. year.gt.(nyears_spinup_bgc-10-2)) then
      ! average over 10 years to remove possible noise
      avg = 1._wp/10._wp
      do j=1,maxj
        do i=1,maxi
          if (ocn%grid%k1(i,j).le.ocn%grid%nk) then
            n = bgc%grid%id_map(i,j)
            cmn%alk_to_ocn_glob = cmn%alk_to_ocn_glob + bgc%bgc_1d(n)%flux%inp_alk*1e3_wp*area(i,j)*cmn%f_ocn(i,j)*sec_year * avg ! kmol/m2/s * mol/kmol * m2 * s/yr = mol/yr
          endif
        enddo
      enddo
    endif

  end subroutine bgc_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ c o 2
  ! Purpose    :  from common grid to co2 model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_co2(cmn,co2)

    type(cmn_class), intent(in) :: cmn
    type(co2_class), intent(inout) :: co2

    integer :: i, j
    logical, save :: firstcall = .true.


    ! ocean carbon fluxes
    if (flag_bgc) then
      co2%dCocn_dt   = cmn%delta_C_ocn   ! kgC/yr
      co2%dC13ocn_dt = cmn%delta_C13_ocn ! kgC13/yr
      co2%dC14ocn_dt = cmn%delta_C14_ocn ! kgC14/yr
    else
      co2%dCocn_dt   = 0._wp 
      co2%dC13ocn_dt = 0._wp 
      co2%dC14ocn_dt = 0._wp 
      co2%dCocn_dt_avg = 0._wp
    endif

    ! land carbon fluxes
    if (flag_lnd) then
      co2%dClnd_dt   = cmn%delta_C_lnd   ! kgC/yr
      co2%dC13lnd_dt = cmn%delta_C13_lnd ! kgC13/yr
      co2%dC14lnd_dt = cmn%delta_C14_lnd ! kgC14/yr
    else
      co2%dClnd_dt   = 0._wp 
      co2%dC13lnd_dt = 0._wp 
      co2%dC14lnd_dt = 0._wp 
    endif

    ! weathering fluxes
    if (flag_lnd) then
      co2%dCweath_dt = 0._wp
      co2%dC13weath_dt = 0._wp
      co2%dC14weath_dt = 0._wp
      do j=1,nj
        do i=1,ni
          if (cmn%f_ocn(i,j).gt.0._wp) then     ! cmn%weath_* fluxes are on the ocean grid!
            co2%dCweath_dt = co2%dCweath_dt &   ! kgC/yr
              + (0.5_wp*cmn%weath_carb(i,j) * area(i,j)*cmn%f_ocn(i,j) &    ! mol C/yr, 0.5 mole of CO2 per mole of HCO3- removed by carbonate weathering
              + cmn%weath_sil(i,j) * area(i,j)*cmn%f_ocn(i,j)) &            ! mol C/yr, one mole of CO2 per mole of HCO3- removed by silicate weathering
              * 12._wp*1.e-3_wp     ! mol C -> kgC
            co2%dC13weath_dt = co2%dC13weath_dt &   ! kgC13/yr
              + (0.5_wp*cmn%weath_carb(i,j)*cmn%c13_c12_atm * area(i,j)*cmn%f_ocn(i,j) &  ! mol C13/yr, only C13 weathering from atmosphere 
              + cmn%weath13_sil(i,j) * area(i,j)*cmn%f_ocn(i,j)) &          ! mol C13/yr, all C13 from silicate weathering comes from the atmosphere 
              * 12._wp*1.e-3_wp     ! mol C -> kgC
            co2%dC14weath_dt = co2%dC14weath_dt &   ! kgC14/yr
              + (cmn%weath14_carb(i,j) * area(i,j)*cmn%f_ocn(i,j) &         ! mol C14/yr, all weathered C14 comes from the atmosphere 
              + cmn%weath14_sil(i,j) * area(i,j)*cmn%f_ocn(i,j)) &          ! mol C14/yr, all weathered C14 comes from the atmosphere 
              * 12._wp*1.e-3_wp     ! mol C -> kgC
          endif
        enddo
      enddo
    else
      co2%dCweath_dt   = 0._wp 
      co2%dC13weath_dt = 0._wp 
      co2%dC14weath_dt = 0._wp 
    endif

    if (firstcall) then
      co2%weath_sil_avg = cmn%weath_sil_avg     ! kgC/yr
      co2%dCvolc_dt_eq = 0.5_wp*cmn%weath_carb_avg + cmn%weath_sil_avg + cmn%Cflx_ocn_avg + cmn%Cflx_lnd_avg       ! kgC/yr
    endif

    ! anthropogenic CH4 emissions to CO2 equivalent, exclude emissions from agriculture to avoid double-counting!
    co2%dCH4_dt = (1._wp-cmn%f_ch4emis_agro)*cmn%ch4_emis * 12._wp/16._wp  ! kgCH4 <- kgC

    ! global annual temperature change for extra emissions 
    if (co2%T_glob.eq.0._wp) co2%T_glob = cmn%t2m_glob_ann  ! initialize 
    co2%dT_glob_cum = co2%dT_glob_cum + (cmn%t2m_glob_ann - co2%T_glob)
    co2%T_glob = cmn%t2m_glob_ann

    firstcall = .false.

    return
    
  end subroutine cmn_to_co2


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c o 2 _ t o _ c m n
  ! Purpose    :  from co2 to common
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine co2_to_cmn(co2,cmn)

    type(co2_class), intent(in) :: co2
    type(cmn_class), intent(inout) :: cmn

    cmn%co2 = co2%co2   ! ppm
    if (ico2_rad.eq.0) cmn%co2_rad = cmn%co2   ! ppm
    cmn%c13_c12_atm = co2%c13_c12
    cmn%c14_c_atm   = co2%c14_c

  end subroutine co2_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ c h 4
  ! Purpose    :  from common grid to ch4 model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_ch4(cmn,ch4)

    type(cmn_class), intent(in) :: cmn
    type(ch4_class), intent(inout) :: ch4


    ! ocean methane fluxes
    ch4%dch4ocn_dt   = 0._wp 

    ! land methane fluxes
    if (flag_lnd) then
      ch4%dch4lnd_dt   = cmn%ch4_flx_lnd   ! kgCH4/yr
    else
      ch4%dch4lnd_dt   = 0._wp 
    endif

    return
    
  end subroutine cmn_to_ch4


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c h 4 _ t o _ c m n
  ! Purpose    :  from ch4 to common
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ch4_to_cmn(ch4,cmn)

    type(ch4_class), intent(in) :: ch4
    type(cmn_class), intent(inout) :: cmn

    cmn%ch4_emis = ch4%dch4emis_dt ! kgCH4/yr
    cmn%f_ch4emis_agro = ch4%f_ch4emis_agro

    cmn%ch4 = ch4%ch4   ! ppb
    if (ich4_rad.eq.0) cmn%ch4_rad = cmn%ch4   ! ppb

  end subroutine ch4_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ s m b
  ! Purpose    :  transfer common grid variables to smb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_smb(cmn,smb_in)

    implicit none

    type(cmn_class) :: cmn
    type(smb_in_class) :: smb_in

    integer :: i, j, n, ns, nc
    real(wp) :: w1, w2
    real(wp) :: z(nsurf_macro)
    real(wp) :: lw(nsurf_macro)
    real(wp), parameter :: gam_lw     = -2.5e-2_wp
    integer, parameter :: nfil = -1


    ! CO2
    smb_in%co2 = cmn%co2
    ! maximum summer insolation at 65N
    smb_in%Smax65N = maxval(cmn%solarm(:,minloc(abs(cmn%grid%lat(1,:)-65._wp))))

    !$omp parallel do private(i,j,n,ns,nc,z,lw,w1,w2)
    do j=1,nj
      do i=1,ni

        smb_in%z_sur(i,j) = cmn%z_sur(i,j)  ! grid cell mean elevation relative to sea level [m]
        smb_in%t2m(i,j)   = sum(cmn%t2m(i,j,:)*cmn%f_stp(i,j,:))   ! 2m temperature at mean surface elevation on coupler grid [K]
        smb_in%tam(i,j)   = cmn%tam(i,j)   ! atmospheric temperature on coupler grid [K] 
        smb_in%ram(i,j)   = cmn%ram(i,j)   ! atmospheric relative humidity on coupler grid [/] 
        smb_in%gam(i,j)   = cmn%gams(i,j)   ! atmospheric lapse rate on coupler grid [K/m] 
        smb_in%tstd(i,j)  = 2._wp ! cmn%tstd  ! standard deviation of daily surface air temperature on coupler grid [K] TODO
        smb_in%prc(i,j)   = cmn%prc(i,j)   ! total precipitation rate on coupler grid [kg/m2/s] 
        smb_in%u700(i,j)  = cmn%u700(i,j)  ! zonal wind component at 700 hPa on coupler grid [m/s]
        smb_in%v700(i,j)  = cmn%v700(i,j)  ! meridional wind component at 700 hPa on coupler grid [m/s]
        smb_in%wind(i,j)  = sum(cmn%wind(i,j,:)*cmn%f_stp(i,j,:)) ! surface wind speed on coupler grid [m/s] grid cell mean or over ice? fixme
        smb_in%cld(i,j)   = cmn%cld(i,j)   ! cloud cover fraction on coupler grid [/]
        smb_in%dust(i,j)  = cmn%dust_dep(i,j) ! dust deposition rate on coupler grid [kg/m2/s]
        smb_in%t_ground(i,j) = cmn%f_lnd(i,j)*(cmn%t_soil(i,j)-T0) &    ! annual mean deepest soil layer temperature on coupler grid [degC]
          + (1._wp-cmn%f_lnd(i,j))*(cmn%t_shelf(i,j)-T0) 

        ! downward longwave radiation at the surface, depends only on surface elevation
        ! pass value for grid cell mean elevation + 'lapse rate' to smb
        ns = 0
        do n=1,nsurf_macro
          ! fluxes over ocean are always available, exclude sea ice (same elevation anyway)
          if (n.eq.i_surf_macro_ocn .or. (n.ne.i_surf_macro_sic .and. cmn%f_astp(i,j,n).gt.0._wp)) then   
            ns = ns+1
            z(ns)      = cmn%z_sur_n(i,j,n)
            lw(ns)     = cmn%cld(i,j)*cmn%lwd_cld(i,j,n) + (1._wp-cmn%cld(i,j))*cmn%lwd_cs(i,j,n)
          endif
        enddo
        if (ns.eq.1) then
          ! longwave radiation available only for one surface type
          smb_in%lwdown(i,j)     = lw(ns)
          ! guess gradient
          smb_in%gam_lw(i,j)     = gam_lw
        else
          nc = 0
          do n=1,ns-1
            if (z(n).lt.cmn%z_sur(i,j) .and. z(n+1).ge.(cmn%z_sur(i,j)-epsilon(1.))) then
              nc = n
              exit
            endif
          enddo
          if (nc.eq.0) then
            ! all surface types have the same elevation
            smb_in%lwdown(i,j)     = lw(1)
            ! guess gradient
            smb_in%gam_lw(i,j)     = gam_lw
          else
            ! compute surface lw at grid-cell mean surface elevation
            w1 = (z(nc+1)-cmn%z_sur(i,j))/(z(nc+1)-z(nc))
            w2 = 1._wp-w1
            smb_in%lwdown(i,j)     = w1*lw(nc)     + w2*lw(nc+1)
            ! compute surface lw 'lapse rate'
            smb_in%gam_lw(i,j)     = (lw(nc+1)-lw(nc))/(z(nc+1)-z(nc))
          endif
        endif
        if (smb_in%gam_lw(i,j).gt.-5e-3) then
          print *
          print *,i,j
          print *,nc
          print *,z(1).lt.cmn%z_sur(i,j)
          print *,z(2).ge.cmn%z_sur(i,j)
          print *,'f_astp',cmn%f_astp(i,j,:)
          print *,'f_stp',cmn%f_stp(i,j,:)
          print *,'z_sur',cmn%z_sur(i,j)
          print *,'z',z(1:ns)
          print *,'lw',lw(1:ns)
          print *,smb_in%lwdown(i,j)
          print *,smb_in%gam_lw(i,j)
        endif

        ! grid-cell mean surface albedoes
        smb_in%alb_vis_dir(i,j) = sum(cmn%alb_vis_dir(i,j,:)*cmn%f_stp(i,j,:))
        smb_in%alb_nir_dir(i,j) = sum(cmn%alb_nir_dir(i,j,:)*cmn%f_stp(i,j,:))
        smb_in%alb_vis_dif(i,j) = sum(cmn%alb_vis_dif(i,j,:)*cmn%f_stp(i,j,:))
        smb_in%alb_nir_dif(i,j) = sum(cmn%alb_nir_dif(i,j,:)*cmn%f_stp(i,j,:))

        smb_in%swd_sur_vis_dir(i,j) = cmn%swd_vis_dir(i,j)
        smb_in%swd_sur_nir_dir(i,j) = cmn%swd_nir_dir(i,j)
        smb_in%swd_sur_vis_dif(i,j) = cmn%swd_vis_dif(i,j)
        smb_in%swd_sur_nir_dif(i,j) = cmn%swd_nir_dif(i,j)

        smb_in%dswd_dalb_vis_dir(i,j) = cmn%dswd_dalb_vis_dir(i,j)  
        smb_in%dswd_dalb_nir_dir(i,j) = cmn%dswd_dalb_nir_dir(i,j)
        smb_in%dswd_dalb_vis_dif(i,j) = cmn%dswd_dalb_vis_dif(i,j)
        smb_in%dswd_dalb_nir_dif(i,j) = cmn%dswd_dalb_nir_dif(i,j)
        smb_in%dswd_dz_nir_dir(i,j)   = cmn%dswd_dz_nir_dir(i,j)  
        smb_in%dswd_dz_nir_dif(i,j)   = cmn%dswd_dz_nir_dif(i,j)  

        smb_in%coszm(i,j) = cmn%coszm(doy,j)
        smb_in%swd_toa(i,j) = cmn%solarm(doy,j)
        smb_in%swd_toa_min(i,j) = cmn%solarmin(doy,j)

      enddo

    enddo
    !$omp end parallel do

    ! filtering 
    !$omp parallel sections
    !$omp section
    call filter_smb(smb_in%z_sur) 
    !$omp section
    call filter_smb(smb_in%t2m)   
    !$omp section
    call filter_smb(smb_in%tstd)  
    !$omp section
    call filter_smb(smb_in%tam)    
    !$omp section
    call filter_smb(smb_in%ram)    
    !$omp section
    call filter_smb(smb_in%gam)    
    !$omp section
    call smooth2(smb_in%prc,1)
    call filter_smb(smb_in%prc)   
    !$omp section
    call filter_smb(smb_in%wind)  
    !$omp section
    call filter_smb(smb_in%cld)   
    !$omp section
    call smooth2(smb_in%dust,1)
    call filter_smb(smb_in%dust)  
    !$omp section
    call filter_smb(smb_in%t_ground) 
    !$omp section
    call filter_smb(smb_in%lwdown)  
    !$omp section
    call filter_smb(smb_in%gam_lw)  
    !$omp section
    call filter_smb(smb_in%alb_vis_dir) 
    !$omp section
    call filter_smb(smb_in%alb_nir_dir) 
    !$omp section
    call filter_smb(smb_in%alb_vis_dif) 
    !$omp section
    call filter_smb(smb_in%alb_nir_dif) 
    !$omp section
    call filter_smb(smb_in%swd_sur_vis_dir) 
    !$omp section
    call filter_smb(smb_in%swd_sur_nir_dir) 
    !$omp section
    call filter_smb(smb_in%swd_sur_vis_dif) 
    !$omp section
    call filter_smb(smb_in%swd_sur_nir_dif) 
    !$omp section
    call filter_smb(smb_in%dswd_dalb_vis_dir) 
    !$omp section
    call filter_smb(smb_in%dswd_dalb_nir_dir) 
    !$omp section
    call filter_smb(smb_in%dswd_dalb_vis_dif) 
    !$omp section
    call filter_smb(smb_in%dswd_dalb_nir_dif) 
    !$omp section
    call filter_smb(smb_in%dswd_dz_nir_dir)   
    !$omp section
    call filter_smb(smb_in%dswd_dz_nir_dif)   
    !$omp section
    call filter_smb(smb_in%u700)
    !$omp section
    call filter_smb(smb_in%v700)
    !$omp end parallel sections


    return

  end subroutine cmn_to_smb


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  f i l t e r _ s m b
  ! Purpose    :  filter smb_in variables
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine filter_smb(var2d)

    implicit none

    real(wp), dimension(:,:), intent(inout) :: var2d

    integer :: j, nfil

    do j=1,nj
      if (j.eq.1 .or. j.eq.nj) then
        nfil = -1
      else
        !nfil = nint(1._wp/cos(lat(j)*pi/180._wp))-1
        nfil = nint(0.5_wp/cos(lat(j)*pi/180._wp))-1
      endif
      call filter1d(var2d(:,j), nfil) 
    enddo

    return

  end subroutine filter_smb


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s m b _ t o _ c m n
  ! Purpose    :  transfer smb variables to common grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_to_cmn(smb,cmn)

    implicit none

    type(smb_class) :: smb
    type(cmn_class) :: cmn

    integer :: i, j, ii, jj, ncells_ice
    logical, save :: firstcall = .true.

    real(wp), allocatable, dimension(:,:,:) :: runoff_ice_i_mon
    real(wp), allocatable, dimension(:,:,:) :: melt_ice_i_mon


    if (firstcall .or. time_eoy_smb) then

      ! aggregate monthly runoff from ice sheets computed in smb to coupler grid

      ! initialize runoff to 0 only where smb domain covers coupler grid (to avoid overwriting possible other ice sheet domains)
      allocate(runoff_ice_i_mon(cmn%grid%G%nx,cmn%grid%G%ny,nmon_year))
      allocate(melt_ice_i_mon(cmn%grid%G%nx,cmn%grid%G%ny,nmon_year))
      do i=1,cmn%grid%G%nx
        do j=1,cmn%grid%G%ny
          if (smb%grid_smb_to_cmn%ncells(i,j)>0) then
            runoff_ice_i_mon(i,j,:) = 0._wp
            melt_ice_i_mon(i,j,:) = 0._wp
          endif
        enddo
      enddo

      ! integrate runoff and icemelt over grid-cell smb ice sheet area
      !!$omp parallel do private(ii,jj,i,j)
      do ii=1,smb%grid%G%nx
        do jj=1,smb%grid%G%ny
          if (smb%mask_ice(ii,jj).eq.1) then
          !if (smb%mask_ice(ii,jj).eq.1 .and. smb%h_ice(ii,jj).gt.h_ice_min) then
            i = smb%grid_smb_to_cmn%i_lowres(ii,jj)
            j = smb%grid_smb_to_cmn%j_lowres(ii,jj)
            runoff_ice_i_mon(i,j,:) = runoff_ice_i_mon(i,j,:) + smb%mon_runoff(ii,jj,:)*smb%grid%area(ii,jj)*1.e6_wp  ! kg/m2/s * m2 -> kg/s 
            melt_ice_i_mon(i,j,:)   = melt_ice_i_mon(i,j,:)   + smb%mon_icemelt(ii,jj,:)*smb%grid%area(ii,jj)*1.e6_wp  ! kg/m2/s * m2 -> kg/s 
          endif
        enddo
      enddo
      !!$omp end parallel do

      ! smooth runoff in time 
      do i=1,cmn%grid%G%nx
        do j=1,cmn%grid%G%ny
          if (smb%grid_smb_to_cmn%ncells(i,j)>0) then
            cmn%runoff_ice_i_mon(i,j,:) = relax_run*cmn%runoff_ice_i_mon(i,j,:) + (1._wp-relax_run)*runoff_ice_i_mon(i,j,:)
            cmn%melt_ice_i_mon(i,j,:)   = relax_run*cmn%melt_ice_i_mon(i,j,:)   + (1._wp-relax_run)*melt_ice_i_mon(i,j,:)
          endif
        enddo
      enddo

      deallocate(runoff_ice_i_mon)
      deallocate(melt_ice_i_mon)

    endif

    ! compute average albedo over ice sheet over climate model coarse grid
    ! initialize to 0 only where smb domain covers coupler grid (to avoid overwriting possible other ice sheet domains)
    do i=1,cmn%grid%G%nx
      do j=1,cmn%grid%G%ny
        if (smb%grid_smb_to_cmn%ncells(i,j)>0) then
          cmn%alb_vis_dir_ice_semi(i,j,doy) = 0._wp
          cmn%alb_vis_dif_ice_semi(i,j,doy) = 0._wp
          cmn%alb_nir_dir_ice_semi(i,j,doy) = 0._wp
          cmn%alb_nir_dif_ice_semi(i,j,doy) = 0._wp
        endif
      enddo
    enddo
    !!$omp parallel do private(ii,jj,i,j,ncells_ice)
    do ii=1,smb%grid%G%nx
      do jj=1,smb%grid%G%ny
        if (smb%mask_ice(ii,jj).eq.1) then
          i = smb%grid_smb_to_cmn%i_lowres(ii,jj)
          j = smb%grid_smb_to_cmn%j_lowres(ii,jj)
          ncells_ice = smb%grid_smb_to_cmn%ncells_ice(i,j)
          if (ncells_ice.gt.0) then
            cmn%alb_vis_dir_ice_semi(i,j,doy) = cmn%alb_vis_dir_ice_semi(i,j,doy) + smb%alb_vis_dir(ii,jj)/ncells_ice
            cmn%alb_vis_dif_ice_semi(i,j,doy) = cmn%alb_vis_dif_ice_semi(i,j,doy) + smb%alb_vis_dif(ii,jj)/ncells_ice
            cmn%alb_nir_dir_ice_semi(i,j,doy) = cmn%alb_nir_dir_ice_semi(i,j,doy) + smb%alb_nir_dir(ii,jj)/ncells_ice
            cmn%alb_nir_dif_ice_semi(i,j,doy) = cmn%alb_nir_dif_ice_semi(i,j,doy) + smb%alb_nir_dif(ii,jj)/ncells_ice
          endif
        endif
      enddo
    enddo
    !!$omp end parallel do

    firstcall = .false.

  return

  end subroutine smb_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  i c e _ t o _ s m b
  ! Purpose    :  transfer ice to smb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_to_smb(ice,smb)

    implicit none

    type(ice_class) :: ice
    type(smb_class) :: smb


    if (time_soy_smb) then
      ! surface elevation
      smb%z_sur = ice%z_sur
      where (smb%z_sur.lt.0._wp) 
        smb%z_sur = 0._Wp
      endwhere
      ! ice thickness
      smb%h_ice = ice%h_ice
      ! ice mask
      ! save old ice mask
      smb%mask_ice_old = smb%mask_ice
      ! derive new mask
      where (ice%h_ice.gt.0._wp)
        smb%mask_ice = 1
      elsewhere
        smb%mask_ice = 0
      endwhere
    endif

    return

  end subroutine ice_to_smb


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g e o _ t o _ s m b
  ! Purpose    :  transfer geo to smb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_to_smb(geo,smb)

    implicit none

    type(geo_class) :: geo
    type(smb_class) :: smb

    type(map_scrip_class) :: maps_geo_to_smb

    real(wp), allocatable :: h_ice(:,:)


    ! no ice model, use ice sheets from bnd
    allocate(h_ice(smb%grid%G%nx,smb%grid%G%ny))
    call map_scrip_init(maps_geo_to_smb,geo%hires%grid,smb%grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)
    ! map surface elevation
    call map_scrip_field(maps_geo_to_smb,"topo",geo%hires%z_topo,smb%z_sur,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5._dp*smb%grid%G%dx,smb%grid%G%dx])
    ! map ice thickness
    call map_scrip_field(maps_geo_to_smb,"h_ice",geo%hires%h_ice,h_ice,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5._dp*smb%grid%G%dx,smb%grid%G%dx])
    where (h_ice<10._wp) h_ice = 0._wp

    where (smb%z_sur.lt.0._wp) 
      smb%z_sur = 0._Wp
    endwhere
    ! generate ice mask
    where (h_ice.gt.0._wp) 
      smb%mask_ice = 1
    elsewhere
      smb%mask_ice = 0
    endwhere
    smb%mask_ice_old = smb%mask_ice
    ! ice thickness
    smb%h_ice = h_ice
    deallocate(h_ice)


    return

  end subroutine geo_to_smb


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s m b _ t o _ i c e
  ! Purpose    :  transfer smb to ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_to_ice(smb,ice)

    implicit none

    type(smb_class) :: smb
    type(ice_class) :: ice

    integer :: i, j, ii, jj, i_f, j_f, nx, ny, n_filter
    real(wp) :: sigma_filter, dx, dy, dist, weigh, sum_weigh
    real(wp), dimension(:,:), allocatable :: z_bed_std

    ! annual surface mass balance
    ice%smb = smb%ann_smb/sec_year / rho_i ! kg/m2 -> m(ice equivalent)/s
    ! annual accumulation (total precipitation - evaporation)
    ice%accum  = max(0._wp,smb%ann_prc-smb%ann_evp)/sec_year / rho_i ! kg/m2 -> m(ice equivalent)/s
    ! annual runoff
    ice%runoff = smb%ann_runoff/sec_year / rho_i ! kg/m2 -> m(ice equivalent)/s
    ! ice surface temperature
    ice%temp_s = smb%t_ice ! degC
    ! ground temperature
    ice%temp_g = smb%t_ground ! degC
    ! standard deviation of surface topography
    ice%z_sur_std = smb%z_sur_std ! m
    ! standard deviation of bedrock topography
    ice%z_bed_std = smb%z_bed_std ! m
    ! filter z_bed_std
    nx = smb%grid%G%nx
    ny = smb%grid%G%ny
    allocate(z_bed_std(1:nx,1:ny))
    z_bed_std = ice%z_bed_std
    dx = smb%grid%G%dx ! km
    dy = smb%grid%G%ny ! km
    sigma_filter = 100._wp/dx   ! half span of filtered area, in grid points
    n_filter     = ceiling(2.0_wp*sigma_filter)
    do j=1,ny 
      do i=1,nx
        sum_weigh = 0.0_wp
        ice%z_bed_std(i,j) = 0._wp
        do ii=-n_filter, n_filter
          do jj=-n_filter, n_filter
            i_f = i+ii
            j_f = j+jj
            if (i_f <  1) i_f = 1
            if (i_f > nx) i_f = nx
            if (j_f <  1) j_f = 1
            if (j_f > ny) j_f = ny
            dist      = sqrt(real(ii,wp)**2+real(jj,wp)**2)
            weigh     = exp(-(dist/sigma_filter)**2)
            sum_weigh = sum_weigh + weigh
            ice%z_bed_std(i,j) = ice%z_bed_std(i,j) + weigh*z_bed_std(i_f,j_f)
          end do
        end do
        ice%z_bed_std(i,j) = ice%z_bed_std(i,j)/sum_weigh
      end do
    end do
    deallocate(z_bed_std)


    return

  end subroutine smb_to_ice


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ b m b
  ! Purpose    :  transfer common grid variables to bmb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_bmb(cmn,bmb)

    implicit none

    type(cmn_class) :: cmn
    type(bmb_class) :: bmb

    integer :: nk_ocn, nk_lake


    nk_ocn = size(cmn%z_ocn_bmb)
    if (.not.allocated(bmb%z_ocn_in)) allocate(bmb%z_ocn_in(nk_ocn))
    if (.not.allocated(bmb%t_ocn_in)) allocate(bmb%t_ocn_in(ni,nj,nk_ocn))
    if (.not.allocated(bmb%s_ocn_in)) allocate(bmb%s_ocn_in(ni,nj,nk_ocn))
    if (.not.allocated(bmb%mask_ocn_in)) allocate(bmb%mask_ocn_in(ni,nj,nk_ocn))

    nk_lake = size(cmn%z_lake_bmb)
    if (.not.allocated(bmb%z_lake_in)) allocate(bmb%z_lake_in(nk_lake))
    if (.not.allocated(bmb%t_lake_in)) allocate(bmb%t_lake_in(ni,nj,nk_lake))
    if (.not.allocated(bmb%s_lake_in)) allocate(bmb%s_lake_in(ni,nj,nk_lake))
    if (.not.allocated(bmb%mask_lake_in)) allocate(bmb%mask_lake_in(ni,nj))

    bmb%z_ocn_in = cmn%z_ocn_bmb
    bmb%t_ocn_in(:,:,:) = cmn%t_ocn_bmb(:,:,:)
    bmb%s_ocn_in(:,:,:) = cmn%s_ocn_bmb(:,:,:)
    bmb%mask_ocn_in(:,:,:) = cmn%mask_ocn_bmb(:,:,:)

    bmb%z_lake_in = cmn%z_lake_bmb
    bmb%t_lake_in(:,:,:) = cmn%t_lake_bmb(:,:,:)
    bmb%s_lake_in(:,:,:) = cmn%s_lake_bmb(:,:,:)
    bmb%mask_lake_in(:,:) = cmn%mask_lake_bmb(:,:)

    return

  end subroutine cmn_to_bmb


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  b m b _ t o _ c m n
  ! Purpose    :  transfer bmb to common grid 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_to_cmn(bmb,cmn)

    implicit none

    type(bmb_class) :: bmb
    type(cmn_class) :: cmn




    return

  end subroutine bmb_to_cmn



  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  i c e _ t o _ b m b
  ! Purpose    :  transfer ice to bmb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_to_bmb(ice,bmb)

    implicit none

    type(ice_class) :: ice
    type(bmb_class) :: bmb


    ! ice base elevation
    bmb%zb = ice%z_base
    ! filtered lithosphere elevation
    bmb%zl_fil = ice%z_bed_fil
    ! ocean/lake mask
    bmb%mask_ocn_lake = ice%mask_ocn_lake
    ! ice shelf mask
    where (ice%h_ice.gt.0._wp .and. ice%h_ice.lt.(-ice%z_bed*rho_sw/rho_i) .and. ice%z_bed.lt.0._wp)
      bmb%mask_ice_shelf = 1
    elsewhere
      bmb%mask_ice_shelf = 0
    endwhere


    return

  end subroutine ice_to_bmb


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g e o _ t o _ b m b
  ! Purpose    :  transfer geo to bmb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_to_bmb(geo,bmb)

    implicit none

    type(geo_class) :: geo
    type(bmb_class) :: bmb

    type(map_scrip_class) :: maps_geo_to_bmb

    real(wp), allocatable :: h_ice(:,:)
    real(wp), allocatable :: z_bed(:,:)


    ! no ice model, use ice sheets from bnd
    allocate(h_ice(bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(z_bed(bmb%grid%G%nx,bmb%grid%G%ny))
    ! initialize map
    call map_scrip_init(maps_geo_to_bmb,geo%hires%grid,bmb%grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)
    ! map ice base elevation
    call map_scrip_field(maps_geo_to_bmb,"zb",geo%hires%z_topo-geo%hires%h_ice,bmb%zb,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5._dp*bmb%grid%G%dx,bmb%grid%G%dx])
    ! map ice thickness and bedrock elevation
    call map_scrip_field(maps_geo_to_bmb,"h_ice",geo%hires%h_ice,h_ice,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5._dp*bmb%grid%G%dx,bmb%grid%G%dx])
    where (h_ice<10._wp) h_ice = 0._wp
      call map_scrip_field(maps_geo_to_bmb,"z_bed",geo%hires%z_bed,z_bed,method="mean",missing_value=-9999._dp, &
        filt_method="none",filt_par=[5._dp*bmb%grid%G%dx,bmb%grid%G%dx])

    ! generate ice shelf mask
    where (h_ice.gt.0._wp .and. h_ice.lt.(-z_bed*rho_sw/rho_i) .and. z_bed.lt.0._wp)
      bmb%mask_ice_shelf = 1
    elsewhere
      bmb%mask_ice_shelf = 0
    endwhere
    bmb%zl_fil = z_bed
    deallocate(h_ice)
    deallocate(z_bed)


    return

  end subroutine geo_to_bmb


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  b m b _ t o _ i c e
  ! Purpose    :  transfer bmb to ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_to_ice(bmb,ice)

    implicit none

    type(bmb_class) :: bmb
    type(ice_class) :: ice


    ! basal mass balance of floating ice
    ice%bmb_shlf = bmb%bmb_ann/sec_year / rho_i ! kg/m2/a -> m(ice equivalent)/s

    ! ocean temperature for small-scale discharge parameterisation
    ice%t_ocn = bmb%t_disc     ! degC
    ! ocean salinity for small-scale discharge parameterisation
    ice%s_ocn = bmb%s_disc     ! psu


    return

  end subroutine bmb_to_ice


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  i c e _ t o _ c m n
  ! Purpose    :  transfer ice variables to common grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_to_cmn(ice,cmn)

    implicit none

    type(ice_class) :: ice
    type(cmn_class) :: cmn

    integer :: i, j, ii, jj

    real(wp), allocatable, dimension(:,:) :: calving_ice_i
    real(wp), allocatable, dimension(:,:) :: bmelt_grd_i
    real(wp), allocatable, dimension(:,:) :: bmelt_flt_i


    ! aggregate calving and basal melt from ice sheets to coupler grid

    ! initialize to 0 only where ice domain covers coupler grid (to avoid overwriting possible other ice sheet domains)
    allocate(calving_ice_i(cmn%grid%G%nx,cmn%grid%G%ny))
    allocate(bmelt_grd_i(cmn%grid%G%nx,cmn%grid%G%ny))
    allocate(bmelt_flt_i(cmn%grid%G%nx,cmn%grid%G%ny))
    do i=1,cmn%grid%G%nx
      do j=1,cmn%grid%G%ny
        if (ice%grid_ice_to_cmn%ncells(i,j)>0) then
          calving_ice_i(i,j) = 0._wp
          bmelt_grd_i(i,j) = 0._wp
          bmelt_flt_i(i,j) = 0._wp
        endif
      enddo
    enddo

    ! integrate over grid cell ice sheet area
    !!$omp parallel do private(ii,jj,i,j)
    do ii=1,ice%grid%G%nx
      do jj=1,ice%grid%G%ny
        i = ice%grid_ice_to_cmn%i_lowres(ii,jj)
        j = ice%grid_ice_to_cmn%j_lowres(ii,jj)
        calving_ice_i(i,j) = calving_ice_i(i,j) + &
          ice%calv(ii,jj)*ice%grid%area(ii,jj)*1.e6_wp*rho_i  ! m(ice)/s * m2 * kg(ice)/m3(ice) = kg(ice)/s = kg(water)/s
        if (ice%h_ice(ii,jj).gt.0._wp .and. ice%h_ice(ii,jj).lt.(-ice%z_bed(ii,jj)*rho_sw/rho_i) .and. ice%z_bed(ii,jj).lt.0._wp) then
          ! floating ice
          bmelt_flt_i(i,j) = bmelt_flt_i(i,j) + &
            ice%Q_b(ii,jj)*ice%grid%area(ii,jj)*1.e6_wp*rho_i  ! m(ice)/s * m2 * kg(ice)/m3(ice) = kg(ice)/s = kg(water)/s
        else
          ! grounded ice
          bmelt_grd_i(i,j) = bmelt_grd_i(i,j) + &
            ice%Q_b(ii,jj)*ice%grid%area(ii,jj)*1.e6_wp*rho_i  ! m(ice)/s * m2 * kg(ice)/m3(ice) = kg(ice)/s = kg(water)/s
        endif
      enddo
    enddo
    !!$omp end parallel do

    ! smooth calving and basal melt in time 
    do i=1,cmn%grid%G%nx
      do j=1,cmn%grid%G%ny
        if (ice%grid_ice_to_cmn%ncells(i,j)>0) then
          cmn%calving_ice_i(i,j) = relax_calv*cmn%calving_ice_i(i,j) + (1._wp-relax_calv)*calving_ice_i(i,j)
          cmn%bmelt_grd_i(i,j) = relax_bmelt*cmn%bmelt_grd_i(i,j) + (1._wp-relax_bmelt)*bmelt_grd_i(i,j)
          cmn%bmelt_flt_i(i,j) = relax_bmelt*cmn%bmelt_flt_i(i,j) + (1._wp-relax_bmelt)*bmelt_flt_i(i,j)
        endif
      enddo
    enddo

    deallocate(calving_ice_i)
    deallocate(bmelt_grd_i)
    deallocate(bmelt_flt_i)

  return

  end subroutine ice_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  i c e _ t o _ g e o
  ! Purpose    :  transfer ice to geo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_to_geo(n_ice_domain,ice,geo)

    implicit none

    integer, intent(in) :: n_ice_domain
    type(ice_class) :: ice(:)
    type(geo_class) :: geo

    integer :: i, j, n
    real(wp) :: V_grounded, V_gr_redu, V_ice_af
    real(wp), save :: V_ice_af_old
    real(wp), allocatable, dimension(:,:) :: mask_ice
    real(wp), allocatable, dimension(:,:) :: mask_ice_geo
    type(map_scrip_class), save, allocatable, dimension(:) :: maps_hice_to_geo
    type(map_scrip_class), save, allocatable, dimension(:) :: maps_mask_to_geo

    logical, save :: firstcall = .true.


    !----------------------------------------------------------------
    ! initialize map
    !----------------------------------------------------------------
    if (firstcall) then
      allocate(maps_hice_to_geo(n_ice_domain))
      allocate(maps_mask_to_geo(n_ice_domain))
      do n=1,n_ice_domain
        call map_scrip_init(maps_hice_to_geo(n),ice(n)%grid,geo%hires%grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)
        call map_scrip_init(maps_mask_to_geo(n),ice(n)%grid,geo%hires%grid,method="nn",fldr="maps",load=.TRUE.,clean=.FALSE.)
      enddo
    endif

    !----------------------------------------------------------------
    ! map from ice domains to geo domain
    !----------------------------------------------------------------
    !!$omp parallel do
    do n=1,n_ice_domain
      ! ice thickness
      call map_scrip_field(maps_hice_to_geo(n),"hice",ice(n)%H_ice,geo%hires%h_ice,method="mean",missing_value=-9999._dp,reset=.false., &
        filt_method="none",filt_par=[5._dp*geo%hires%grid%G%dx,geo%hires%grid%G%dx])
      allocate(mask_ice(ice(n)%grid%G%nx,ice(n)%grid%G%ny))
      allocate(mask_ice_geo(geo%hires%grid%G%nx,geo%hires%grid%G%ny))
      where (ice(n)%H_ice>h_ice_min) 
        mask_ice = 1.
      elsewhere
        mask_ice = 0.
      endwhere
      mask_ice_geo = 1.
      call map_scrip_field(maps_hice_to_geo(n),"mask",mask_ice,mask_ice_geo,method="mean",missing_value=-9999._dp,reset=.false., &
        filt_method="none",filt_par=[5._dp*geo%hires%grid%G%dx,geo%hires%grid%G%dx])
      where (mask_ice_geo<0.5) geo%hires%h_ice = 0._wp
      !mask_ice_geo = 0._wp
      !call map_scrip_field(maps_hice_to_geo(n),"mask",mask_ice,mask_ice_geo,method="mean",missing_value=-9999._wp,reset=.false., &
      !  filt_method="none",filt_par=[5.0*geo%hires%grid%G%dx,geo%hires%grid%G%dx])
      !where (mask_ice_geo>0.5_wp) mask_ice_geo = 1._wp
      deallocate(mask_ice)
      deallocate(mask_ice_geo)
      !print *,'ice vol ice',sum(ice(n)%H_ice(:,:)*ice(n)%grid%area(:,:)) * 1.e6_wp ! m3
      !print *,'ice vol geo',sum(geo%hires%h_ice(:,:)*geo%hires%grid%area(:,:), mask=mask_ice_geo==1._wp) * 1.e6_wp ! m3

    enddo
    !!$omp end parallel do

!    where (geo%hires%grid%lat<45._wp .and. geo%hires%grid%lat>0._wp) 
!      geo%hires%h_ice = 0._wp
!    endwhere

    !----------------------------------------------------------------
    ! compute sea level change from ice volume changes
    !----------------------------------------------------------------

    ! compute total ice volume of all domains
    V_grounded = 0._wp
    V_gr_redu = 0._wp

    do n=1,n_ice_domain ! loop over all ice domains

      do i=1, ice(n)%grid%G%nx
        do j=1, ice(n)%grid%G%ny
          if (ice(n)%H_ice(i,j).gt.0._wp .and. ice(n)%H_ice(i,j).gt.(-ice(n)%z_bed(i,j)*rho_sw/rho_i)) then   ! grounded ice
            ! total grounded ice
            V_grounded = V_grounded + ice(n)%H_ice(i,j) * ice(n)%grid%area(i,j) * 1.e6_wp ! m3
            ! grounded ice below floatation
            V_gr_redu  = V_gr_redu + rho_sw/rho_i*max((ice(n)%z_sl(i,j)-ice(n)%z_bed(i,j)),0.0_wp) &
                                                                  *ice(n)%grid%area(i,j) * 1.e6_wp ! m3
          endif
        enddo
      enddo

    enddo

    ! total ice volume above floatation (which contributes to sea level)
    V_ice_af    = V_grounded - V_gr_redu

    if (firstcall) then
      if (geo_restart) then
        V_ice_af_old = geo%V_ice_af
      else
        V_ice_af_old = V_ice_af
      endif
    endif

    geo%V_ice_af = V_ice_af

    ! change in sea level (m)
    geo%d_sea_level = -(V_ice_af-V_ice_af_old)*(rho_i/rho_w)/geo%ocn_area_tot     ! m^3 ice equiv./m^2 -> m water equiv. -> sea level equivalent

    ! save ice volume above floatation
    V_ice_af_old = V_ice_af

    if (firstcall) firstcall = .false.

    return

  end subroutine ice_to_geo


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g e o _ t o _ i c e
  ! Purpose    :  transfer geo to ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_to_ice(n_ice_domain,geo,ice)

    implicit none

    integer, intent(in) :: n_ice_domain
    type(geo_class) :: geo
    type(ice_class) :: ice(:)

    integer :: n 
    type(map_scrip_class), save, allocatable, dimension(:) :: maps_geo_to_ice
    type(map_scrip_class), save, allocatable, dimension(:) :: maps_geo_to_ice_nn

    logical, save :: firstcall = .true.


    !----------------------------------------------------------------
    ! initialize map
    !----------------------------------------------------------------
    if (firstcall) then
      allocate(maps_geo_to_ice(n_ice_domain))
      allocate(maps_geo_to_ice_nn(n_ice_domain))
      do n=1,n_ice_domain
        call map_scrip_init(maps_geo_to_ice(n),geo%hires%grid,ice(n)%grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)
        call map_scrip_init(maps_geo_to_ice_nn(n),geo%hires%grid,ice(n)%grid,method="nn",fldr="maps",load=.TRUE.,clean=.FALSE.)
      enddo
    endif

    !----------------------------------------------------------------
    ! map from geo to ice domains
    !----------------------------------------------------------------
    do n=1,n_ice_domain

      !$omp parallel sections
      !$omp section
      ! map ocean/lake mask to ice grid
      call map_scrip_field(maps_geo_to_ice_nn(n),"mask_ocn_lake",geo%hires%mask_ocn_lake,ice(n)%mask_ocn_lake,method="mean",missing_value=-9999._dp, &
        filt_method="none",filt_par=[100._dp,ice(n)%grid%G%dx])
      !$omp section
      ! map sea (lake) level to ice grid, sea level over ocean is ==0
      call map_scrip_field(maps_geo_to_ice(n),"z_sl",geo%hires%z_sea_lake,ice(n)%z_sl,method="mean",missing_value=-9999._dp, &
        filt_method="gaussian",filt_par=[100._dp,ice(n)%grid%G%dx])
      !$omp section
      ! map lithosphere elevation relative to contemporary sea level (geoid)
      call map_scrip_field(maps_geo_to_ice(n),"zbed",geo%hires%z_bed,ice(n)%z_bed,method="mean",missing_value=-9999._dp, &
        filt_method="none",filt_par=[5._dp*ice(n)%grid%G%dx,ice(n)%grid%G%dx])
      !$omp section
      ! map lithosphere elevation relative to contemporary sea level (geoid), filtered
      call map_scrip_field(maps_geo_to_ice(n),"zbed",geo%hires%z_bed,ice(n)%z_bed_fil,method="mean",missing_value=-9999._dp, &
        filt_method="gaussian",filt_par=[100._dp,ice(n)%grid%G%dx])
      !$omp section
      ! map geothermal heat flux to ice grid
      call map_scrip_field(maps_geo_to_ice(n),"q_geo",geo%hires%q_geo,ice(n)%q_geo,method="mean",missing_value=-9999._dp, &
        filt_method="gaussian",filt_par=[100._dp,ice(n)%grid%G%dx])
      !$omp section
      ! map sediment mask to ice grid
      call map_scrip_field(maps_geo_to_ice(n),"h_sed",geo%hires%h_sed,ice(n)%H_sed,method="mean",missing_value=-9999._dp, &
        filt_method="none",filt_par=[100._dp,ice(n)%grid%G%dx])
      !$omp end parallel sections

    enddo

    if (firstcall) firstcall = .false.

    return

  end subroutine geo_to_ice


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  b n d _ t o _ g e o
  ! Purpose    :  transfer bnd to geo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bnd_to_geo(bnd,geo)

    implicit none

    type(bnd_class) :: bnd
    type(geo_class) :: geo

    real(wp), allocatable, dimension(:,:) :: mask_ice
    real(wp), allocatable, dimension(:,:) :: mask_ice_geo
    type(map_scrip_class), save :: maps_bndice_to_geo
    type(map_scrip_class), save :: maps_bndgeo_to_geo

    real(wp), dimension(:,:), allocatable, save :: z_bed_anom

    logical, save :: firstcall = .true.
    logical, save :: l_samegrid_geo_bndgeo
    logical, save :: l_samegrid_geo_bndice


    if (firstcall) then

      if (.not.flag_ice .or. ifake_ice.eq.1) then
        if (bnd%ice%grid%name==geo%hires%grid%name) then
          ! grids are equal
          l_samegrid_geo_bndice = .true.
        else
          ! grids are different
          l_samegrid_geo_bndice = .false.
          ! initialize map
          call map_scrip_init(maps_bndice_to_geo,bnd%ice%grid,geo%hires%grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)
        endif
      endif

      if (.not.flag_geo) then
        allocate(z_bed_anom(geo%hires%grid%G%nx,geo%hires%grid%G%ny))
        if (bnd%geo%grid%name==geo%hires%grid%name) then
          ! grids are equal
          l_samegrid_geo_bndgeo = .true.
        else
          ! grids are different
          l_samegrid_geo_bndgeo = .false.
          ! initialize map
          call map_scrip_init(maps_bndgeo_to_geo,bnd%geo%grid,geo%hires%grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)
        endif
      endif

    endif

    ! if required, get ice thickness from bnd
    if (.not.flag_ice .or. ifake_ice.eq.1) then
      if (l_samegrid_geo_bndice) then
        geo%hires%h_ice = bnd%ice%h_ice
      else
        call map_scrip_field(maps_bndice_to_geo,"hice",bnd%ice%h_ice,geo%hires%h_ice,method="mean",missing_value=-9999._dp,reset=.false., &
          filt_method="none",filt_par=[5._dp*geo%hires%grid%G%dx,geo%hires%grid%G%dx])
        allocate(mask_ice(bnd%ice%grid%G%nx,bnd%ice%grid%G%ny))
        allocate(mask_ice_geo(geo%hires%grid%G%nx,geo%hires%grid%G%ny))
        where (bnd%ice%h_ice>h_ice_min) 
          mask_ice = 1.
        elsewhere
          mask_ice = 0.
        endwhere
        mask_ice_geo = 1.
        call map_scrip_field(maps_bndice_to_geo,"mask",mask_ice,mask_ice_geo,method="mean",missing_value=-9999._dp,reset=.false., &
          filt_method="none",filt_par=[5._dp*geo%hires%grid%G%dx,geo%hires%grid%G%dx])
        where (mask_ice_geo<0.5) geo%hires%h_ice = 0._wp
        deallocate(mask_ice)
        deallocate(mask_ice_geo)
      endif
    endif

    ! if solid Earth model not active, derive bedrock elevation anomaly from bnd
    if (.not.flag_geo) then
      if (l_samegrid_geo_bndgeo) then
        geo%hires%z_bed = bnd%geo%z_bed
        bnd%geo%z_bed_ref = bnd%geo%z_bed
      else
        ! map bedrock topography to geo
        call map_scrip_field(maps_bndgeo_to_geo,"zbed",bnd%geo%z_bed-bnd%geo%z_bed_ref,z_bed_anom,method="mean",missing_value=-9999._dp, &
          filt_method="none",filt_par=[5._dp*geo%hires%grid%G%dx,geo%hires%grid%G%dx])
        geo%hires%z_bed = geo%hires%z_bed_ref + z_bed_anom 
      endif
    endif

    if (firstcall) firstcall = .false.

    return

  end subroutine bnd_to_geo


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c m n _ t o _ g e o
  ! Purpose    :  transfer cmn to geo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_to_geo(cmn,geo)

    implicit none

    type(cmn_class) :: cmn
    type(geo_class) :: geo

 
    ! sea level
    if (.not.flag_geo) then
      geo%sea_level = cmn%sea_level
    endif

    if (flag_lnd .and. flag_lakes) then
      ! update volume of lakes 
      geo%lake%vol = cmn%lake%vol
      ! lake depth changes for diagnostic
      geo%lake%dh_p_e = cmn%lake%dh_p_e
      geo%lake%dh_run = cmn%lake%dh_run
      geo%lake%runoff_diag = cmn%lake%runoff_diag
    endif
    

    return

  end subroutine cmn_to_geo


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g e o _ t o _ c m n
  ! Purpose    :  transfer geo to cmn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_to_cmn(geo,cmn)

    implicit none

    type(geo_class), intent(in) :: geo
    type(cmn_class) :: cmn

    integer :: k


    ! sea level
    cmn%sea_level = geo%sea_level

    ! Bering Strait cross-sectional area
    cmn%A_bering = geo%A_bering

    ! ocean fraction of grid cell, including floating ice
    cmn%f_ocn   = geo%f_ocn
    ! ocean fraction of grid cell, excluding floating ice
    cmn%f_ocn2  = geo%f_ocn2
    ! update ocean mask
    where (cmn%f_ocn.gt.0._wp)
      cmn%mask_ocn = 1
    elsewhere
      cmn%mask_ocn = 0
    endwhere
    cmn%mask_coast = geo%mask_coast
    cmn%mask_coast2 = geo%mask_coast2
    ! total ocean volume
    cmn%ocn_vol_tot = geo%ocn_vol_tot

    ! reference land fraction for present day sea level
    cmn%f_lnd0 = geo%f_lnd0

    ! land fraction of grid cell
    cmn%f_lnd   = geo%f_lnd

    ! ice sheet fraction of grid cell
    cmn%f_ice     = geo%f_ice
    cmn%f_ice_grd = geo%f_ice_grd
    cmn%f_ice_flt = geo%f_ice_flt

    ! lake fraction
    cmn%f_lake = geo%f_lake

    ! update land mask
    where (cmn%f_lnd+geo%f_ice_flt.gt.0._wp)
      cmn%mask_lnd = 1
    elsewhere
      cmn%mask_lnd = 0
    endwhere

    ! runoff directions
    cmn%i_runoff = geo%i_runoff  
    cmn%j_runoff = geo%j_runoff
    cmn%i_runoff_veg = geo%i_runoff_veg  
    cmn%j_runoff_veg = geo%j_runoff_veg
    cmn%i_runoff_ice = geo%i_runoff_ice 
    cmn%j_runoff_ice = geo%j_runoff_ice

    ! number of lakes
    cmn%n_lakes = geo%n_lakes

    if (flag_lakes) then
      ! lake type
      if (allocated(cmn%lake)) deallocate(cmn%lake)
      allocate(cmn%lake(cmn%n_lakes))
      cmn%lake = geo%lake
      ! save lake volume
      cmn%lake%vol_old =  cmn%lake%vol
      ! individual lake fractions
      if (allocated(cmn%f_lake_n)) deallocate(cmn%f_lake_n)
      allocate(cmn%f_lake_n(cmn%n_lakes,ni,nj))
      cmn%f_lake_n = geo%f_lake_n
    endif

    ! drainage fractions to lakes and ocean
    if (allocated(cmn%f_drain_veg)) deallocate(cmn%f_drain_veg)
    if (allocated(cmn%f_drain_ice)) deallocate(cmn%f_drain_ice)
    allocate(cmn%f_drain_veg(0:cmn%n_lakes,ni,nj))
    allocate(cmn%f_drain_ice(0:cmn%n_lakes,ni,nj))
    cmn%f_drain_veg = geo%f_drain_veg
    cmn%f_drain_ice = geo%f_drain_ice

    ! atlantic catchment area
    cmn%idivide_pac_atl = geo%idivide_pac_atl
    cmn%idivide_atl_indpac = geo%idivide_atl_indpac

    ! surface elevation on common grid
    cmn%z_sur = geo%z_sur
    cmn%z_ocn_max = geo%z_ocn_max_q
    cmn%z_veg = geo%z_veg
    cmn%z_veg_min = geo%z_veg_min
    cmn%z_veg_max = geo%z_veg_max
    cmn%z_ice = geo%z_ice
    cmn%z_lake = geo%z_lake
    cmn%z_sur_n(:,:,i_surf_macro_ocn) = 0._wp
    cmn%z_sur_n(:,:,i_surf_macro_sic) = 0._wp
    cmn%z_sur_n(:,:,i_surf_macro_lnd) = geo%z_veg
    cmn%z_sur_n(:,:,i_surf_macro_ice) = geo%z_ice
    cmn%z_sur_n(:,:,i_surf_macro_lake) = geo%z_lake

    cmn%z_sur_smooth_std = geo%z_sur_smooth_std
    cmn%z_veg_std = geo%z_veg_std

    cmn%coast_nbr   = geo%coast_nbr
    cmn%i_coast_nbr = geo%i_coast_nbr
    cmn%j_coast_nbr = geo%j_coast_nbr

    ! for corals
    if (flag_bgc) then
      do k=-250,50
        where (cmn%f_ocn.gt.0._wp) 
          cmn%coral_f_area(:,:,k) = geo%coral_f_area(:,:,k) / cmn%f_ocn
          cmn%coral_f_topo(:,:,k) = geo%coral_f_topo(:,:,k)
        elsewhere
          cmn%coral_f_area(:,:,k) = 0._wp
          cmn%coral_f_topo(:,:,k) = 0._wp
        endwhere
      enddo
    endif

    ! geothermal heat flux
    cmn%q_geo = geo%q_geo

    return

  end subroutine geo_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  b n d _ t o _ c m n
  ! Purpose    :  transfer boundary data to common grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bnd_to_cmn(bnd,cmn,ini)

    implicit none

    type(bnd_class) :: bnd
    type(cmn_class) :: cmn
    integer, optional :: ini

    integer :: i, j, ii, jj, n, k_bmb
    real(wp) :: tacl, ftemp, frcold
    real(wp), parameter :: tshift = 5._wp

    real(wp), allocatable, dimension(:,:,:) :: melt_ice_i_mon


    if (time_soy_bnd .or. present(ini)) then 

      ! sea level
      if (.not.flag_geo) then
        cmn%sea_level = bnd%sea_level
      endif

      ! orbital parameters
      cmn%eccentricity = bnd%eccentricity
      cmn%precession   = bnd%precession
      cmn%obliquity    = bnd%obliquity

      cmn%solar  = bnd%solar
      cmn%solarm = bnd%solarm
      cmn%solarmin = bnd%solarmin
      cmn%solarmax = bnd%solarmax
      cmn%cosz   = bnd%cosz
      cmn%coszm  = bnd%coszm
      cmn%daylength = bnd%daylength

      if (.not.flag_co2) then
        cmn%co2 = bnd%co2
        if (ico2_rad.eq.0) then
          cmn%co2_rad = cmn%co2
        endif
        cmn%c13_c12_atm = (bnd%d13c_atm/1000._wp+1._wp)*c13_c12_std      
        cmn%c14_c_atm = (bnd%D14c_atm/1000._wp+1._wp)*c14_c_std / ((0.975_wp/(1._wp+bnd%d13c_atm/1000._wp))**2)
      endif
      if (ico2_rad.ne.0) then
        cmn%co2_rad = bnd%co2_rad
      endif

      if (.not.flag_ch4) then
        cmn%ch4 = bnd%ch4
        if (ich4_rad.eq.0) then
          cmn%ch4_rad = cmn%ch4
        endif
      endif
      if (ich4_rad.ne.0) then
        cmn%ch4_rad = bnd%ch4_rad
      endif

      cmn%n2o = bnd%n2o

      cmn%cfc11 = bnd%cfc11
      cmn%cfc12 = bnd%cfc12

      cmn%so4(:,:) = bnd%so4(:,:)

      cmn%o3_pl(:) = bnd%o3_pl(:)
      cmn%o3(:,:,:) = bnd%o3(:,:,:)

      ! land use state
      cmn%f_crop(:,:)    = bnd%f_crop(:,:)
      cmn%f_pasture(:,:) = bnd%f_pasture(:,:)

      ! disturbance rate
      cmn%disturbance(:,:,:) = bnd%disturbance(:,:,:)

      if (ifake_ice.eq.1) then

        ! derive icemelt flux going into runoff and accumulation to be substracted from calving in the land model
        ! from changes in thickness of the prescribed ice sheets
        allocate(melt_ice_i_mon(cmn%grid%G%nx,cmn%grid%G%ny,nmon_year))
        melt_ice_i_mon(:,:,:) = 0._wp
        cmn%acc_ice_i_mon(:,:,:) = 0._wp
        do ii=1,bnd%ice%grid%G%nx
          do jj=1,bnd%ice%grid%G%ny
            i = bnd%ice%grid_ice_to_cmn%i_lowres(ii,jj)
            j = bnd%ice%grid_ice_to_cmn%j_lowres(ii,jj)
            if (bnd%ice%h_ice(ii,jj).gt.0._wp) then
              melt_ice_i_mon(i,j,:) = melt_ice_i_mon(i,j,:) + &
                max(0._wp,(-bnd%ice%dh_ice_dt(ii,jj))*rho_i*bnd%ice%grid%area(ii,jj))  ! m /s * kg/m3 * m2 = kg/s
              cmn%acc_ice_i_mon(i,j,:) = cmn%acc_ice_i_mon(i,j,:) + &
                max(0._wp,(bnd%ice%dh_ice_dt(ii,jj))*rho_i*bnd%ice%grid%area(ii,jj))  ! m /s * kg/m3 * m2 = kg/s
            endif
          enddo
        enddo

        ! smooth in time 
        do i=1,cmn%grid%G%nx
          do j=1,cmn%grid%G%ny
            cmn%melt_ice_i_mon(i,j,:) = relax_run*cmn%melt_ice_i_mon(i,j,:) + (1._wp-relax_run)*melt_ice_i_mon(i,j,:)
          enddo
        enddo

        deallocate(melt_ice_i_mon)

      endif

    endif

    ! ocean relaxation setup, top layer temperature and salinity to restore to
    if (flag_ocn .and. ocn_restore_sal) then
      cmn%sss_dat = bnd%ocn%s1l ! psu
    endif
    if (flag_ocn .and. ocn_restore_temp) then
      cmn%sst_dat = bnd%ocn%t1l ! K
    endif

    ! if no ocean
    if (.not. flag_ocn .and. .not.l_spinup_cc) then
      if (flag_sic) then
        cmn%sst = bnd%ocn%sst ! K, sea surface temperature
        cmn%sss = bnd%ocn%sss ! psu, sea surface salinity
        cmn%uo1 = bnd%ocn%uo1 ! m/s, zonal component of surface ocean current on u-grid
        cmn%vo1 = bnd%ocn%vo1 ! m/s, meridional component of surface ocean current on v-grid
      endif
      if (flag_lnd .and. .not.l_spinup_cc) then
        cmn%t_shelf = bnd%ocn%t_shelf  ! K         
      endif
      if (flag_bmb) then
        k_bmb = minloc(abs(bnd%ocn%depth-1000._wp),1)
        if (.not.allocated(cmn%z_ocn_bmb)) allocate(cmn%z_ocn_bmb(k_bmb))
        if (.not.allocated(cmn%t_ocn_bmb)) allocate(cmn%t_ocn_bmb(ni,nj,k_bmb))
        if (.not.allocated(cmn%s_ocn_bmb)) allocate(cmn%s_ocn_bmb(ni,nj,k_bmb))
        cmn%z_ocn_bmb = bnd%ocn%depth(1:k_bmb)
        cmn%t_ocn_bmb = bnd%ocn%t(:,:,1:k_bmb)  ! K         
        cmn%s_ocn_bmb = bnd%ocn%s(:,:,1:k_bmb)  ! psu
      endif
    endif

    ! if no sea ice
    if (.not. flag_sic .and. .not.l_spinup_cc) then
      cmn%f_sic = bnd%sic%f_sic
    endif

    ! if no atmosphere
    if (.not. flag_atm .and. .not.l_spinup_cc) then
      do n=1,nsurf_macro
        cmn%taux(:,:,n) = bnd%atm%taux(:,:)
        cmn%tauy(:,:,n) = bnd%atm%tauy(:,:)
      enddo
      cmn%tauxo(:,:) = bnd%atm%taux(:,:)
      cmn%tauyo(:,:) = bnd%atm%tauy(:,:)
      cmn%cld      = bnd%atm%cld
      cmn%slp      = bnd%atm%slp
      cmn%swd  = bnd%atm%swdown
      cmn%swnet(:,:,i_surf_macro_sic) = bnd%atm%swdown * (1._wp-0.8_wp)  ! assuming fixed sea ice albedo
      cmn%swnet(:,:,i_surf_macro_ocn) = bnd%atm%swdown * (1._wp-0.06_wp) ! assuming fixed ocean albedo
      cmn%usur = bnd%atm%usur
      cmn%vsur = bnd%atm%vsur
      cmn%u700 = bnd%atm%usur
      cmn%v700 = bnd%atm%vsur
      cmn%tam = bnd%atm%tair
      do n=1,nsurf
        cmn%ps(:,:,n)   = bnd%atm%pressure
        cmn%t2m(:,:,n)  = bnd%atm%tair
        cmn%q2m(:,:,n)  = bnd%atm%qair
        cmn%lwd(:,:,n)  = bnd%atm%lwdown
        cmn%wind(:,:,n)     = bnd%atm%wind
      enddo
      do n=1,nsurf_macro
        cmn%lwd_cld(:,:,n) = bnd%atm%lwdown
        cmn%lwd_cs(:,:,n)  = bnd%atm%lwdown
      enddo
      cmn%swd_vis_dir(:,:)  = bnd%atm%swdown_c    ! clear sky
      cmn%swd_nir_dir(:,:)  = bnd%atm%swdown_c    ! clear sky
      cmn%swd_vis_dif(:,:)  = bnd%atm%swdown_d    ! diffuse
      cmn%swd_nir_dif(:,:)  = bnd%atm%swdown_d    ! diffuse
      cmn%dswd_dalb_vis_dir = 0._wp 
      cmn%dswd_dalb_nir_dir = 0._wp
      cmn%dswd_dalb_vis_dif = 0._wp
      cmn%dswd_dalb_nir_dif = 0._wp
      cmn%dswd_dz_nir_dir   = 0._wp
      cmn%dswd_dz_nir_dif   = 0._wp  
      do n=1,nsurf
        cmn%alb_vis_dir(:,:,n) = 0.1_wp ! dummy value for SEMI  
        cmn%alb_nir_dir(:,:,n) = 0.1_wp ! dummy value for SEMI
        cmn%alb_vis_dif(:,:,n) = 0.1_wp ! dummy value for SEMI
        cmn%alb_nir_dif(:,:,n) = 0.1_wp ! dummy value for SEMI
      enddo

      do j=1,nj
        do i=1,ni
          if (cmn%t2m(i,j,1).gt.T0) then
            cmn%ram(i,j) = cmn%q2m(i,j,1) / q_sat_w(cmn%t2m(i,j,1),cmn%ps(i,j,1)) ! fixme
          else
            cmn%ram(i,j) = cmn%q2m(i,j,1) / q_sat_i(cmn%t2m(i,j,1),cmn%ps(i,j,1))
          endif
        enddo
      enddo
      cmn%t2m_min_mon = bnd%atm%tair_min_mon
      if (prc_forcing.eq.0) then
        do j=1,nj
          do i=1,ni
            cmn%rain(i,j,:) = bnd%atm%rain(i,j)
            cmn%snow(i,j,:) = bnd%atm%snow(i,j)
            cmn%prc(i,j)    = bnd%atm%rain(i,j) + bnd%atm%snow(i,j)
          enddo
        enddo
      else if (prc_forcing.eq.1) then
        cmn%prc  = bnd%atm%prc
        ! split total precipitation between rain and snow based on atmospheric temperature, from CLIMBER-2, check with data again!
        do j=1,nj
          do i=1,ni
            tacl = bnd%atm%tair(i,j)-T0-1._wp*tshift
            ftemp = cos(pi*tacl/(2._wp*tshift))
            frcold = 0.5_wp*(1._wp-ftemp)
            if (tacl.gt.0._wp) frcold = 0._wp
            if (tacl.lt.(-2._wp*tshift)) frcold = 1._wp
            cmn%rain(i,j,:) = (1._wp-frcold)*bnd%atm%prc(i,j)
            cmn%snow(i,j,:) = frcold*bnd%atm%prc(i,j)
          enddo
        enddo
      endif
    endif

    ! if no dust model
    if (.not. flag_dust) then
      cmn%dust_dep  = bnd%dust%dust_dep
    endif

    if (flag_atm .and. atm_fix_tau) then
      cmn%taux_dat(:,:) = bnd%atm%taux(:,:)
      cmn%tauy_dat(:,:) = bnd%atm%tauy(:,:)
    endif

   return

  end subroutine bnd_to_cmn


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  l a k e s _ u p d a t e 
  ! Purpose    :  update lakes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lakes_update(cmn)

    implicit none

    type(cmn_class) :: cmn

    integer :: i, j, n
    logical, parameter :: flag_print = .false.
    integer, parameter :: n_lake_print = 6 
    real(wp) :: lake_n_cell_area, lake_area_scale


    if (time_soy_lnd) then
      cmn%lake(:)%dh_p_e = 0._wp
      cmn%lake(:)%dh_run = 0._wp
      cmn%lake(:)%runoff_diag = 0._wp
    endif

    ! update lake volume using surface lake water balance and runoff from land and ice to lake

    if (flag_print) then
      print *
      print *,'lake ',n_lake_print
      print *,'vol old',cmn%lake(n_lake_print)%vol_old*1e-6
      print *,'vol before',cmn%lake(n_lake_print)%vol*1e-6
      print *,'vol min',cmn%lake(n_lake_print)%vol_min*1e-6
      print *,'vol pot',cmn%lake(n_lake_print)%vol_pot*1e-6
      print *,'area ',cmn%lake(n_lake_print)%area*1e-6
      print *,'area_min ',cmn%lake(n_lake_print)%area_min*1e-6
      print *,'area_pot ',cmn%lake(n_lake_print)%area_pot*1e-6
    endif

    !$omp parallel do private(n,i,j,lake_area_scale,lake_n_cell_area)
    do n=1,cmn%n_lakes
      do j=1,nj
        do i=1,ni
!          if (flag_print) then
!            if (n.eq.n_lake_print .and. cmn%f_lake_n(n,i,j).gt.0._wp) then
!              print *,'surface water balance'
!              print *,'i,j',i,j
!              print *,'f_lake',cmn%f_lake_n(n,i,j)
!              print *,'lake P-E (kg/m2)',cmn%lake_p_e(i,j)*dt_lnd
!              print *,'dV (m3) = ',area(i,j)*cmn%f_lake_n(n,i,j)*cmn%lake_p_e(i,j) * dt_lnd / rho_w * 1e-6
!            endif
!            if (n.eq.n_lake_print .and. cmn%f_drain_veg(n,i,j).gt.0._wp) then
!              print *,'runoff from land'
!              print *,'i,j',i,j
!              print *,'f_drain_veg',cmn%f_drain_veg(n,i,j)
!              print *,'runoff (kg/m2)',cmn%runoff_veg(i,j)*dt_lnd
!              print *,'dV (m3) = ',cmn%runoff_veg(i,j) * dt_lnd * cmn%f_drain_veg(n,i,j) / rho_w * 1e-6
!            endif
!          endif
          ! contribution of runoff from land in lake catchment to lake volume changes
          cmn%lake(n)%vol = cmn%lake(n)%vol &
            + cmn%runoff_veg(i,j) * dt_lnd*real(n_accel,wp) * cmn%f_drain_veg(n,i,j) / rho_w ! kg/s * s * m3/kg = m3
          ! contribution of runoff from ice in lake catchment to lake volume changes
          cmn%lake(n)%vol = cmn%lake(n)%vol &
            + cmn%runoff_ice(i,j) * dt_lnd*real(n_accel,wp) * cmn%f_drain_ice(n,i,j) / rho_w ! kg/s * s * m3/kg = m3
          ! contribution of surface water balance to lake volume changes
          ! scaling factor to account for changing lake area as a result of net water balance
          if (cmn%lake(n)%vol_old.gt.0._wp) then
            lake_area_scale = cmn%lake(n)%vol/cmn%lake(n)%vol_old
            lake_area_scale = min(2._wp,lake_area_scale)
          else
            lake_area_scale = 1._wp
          endif
          lake_n_cell_area = lake_area_scale*area(i,j)*cmn%f_lake_n(n,i,j)    ! area of lake n in cell (i,j)
          !lake_n_cell_area = area(i,j)*cmn%f_lake_n(n,i,j)    ! area of lake n in cell (i,j)
          cmn%lake(n)%vol = cmn%lake(n)%vol &
            + cmn%lake_p_e(i,j) * dt_lnd*real(n_accel,wp) * lake_n_cell_area / rho_w ! kg/m2/s * s * m2 * m3/kg = m3   

          if (cmn%lake(n)%area.gt.0._wp) then
            ! average lake depth change from net surface water balance, diagnostic only 
            cmn%lake(n)%dh_p_e = cmn%lake(n)%dh_p_e &
              + cmn%lake_p_e(i,j) * lake_n_cell_area/cmn%lake(n)%area * dt_lnd*real(n_accel,wp) / rho_w ! kg/m2/s * s * m3/kg = m
            ! lake depth change from runoff input, diagnostic only 
            cmn%lake(n)%dh_run = cmn%lake(n)%dh_run &
              + cmn%runoff_veg(i,j) * cmn%f_drain_veg(n,i,j)/cmn%lake(n)%area * dt_lnd*real(n_accel,wp) / rho_w & ! kg/s / m2 * s * m3/kg = m
              + cmn%runoff_ice(i,j) * cmn%f_drain_ice(n,i,j)/cmn%lake(n)%area * dt_lnd*real(n_accel,wp) / rho_w   ! kg/s / m2 * s * m3/kg = m
!            if (flag_print) then
!              if (n.eq.n_lake_print .and. lake_n_cell_area.gt.0._wp) then
!                print *
!                print *,n,cmn%f_lake_n(n,i,j)
!                print *,cmn%lake_p_e(i,j) * lake_n_cell_area/cmn%lake(n)%area * dt_lnd / rho_w
!                print *,cmn%runoff_veg(i,j) * cmn%f_drain_veg(n,i,j)/cmn%lake(n)%area * dt_lnd / rho_w
!                print *,cmn%runoff_ice(i,j) * cmn%f_drain_ice(n,i,j)/cmn%lake(n)%area * dt_lnd / rho_w
!                print *,cmn%f_drain_veg(n,i,j)
!                print *,lake_n_cell_area,cmn%lake(n)%area
!                print *,doy,cmn%lake(n)%dh_p_e,cmn%lake(n)%dh_run
!              endif
!            endif
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do

    do n=1,cmn%n_lakes
      ! runoff from lakes to the ocean when the lake is full
      cmn%lake(n)%runoff = max(0._wp,cmn%lake(n)%vol-cmn%lake(n)%vol_pot) * rho_w / (dt_lnd*real(n_accel,wp))   ! m3 * kg/m3 / s = kg/s
      ! 'negative runoff' from ocean to lakes if lake is below minimum volume
      cmn%lake(n)%runoff = cmn%lake(n)%runoff + min(0._wp,cmn%lake(n)%vol-cmn%lake(n)%vol_min) * rho_w / (dt_lnd*real(n_accel,wp))   ! m3 * kg/m3 / s = kg/s
      ! add runoff originating from lake mask changes computed in geo 
      cmn%lake(n)%runoff = cmn%lake(n)%runoff + cmn%lake(n)%runoff_geo    ! kg/s
      ! annualy integrated runoff from lake, diagnostic only 
      cmn%lake(n)%runoff_diag = cmn%lake(n)%runoff_diag &
        + cmn%lake(n)%runoff/rho_w*1.e-6_wp * dt_lnd/sec_year  ! kg/s * m3/kg * mln m3/m3 * s = mln m3
      ! 'calving' from lakes to the ocean, TODO
      cmn%lake(n)%calving = 0._wp
      ! limit lake volume to be lower than potential volume
      cmn%lake(n)%vol = min(cmn%lake(n)%vol,cmn%lake(n)%vol_pot)
      ! limit lake volume to be larger than minimum volume
      cmn%lake(n)%vol = max(cmn%lake(n)%vol,cmn%lake(n)%vol_min)
    enddo

    if (flag_print) then
      if (cmn%lake(n_lake_print)%vol_old.gt.0._wp) then
        lake_area_scale = cmn%lake(n_lake_print)%vol/cmn%lake(n_lake_print)%vol_old
        lake_area_scale = min(2._wp,lake_area_scale)
      else
        lake_area_scale = 1._wp
      endif
      print *,'lake_area_scale',lake_area_scale
      print *,'dh_p_e',cmn%lake(n_lake_print)%dh_p_e
      print *,'dh_run',cmn%lake(n_lake_print)%dh_run
      print *,'runoff',cmn%lake(n_lake_print)%runoff/rho_w*1.e-6
      print *,'runoff_diag',cmn%lake(n_lake_print)%runoff_diag
      print *,'vol after',cmn%lake(n_lake_print)%vol*1e-6
      print *,'vol pot',cmn%lake(n_lake_print)%vol_pot*1e-6
    endif

!    if (time_soy_lnd .or. time_eoy_lnd) then
!    print*,'lake #, i, j, vol_min, vol_pot, vol, runoff, runoff_geo'
!    write(6,'(3i7,5f12.4)')(n,cmn%lake(n)%hires%i,cmn%lake(n)%hires%j,cmn%lake(n)%vol_min*1e-10,cmn%lake(n)%vol_pot*1e-10,cmn%lake(n)%vol*1e-10,cmn%lake(n)%runoff/rho_w*1e-6,cmn%lake(n)%runoff_geo/rho_w*1e-6,n=1,cmn%n_lakes)
!  endif

   return

  end subroutine lakes_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  r u n o f f _ m e r g e
  ! Purpose    :  merge runoff computed by different models 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine runoff_merge(cmn)

    implicit none

    type(cmn_class) :: cmn

    integer :: i, j


    do j=1,nj
      do i=1,ni

        if (cmn%mask_smb(i,j).eq.0) then
          ! grid cell not covered by smb model domain(s) 
          ! melt from prescribed ice sheet changes (if applicable)
          cmn%melt_ice(i,j) = scale_runoff_ice*cmn%melt_ice_i_mon(i,j,mon)
          ! use runoff computed in the land model + melt from prescribed ice sheet changes (if applicable)
          cmn%runoff_ice(i,j) = cmn%runoff_ice_l(i,j) + cmn%melt_ice(i,j)
        else
          ! grid cell covered by smb model domain(s), use monthly runoff and icemelt computed by smb 
          cmn%runoff_ice(i,j) = scale_runoff_ice*cmn%runoff_ice_i_mon(i,j,mon)
          cmn%melt_ice(i,j)   = scale_runoff_ice*cmn%melt_ice_i_mon(i,j,mon)  
        endif

        if (cmn%mask_ice(i,j).eq.0) then
          ! grid cell not covered by ice sheet model domain(s), use calving computed in the land model, no basal melt
          ! substract ice accumulation from prescribed ice sheet changes (if applicable)
          !cmn%calving_ice(i,j) = cmn%calving_ice_l(i,j)
          cmn%calving_ice(i,j) = max(0._wp, cmn%calving_ice_l(i,j) - scale_runoff_ice*cmn%acc_ice_i_mon(i,j,mon)) 
          cmn%bmelt_grd(i,j) = 0._wp
          cmn%bmelt_flt(i,j) = 0._wp
        else 
          ! grid cell covered by ice model domain(s), use calving and basal melt computed by ice sheet model
          cmn%calving_ice(i,j) = cmn%calving_ice_i(i,j)
          cmn%bmelt_grd(i,j) = cmn%bmelt_grd_i(i,j)
          cmn%bmelt_flt(i,j) = cmn%bmelt_flt_i(i,j)
        endif

      enddo
    enddo

   return

  end subroutine runoff_merge


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  r u n o f f _ t o _ o c n
  ! Purpose    :  transfer runoff to ocean domain
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine runoff_to_ocn(cmn)

    implicit none

    type(cmn_class) :: cmn

    integer :: i, j, ir, jr, ii, jj, n, nr, nbr
    real(wp) :: area_nbr


    ! initialize
    cmn%runoff_veg_o  = 0._wp
    cmn%calving_veg_o = 0._wp
    cmn%runoff_ice_o  = 0._wp
    cmn%melt_ice_o    = 0._wp
    cmn%calving_ice_o = 0._wp
    cmn%runoff_lake_o  = 0._wp
    cmn%calving_lake_o = 0._wp
    cmn%bmelt_grd_o    = 0._wp
    cmn%bmelt_flt_o    = 0._wp

    do j=1,nj
      do i=1,ni

        ! runoff and calving from ice-free and lake-free land 
        if (cmn%f_drain_veg(0,i,j)>0._wp) then
          ! index of runoff destination cells (i,j) -> (ir,jr)
          ir = cmn%i_runoff_veg(i,j)
          jr = cmn%j_runoff_veg(i,j)
          ! runoff
          ! compute total area of neighbors where runoff is distributed
          area_nbr = 0._wp
          nbr = min(cmn%coast_nbr(ir,jr),n_cells_dist_runoff)
          do n=1,nbr
            ii = cmn%i_coast_nbr(ir,jr,n)
            jj = cmn%j_coast_nbr(ir,jr,n)
            area_nbr = area_nbr + area(ii,jj)*cmn%f_ocn(ii,jj)  ! m2
          enddo
          ! transfer runoff to ocean domain
          do n=1,nbr
            ii = cmn%i_coast_nbr(ir,jr,n)
            jj = cmn%j_coast_nbr(ir,jr,n)
            cmn%runoff_veg_o(ii,jj) = cmn%runoff_veg_o(ii,jj) &
              + cmn%runoff_veg(i,j) * cmn%f_drain_veg(0,i,j) / area_nbr     ! kg/m2/s
          enddo
          ! 'calving'
          ! compute total area of neighbors where calving is distributed
          area_nbr = 0._wp
          nbr = min(cmn%coast_nbr(ir,jr),n_cells_dist_calving)
          do n=1,nbr
            ii = cmn%i_coast_nbr(ir,jr,n)
            jj = cmn%j_coast_nbr(ir,jr,n)
            area_nbr = area_nbr + area(ii,jj)*cmn%f_ocn(ii,jj)  ! m2
          enddo
          ! transfer calving to ocean domain
          do n=1,nbr
            ii = cmn%i_coast_nbr(ir,jr,n)
            jj = cmn%j_coast_nbr(ir,jr,n)
            cmn%calving_veg_o(ii,jj) = cmn%calving_veg_o(ii,jj) &
              + cmn%calving_veg(i,j) * cmn%f_drain_veg(0,i,j) / area_nbr    ! kg/m2/s
          enddo
        endif

        ! runoff and calving from ice-sheets
        if (cmn%f_drain_ice(0,i,j)>0._wp) then
          ! index of runoff destination cells (i,j) -> (ir,jr)
          ir = cmn%i_runoff_ice(i,j)
          jr = cmn%j_runoff_ice(i,j)
          ! runoff and basal melt
          ! compute total area of neighbors where runoff and basal melt are distributed
          area_nbr = 0._wp
          nbr = min(cmn%coast_nbr(ir,jr),n_cells_dist_runoff)
          do n=1,nbr
            ii = cmn%i_coast_nbr(ir,jr,n)
            jj = cmn%j_coast_nbr(ir,jr,n)
            area_nbr = area_nbr + area(ii,jj)*cmn%f_ocn(ii,jj)  ! m2
          enddo
          ! transfer runoff and basal melt to ocean domain
          do n=1,nbr
            ii = cmn%i_coast_nbr(ir,jr,n)
            jj = cmn%j_coast_nbr(ir,jr,n)
            cmn%runoff_ice_o(ii,jj) = cmn%runoff_ice_o(ii,jj) &
              + cmn%runoff_ice(i,j) * cmn%f_drain_ice(0,i,j) / area_nbr     ! kg/m2/s
            cmn%melt_ice_o(ii,jj)   = cmn%melt_ice_o(ii,jj) &
              + cmn%melt_ice(i,j) * cmn%f_drain_ice(0,i,j) / area_nbr     ! kg/m2/s
            cmn%bmelt_grd_o(ii,jj) = cmn%bmelt_grd_o(ii,jj) &
              + cmn%bmelt_grd(i,j) * cmn%f_drain_ice(0,i,j) / area_nbr    ! kg/m2/s
            cmn%bmelt_flt_o(ii,jj) = cmn%bmelt_flt_o(ii,jj) &
              + cmn%bmelt_flt(i,j) * cmn%f_drain_ice(0,i,j) / area_nbr    ! kg/m2/s
          enddo
          ! calving 
          ! compute total area of neighbors where calving is distributed
          area_nbr = 0._wp
          nbr = min(cmn%coast_nbr(ir,jr),n_cells_dist_calving)
          do n=1,nbr
            ii = cmn%i_coast_nbr(ir,jr,n)
            jj = cmn%j_coast_nbr(ir,jr,n)
            area_nbr = area_nbr + area(ii,jj)*cmn%f_ocn(ii,jj)  ! m2
          enddo
          ! transfer calving to ocean domain
          do n=1,nbr
            ii = cmn%i_coast_nbr(ir,jr,n)
            jj = cmn%j_coast_nbr(ir,jr,n)
            cmn%calving_ice_o(ii,jj) = cmn%calving_ice_o(ii,jj) &
              + cmn%calving_ice(i,j) * cmn%f_drain_ice(0,i,j) / area_nbr    ! kg/m2/s
          enddo
        endif
      enddo
    enddo

    ! runoff and 'calving' from lakes
    if (flag_lakes) then
      do n=1,cmn%n_lakes
        ! index of runoff destination cells (lake n) -> (ir,jr)
        ir = cmn%lake(n)%i_runoff
        jr = cmn%lake(n)%j_runoff
        ! compute total area of neighbors where runoff is distributed
        area_nbr = 0._wp
        nbr = min(cmn%coast_nbr(ir,jr),n_cells_dist_runoff)
        do nr=1,nbr
          ii = cmn%i_coast_nbr(ir,jr,nr)
          jj = cmn%j_coast_nbr(ir,jr,nr)
          area_nbr = area_nbr + area(ii,jj)*cmn%f_ocn(ii,jj)  ! m2
        enddo
        ! transfer runoff and calving to ocean domain
        do nr=1,nbr
          ii = cmn%i_coast_nbr(ir,jr,nr)
          jj = cmn%j_coast_nbr(ir,jr,nr)
          cmn%runoff_lake_o(ii,jj) = cmn%runoff_lake_o(ii,jj) &
            + cmn%lake(n)%runoff / area_nbr     ! kg/m2/s
          cmn%calving_lake_o(ii,jj) = cmn%calving_lake_o(ii,jj) &
            + cmn%lake(n)%calving / area_nbr    ! kg/m2/s
          !cmn%bmelt_grd_o(ii,jj) = cmn%bmelt_grd_o(ii,jj) &
          !  + cmn%bmelt_grd(i,j) / area_nbr    ! kg/m2/s
          !cmn%bmelt_flt_o(ii,jj) = cmn%bmelt_flt_o(ii,jj) &
          !  + cmn%bmelt_flt(i,j) / area_nbr    ! kg/m2/s
        enddo
      enddo
    endif
    
    cmn%runoff_o  = cmn%runoff_veg_o  + cmn%runoff_ice_o  + cmn%runoff_lake_o
    cmn%calving_o = cmn%calving_veg_o + cmn%calving_ice_o + cmn%calving_lake_o
    cmn%bmelt_o   = cmn%bmelt_grd_o   + cmn%bmelt_flt_o

    ! annual mean runoff 
    cmn%runoff_o_ann = cmn%runoff_o_ann + cmn%runoff_o/real(nday_year,wp)


   return

  end subroutine runoff_to_ocn


  subroutine aqua_init(cmn)

    implicit none

    type(cmn_class) :: cmn

    if (l_aqua_slab) then
      call nc_read(trim(restart_in_dir)//"/aqua_restart.nc","qflux",cmn%qflux) 
    endif

    return

  end subroutine aqua_init

  subroutine aqua_end(cmn)

    implicit none

    type(cmn_class) :: cmn
    character (len=256) :: fnm

    if (.not.l_aqua_slab) then
      fnm = trim(restart_in_dir)//"/aqua_restart.nc" 
      call nc_create(fnm)
      call nc_write_dim(fnm,"lon",x=lon,axis="x")
      call nc_write_dim(fnm,"lat",x=lat,axis="y")
      call nc_write(fnm,"qflux",cmn%flx_ocn,dims=["lon","lat"],long_name="",units="")
    endif

    return

  end subroutine aqua_end

  subroutine aquaplanet(cmn)

    implicit none

    type(cmn_class) :: cmn

    integer :: i, j
    real(wp) :: sst, fi, phi
    real(wp), parameter :: fim = pi/3._wp
    real(wp), parameter :: cslab = cap_w*1000.*30._wp ! 30 m water equivalent


    cmn%f_lnd = 0._wp 
    cmn%f_ocn = 1._wp
    cmn%z_sur_smooth_std = 0._wp

    cmn%idivide_pac_atl = 1
    cmn%idivide_atl_indpac = 1

    cmn%z_veg = 0._wp 
    cmn%z_ice = 0._wp 
    cmn%z_lake = 0._wp 
    cmn%z_sur = 0._wp 

    do i=1,ni
      do j=1,nj

        cmn%f_stp(i,j,:) = 0._wp
        cmn%f_stp(i,j,1) = 1._wp

        cmn%alb_vis_dir(i,j,1) = 0.07_wp 
        cmn%alb_vis_dif(i,j,1) = 0.07_wp 
        cmn%alb_nir_dir(i,j,1) = 0.07_wp 
        cmn%alb_nir_dif(i,j,1) = 0.07_wp 

        if (l_aqua_slab) then
          ! slab ocean
          cmn%t_skin(i,j,1) = cmn%t_skin(i,j,1) + (cmn%flx_ocn(i,j) - cmn%qflux(i,j))/cslab*dt_atm
        else
          ! prescribed zonally uniform SSTs
          phi = lat(j)*pi/180._wp
          fi=phi*pi/(2._wp*fim)
          if (abs(phi).lt.fim)then
            sst=27._wp*0.5_wp*(2._wp-sin(fi)**4-sin(fi)**2)
          else
            sst=0._wp
          endif
          cmn%t_skin(i,j,1)=sst+T0         
        endif

        ! compute surface fluxes
        cmn%sh(i,j,1) = 1.3e-3*cap_a*1.3_wp*cmn%wind(i,j,1)*(cmn%t_skin(i,j,1)-cmn%t2m(i,j,1))
        cmn%evp(i,j,1) = 1.3e-3*1.3_wp*cmn%wind(i,j,1)*(fqsat(cmn%t_skin(i,j,1),p0)-cmn%q2m(i,j,1))
        cmn%evp(i,j,1) = max(0._wp,cmn%evp(i,j,1))
        cmn%lh(i,j,1) = Le*cmn%evp(i,j,1) 
        cmn%lwu(i,j,1) = sigma*cmn%t_skin(i,j,1)**4

        ! 'ocean' heat flux, needed for slab ocean setup
        cmn%flx_ocn(i,j) = cmn%swnet(i,j,1) + cmn%lwd(i,j,1) - cmn%lwu(i,j,1) - cmn%sh(i,j,1) - cmn%lh(i,j,1) - cmn%snow(i,j,1)*Lf

      enddo
    enddo

    return

  end subroutine aquaplanet


  subroutine cmn_init(cmn)

    implicit none

    type(cmn_class) :: cmn

   
    ! Common grid definition
    call grid_init(cmn%grid,name="CMN-5x5",mtype="latlon",units="degrees",x=real(lon,dp),y=real(lat,dp))

    call cmn_alloc(cmn)

    ! get weights for interpolation from monthly to daily
    call monthly2daily(m0,m1,wtm0,wtm1)

    ! initalize all variables to zero
    cmn%mask_ocn = 0
    cmn%mask_lnd = 0
    cmn%mask_smb = 0
    cmn%mask_ice = 0
    cmn%f_ocn = 0._wp
    cmn%f_sic = 0._wp
    cmn%f_lnd0 = 0._wp
    cmn%f_lnd = 0._wp
    cmn%f_ice = 0._wp
    cmn%f_ice_grd = 0._wp
    cmn%f_ice_flt = 0._wp
    cmn%f_lake= 0._wp
    cmn%z_sur = 0._wp
    cmn%z_sur_n = 0._wp
    cmn%z_ocn_max = 0._wp
    cmn%z_veg = 0._wp
    cmn%z_veg_min = 0._wp
    cmn%z_veg_max = 0._wp
    cmn%z_ice = 0._wp
    cmn%z_lake = 0._wp
    cmn%z_sur_smooth_std = 0._wp
    cmn%z_veg_std = 0._wp
    cmn%f_stp = 1._wp/real(nsurf,wp)
    cmn%f_astp = 0._wp
    cmn%sst = 0._wp
    cmn%sss = 0._wp
    cmn%sss_dat = 0._wp
    cmn%sst_dat = 0._wp
    cmn%taux_dat = 0._wp
    cmn%tauy_dat = 0._wp
    cmn%flx_ocn = 0._wp
    cmn%qflux = 0._wp
    cmn%p_e_sic_ocn = 0._wp
    cmn%fw_brines = 0._wp
    cmn%taux = 0._wp
    cmn%tauy = 0._wp
    cmn%uo1 = 0._wp
    cmn%vo1 = 0._wp
    cmn%ssh = 0._wp
    cmn%t_skin= T0
    cmn%t2m = T0 
    cmn%q2m = 0._wp
    cmn%tam = T0 
    cmn%ram = 0._wp
    cmn%gams = 0._wp
    cmn%t2m_mon = 0._wp
    cmn%t2m_mon_lnd = 0._wp
    cmn%t2m_min_mon = T0+20._wp
    cmn%t_soil = T0 
    cmn%t_shelf = T0 
    cmn%alb_vis_dir = 0._wp
    cmn%alb_vis_dif = 0._wp
    cmn%alb_nir_dir = 0._wp
    cmn%alb_nir_dif = 0._wp
    cmn%alb_vis_dir_ice_semi = 0._wp
    cmn%alb_vis_dif_ice_semi = 0._wp
    cmn%alb_nir_dir_ice_semi = 0._wp
    cmn%alb_nir_dif_ice_semi = 0._wp
    cmn%prc  = 0._wp
    cmn%rain = 0._wp
    cmn%snow = 0._wp
    cmn%cld = 0._wp
    cmn%ps = 0._wp
    cmn%slp = 0._wp
    cmn%usur = 0._wp
    cmn%vsur = 0._wp
    cmn%u700 = 0._wp
    cmn%v700 = 0._wp
    cmn%wind = 0._wp
    cmn%z0m = 0._wp
    cmn%swd = 0._wp
    cmn%swd_vis_dir = 0._wp
    cmn%swd_nir_dir = 0._wp
    cmn%swd_vis_dif = 0._wp
    cmn%swd_nir_dif = 0._wp
    cmn%dswd_dalb_vis_dir = 0._wp
    cmn%dswd_dalb_nir_dir = 0._wp
    cmn%dswd_dalb_vis_dif = 0._wp
    cmn%dswd_dalb_nir_dif = 0._wp
    cmn%dswd_dz_nir_dir = 0._wp
    cmn%dswd_dz_nir_dif = 0._wp
    cmn%lwd = 0._wp
    cmn%swnet = 0._wp
    cmn%lwu = 0._wp
    cmn%sh = 0._wp
    cmn%lh = 0._wp
    cmn%evp = 0._wp
    cmn%runoff_veg = 0._wp
    cmn%calving_veg = 0._wp
    cmn%runoff_ice = 0._wp
    cmn%melt_ice = 0._wp
    cmn%calving_ice = 0._wp
    cmn%bmelt_grd = 0._wp
    cmn%bmelt_flt = 0._wp
    cmn%runoff_ice_l = 0._wp
    cmn%calving_ice_l = 0._wp
    cmn%melt_ice_i_mon = 0._wp
    cmn%acc_ice_i_mon = 0._wp
    cmn%runoff_ice_i_mon = 0._wp
    cmn%calving_ice_i = 0._wp
    cmn%bmelt_grd_i = 0._wp
    cmn%bmelt_flt_i = 0._wp
    cmn%runoff_o = 0._wp
    cmn%runoff_o_ann = 0._wp
    cmn%calving_o = 0._wp
    cmn%bmelt_o = 0._wp
    cmn%bmelt_grd_o = 0._wp
    cmn%bmelt_flt_o = 0._wp
    cmn%runoff_veg_o = 0._wp
    cmn%calving_veg_o = 0._wp
    cmn%runoff_ice_o = 0._wp
    cmn%melt_ice_o = 0._wp
    cmn%calving_ice_o = 0._wp
    cmn%runoff_lake_o = 0._wp
    cmn%calving_lake_o = 0._wp
    cmn%lake_p_e = 0._wp
    cmn%f_ice_lake = 0._wp
    cmn%i_runoff = 0
    cmn%j_runoff = 0
    cmn%i_runoff_ice = 0
    cmn%j_runoff_ice = 0
    cmn%i_runoff_veg = 0
    cmn%j_runoff_veg = 0
    cmn%cosz = 0._wp
    cmn%coszm = 0._wp
    cmn%solar = 0._wp
    cmn%solarm = 0._wp
    cmn%solarmin = 0._wp
    cmn%solarmax = 0._wp
    cmn%daylength = 0._wp
    cmn%dust_emis = 0._wp
    cmn%dust_dep = 0._wp
    cmn%delta_C_ocn_2d = 0._wp
    cmn%delta_C_lnd_2d = 0._wp
    cmn%doc_export = 0._wp
    cmn%doc13_export = 0._wp
    cmn%doc14_export = 0._wp
    cmn%poc_export = 0._wp
    cmn%poc13_export = 0._wp
    cmn%poc14_export = 0._wp
    cmn%weath_carb = 0._wp
    cmn%weath_sil = 0._wp
    cmn%weath13_carb = 0._wp
    cmn%weath13_sil = 0._wp
    cmn%weath14_carb = 0._wp
    cmn%weath14_sil = 0._wp
    cmn%f_carb = 1._wp
    cmn%alk_to_ocn_glob = 0._wp
    cmn%alk_from_lnd_glob = 0._wp
    

   return

  end subroutine cmn_init


  subroutine cmn_alloc(cmn)

    implicit none

    type(cmn_class) :: cmn


    allocate(cmn%delta_C_lnd_2d(ni,nj)) 
    allocate(cmn%delta_C13_lnd_2d(ni,nj)) 
    allocate(cmn%delta_C14_lnd_2d(ni,nj)) 
    allocate(cmn%delta_C_ocn_2d(ni,nj))  
    allocate(cmn%delta_C13_ocn_2d(ni,nj)) 
    allocate(cmn%delta_C14_ocn_2d(ni,nj)) 
    allocate(cmn%so4(ni,nj))
    allocate(cmn%o3_pl(11))
    allocate(cmn%o3(ni,nj,11))
    allocate(cmn%mask_ocn(ni,nj))
    allocate(cmn%mask_lnd(ni,nj))
    allocate(cmn%mask_smb(ni,nj))
    allocate(cmn%mask_ice(ni,nj))
    allocate(cmn%mask_coast(ni,nj))
    allocate(cmn%mask_coast2(ni,nj))
    allocate(cmn%f_ocn(ni,nj))
    allocate(cmn%f_ocn2(ni,nj))
    allocate(cmn%f_sic(ni,nj))
    allocate(cmn%f_lnd0(ni,nj))
    allocate(cmn%f_lnd(ni,nj))
    allocate(cmn%f_ice(ni,nj))    
    allocate(cmn%f_ice_grd(ni,nj))    
    allocate(cmn%f_ice_flt(ni,nj))    
    allocate(cmn%f_lake(ni,nj))    
    allocate(cmn%f_crop(ni,nj))    
    allocate(cmn%f_pasture(ni,nj))    
    allocate(cmn%disturbance(5,ni,nj))    
    allocate(cmn%z_sur(ni,nj))    
    allocate(cmn%z_sur_n(ni,nj,nsurf_macro))    
    allocate(cmn%z_ocn_max(ni,nj))    
    allocate(cmn%z_veg(ni,nj))    
    allocate(cmn%z_veg_min(ni,nj))    
    allocate(cmn%z_veg_max(ni,nj))    
    allocate(cmn%z_ice(ni,nj))    
    allocate(cmn%z_lake(ni,nj))    
    allocate(cmn%z_sur_smooth_std(ni,nj))    
    allocate(cmn%z_veg_std(ni,nj))    
    allocate(cmn%f_stp(ni,nj,nsurf))
    allocate(cmn%f_astp(ni,nj,nsurf_macro))
    allocate(cmn%coral_f_area(ni,nj,-250:50))
    allocate(cmn%coral_f_topo(ni,nj,-250:50))
    allocate(cmn%q_geo(ni,nj))
    allocate(cmn%sst(ni,nj))
    allocate(cmn%sss(ni,nj))
    allocate(cmn%sss_dat(ni,nj))
    allocate(cmn%sst_dat(ni,nj))
    allocate(cmn%taux_dat(ni,nj))
    allocate(cmn%tauy_dat(ni,nj))
    allocate(cmn%rbatm_dat(ni,nj))
    allocate(cmn%cld_dat(ni,nj))
    allocate(cmn%cld_day_dat(ni,nj))
    allocate(cmn%slp_dat(ni,nj))
    allocate(cmn%tsl_dat(ni,nj))
    allocate(cmn%htrop_dat(ni,nj))
    allocate(cmn%flx_ocn(ni,nj))
    allocate(cmn%qflux(ni,nj))
    allocate(cmn%p_e_sic_ocn(ni,nj))
    allocate(cmn%fw_brines(ni,nj))
    allocate(cmn%t_shelf(ni,nj))
    allocate(cmn%usur(ni,nj))
    allocate(cmn%vsur(ni,nj))
    allocate(cmn%u700(ni,nj))
    allocate(cmn%v700(ni,nj))
    allocate(cmn%wind(ni,nj,nsurf))
    allocate(cmn%taux(ni,nj,nsurf_macro))
    allocate(cmn%tauy(ni,nj,nsurf_macro))
    allocate(cmn%tauxo(ni,nj))
    allocate(cmn%tauyo(ni,nj))
    allocate(cmn%uo1(ni,nj))
    allocate(cmn%vo1(ni,nj))
    allocate(cmn%ssh(ni,nj))
    allocate(cmn%t_skin(ni,nj,nsurf))
    allocate(cmn%t2m(ni,nj,nsurf))
    allocate(cmn%q2m(ni,nj,nsurf))
    allocate(cmn%tam(ni,nj))
    allocate(cmn%ram(ni,nj))
    allocate(cmn%gams(ni,nj))
    allocate(cmn%t2m_min_mon(ni,nj))
    allocate(cmn%t2m_mon(ni,nj,nmon_year))
    allocate(cmn%t2m_mon_lnd(ni,nj,nmon_year))
    allocate(cmn%t_soil(ni,nj))
    allocate(cmn%alb_vis_dir(ni,nj,nsurf))
    allocate(cmn%alb_vis_dif(ni,nj,nsurf))
    allocate(cmn%alb_nir_dir(ni,nj,nsurf))
    allocate(cmn%alb_nir_dif(ni,nj,nsurf))
    allocate(cmn%alb_vis_dir_ice_semi(ni,nj,nday_year))
    allocate(cmn%alb_vis_dif_ice_semi(ni,nj,nday_year))
    allocate(cmn%alb_nir_dir_ice_semi(ni,nj,nday_year))
    allocate(cmn%alb_nir_dif_ice_semi(ni,nj,nday_year))
    allocate(cmn%prc(ni,nj))
    allocate(cmn%rain(ni,nj,nsurf))
    allocate(cmn%snow(ni,nj,nsurf))
    allocate(cmn%cld(ni,nj))
    allocate(cmn%ps(ni,nj,nsurf))
    allocate(cmn%slp(ni,nj))
    allocate(cmn%z0m(ni,nj,nsurf))
    allocate(cmn%swnet(ni,nj,nsurf))
    allocate(cmn%swd(ni,nj))
    allocate(cmn%swd_vis_dir(ni,nj))
    allocate(cmn%swd_nir_dir(ni,nj))
    allocate(cmn%swd_vis_dif(ni,nj))
    allocate(cmn%swd_nir_dif(ni,nj))
    allocate(cmn%dswd_dalb_vis_dir(ni,nj))
    allocate(cmn%dswd_dalb_nir_dir(ni,nj))
    allocate(cmn%dswd_dalb_vis_dif(ni,nj))
    allocate(cmn%dswd_dalb_nir_dif(ni,nj))
    allocate(cmn%dswd_dz_nir_dir(ni,nj))
    allocate(cmn%dswd_dz_nir_dif(ni,nj))
    allocate(cmn%lwd(ni,nj,nsurf))
    allocate(cmn%lwd_cs(ni,nj,nsurf_macro))
    allocate(cmn%lwd_cld(ni,nj,nsurf_macro))
    allocate(cmn%lwu(ni,nj,nsurf))
    allocate(cmn%sh(ni,nj,nsurf))
    allocate(cmn%lh(ni,nj,nsurf))
    allocate(cmn%evp(ni,nj,nsurf))
    allocate(cmn%runoff_veg(ni,nj)) 
    allocate(cmn%calving_veg(ni,nj))
    allocate(cmn%runoff_ice(ni,nj)) 
    allocate(cmn%melt_ice(ni,nj)) 
    allocate(cmn%calving_ice(ni,nj))
    allocate(cmn%bmelt_grd(ni,nj))
    allocate(cmn%bmelt_flt(ni,nj))
    allocate(cmn%runoff_ice_l(ni,nj))
    allocate(cmn%calving_ice_l(ni,nj))
    allocate(cmn%melt_ice_i_mon(ni,nj,nmon_year)) 
    allocate(cmn%acc_ice_i_mon(ni,nj,nmon_year)) 
    allocate(cmn%runoff_ice_i_mon(ni,nj,nmon_year)) 
    allocate(cmn%calving_ice_i(ni,nj))
    allocate(cmn%bmelt_grd_i(ni,nj))
    allocate(cmn%bmelt_flt_i(ni,nj))
    allocate(cmn%runoff_o(ni,nj)) 
    allocate(cmn%runoff_o_ann(ni,nj)) 
    allocate(cmn%calving_o(ni,nj)) 
    allocate(cmn%bmelt_o(ni,nj))
    allocate(cmn%bmelt_grd_o(ni,nj))
    allocate(cmn%bmelt_flt_o(ni,nj))
    allocate(cmn%runoff_veg_o(ni,nj)) 
    allocate(cmn%calving_veg_o(ni,nj))
    allocate(cmn%runoff_ice_o(ni,nj)) 
    allocate(cmn%melt_ice_o(ni,nj)) 
    allocate(cmn%calving_ice_o(ni,nj)) 
    allocate(cmn%runoff_lake_o(ni,nj)) 
    allocate(cmn%calving_lake_o(ni,nj)) 
    allocate(cmn%lake_p_e(ni,nj))
    allocate(cmn%h_lake(ni,nj))
    allocate(cmn%f_ice_lake(ni,nj))
    allocate(cmn%i_runoff(ni,nj))
    allocate(cmn%j_runoff(ni,nj))
    allocate(cmn%i_runoff_veg(ni,nj))
    allocate(cmn%j_runoff_veg(ni,nj))
    allocate(cmn%i_runoff_ice(ni,nj))
    allocate(cmn%j_runoff_ice(ni,nj))
    allocate(cmn%idivide_pac_atl(nj))
    allocate(cmn%idivide_atl_indpac(nj))
    allocate(cmn%cosz(nday_year,24,nj))
    allocate(cmn%coszm(nday_year,nj))
    allocate(cmn%solar(nday_year,24,nj))
    allocate(cmn%solarm(nday_year,nj))
    allocate(cmn%solarmin(nday_year,nj))
    allocate(cmn%solarmax(nday_year,nj))
    allocate(cmn%daylength(nday_year,nj))
    allocate(cmn%dust_emis(ni,nj))
    allocate(cmn%dust_dep(ni,nj))
    allocate(cmn%doc_export(ni,nj))
    allocate(cmn%doc13_export(ni,nj))
    allocate(cmn%doc14_export(ni,nj))
    allocate(cmn%poc_export(ni,nj))
    allocate(cmn%poc13_export(ni,nj))
    allocate(cmn%poc14_export(ni,nj))
    allocate(cmn%weath_carb(ni,nj))
    allocate(cmn%weath_sil(ni,nj))
    allocate(cmn%weath13_carb(ni,nj))
    allocate(cmn%weath13_sil(ni,nj))
    allocate(cmn%weath14_carb(ni,nj))
    allocate(cmn%weath14_sil(ni,nj))
    allocate(cmn%f_carb(ni,nj))

   return

  end subroutine cmn_alloc

end module coupler
 
