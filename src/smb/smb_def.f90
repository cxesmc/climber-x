!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ d e f
!
!  Purpose : definition of surface mass balance model class 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
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
module smb_def

    use precision, only : wp
    use coord, only : grid_class
    use coord, only : map_class, map_scrip_class
    use timer, only : nmon_year, nday_year

    implicit none

    type smb_in_class
      ! grid definitions
      type(grid_class) :: grid    !! coupler grid object
      real(wp) :: co2         ! atmospheric CO2 concentration [ppm]
      real(wp) :: Smax65N     ! maximum summer insolation at 65N [W/m2]
      ! input variables on coupler grid
      real(wp), dimension(:,:), allocatable :: z_sur  !! grid cell mean elevation [m]
      real(wp), dimension(:,:), allocatable :: t2m    !! surface air temperature at grid cell mean elevation [K]
      real(wp), dimension(:,:,:), allocatable :: t2m_bias    !! surface air temperature bias at present-day [K]
      real(wp), dimension(:,:,:), allocatable :: prc_bias    !! precipitation bias at present-day (ratio model/obs) [1]
      real(wp), dimension(:,:), allocatable :: tam    !! atmospheric temperature [K] 
      real(wp), dimension(:,:), allocatable :: ram    !! atmospheric relative humidity [1] 
      real(wp), dimension(:,:), allocatable :: gam    !! atmospheric lapse rate at the surface [K/m] 
      real(wp), dimension(:,:), allocatable :: tstd   !! standard deviation on daily 2m temperature [K] 
      real(wp), dimension(:,:), allocatable :: prc    !! total precipitation rate [kg/m2/s] 
      real(wp), dimension(:,:), allocatable :: u700   !! zonal wind component at 700 hPa [m/s]
      real(wp), dimension(:,:), allocatable :: v700   !! meridional wind component at 700 hPa [m/s]
      real(wp), dimension(:,:), allocatable :: wind   !! surface wind speed [m/s]
      real(wp), dimension(:,:), allocatable :: cld    !! cloud cover fraction [/]
      real(wp), dimension(:,:), allocatable :: dust   !! dust deposition rate [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: swd_toa   !! downward shortwave radiation at TOA [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_toa_min   !! minimum diurnal downward shortwave radiation at TOA [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_vis_dir   !! downward shortwave visible radiation at surface, clear sky [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_nir_dir   !! downward shortwave near-infrared radiation at surface, clear sky [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_vis_dif   !! downward shortwave visible radiation at surface, cloudy [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_nir_dif   !! downward shortwave near-infrared radiation at surface, cloudy [W/m2]
      real(wp), dimension(:,:), allocatable :: dswd_dalb_vis_dir
      real(wp), dimension(:,:), allocatable :: dswd_dalb_nir_dir
      real(wp), dimension(:,:), allocatable :: dswd_dalb_vis_dif
      real(wp), dimension(:,:), allocatable :: dswd_dalb_nir_dif
      real(wp), dimension(:,:), allocatable :: dswd_dz_nir_dir
      real(wp), dimension(:,:), allocatable :: dswd_dz_nir_dif
      real(wp), dimension(:,:), allocatable :: alb_vis_dir   !! surface albedo for visible radiation at surface, clear sky []
      real(wp), dimension(:,:), allocatable :: alb_nir_dir   !! surface albedo for near-infrared radiation at surface, clear sky []
      real(wp), dimension(:,:), allocatable :: alb_vis_dif   !! surface albedo for visible radiation at surface, cloudy []
      real(wp), dimension(:,:), allocatable :: alb_nir_dif   !! surface albedo for near-infrared radiation at surface, cloudy []
      real(wp), dimension(:,:), allocatable :: lwdown   !! downward surface longwave radiation at grid-cell mean elevation [W/m2]
      real(wp), dimension(:,:), allocatable :: gam_lw   !! downward surface longwave radiation 'lapse rate' [W/m2]
      real(wp), dimension(:,:), allocatable :: t_ground  !! annual mean ground temperature (soil over land and bottom water over ocean) [degC]
      real(wp), dimension(:,:), allocatable :: coszm  !! daily mean cosine of solar zenith angle []
    end type

    type ts_out
      real(wp) :: Aice
      real(wp) :: smb
      real(wp) :: prc
      real(wp) :: snow
      real(wp) :: refreezing
      real(wp) :: evp
      real(wp) :: melt
      real(wp) :: run
      real(wp) :: smb_avg
      real(wp) :: prc_avg
      real(wp) :: snow_avg
      real(wp) :: refreezing_avg
      real(wp) :: evp_avg
      real(wp) :: melt_avg
      real(wp) :: run_avg
    end type

    type s_out
      ! interpolated variables
      real(wp), dimension(:,:), allocatable :: z_sur_i  !! mean surface elevation for atmosphere interpolated to ice sheet grid [m]
      real(wp), dimension(:,:), allocatable :: t2m_i    !! surface air temperature at mean surface elevation interpolated to ice sheet grid [K]
      real(wp), dimension(:,:), allocatable :: t2m_bias_i    !! surface air temperature bias interpolated to ice sheet grid [K]
      real(wp), dimension(:,:), allocatable :: prc_bias_i    !! precipitation bias at present-day (ratio model/obs) interpolated to ice sheet grid [1]
      real(wp), dimension(:,:), allocatable :: tam_i     !! atmospheric temperature interpolated to ice sheet grid [K] 
      real(wp), dimension(:,:), allocatable :: ram_i     !! atmospheric relative humidity interpolated to ice sheet grid [1] 
      real(wp), dimension(:,:), allocatable :: gam_i     !! atmospheric temperature lapse rate interpolated to ice sheet grid [K/m] 
      real(wp), dimension(:,:), allocatable :: tstd_i    !! standard deviation of daily surface air/skin temperature interpolated to ice sheet grid [K]
      real(wp), dimension(:,:), allocatable :: prc_i    !! total precipitation rate interpolated to ice sheet grid [kg/m2/s] 
      real(wp), dimension(:,:), allocatable :: u700_i      !! zonal wind component at 700 hPa interpolated to ice sheet grid [m/s]
      real(wp), dimension(:,:), allocatable :: v700_i      !! meridional wind component at 700 hPa interpolated to ice sheet grid [m/s]
      real(wp), dimension(:,:), allocatable :: wind_i   !! surface wind speed interpolated to ice sheet grid [m/s]
      real(wp), dimension(:,:), allocatable :: cld_i    !! cloud cover fraction interpolated to ice sheet grid [/]
      real(wp), dimension(:,:), allocatable :: dust_i   !! dust deposition rate interpolated to ice sheet grid [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: swdown_i   !! downward shortwave radiation at the surface interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_toa_i   !! downward shortwave radiation at top of atmosphere interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_toa_min_i   !! minimum diurnal downward shortwave radiation at top of atmosphere interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_vis_dir_i   !! downward surface shortwave visible radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_nir_dir_i   !! downward surface shortwave near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_vis_dif_i   !! downward surface shortwave visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_nir_dif_i   !! downward surface shortwave near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: dswd_dalb_vis_dir_i  
      real(wp), dimension(:,:), allocatable :: dswd_dalb_nir_dir_i
      real(wp), dimension(:,:), allocatable :: dswd_dalb_vis_dif_i
      real(wp), dimension(:,:), allocatable :: dswd_dalb_nir_dif_i
      real(wp), dimension(:,:), allocatable :: dswd_dz_nir_dir_i
      real(wp), dimension(:,:), allocatable :: dswd_dz_nir_dif_i
      real(wp), dimension(:,:), allocatable :: alb_vis_dir_i   !! surface albedo for visible radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_nir_dir_i   !! surface albedo for near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_vis_dif_i   !! surface albedo for visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_nir_dif_i   !! surface albedo for near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: coszm_i  !! daily mean cosine of solar zenith angle interpolated to ice sheet grid []
      real(wp), dimension(:,:), allocatable :: lwdown_i    !! downward surface longwave radiation interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: gam_lw_i    !! 'lapse rate' of downward surface longwave radiation interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: t_ground_i  !! annual mean ground temperature (soil over land and bottom water over ocean) interpolated to ice sheet grid [degC]
      ! smb variables
      real(wp), dimension(:,:), allocatable :: area    !! grid cell area [km2]
      integer, dimension(:,:), allocatable :: mask_smb   !! mask where smb is computed []
      integer, dimension(:,:), allocatable :: mask_ice     !! ice mask on ice sheet grid [/]
      integer, dimension(:,:), allocatable :: mask_maxice  !! mask of maximum ice extent []
      integer, dimension(:,:), allocatable :: mask_margin    !! ice margin mask on ice sheet grid [/]
      real(wp), dimension(:,:), allocatable :: f_ice !! ice cover fraction [1]
      real(wp), dimension(:,:), allocatable :: h_ice !! ice thickness [m]
      real(wp), dimension(:,:), allocatable :: z_sur     !! surface elevation on ice sheet grid [m]
      real(wp), dimension(:,:), allocatable :: z_sur_eff     !! effective surface elevation on ice sheet grid [m]
      real(wp), dimension(:,:), allocatable :: z_sur_fil     !! filtered surface elevation on ice sheet grid [m]
      real(wp), dimension(:,:), allocatable :: z_sur_std     !! sub-grid standard deviation of surface elevation [m]
      real(wp), dimension(:,:), allocatable :: dz_dx_sur  !! surface x-slope on ice sheet grid [m/m]
      real(wp), dimension(:,:), allocatable :: dz_dy_sur  !! surface y-slope on ice sheet grid [m/m]
      real(wp), dimension(:,:), allocatable :: dz_sur     !! surface slope on ice sheet grid [m/m]
      real(wp), dimension(:,:), allocatable :: tam      !! atmospheric temperature at ice sheet elevation [K]
      real(wp), dimension(:,:), allocatable :: t2m      !! 2m temperature at ice sheet elevation [K]
      real(wp), dimension(:,:), allocatable :: t_skin   !! skin temperature [K]
      real(wp), dimension(:,:), allocatable :: tskin_tam   !! skin - atm temperature [K]
      real(wp), dimension(:,:), allocatable :: t_skin_amp   !! amplitude of diurnal cycle of skin temperature [K]
      real(wp), dimension(:,:,:), allocatable :: t_prof  !! snow/ice temperature [K]
      real(wp), dimension(:,:), allocatable :: q2m    !! surface air specific humidity [kg/kg]
      real(wp), dimension(:,:), allocatable :: pressure    !! surface pressure [Pa]
      real(wp), dimension(:,:), allocatable :: mask_snow    !! snow mask [/]
      real(wp), dimension(:,:), allocatable :: f_snow !! snow cover fraction [1]
      real(wp), dimension(:,:), allocatable :: h_snow !! snow thickness [m]
      real(wp), dimension(:,:), allocatable :: w_snow !! snow water equivalent [kg/m2]
      real(wp), dimension(:,:), allocatable :: w_snow_max !! max snow water equivalent [kg/m2]
      real(wp), dimension(:,:), allocatable :: snowmelt !! snow melt [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: icemelt !! ice melt [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: melt !! ice + snow melt [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: refreezing !! refreezing [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: f_rfz_to_snow !! fraction of refreezing going to snow [1]
      real(wp), dimension(:,:), allocatable :: runoff !! runoff [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: albedo  !! allwave surface albedo [/]
      real(wp), dimension(:,:), allocatable :: alb_bg  !! allwave background albedo [/]
      real(wp), dimension(:,:), allocatable :: alb_ice  !! allwave ice albedo [/]
      real(wp), dimension(:,:), allocatable :: alb_vis_dir   !! surface albedo for visible radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_nir_dir   !! surface albedo for near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_vis_dif   !! surface albedo for visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_nir_dif   !! surface albedo for near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_snow_vis_dir   !! snow_albedo for visible radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_snow_nir_dir   !! snow_albedo for near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_snow_vis_dif   !! snow_albedo for visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_snow_nir_dif   !! snow_albedo for near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: snow_grain  !! snow grain size [um]
      real(wp), dimension(:,:), allocatable :: dust_con  !! dust concentration in snow [ppmw]
      real(wp), dimension(:,:), allocatable :: cld  !! cloud cover fraction [/]
      real(wp), dimension(:,:), allocatable :: swnet  !! net surface shortwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: swnet_min  !! minimum diurnal net surface shortwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: swdown  !! downward surface shortwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: lwdown !! downward surface longwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: r_a  !! aerodynamic resistance [s/m]
      real(wp), dimension(:,:), allocatable :: flx_g  !! ground heat flux [W/m2]
      real(wp), dimension(:,:), allocatable :: dflxg_dT  !! derivative of ground heat flux wrt top layer temperature [W/m2/K]
      real(wp), dimension(:,:), allocatable :: flx_melt !! heat flux used to melt snow [W/m2]
      real(wp), dimension(:,:), allocatable :: flx_sh  !! sensible heat flux [W/m2]
      real(wp), dimension(:,:), allocatable :: flx_lwu  !! upwelling longwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: flx_lh  !! latent heat flux [W/m2]
      real(wp), dimension(:,:), allocatable :: evp   !! evaporation rate [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: prc   !! precipitation rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: f_ele  !! precipitation rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: f_wind !! precipitation rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: rain   !! rainfall rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: snow   !! snowfall rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: u700   !! zonal wind component at 700 hPa, downscaled [m/s]
      real(wp), dimension(:,:), allocatable :: v700   !! meridional wind component at 700 hPa, downscaled [m/s]
      real(wp), dimension(:,:), allocatable :: wind   !! surface wind speed, downscaled [m/s]
      real(wp), dimension(:,:), allocatable :: pdd    !! positive degree days [degC]
      ! annually cumulated variables
      real(wp), dimension(:,:), allocatable :: ann_smb        !! annual surface mass balance [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_smb_ice    !! annual surface mass balance, for ice cells only [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_smb_noice  !! annual surface mass balance, for ice-free cells only [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_smb_pdd    !! annual surface mass using PDD scheme [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_prc        !! annual precipitation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_prc_ice        !! annual precipitation, for ice cells only [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_snow       !! annual snowfall [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_snow_ice       !! annual snowfall, for ice cells only [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_ablation   !! annual surface ablation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_ablation_ice   !! annual surface ablation, for ice cells only [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_melt       !! annual surface melt [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_melt_ice       !! annual surface melt, for ice cells only [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_evp        !! annual evaporation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_evp_ice        !! annual evaporation, for ice cells only [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_runoff     !! annual runoff [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_runoff_ice     !! annual runoff, for ice cells only [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_refreezing !! annual refreezing [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_refreezing_ice !! annual refreezing, for ice cells only [kg/m2]
      ! output variables
      real(wp), dimension(:,:), allocatable :: t_ice      !! annual mean top ice temperature on ice sheet grid [degC]
      real(wp), dimension(:,:), allocatable :: t_ground   !! ground temperature on ice sheet grid [degC]
    end type

    type simple_smb_class
      real(wp), dimension(:,:), allocatable :: pdd
      real(wp), dimension(:,:), allocatable :: snow_cum
      real(wp), dimension(:,:), allocatable :: rain_cum
      real(wp), dimension(:,:), allocatable :: t2m_cum
      real(wp), dimension(:,:), allocatable :: melt
      real(wp), dimension(:,:), allocatable :: melt_star
      real(wp), dimension(:,:), allocatable :: runoff
      real(wp), dimension(:,:), allocatable :: smb
      real(wp), dimension(:,:), allocatable :: t_ice
    end type

    type grid_smb_to_cmn_type
      integer, dimension(:,:), allocatable :: i_lowres 
      integer, dimension(:,:), allocatable :: j_lowres
      integer, dimension(:,:), allocatable :: ncells
      integer, dimension(:,:), allocatable :: ncells_ice
    end type

    type smb_class

      ! grid definition
      type(grid_class) :: grid
      type(grid_class) :: grid_latlon
      type(grid_smb_to_cmn_type) :: grid_smb_to_cmn
      type(map_class) :: map_cmn_to_ice
      type(map_class) :: map_to_latlon
      type(map_class) :: map_from_latlon
      type(map_scrip_class) :: maps_cmn_to_ice
      type(map_scrip_class) :: maps_to_latlon
      type(map_scrip_class) :: maps_from_latlon

      type(simple_smb_class) :: simple

      integer :: ncells                                        !! total number of active smb columns
      integer, dimension(:), allocatable :: idx_cell_active    !! index of active smb cells
      integer, dimension(:, :), allocatable :: ij_1d           !! n --> i,j  (2 x ncol)
      integer, dimension(:, :), allocatable :: id_map          !! i,j --> n  (ni, nj)

      real(wp) :: dTvar       ! imposed temperature variability [K]
      real(wp) :: dTvar_ann       ! imposed temperature variability [K]
      real(wp) :: dTvar_day       ! imposed temperature variability [K]

      ! input variables from ice sheet (matches smb grid)
      integer, dimension(:,:), allocatable :: mask_ice   !! ice mask on ice sheet grid []
      integer, dimension(:,:), allocatable :: mask_ice_old   !! old ice mask on ice sheet grid []
      integer, dimension(:,:), allocatable :: mask_smb   !! mask where smb is computed []
      integer, dimension(:,:), allocatable :: mask_smb_tmp   !! mask where smb is computed []
      integer, dimension(:,:), allocatable :: mask_maxice  !! mask of maximum ice extent []
      integer, dimension(:,:), allocatable :: mask_margin    !! ice margin mask on ice sheet grid [/]
      real(wp), dimension(:,:), allocatable :: z_sur      !! surface elevation on ice sheet grid [m]
      real(wp), dimension(:,:), allocatable :: z_sur_eff     !! effective surface elevation on ice sheet grid [m]
      real(wp), dimension(:,:), allocatable :: z_sur_fil     !! filtered surface elevation on ice sheet grid [m]
      real(wp), dimension(:,:), allocatable :: z_sur_std     !! sub-grid standard deviation of surface elevation [m]
      real(wp), dimension(:,:), allocatable :: z_bed_std     !! sub-grid standard deviation of bedrock elevation [m]
      real(wp), dimension(:,:), allocatable :: h_ice      !! ice thickness on ice sheet grid [m]

      ! interpolated fields
      real(wp), dimension(:,:), allocatable :: z_sur_i  !! mean surface elevation for atmosphere interpolated to ice sheet grid [m]
      real(wp), dimension(:,:), allocatable :: t2m_i    !! surface air temperature at mean surface elevation interpolated to ice sheet grid [K]
      real(wp), dimension(:,:,:), allocatable :: t2m_bias_i    !! surface air temperature bias interpolated to ice sheet grid [K]
      real(wp), dimension(:,:,:), allocatable :: prc_bias_i    !! precipitation bias at present-day (ratio model/obs) interpolated to ice sheet grid [1]
      real(wp), dimension(:,:), allocatable :: tam_i    !! atmospheric temperature interpolated to ice sheet grid [K] 
      real(wp), dimension(:,:), allocatable :: ram_i    !! atmospheric relative humidity interpolated to ice sheet grid [1] 
      real(wp), dimension(:,:), allocatable :: gam_i     !! atmospheric temperature lapse rate interpolated to ice sheet grid [K/m] 
      real(wp), dimension(:,:), allocatable :: tstd_i   !! standard deviation of daily surface air/skin temperature interpolated to ice sheet grid [K]
      real(wp), dimension(:,:), allocatable :: prc_i    !! total precipitation rate interpolated to ice sheet grid [kg/m2/s] 
      real(wp), dimension(:,:), allocatable :: u700_i      !! zonal wind component at 700 hPa interpolated to ice sheet grid [m/s]
      real(wp), dimension(:,:), allocatable :: v700_i      !! meridional wind component at 700 hPa interpolated to ice sheet grid [m/s]
      real(wp), dimension(:,:), allocatable :: wind_i   !! surface wind speed interpolated to ice sheet grid [m/s]
      real(wp), dimension(:,:), allocatable :: cld_i    !! cloud cover fraction interpolated to ice sheet grid [/]
      real(wp), dimension(:,:), allocatable :: dust_i   !! dust deposition rate interpolated to ice sheet grid [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: swd_toa_i   !! downward shortwave radiation at top of atmosphere interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_toa_min_i   !! minimum diurnal downward shortwave radiation at top of atmosphere interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_vis_dir_i   !! downward surface shortwave visible radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_nir_dir_i   !! downward surface shortwave near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_vis_dif_i   !! downward surface shortwave visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: swd_sur_nir_dif_i   !! downward surface shortwave near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: dswd_dalb_vis_dir_i  
      real(wp), dimension(:,:), allocatable :: dswd_dalb_nir_dir_i
      real(wp), dimension(:,:), allocatable :: dswd_dalb_vis_dif_i
      real(wp), dimension(:,:), allocatable :: dswd_dalb_nir_dif_i
      real(wp), dimension(:,:), allocatable :: dswd_dz_nir_dir_i
      real(wp), dimension(:,:), allocatable :: dswd_dz_nir_dif_i
      real(wp), dimension(:,:), allocatable :: alb_vis_dir_i   !! surface albedo for visible radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_nir_dir_i   !! surface albedo for near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_vis_dif_i   !! surface albedo for visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_nir_dif_i   !! surface albedo for near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: coszm_i  !! daily mean cosine of solar zenith angle interpolated to ice sheet grid []
      real(wp), dimension(:,:), allocatable :: lwdown_i    !! downward surface longwave radiation interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: gam_lw_i    !! 'lapse rate' of downward surface longwave radiation interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: t_ground_i  !! annual mean ground temperature (soil over land and bottom water over ocean) interpolated to ice sheet grid [degC]

      ! internal variables
      real(wp), dimension(:,:), allocatable :: dz_dx_sur  !! ice surface x-slope on ice sheet grid [/]
      real(wp), dimension(:,:), allocatable :: dz_dy_sur  !! ice surface y-slope on ice sheet grid [/]
      real(wp), dimension(:,:), allocatable :: dz_sur     !! ice surface slope on ice sheet grid [/]
      real(wp), dimension(:,:), allocatable :: tam      !! atmospheric temperature at ice sheet elevation [K]
      real(wp), dimension(:,:), allocatable :: t2m      !! 2m temperature at ice sheet elevation [K]
      real(wp), dimension(:,:), allocatable :: t_skin   !! skin temperature [K]
      real(wp), dimension(:,:), allocatable :: t_skin_amp   !! amplitude of diurnal cycle of skin temperature [K]
      real(wp), dimension(:,:), allocatable :: t_skin_old   !! old skin temperature [K]
      real(wp), dimension(:,:,:), allocatable :: t_prof  !! snow/ice temperature [K]
      real(wp), dimension(:,:,:), allocatable :: t_prof_old  !! old snow/ice temperature [K]
      real(wp), dimension(:,:), allocatable :: q2m    !! surface air specific humidity [kg/kg]
      real(wp), dimension(:,:), allocatable :: pressure    !! surface pressure [Pa]
      integer, dimension(:,:), allocatable :: mask_snow    !! snow mask [/]
      real(wp), dimension(:,:), allocatable :: f_snow !! snow cover fraction [1]
      real(wp), dimension(:,:), allocatable :: h_snow !! snow thickness [m]
      real(wp), dimension(:,:), allocatable :: w_snow !! snow water equivalent [kg/m2]
      real(wp), dimension(:,:), allocatable :: w_snow_old !! old snow water equivalent [kg/m2]
      real(wp), dimension(:,:), allocatable :: w_snow_max !! max snow water equivalent [kg/m2]
      real(wp), dimension(:,:), allocatable :: snowmelt !! snow melt [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: icemelt !! ice melt [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: melt !! ice + snow melt [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: refreezing !! refreezing [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: refreezing_sum !! refreezing summed up over the year [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: f_rfz_to_snow !! fraction of refreezing going to snow [1]
      real(wp), dimension(:,:), allocatable :: runoff  !! runoff [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: cod  !! cloud optical depth [/]
      real(wp), dimension(:,:), allocatable :: albedo  !! allwave surface albedo [/]
      real(wp), dimension(:,:), allocatable :: alb_bg  !! allwave background albedo [/]
      real(wp), dimension(:,:), allocatable :: f_ice !! ice cover fraction [1]
      real(wp), dimension(:,:), allocatable :: f_ice_old !! old ice cover fraction [1]
      real(wp), dimension(:,:), allocatable :: dt_snowfree !! cumulated time over the year when snowfree [s]
      real(wp), dimension(:,:), allocatable :: alb_ice  !! allwave ice albedo [/]
      real(wp), dimension(:,:), allocatable :: alb_vis_dir   !! surface albedo for visible radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_nir_dir   !! surface albedo for near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_vis_dif   !! surface albedo for visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_nir_dif   !! surface albedo for near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_snow_vis_dir   !! snow_albedo for visible radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_snow_nir_dir   !! snow_albedo for near-infrared radiation for clear sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_snow_vis_dif   !! snow_albedo for visible radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: alb_snow_nir_dif   !! snow_albedo for near-infrared radiation for cloudy sky interpolated to ice sheet grid [W/m2]
      real(wp), dimension(:,:), allocatable :: snow_grain  !! snow grain size [um]
      real(wp), dimension(:,:), allocatable :: dust_con  !! dust concentration in snow [kg/kg]
      real(wp), dimension(:,:), allocatable :: cld  !! cloud cover fraction [/]
      real(wp), dimension(:,:), allocatable :: swdown  !! downward surface shortwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: swnet  !! net surface shortwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: swnet_min  !! minimum diurnal net surface shortwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: lwdown !! downward surface longwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: r_a  !! aerodynamic resistance [s/m]
      real(wp), dimension(:,:), allocatable :: flx_g  !! ground heat flux [W/m2]
      real(wp), dimension(:,:), allocatable :: dflxg_dT  !! derivative of ground heat flux wrt top layer temperature [W/m2/K]
      real(wp), dimension(:,:), allocatable :: flx_melt !! heat flux used to melt snow [W/m2]
      real(wp), dimension(:,:), allocatable :: flx_sh  !! sensible heat flux [W/m2]
      real(wp), dimension(:,:), allocatable :: flx_lwu  !! upwelling longwave radiation [W/m2]
      real(wp), dimension(:,:), allocatable :: flx_lh  !! latent heat flux [W/m2]
      real(wp), dimension(:,:), allocatable :: evp   !! evaporation rate [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: prc   !! precipitation rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: f_ele  !! precipitation rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: f_wind  !! precipitation rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: rain   !! rainfall rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: snow   !! snowfall rate, downscaled [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: u700   !! zonal wind component at 700 hPa, downscaled [m/s]
      real(wp), dimension(:,:), allocatable :: v700   !! meridional wind component at 700 hPa, downscaled [m/s]
      real(wp), dimension(:,:), allocatable :: wind   !! surface wind speed, downscaled [m/s]

      real(wp), dimension(:,:), allocatable :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
      real(wp), dimension(:,:), allocatable :: f_sh, f_e, f_lh, f_lw, qsat, dqsatdT

      ! annually cumulated variables
      real(wp), dimension(:,:), allocatable :: ann_smb        !! annual surface mass balance [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_smb_ref    !! annual reference surface mass balance [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_smb_cx_ref !! annual CLIMBER-X reference surface mass balance [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_smb_avg_ref
      real(wp), dimension(:,:), allocatable :: ann_prc        !! annual precipitation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_prc_ref    !! annual reference precipitation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_prc_cx_ref !! annual CLIMBER-X reference precipitation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_prc_avg_ref
      real(wp), dimension(:,:), allocatable :: ann_snow       !! annual snowfall [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_ablation   !! annual surface ablation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_melt       !! annual surface melt [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_icemelt    !! annual ice melt [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_evp        !! annual evaporation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_evp_ref    !! annual reference evaporation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_evp_cx_ref !! annual CLIMBER-X reference evaporation [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_evp_avg_ref
      real(wp), dimension(:,:), allocatable :: ann_runoff     !! annual runoff [kg/m2]
      real(wp), dimension(:,:), allocatable :: ann_refreezing !! annual refreezing [kg/m2]

      ! output variables
      real(wp), dimension(:,:,:), allocatable :: mon_runoff     !! monthly runoff on ice sheet grid [kg/m2/s]
      real(wp), dimension(:,:,:), allocatable :: mon_icemelt    !! monthly icemelt on ice sheet grid [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: t_ice       !! annual mean top ice temperature on ice sheet grid [degC]
      real(wp), dimension(:,:), allocatable :: t_ground    !! ground temperature (soil over land and bottom water over ocean) on ice sheet grid [degC]

      ! output types
      integer :: nout
      type(ts_out), allocatable :: ann_ts(:)
      type(ts_out) :: mon_ts(nmon_year)
      type(s_out) :: day_s(nday_year), mon_s(nmon_year), ann_s

    end type

    private
    public :: smb_class, smb_in_class, ts_out, s_out

end module
