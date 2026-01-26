!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ p a r a m s
!
!  Purpose : SMB model parameters
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
module smb_params

  use precision, only : wp, dp
  use ncio
  use nml
  use timer, only : dt_smb, sec_day, sec_year, n_year_smb
  use control, only : out_dir
  use constants, only : T0, rho_w, rho_i

  implicit none

  real(wp) :: dt, rdt

  integer :: i_smb

  logical :: l_maxice                ! apply large negative SMB (-4 m/yr) outside ice area specified by mask_maxice_file?
  character (len=256) :: mask_maxice_file ! maximum ice extent mask file
  real(wp) :: smb_maxice    ! smb to be applied outside of allowed ice mask

  real(wp) :: smb_crit_mask
  integer :: n_smb_mask_ext

  logical :: l_regional_climate_forcing

  logical :: l_smb_bias_corr 
  character (len=256) :: smb_ref_file 
  character (len=256) :: smb_cx_ref_file 
  integer :: year_ini_smb_ref 
  integer :: year_end_smb_ref 

  integer :: i_t2m_bias_corr
  integer :: i_prc_bias_corr
  character (len=256) :: bias_corr_file
  real(wp) :: t2m_bias_scale_fac
  real(wp) :: t2m_bias_corr_uniform
  real(wp) :: dLW_dT_fac

  real(wp) :: d_cld

  real(wp) :: wind_ele_fac

  logical :: l_Tvar_ann
  real(wp) :: Tvar_ann_amp 
  integer :: Tvar_ann_period
  logical :: l_Tvar_day
  real(wp) :: Tvar_day_amp 
  integer :: Tvar_day_period

  character(len=256) :: map_method
  real(dp) :: filt_sigma

  real(wp), parameter :: p0 = 1010.e2_wp    !! reference sea level pressure [Pa], fixme, should be passed by coupler!!

  real(wp), parameter :: h_atm = 8600._wp !! atmosphere scale height [m], fixme, not consistent with atmosphere!

  logical :: l_diurnal_cycle
  real(wp) :: tstd_scale

  logical :: l_neutral
  logical :: l_dew

  integer :: nday_update_climate

  logical :: l_write_timer

  logical :: l_monthly_output
  logical :: l_daily_output

  integer :: i_gamma=1
  real(wp) :: gamma

  integer :: i_z_sur_eff
  real(wp) :: alpha_zstd

  ! PDD parameters
  type pdd_par_type
    integer :: i_ablation
    real(wp) :: s_stat  !! Standard deviation of the air temperature, in deg C
    real(wp) :: beta1   !! Degree-day factor for snow, in (mm WE)/(d*deg C)
    real(wp) :: beta2   !! Degree-day factor for ice, in (mm WE)/(d*deg C)
    real(wp) :: Pmax    !! Saturation factor for the formation of superimposed ice
    real(wp) :: mu      !! Firn-warming correction, in (d*deg C)/(mm WE)
  end type pdd_par_type
  type(pdd_par_type) :: pdd_par

  ! snow parameters
  type snow_par_type
    logical :: lsnow_dust 
    logical :: lsnow_aging 
    integer :: isnow_albedo 
    real(wp) :: dalb_snow_vis
    real(wp) :: dalb_snow_nir
    real(wp) :: w_snow_crit = 5._wp    ! kg/m2, minimum snow water equivalent for explicit snow layer
    real(wp) :: snow_grain_fresh ! um, fresh snow grain size
    real(wp) :: snow_grain_old = 1000._wp ! um, old snow grain size
    real(wp) :: f_age_t      ! temperature parameter for snow aging
    real(wp) :: dT_age
    real(wp) :: snow_0      ! kg/m2/s, critical snowfall rate for aging
    real(wp) :: snow_1      ! kg/m2/s, minimum snowfall rate for aging
    real(wp) :: k_sigma_orog
    real(wp) :: sigma_orog_crit
    real(wp) :: c_fsnow
    real(wp) :: c_fsnow_orog
    logical :: l_fsnow_orog
    real(wp) :: rho
    real(wp) :: lambda 
    real(wp) :: w_snow_max  = 1000._wp ! kg/m2, maximum snow water equivalent
    real(wp) :: w_snow_dust
    real(wp) :: dust_con_scale
    integer :: i_rfz
    real(wp) :: porosity 
    real(wp) :: f_rfz_max   ! parameter to control maxium fraction of refreezing 
    real(wp) :: f_rfz_to_snow_max   ! parameter to control fraction of refreezing going into snow
    real(wp) :: wsnow_crit_rfz   ! parameter to control fraction of refreezing going into snow
  end type
  type(snow_par_type) :: snow_par

  type prc_par_type
    logical :: l_elevation_corr
    real(wp) :: z_sur_crit_fele
    real(wp) :: z_sur_high_fele
    real(wp) :: dT_rain_snow
    real(wp) :: dP_dT
    logical :: l_slope_effect
    integer :: iwind_synoptic
    real(wp) :: topo_filter_width
    real(wp) :: wind_factor
    real(wp) :: wind_mod_factor
    real(wp) :: f_wind_max
  end type
  type(prc_par_type) :: prc_par

  ! surface parameters
  type surf_par_type

    integer :: i_f_ice
    real(wp) :: h_ice_crit
    real(wp) :: c_fice
    real(wp) :: z_sur_std_crit

    integer :: i_alb_ice
    integer :: i_alb_ice_margin

    real(wp) :: alb_soil
    real(wp) :: alb_firn
    real(wp) :: alb_ice_const
    real(wp) :: alb_ice_margin_const
    real(wp) :: alb_ice_clean 
    real(wp) :: alb_ice_dirty
    real(wp) :: tau_alb_ice_dirty
    real(wp) :: w_firn

    real(wp) :: alb_vis_dif_snow_new = 0.99
    real(wp) :: alb_nir_dif_snow_new = 0.65

    real(wp) :: d_alb_age_vis = 0.05_wp
    real(wp) :: d_alb_age_nir = 0.25_wp

    real(wp) :: z0m_ice
    real(wp) :: z0m_snow = 0.0024_wp

    real(wp) :: zm_to_zh = exp(-2._wp)  ! LandLad, Garratt 1992

    real(wp) :: emissivity_ice  = 0.99_wp 
    real(wp) :: emissivity_snow = 0.99_wp
    real(wp) :: emissivity_bare = 0.95_wp
  end type
  type(surf_par_type) :: surf_par

  integer :: ifake_atm_hires
  character (len=256) :: fake_atm_hires_const_file, fake_atm_hires_var_file
  integer :: prc_forcing
  integer :: wind_forcing

contains

    subroutine smb_params_init

    implicit none


    ! time step
    dt = dt_smb
    rdt = 1._wp/dt

    ! read smb parameter file
    call smb_par_load(trim(out_dir)//"/smb_par.nml")


    return

    end subroutine smb_params_init


subroutine smb_par_load(filename)

    implicit none

    character (len=*) :: filename

    real(wp) :: h_firn
    real(wp), parameter :: rho_firn = 600._wp   ! kg/m3


    ! Read parameters from file
    write(*,*) "smb parameters ==========="
    call nml_read(filename,"smb_par","i_ablation",pdd_par%i_ablation)
    call nml_read(filename,"smb_par","s_stat ",pdd_par%s_stat )
    call nml_read(filename,"smb_par","beta1  ",pdd_par%beta1)
    call nml_read(filename,"smb_par","beta2  ",pdd_par%beta2)
    call nml_read(filename,"smb_par","Pmax   ",pdd_par%Pmax   )
    call nml_read(filename,"smb_par","mu     ",pdd_par%mu   )
    pdd_par%beta1  = pdd_par%beta1 *(0.001_wp/86400.0_wp)*(rho_w/rho_i) ! (mm WE)/(d*deg C) --> (m IE)/(s*deg C)
    pdd_par%beta2  = pdd_par%beta2 *(0.001_wp/86400.0_wp)*(rho_w/rho_i) ! (mm WE)/(d*deg C) --> (m IE)/(s*deg C)
    pdd_par%mu     = pdd_par%mu    *(1000.0_wp*86400.0_wp)*(rho_i/rho_w) ! (d*deg C)/(mm WE) --> (s*deg C)/(m IE)

    call nml_read(filename,"smb_par","i_gamma ",i_gamma )
    call nml_read(filename,"smb_par","gamma ",gamma )
    gamma = gamma*1e-3_wp   ! k/m

    call nml_read(filename,"smb_par","i_z_sur_eff",i_z_sur_eff)
    call nml_read(filename,"smb_par","alpha_zstd",alpha_zstd )

    call nml_read(filename,"smb_par","i_smb",i_smb)

    call nml_read(filename,"smb_par","l_maxice",l_maxice)
    call nml_read(filename,"smb_par","mask_maxice_file",mask_maxice_file)
    call nml_read(filename,"smb_par","smb_maxice",smb_maxice)

    call nml_read(filename,"smb_par","smb_crit_mask",smb_crit_mask)
    call nml_read(filename,"smb_par","n_smb_mask_ext",n_smb_mask_ext)

    call nml_read(filename,"smb_par","l_regional_climate_forcing",l_regional_climate_forcing)
    call nml_read(filename,"smb_par","map_method",map_method)
    call nml_read(filename,"smb_par","filt_sigma",filt_sigma)

    call nml_read(filename,"smb_par","l_smb_bias_corr",l_smb_bias_corr)
    call nml_read(filename,"smb_par","smb_ref_file",smb_ref_file)
    call nml_read(filename,"smb_par","smb_cx_ref_file",smb_cx_ref_file)
    call nml_read(filename,"smb_par","year_ini_smb_ref",year_ini_smb_ref)
    call nml_read(filename,"smb_par","year_end_smb_ref",year_end_smb_ref)
    year_ini_smb_ref = year_ini_smb_ref - 2000
    year_end_smb_ref = year_end_smb_ref - 2000

    call nml_read(filename,"smb_par","i_t2m_bias_corr",i_t2m_bias_corr)
    call nml_read(filename,"smb_par","i_prc_bias_corr",i_prc_bias_corr)
    call nml_read(filename,"smb_par","bias_corr_file",bias_corr_file)
    call nml_read(filename,"smb_par","t2m_bias_scale_fac",t2m_bias_scale_fac)
    call nml_read(filename,"smb_par","t2m_bias_corr_uniform",t2m_bias_corr_uniform)
    call nml_read(filename,"smb_par","dLW_dT_fac",dLW_dT_fac)

    call nml_read(filename,"smb_par","d_cld",d_cld)

    call nml_read(filename,"smb_par","wind_ele_fac",wind_ele_fac)

    call nml_read(filename,"smb_par","l_Tvar_ann",l_Tvar_ann)
    call nml_read(filename,"smb_par","Tvar_ann_amp",Tvar_ann_amp)
    call nml_read(filename,"smb_par","Tvar_ann_period",Tvar_ann_period)
    Tvar_ann_period = nint(real(Tvar_ann_period,wp)/real(n_year_smb,wp))*n_year_smb
    call nml_read(filename,"smb_par","l_Tvar_day",l_Tvar_day)
    call nml_read(filename,"smb_par","Tvar_day_amp",Tvar_day_amp)
    call nml_read(filename,"smb_par","Tvar_day_period",Tvar_day_period)

    call nml_read(filename,"smb_par","nday_update_climate" ,nday_update_climate )

    call nml_read(filename,"smb_par","l_diurnal_cycle" ,l_diurnal_cycle )
    call nml_read(filename,"smb_par","tstd_scale" ,tstd_scale )

    call nml_read(filename,"smb_par","l_neutral" ,l_neutral)
    call nml_read(filename,"smb_par","l_dew" ,l_dew)

    call nml_read(filename,"smb_par","rho_snow" ,snow_par%rho )
    call nml_read(filename,"smb_par","lambda_snow" ,snow_par%lambda)
    call nml_read(filename,"smb_par","lsnow_aging" ,snow_par%lsnow_aging )
    call nml_read(filename,"smb_par","lsnow_dust" ,snow_par%lsnow_dust )
    call nml_read(filename,"smb_par","f_age_t" ,snow_par%f_age_t)
    call nml_read(filename,"smb_par","dT_age" ,snow_par%dT_age)
    call nml_read(filename,"smb_par","snow_0" ,snow_par%snow_0 )
    snow_par%snow_0 = snow_par%snow_0/sec_day
    call nml_read(filename,"smb_par","snow_1" ,snow_par%snow_1 )
    call nml_read(filename,"smb_par","snow_grain_fresh" ,snow_par%snow_grain_fresh )
    call nml_read(filename,"smb_par","isnow_albedo" ,snow_par%isnow_albedo )
    call nml_read(filename,"smb_par","dalb_snow_vis" ,snow_par%dalb_snow_vis )
    call nml_read(filename,"smb_par","dalb_snow_nir" ,snow_par%dalb_snow_nir )
    call nml_read(filename,"smb_par","k_sigma_orog" ,snow_par%k_sigma_orog )
    call nml_read(filename,"smb_par","sigma_orog_crit" ,snow_par%sigma_orog_crit )
    call nml_read(filename,"smb_par","c_fsnow" ,snow_par%c_fsnow)
    call nml_read(filename,"smb_par","c_fsnow_orog" ,snow_par%c_fsnow_orog)
    call nml_read(filename,"smb_par","l_fsnow_orog" ,snow_par%l_fsnow_orog)
    call nml_read(filename,"smb_par","w_snow_max" ,snow_par%w_snow_max )
    call nml_read(filename,"smb_par","w_snow_dust" ,snow_par%w_snow_dust )
    call nml_read(filename,"smb_par","dust_con_scale" ,snow_par%dust_con_scale)
    call nml_read(filename,"smb_par","i_rfz" ,snow_par%i_rfz)
    call nml_read(filename,"smb_par","porosity" ,snow_par%porosity)
    call nml_read(filename,"smb_par","f_rfz_max" ,snow_par%f_rfz_max )
    call nml_read(filename,"smb_par","f_rfz_to_snow_max" ,snow_par%f_rfz_to_snow_max )
    call nml_read(filename,"smb_par","wsnow_crit_rfz" ,snow_par%wsnow_crit_rfz)

    call nml_read(filename,"smb_par","i_f_ice" ,surf_par%i_f_ice)
    call nml_read(filename,"smb_par","h_ice_crit" ,surf_par%h_ice_crit)
    call nml_read(filename,"smb_par","c_fice" ,surf_par%c_fice)
    call nml_read(filename,"smb_par","z_sur_std_crit" ,surf_par%z_sur_std_crit)
    call nml_read(filename,"smb_par","z0m_ice" ,surf_par%z0m_ice)
    call nml_read(filename,"smb_par","i_alb_ice" ,surf_par%i_alb_ice)
    call nml_read(filename,"smb_par","i_alb_ice_margin" ,surf_par%i_alb_ice_margin)
    call nml_read(filename,"smb_par","alb_soil" ,surf_par%alb_soil)
    call nml_read(filename,"smb_par","alb_firn" ,surf_par%alb_firn)
    call nml_read(filename,"smb_par","alb_ice_const" ,surf_par%alb_ice_const)
    call nml_read(filename,"smb_par","alb_ice_margin_const" ,surf_par%alb_ice_margin_const)
    call nml_read(filename,"smb_par","alb_ice_clean" ,surf_par%alb_ice_clean)
    call nml_read(filename,"smb_par","alb_ice_dirty" ,surf_par%alb_ice_dirty)
    call nml_read(filename,"smb_par","tau_alb_ice_dirty" ,surf_par%tau_alb_ice_dirty)
    surf_par%tau_alb_ice_dirty = surf_par%tau_alb_ice_dirty * sec_year  ! yr -> s
    call nml_read(filename,"smb_par","h_firn" ,h_firn)
    surf_par%w_firn = h_firn*rho_firn

    call nml_read(filename,"smb_par","l_slope_effect" ,prc_par%l_slope_effect)
    call nml_read(filename,"smb_par","l_elevation_corr" ,prc_par%l_elevation_corr)
    call nml_read(filename,"smb_par","z_sur_crit_fele" ,prc_par%z_sur_crit_fele)
    call nml_read(filename,"smb_par","z_sur_high_fele" ,prc_par%z_sur_high_fele)
    call nml_read(filename,"smb_par","dT_rain_snow" ,prc_par%dT_rain_snow)
    call nml_read(filename,"smb_par","dP_dT" ,prc_par%dP_dT)
    call nml_read(filename,"smb_par","iwind_synoptic" ,prc_par%iwind_synoptic )
    call nml_read(filename,"smb_par","topo_filter_width" ,prc_par%topo_filter_width )
    call nml_read(filename,"smb_par","wind_factor" ,prc_par%wind_factor )
    call nml_read(filename,"smb_par","wind_mod_factor" ,prc_par%wind_mod_factor )
    call nml_read(filename,"smb_par","f_wind_max" ,prc_par%f_wind_max)

    call nml_read(filename,"smb_par","ifake_atm_hires",ifake_atm_hires)
    call nml_read(filename,"smb_par","fake_atm_hires_const_file",fake_atm_hires_const_file)
    call nml_read(filename,"smb_par","fake_atm_hires_var_file",fake_atm_hires_var_file)
    call nml_read(filename,"smb_par","prc_forcing",prc_forcing)
    call nml_read(filename,"smb_par","wind_forcing",wind_forcing)

    call nml_read(filename,"smb_par","l_write_timer",l_write_timer)

    call nml_read(filename,"smb_par","l_monthly_output",l_monthly_output)
    call nml_read(filename,"smb_par","l_daily_output",l_daily_output)

   return

end subroutine smb_par_load


end module smb_params
