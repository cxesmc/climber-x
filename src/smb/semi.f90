!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s e m i _ m 
!
!  Purpose : Main SEMIX surface energy and mass balance model 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit, Reihard Calov and Andrey Ganopolski
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
module semi_m

use precision, only : wp
use constants, only : fqsat, q_sat_i
use smb_params, only : dt, i_gamma, gamma, l_regional_climate_forcing, prc_par, snow_par
use downscaling_mod, only : wind_downscaling, prc_downscaling, rad_downscaling
use smb_surface_par, only : surface_albedo, resistance
use smb_ebal_mod, only : ebal, update_tskin 
use smb_temp_mod, only : smb_temp 
use snow_mod, only : snow_update

private
public :: semi

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s e m i
  !   Purpose    :  surface energy and mass balance interface
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine semi(i, j, f_ice, alb_ice, z_sur, z_sur_i, z_sur_std, dz_dx_sur, dz_dy_sur, dz_sur, f_ele, & ! in
      tam_i, t2m_bias_i, dTvar, gam_i, tstd_i, ram_i, pressure, u700_i, v700_i, wind_i, prc_i, prc_bias_i, &    ! in
      alb_vis_dir_i, alb_nir_dir_i, alb_vis_dif_i, alb_nir_dif_i, &     ! in
      swd_sur_vis_dir_i, swd_sur_nir_dir_i, swd_sur_vis_dif_i, swd_sur_nir_dif_i, &     ! in
      dswd_dalb_vis_dir_i, dswd_dalb_nir_dir_i, dswd_dalb_vis_dif_i, dswd_dalb_nir_dif_i, & ! in
      dswd_dz_nir_dir_i, dswd_dz_nir_dif_i, dust_i, coszm_i, &  ! in
      swd_toa_i, swd_toa_min_i, cld_i, lwdown_i, gam_lw_i, & ! in
      tam, t2m, t_skin, t_skin_old, t_skin_amp, t_prof, t_prof_old, & ! out
      q2m, qsat, dqsatdT, r_a, &   ! out 
      u700, v700, wind, snow, rain, prc, f_wind, &   ! out
      mask_snow, f_snow, h_snow, w_snow, w_snow_old, w_snow_max, &    ! out
      snow_grain, dust_con, & ! out
      alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif, &   ! out
      alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif, &   ! out
      dt_snowfree, alb_bg, albedo, cld, swnet, swnet_min, swdown, lwdown, &  ! out
      flx_g, dflxg_dT, flx_melt, flx_sh, flx_lwu, flx_lh, evp, &   ! out
      num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &  ! out
      f_sh, f_e, f_lh, f_lw, & ! out
      snowmelt, icemelt, refreezing, refreezing_sum, f_rfz_to_snow)   ! out

  implicit none

  integer, intent(in)  :: i, j
  real(wp), intent(in) :: f_ice
  real(wp), intent(in) :: z_sur
  real(wp), intent(in) :: z_sur_i
  real(wp), intent(in) :: z_sur_std
  real(wp), intent(in) :: dz_dx_sur
  real(wp), intent(in) :: dz_dy_sur
  real(wp), intent(in) :: dz_sur
  real(wp), intent(in) :: f_ele
  real(wp), intent(in) :: tam_i
  real(wp), intent(in) :: t2m_bias_i
  real(wp), intent(in) :: dTvar
  real(wp), intent(in) :: gam_i
  real(wp), intent(in) :: tstd_i
  real(wp), intent(in) :: ram_i
  real(wp), intent(in) :: pressure
  real(wp), intent(in) :: u700_i
  real(wp), intent(in) :: v700_i
  real(wp), intent(in) :: wind_i
  real(wp), intent(in) :: prc_i
  real(wp), intent(in) :: prc_bias_i
  real(wp), intent(in) :: alb_vis_dir_i
  real(wp), intent(in) :: alb_nir_dir_i
  real(wp), intent(in) :: alb_vis_dif_i
  real(wp), intent(in) :: alb_nir_dif_i
  real(wp), intent(in) :: swd_sur_vis_dir_i
  real(wp), intent(in) :: swd_sur_nir_dir_i
  real(wp), intent(in) :: swd_sur_vis_dif_i
  real(wp), intent(in) :: swd_sur_nir_dif_i
  real(wp), intent(in) :: dswd_dalb_vis_dir_i
  real(wp), intent(in) :: dswd_dalb_nir_dir_i
  real(wp), intent(in) :: dswd_dalb_vis_dif_i
  real(wp), intent(in) :: dswd_dalb_nir_dif_i
  real(wp), intent(in) :: dswd_dz_nir_dir_i
  real(wp), intent(in) :: dswd_dz_nir_dif_i
  real(wp), intent(in) :: dust_i
  real(wp), intent(in) :: coszm_i
  real(wp), intent(in) :: swd_toa_i
  real(wp), intent(in) :: swd_toa_min_i
  real(wp), intent(in) :: cld_i
  real(wp), intent(in) :: lwdown_i
  real(wp), intent(in) :: gam_lw_i

  real(wp), intent(inout) :: tam
  real(wp), intent(inout) :: t2m
  real(wp), intent(inout) :: t_skin
  real(wp), intent(inout) :: t_skin_old
  real(wp), intent(inout) :: t_skin_amp
  real(wp), intent(inout) :: t_prof(0:)
  real(wp), intent(inout) :: t_prof_old(0:)
  real(wp), intent(inout) :: q2m
  real(wp), intent(inout) :: qsat
  real(wp), intent(inout) :: dqsatdT
  real(wp), intent(inout) :: r_a 
  real(wp), intent(inout) :: u700
  real(wp), intent(inout) :: v700
  real(wp), intent(inout) :: wind
  real(wp), intent(inout) :: snow
  real(wp), intent(inout) :: rain
  real(wp), intent(inout) :: prc
  real(wp), intent(inout) :: f_wind
  integer,  intent(inout)  :: mask_snow
  real(wp), intent(inout) :: f_snow
  real(wp), intent(inout) :: h_snow
  real(wp), intent(inout) :: w_snow
  real(wp), intent(inout) :: w_snow_old
  real(wp), intent(inout) :: w_snow_max
  real(wp), intent(inout) :: snow_grain
  real(wp), intent(inout) :: dust_con
  real(wp), intent(inout) :: alb_vis_dir
  real(wp), intent(inout) :: alb_nir_dir
  real(wp), intent(inout) :: alb_vis_dif
  real(wp), intent(inout) :: alb_nir_dif
  real(wp), intent(inout) :: alb_snow_vis_dir
  real(wp), intent(inout) :: alb_snow_nir_dir
  real(wp), intent(inout) :: alb_snow_vis_dif
  real(wp), intent(inout) :: alb_snow_nir_dif
  real(wp), intent(inout) :: alb_ice
  real(wp), intent(inout) :: dt_snowfree 
  real(wp), intent(inout) :: alb_bg
  real(wp), intent(inout) :: albedo
  real(wp), intent(inout) :: cld
  real(wp), intent(inout) :: swnet
  real(wp), intent(inout) :: swnet_min
  real(wp), intent(inout) :: swdown
  real(wp), intent(inout) :: lwdown
  real(wp), intent(inout) :: flx_g
  real(wp), intent(inout) :: dflxg_dT
  real(wp), intent(inout) :: flx_melt
  real(wp), intent(inout) :: flx_sh
  real(wp), intent(inout) :: flx_lwu
  real(wp), intent(inout) :: flx_lh
  real(wp), intent(inout) :: evp
  real(wp), intent(inout) :: num_lh
  real(wp), intent(inout) :: num_sh
  real(wp), intent(inout) :: num_sw
  real(wp), intent(inout) :: num_lw
  real(wp), intent(inout) :: denom_lh
  real(wp), intent(inout) :: denom_sh
  real(wp), intent(inout) :: denom_lw
  real(wp), intent(inout) :: f_sh
  real(wp), intent(inout) :: f_e
  real(wp), intent(inout) :: f_lh
  real(wp), intent(inout) :: f_lw
  real(wp), intent(inout) :: snowmelt
  real(wp), intent(inout) :: icemelt
  real(wp), intent(inout) :: refreezing
  real(wp), intent(inout) :: refreezing_sum
  real(wp), intent(inout) :: f_rfz_to_snow

  real(wp) :: ram, qam, qsat_skin, r2m, r_skin


  ! +++++++++++++++++++++++++++++++++
  ! downscaling

  if (.not.l_regional_climate_forcing) then
    ! free atmosphere temperature at ice sheet elevation
    if (i_gamma.eq.1) then
      tam = tam_i + gamma*(z_sur_i-z_sur) - t2m_bias_i 
    else if (i_gamma.eq.2) then
      tam = tam_i + gam_i*(z_sur_i-z_sur) - t2m_bias_i 
    endif
    ! add artificial interannual temperature variability 
    tam = tam + dTvar
    ! surface air temperature, average of atmospheric temperature and skin temperature
    t2m = 0.5_wp*(tam+t_skin) 
    ! surface air specific humidity, relative humidity as average between ram and rskina 
    qsat_skin = q_sat_i(t_skin,pressure)
    ram = ram_i
    qam = ram*fqsat(tam,pressure)
    r_skin = min(qsat_skin,qam)/qsat_skin
    r2m = 0.5_wp*(ram+r_skin) 
    q2m = r2m*q_sat_i(t2m,pressure)
    ! wind downscaling 
    call wind_downscaling(u700_i, v700_i, wind_i, z_sur, z_sur_i, & ! in
      u700, v700, wind) ! out
    ! precipitation downscaling
    call prc_downscaling(t2m, prc_i, prc_bias_i, u700, v700, wind_i, dz_dx_sur, dz_dy_sur, dz_sur, f_ele, &    ! in
      snow, rain, prc, f_wind)  ! out
  endif

  ! surface albedo (needed for radiation downscaling)
  call surface_albedo(z_sur_std, f_ice, alb_ice, h_snow, w_snow, w_snow_max, &
    snow, &
    t_skin, dust_i, coszm_i, &
    snow_grain, dust_con, f_snow, dt_snowfree, alb_bg, &
    alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif, &
    alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)

  if (.not.l_regional_climate_forcing) then
    ! radiation downscaling
    call rad_downscaling(z_sur, z_sur_i, t2m, t2m_bias_i, &
      alb_vis_dir_i, alb_nir_dir_i, alb_vis_dif_i, alb_nir_dif_i, &
      swd_sur_vis_dir_i, swd_sur_nir_dir_i, swd_sur_vis_dif_i, swd_sur_nir_dif_i, &
      dswd_dalb_vis_dir_i, dswd_dalb_nir_dir_i, dswd_dalb_vis_dif_i, dswd_dalb_nir_dif_i, &
      dswd_dz_nir_dir_i, dswd_dz_nir_dif_i, &
      swd_toa_i, swd_toa_min_i, cld_i, &
      alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif, &
      lwdown_i, gam_lw_i, &
      cld, albedo, swnet, swnet_min, swdown, lwdown,i,j)  ! out
  else
    swnet = swdown*(1._wp-albedo) 
    if (swd_toa_i>0._wp) then
      swnet_min = swnet * swd_toa_min_i/swd_toa_i
    else
      swnet_min = 0._wp
    endif
  endif

  ! +++++++++++++++++++++++++++++
  ! aerodynamic resistance 
  call resistance(h_snow, t2m, t_skin, wind, &  ! in
    r_a)    ! out

  ! +++++++++++++++++++++++++++++
  ! add snowfall to snow layer, so that it can be melted in smb_temp
  w_snow_old = w_snow
  w_snow = w_snow + snow*dt
  h_snow = w_snow / snow_par%rho  ! m

  ! +++++++++++++++++++++++++++++
  ! surface energy balance
  call ebal(mask_snow, h_snow, & ! in
    t_prof(:), t2m, tstd_i, q2m, pressure, &                            ! in
    swnet, swnet_min, lwdown, r_a, &                         ! in
    t_skin, &                                                           ! inout
    t_skin_old, t_skin_amp, flx_g, dflxg_dT, flx_melt, &                ! out
    num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &     ! out
    f_sh, f_e, f_lh, f_lw, qsat, dqsatdT,i,j)                           ! out

  ! +++++++++++++++++++++++++++++++
  ! snow/ice/soil thermodynamics
  call smb_temp(mask_snow, h_snow, rain, &
    flx_g, dflxg_dT, flx_melt, refreezing_sum, &
    t_prof(:), w_snow, snowmelt, icemelt, refreezing, f_rfz_to_snow, &
    t_prof_old(:), i,j)

  refreezing_sum = refreezing_sum + refreezing*dt 
  if (refreezing.eq.0._wp) then
    refreezing_sum = refreezing_sum - snow*dt
    refreezing_sum = max(0._wp,refreezing_sum)
  endif
  
  ! +++++++++++++++++++++++++++++++
  ! update skin temperature
  call update_tskin(mask_snow, t_skin_old, dflxg_dT, t2m, q2m, &
    swnet, lwdown, &
    t_prof(:), t_prof_old(:), flx_melt, & 
    num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
    f_sh, f_e, f_lh, f_lw, qsat, dqsatdT, &
    t_skin, flx_g, flx_sh, flx_lwu, flx_lh, evp,i,j)

  ! +++++++++++++++++++++++++++++++
  ! update snow layer
  call snow_update(mask_snow, evp, &
    w_snow, w_snow_old, w_snow_max, &
    h_snow)


  return

end subroutine semi

end module semi_m
