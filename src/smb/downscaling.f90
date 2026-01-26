!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d o w n s c a l i n g _ m o d 
!
!  Purpose : downscaling of climate fields
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
module downscaling_mod

  use precision, only : wp
  use constants, only : pi, T0, frac_vu
  use smb_params, only : prc_par, wind_ele_fac, dLW_dT_fac, d_cld

  implicit none

  private
  public :: wind_downscaling, prc_downscaling, rad_downscaling

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  w i n d _ d o w n s c a l i n g
  !   Purpose    :  downscaling of wind
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine wind_downscaling(u700_i,v700_i,wind_i,z_sur,z_sur_i, &
                             u700,v700, wind)

  implicit none

  real(wp), intent(in) :: u700_i, v700_i
  real(wp), intent(in) :: wind_i
  real(wp), intent(in) :: z_sur, z_sur_i

  real(wp), intent(out) :: u700, v700, wind


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! effective wind used to compute surface energy fluxes

  ! increase with elevation
  wind = wind_i + wind_ele_fac*max(z_sur-z_sur_i,0._wp)  ! todo, check factor with CORDEX data

  ! +++++++++++++++++++++++++++++++++++++++++++++++
  ! effective boundary-layer wind for precipitation downscaling

  u700 = u700_i  ! 700 hPa
  v700 = v700_i  ! 700 hPa

  return

  end subroutine wind_downscaling


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p r c _ d o w n s c a l i n g
  !   Purpose    :  downscaling of precipitation 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine prc_downscaling(t2m,prc_i,prc_bias_i,u700,v700,wind,dz_dx_sur,dz_dy_sur,dz_sur,f_ele, &
                            snow,rain,prc,f_wind)

  implicit none

  real(wp), intent(in) :: t2m, prc_i, prc_bias_i, u700, v700, wind
  real(wp), intent(in) :: dz_dx_sur, dz_dy_sur, dz_sur, f_ele
  real(wp), intent(out) :: snow, rain, prc, f_wind

  real(wp) :: frsnw, w


  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Temperature dependence for snowfall fraction of precipitation

  if (t2m.lt.T0-prc_par%dT_rain_snow) frsnw=1._wp
  if (t2m.gt.T0+prc_par%dT_rain_snow) frsnw=0._wp
  if (t2m.ge.T0-prc_par%dT_rain_snow .and. t2m.le.T0+prc_par%dT_rain_snow) frsnw=(T0+prc_par%dT_rain_snow-t2m)/(2._wp*prc_par%dT_rain_snow)
  
  ! ++++++++++++++++++++++++
  ! wind slope effect

  if (prc_par%l_slope_effect) then

    if (prc_par%iwind_synoptic.eq.1) then
      w = u700*dz_dx_sur + v700*dz_dy_sur + prc_par%wind_mod_factor*wind*dz_sur
    else if (prc_par%iwind_synoptic.eq.2) then
      w = u700*dz_dx_sur + v700*dz_dy_sur + prc_par%wind_mod_factor*(wind-sqrt(u700**2+v700**2))*dz_sur
    else if (prc_par%iwind_synoptic.eq.3) then
      w = prc_par%wind_mod_factor*wind*dz_sur
    else if (prc_par%iwind_synoptic.eq.4) then
      w = prc_par%wind_mod_factor*(wind-sqrt(u700**2+v700**2))*dz_sur
    endif
    if (dz_sur.gt.1.e-6_wp) then
      f_wind = 1._wp+prc_par%wind_factor*w
    else 
      f_wind = 1._wp
    endif 
    ! limit wind factor, limits derived from CORDEX model data
    f_wind = min(f_wind,prc_par%f_wind_max)    
    f_wind = max(f_wind,1._wp/prc_par%f_wind_max) 

  else

    f_wind = 1._wp

  endif
      
  ! +++++++++++++++++++++++++++++++++++++++++++++++
  ! final corrected precipitation and snowfall rate

  prc = prc_i/prc_bias_i * f_ele * f_wind  ! kg/m2/s

  snow = prc * frsnw
  rain = prc - snow 


  return

  end subroutine prc_downscaling


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r a d _ d o w n s c a l i n g
  !   Purpose    :  downscaling of sw and lw radiation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rad_downscaling(z_sur, z_sur_i, t2m, t2m_bias_i, &
      alb_vis_dir_i, alb_nir_dir_i, alb_vis_dif_i, alb_nir_dif_i, &
      swd_sur_vis_dir_i, swd_sur_nir_dir_i, swd_sur_vis_dif_i, swd_sur_nir_dif_i, &
      dswd_dalb_vis_dir_i, dswd_dalb_nir_dir_i, dswd_dalb_vis_dif_i, dswd_dalb_nir_dif_i, &
      dswd_dz_nir_dir_i, dswd_dz_nir_dif_i, &
      swd_toa_i, swd_toa_min_i, cld_i, &
      alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif, &
      lwdown_i, gam_lw_i, &
      cld, albedo, swnet, swnet_min, swdown, lwdown, i,j)

  implicit none

  real(wp), intent(in) :: z_sur   !! surface elevation [m]
  real(wp), intent(in) :: z_sur_i !! mean surface elevation for atmosphere interpolated to ice sheet grid [m]
  real(wp), intent(in) :: t2m !! 2m temperature [K]
  real(wp), intent(in) :: t2m_bias_i !! 2m temperature bias at present-day interpolated to ice sheet grid [K]
  real(wp), intent(in) :: alb_vis_dir_i    !! surface visible clear sky albedo []
  real(wp), intent(in) :: alb_nir_dir_i    !! surface near-infrared clear sky albedo []
  real(wp), intent(in) :: alb_vis_dif_i    !! surface visible cloudy sky albedo []
  real(wp), intent(in) :: alb_nir_dif_i    !! surface near-infrared cloudy sky albedo []
  real(wp), intent(in) :: swd_sur_vis_dir_i    !! downward surface visible clear sky radiation at elevation z_sur_i [W/m2]
  real(wp), intent(in) :: swd_sur_nir_dir_i    !! downward surface near-infrared clear sky radiation at elevation z_sur_i [W/m2]
  real(wp), intent(in) :: swd_sur_vis_dif_i    !! downward surface visible cloudy sky radiation at elevation z_sur_i [W/m2]
  real(wp), intent(in) :: swd_sur_nir_dif_i    !! downward surface near-infrared cloudy sky radiation at elevation z_sur_i [W/m2]
  real(wp), intent(in) :: dswd_dalb_vis_dir_i    
  real(wp), intent(in) :: dswd_dalb_nir_dir_i    
  real(wp), intent(in) :: dswd_dalb_vis_dif_i    
  real(wp), intent(in) :: dswd_dalb_nir_dif_i    
  real(wp), intent(in) :: dswd_dz_nir_dir_i    
  real(wp), intent(in) :: dswd_dz_nir_dif_i    
  real(wp), intent(in) :: swd_toa_i    !! downward radiation at top of atmosphere [W/m2]
  real(wp), intent(in) :: swd_toa_min_i    !! minimum diurnal downward radiation at top of atmosphere [W/m2]
  real(wp), intent(in) :: alb_vis_dir    !! surface visible clear sky albedo []
  real(wp), intent(in) :: alb_nir_dir    !! surface near-infrared clear sky albedo []
  real(wp), intent(in) :: alb_vis_dif    !! surface visible cloudy sky albedo []
  real(wp), intent(in) :: alb_nir_dif    !! surface near-infrared cloudy sky albedo []
  real(wp), intent(in) :: cld_i   !! cloud cover fraction [/]
  real(wp), intent(in) :: lwdown_i  !! downward longwave radiation [W/m2] 
  real(wp), intent(in) :: gam_lw_i  !! downward longwave radiation [W/m2] 

  real(wp), intent(out) :: cld    !! cloud cover fraction [/]
  real(wp), intent(out) :: albedo    !! surface albedo [/]
  real(wp), intent(out) :: swnet    !! net surface shortwave radiation at ice sheet elevation [W/m2]
  real(wp), intent(out) :: swnet_min    !! minimum diurnal net surface shortwave radiation at ice sheet elevation [W/m2]
  real(wp), intent(out) :: swdown   !! downward surface shortwave radiation at ice sheet elevation [W/m2]
  real(wp), intent(out) :: lwdown   !! downward surface longwave radiation at ice sheet elevation [W/m2]
      
  integer :: i,j
  real(wp) :: swd_sur_vis_dir, swd_sur_nir_dir, swd_sur_vis_dif, swd_sur_nir_dif


  ! +++++++++++++++++++++++
  ! cloud fraction

  cld = min(1._wp, cld_i + d_cld)

  ! +++++++++++++++++++++++
  ! shortwave radiation

  if (swd_toa_i.gt.0.1_wp) then

    swd_sur_vis_dir = swd_sur_vis_dir_i + dswd_dalb_vis_dir_i*(alb_vis_dir-alb_vis_dir_i) 
    swd_sur_nir_dir = swd_sur_nir_dir_i + dswd_dalb_nir_dir_i*(alb_nir_dir-alb_nir_dir_i) + dswd_dz_nir_dir_i*(z_sur-z_sur_i)
    swd_sur_vis_dif = swd_sur_vis_dif_i + dswd_dalb_vis_dif_i*(alb_vis_dif-alb_vis_dif_i) 
    swd_sur_nir_dif = swd_sur_nir_dif_i + dswd_dalb_nir_dif_i*(alb_nir_dif-alb_nir_dif_i) + dswd_dz_nir_dif_i*(z_sur-z_sur_i)

    swdown = (1._wp-cld)*(frac_vu*swd_sur_vis_dir+(1._wp-frac_vu)*swd_sur_nir_dir) &
           + cld * (frac_vu*swd_sur_vis_dif+(1._wp-frac_vu)*swd_sur_nir_dif)

    swnet  = (1._wp-cld)*(frac_vu*swd_sur_vis_dir*(1._wp-alb_vis_dir)+(1._wp-frac_vu)*swd_sur_nir_dir*(1._wp-alb_nir_dir)) &
           + cld * (frac_vu*swd_sur_vis_dif*(1._wp-alb_vis_dif)+(1._wp-frac_vu)*swd_sur_nir_dif*(1._wp-alb_nir_dif))

    ! estimate of minimum net surface solar radiation for amplitude of diurnal cycle
    swnet_min = swnet * swd_toa_min_i/ swd_toa_i

    !if (i.eq.80 .and. j.eq.170) then
    !  print *
    !  print *,'z_sur',z_sur
    !  print *,'INTERPOLATED: SW TOA,SW dir,SW diff',swd_toa_i,swd_sur_dir_i,swd_sur_dif_i
    !  print *,'DOWNSCALED  : SW TOA,SW dir,SW diff',swd_toa_i,swd_dir,swd_dif
    !  print *,'DOWNSCALED  : SW',swdown
    !  print *,'cld',cld
    !endif

    ! Effective surface albedo, diagnostic only
    if (swdown.gt.0._wp) then
      albedo = 1._wp-swnet/swdown
    else
      albedo = 0._wp
    endif
    albedo = max(albedo,0._wp)
    albedo = min(albedo,1._wp)

  else
    ! polar night

    swdown = 0._wp
    swnet = 0._wp
    swnet_min = 0._wp
    albedo = 0._wp

  endif  


  ! ++++++++++++++++++++++++++++
  ! Longwave radiation

  ! Vertical interpolation of down longwave radiation at the surface
  lwdown = lwdown_i+gam_lw_i*(z_sur-z_sur_i)
  ! correction for near surface temperature bias
  ! based on quadratic fit of LWdown vs T2m in ERA-Interim 
  lwdown = lwdown - dLW_dT_fac*(-5._wp+2._wp*0.017_wp*t2m)*t2m_bias_i


  return

  end subroutine rad_downscaling

end module downscaling_mod

