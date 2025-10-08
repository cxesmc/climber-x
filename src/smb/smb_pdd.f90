!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ p d d _ m
!
!  Purpose : PDD model
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Reinhard Calov
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
module smb_pdd_m

use precision, only : wp
use constants, only : rho_w, rho_i, T0, pi
use timer, only : sec_year, sec_mon, time_eoy_smb
use smb_params, only : prc_par, pdd_par, gamma
use downscaling_mod, only : wind_downscaling, prc_downscaling

real(wp), parameter :: dt_simple = sec_mon 

real(wp), parameter :: inv_sqrt2pi = 1.0_wp/sqrt(2.0_wp*pi)
real(wp), parameter :: inv_sqrt2   = 1.0_wp/sqrt(2.0_wp)


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ p d d
  !   Purpose    :  surface mass balance scheme based on PDD
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine smb_pdd(z_sur, z_sur_i, dz_dx_sur, dz_dy_sur, dz_sur, f_ele, &   
          t2m_i, t2m_bias_i, u700_i, v700_i, wind_i, prc_i, prc_bias_i, &    
          t2m, u700, v700, wind, snow, rain, prc, f_wind, &   
          pdd, t2m_cum, snow_cum, rain_cum,&
          smb, melt, melt_star, runoff, t_ice)

  implicit none

  real(wp), intent(in) :: z_sur
  real(wp), intent(in) :: z_sur_i
  real(wp), intent(in) :: dz_dx_sur
  real(wp), intent(in) :: dz_dy_sur
  real(wp), intent(in) :: dz_sur
  real(wp), intent(in) :: f_ele
  real(wp), intent(in) :: t2m_i
  real(wp), intent(in) :: t2m_bias_i
  real(wp), intent(in) :: u700_i
  real(wp), intent(in) :: v700_i
  real(wp), intent(in) :: wind_i
  real(wp), intent(in) :: prc_i
  real(wp), intent(in) :: prc_bias_i
  real(wp), intent(out) :: t2m
  real(wp), intent(out) :: u700
  real(wp), intent(out) :: v700
  real(wp), intent(out) :: wind
  real(wp), intent(out) :: snow
  real(wp), intent(out) :: rain
  real(wp), intent(out) :: prc
  real(wp), intent(out) :: f_wind

  real(wp), intent(inout) :: pdd, t2m_cum, snow_cum, rain_cum
  real(wp), intent(out) :: smb, melt, melt_star, runoff, t_ice

  real(wp) :: t
  real(wp) :: snow_ann, rain_ann


  !-------------------------------------
  ! 2m temperature at ice sheet elevation
  t2m = t2m_i + gamma*(z_sur_i-z_sur) - t2m_bias_i

  !-------------------------------------
  ! wind downscaling 
  call wind_downscaling(u700_i, v700_i, wind_i, &  ! in
    z_sur, z_sur_i, &  ! in
    u700, v700, wind)   ! out

  !-------------------------------------
  ! precipitation downscaling
  call prc_downscaling(t2m, prc_i, prc_bias_i, u700, v700, wind, &    ! in
    dz_dx_sur, dz_dy_sur, dz_sur, f_ele, &   ! in
    snow, rain, prc, f_wind) ! out

  ! cumulate over the year
  snow_cum = snow_cum + snow * dt_simple  ! kg/m2
  rain_cum = rain_cum + rain * dt_simple  ! kg/m2

  ! PDD
  t = t2m - T0 ! degC
  pdd = pdd &   ! positive degree days   (deg C)
    + ( pdd_par%s_stat*inv_sqrt2pi*exp(-0.5_wp*(t/pdd_par%s_stat)**2) &
    + 0.5_wp*t * erfc(-t/pdd_par%s_stat*inv_sqrt2) ) &
    * dt_simple/sec_year 
  
  ! cumulate temperature
  t2m_cum = t2m_cum + t*dt_simple


  ! Formation rate of superimposed ice (melt_star), melt rate (melt)
  ! and runoff rate (runoff)

  if (time_eoy_smb) then

    snow_ann = snow_cum / sec_year / rho_i     ! kg/m2/yr * yr/s * m3/kg = m IE/s
    rain_ann = rain_cum / sec_year / rho_i     ! kg/m2/yr * yr/s * m3/kg = m IE/s

    if (pdd_par%i_ablation==1) then
      ! Ablation parameterized by positive-degree-day (PDD) method.
      ! Rainfall assumed to run off instantaneously.

      if ((pdd_par%beta1*pdd) <= (pdd_par%Pmax*snow_ann)) then
        melt_star = pdd_par%beta1*pdd  ! m IE/s/degC * degC = m IE/s
        melt      = 0.0_wp
        runoff    = melt+rain_ann
      else
        melt_star = pdd_par%Pmax*snow_ann
        melt      = pdd_par%beta2*(pdd-melt_star/pdd_par%beta1)
        runoff    = melt+rain_ann
      end if

    else if (pdd_par%i_ablation==2) then
      ! Ablation parameterized ! by positive-degree-day (PDD) method.
      ! Rainfall assumed to contribute to formation of superimposed ice.

      if ( rain<= (pdd_par%Pmax*snow_ann) ) then

        if ( (rain_ann+pdd_par%beta1*pdd) <= (pdd_par%Pmax*snow_ann) ) then
          melt_star = rain_ann+pdd_par%beta1*pdd
          melt      = 0.0_wp
          runoff    = melt
        else
          melt_star = pdd_par%Pmax*snow_ann
          melt      = pdd_par%beta2*(pdd-(melt_star-rain_ann)/pdd_par%beta1)
          runoff    = melt
        end if

      else

        melt_star = pdd_par%Pmax*snow_ann
        melt      = pdd_par%beta2*pdd
        runoff    = melt + rain_ann-pdd_par%Pmax*snow_ann

      end if

    endif

    smb = (snow_ann+rain_ann) - runoff     ! m IE/s
    smb = smb * rho_i     ! kg/m2/s

    ! Ice-surface temperature (10-m firn temperature),
    ! including empirical firn-warming correction due to
    ! refreezing meltwater when superimposed ice is formed
    if (melt_star >= melt) then
      t_ice = t2m_cum/sec_year + pdd_par%mu*(melt_star-melt)
    else
      t_ice = t2m_cum/sec_year
    endif

    ! convert for output
    melt      = melt*rho_i   ! kg/m2/s
    melt_star = melt_star*rho_i   ! kg/m2/s
    runoff    = runoff*rho_i   ! kg/m2/s

  endif

end subroutine smb_pdd

end module smb_pdd_m
